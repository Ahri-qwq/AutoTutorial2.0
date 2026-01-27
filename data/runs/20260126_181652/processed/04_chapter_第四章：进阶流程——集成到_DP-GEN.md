# 第四章：进阶流程——集成到 DP-GEN

在前面的章节中，我们已经掌握了如何单独运行 ABACUS 进行单点能计算和分子动力学模拟。然而，在构建机器学习势函数（Machine Learning Potential, MLP）的大规模生产环境中，手工提交任务、收集数据、转换格式的效率极低且容易出错。

本章将介绍如何将 ABACUS 作为高精度的**标注引擎（Labeling Engine）**嵌入到 **DP-GEN**（Deep Potential Generator）的并发主动学习工作流中。

---

## 4.1 DP-GEN 中的 ABACUS 接口配置

DP-GEN 的核心控制文件通常命名为 `param.json`。为了让 DP-GEN 识别并正确调用 ABACUS，我们需要在 `fp_params`（First Principles parameters）部分进行特定配置。

### 4.1.1 核心参数解析

以下是一个典型的 `param.json` 片段，展示了如何指定 ABACUS 为计算引擎：

```json
{
    "type_map": ["Si", "C"],
    "mass_map": [28.0855, 12.0107],
    "fp_style": "abacus",
    "fp_task_max": 20,
    "fp_task_min": 5,
    "fp_pp_path": "./PP_ORB",
    "fp_orb_path": "./PP_ORB",
    "fp_params": {
        "ecutwfc": 100,
        "scf_thr": 1e-6,
        "basis_type": "lcao",
        "smearing_method": "gauss",
        "smearing_sigma": 0.002,
        "mixing_type": "pulay",
        "mixing_beta": 0.3,
        "cal_force": 1,
        "cal_stress": 1
    }
}
```

**关键参数详解：**

1.  **`fp_style`: `"abacus"`**
    *   这是告诉 DP-GEN 此时的第一性原理计算引擎是 ABACUS。
2.  **`fp_pp_path` / `fp_orb_path`**
    *   指定赝势（Pseudopotentials）和轨道（Orbitals）文件所在的**绝对路径**或相对路径。DP-GEN 在分发任务时会将这些文件链接到计算目录中。
3.  **`fp_params` 字典**
    *   这里定义的参数会被直接写入 ABACUS 的 `INPUT` 文件中。
    *   **必须包含**：
        *   `cal_force`: `1` —— 计算原子受力（Force），这是训练势函数的核心数据。
        *   `cal_stress`: `1` —— 计算维里应力（Virial Stress）。**这是初学者最容易遗漏的参数**。如果缺失，训练过程中的 Virial Loss 将无法计算，导致模型对晶胞体积变化的预测能力极差。

### 4.1.2 迭代轮次（Stages）

在 `param.json` 的外层，通常还需要定义迭代流程：

```json
"model_devi_jobs": [
    {"sys_idx": [0], "temps": [300], "press": [1.0], "trj_freq": 10, "nsteps": 1000, "ensemble": "nvt", "_idx": "00"},
    {"sys_idx": [0], "temps": [600], "press": [1.0], "trj_freq": 10, "nsteps": 2000, "ensemble": "nvt", "_idx": "01"}
]
```

*   **注意**：这里的 `model_devi_jobs` 定义的是**探索（Exploration）**阶段的参数。通常这个阶段使用 LAMMPS 运行 DPMD 模型。
*   当模型发现“不准确”的构型后，DP-GEN 会自动将其筛选出来，生成 ABACUS 的输入文件（`STRU`, `INPUT`, `KPT` 等），进入 **02.fp (Labeling)** 阶段进行高精度计算。

---

## 4.2 目录结构与文件依赖

自动化流程对文件结构极其敏感。为了确保 ABACUS 在远程集群上能顺利运行，除了 DP-GEN 自动生成的 `STRU` 和 `INPUT` 外，你必须提前准备好静态依赖文件。

### 4.2.1 推荐的目录树

```text
Work_Dir/
├── param.json          # DP-GEN 主配置文件
├── machine.json        # 机器配置（定义如何提交 Slurm/PBS 任务）
├── PP_ORB/             # 存放物理势文件的目录
│   ├── Si_ONCV_PBE-1.0.upf
│   ├── Si_orb_6-31G_1.2mm.orb
│   ├── C_ONCV_PBE-1.0.upf
│   └── C_orb_6-31G_1.2mm.orb
└── init_data/          # 初始训练数据（可选）
    ├── ABACUS_MD_01/
    │   ├── INPUT
    │   ├── STRU
    │   ├── KPT
    │   └── ...
```

### 4.2.2 关键依赖文件说明

1.  **赝势与轨道文件 (`*.upf`, `*.orb`)**
    *   **来源**：必须与 `STRU` 文件中的元素名称严格对应。
    *   **配置**：在 `param.json` 中通过 `fp_pp_path` 指定目录。
    *   **注意**：ABACUS 的 LCAO 基组计算强依赖于 `.orb` 文件，请确保它们与 `.upf` 文件是配套生成的。

2.  **K点设置 (`KPT`)**
    *   **生成机制**：DP-GEN 通常允许在 `param.json` 的 `fp_params` 中指定 `kspacing`（K点间距），从而自动为每个结构生成 `KPT` 文件；或者你可以提供一个通用的模板。
    *   **重要性**：对于金属体系，K点采样密度直接影响能量和力的准确性。

3.  **机器配置 (`machine.json`)**
    *   该文件定义了 ABACUS 具体的执行命令。例如：
    ```json
    "command": "mpirun -np 32 abacus"
    ```
    *   **环境加载**：务必在 `machine.json` 的 `context_type` 或 `prepend_script` 中写入 `module load abacus` 或 `source /path/to/abacus_env.sh`，确保计算节点能找到 ABACUS 可执行程序。

---

## 附录：常见陷阱与调试指南

在将 ABACUS 集成到 DP-GEN 的过程中，以下四个陷阱占据了用户报错的 80%，请务必逐一核查。

### 1. 极易混淆点：DFT 采样 vs DPMD 推理

这是概念上最容易混淆的地方，直接导致计算毫无意义。

*   **场景 A：生成训练数据 (Labeling / FP step)**
    *   **目的**：利用量子力学计算真实的能量、力和应力，作为“标准答案”。
    *   **设置**：`INPUT` 文件中 **`esolver_type` 必须是 `ks`** (默认值) 或显式不写。
    *   **绝对禁止**：在此阶段设置 `esolver_type dp`。如果你这样做了，等于是在用一个未训练好的模型去教它自己（循环论证），产生的数据是垃圾数据。

*   **场景 B：使用模型进行探索 (Exploration / Model Devi step)**
    *   **目的**：利用现有的 DP 模型快速跑 MD，寻找高不确定性构型。
    *   **设置**：通常由 LAMMPS 完成。如果必须用 ABACUS 做推理，此时 `INPUT` 中才设置 `esolver_type dp`。

### 2. 物理量缺失 (Missing Physical Quantities)

*   **症状**：DP-GEN 流程跑通了，但 `dp train` 阶段报错，提示 Virial 维度不匹配或数据全为零。
*   **原因**：ABACUS 默认 **不计算** 应力张量。
*   **解决方案**：
    *   在 `param.json` 的 `fp_params` 中，**必须显式添加**：
        ```text
        cal_stress 1
        cal_force 1
        ```
    *   对于 LCAO 基组，计算应力会显著增加计算耗时，但这是训练 NPT 系综适用势函数的必要代价。

### 3. 数据格式陷阱 (`dpdata` 兼容性)

*   **问题**：ABACUS 版本迭代较快（如 3.2.x 到 3.5.x），输出日志格式可能微调，导致 `dpdata` 无法正则匹配数据。
*   **风险提示**：
    *   DP-GEN 依赖 `dpdata` 库来解析 ABACUS 的输出。
    *   如果遇到 `KeyError` 或解析失败，请优先检查 `dpdata` 是否为最新版。
    *   建议查阅 `dpdata` 官方文档中关于 `abacus/md` 或 `abacus/scf` 格式的说明，确认你的 ABACUS 输出文件名（如 `Running_MD` 或 `MD_dump`）是否符合预期。

### 4. 数据集划分策略

*   **操作建议**：虽然 DP-GEN 的 `00.train` 环节可以自动划分验证集，但在准备初始数据集（`init_data`）时，建议手动进行划分。
*   **最佳实践**：
    *   引用资料建议，将数据集划分为 **训练集 (Training)** 和 **验证集 (Validation)**。
    *   例如：对于一条包含 1000 帧的 AIMD 轨迹，随机抽取 100 帧放入验证集目录，其余 900 帧作为训练集。这能让你在模型训练初期就真实地评估模型的泛化能力，防止过拟合。

### 5. 收敛性检查 (Sanity Check)

*   在将 ABACUS 计算出的数据喂给神经网络之前，建议编写一个简单的 Python 脚本扫描所有 `OUT.*/running_scf.log`。
*   **剔除标准**：如果某一步 SCF 未收敛（`converge` 标志为 false），或者能量出现非物理的巨大跳变（例如单步跳变 > 10 eV/atom），应将其从训练集中剔除，以免污染模型。