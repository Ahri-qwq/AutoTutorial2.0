根据提供的参考资料和 ABACUS 与 DeepMD 的工作流逻辑，以下是关于 **DeepMD 训练数据生成流程** 的结构化元数据报告。

## 1. 物理本质 (Physics Concepts)

- **核心物理概念**:
    - **势能面采样 (PES Sampling)**: 通过第一性原理分子动力学 (AIMD) 或结构弛豫，在构型空间中采集具有代表性的原子位置、能量、受力和维里 (Virial) 数据。
    - **监督学习 (Supervised Learning)**: 将 DFT 计算得到的高精度数据（Label）作为“真值”，训练神经网络势函数（NNP），使其能以经典力学的成本预测量子力学精度的势能面。
- **该计算解决什么科学问题？**:
    - 解决传统 DFT 计算成本高、无法进行大尺度或长时间模拟的问题。
    - 解决经验力场精度不足、难以描述复杂化学反应或断键过程的问题。
    - 生成用于训练 Deep Potential (DP) 模型的高质量数据集。

## 2. 关键输入参数 (Key Parameters)

此部分主要关注**生成训练数据**所需的 ABACUS DFT 计算设置，以及数据处理工具的参数。

### 2.1 ABACUS INPUT (用于生成 DFT 原始数据)
为了生成可用于 DeepMD 训练的数据，ABACUS 必须运行在 DFT 模式（KS 轨道），并输出力和应力。

| 参数名 | 推荐值 | 物理意义 |
| :--- | :--- | :--- |
| `calculation` | `md` 或 `scf` | **计算类型**。通常使用 AIMD (`md`) 采样动态构型，或使用 `scf` 计算静态构型的单点能。 |
| `esolver_type` | `ks` (默认) | **能量求解器**。**注意**：生成训练数据必须使用 DFT (`ks`)，严禁设置为 `dp`（那是使用模型，见资料 2/5/7）。 |
| `cal_force` | `1` | **计算力**。训练 DP 模型必须包含原子受力信息。 |
| `cal_stress` | `1` | **计算应力**。训练 DP 模型通常需要维里 (Virial) 信息，特别是涉及晶胞变化的模拟。 |
| `md_dumpfreq` | `1` | **MD 输出频率**。决定了采样密度，通常每步都输出以最大化数据利用率。 |
| `out_stru` | `1` (需确认) | **输出结构**。确保每步的结构信息被记录（资料未明确提及参数名，但逻辑上必须有结构文件输出，通常包含在 MD 输出流中）。 |

### 2.2 Python / dpdata (数据转换)
使用 `dpdata` 将 ABACUS 数据转换为 DeepMD 格式。

| 参数/方法 | 关键设置 | 说明 |
| :--- | :--- | :--- |
| `dpdata.LabeledSystem` | `fmt='abacus/md'` | **读取格式**。指定读取 ABACUS 的 MD 轨迹格式（资料 1）。 |
| `to_deepmd_npy` | 路径字符串 | **输出方法**。将数据导出为 DeepMD 训练所需的 NumPy 压缩格式。 |

### 2.3 DP-GEN param.json (如果使用主动学习)
如果使用 DP-GEN 自动化流程（资料 3）。

| 参数名 | 设置 | 说明 |
| :--- | :--- | :--- |
| `init_fp_style` | `"ABACUS"` | **第一性原理代码接口**。指定使用 ABACUS 进行 Labeling 步骤。 |
| `stages` | `[1,2,3,4]` | **迭代阶段**。定义训练-探索-标注的循环。 |

## 3. 体系与接口配置 (System & Interfaces)

- **结构 (STRU)**:
    - 需要标准的 ABACUS `STRU` 文件。
    - 在 DP-GEN 流程中，需要为每个构型建立单独文件夹（如 `init/3C`, `init/2H`）。
- **外部接口**:
    - **dpdata**: 必选。用于数据清洗、格式转换（ABACUS -> DeepMD）。
    - **DeepMD-kit**: 用于后续的模型训练（输入为 dpdata 的输出）。
    - **DP-GEN**: 可选。用于自动化并发管理 "训练 -> 探索 (DPMD) -> 标注 (ABACUS DFT)" 的闭环。
- **接口注意事项**:
    - **格式对应**: `dpdata` 读取 ABACUS MD 数据时，需指向包含 `MD_dump` 或相关输出的目录（资料 1 示例为 `./LiCl_DP_Tutorial_Example/chapter2/abacus_md`）。
    - **文件依赖**: 若使用 DP-GEN，除了 `STRU`，还需要准备 `*.orb` (轨道), `*.upf` (赝势), `KPT`, `machine.json` 和 `param.json`。

## 4. 教程编写特殊指令 (Special Instructions for Writer)

- **Critical (关键区分)**:
    - **极易混淆点**: 资料中同时出现了“生成训练数据”（资料 1, 3）和“运行 DPMD”（资料 2, 5, 7）。
    - **撰写要求**: 必须明确区分 **DFT 采样** 和 **DPMD 推理**。
        - **生成数据时**: `esolver_type` 必须是 DFT 相关（默认 `ks`），**绝对不能**是 `dp`。
        - **使用模型时**: `esolver_type` 才是 `dp`。
    - 本 Topic 聚焦于“训练数据生成”，因此重点应放在如何跑通 DFT MD 以及如何设置 `cal_force=1`, `cal_stress=1`，随后如何用 `dpdata` 转换。
- **数据划分策略**:
    - 引用资料 1 和 6，说明通常需要将数据集划分为训练集 (Training) 和验证集 (Validation)。例如：随机抽取 100 帧作为验证集，其余作为训练集。
- **风险提示**:
    - 资料中未明确提及 ABACUS MD 输出的具体文件名（如 `MD_dump` 或 `Running_MD` 等）在 `dpdata` 中的具体对应关系，建议提醒读者查阅 `dpdata` 官方文档中关于 `abacus/md` 格式的详细说明，确保目录结构正确。

## 5. 常见报错与注意事项 (Pitfalls)

- **数据格式不匹配**: `dpdata` 对 ABACUS 输出文件的版本或命名敏感。如果 ABACUS 版本更新导致输出格式微调，`dpdata` 可能报错。需确保两者版本兼容。
- **缺少关键数据**: 如果 DFT 计算时忘记设置 `cal_stress 1`，生成的训练数据将不包含维里信息，导致训练出的模型无法正确预测压强或晶格常数。
- **KPT 设置**: DFT 计算（生成数据）必须设置 `KPT`，而 DPMD（使用模型）不需要 `KPT`。初学者容易在切换两种模式时混淆。
- **收敛性**: 在高温 AIMD 采样时，需注意 SCF 的收敛性，未收敛的步数产生的 Force/Energy 是噪声，会污染训练集，应在 `dpdata` 处理时剔除或在计算时优化参数。