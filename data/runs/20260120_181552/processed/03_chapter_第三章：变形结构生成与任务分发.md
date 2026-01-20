# 第三章：变形结构生成与任务分发

在上一章中，我们完成了弹性常数计算的理论准备和环境配置。从本章开始，我们将进入核心实战环节。

计算弹性常数的本质，是观测晶体在受到微小形变（Strain）时，其能量或应力（Stress）的响应。为了获得完整的刚度矩阵，我们需要人为地构造一系列“变形晶胞”，算出它们的应力状态，最后通过拟合得到弹性常数。

本章将指导你如何生成这套标准的变形结构，并设计高效的批量计算工作流。

---

## 3.1 变形构型生成 (Pre-processing)

在 ABACUS 的弹性常数计算工作流中，我们通常采用 **“应力-应变法”（Stress-Strain Method）**。这要求我们对晶胞施加特定的拉伸或剪切应变，并计算对应的柯西应力（Cauchy Stress）。

### 3.1.1 核心逻辑：为什么是 24 种构型？

在自动化脚本（如 `gene_dfm.py`）的背后，遵循着一套严格的物理逻辑。为了求解完整的弹性张量，我们通常会生成 **24 个变形结构**。这个数字并非随机选择，而是基于以下考量：

1.  **6 种独立应变模式**：
    根据 Voigt 记号，三维晶体的应变张量有 6 个独立分量（$\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \varepsilon_{yz}, \varepsilon_{xz}, \varepsilon_{xy}$）。为了解出刚度矩阵的所有分量，我们需要分别激活这 6 种变形模式。

2.  **4 种应变幅值（Amplitudes）**：
    对于每种变形模式，我们通常取 4 个不同的应变点，例如：
    $$ \delta = \{ +0.01, -0.01, +0.005, -0.005 \} $$
    即 $\pm 1\%$ 和 $\pm 0.5\%$ 的形变。

    > **教授视点 (Professor's Insight)**：
    > **为什么要取 4 个点而不是 1 个点？**
    > 理论上，胡克定律（$ \sigma = C \varepsilon $）是线性的，一个点似乎就够了。但在 DFT 计算中，我们需要验证结果是否处于“线性响应区”。通过 4 个点进行线性拟合，可以有效消除数值噪声，并检测是否意外进入了非线性区（Anharmonicity）。

### 3.1.2 使用脚本生成结构

假设你已经有了一个经过**极高精度弛豫**的原始结构文件 `STRU_original`。我们将使用 Python 脚本（这里以通用的 `gene_dfm.py` 为例）来生成变形结构。

脚本的核心操作是将应变矩阵 $\boldsymbol{\varepsilon}$ 作用于原始晶格矢量 $\mathbf{R}$：
$$ \mathbf{R}' = \mathbf{R} (\mathbf{I} + \boldsymbol{\varepsilon}) $$

执行生成命令通常如下（示意）：
```bash
python gene_dfm.py --input STRU_original --max_strain 0.01 --points 4
```

执行后，你将获得如下目录结构：

```text
/work_dir/
├── STRU_original
├── gene_dfm.py
└── relax/                  <-- 生成的任务目录
    ├── 00_strain_mode1_p0.01/
    │   └── STRU            <-- 施加了 +1% xx应变的结构
    ├── 01_strain_mode1_m0.01/
    │   └── STRU            <-- 施加了 -1% xx应变的结构
    ├── ...
    └── 23_strain_mode6_m0.005/
        └── STRU
```

### 3.1.3 关键权衡：应变幅值的选择

在生成结构时，最大应变幅值（Maximum Strain Amplitude）的选择是一门艺术：

*   **太小（< 0.001）**：DFT 计算存在数值噪声（Numerical Noise）。如果应变引起的能量/应力变化小于计算误差，结果将完全不可信。
*   **太大（> 0.02）**：晶体将超出胡克定律的线性弹性区，高阶非线性项（三阶弹性常数）开始占主导，导致拟合出的二阶弹性常数不准确。
*   **推荐值**：对于大多数硬度正常的无机晶体，**0.005 (0.5%) 到 0.01 (1.0%)** 是最佳的“甜点区”。对于极软材料（如范德华晶体），可能需要适当增大；对于极硬材料，可适当减小。

---

## 3.2 批量任务管理与 INPUT 设置

结构生成后，我们需要为这 24 个任务准备 `INPUT` 文件并提交计算。

### 3.2.1 目录结构设计

为了便于后期数据处理，建议保持整洁的目录结构。一个标准的实战目录如下：

```text
elastic_calc/
├── 1_generation/
│   ├── STRU_original
│   └── gene_dfm.py
├── 2_running/          <-- 所有的计算在这里运行
│   ├── run_task.sh     <-- 批量提交脚本
│   ├── INPUT_template  <-- INPUT 模板文件
│   ├── KPT             <-- K点文件（所有任务共用）
│   ├── pseudo/         <-- 赝势目录
│   ├── task_00/
│   ├── task_01/
│   └── ...
└── 3_postprocessing/
```

### 3.2.2 核心配置文件：INPUT (Critical)

这是本章**最关键**的部分。在计算变形结构的应力时，`INPUT` 文件的设置有极高的排他性要求。

以下是一个标准的 `INPUT` 模板：

```bash
INPUT_PARAMETERS
# 1. 核心计算参数
suffix              elastic_task
calculation         relax       # <--- 极其重要！！见下文解释
cal_stress          1           # <--- 必须开启，否则无法计算应力张量

# 2. 精度控制 (必须与原始结构优化时保持一致或更高)
ecutwfc             80          # 推荐：根据赝势建议值设置
scf_thr             1.0e-8      # 电子自洽收敛精度
force_thr_ev        0.001       # 离子弛豫力的收敛精度

# 3. 涂抹与基组
smearing_method     gauss
smearing_sigma      0.05
basis_type          lcao        # 或 pw

# 4. 对称性设置 (风险控制)
symmetry            0           # <--- 建议关闭，见下文“风险提示”
```

#### 🚨 核心强调：calculation 参数
**绝对禁止**将 `calculation` 设置为 `cell-relax`。

*   **错误做法 (`cell-relax`)**：ABACUS 会优化晶胞矢量。这意味着程序会发现你施加的应变增加了系统能量，于是它会通过优化晶胞把应变“释放”掉，让晶格回到未变形的平衡态。结果：你计算的是平衡态的应力（接近0），得不到弹性常数。
*   **正确做法 (`relax`)**：固定晶胞形状（Fixed Cell），只优化原子位置（Relax Atoms）。这允许原子在变形后的晶胞内寻找新的平衡位置（Internal Relaxation），这是计算“内弹性常数”所必需的，同时保持了我们施加的外部应变。

### 3.2.3 批量提交脚本 (`run_task.sh`)

不要手动去 24 个文件夹里提交任务。编写一个简单的 Shell 脚本来分发任务。

```bash
#!/bin/bash

# 定义 ABACUS 可执行程序路径
ABACUS_EXEC="mpirun -np 32 abacus"
CURRENT_DIR=$(pwd)

# 循环遍历所有任务目录
for task_dir in task_*; do
    if [ -d "$task_dir" ]; then
        cd "$task_dir"
        echo "Processing $task_dir ..."

        # 1. 准备文件 (如果生成脚本未处理，需在此处复制)
        # cp ../INPUT_template INPUT
        # cp ../KPT KPT
        # ln -sf ../pseudo .

        # 2. 提交任务
        # 注意：实际生产环境中，这里通常是 sbatch submit.slurm
        # 下面是直接运行的示例：
        $ABACUS_EXEC > output.log 2>&1
        
        # 3. 简单检查是否正常结束
        if grep -q "Job is finished" output.log; then
            echo "  -> Done."
        else
            echo "  -> FAILED!"
        fi

        cd "$CURRENT_DIR"
    fi
done
```

---

## 3.3 专家锦囊与避坑指南

作为一名资深开发者，在此特别提醒几个容易导致全盘返工的细节：

### 1. 前置优化的精度陷阱
在执行 `gene_dfm.py` 之前，你的 `STRU_original` 必须处于**完美的零应力状态**。
*   **要求**：剩余应力（Residual Stress）应小于 0.1 kbar（甚至更低）。
*   **后果**：如果原始结构带有 2 kbar 的剩余应力，这个背景噪声会直接叠加到所有变形结构中，导致最终拟合的弹性常数出现严重偏差（尤其是 $C_{11}, C_{22}$ 等对角项）。

### 2. 对称性破缺的风险 (`symmetry 0`)
当对高对称性晶体（如立方晶系）施加剪切应变（如 $\varepsilon_{xy}$）时，晶体的对称性会显著降低（变成三斜晶系）。
*   **现象**：有时 ABACUS 的对称性分析模块在处理这些轻微变形的低对称结构时，可能会因为容差（tolerance）问题报错，或者强制施加了错误的对称性约束。
*   **对策**：为了安全起见，建议在计算变形结构时，在 `INPUT` 中显式设置 **`symmetry 0`**（关闭对称性分析）。虽然这会略微增加计算量（不再利用不可约 k 点），但能确保物理结果的绝对正确性。

### 3. 单位与后处理
ABACUS 的输出文件（如 `OUT.${suffix}/running_cell-relax.log` 或 `running_relax.log`）中，应力单位通常是 **Ry/Bohr³** 或 **kbar**（取决于版本和设置）。
*   **提示**：虽然本章只涉及计算，但请留意，最终的 Python 后处理脚本需要正确识别这些单位并转换为 **GPa**。通常脚本会自动处理，但如果你手动提取数据，请务必检查单位换算：
    $$ 1 \text{ Ry/Bohr}^3 \approx 14710.5 \text{ GPa} $$
    $$ 1 \text{ kbar} = 0.1 \text{ GPa} $$

---
**下一章预告**：
计算任务完成后，我们将进入第四章《数据提取与弹性张量拟合》，届时我们将学习如何从这 24 个复杂的输出文件中提取应力数据，并利用最小二乘法解出刚度矩阵。