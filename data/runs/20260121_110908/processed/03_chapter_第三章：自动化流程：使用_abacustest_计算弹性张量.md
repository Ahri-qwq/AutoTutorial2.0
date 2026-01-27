# 第三章：自动化流程：使用 abacustest 计算弹性张量

这里是《ABACUS 实战教程》第三章的完整内容。本章严格遵循您的专家级要求，以 `abacustest` 自动化工作流为核心，明确区分了结构优化与变形计算两个关键阶段，并对物理原理进行了深度解析。

---

# 第三章：自动化流程：使用 abacustest 计算弹性张量

在材料科学中，弹性张量（Elastic Tensor, $C_{ij}$）是描述材料力学性质的指纹。它不仅决定了材料的硬度（体模量、剪切模量），还是判断晶体结构力学稳定性的关键判据（波恩稳定性准则）。

计算弹性张量的传统方法是“应力-应变法”（Stress-Strain Method）：手动对晶胞施加一系列微小形变，计算对应的应力，最后通过胡克定律拟合。这一过程涉及数十个算例的生成、管理和数据提取，手工操作极易出错。

本章将介绍如何使用 ABACUS 官方推荐的辅助工具 **`abacustest`**，一键完成从变形结构生成到结果拟合的全自动化工作流。

## 3.1 流程概览与环境配置

### 3.1.1 核心工作流
计算弹性常数必须严格遵循以下两个阶段，切勿混淆：

1.  **阶段一：原始结构优化 (Zero-strain Optimization)**
    *   **目标**：消除晶胞内部的残余应力。
    *   **方法**：使用 `calculation cell-relax`（变胞优化）。
    *   **重要性**：如果初始结构存在残余应力（例如 > 1 kBar），后续施加微小应变（如 1%）引起的应力变化将被由于未弛豫导致的误差掩盖，导致结果完全错误。

2.  **阶段二：变形结构计算 (Deformation & Stress Calculation)**
    *   **目标**：计算施加应变后的应力响应。
    *   **方法**：使用 `abacustest` 自动生成变形结构，进行 `calculation relax`（定胞动原子）或 `calculation scf`（定胞定原子）计算。

### 3.1.2 环境准备
`abacustest` 的弹性计算模块依赖于 `pymatgen` 进行晶体操作。请确保你的 Python 环境中已安装以下库：

```bash
pip install abacustest pymatgen ase
```

---

## 3.2 步骤一：生成变形结构 (Deformation Generation)

假设你已经完成了**阶段一**，获得了一个充分弛豫的结构文件（例如 `STRU_relaxed`）。现在我们将进入**阶段二**。

我们将使用 `abacustest` 的库功能来生成一系列施加了正应变（Normal strain）和剪切应变（Shear strain）的子任务文件夹。

### 3.2.1 命令行操作
在包含 `STRU_relaxed` 和势函数文件的目录下，运行以下命令：

```bash
# 假设你的优化后结构文件名为 STRU
# 准备生成弹性计算任务，指定任务类型为 elastic
abacustest lib -t elastic \
    --stru STRU \
    --orb <你的轨道文件>.orb \
    --pp <你的赝势文件>.upf \
    --norm 0.01 \
    --shear 0.01 \
    --kpt <KPT文件> 
```

### 3.2.2 关键参数详解
*   **`--norm 0.01`**: 设置最大正应变幅度为 1% (0.01)。工具通常会生成 $\pm 0.005$ 和 $\pm 0.01$ 四个点。
*   **`--shear 0.01`**: 设置最大剪切应变幅度为 1%。
*   **`--norelax` (高风险参数)**: 
    *   **默认情况（不加此参数）**: 生成的 INPUT 文件中 `calculation` 为 `relax`。这意味着在变形后的晶胞中，允许**原子内部坐标弛豫**（Relaxed Ions）。这对应于**绝热弹性常数**，通常与实验测量的低频弹性模量一致（推荐）。
    *   **加上 `--norelax`**: 生成的 INPUT 文件中 `calculation` 为 `scf`。这意味着原子坐标被强制固定在仿射变换后的位置（Clamped Ions）。这对应于**高频极限**下的弹性常数。

### 3.2.3 生成的目录结构
命令执行后，当前目录下会出现名为 `elastic_result`（或类似命名）的文件夹，内部包含如下结构：
*   `run_00/`: 原始结构的计算任务（用于校验残余应力）。
*   `run_01/` ~ `run_24/`: 对应不同应变模式（Strain Patterns）的子任务。
    *   每个子文件夹内都包含独立的 `INPUT`, `STRU`, `KPT` 及势函数链接。

---

## 3.3 步骤二：批量提交与应力计算

在提交计算之前，作为高阶用户，你需要检查生成的 `INPUT` 文件，确保物理参数设置正确。

### 3.3.1 检查 INPUT 核心参数
进入任意一个子文件夹（如 `run_01`），检查 `INPUT` 文件。以下参数对于弹性计算至关重要：

```bash
# INPUT file snippet

# 1. 计算模式 (由是否使用 --norelax 决定)
calculation     relax   # 推荐: 允许原子弛豫 (Relaxed ions)
# calculation   scf     # 仅当使用了 --norelax 时 (Clamped ions)

# 2. 应力计算 (必须开启!)
cal_stress      1       # 核心参数：计算应力张量

# 3. 对称性设置
symmetry        0       # 强烈推荐：关闭对称性
                        # 原因：施加形变后晶体对称性降低。
                        # 强制开启对称性可能导致 ABACUS 错误地将应力张量对称化，
                        # 或者限制原子的弛豫方向，导致结果失真。

# 4. 精度控制
ecutwfc         100     # 需与阶段一保持一致或更高
scf_thr         1e-8    # 严格的收敛标准有助于获得精确的应力值
force_thr_ev    0.001   # 如果是 relax，力收敛标准要足够严
```

### 3.3.2 批量提交任务
使用 `abacustest` 的提交功能或你自己的调度脚本批量运行这些任务。

```bash
# 使用 abacustest 提交 (假设你配置好了 machine.json)
abacustest submit -p elastic_result/ -c machine.json
```

或者，如果你习惯使用 Shell 脚本：
```bash
for dir in elastic_result/run_*; do
    cd $dir
    mpirun -np 16 abacus > abacus.log &  # 示例命令
    cd ../..
done
```

---

## 3.4 步骤三：数据后处理与拟合

当所有子任务计算完成后（即所有目录下都生成了 `OUT.` 文件夹），我们使用 `abacustest` 进行数据收集和拟合。

### 3.4.1 执行后处理
在包含 `elastic_result` 的父目录下运行：

```bash
# 收集数据并拟合
abacustest post -t elastic -d elastic_result/
```

### 3.4.2 结果分析与验证
程序将输出 Voigt 标记法下的 $6 \times 6$ 弹性刚度矩阵 $C_{ij}$ 以及模量性质。

**示例输出 (Si 晶体)**:
```text
Elastic Tensor Cij (GPa):
      155.5    58.2    58.2     0.0     0.0     0.0
       58.2   155.5    58.2     0.0     0.0     0.0
       58.2    58.2   155.5     0.0     0.0     0.0
        0.0     0.0     0.0    76.2     0.0     0.0
        0.0     0.0     0.0     0.0    76.2     0.0
        0.0     0.0     0.0     0.0     0.0    76.2

Bulk Modulus (Voigt): 90.6 GPa
Shear Modulus (Voigt): 68.1 GPa
```

**专家验证清单 (Checklist)**：
1.  **对称性检查**：
    *   对于立方晶系（如 Si, Cu），必须满足 $C_{11}=C_{22}=C_{33}$，$C_{12}=C_{13}=C_{23}$，$C_{44}=C_{55}=C_{66}$。
    *   非对角项（如 $C_{14}$）应接近于 0（通常 < 1-2 GPa 的波动是允许的数值误差）。
2.  **正定性检查**：
    *   确保特征值均为正值，否则结构力学不稳定。
3.  **结果对比**：
    *   将结果与 Materials Project 或参考文献对比。
    *   **风险提示**：Materials Project 数据库通常会将晶体旋转到标准取向（Standard Orientation）。如果你的输入结构是旋转过的（例如 $z$ 轴沿 [111] 方向），你计算出的 $C_{ij}$ 矩阵将是旋转后的结果，与标准数据库直接对比会由于坐标系不同而看似“错误”。请务必确认你的 `STRU` 晶格矢量方向。

---

## 附录：Workflow B (Python 脚本模式)
> **Note**: 这是早期的手动脚本方法，仅供参考，现代工作流推荐使用上述的 `abacustest`。

在 `abacustest` 集成之前，用户通常使用 `gene_dfm.py`（基于 pymatgen 生成变形）和 `post_process.py` 进行计算。这要求用户手动修改 Python 脚本中的应变参数，并手动编写 Shell 脚本来管理文件夹。虽然原理与 Workflow A 相同，但维护成本较高且容易出错。