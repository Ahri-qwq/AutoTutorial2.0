# 第一章：计算准备与物理基础

在探索材料的电子结构之前，我们必须建立一个坚实的物理基石。能带结构（Band Structure）本质上是电子在周期性势场中的能量色散关系，而这个势场由原子核位置和基态电子密度共同决定。

本章将指导你完成能带计算最关键的前置步骤：构建合理的晶体结构文件，并进行高精度的基态自洽（SCF）计算。

> **专家提示**：能带计算遵循 "Garbage In, Garbage Out" 原则。如果基态电荷密度（Charge Density）未收敛或结构不合理，后续无论是传统的 NSCF 计算还是基于紧束缚模型的 PYATB 分析，其结果都将毫无物理意义。

---

## 1.1 晶胞选择与结构文件 (STRU)

### 1.1.1 原胞 (Primitive Cell) vs. 超胞 (Supercell)
在进行能带计算时，晶胞的选择直接决定了布里渊区（Brillouin Zone）的形状和能带的折叠方式。

*   **标准能带计算**：请务必使用**原胞 (Primitive Cell)**。这是因为能带图通常是沿着原胞布里渊区的高对称线（High Symmetry Lines）绘制的。如果使用超胞，能带会发生折叠（Band Folding），导致本征值数量增加且难以辨认物理图像，除非你明确需要研究缺陷态或进行能带反折叠（Unfolding）分析。
*   **PYATB/紧束缚路线**：虽然 PYATB 支持反折叠功能，但在基础学习阶段，依然建议从原胞开始，以确保物理图像的直观性。

### 1.1.2 详解 STRU 文件
`STRU` 文件定义了晶格常数、晶格矢量、原子种类及位置。它是 ABACUS 识别物理系统的核心。

以下是一个典型的 `STRU` 文件模板（以硅为例）：

```bash
ATOMIC_SPECIES
Si  28.0855  Si_ONCV_PBE-1.0.upf  Si_Mz1.3-5sp.orb
# 格式说明：
# 元素符号  原子质量  赝势文件名  (LCAO模式下必需)数值原子轨道文件名

LATTICE_CONSTANT
10.2631  # 晶格常数缩放因子，单位通常为 Bohr (1 Bohr ≈ 0.529 Angstrom)

LATTICE_VECTORS
0.5 0.5 0.0  # 矢量 a1
0.5 0.0 0.5  # 矢量 a2
0.0 0.5 0.5  # 矢量 a3

ATOMIC_POSITIONS
Direct  # 坐标类型：Direct (分数坐标) 或 Cartesian (笛卡尔坐标，单位为晶格常数)
Si      # 元素名称，需与 ATOMIC_SPECIES 对应
0.0     # 初始磁矩 (m_start)，非磁性计算设为 0.0
2       # 该元素原子数量
0.00 0.00 0.00 1 1 1  # 原子1坐标 及 移动约束(1=允许移动, 0=固定)
0.25 0.25 0.25 1 1 1  # 原子2坐标
```

### 1.1.3 LCAO 模式下的特殊依赖
如果你计划使用 **PYATB** 进行后处理，或者在 ABACUS 中使用 `basis_type lcao`（原子轨道线性组合）进行计算，`ATOMIC_SPECIES` 部分必须严格包含**数值原子轨道文件（.orb）**。

*   **匹配性原则**：`.orb` 文件必须与赝势文件（`.upf` / `.vps`）在交换关联泛函（如 PBE）和截断半径上保持物理上的一致性。
*   **文件路径**：ABACUS 会在 `INPUT` 文件中指定的 `pseudo_dir` 和 `orbital_dir` 目录下寻找这些文件。

---

## 1.2 基态自洽计算 (SCF)

所有能带计算的第一步都是通过自洽场（Self-Consistent Field, SCF）迭代求解 Kohn-Sham 方程，以获得收敛的基态电子密度。

### 1.2.1 区分两条计算路线
在开始 SCF 之前，你必须明确你的后续目标，因为这决定了 `INPUT` 文件的设置策略：

1.  **路线 A：传统 NSCF 路线**
    *   **目标**：获得高精度的电荷密度文件（如 `SPIN1_CHG.cube` 或密度目录）。
    *   **后续**：读取该密度，进行非自洽（NSCF）计算，K 点取沿高对称线的 Line 模式。
    *   **适用场景**：标准 DFT 能带计算，使用平面波（PW）或 LCAO 基组均可。

2.  **路线 B：PYATB 紧束缚路线**
    *   **目标**：获得哈密顿量（Hamiltonian）和重叠矩阵（Overlap Matrix）等紧束缚矩阵信息。
    *   **后续**：运行 PYATB 程序处理矩阵，得到能带、Fat band 或进行反折叠。
    *   **适用场景**：需要高效处理大规模体系、投影能带或类似 Wannier 函数的分析，**必须使用 LCAO 基组**。

### 1.2.2 准备 INPUT 文件 (SCF 阶段)
为了工程管理的规范性，建议将 SCF 的输入文件命名为 `INPUT_scf`，计算时再复制为 `INPUT`。

**通用核心参数 (INPUT_scf)**：

```bash
INPUT_PARAMETERS
# --- 基础设置 ---
suffix          Si_scf      # 输出文件后缀，用于区分任务
calculation     scf         # 计算类型：自洽计算
basis_type      lcao        # 基组类型：lcao (PYATB必需) 或 pw (传统路线常用)
ntype           1           # 元素种类数量

# --- 精度控制 ---
ecutwfc         50          # 平面波截断能 (Ry)，PW 模式下需严格收敛测试
scf_thr         1.0e-7      # 自洽收敛阈值。能带计算建议设为 1e-7 或更高精度

# --- 路径设置 ---
pseudo_dir      ./PP_ORB/   # 赝势文件夹路径
orbital_dir     ./PP_ORB/   # 轨道文件夹路径 (LCAO模式需要)

# --- 结构优化 (可选，但推荐先做) ---
# relax_nmax    100         # 如果需要先优化结构，开启此项并设 calculation 为 relax
```

> **风险提示 (针对 PYATB 路线)**：
> 如果你的目标是路线 B (PYATB)，单纯的 SCF 计算是不够的。你需要让 ABACUS 输出紧束缚矩阵（HR, SR 等）。
> **请务必查阅 ABACUS 官方文档中关于 "Tight-Binding Interface" 或 "Matrices Output" 的章节**，在 `INPUT` 文件中添加开启矩阵输出的特定参数（例如涉及 `out_mat_...` 类的参数）。由于版本更新，具体参数名请以官方文档为准。

### 1.2.3 准备 KPT 文件 (SCF 阶段)
**特别注意**：在 SCF 阶段，必须使用均匀网格（Monkhorst-Pack）覆盖整个布里渊区，以准确积分电荷密度。**切勿使用 Line 模式！**

`KPT` 文件示例（SCF 阶段）：
```bash
K_POINTS
0
Gamma           # 生成模式：Gamma Centered 或 Monkhorst-Pack
8 8 8 0 0 0     # 8x8x8 的均匀网格，后三位通常为 0 0 0 (Shift)
```

### 1.2.4 执行计算与工作流管理
推荐使用脚本化的工作流来避免手动修改文件的错误。

**最佳实践脚本逻辑**：
```bash
# 1. 准备阶段
cp INPUT_scf INPUT    # 将准备好的 SCF 输入文件复制为标准输入名
cp KPT_scf KPT        # 确保使用均匀网格 K 点

# 2. 运行 ABACUS (根据实际环境调用 mpirun 或 yhrun)
mpirun -np 16 abacus

# 3. 检查收敛
# 检查输出日志中是否出现 "Convergence achieved"
```

### 1.2.5 关键产物检查
计算完成后，根据你的路线检查输出：
*   **传统路线**：检查 `OUT.suffix/` 目录下是否生成了电荷密度文件（通常是 `SPIN1_CHG.cube` 或存储密度的二进制文件夹）。这是下一步 NSCF 计算的“燃料”。
*   **PYATB 路线**：检查是否生成了矩阵文件（通常位于输出目录下的矩阵子文件夹中）。

完成这一步后，我们就拥有了准确的基态电子密度，可以根据第二章的指引，分流进入“传统 NSCF 能带计算”或“PYATB 紧束缚能带计算”。