# ABACUS 电子能带结构计算与进阶分析实战指南

## 前言

### 欢迎开启 ABACUS 探索之旅
欢迎阅读本教程。ABACUS（Atomic-orbital Based Ab-initio Computation Software）作为一款优秀的国产开源电子结构计算软件，在处理大规模体系及高效计算方面展现出显著优势，尤其是其基于数值原子轨道的线性组合（LCAO）方案，在计算效率与精度之间取得了极佳的平衡。虽然对于初学者而言，其输入参数的精细调节具有一定的学习曲线，但一旦掌握，它将成为你探索微观物质世界的利器。

### 学习路线图 (Roadmap)
为了帮助读者系统性地掌握能带计算，本书精心设计了由浅入深的三阶学习路径：
1.  **第一章：物理基石**。重点讲解结构构建与自洽场（SCF）计算。这是所有电子结构分析的起点，旨在确保你的计算输入“言之有物”。
2.  **第二章：标准路径**。深入探讨非自洽场（NSCF）方法。这是最符合物理直觉的计算方式，用于获取布里渊区高对称路径上的本征能量。
3.  **第三章：进阶分析**。引入 PYATB 与紧束缚模型。针对复杂体系（如掺杂、界面），讲解如何进行轨道成分（Fat Band）分析及超胞能带反折叠，解决传统方法难以应对的痛点。

### 知识体系定位 (Knowledge Graph)
在整个 ABACUS/DFT 知识体系中，本教程处于核心地位。它承接了“结构优化”这一前置环节，并为后续的态密度（DOS）分析、有效质量提取、光学性质计算以及输运性质研究提供了关键的电子结构信息。它是从“几何结构”跨越到“物理性能”的必经桥梁。

### 前置要求 (Prerequisites)
在开始本教程之前，我们建议你已具备以下基础：
- **物理基础**：了解固体物理基础概念（如倒空间、布里渊区、费米面）。
- **理论背景**：对密度泛函理论（DFT）及 Kohn-Sham 方程有基本认知。
- **操作技能**：熟悉 Linux 命令行基础操作，并已成功安装 ABACUS 及其相关依赖工具（如 PYATB）。

---

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

# 第二章：标准能带计算路线 (NSCF Method)

在上一章中，我们已经掌握了如何通过自洽场（SCF）计算获得体系的基态电子密度。本章将进入电子结构分析的核心环节——**能带结构（Band Structure）计算**。

在 ABACUS 中，计算能带主要有两种路线：
1.  **标准 NSCF 路线（本章重点）**：基于 SCF 收敛的电荷密度，固定哈密顿量，沿着高对称 K 点路径进行非自洽（Non-Self-Consistent Field）计算，直接求解本征值。这是最物理、最标准的做法。
2.  **Tight-Binding/Wannier 路线（后续章节）**：基于 SCF 输出的哈密顿量矩阵（HR/SR），使用后处理工具（如 PYATB）进行紧束缚模型插值。

本章我们将聚焦于**第一种路线**，即“分步计算法”：SCF $\rightarrow$ NSCF。

---

## 2.1 高对称 K 点路径设置

能带图本质上是电子能量 $E$ 随波矢量 $k$ 变化的色散关系。为了直观展示能带，我们需要在布里渊区（Brillouin Zone）中选取具有代表性的高对称点（如 $\Gamma, X, M, K$ 等）并连成路径。

### 2.1.1 理解 KPT 文件模式
在 SCF 计算中，我们通常使用 `Gamma` 或 `mp`（Monkhorst-Pack）模式来均匀采样布里渊区。而在能带计算（NSCF）中，必须将 `KPT` 文件切换为 **Line** 模式。

### 2.1.2 使用辅助工具生成路径
虽然可以手动查表，但强烈建议使用 `Atomkit`、`SeeK-path` 或 `VASPKIT` 等工具自动生成标准路径，以确保不遗漏关键的高对称点。

假设我们通过工具获得了硅（Si）的 FCC 布里渊区路径：$L \rightarrow \Gamma \rightarrow X \rightarrow U \rightarrow K \rightarrow \Gamma$。

### 2.1.3 标准 Line 模式 KPT 示例
一个典型的用于能带计算的 `KPT` 文件内容如下：

```fortran
K_POINTS
20  # 每个线段（Segment）上的插值点数
Line
6   # 高对称点的总数
0.500 0.500 0.500 1.0  # L
0.000 0.000 0.000 1.0  # Gamma
0.500 0.000 0.500 1.0  # X
0.625 0.250 0.625 1.0  # U
0.375 0.375 0.750 1.0  # K
0.000 0.000 0.000 1.0  # Gamma
```

**参数详解：**
*   **第二行 (20)**:  表示每两个高对称点之间插入 20 个 K 点。数值越大，能带图越平滑，但计算量越大。
*   **第三行 (Line)**: **关键模式**。告诉 ABACUS 这是一个路径计算，而不是均匀网格。
*   **第四行 (6)**:  接下来定义的高对称点数量。

---

## 2.2 非自洽计算 (NSCF) 参数配置

NSCF 计算的核心逻辑是：**读取**上一步 SCF 得到的电荷密度（`SPIN1_CHG.cube` 或密度矩阵），构建势函数 $V_{eff}$，然后**不再更新**电荷密度，仅对新的 K 点（路径上的点）进行本征值求解。

我们需要修改 `INPUT` 文件来实现这一逻辑。

### 2.2.1 关键参数清单

以下是 `INPUT` 文件中必须调整的核心参数：

| 参数名 | 推荐值/设定值 | 物理含义与解释 |
| :--- | :--- | :--- |
| **`calculation`** | `nscf` | **核心开关**。指定计算类型为非自洽计算。 |
| **`init_chg`** | `file` | **核心开关**。指示 ABACUS 从文件中读取电荷密度，而不是从原子叠加开始初始化。这是连接 SCF 和 NSCF 的桥梁。 |
| **`symmetry`** | `0` | **必须关闭**。在 Line 模式下，K 点路径通常不具备完整的晶体对称性。开启对称性可能导致 K 点被错误折叠或约化，导致能带路径断裂。 |
| `out_band` | `1` | 输出能带本征值文件（通常在输出目录的 `BANDS_1.dat` 等文件中）。 |
| `out_bandgap` | `1` | 在日志文件（`running_nscf.log`）中自动分析并输出带隙信息（VBM, CBM, Gap）。 |
| `scf_nmax` | `1` (或小值) | 对于 LCAO 基组，NSCF 本质上是一次对角化过程，不需要像 SCF 那样迭代混合电荷密度。 |

### 2.2.2 INPUT 文件示例 (NSCF)

```bash
INPUT_PARAMETERS
# 1. General
suffix          Si_band      # 输出目录后缀
ntype           1
symmetry        0            # [CRITICAL] 必须关闭对称性

# 2. Calculation Mode
calculation     nscf         # [CRITICAL] 非自洽模式
init_chg        file         # [CRITICAL] 读取电荷密度文件
dft_functional  pbe

# 3. Basis & Convergence
basis_type      lcao
ecutwfc         60           # 必须与 SCF 保持一致
scf_thr         1e-7         # 对角化收敛精度

# 4. Output
out_band        1            # 输出能带数据
out_bandgap     1            # 输出带隙分析
```

---

## 2.3 数据流与作业提交规范

为了保证计算的严谨性，建议采用**分步目录管理**或**分步脚本控制**的方式。切勿在同一个目录下频繁覆盖文件，导致数据混乱。

### 2.3.1 标准操作流程 (Workflow)

**第一步：自洽计算 (SCF)**
1.  准备 `INPUT` (设 `calculation scf`, `symmetry 1`)。
2.  准备 `KPT` (使用 Gamma 或 mp 模式)。
3.  运行 ABACUS。
4.  **检查收敛**：确保 SCF 达到 `scf_thr` 收敛标准。
5.  **确认产物**：检查输出目录（默认为 `OUT.suffix`）下是否生成了电荷密度文件。
    *   对于 LCAO/PW 体系，通常是 `SPIN1_CHG.cube`。

**第二步：非自洽计算 (NSCF)**
1.  **准备环境**：
    *   保留 SCF 的输出目录。
    *   修改 `INPUT` 为 NSCF 模式（见 2.2 节）。
    *   替换 `KPT` 为 Line 模式（见 2.1 节）。
2.  **链接电荷密度**：
    *   ABACUS 默认会在 `OUT.suffix` 目录下寻找电荷密度文件。
    *   **注意**：如果你更改了 NSCF 的 `suffix`（例如从 `Si_scf` 改为 `Si_band`），ABACUS 会去 `OUT.Si_band` 下找电荷密度，但这显然是空的。
    *   **解决方案 A (推荐)**：保持 `suffix` 不变，或者手动将 SCF 的 `OUT.Si_scf` 目录下的电荷文件复制到当前工作区或指定目录。
    *   **解决方案 B**：使用 `read_file_dir` 参数显式指定 SCF 结果所在的路径。

### 2.3.2 实战脚本示例

以下是一个推荐的 Shell 脚本逻辑，展示了如何在一个脚本中串联两步计算：

```bash
#!/bin/bash

# --- 阶段 1: SCF 计算 ---
echo "Start SCF Calculation..."
# 准备 SCF 的输入文件
cp INPUT_scf INPUT
cp KPT_scf KPT

# 运行 ABACUS (假设可执行文件名为 abacus)
mpirun -np 16 abacus > log_scf.txt

# 检查是否生成电荷密度 (假设 suffix 为 ABACUS)
if [ ! -f "OUT.ABACUS/SPIN1_CHG.cube" ]; then
    echo "Error: SCF failed or charge density not found!"
    exit 1
fi

# --- 阶段 2: NSCF 计算 ---
echo "Start NSCF Calculation..."
# 准备 NSCF 的输入文件
cp INPUT_nscf INPUT
cp KPT_line KPT

# 关键：确保 NSCF 能读到 SCF 的电荷
# 方法：INPUT_nscf 中设置 init_chg file，且 suffix 保持一致
# 或者显式指定读取路径 (如果支持)

mpirun -np 16 abacus > log_nscf.txt

echo "All Done. Check OUT.ABACUS/BANDS_1.dat"
```

### 2.3.3 结果验证

计算完成后，进入输出目录（如 `OUT.ABACUS`），你应该关注以下文件：
*   **`running_nscf.log`**: 
    *   搜索 "Fermi Energy" 确认费米能级。
    *   搜索 "Band Gap" 查看程序自动计算的带隙。
*   **`BANDS_1.dat`**: 
    *   这是原始能带数据文件，包含 K 点坐标和对应的本征值。
    *   可以使用 `abacus-plot` (Python工具) 或手动编写脚本将其可视化。

> **专家提示 (Expert Tip)**: 
> 如果在 `running_nscf.log` 中看到 "Charge density file not found" 错误，请立即检查：
> 1. 上一步 SCF 是否真的成功结束？
> 2. `INPUT` 中的 `suffix` 是否与上一步一致？
> 3. 是否误删了 `OUT.suffix` 目录？

---

**下一章预告**：
掌握了标准的 NSCF 流程后，我们将学习如何处理更复杂的体系。第三章将介绍如何利用 **PYATB** 结合 ABACUS 的 Tight-Binding 接口进行能带反折叠（Band Unfolding）和投影能带（Fat Band）分析。

# 第三章：进阶能带分析路线 (PYATB & Tight-Binding)

在上一章中，我们通过传统的 **NSCF（非自洽场）** 路线，利用 `calculation = nscf` 和 Line 模式的 K 点，成功绘制了基础能带图。这种方法对于完美的晶体结构非常有效，能直接给出本征值。

然而，当我们面对更复杂的物理场景时，传统 NSCF 路线会显得力不从心：
1.  **轨道成分分析**：如果你想知道某条能带是由 $d_{xy}$ 轨道贡献还是 $p_z$ 轨道贡献（即 Fat Band 投影能带）。
2.  **大体系反折叠**：当你计算掺杂、缺陷或界面体系时，必须使用超胞（Supercell）。超胞会导致布里渊区折叠，产生极其密集且难以辨认的“意大利面”能带。此时需要 **能带反折叠 (Band Unfolding)** 技术将其还原回原胞的布里渊区。

为了解决这些问题，我们需要引入 **PYATB (Python Ab-initio Tight-Binding)** 路线。

> **核心区分 (Critical Distinction)**：
> *   **传统 NSCF 路线**：`SCF` -> `NSCF` (读取电荷密度) -> 获得本征值。
> *   **PYATB 紧束缚路线**：`SCF` (开启矩阵输出) -> `PYATB` (后处理) -> 获得能带/反折叠/投影。
>
> **注意**：本章的计算流程与第二章完全独立，请勿混淆。

---

## 3.1 构建紧束缚哈密顿量 (ABACUS 端)

PYATB 的工作原理是基于 ABACUS 计算出的**局域轨道（LCAO）哈密顿量矩阵**。因此，第一步是在 ABACUS 中进行一次 SCF 计算，并命令其将这些矩阵写出到磁盘。

### 3.1.1 准备输入文件

我们需要进行一次标准的自洽计算（`calculation = scf`），但在 `INPUT` 文件中需要增加特定的参数以输出矩阵。

**关键 `INPUT` 参数设置**：

```bash
INPUT_PARAMETERS
# ... 基础参数 (ecutwfc, scf_nmax 等) ...

# 1. 必须使用 LCAO 基组
basis_type      lcao

# 2. 计算模式为自洽计算
calculation     scf

# 3. 开启紧束缚矩阵输出 (关键步骤)
# 注意：不同版本的 ABACUS 参数名可能略有不同，请以官方文档 "Tight-Binding Matrices" 章节为准
out_mat_hs2     1        # 输出稀疏格式的哈密顿量(H)和重叠矩阵(S)
out_mat_r       1        # 输出位置算符偶极矩阵(rR)，用于计算 Berry Phase 或光学性质

# 4. 推荐设置
symmetry        1        # 生成矩阵时通常保持对称性，PYATB 会处理
cal_force       0        # 能带分析通常不需要算力
cal_stress      0        # 不需要算应力
```

### 3.1.2 运行计算与输出检查

运行 ABACUS 可执行文件：
```bash
mpirun -np 16 abacus | tee out.log
```

计算结束后，检查输出目录（默认为 `OUT.ABACUS/`），你必须找到以下关键文件，它们是连接 ABACUS 与 PYATB 的桥梁：

1.  **`data-HR-sparse_SPIN0.csr`** (或类似命名): 哈密顿量矩阵 $H(R)$。
2.  **`data-SR-sparse_SPIN0.csr`**: 重叠矩阵 $S(R)$。
3.  **`data-rR-sparse_SPIN0.csr`**: 偶极矩阵 $\vec{r}(R)$。
4.  **`STRU`**: 结构文件（PYATB 需要读取晶格信息）。
5.  **`*.orb`**: 原子轨道文件（**极度重要**，PYATB 需要读取轨道信息来构建基组，请确保该文件在当前目录或路径正确）。

> **风险提示**：如果未生成 `.csr` 文件，请检查 `out_mat_hs2` 参数是否被正确识别，或查阅你所使用的 ABACUS 版本的官方文档中关于 "Tight Binding Interface" 的说明。

---

## 3.2 投影能带 (Fat Band) 计算

得到矩阵后，我们脱离 ABACUS 主程序，转而使用 Python 脚本调用 PYATB 库进行后处理。Fat Band 可以直观地展示不同原子轨道对能带的贡献权重。

### 3.2.1 准备 PYATB 脚本

创建一个 Python 脚本（例如 `run_fatband.py`）。以下是核心逻辑的伪代码示例：

```python
from pyatb import TB_model, Solver

# 1. 初始化 Tight-Binding 模型
# 需要指定 ABACUS 的输出目录和轨道文件
model = TB_model.from_abacus(
    path='./OUT.ABACUS',           # ABACUS 输出目录
    binary=False,                  # 是否为二进制格式，csr 通常为文本
    orb_files=['Si_gga_8au_60Ry_2s2p1d.orb'] # 必须与 STRU 中使用的轨道文件一致
)

# 2. 设置求解器
solver = Solver(model)

# 3. 设置 K 点路径 (高对称点)
# 坐标需参考晶体结构，例如面心立方
kpath = [
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.0, 0.5],  # X
    # ... 添加更多点
]
labels = ['G', 'X']
solver.set_kpath(kpath, npts=100, labels=labels)

# 4. 计算能带与轨道投影
# projection=True 开启 Fat Band 计算
results = solver.solve_bands(projection=True)

# 5. 绘图与保存
# PYATB 通常提供内置绘图功能，或提取 results 数据自行绘制
# results['band'] 包含本征值
# results['proj'] 包含投影权重
print("Fat band calculation done.")
```

### 3.2.2 结果分析
运行脚本后，你将得到包含轨道权重的能带数据。在绘图中，能带的粗细（或颜色）代表了特定轨道（如 Fe 的 $3d$ 轨道）的贡献。这对于分析磁性材料、强关联体系的电子结构至关重要。

---

## 3.3 能带反折叠 (Band Unfolding)

当你计算一个 $2\times2\times2$ 的超胞（例如为了引入一个空位缺陷）时，原本清晰的能带会被折叠进变小了 8 倍的布里渊区中，变得杂乱无章。Band Unfolding 技术可以将这些能带“展开”回原胞的布里渊区，使其与实验 ARPES 结果可比。

### 3.3.1 物理图像
*   **Supercell (SC)**: 实际计算的体系（大实空间，小倒空间）。
*   **Primitive Cell (PC)**: 理想的参考体系（小实空间，大倒空间）。
*   **目标**: 计算 SC 的波函数在 PC 波函数基组上的投影权重（Spectral Weight）。

### 3.3.2 操作流程

1.  **ABACUS 计算**: 对**超胞**结构进行 SCF 计算，并按照 3.1 节的方法输出 HR/SR 矩阵。
2.  **PYATB 设置**: 在 Python 脚本中，除了加载模型，还需要定义从原胞到超胞的变换矩阵。

```python
# ... (加载模型步骤同上，注意这里加载的是超胞的矩阵) ...

# 定义超胞变换矩阵 (Supercell Matrix)
# 例如 2x2x2 超胞
sc_matrix = [
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2]
]

# 设置反折叠求解
# 注意：K 点路径应该是针对“原胞”布里渊区的高对称点
solver.set_kpath(kpath_primitive, npts=100)

# 运行反折叠
# unfolding=True 激活反折叠算法
unfold_data = solver.solve_unfolding(sc_matrix=sc_matrix)

# 绘图
# 结果通常是 (E, k) 空间的谱权重图 (Spectral Function)
# 颜色深的地方代表“有效能带”，模糊的地方代表对称性破缺导致的展宽
```

---

## 附录：常见问题与避坑指南

### 1. 对称性陷阱 (Symmetry Trap)
*   **现象**: PYATB 报错提示矩阵维度不匹配或 K 点映射错误。
*   **原因**: 在进行能带反折叠时，如果 ABACUS 计算超胞时自动通过对称性简化了 K 点或调整了结构，可能导致映射失败。
*   **对策**: 虽然 SCF 计算通常开启对称性，但在处理复杂的反折叠任务时，如果遇到问题，尝试在 ABACUS `INPUT` 中设置 `symmetry 0` 以关闭对称性操作，确保输出的矩阵完全对应输入的 `STRU` 几何。

### 2. 文件依赖 (File Dependencies)
*   **现象**: PYATB 运行时报错 `FileNotFound` 或 `Orbital mismatch`。
*   **原因**: 很多人只拷贝了 `OUT.ABACUS` 文件夹，却忘记了 `.orb` 轨道文件。
*   **对策**: PYATB 需要读取 `.orb` 文件中的截断半径和轨道形状来计算重叠积分修正。**务必确保 `.orb` 文件与 Python 脚本可见，且文件名与 `STRU` 中定义的一模一样。**

### 3. 流程中断 (Workflow Break)
*   **误区**: 试图用 `nscf` 的结果去跑 PYATB。
*   **纠正**: PYATB 需要的是实空间哈密顿量矩阵（HR）。这通常是在 `scf` 计算结束时通过 `out_mat_hs2` 输出的。`nscf` 主要用于在倒空间计算本征值，通常不输出实空间稀疏矩阵。请务必分清这两个流程。

### 4. 物理图像：原胞 vs 超胞
*   **提示**: 在做反折叠时，你的 K 点路径（`kpath`）必须是写在**原胞**的倒格矢基底上的，而不是超胞的。否则你画出来的能带图横坐标会完全错误。

---

## 附录：进阶学习指南

### 延伸阅读 (Extended Reading)
完成本教程后，你已掌握了能带计算的核心技能。若想在计算材料学领域进一步深造，建议探索以下主题：
1.  **态密度 (DOS) 计算**：能带图展示了能量随动量的变化，而 DOS 则反映了能量状态的分布密度，两者结合可提供完整的电子结构图景。
2.  **结构弛豫 (Relaxation)**：本教程假设你已拥有合理的结构。在实际科研中，学会如何使用 ABACUS 进行离子步优化以获得最低能量构型是必不可少的。
3.  **有效质量与迁移率**：利用能带顶或带底的曲率，可以进一步计算载流子的有效质量，进而预测材料的电学输运特性。
4.  **Wannier 函数与拓扑分析**：对于希望研究拓扑绝缘体或反常霍尔效应的读者，学习如何通过 ABACUS 接口构建 Wannier 函数是下一步的重点。

### 故障排除 (Troubleshooting)
计算过程中遇到问题是常态，请参考以下通用建议：
- **收敛性检查**：如果 SCF 不收敛，请尝试调整混合参数（Mixing）或增加 K 点密度。记住：Garbage In, Garbage Out。
- **对称性问题**：在 NSCF 计算中，确保 K 点路径经过体系的高对称点。可以使用工具如 Seek-path 辅助确定。
- **基组与赝势**：确保使用的数值原子轨道（Basis）和赝势（Pseudo）在目标体系下经过了充分验证。

### 官方资源
- **ABACUS 官方文档**：包含最详尽的参数说明与技术细节。
- **Bohrium 社区**：提供丰富的实战案例与在线计算环境。
- **GitHub 仓库**：关注 ABACUS 与 PYATB 的更新动态，积极参与社区讨论。

愿你在计算材料学的道路上不断突破，探索未知！
