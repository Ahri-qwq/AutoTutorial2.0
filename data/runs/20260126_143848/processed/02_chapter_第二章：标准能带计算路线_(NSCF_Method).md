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