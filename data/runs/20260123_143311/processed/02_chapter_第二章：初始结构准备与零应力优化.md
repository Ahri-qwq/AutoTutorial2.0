# 第二章：初始结构准备与零应力优化

> **本章核心逻辑**：“Garbage in, garbage out”。
> 弹性常数（Elastic Constants）本质上是能量对格点应变的二阶导数，或者是应力对应变的一阶导数。如果初始结构的残余应力（Residual Stress）未被消除，后续施加微小应变（通常仅为 $\pm 0.5\% \sim 1\%$）引起的应力变化就会被初始误差掩盖，导致计算出的弹性模量出现严重偏差，甚至出现非物理的结果（如负模量）。
>
> 本章将指导你如何获得一个完美的、无残余应力的初始构型，这是后续所有计算的基石。

---

## 2.1 晶胞选择与标准化

在计算弹性常数之前，第一步是选择合适的晶胞（Unit Cell）。对于晶体材料，我们通常面临 **原胞 (Primitive Cell)** 和 **惯用胞 (Conventional Cell)** 的选择。

### 2.1.1 原胞 vs. 惯用胞
*   **原胞 (Primitive Cell)**: 包含原子数最少，计算效率最高。但其晶格矢量通常不与笛卡尔坐标轴（x, y, z）平行。
*   **惯用胞 (Conventional Cell)**: 原子数较多，计算成本较高。但其晶格矢量通常与笛卡尔坐标轴对齐，且更能直观反映晶体的对称性（如立方晶系）。

**教授建议**: 
在计算弹性常数时，**强烈建议使用惯用胞**，特别是对于立方（Cubic）、四方（Tetragonal）和正交（Orthorhombic）晶系。
原因如下：
1.  **物理定义的直观性**: 弹性刚度张量 $C_{ij}$ 通常是基于笛卡尔坐标系定义的。如果使用原胞，你得到的 $C_{ij}$ 也是基于原胞基矢的，需要进行复杂的张量旋转才能与实验值或数据库（如 Materials Project）对比。
2.  **IEEE 标准**: Materials Project 数据库通常将晶体旋转至 IEEE 176-1987 标准规定的取向。使用与坐标轴对齐的惯用胞可以天然满足这一标准。

### 2.1.2 实战：生成标准化的 Si 结构
我们将使用 `abacustest` 及其依赖的 ASE (Atomic Simulation Environment) 来生成标准的硅（Si）惯用胞结构。

创建一个 Python 脚本 `generate_si.py`：

```python
from ase.build import bulk
from ase.io import write

# 生成 Si 的惯用胞 (cubic=True)
# a=5.43 是初始猜测值，后续我们会通过 ABACUS 优化它
atoms = bulk('Si', 'diamond', a=5.43, cubic=True)

# 保存为 CIF 格式，供 abacustest 使用
write('Si.cif', atoms)
print("Si.cif (Conventional Cell) has been generated.")
```

运行该脚本后，你将获得一个晶格矢量与 x, y, z 轴平行的 `Si.cif` 文件。

---

## 2.2 消除残余应力 (高精度弛豫)

获得初始结构后，必须使用 ABACUS 进行**变胞优化 (Variable-Cell Relaxation)**。我们的目标是将晶胞的残余应力降至忽略不计的水平（通常建议 $< 1 \text{ kBar}$）。

### 2.2.1 使用 abacustest 准备任务
`abacustest` 封装了繁琐的输入文件准备过程。我们将使用 `lib-prepare` 命令来生成优化任务。

**⚠️ 风险提示**: 
`abacustest` 的准备命令会**直接删除**目录下已有的同名文件夹（如 `Si/`），请务必在操作前备份重要数据。

在终端执行以下命令：

```bash
# 假设你已经有了 Si.cif 和相应的赝势/轨道文件
# 环境变量 ABACUS_PP_PATH 和 ABACUS_ORB_PATH 需提前设置好

abacustest lib-prepare \
    -f Si.cif \
    --ftype cif \
    --jtype cell-relax \
    --lcao \
    --folder-syntax Si_opt
```

*   `-f Si.cif`: 指定输入结构。
*   `--jtype cell-relax`: **关键参数**。指定任务类型为晶胞弛豫，这会自动设置 `INPUT` 中的 `calculation` 参数。
*   `--lcao`: 使用 LCAO 基组（效率更高，适合高通量）。如需使用平面波，请去掉此参数并调整相关设置。

### 2.2.2 核心参数解析 (INPUT)
进入生成的 `Si_opt` 文件夹，检查 `INPUT` 文件。为了保证弹性常数的准确性，我们需要手动检查或微调以下关键参数：

```plaintext
INPUT_PARAMETERS
# ... (省略基础参数)

# --- 核心计算控制 ---
calculation     cell-relax  # 进行变胞优化，同时优化晶格矢量和原子位置
basis_type      lcao        # 基组类型

# --- 力与应力计算 (必须开启) ---
cal_force       1           # 计算原子力
cal_stress      1           # 计算晶格应力 (对 cell-relax 默认开启，但显式写出是个好习惯)

# --- 收敛标准 (至关重要) ---
# 建议将应力收敛标准设置得比默认值更严格
stress_thr      1.0         # 单位: kBar (1 kBar = 0.1 GPa)。
                            # 默认值通常为 10 kBar 左右，对于弹性计算这太大了。
                            # 建议设为 1.0 甚至 0.1。

force_thr_ev    0.001       # 单位: eV/Angstrom。
                            # 原子受力的收敛标准，建议 1e-3 或更小。

# --- 迭代控制 ---
relax_nmax      100         # 最大离子步数，防止不收敛时无限计算
out_stru        1           # 输出优化后的结构文件 (STRU_ION_D)
```

**专家解读**:
*   **`stress_thr`**: 这是本章最重要的参数。如果残余应力是 5 kBar (0.5 GPa)，而你后续施加 0.5% 应变产生的应力变化可能只有几 GPa，那么 0.5 GPa 的基线误差是不可接受的。
*   **`out_stru`**: 确保设置为 1。优化完成后，ABACUS 会生成 `STRU_ION_D` 文件，这就是我们下一章需要的“零应力完美结构”。

### 2.2.3 提交与检查
提交计算（根据你的集群环境使用 `sbatch` 或直接运行）。计算完成后，检查 `OUT.Si_opt/running_cell-relax.log` 的最后几行，确认：
1.  任务是否正常结束 (`DONE`).
2.  最终的应力值是否小于 `stress_thr`。

---

## 2.3 进阶理解：原子内部弛豫 (Internal Relaxation)

在进入下一章之前，必须澄清一个容易混淆的物理概念，这直接关系到你计算的是哪种弹性常数。

### 2.3.1 什么是内部弛豫？
当你对晶胞施加一个宏观应变（例如将晶格拉长 1%）时：
1.  **均匀形变 (Affine Deformation)**: 晶格矢量发生变化，原子根据分数坐标随之均匀移动。
2.  **内部弛豫 (Internal Relaxation)**: 在新的晶格形状下，原子原本的位置可能不再是受力平衡点（特别是对于原胞内包含多个原子的复式晶格，如 Si 的金刚石结构）。原子会受到微小的力，需要移动到一个新的能量最低点。

### 2.3.2 Clamped-ion vs. Relaxed-ion
*   **Clamped-ion (无弛豫) 弹性常数**: 
    *   对应过程：施加应变 -> **固定原子分数坐标** -> 计算应力。
    *   物理意义：对应极高频响应（原子来不及移动）。
    *   ABACUS 操作：在后续形变计算中使用 `--norelax` (即 `calculation = scf`)。
    *   **结果特征**: 通常数值**偏大**，不符合常规实验值。

*   **Relaxed-ion (完全弛豫) 弹性常数**:
    *   对应过程：施加应变 -> **固定晶格，优化原子位置** -> 计算应力。
    *   物理意义：对应静态或低频响应（实验测量的通常是这个）。
    *   ABACUS 操作：在后续形变计算中使用 `calculation = relax` (固定晶胞优化原子)。

**教程策略**: 
在本教程后续的第三章中，我们将默认计算 **Relaxed-ion** 弹性常数。
但在本章（第二章）的准备阶段，我们做的是 `cell-relax`，这是为了让**初始态**既没有宏观应力，原子也处于受力平衡位置。

---

## 2.4 对称性提示

对于我们使用的硅（立方晶系），理论上只有 3 个独立的弹性常数：
*   $C_{11}$ (纵向压缩)
*   $C_{12}$ (横向膨胀)
*   $C_{44}$ (剪切)

虽然 ABACUS 会计算出完整的 $6 \times 6$ 矩阵，但你会发现：
*   $C_{11} \approx C_{22} \approx C_{33}$
*   $C_{12} \approx C_{13} \approx C_{23}$
*   $C_{44} \approx C_{55} \approx C_{66}$
*   其他分量应接近于 0。

如果你的计算结果严重偏离这个规律（例如 $C_{11}$ 和 $C_{22}$ 相差巨大），通常意味着：
1.  **初始结构未充分优化**（残余应力过大）。
2.  k 点网格密度不够，破坏了对称性。
3.  基组精度不足（LCAO 截断半径过小或 PW 截断能过低）。

---

## 本章小结
1.  **结构选择**: 优先使用**惯用胞**，并确保晶格矢量与坐标轴对齐。
2.  **零应力优化**: 使用 `calculation = cell-relax`，并将 `stress_thr` 设为 `1.0` kBar 或更低。
3.  **产物**: 获得优化后的 `STRU_ION_D`，将其重命名并保存，作为第三章形变计算的唯一输入源。