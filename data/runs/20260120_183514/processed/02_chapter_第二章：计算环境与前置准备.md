# 第二章：计算环境与前置准备

在进行弹性常数（Elastic Constants）这一高阶物理性质计算之前，构建一个稳定、精确的计算流水线至关重要。弹性常数的计算本质上是对晶体施加微小形变后，通过响应的应力或能量变化来拟合二阶导数。这一过程对初始结构的平衡状态异常敏感——**任何未消除的残余应力（Residual Stress）都会像噪声一样叠加在物理响应上，导致计算结果严重失真**。

本章将指导你完成两个核心任务：配置基于 ABACUS 和 Python 的混合工具链，以及获得一个完美的“零应力”金刚石（Diamond）初始结构。

---

## 2.1 工具链准备

弹性常数计算是一个典型的“生成-计算-后处理”复合工作流。ABACUS 负责核心的量子力学计算（DFT），而结构的生成与数据的拟合则由 Python 脚本生态完成。

### 2.1.1 核心软件环境
请确保你的计算集群或工作站已正确安装以下组件：

1.  **ABACUS 主程序**: 建议使用 v3.0 及以上版本（支持更完善的 `cal_stress` 功能）。
    *   *注意*: 本教程示例基于 `PW`（平面波）或 `LCAO`（原子轨道）基组均可，但考虑到金刚石属于硬材料，且弹性常数对基组完备性敏感，**推荐在初始学习阶段使用平面波（PW）基组**以减少基组重叠误差的影响。
2.  **Python 环境**: 推荐 Python 3.8+。
3.  **Pymatgen**: 这是一个强大的材料分析库，我们将利用它来处理晶体结构的张量变换。
    ```bash
    # 在终端中执行安装
    pip install pymatgen numpy
    ```

### 2.1.2 辅助脚本获取
本教程采用 **Stress-Strain（应力-应变）** 方法进行计算。为了简化操作，我们提供了两个非 ABACUS 内置的专用 Python 脚本。请从 ABACUS 官方文档仓库或本教程附件中下载：

*   **`gene_dfm.py`**: 用于生成形变结构。它会读取平衡态结构，根据指定的应变模式（Strain Modes）和应变大小（$\delta$），生成一系列发生畸变的结构文件。
*   **`compute_dfm.py`**: 用于后处理。它会读取所有形变任务的应力输出，通过线性回归拟合胡克定律（Hooke's Law），计算出刚度矩阵（Stiffness Matrix）。

### 2.1.3 推荐的目录结构
为了避免文件混乱，建议按照以下树状结构组织你的工作目录（以 Diamond 为例）：

```text
Diamond_Elastic/
├── 00_relax/          # 本章重点：用于初始结构几何优化
│   ├── INPUT
│   ├── STRU
│   └── KPT
├── 01_elastic_run/    # 下一章内容：用于存放形变任务
│   ├── gene_dfm.py    # 放入辅助脚本
│   ├── compute_dfm.py # 放入辅助脚本
│   └── ...
└── pseudo/            # 存放碳原子的赝势文件 (C.pbe-n-kjpaw_psl.1.0.0.UPF 等)
```

---

## 2.2 初始结构的几何优化 (Pre-relaxation)

这是整个计算流程中**最关键**的一步。我们需要对原始的金刚石晶胞进行变胞弛豫（Variable-Cell Relaxation），使其内部应力降至最低（理想目标是 $< 0.1$ kBar）。

### 2.2.1 准备 STRU 文件 (Diamond 8原子超胞)
虽然金刚石的原胞仅包含 2 个原子，但在实际计算中，为了方便施加不同方向的应变，我们通常使用包含 8 个原子的常规立方晶胞。

**文件路径**: `Diamond_Elastic/00_relax/STRU`

```text
ATOMIC_SPECIES
C 12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF

LATTICE_CONSTANT
3.567  # 初始猜测值，稍后ABACUS会自动优化它

LATTICE_VECTORS
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0

ATOMIC_POSITIONS
Direct
C # Label
0.00 0.00 0.00 1 1 1
0.00 0.50 0.50 1 1 1
0.50 0.00 0.50 1 1 1
0.50 0.50 0.00 1 1 1
0.25 0.25 0.25 1 1 1
0.25 0.75 0.75 1 1 1
0.75 0.25 0.75 1 1 1
0.75 0.75 0.25 1 1 1
```
*注：`1 1 1` 表示允许原子在三个方向上移动。*

### 2.2.2 编写 INPUT 文件
在此阶段，我们必须开启 `cell-relax`（或 `vc-relax`），允许 ABACUS 调整晶格常数以释放应力。

**文件路径**: `Diamond_Elastic/00_relax/INPUT`

```text
INPUT_PARAMETERS
# 1. General Parameters
suffix          diamond_relax
calculation     cell-relax    # 核心参数：变胞弛豫，优化晶格矢量和原子位置
symmetry        1             # 保持晶体对称性，提高计算效率

# 2. Basis Set (推荐使用 PW 基组进行高精度 Stress 计算)
basis_type      pw
ecutwfc         60            # 60 Ry (~816 eV)，硬材料需要较高的截断能

# 3. Iteration & Convergence
scf_thr         1e-8          # 自洽场收敛精度，弹性常数计算建议设高
scf_nmax        100
relax_nmax      100           # 几何优化最大步数

# 4. Stress & Force (风险提示：必须开启)
cal_force       1             # 计算原子受力
cal_stress      1             # 核心参数：计算应力张量。这是弹性常数计算的充要条件！

# 5. Smearing (绝缘体/半导体设置)
smearing_method gaussian
smearing_sigma  0.01
```

> **专家提示 (Expert Tip)**: 
> *   **`cal_stress 1`**: 无论是当前的预处理还是后续的形变计算，此参数必须始终开启。如果关闭，ABACUS 将不会输出应力张量，后续的 `compute_dfm.py` 将无法提取数据。
> *   **`calculation cell-relax`**: 仅在本步骤（00_relax）使用。

### 2.2.3 运行与检查
提交任务运行完成后，检查输出文件 `OUT.diamond_relax/`。

1.  **检查收敛性**: 确认任务正常结束。
2.  **提取优化结构**: ABACUS 会将优化后的结构写入 `OUT.diamond_relax/STRU_ION_D`。
    *   **操作**: 将此 `STRU_ION_D` 复制到 `01_elastic_run/` 目录下，并重命名为 `STRU`。
    *   **重要**: 这个新的 `STRU` 将作为后续生成所有形变结构的**唯一母本**。

---

## 2.3 核心逻辑预警：从弛豫到形变

在进入下一章之前，请务必理解以下**计算逻辑的切换**，这是初学者最容易犯错的地方：

| 阶段 | 任务目标 | INPUT `calculation` 设置 | 物理含义 |
| :--- | :--- | :--- | :--- |
| **阶段一 (本章)** | **消除内应力** | `cell-relax` (或 vc-relax) | **改变晶格**，寻找能量最低点，使 $\sigma_{total} \to 0$。 |
| **阶段二 (下一章)** | **计算弹性响应** | **`relax`** (固定晶胞) | **固定晶格**（保持人为施加的应变），仅弛豫原子位置。 |

**危险操作提示**: 
如果在生成了形变结构（如拉伸了 1%）后，错误地使用了 `calculation cell-relax`，ABACUS 会自动优化晶格，把我们辛苦施加的 1% 应变“优化”掉，使其变回完美的金刚石结构。这样计算出的应力将接近于零，导致弹性常数计算结果为零或完全错误。

**结论**: 
1.  本章产生的 `STRU` 是完美的、无应力的。
2.  下一章我们将使用脚本对这个 `STRU` 施加 $\delta \in \{-0.01, -0.005, 0.005, 0.01\}$ 的应变。
3.  在计算那些形变结构时，**严禁**改变晶胞大小。