# 第一章：理论基础与工具准备

在正式开始“跑数据”之前，我们必须清楚两件事：**我们在算什么（物理图像）**，以及**我们用什么算（工具分工）**。

本章将带你快速理解声子计算的核心——有限位移法，并手把手搭建基于 ABACUS 和 Phonopy 的计算环境。

## 1.1 晶格动力学与有限位移法

### 1.1.1 什么是声子？
在晶体材料中，原子并非静止不动，而是在平衡位置附近做微小的热振动。这种集体振动的能量量子化描述就是**声子（Phonon）**。声子谱（Phonon Spectrum）或色散关系（Dispersion Relation）描述了声子频率 $\omega$ 与波矢量 $\mathbf{q}$ 之间的关系，它是理解材料热导率、热容、超导电性以及相稳定性等宏观性质的钥匙。

### 1.1.2 为什么需要“有限位移法”？
为了得到声子谱，我们需要计算**力常数矩阵（Force Constants）**。在简谐近似（Harmonic Approximation）下，力常数对应于体系势能面在平衡位置的二阶导数。

**有限位移法（Finite Displacement Method, FDM）**，又称“冻结声子法”，是目前最直观、最鲁棒的计算方法。其核心逻辑如下：

1.  **扩胞（Supercell）**：建立一个足够大的超胞。这是为了截断原子间的长程相互作用，保证周期性边界条件下的力常数计算准确。
2.  **微扰（Displacement）**：人为地将超胞中的原子偏离平衡位置一个微小距离（通常为 0.01 Å 量级）。这一步破坏了原有的晶体对称性。
3.  **算力（Force Calculation）**：使用高精度的 DFT 引擎（这里是 **ABACUS**）计算微扰结构中所有原子的受力（Hellmann-Feynman Force）。
4.  **重构（Reconstruction）**：根据力与位移的关系（$F = -kx$ 的矩阵形式），反推力常数矩阵，进而构建动力学矩阵（Dynamical Matrix），最终对角化求解得到声子频率。

### 1.1.3 ABACUS 与 Phonopy 的分工
在这个工作流中，两个软件各司其职，缺一不可：

*   **Phonopy**（大脑）：
    *   负责根据晶体对称性，生成最少数量的必要微扰结构。
    *   负责收集 DFT 计算得到的力数据。
    *   负责计算力常数、声子谱、热力学性质并绘图。
*   **ABACUS**（引擎）：
    *   负责对 Phonopy 生成的每一个微扰结构进行自洽计算（SCF）。
    *   输出高精度的原子受力（Force）。

---

## 1.2 工具链环境搭建与案例准备

本教程将以经典的 **FCC 铝（Al）** 为例，演示完整的计算流程。

### 1.2.1 安装 Phonopy
Phonopy 是一个基于 Python 的开源软件。推荐使用 `pip` 或 `conda` 进行安装。

在终端执行以下命令（确保你已安装 Python 3）：
```bash
# 方法 A: 使用 pip (推荐)
pip install phonopy

# 方法 B: 从源码安装
git clone https://github.com/phonopy/phonopy.git
cd phonopy
python3 setup.py install
```
安装完成后，输入 `phonopy --version` 检查是否成功。

### 1.2.2 准备 ABACUS 环境
请确保你已在计算集群或本地机器上安装了 ABACUS（建议版本 3.2.x 及以上）。
*   你需要知道 ABACUS 可执行文件的路径（例如 `abacus` 或 `mpirun -np 4 abacus`）。
*   本教程使用 LCAO 基组进行计算，请确保已准备好相应的赝势文件（`.upf`）和轨道文件（`.orb`）。

### 1.2.3 案例文件结构
为了保持工作流清晰，建议按照以下目录结构组织文件。我们将以 `1_Al` 作为项目根目录。

```text
1_Al/
├── psp/                        # 存放赝势和轨道文件
│   ├── Al_ONCV_PBE-1.0.upf
│   └── Al_gga_7au_100Ry_4s4p1d.orb
├── STRU                        # 初始晶体结构文件
└── INPUT                       # ABACUS 输入参数模版
```

#### 1. 结构文件 (`STRU`)
这是 FCC Al 的原胞结构，也是我们计算的起点。
**文件名**: `STRU`
```abacus
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf upf201

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.88972612546  # 这里的常数通常用于将 Angstrom 转换为 Bohr

LATTICE_VECTORS
4.03459549706 0 0 #latvec1
0 4.03459549706 0 #latvec2
0 0 4.03459549706 #latvec3

ATOMIC_POSITIONS
Direct

Al #label
0 #magnetism
4 #number of atoms
0  0  0  m  0  0  0
0.5  0.5  0  m  0  0  0
0.5  0  0.5  m  0  0  0
0  0.5  0.5  m  0  0  0
```
> **注意**：此 `STRU` 文件中的原子位置必须是经过充分结构优化（Relaxation）后的位置。如果结构处于高能非平衡态，计算出的声子谱可能会出现虚频。

#### 2. 输入参数文件 (`INPUT`)
这是用于计算受力的核心参数文件。
**文件名**: `INPUT`
```abacus
INPUT_PARAMETERS
#Parameters (1.General)
suffix          Al-fcc
calculation     scf        # 必须进行自洽计算
esolver_type    ksdft
symmetry        1          # 开启对称性分析
pseudo_dir      ./psp      # 指定赝势目录
orbital_dir     ./psp      # 指定轨道目录
cal_stress      1          # 计算应力（推荐开启）
cal_force       1          # 【关键】必须开启力的计算！

#Parameters (2.Iteration)
ecutwfc         100        # 平面波截断能 (Ry)
scf_thr         1e-7       # 自洽收敛阈值 (推荐 1e-7 或更严)
scf_nmax        50

#Parameters (3.Basis)
basis_type      lcao       # 使用原子轨道基组
gamma_only      0

#Parameters (4.Smearing)
smearing_method mp         # 金属体系推荐 Methfessel-Paxton
smearing_sigma  0.015      # 展宽参数 (Ry)

#Parameters (5.Mixing)
mixing_type     pulay
mixing_beta     0.7
mixing_gg0      1.5
```

### 1.2.4 关键参数解析
在 `INPUT` 文件中，有几个参数对于声子计算至关重要：

1.  **`calculation scf`**: 我们只需要进行电子自洽计算，不需要离子移动（`relax` 或 `md`），因为原子位置是由 Phonopy 固定好的。
2.  **`cal_force 1`**: **这是最核心的参数**。它指示 ABACUS 计算并输出原子受力。如果没有这一行，Phonopy 将无法获取数据。
3.  **`scf_thr 1e-7`**: 声子计算对力的精度要求较高，建议将自洽收敛标准设置得比常规能量计算更严格（如 `1e-7` 或 `1e-8`）。
4.  **`stru_file`** (后续使用): 在后续章节实际运行微扰结构计算时，我们会用到 `stru_file` 参数来指定不同的结构文件（如 `STRU-001`），这一点在下一章会详细说明。

---
**本章小结**：
你已经理解了有限位移法的物理原理，并准备好了 ABACUS + Phonopy 的计算环境。手中的 `STRU` 和 `INPUT` 文件是通往声子谱的第一把钥匙。下一章，我们将正式进入实战，生成超胞并计算受力。