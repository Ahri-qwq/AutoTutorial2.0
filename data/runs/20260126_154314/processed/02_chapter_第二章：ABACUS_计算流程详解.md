# 第二章：ABACUS 计算流程详解

在第一章中，我们已经完成了软件的安装与环境配置。本章将进入实战核心，以**铁电材料 PbTiO$_3$ 的自发极化计算**为例，详细演示基于 ABACUS LCAO 基组的完整工作流。

计算铁电极化（Berry Phase 方法）是很多初学者的“噩梦”，因为它不仅涉及复杂的物理概念（现代极化理论），在操作上还必须严格遵循 **“先 SCF 后 NSCF”** 的分步流程。本章将带你避开这些常见的坑。

---

## 2.1 体系建模与对称性检查

计算自发极化（Spontaneous Polarization）的物理本质是计算两个状态之间的极化差值：**中心反演对称相（参考相）** 与 **铁电相（畸变相）**。

**核心原则**：你用于 Berry Phase 计算的结构，**必须**是破坏了中心反演对称性的。如果输入的结构具有中心反演对称性（例如理想的立方钙钛矿结构），理论计算出的电子极化值模数应为 0。

### 2.1.1 准备 STRU 文件
我们以四方相 PbTiO$_3$ 为例，这是一个经典的铁电体，Ti 原子沿 z 轴发生了位移。

**文件名**: `STRU`

```abacus
ATOMIC_SPECIES
Pb 207.2 Pb_ONCV_PBE-1.0.upf Pb_gga_7au_100Ry_4s2p2d.orb
Ti 47.87 Ti_ONCV_PBE-1.0.upf Ti_gga_7au_100Ry_4s2p2d.orb
O  15.99 O_ONCV_PBE-1.0.upf  O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
7.37  // 约为 3.90 Angstrom converted to Bohr (1 Ang = 1.8897 Bohr)

LATTICE_VECTORS
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.060000 // c/a > 1，四方拉伸

ATOMIC_POSITIONS
Direct  // 使用分数坐标

Pb // Pb at corner
0.000000 0.000000 0.000000 1 1 1

Ti // Ti at center, displaced along z
0.500000 0.500000 0.540000 1 1 1 

O // O1
0.500000 0.500000 0.110000 1 1 1

O // O2
0.500000 0.000000 0.610000 1 1 1

O // O3
0.000000 0.500000 0.610000 1 1 1
```

> **专家提示**：
> 1. 请确保 `*.upf` (赝势) 和 `*.orb` (轨道) 文件名与你目录下的实际文件一致。
> 2. 上述坐标中，Ti 和 O 均相对于理想位置发生了 z 方向的位移，这是产生极化的来源。

---

## 2.2 第一步：自洽计算 (SCF)

Berry Phase 计算属于非自洽计算（NSCF），它依赖于基态的布洛赫波函数。因此，**必须**先进行一次标准的自洽场（SCF）计算来获取收敛的电荷密度。

### 2.2.1 准备 INPUT 文件 (SCF)
创建一个名为 `SCF` 的文件夹，将 `STRU` 和 `KPT` 以及赝势/轨道文件放入其中。

**文件名**: `INPUT`

```abacus
INPUT_PARAMETERS
# ---------------------------
# General Parameters
# ---------------------------
suffix          pto_scf      # 输出文件后缀
calculation     scf          # 核心参数：执行自洽计算
basis_type      lcao         # 使用原子轨道基组
symmetry        1            # 开启对称性分析

# ---------------------------
# Electronic Structure
# ---------------------------
ecutwfc         100          # 平面波截断能 (Ry)，用于构建格点积分
scf_nmax        100          # 最大迭代步数
scf_thr         1e-7         # 收敛精度，建议比默认值高

# ---------------------------
# LCAO Specific
# ---------------------------
mixing_type     pulay        # 电荷密度混合方法
mixing_beta     0.7          # 混合因子

# ---------------------------
# I/O
# ---------------------------
out_chg         1            # 核心参数：必须输出电荷密度文件！
cal_force       0            # 极化计算通常固定离子位置
cal_stress      0
```

### 2.2.2 准备 KPT 文件
对于绝缘体，K 点网格不需要极其密集，但为了保证 Berry Phase 的精度，建议适当加密。

**文件名**: `KPT`
```abacus
K_POINTS
0
Gamma
6 6 6 0 0 0
```

### 2.2.3 运行计算
执行 ABACUS（假设使用 MPI 并行）：
```bash
mpirun -n 4 abacus
```
计算完成后，检查 `OUT.pto_scf` 目录，确保生成了电荷密度文件（通常为 `SPIN1_CHG.cube` 或二进制密度文件，取决于版本和设置）。这是下一步的“燃料”。

---

## 2.3 第二步：非自洽 Berry Phase 计算 (NSCF)

现在我们进入关键步骤。我们需要读取上一步的电荷密度，进行一次非自洽计算来求解 Berry Phase。

### 2.3.1 准备工作
1. 新建文件夹 `Berry_Phase`。
2. 将 `SCF` 步骤中的 `STRU`, `KPT`, `*.upf`, `*.orb` 复制过来。
3. **关键操作**：将 `SCF` 计算生成的 `OUT.pto_scf` 文件夹路径记录下来，或者将其复制到当前目录下。

### 2.3.2 准备 INPUT 文件 (NSCF)

**文件名**: `INPUT`

```abacus
INPUT_PARAMETERS
# ---------------------------
# General
# ---------------------------
suffix          pto_berry
calculation     nscf         # 核心参数：非自洽计算
basis_type      lcao

# ---------------------------
# Restart & Initialization
# ---------------------------
init_chg        file         # 核心参数：从文件读取电荷密度
read_file_dir   ../SCF/      # 指向第一步 SCF 的输出目录(包含 OUT.pto_scf)

# ---------------------------
# Berry Phase Settings
# ---------------------------
berry_phase     1            # 开启 Berry Phase 计算
# g_evec        1            # 部分版本可能需要显式开启特征向量计算

# ---------------------------
# Risk Warning / Direction
# ---------------------------
# 注意：ABACUS 计算 Berry Phase 的具体方向控制参数（如 berry_phase_dir）
# 可能随版本更新而变化。
# 请务必查阅官方文档中关于 "Berry Phase" 的最新说明，
# 确认是否需要通过 KPT 路径或特定参数指定极化方向。
```

> **重要参数解析**：
> *   `calculation nscf`: 告诉 ABACUS 不要更新电荷密度，仅求解波函数。
> *   `init_chg file`: 强制程序读取外部电荷密度，而不是从原子叠加密度开始。
> *   `berry_phase 1`: 激活现代极化理论计算模块。ABACUS 的 LCAO 模块在计算 Berry Curvature 方面具有天然优势，通常能以较低的计算成本获得与有限差分法（Wannier90）相当的精度。

### 2.3.3 运行与输出
运行计算后，查看输出文件（通常在日志文件或特定的 `*.berry` 文件中）。你将看到类似于以下的输出：

```text
Polarization along direction 1: ...
Polarization along direction 2: ...
Polarization along direction 3: ...
Total Polarization (C/m^2): ...
```

---

## 2.4 结果分析：多值性陷阱 (Modulo 2π)

当你拿到计算结果时，千万不要直接使用数值，必须理解 **Berry Phase 的多值性**。

### 2.4.1 极化量子 (Polarization Quantum)
Berry Phase 定义的极化强度 $P$ 并不是一个单值，而是一个晶格（Lattice of values）：
$$ P = P_0 + n \cdot \frac{e \mathbf{R}}{\Omega} $$
其中 $\frac{e \mathbf{R}}{\Omega}$ 被称为**极化量子**。这意味着计算出的极化值在相差整数倍个极化量子时是物理等价的。

### 2.4.2 如何计算自发极化
自发极化 $P_s$ 定义为：
$$ P_s = P_{\text{ferro}} - P_{\text{para}} $$
其中 $P_{\text{ferro}}$ 是你刚刚计算的铁电相极化，$P_{\text{para}}$ 是参考顺电相（中心对称结构）的极化。

**操作步骤**：
1. 计算铁电相的 $P_{\text{ferro}}$（包含电子部分和离子部分）。
2. 计算顺电相的 $P_{\text{para}}$（通常理论为 0，但建议计算以消除系统误差）。
3. 相减时，如果结果大得离谱，很可能是因为两个值处于不同的“分支”上。你需要加上或减去整数倍的极化量子，使得 $P_s$ 落在合理的物理范围内（通常对于钙钛矿，极化值在 0.1 ~ 1.0 C/m$^2$ 之间）。

---

## 2.5 本章小结

本章我们演示了 ABACUS 中最经典的 **SCF -> NSCF** 两步工作流。
1. **建模**：确保结构破坏中心反演对称性。
2. **SCF**：计算并保存基态电荷密度 (`out_chg 1`)。
3. **NSCF**：读取电荷密度 (`init_chg file`) 并开启 `berry_phase 1`。

掌握了这一流程，你就掌握了 ABACUS 进行大部分电子性质分析（如能带、态密度、光学性质）的通用模板。下一章，我们将深入探讨如何利用 ABACUS 进行结构优化与动力学模拟。