# 第二章：基态自洽计算与矩阵生成 (ABACUS)

在进行光学性质计算（如介电函数、光导率）之前，我们必须首先获得高质量的电子基态波函数。与常规的能带计算或结构优化不同，结合 pyatb 进行光学计算时，ABACUS 的任务不仅仅是输出能量，更重要的是生成描述电子跃迁所需的**稀疏矩阵文件**（哈密顿量、重叠矩阵、位置算符矩阵）。

本章将指导你配置 ABACUS 的输入文件，运行 SCF 计算，并正确提取后续步骤所需的关键数据。

## 2.1 结构建模与输入文件准备

我们将以晶态 $\mathbf{SiO_2}$ 为例。为了计算光学性质，我们需要准备标准的 ABACUS 四大输入文件：`STRU`、`KPT`、`INPUT` 以及赝势/轨道文件。

### 2.1.1 晶体结构 (STRU)
在 `STRU` 文件中，我们定义晶格常数、晶格矢量和原子位置。

> **⚠️ 专家提示 (Critical)**：
> 后续使用 pyatb 时，**必须手动将 `STRU` 中的晶格信息复制到 pyatb 的输入文件中**。请务必注意长度单位。ABACUS 默认使用 **Bohr**，而许多后处理工具可能默认使用 Angstrom。本教程案例统一使用 **Bohr**。

**文件示例：`STRU`**
```abacus
ATOMIC_SPECIES
Si 28.086  Si_ONCV_PBE-1.0.upf
O 15.999   O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_7au_100Ry_2s2p1d.orb
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897261246257702  # 1.8897 Bohr approx 1 Angstrom

LATTICE_VECTORS
7.1199998856 0.0 0.0 
0.0 7.1199998856 0.0 
0.0 0.0 7.1199998856 

ATOMIC_POSITIONS
Cartesian
Si
0.0
8
0.000000000000 0.000000000000 0.000000000000 1 1 1
... (其余原子坐标略)
```

### 2.1.2 K 点设置 (KPT)
光学性质的计算涉及布里渊区积分。为了获得平滑的介电函数曲线，**通常需要比常规 SCF 计算更密的 K 点网格**。虽然本教程演示使用 `6 6 6`，但在实际科研中，可能需要 `20 20 20` 或更高。

**文件示例：`KPT`**
```abacus
K_POINTS
0
Gamma
6 6 6 0 0 0
```

## 2.2 关键控制参数设置 (INPUT)

这是本章的核心。为了让 ABACUS 输出 pyatb 所需的矩阵，我们需要在 `INPUT` 文件中添加特定的控制参数。

### 2.2.1 基础参数
*   `calculation`: 设置为 `scf` (自洽场计算)。
*   `basis_type`: 必须设置为 `lcao` (线性组合原子轨道)，因为 pyatb 是基于紧束缚模型（Tight-Binding）原理工作的，依赖于局域基组生成的稀疏矩阵。

### 2.2.2 矩阵输出参数 (核心)
请在 `INPUT` 文件中确保包含以下参数：

*   **`out_mat_hs2`**: 设置为 `1`。
    *   **含义**: 输出稀疏格式（CSR）的哈密顿量矩阵 ($H$) 和重叠矩阵 ($S$)。这是计算能带本征值的基础。
*   **`out_mat_r`**: 设置为 `1`。
    *   **含义**: 输出位置算符矩阵（偶极矩阵）。这是计算电子跃迁几率（速度矩阵元）的必要条件。
*   **`out_chg`**: 设置为 `1`。
    *   **含义**: 输出电荷密度文件。虽然 pyatb 直接读取矩阵，但保存电荷密度便于后续分析或作为非自洽计算的起点。

**文件示例：`INPUT`**
```abacus
INPUT_PARAMETERS

# --- 基础设置 ---
suffix          silica
calculation     scf
basis_type      lcao
ks_solver       genelpa  # 或其他对角化求解器
symmetry        0        # 建议关闭对称性以避免矩阵生成中的潜在旋转问题

# --- 电子结构参数 ---
ecutwfc         100      # Rydberg
scf_thr         1e-8     # 精度要求较高
smearing_method gaussian
smearing_sigma  0.01

# --- 核心：矩阵输出 ---
out_chg         1        # 输出电荷密度
out_mat_hs2     1        # 输出 H 和 S 稀疏矩阵 (Critical for pyatb)
out_mat_r       1        # 输出 r 矩阵 (Critical for Optics)
```

## 2.3 运行计算与结果检查

准备好所有文件后，提交 ABACUS 任务。

```bash
# 设置线程数并运行 (示例使用 MPI)
export OMP_NUM_THREADS=1
mpirun -np 16 abacus
```

### 2.3.1 检查输出日志与提取 Fermi Energy
计算完成后，进入输出目录（例如 `OUT.silica`）。我们需要从日志文件 `running_scf.log` 中提取两个至关重要的参数，用于下一章 pyatb 的输入。

1.  **检查收敛性**: 确认日志末尾有 `TOTAL Time`，且 SCF 迭代达到收敛精度。
2.  **提取占据带数目 (Occupied Bands)**。
3.  **提取费米能级 (Fermi Energy)**。

> **🛑 关键步骤 (Data Transfer)**：
> pyatb **不会**自动读取 ABACUS 的费米能级。你必须手动查找并填入。如果填错，会导致能带占据数计算错误，进而导致光谱结果完全错误（例如绝缘体变成了金属）。

**操作指令：**
```bash
grep 'occupied bands' OUT.silica/running_scf.log
grep 'E_Fermi' OUT.silica/running_scf.log | tail -n 1
```

**预期输出示例：**
```text
occupied bands = 64
E_Fermi            0.4070749075         5.5385382545
```
*   注意：`E_Fermi` 通常有两列数值。第一列单位通常为 **Ry** (Rydberg)，第二列为 **eV**。
*   **记录**: 费米能级 $E_f \approx 5.5385$ eV。我们在下一章配置 pyatb 时将使用此数值。

### 2.3.2 确认矩阵文件
检查 `OUT.silica` 目录下是否生成了以下文件：
*   `data-HR-sparse_SPIN0.csr` (哈密顿量矩阵)
*   `data-SR-sparse_SPIN0.csr` (重叠矩阵)
*   `data-rR-sparse.csr` (位置/偶极矩阵)

如果这些文件缺失，请检查 `INPUT` 中 `out_mat_*` 参数是否设置正确。

### 2.3.3 数据转移
由于 ABACUS 的输出在子文件夹中，而 pyatb 通常在独立目录运行，建议将矩阵文件复制到 pyatb 的工作目录：

```bash
mkdir -p pyatb_OpticConductivity
cp OUT.silica/data-*-sparse*.csr pyatb_OpticConductivity/
# 同时也建议备份 INPUT 和 STRU 以便查阅晶格信息
cp INPUT STRU pyatb_OpticConductivity/
```

---

## 2.4 预备知识：光学性质单位换算逻辑

在进入下一章使用 pyatb 计算之前，我们需要理解后续 Python 后处理脚本中将涉及的单位换算，特别是**吸收系数 (Absorption Coefficient)**。

pyatb 计算出的介电函数 $\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)$ 是无量纲的。但在计算吸收系数 $\alpha(\omega)$ 时，公式如下：

$$
\alpha(\omega) = \frac{\omega}{c} \sqrt{2 \left( \sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1 \right)}
$$

在编写 Python 脚本时，我们需要注意以下物理常数的单位匹配：

1.  **能量与频率**:
    *   通常我们使用能量 $E$ (单位 eV) 作为横坐标。
    *   频率 $\omega = E / \hbar$。
    *   $\hbar \approx 4.1357 \times 10^{-15} \text{ eV}\cdot\text{s}$。
    *   因此，$\omega (\text{s}^{-1}) \approx E (\text{eV}) \times 2.418 \times 10^{14}$。

2.  **光速 $c$ 与 吸收系数单位**:
    *   吸收系数的常用单位是 $\text{cm}^{-1}$。
    *   为了得到 $\text{cm}^{-1}$，光速 $c$ 必须使用 **cm/s** 为单位。
    *   $c \approx 3.0 \times 10^{10} \text{ cm/s}$。
    *   *注意*：如果代码中使用 $c \approx 3.0 \times 10^8 \text{ m/s}$，计算出的吸收系数单位将是 $\text{m}^{-1}$，数值上会相差 100 倍。

**Python 代码逻辑预演：**
```python
# 这里的常数是为了将 eV 转换为 s^-1，并将结果归一化到 cm^-1
# 2.418e14 来自 1/hbar
# 3.0e8 是光速 (m/s)。如果除以 m/s，结果是 m^-1。
# 若要得到 cm^-1，需再除以 100，或者使用 3.0e10 cm/s。
# 请在下一章的脚本中仔细检查此处的数量级。
absorption_coefficient = ... * (2.418*10**14) / (3.0*10**8) 
```

至此，我们已经完成了 ABACUS 部分的所有准备工作，生成了核心矩阵数据并提取了费米能级。下一章，我们将使用这些数据驱动 pyatb 进行线性光学性质计算。