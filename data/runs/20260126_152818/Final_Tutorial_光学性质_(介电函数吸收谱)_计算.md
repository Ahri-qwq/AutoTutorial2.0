# ABACUS 第一性原理计算实战：从基态电子结构到线性光学性质

## 前言

**欢迎来到《ABACUS 第一性原理计算实战》系列教程。**

在材料科学的研究中，预测材料对光的响应能力——即光学性质——是连接理论模拟与实验表征的核心桥梁。ABACUS（原子算筹）作为国产高性能第一性原理计算软件，凭借其独具特色的数值原子轨道（LCAO）基组，在处理大规模体系和复杂电子结构时展现出极高的效率。本教程聚焦于 ABACUS 与 pyatb（Python Ab-initio Tight-Binding）工具包的协同应用，旨在为读者提供一套从底层理论到实际操作的完整线性光学计算方案。

### 学习路线图
本教程采用循序渐进的结构，确保读者能够平滑掌握计算全流程：
1.  **理论筑基**（第一章）：深入理解 Kubo-Greenwood 公式与介电函数的物理本质，明确 ABACUS 与 pyatb 的分工协作关系。
2.  **基态准备**（第二章）：学习如何配置 ABACUS 进行自洽场（SCF）计算，重点掌握生成哈密顿量、重叠矩阵及位置算符矩阵的关键参数。
3.  **模型构建**（第三章）：利用 pyatb 模块提取矩阵元，构建紧束缚模型，并进行非自洽的光学响应计算。
4.  **数据转化**（第四章）：通过后处理脚本将抽象的介电函数转化为折射率、吸收系数等实验可观测物理量，并攻克初学者最难跨越的单位换算关卡。

### 知识体系定位
在 ABACUS 的知识图谱中，本教程处于“激发态与响应性质”模块。相较于实时演化含时密度泛函理论（rt-TDDFT）等全电子动力学方法，基于 pyatb 的线性响应计算在计算效率与物理图像的清晰度之间取得了优秀的平衡。它不仅是研究半导体、绝缘体光学特性的利器，也是深入理解材料电子跃迁机制的进阶必修课。

### 前置要求
在开始本教程的学习之前，我们建议读者具备以下基础：
*   **量子力学与固体物理基础**：了解能带结构、费米能级及电磁波与物质相互作用的基本概念。
*   **ABACUS 基础操作**：熟悉 `INPUT`、`STRU`、`KPT` 等核心输入文件的配置。
*   **Python 环境基础**：具备基本的 Python 脚本运行与数据绘图能力（如使用 Matplotlib）。

---

# 第一章：理论基础与计算原理

欢迎来到《ABACUS 实战教程》。在开始敲击键盘运行计算之前，我们需要先建立坚实的理论地基。

计算材料的光学性质（如介电函数、吸收光谱），本质上是在回答一个问题：**当光（电磁波）照射到材料上时，材料内部的电子如何响应？**

ABACUS 结合 pyatb（Python Ab-initio Tight-Binding）工具包，提供了一套高效的解决方案：首先利用 ABACUS 基于数值原子轨道（LCAO）进行自洽计算，构建实空间的紧束缚模型；然后利用 pyatb 基于 Kubo-Greenwood 公式计算线性光学响应。

本章将带你深入理解这一流程背后的物理图像与核心参数。

---

## 1.1 线性光学响应与介电函数

在宏观层面，材料的光学性质通常由**复介电函数** $\epsilon(\omega)$ 来描述：

$$
\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)
$$

这两个分量分别对应着不同的物理过程：

*   **虚部 $\epsilon_2(\omega)$（吸收）**：描述了材料对光能量的耗散。在微观上，这对应于电子吸收光子能量，从占据态（价带）跃迁到非占据态（导带）的过程（Interband transition）。这是我们计算的**起点**。
*   **实部 $\epsilon_1(\omega)$（色散）**：描述了光在材料中传播时的极化响应，直接关联折射率。

### Kramers-Kronig 关系
在计算流程中，我们无法直接同时算出实部和虚部。通常的做法是：
1.  利用电子跃迁概率直接计算**虚部** $\epsilon_2(\omega)$。
2.  利用 **Kramers-Kronig (KK) 关系**，通过虚部积分推导出**实部** $\epsilon_1(\omega)$：

$$
\epsilon_1(\omega) = 1 + \frac{2}{\pi} \mathcal{P} \int_0^\infty \frac{\omega' \epsilon_2(\omega')}{\omega'^2 - \omega^2} d\omega'
$$

> **实战提示**：因为实部依赖于对全频段虚部的积分，所以计算 $\epsilon_2(\omega)$ 时，频率范围（$\omega$）需要覆盖足够宽的能量区间，以保证实部在低能区的准确性。

---

## 1.2 基于原子轨道基组的紧束缚模型

ABACUS 的核心优势之一是支持数值原子轨道（LCAO）基组。与平面波基组不同，LCAO 天然适合构建**紧束缚（Tight-Binding）模型**，这使得后续的光学性质计算非常高效。

### 核心矩阵的构建
在 ABACUS 完成自洽场（SCF）计算后，我们需要将电子结构信息“转录”为实空间（Real Space）的矩阵形式，供 pyatb 读取。这涉及三个关键矩阵：

1.  **哈密顿矩阵 $H(R)$**：描述电子在轨道间的跳跃能量。
2.  **重叠矩阵 $S(R)$**：描述原子轨道之间的非正交性（$\langle \phi_i | \phi_j \rangle \neq \delta_{ij}$）。
3.  **位置/偶极矩阵 $r(R)$**：描述电子在轨道间跃迁时的偶极矩，这是计算光吸收强度的核心。

### ABACUS 参数设置 (INPUT)
为了输出上述矩阵，必须在 ABACUS 的 `INPUT` 文件中进行如下设置：

```fortran
INPUT_PARAMETERS
{
    # 基础设置
    basis_type      lcao    # 必须使用 LCAO 基组
    calculation     scf     # 进行自洽计算

    # 关键输出控制 (供 pyatb 使用)
    out_mat_hs2     1       # 输出二阶稀疏格式的 H(R) 和 S(R) 矩阵
    out_mat_r       1       # 输出位置矩阵 r(R)
    out_chg         1       # 建议开启，输出电荷密度以便查错或重启
}
```

**文件流转说明**：
运行 ABACUS 后，会在输出目录（如 `OUT.suffix/`）下生成以下文件：
*   `data-HR-sparse_SPIN0.csr` (哈密顿量)
*   `data-SR-sparse_SPIN0.csr` (重叠矩阵)
*   `data-rR-sparse.csr` (位置矩阵)

> **Data Transfer**: 在运行 pyatb 之前，你必须将这些 `data-*` 文件从 ABACUS 的输出子目录复制到 pyatb 的工作目录中：
> `cp OUT*/data* ./pyatb_directory`

---

## 1.3 Kubo-Greenwood 公式与 pyatb 核心逻辑

pyatb 读取上述矩阵后，利用 **Kubo-Greenwood 公式** 计算光电导率 $\sigma(\omega)$。

### 物理原理
光电导率 $\sigma(\omega)$ 与介电函数的关系为：
$$
\epsilon(\omega) = 1 + \frac{i\sigma(\omega)}{\omega\epsilon_0}
$$

计算的核心在于**速度矩阵元**（Velocity Matrix Elements），它描述了电子在能带 $n$ 和 $m$ 之间跃迁的概率幅。pyatb 通过偶极矩阵 $r(R)$ 和哈密顿量 $H(R)$ 自动处理这一过程。

### 关键输入参数：Fermi Energy (费米能级)
这是新手最容易犯错的地方。Kubo 公式中包含费米-狄拉克分布函数 $f_{n\mathbf{k}}$，它决定了哪些态是占据的（电子源），哪些是空的（跃迁终态）。

**必须手动传递费米能级**：
pyatb 不会自动从 ABACUS 的二进制输出中读取费米能级，你必须：
1.  查看 ABACUS 的日志文件 `OUT.*/running_scf.log`。
2.  搜索关键词 `E_Fermi`。
3.  将该数值精确填入 pyatb 的 `Input` 文件中。

**ABACUS Log 示例**:
```text
E_Fermi            0.4070749075         5.5385382545
```
*(注意：通常第二列是 Rydberg 单位，第三列是 eV 单位，请根据 pyatb Input 中的单位设置选择)*

**pyatb Input 示例**:
```text
INPUT_PARAMETERS
{
    fermi_energy        5.5385382545  # 必须与 SCF 结果一致！
    fermi_energy_unit   eV
    ...
}
```
> **警告**：如果费米能级填错，会导致占据态判断错误（例如把导带当成价带），计算出的光谱将完全错误且无物理意义。

### 关键输入参数：晶格信息 (Lattice)
pyatb 同样需要手动输入晶格常数和矢量。请直接复制 ABACUS `INPUT` 或 `STRU` 文件中的信息。

**特别注意单位**：
*   ABACUS 的 `STRU` 文件中，晶格常数通常以 **Bohr** 为单位。
*   在 pyatb 的 `Input` 中，务必指定 `lattice_constant_unit Bohr` 以保持一致。

---

## 1.4 后处理中的单位陷阱：吸收系数

计算得到介电函数 $\epsilon_1, \epsilon_2$ 后，我们通常需要计算**吸收系数 (Absorption Coefficient, $\alpha$)**。

原始物理公式为：
$$
\alpha(\omega) = \frac{2\omega}{c} \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1}{2}}
$$
或者更直观的形式：$\alpha(\omega) \approx \frac{\omega \epsilon_2}{n c}$ （在弱吸收近似下）。

**单位换算的痛点**：
*   计算输出的 $\omega$ 通常是能量单位 (**eV**)。
*   光速 $c$ 是 SI 单位 (**m/s**)。
*   吸收系数 $\alpha$ 的标准单位通常是 **$cm^{-1}$**。

为了得到正确的 $cm^{-1}$ 量级，必须在 Python 后处理脚本中进行如下转换：

1.  **常数准备**：
    *   $\hbar \approx 4.1357 \times 10^{-15} \text{ eV}\cdot\text{s}$
    *   光速 $c \approx 3 \times 10^{10} \text{ cm/s}$ (注意这里用 cm/s)

2.  **频率转换**：
    将能量 $E$ (eV) 转换为角频率 $\omega$ (rad/s)：
    $$ \omega [\text{s}^{-1}] = \frac{E [\text{eV}]}{\hbar [\text{eV}\cdot\text{s}]} $$

3.  **代码实现逻辑**：
    ```python
    # 伪代码示例
    # energy_ev 是从文件读取的能量轴
    # epsilon_1, epsilon_2 是介电函数实部和虚部
    
    # 1. 计算消光系数 kappa
    kappa = np.sqrt((np.sqrt(epsilon_1**2 + epsilon_2**2) - epsilon_1) / 2)
    
    # 2. 定义常数
    hbar_ev_s = 4.135667696e-15  # Planck constant in eV*s
    c_cm_s = 2.99792458e10       # Speed of light in cm/s
    
    # 3. 计算角频率 omega (s^-1)
    omega = energy_ev / hbar_ev_s
    
    # 4. 计算吸收系数 alpha (cm^-1)
    # 公式: alpha = 2 * omega * kappa / c
    alpha = 2 * omega * kappa / c_cm_s
    ```

如果不进行这一步严谨的单位换算，得到的吸收系数值可能会相差十几个数量级，导致无法与实验数据（如 UV-Vis 光谱）对比。

---

**本章小结**：
我们已经理解了从 DFT 到 LCAO 紧束缚模型，再到介电函数的完整数据流。关键在于：
1.  ABACUS `INPUT` 中开启 `out_mat_*`。
2.  pyatb `Input` 中精确匹配 `fermi_energy` 和晶格参数。
3.  后处理时严谨对待物理常数的单位转换。

下一章，我们将进入实战环节，配置环境并运行第一个二氧化硅（SiO$_2$）的光学性质计算案例。

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

# 第三章：介电函数计算 (pyatb)

在上一章中，我们通过 ABACUS 的自洽场（SCF）计算，成功得到了体系的哈密顿量矩阵（HR）、重叠矩阵（SR）和位置算符偶极矩阵（rR）。这些矩阵是构建紧束缚模型（Tight-Binding Model）的基石。

本章我们将进入 `pyatb` 模块，利用这些矩阵进行非自洽的光学响应计算。本章的核心在于**数据的一致性**——即如何确保 `pyatb` 的输入参数（特别是费米能级和晶格信息）与 ABACUS 的 SCF 计算结果严格对应。

## 3.1 数据准备与环境配置

`pyatb` 是一个独立的后处理模块，通常建议在一个新的工作目录中运行，以保持文件结构清晰。我们需要将 ABACUS 计算生成的稀疏矩阵文件转移到该目录。

### 3.1.1 建立工作目录

假设你的 ABACUS SCF 计算目录为 `silica_PrimaryCell`，且输出结果位于 `OUT.silica` 文件夹中。

```bash
# 在案例根目录下
mkdir pyatb_OpticConductivity
cd pyatb_OpticConductivity
```

### 3.1.2 转移核心矩阵文件

ABACUS 在设置了 `out_mat_hs2 1` 和 `out_mat_r 1` 后，会在输出目录中生成以下关键文件，请务必将其复制到当前目录：

*   `data-HR-sparse_SPIN0.csr`: 哈密顿量矩阵（稀疏格式）
*   `data-SR-sparse_SPIN0.csr`: 重叠矩阵（稀疏格式）
*   `data-rR-sparse.csr`: 偶极矩阵（用于计算跃迁矩阵元）

**执行命令：**
```bash
# 注意：请根据你的实际输出文件夹名称调整路径（如 OUT.silica 或 OUT.suffix）
cp ../OUT.silica/data-HR-sparse_SPIN0.csr .
cp ../OUT.silica/data-SR-sparse_SPIN0.csr .
cp ../OUT.silica/data-rR-sparse.csr .
```

> **注意**：如果你的计算是自旋极化的（`nspin 2`），你可能还会看到 `SPIN1` 的文件，需要根据研究需求一并复制。本教程以非磁性 SiO2 为例，仅涉及 `SPIN0`。

---

## 3.2 编写 pyatb 输入文件

`pyatb` 的输入文件通常命名为 `Input`。这个文件分为三个部分：基本参数、晶格信息、以及光学计算参数。

### 3.2.1 关键参数获取（Critical）

在编写 `Input` 之前，你必须从 ABACUS 的输出日志中提取两个至关重要的物理量。如果这一步出错，后续计算的光谱将完全错误。

**1. 费米能级 (Fermi Energy)**
`pyatb` 不会重新计算费米能级，必须手动指定。请查看 SCF 计算的日志文件 `running_scf.log`。

```bash
# 在 ABACUS 输出目录执行
grep "E_Fermi" running_scf.log | tail -n 1
```
*输出示例：*
```text
E_Fermi        0.4070749075         5.5385382545
```
这里第二列通常是 Rydberg 单位，**第三列是 eV 单位**。我们在 `pyatb` 中通常使用 eV。

**2. 占据带数目 (Occupied Bands)**
同理，确认体系有多少条占据能带：
```bash
grep "occupied bands" running_scf.log
```

### 3.2.2 编写 Input 文件

新建文件 `Input` 并填入以下内容。请务必根据你的实际 SCF 结果修改 `fermi_energy` 和 `LATTICE` 部分。

```text
INPUT_PARAMETERS
{
    nspin           1
    package         ABACUS
    
    # [CRITICAL] 必须与 SCF log 中的最终 E_Fermi (eV) 完全一致
    fermi_energy    5.5385382545
    fermi_energy_unit eV
    
    # 指定矩阵文件路径
    HR_route        data-HR-sparse_SPIN0.csr
    SR_route        data-SR-sparse_SPIN0.csr
    rR_route        data-rR-sparse.csr
    
    # ABACUS 输出矩阵的默认单位
    HR_unit         Ry
    rR_unit         Bohr
}

LATTICE
{
    # [CRITICAL] 必须与 ABACUS INPUT/STRU 文件一致
    # 注意单位：这里使用的是 Bohr
    lattice_constant 1.8897261246257702
    lattice_constant_unit Bohr
    
    lattice_vector
    7.1199998856 0.0 0.0 
    0.0 7.1199998856 0.0 
    0.0 0.0 7.1199998856 
}

OPTICAL_CONDUCTIVITY
{
    # [CRITICAL] 占据带数目
    occ_band        64
    
    # 能量范围 (eV): 从 0 到 30 eV
    omega           0 30
    # 能量步长 (eV)
    domega          0.01
    # 展宽因子 (eV)，通常取 0.05 - 0.1
    eta             0.1
    
    # 光学性质计算通常需要比 SCF 更密的 K 点网格
    grid            20 20 20
}
```

**参数详解：**
*   **`fermi_energy`**: 决定了电子的占据分布（Fermi-Dirac 分布）。如果填错，会导致本来空的带被占据，或者占据带变空，直接导致错误的光学跃迁。
*   **`LATTICE`**: `pyatb` 需要利用晶格矢量进行倒空间积分。建议直接从 ABACUS 的 `OUT.*/STRU_READIN_ADJUST.cif` 或原始 `INPUT` 中复制，**务必确认单位是 Bohr 还是 Angstrom**。
*   **`grid`**: 光学性质计算涉及布里渊区积分，对 K 点密度收敛要求较高。如果光谱出现非物理的剧烈震荡，尝试增加此处的网格密度（如 `30 30 30`）。

---

## 3.3 运行 pyatb 与输出解析

### 3.3.1 运行程序

配置完成后，使用 MPI 并行运行 `pyatb`。计算时间取决于 K 点网格密度和能带数量。

```bash
export OMP_NUM_THREADS=1
mpirun -np 16 pyatb
```

运行过程中，终端可能不会有大量输出。程序运行完毕后，会在当前目录下生成 `Out/Optical_Conductivity` 文件夹。

### 3.3.2 输出文件解析

进入输出目录：
```bash
cd Out/Optical_Conductivity
ls
```

主要关注以下两个文件：
1.  **`dielectric_function_real_part.dat`**: 介电函数实部 $\epsilon_1(\omega)$
2.  **`dielectric_function_imag_part.dat`**: 介电函数虚部 $\epsilon_2(\omega)$

**数据格式：**
文件第一列为光子能量 $\hbar\omega$ (eV)，后续列为介电张量的分量。对于正交晶系，主要关注对角项 `xx`, `yy`, `zz`。

```text
# omega(eV)      xx            xy            xz ...
   0.00000      1.896e+00     ...
   0.01000      1.896e+00     ...
```

---

## 3.4 后处理：计算吸收系数与单位换算

得到介电函数的实部 $\epsilon_1$ 和虚部 $\epsilon_2$ 后，我们可以推导出折射率 $n$、消光系数 $\kappa$ 和吸收系数 $\alpha$。

其中，**吸收系数 $\alpha(\omega)$ 的单位换算**是初学者最容易出错的地方。

### 3.4.1 物理公式

吸收系数的定义为：
$$ \alpha(\omega) = \frac{2\omega}{c} \kappa(\omega) = \frac{2\omega}{c} \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1}{2}} $$

或者写为：
$$ \alpha(\omega) = \frac{\sqrt{2}\omega}{c} \sqrt{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1} $$

### 3.4.2 Python 脚本中的单位陷阱

在数值计算中，我们需要统一单位。通常输出的吸收系数单位为 $\text{cm}^{-1}$。

*   **能量 $\omega$**: 数据文件中的第一列是能量 $E$ (eV)。我们需要将其转换为频率或角频率，或者直接在公式中带入常数。
    *   关系式：$E = \hbar \omega \implies \omega = E / \hbar$
*   **光速 $c$**: 常用单位是 $\text{m/s}$，但为了得到 $\text{cm}^{-1}$，建议使用 $\text{cm/s}$。
    *   $c \approx 2.9979 \times 10^{10} \text{ cm/s}$

**推荐的 Python 处理逻辑：**

```python
import numpy as np

# 1. 读取数据 (假设已加载为 numpy 数组)
# energy_eV: 第一列能量
# eps_real: 实部 (xx 分量)
# eps_imag: 虚部 (xx 分量)

# 2. 物理常数
hbar_eVs = 6.582119569e-16      # Planck constant in eV*s
c_cm_s   = 2.99792458e10        # Speed of light in cm/s

# 3. 将能量(eV)转换为角频率 omega (rad/s)
omega = energy_eV / hbar_eVs

# 4. 计算模长
eps_abs = np.sqrt(eps_real**2 + eps_imag**2)

# 5. 计算吸收系数 alpha (cm^-1)
# 公式: alpha = (omega / c) * sqrt(2) * sqrt( |eps| - eps_real )
# 注意：有些文献公式略有不同，本质是 2*omega*kappa/c
alpha = (omega / c_cm_s) * np.sqrt(2) * np.sqrt(eps_abs - eps_real)

# 6. 绘图 (略)
```

> **专家提示**：在某些简化的脚本中（如本章参考案例），可能会直接使用数值因子 `2.418e14` 将 eV 转换为频率，并使用 $c=3.0\times10^8$ m/s。如果这样做，计算出的 $\alpha$ 单位是 $\text{m}^{-1}$。**为了得到标准的 $\text{cm}^{-1}$，结果必须除以 100**。请务必检查你的数量级，典型的半导体/绝缘体带边吸收系数在 $10^4 \sim 10^5 \text{ cm}^{-1}$ 量级。

# 第四章：后处理与线性光学性质分析

在上一章中，我们已经成功运行了 ABACUS 的 SCF 计算，并利用 `pyatb` 模块获得了介电函数（Dielectric Function）的原始数据。然而，原始的介电函数 $\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)$ 并不直观，无法直接与实验测量的光谱（如吸收光谱、反射光谱）进行对比。

本章将重点介绍如何通过 Python 脚本对数据进行后处理，将介电函数转化为折射率、消光系数、吸收系数等实验可观测的物理量。特别是**吸收系数的单位换算**，这是初学者最容易出错的环节。

## 4.0 前置检查：数据的正确传递

在开始画图之前，我们必须确保 `pyatb` 计算所依据的物理图像是正确的。光学性质计算对电子结构极为敏感，以下三个步骤是**绝对关键**的（Critical），任何一步出错都会导致后续分析毫无意义。

### 1. 数据的转移
ABACUS 的计算结果通常存储在 `OUT.${suffix}/` 目录下，而 `pyatb` 通常在当前目录或独立目录运行。你需要将紧束缚矩阵文件复制过来：

```bash
# 假设你在案例根目录
# 复制 ABACUS 生成的哈密顿量(HR)、重叠矩阵(SR)和位置矩阵(rR)
cp OUT.silica/data-HR-sparse_SPIN0.csr ./pyatb_OpticConductivity/
cp OUT.silica/data-SR-sparse_SPIN0.csr ./pyatb_OpticConductivity/
cp OUT.silica/data-rR-sparse.csr       ./pyatb_OpticConductivity/
```

### 2. 费米能级 (Fermi Energy) 的精确传递
`pyatb` 不会自动读取 ABACUS 的费米能级。如果费米能级填错，会导致能带占据数（Occupation）错误，进而导致带间跃迁计算完全错误（例如将绝缘体算成金属）。

*   **步骤 A**: 查看 ABACUS 的 SCF 日志：
    ```bash
    grep "E_Fermi" OUT.silica/running_scf.log | tail -n 1
    # 输出示例: E_Fermi        0.4070749075         5.5385382545
    ```
    这里第二列通常是 Ry 单位，**第三列是 eV 单位**。

*   **步骤 B**: 将 eV 单位的数值填入 `pyatb` 的 `Input` 文件：
    ```text
    INPUT_PARAMETERS
    {
        ...
        fermi_energy 5.5385382545  <-- 必须精确匹配
        fermi_energy_unit eV
        ...
    }
    ```

### 3. 晶格信息的对齐
`pyatb` 需要显式的晶格常数和矢量。请务必检查 ABACUS `INPUT` 或 `STRU` 文件中的单位（Bohr 或 Angstrom）。本案例使用 **Bohr**：

```text
LATTICE
{
    lattice_constant 1.8897261246257702
    lattice_constant_unit Bohr
    lattice_vector
    7.1199998856 0.0 0.0 
    0.0 7.1199998856 0.0 
    0.0 0.0 7.1199998856 
}
```

---

## 4.1 介电函数的可视化

运行 `pyatb` 后，在输出目录（如 `Out/Optical_Conductivity`）下会生成两个核心文件：
- `dielectric_function_real_part.dat` ($\epsilon_1$)
- `dielectric_function_imag_part.dat` ($\epsilon_2$)

文件格式通常包含能量（eV）和介电张量的各个分量（xx, xy, xz, ..., zz）。对于各向同性体系或立方晶系（如本例中的 $\beta$-Cristobalite），$xx=yy=zz$；对于各向异性体系，则需要分别分析。

以下 Python 脚本用于读取并绘制介电函数的实部和虚部：

```python
import numpy as np
import matplotlib.pyplot as plt

# --- 文件路径配置 ---
path = './Out/Optical_Conductivity/'
file_real = path + 'dielectric_function_real_part.dat'
file_imag = path + 'dielectric_function_imag_part.dat'

# --- 读取数据函数 ---
def read_dielectric(filename):
    with open(filename, 'r') as f:
        # 跳过第一行表头
        lines = f.readlines()[1:]
    
    energy = []
    avg_val = [] # 存储各向同性平均值 (xx+yy+zz)/3
    
    for line in lines:
        data = [float(x) for x in line.split()]
        if len(data) < 10: continue
        
        eng = data[0] # 第一列是能量 (eV)
        xx, yy, zz = data[1], data[5], data[9] # 提取对角分量
        
        energy.append(eng)
        avg_val.append((xx + yy + zz) / 3.0)
        
    return np.array(energy), np.array(avg_val)

# --- 读取数据 ---
energy, eps_real = read_dielectric(file_real)
_, eps_imag = read_dielectric(file_imag)

# 计算复介电函数的模
eps_tot = np.sqrt(eps_real**2 + eps_imag**2)

# --- 绘图 ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(energy, eps_real, label=r'Real Part ($\epsilon_1$)', color='red', linewidth=2)
ax.plot(energy, eps_imag, label=r'Imaginary Part ($\epsilon_2$)', color='green', linewidth=2)
# ax.plot(energy, eps_tot, label='Total', color='black', linestyle='--')

ax.set_xlabel('Energy (eV)', fontsize=14)
ax.set_ylabel('Dielectric Function', fontsize=14)
ax.set_title(r'Dielectric Function of SiO$_2$', fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

**物理分析提示**：
*   **实部 $\epsilon_1(0)$**：零频极限对应材料的静态介电常数。
*   **虚部 $\epsilon_2(\omega)$**：对应光的吸收。第一个峰值的起始位置通常对应材料的光学带隙（Optical Gap）。

---

## 4.2 导出量的物理推导与计算

有了 $\epsilon_1(\omega)$ 和 $\epsilon_2(\omega)$，我们可以利用经典电磁学关系推导其他光学性质。

### 核心公式
1.  **复折射率** $N(\omega) = n(\omega) + i\kappa(\omega)$，其中：
    *   **折射率 (Refractive Index)** $n(\omega)$:
        $$ n(\omega) = \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} + \epsilon_1}{2}} $$
    *   **消光系数 (Extinction Coefficient)** $\kappa(\omega)$:
        $$ \kappa(\omega) = \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1}{2}} $$

2.  **反射率 (Reflectivity)** $R(\omega)$ (正入射):
    $$ R(\omega) = \frac{(n-1)^2 + \kappa^2}{(n+1)^2 + \kappa^2} $$

3.  **能量损失函数 (Energy Loss Function)** $L(\omega)$:
    常用于分析电子能量损失谱 (EELS) 和等离激元 (Plasmon) 峰。
    $$ L(\omega) = \text{Im}\left(\frac{-1}{\epsilon(\omega)}\right) = \frac{\epsilon_2}{\epsilon_1^2 + \epsilon_2^2} $$

### Python 实现
```python
# 基于 4.1 中读取的 eps_real 和 eps_imag 数组

# 1. 折射率 n
n_idx = np.sqrt((np.sqrt(eps_real**2 + eps_imag**2) + eps_real) / 2)

# 2. 消光系数 kappa
kappa = np.sqrt((np.sqrt(eps_real**2 + eps_imag**2) - eps_real) / 2)

# 3. 反射率 R
reflectivity = ((n_idx - 1)**2 + kappa**2) / ((n_idx + 1)**2 + kappa**2)

# 4. 能量损失函数 L
loss_function = eps_imag / (eps_real**2 + eps_imag**2)
```

---

## 4.3 吸收系数的单位换算（难点）

吸收系数 $\alpha(\omega)$ 描述了光强随穿透深度衰减的快慢，是实验中最常测量的量。其标准定义为：
$$ \alpha(\omega) = \frac{2\omega}{c} \kappa(\omega) $$

### 单位陷阱
在计算中，我们直接得到的变量单位如下：
*   $\omega$ (能量轴): 单位是 **eV**。
*   $c$ (光速): 物理常数，单位通常是 $m/s$。
*   目标 $\alpha$: 实验常用单位是 **$cm^{-1}$**。

如果我们直接把 eV 数值乘进去，结果的数量级会完全错误。我们需要进行严格的单位制转换。

### 推导过程
1.  利用关系 $E = \hbar\omega$，将频率 $\omega$ 替换为能量 $E$：
    $$ \omega = \frac{E}{\hbar} $$
    代入公式得：
    $$ \alpha(E) = \frac{2E}{\hbar c} \kappa(E) $$

2.  **常数处理**：我们需要 $1/(\hbar c)$ 的数值，且要保证最终单位兼容 $cm^{-1}$。
    *   光速 $c \approx 2.9979 \times 10^{10} \text{ cm/s}$ (注意这里用 cm/s 以匹配目标单位)。
    *   约化普朗克常数 $\hbar \approx 6.5821 \times 10^{-16} \text{ eV}\cdot\text{s}$。
    *   组合常数项：
        $$ \frac{1}{\hbar c} \approx \frac{1}{(6.5821 \times 10^{-16} \text{ eV}\cdot\text{s}) \times (2.9979 \times 10^{10} \text{ cm/s})} $$
    *   或者利用 $\hbar$ 的倒数关系：$1/\hbar \approx 1.519 \times 10^{15} \text{ s}^{-1}/\text{eV}$。

    在实际代码中，为了保持物理清晰，我们可以显式写出转换因子。

### Python 代码实现
```python
# --- 物理常数 ---
# 普朗克常数 hbar (eV·s)
hbar_eVs = 6.582119569e-16 
# 光速 c (cm/s) -> 转换为 cm/s 是为了直接得到 cm^-1 的结果
c_cms = 2.99792458e10 

# --- 计算吸收系数 ---
# 公式: alpha = (2 * E / hbar) * kappa / c
# 单位分析: [eV] / [eV·s] * [1] / [cm/s] = [1/s] * [s/cm] = [1/cm]
absorption_coeff = (2 * energy / hbar_eVs) * kappa / c_cms

# --- 绘图展示 ---
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# 折射率
axs[0, 0].plot(energy, n_idx, 'b')
axs[0, 0].set_ylabel('Refractive Index n')
axs[0, 0].set_title('Refractive Index')

# 消光系数
axs[0, 1].plot(energy, kappa, 'm')
axs[0, 1].set_ylabel('Extinction Coefficient $\kappa$')
axs[0, 1].set_title('Extinction Coefficient')

# 吸收系数 (使用对数坐标展示数量级)
axs[1, 0].plot(energy, absorption_coeff, 'y')
axs[1, 0].set_ylabel(r'Absorption Coeff $\alpha$ ($cm^{-1}$)')
axs[1, 0].set_yscale('log') # 吸收系数跨度大，推荐用对数坐标
axs[1, 0].set_title('Absorption Coefficient')

# 反射率
axs[1, 1].plot(energy, reflectivity, 'r')
axs[1, 1].set_ylabel('Reflectivity R')
axs[1, 1].set_title('Reflectivity')

for ax in axs.flat:
    ax.set_xlabel('Energy (eV)')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

通过上述代码，你就可以得到与实验文献中单位一致（$cm^{-1}$）的吸收光谱，通常在 $10^4 \sim 10^6$ 量级。

---

## 附录：常见问题与进阶建议

### 1. K 点收敛性 (K-point Convergence)
**问题**：为什么计算光学性质时，K 点密度要求比计算总能量（SCF）时高得多？
**解答**：介电函数的虚部涉及布里渊区内的积分 $\sum_{k} \delta(E_{ck} - E_{vk} - \hbar\omega)$。为了精确捕捉能带之间的共振跃迁（即 $\delta$ 函数的峰），需要非常密集的 K 点网格来采样能级差。如果 K 点不足，光谱会出现非物理的锯齿状震荡。

### 2. 空带数量 (Number of Empty Bands)
**问题**：高能区的光谱为什么看起来不准？
**解答**：`pyatb` 的 `Input` 中 `occ_band` 参数（或 ABACUS 的 `nbands`）决定了参与跃迁的最高能级。根据求和公式，高能光子激发的电子会跃迁到更高的空带。如果你计算的空带数量不足，高能区的跃迁通道被截断，导致高能端介电函数趋于零。

### 3. 费米能级错误的影响
**现象**：绝缘体算出了巨大的零频光导率（Drude peak）。
**原因**：如果在 `Input` 中 `fermi_energy` 填错（例如填得太高进入了导带），程序会误认为体系是金属，从而产生带内跃迁（Intraband transition）贡献。对于半导体/绝缘体，务必确保费米能级处于带隙之中。

---

## 附录：进阶学习指南

恭喜你完成了本教程的学习！掌握线性光学性质的计算只是探索材料物理世界的起点。为了帮助你进一步深造，我们准备了以下建议：

### 1. 进阶探索主题
*   **非线性光学性质**：本教程侧重于线性响应。若要研究二倍频（SHG）或多光子吸收，需关注更高级的微扰理论方法。
*   **实时演化 TDDFT (rt-TDDFT)**：对于强场下的超快动力学过程，ABACUS 提供的 rt-TDDFT 功能（参考资料5）是更合适的选择。它可以模拟材料在任意电场下的实时电子演化。
*   **激子效应 (GW/BSE)**：pyatb 目前主要处理独立粒子近似下的跃迁。若要考虑电子-空穴相互作用（激子），建议关注更高阶的准粒子修正方法。

### 2. 通用调试建议 (Troubleshooting)
在计算过程中，如果遇到结果异常，请检查以下几点：
*   **收敛性检查**：光学性质对 K 点密度非常敏感。如果光谱曲线抖动剧烈，请尝试增加 K 点采样。
*   **矩阵一致性**：确保 pyatb 读取的费米能级与 ABACUS SCF 计算输出的结果严格一致，否则会导致占据态判断错误。
*   **单位陷阱**：在计算吸收系数时，务必核实能量单位（eV）与长度单位（cm）的转换关系。请参考第四章提供的物理常数换算公式。

### 3. 官方资源与社区
*   **ABACUS 官方文档**：获取最新的输入参数说明与功能更新。
*   **Bohrium 算力平台**：本教程推荐在 Bohrium 平台上运行示例，利用其预装的 ABACUS 与 pyatb 环境快速复现计算结果。
*   **开源社区交流**：通过 GitHub 或开发者社区反馈使用中的 Bug，参与到国产科学软件的建设中。

**科学研究的魅力在于不断探索未知。希望本教程能成为你科研路上的坚实垫脚石，祝你在材料计算领域取得丰硕成果！**
