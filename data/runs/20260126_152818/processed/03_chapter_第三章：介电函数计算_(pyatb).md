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