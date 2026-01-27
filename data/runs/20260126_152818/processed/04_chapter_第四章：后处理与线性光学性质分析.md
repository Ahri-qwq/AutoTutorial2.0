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