# 第四章：后处理与光谱分析

在上一章中，我们通过 `pyatb` 完成了介电函数的计算，得到了 `dielectric_function_real_part.dat` 和 `dielectric_function_imag_part.dat` 等原始数据文件。然而，这些仅仅是“中间数据”。

作为材料科学家，我们最终关心的是**实验可观测**的物理量，如**吸收系数（Absorption Coefficient）**、**折射率（Refractive Index）**和**反射率（Reflectivity）**。此外，原始输出通常基于原子单位或 eV，而实验光谱常使用 $nm$ 或 $cm^{-1}$。

本章将带你完成从“计算数据”到“实验光谱”的最后一步跨越。我们将使用 Python 进行张量处理、物理公式推导以及单位换算。

---

## 4.1 数据读取与张量平均

ABACUS 和 pyatb 计算得到的介电函数是一个 $3 \times 3$ 的张量（Tensor），包含 $\epsilon_{xx}, \epsilon_{xy}, \dots, \epsilon_{zz}$ 等分量。

### 4.1.1 为什么需要张量平均？
- **单晶/各向异性材料**：需要分别分析 $xx$（沿 x 轴极化）或 $zz$（沿 z 轴极化）分量。
- **多晶/无序/粉末样品**：实验测量通常是各向同性的平均效果。此时我们需要对介电张量的**对角元**取算术平均：
  $$ \epsilon_{\text{iso}}(\omega) = \frac{\epsilon_{xx}(\omega) + \epsilon_{yy}(\omega) + \epsilon_{zz}(\omega)}{3} $$

### 4.1.2 Python 读取脚本
我们将使用 `numpy` 读取数据。请确保你的工作目录下包含 `pyatb` 输出的 `.dat` 文件。

```python
import numpy as np

# 文件路径（请根据实际情况修改）
real_file = 'Out/Optical_Conductivity/dielectric_function_real_part.dat'
imag_file = 'Out/Optical_Conductivity/dielectric_function_imag_part.dat'

def read_dielectric_data(filename):
    """
    读取 pyatb 输出文件，返回能量轴和各向同性平均后的介电函数
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 跳过表头，读取数据
    data = []
    energy = []
    
    for line in lines[1:]: # Skip header
        parts = line.split()
        if len(parts) < 10: continue
        
        # 解析数据列
        # 格式: omega(eV) xx xy xz yx yy yz zx zy zz
        w = float(parts[0])      # 能量 (eV)
        xx = float(parts[1])
        yy = float(parts[5])
        zz = float(parts[9])
        
        # 计算各向同性平均
        avg_val = (xx + yy + zz) / 3.0
        
        energy.append(w)
        data.append(avg_val)
        
    return np.array(energy), np.array(data)

# 读取实部 (epsilon_1) 和 虚部 (epsilon_2)
energy, eps_1 = read_dielectric_data(real_file)
_, eps_2 = read_dielectric_data(imag_file)

print(f"数据读取完成。能量范围: {energy[0]} - {energy[-1]} eV")
```

---

## 4.2 吸收系数计算与单位换算 (核心难点)

这是本章最容易出错的部分。直接将计算出的介电函数代入公式得到的数值，往往与实验值相差十几个数量级，原因就在于**单位制**。

### 4.2.1 物理公式推导
吸收系数 $\alpha(\omega)$ 定义为光强在介质中传播单位距离后的衰减率。它与消光系数 $\kappa(\omega)$ 的关系为：
$$ \alpha(\omega) = \frac{2\omega}{c} \kappa(\omega) $$
其中消光系数 $\kappa$ 由介电函数实部 $\epsilon_1$ 和虚部 $\epsilon_2$ 决定：
$$ \kappa(\omega) = \sqrt{ \frac{ \sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1 }{2} } $$

合并后，我们得到计算用的核心公式：
$$ \alpha(\omega) = \frac{\sqrt{2}\omega}{c} \sqrt{ \sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1 } $$

### 4.2.2 单位换算逻辑
ABACUS/pyatb 输出的能量轴 $E$ 单位是 **eV**，而实验常用的吸收系数单位是 **$\text{cm}^{-1}$**。

1.  **频率 $\omega$ 的转换**:
    我们有 $E = \hbar\omega$，所以 $\omega = E / \hbar$。
    *   普朗克常数 $\hbar \approx 4.135667 \times 10^{-15} \text{ eV}\cdot\text{s}$
    *   所以 $\omega [\text{s}^{-1}] = E [\text{eV}] / (4.1357 \times 10^{-15})$

2.  **光速 $c$ 的转换**:
    为了得到 $\text{cm}^{-1}$，光速必须使用 $\text{cm/s}$ 单位。
    *   $c \approx 2.9979 \times 10^{8} \text{ m/s} = 2.9979 \times 10^{10} \text{ cm/s}$

3.  **最终系数**:
    将常数提取出来：
    $$ \text{Factor} = \frac{1}{\hbar \cdot c} \approx \frac{1}{4.1357 \times 10^{-15} \times 2.9979 \times 10^{10}} \approx 8065.5 $$
    这意味着：$1 \text{ eV}$ 的能量对应的波数约为 $8065.5 \text{ cm}^{-1}$。

### 4.2.3 计算脚本实现

```python
def calculate_optical_properties(energy_ev, e1, e2):
    """
    计算线性光学性质
    输入:
        energy_ev: 光子能量 (eV)
        e1: 介电函数实部
        e2: 介电函数虚部
    输出:
        n: 折射率
        k: 消光系数
        alpha: 吸收系数 (cm^-1)
        L: 能量损失函数
    """
    # 1. 中间变量模长
    # epsilon_complex = e1 + 1j * e2
    # modulus = np.abs(epsilon_complex) = sqrt(e1^2 + e2^2)
    modulus = np.sqrt(e1**2 + e2**2)
    
    # 2. 折射率 n
    n = np.sqrt( (modulus + e1) / 2.0 )
    
    # 3. 消光系数 k (kappa)
    k = np.sqrt( (modulus - e1) / 2.0 )
    
    # 4. 吸收系数 alpha (cm^-1)
    # 常数定义
    HBAR_eVs = 4.135667696e-15  # Planck constant in eV*s
    C_cm_s   = 2.99792458e10    # Light speed in cm/s
    
    # 将能量 eV 转换为角频率 omega (rad/s)
    omega = energy_ev / HBAR_eVs
    
    # 公式: alpha = 2 * omega * k / c
    alpha = 2.0 * omega * k / C_cm_s
    
    # 5. 能量损失函数 L (EELS)
    L = e2 / (e1**2 + e2**2)
    
    return n, k, alpha, L

# 执行计算
n, k, alpha, L = calculate_optical_properties(energy, eps_1, eps_2)
```

---

## 4.3 多物理量可视化与分析

最后，我们将计算结果绘制成图。这里我们重点关注**吸收系数**和**能量损失函数**，因为它们最常用于与实验光谱（UV-Vis 吸收谱和 EELS 谱）进行对比。

```python
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Calculated Optical Properties of SiO$_2$', fontsize=16)

# 1. 折射率 (Refractive Index)
axs[0, 0].plot(energy, n, 'b-', linewidth=2)
axs[0, 0].set_ylabel('Refractive Index $n$', fontsize=12)
axs[0, 0].set_xlabel('Energy (eV)', fontsize=12)
axs[0, 0].grid(True, alpha=0.3)

# 2. 消光系数 (Extinction Coefficient)
axs[0, 1].plot(energy, k, 'm-', linewidth=2)
axs[0, 1].set_ylabel('Extinction Coefficient $\kappa$', fontsize=12)
axs[0, 1].set_xlabel('Energy (eV)', fontsize=12)
axs[0, 1].grid(True, alpha=0.3)

# 3. 吸收系数 (Absorption Coefficient)
# 注意：通常使用对数坐标展示，因为吸收边变化剧烈
axs[1, 0].plot(energy, alpha, 'r-', linewidth=2)
axs[1, 0].set_ylabel(r'Absorption Coeff. $\alpha$ (cm$^{-1}$)', fontsize=12)
axs[1, 0].set_xlabel('Energy (eV)', fontsize=12)
axs[1, 0].set_yscale('log') # 对数坐标
axs[1, 0].set_ylim(1e2, 1e6) # 根据实际情况调整范围
axs[1, 0].grid(True, alpha=0.3)

# 4. 能量损失函数 (Energy Loss Function)
axs[1, 1].plot(energy, L, 'g-', linewidth=2)
axs[1, 1].set_ylabel('Energy Loss Function $L(\omega)$', fontsize=12)
axs[1, 1].set_xlabel('Energy (eV)', fontsize=12)
axs[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

### 结果分析指南
- **吸收边 (Absorption Edge)**: 在吸收系数图中，数值急剧上升的能量点通常对应材料的光学带隙（Optical Bandgap）。
- **峰位对应**: 能量损失函数中的峰通常对应等离激元激发（Plasmon resonance），而介电函数虚部 $\epsilon_2$ 的峰则对应带间跃迁（Interband transition）。

---

## 附录：常见问题与进阶建议 (Troubleshooting)

在实际操作中，90% 的错误都发生在 `pyatb` 的输入文件配置上。以下是三个致命陷阱：

### 1. 费米能级陷阱 (The Fermi Energy Trap)
**现象**: 计算出的光谱形状极其怪异，或者完全没有吸收峰。
**原因**: `pyatb` 目前版本**不会**自动从 `.csr` 或 `.orb` 文件中读取费米能级。如果 `Input` 文件中的 `fermi_energy` 填错（例如填了 0 或默认值），程序会错误判断占据态和空态，导致跃迁计算完全错误。
**解决方案**:
必须手动查看 ABACUS 的 SCF 输出日志 `running_scf.log` 或 `OUT.*/running_scf.log`。

```bash
# 在终端运行
grep "E_Fermi" OUT.silica/running_scf.log | tail -n 1
```
输出示例：
```text
E_Fermi                0.4070749075         5.5385382545
```
请务必取**最后一列**（单位通常为 eV）的数值，并精确填入 `pyatb` 的 `Input` 文件中：
```text
fermi_energy 5.5385382545
fermi_energy_unit eV
```

### 2. 占据带数量 (Occupied Bands)
**现象**: 程序报错或计算结果为零。
**原因**: `occ_band` 参数必须准确对应体系的价带数量。
**解决方案**:
同样查看日志：
```bash
grep "occupied bands" OUT.silica/running_scf.log
```
输出示例：
```text
occupied bands = 64
```
在 `Input` 中填入：
```text
occ_band 64
```

### 3. 晶格信息缺失风险
**风险**: `pyatb` 并不总是能从稀疏矩阵文件中重建晶格信息。
**建议**: 永远在 `pyatb` 的 `Input` 文件中显式包含 `LATTICE` 块，内容直接复制自 ABACUS 的 `STRU` 文件（注意单位转换，ABACUS `STRU` 默认是 Bohr，pyatb 也支持 Bohr，但需保持一致）。

### 4. K 点收敛性
**原理**: 光学性质计算涉及布里渊区内的积分。相比于电荷密度（SCF），光学矩阵元随 k 点变化更为剧烈。
**建议**: 光学性质计算所需的 K 点密度通常需要比 SCF 计算密 **2-3 倍**。如果光谱中有大量非物理的毛刺或尖峰，首先尝试增加 `OPTICAL_CONDUCTIVITY` 块中的 `grid` 参数（例如从 `20 20 20` 增加到 `40 40 40`）。