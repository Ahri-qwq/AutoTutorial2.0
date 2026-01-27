# 第三章：介电函数计算 (pyatb)

在上一章中，我们通过 ABACUS 完成了自洽计算（SCF），并成功输出了构建紧束缚模型所需的稀疏矩阵文件（`data-HR-sparse_*.csr`, `data-SR-sparse_*.csr`, `data-rR-sparse.csr`）。

本章将进入**两步走**流程的第二步：使用 **pyatb** 读取这些矩阵，计算介电函数的虚部（$\epsilon_2$），并通过 Kramers-Kronig 变换得到实部（$\epsilon_1$）。这是连接微观电子结构与宏观光学性质的关键桥梁。

> **⚠️ 核心提示**：pyatb 是一个独立的后处理程序。它不进行从头算（ab initio）求解，而是基于 ABACUS 产生的哈密顿量和偶极矩阵进行线性响应计算。因此，**参数的一致性**（特别是费米能级和晶格常数）至关重要。

---

## 3.1 目录结构与数据迁移

首先，我们需要建立一个独立的工作目录来进行光学性质计算，以保持文件条理清晰。

### 3.1.1 建立工作目录
假设你的 ABACUS 计算目录为 `silica_PrimaryCell`，且输出文件位于 `OUT.silica` 中。我们在同级目录下创建一个名为 `pyatb_OpticConductivity` 的文件夹。

```bash
# 在 silica_PrimaryCell 目录下
mkdir pyatb_OpticConductivity
cd pyatb_OpticConductivity
```

### 3.1.2 链接或复制关键矩阵文件
pyatb 运行必须依赖 ABACUS 生成的三个核心矩阵文件。请执行以下命令将它们复制到当前目录：

```bash
# 复制哈密顿量矩阵 (HR)
cp ../OUT.silica/data-HR-sparse_SPIN0.csr .

# 复制重叠矩阵 (SR)
cp ../OUT.silica/data-SR-sparse_SPIN0.csr .

# 复制位置算符偶极矩阵 (rR)
cp ../OUT.silica/data-rR-sparse.csr .
```

> **风险提示**：如果你的 `OUT.silica` 目录下没有这些文件，说明你在上一章的 `INPUT` 文件中未正确设置 `out_mat_hs2` 和 `out_mat_r` 为 `1`。请返回上一章重新运行 ABACUS 计算。

---

## 3.2 pyatb 输入文件编写

这是本章最容易出错的环节。我们需要编写 pyatb 的配置文件 `Input`（注意首字母大写）。该文件包含三个主要部分：参数设置、晶格信息和功能模块。

### 3.2.1 关键参数提取 (Pre-flight Check)
在编写 `Input` 之前，必须从 ABACUS 的日志文件 `running_scf.log` 中提取两个关键数值：**费米能级 (E_Fermi)** 和 **占据带数 (occupied bands)**。

请在终端运行以下命令进行查询：

```bash
# 查询占据带数
grep 'occupied bands' ../OUT.silica/running_scf.log
# 输出示例: occupied bands = 64

# 查询费米能级 (取最后一步迭代的值)
grep 'E_Fermi' ../OUT.silica/running_scf.log | tail -n 1
# 输出示例: E_Fermi                0.4070749075         5.5385382545
```
*注意：ABACUS 输出的费米能级通常有两列，第二列单位为 eV。请记录下 **5.5385382545** (eV) 和 **64** (bands)。*

### 3.2.2 编写 Input 文件
新建名为 `Input` 的文件，填入以下内容。请务必根据上一步提取的数值修改 `fermi_energy` 和 `occ_band`。

```plaintext
INPUT_PARAMETERS
{
    # 基础参数
    nspin           1
    package         ABACUS
    
    # --- 关键参数：必须与 SCF 结果严格一致 ---
    fermi_energy    5.5385382545   # 填入 grep 得到的 eV 数值
    fermi_energy_unit eV
    # -------------------------------------

    # 矩阵文件路径
    HR_route        data-HR-sparse_SPIN0.csr
    SR_route        data-SR-sparse_SPIN0.csr
    rR_route        data-rR-sparse.csr
    
    # 矩阵单位 (ABACUS 默认输出单位)
    HR_unit         Ry
    rR_unit         Bohr
}

LATTICE
{
    # ⚠️ 风险提示：pyatb 不会自动读取 STRU，必须手动填写
    # 这里的数据来自 ABACUS 的 STRU 文件或 running_scf.log
    lattice_constant 1.8897261246257702
    lattice_constant_unit Bohr
    lattice_vector
    7.1199998856 0.0 0.0 
    0.0 7.1199998856 0.0 
    0.0 0.0 7.1199998856 
}

OPTICAL_CONDUCTIVITY
{
    # --- 关键参数：必须与 SCF 结果严格一致 ---
    occ_band        64             # 填入 grep 得到的 occupied bands
    # -------------------------------------

    # 光谱能量范围 (eV)
    omega           0 30           # 计算 0 到 30 eV 的光谱
    domega          0.01           # 能量步长
    
    # 展宽因子 (Broadening)
    eta             0.1            # 洛伦兹展宽，模拟有限寿命效应
    
    # K 点网格
    # 建议比 SCF 计算的网格更密，以获得平滑的光谱
    grid            20 20 20       
}
```

**参数详解**：
*   `eta (0.1)`: 这是一个唯象参数（phenomenological parameter），用于避免分母为零并发散，同时也模拟了电子的有限寿命。值越小，光谱峰越尖锐；值越大，光谱越平滑但特征可能模糊。
*   `grid`: 线性光学性质计算通常需要比 SCF 更密的 K 点网格才能收敛。

---

## 3.3 运行计算与原始数据解析

### 3.3.1 提交任务
配置完成后，使用 MPI 并行运行 pyatb。

```bash
# 推荐使用与 ABACUS 相同的核数，例如 16 核
export OMP_NUM_THREADS=1 
mpirun -np 16 pyatb
```

计算过程通常没有屏幕输出，请查看生成的 `Out/` 目录。

### 3.3.2 监控运行
可以通过查看日志文件确认进度：
```bash
tail -f Out/running.log
```
当看到 `Finish Time` 字样时，表示计算结束。

### 3.3.3 输出文件解析
计算完成后，进入 `Out/Optical_Conductivity` 目录，你会看到以下核心文件：

1.  **`dielectric_function_imag_part.dat`**: 介电函数虚部 $\epsilon_2(\omega)$。这是 pyatb 直接计算的物理量，反映了材料的**光吸收**特性。
2.  **`dielectric_function_real_part.dat`**: 介电函数实部 $\epsilon_1(\omega)$。这是通过 Kramers-Kronig 关系从虚部推导出来的，反映了材料的**折射**特性。

文件格式如下（以虚部为例）：
```text
# omega(eV)      xx           xy           xz           ...
   0.00000   0.000000e+00 0.000000e+00 ...
   0.01000   9.564624e-06 -6.324391e-15 ...
```
*   **第 1 列**: 光子能量 $\hbar\omega$ (eV)。
*   **第 2-10 列**: 介电张量的各个分量 ($\epsilon_{xx}, \epsilon_{xy}, \dots$)。对于各向同性材料或立方晶系（如本例中的 $\beta$-Cristobalite），通常关注 $xx, yy, zz$ 分量的平均值。

---

## 3.4 后处理：计算吸收系数与折射率

得到 $\epsilon_1(\omega)$ 和 $\epsilon_2(\omega)$ 后，我们可以推导出所有其他线性光学性质。

### 3.4.1 物理原理与单位换算
实验上常用的**吸收系数 (Absorption Coefficient, $\alpha$)** 单位是 $\text{cm}^{-1}$。但计算输出是原子单位或 eV。我们需要进行严格的单位换算。

公式如下：
$$ \alpha(\omega) = \frac{2\omega}{c} \kappa(\omega) = \frac{\omega}{c} \sqrt{2\left(\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1\right)} $$

**单位换算逻辑**：
1.  **能量转频率**: $\omega [\text{s}^{-1}] = E [\text{eV}] / \hbar [\text{eV}\cdot\text{s}]$
    *   $\hbar \approx 4.1357 \times 10^{-15} \text{ eV}\cdot\text{s}$
    *   换算因子 $C_{\omega} \approx 2.418 \times 10^{14} \text{ s}^{-1}/\text{eV}$
2.  **光速**: 为匹配 $\text{cm}^{-1}$，光速取 $c \approx 3.0 \times 10^{10} \text{ cm/s}$。

### 3.4.2 Python 后处理脚本
创建一个 Python 脚本 `plot_optics.py`，完成数据读取、物理量计算和绘图。

```python
import numpy as np
import matplotlib.pyplot as plt

# 1. 读取数据函数
def read_dielectric(filename):
    data = np.loadtxt(filename, skiprows=1)
    energy = data[:, 0]
    # 取 xx, yy, zz 的平均值作为各向同性近似
    # xx在第1列(索引1), yy在第5列(索引5), zz在第9列(索引9)
    # 注意：文件列索引从0开始，0是energy
    val_xx = data[:, 1]
    val_yy = data[:, 5]
    val_zz = data[:, 9]
    avg_val = (val_xx + val_yy + val_zz) / 3.0
    return energy, avg_val

# 2. 加载实部和虚部
file_imag = 'Out/Optical_Conductivity/dielectric_function_imag_part.dat'
file_real = 'Out/Optical_Conductivity/dielectric_function_real_part.dat'

energy, epsilon2 = read_dielectric(file_imag)
_, epsilon1 = read_dielectric(file_real)

# 3. 计算光学性质
# 模长 |epsilon|
epsilon_mod = np.sqrt(epsilon1**2 + epsilon2**2)

# 折射率 n (Refractive Index)
n = np.sqrt((epsilon_mod + epsilon1) / 2.0)

# 消光系数 kappa (Extinction Coefficient)
kappa = np.sqrt((epsilon_mod - epsilon1) / 2.0)

# 吸收系数 alpha (Absorption Coefficient) 单位: cm^-1
# 转换常数推导: 
# omega(s^-1) = energy(eV) * 2.418e14
# c(cm/s) = 3.0e10
# factor = 2 * 2.418e14 / 3.0e10 ≈ 1.612e4
const_factor = 1.612e4 
alpha = const_factor * energy * kappa

# 4. 绘图
fig, ax = plt.subplots(2, 2, figsize=(12, 10))

# 介电函数
ax[0, 0].plot(energy, epsilon1, label=r'$\epsilon_1$ (Real)', color='blue')
ax[0, 0].plot(energy, epsilon2, label=r'$\epsilon_2$ (Imag)', color='orange')
ax[0, 0].set_title("Dielectric Function")
ax[0, 0].legend()

# 折射率
ax[0, 1].plot(energy, n, color='green')
ax[0, 1].set_title("Refractive Index (n)")

# 消光系数
ax[1, 0].plot(energy, kappa, color='purple')
ax[1, 0].set_title(r"Extinction Coefficient ($\kappa$)")

# 吸收系数
ax[1, 1].plot(energy, alpha, color='red')
ax[1, 1].set_title(r"Absorption Coefficient ($\alpha$)")
ax[1, 1].set_ylabel(r'$\alpha$ (cm$^{-1}$)')
ax[1, 1].set_yscale('log') # 吸收系数通常跨度很大，建议用对数坐标

for a in ax.flat:
    a.set_xlabel("Energy (eV)")
    a.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("optical_properties.png", dpi=300)
print("绘图完成：optical_properties.png")
```

运行此脚本后，你将获得一张包含介电函数、折射率、消光系数和吸收系数的综合图表。至此，你已经完成了从 DFT 计算到宏观光学性质预测的完整流程。