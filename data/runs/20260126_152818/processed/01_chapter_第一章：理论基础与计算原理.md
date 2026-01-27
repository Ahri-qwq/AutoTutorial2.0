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