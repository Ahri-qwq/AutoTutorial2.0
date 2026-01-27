# 第四章：实战案例 B——体心立方铁的自旋 ELF 分析

在上一章中，我们通过水分子了解了电子局域函数（ELF）的基本计算流程。然而，在实际的材料研究中，磁性体系占据了半壁江山。对于磁性材料，电子的自旋极化效应会导致不同自旋通道的电子表现出不同的局域化特征。

本章我们将以经典的磁性材料——**体心立方铁（BCC Fe）**为例，深入讲解如何在 ABACUS 中开启自旋极化 ELF 计算，并重点解析如何分别处理和分析自旋向上（Spin-up）与自旋向下（Spin-down）的电子局域性。

---

## 4.1 自旋极化计算设置

对于磁性体系，ELF 的计算核心在于正确设置自旋参数（`nspin`）并理解输出文件的变化。

### 4.1.1 输入文件准备

我们需要准备标准的 DFT 计算文件（`INPUT`, `STRU`, `KPT`, 赝势文件）。以下是针对 BCC Fe 的关键设置：

**1. 结构文件 (STRU)**
体心立方铁的晶格常数约为 2.86 Å。

```text
ATOMIC_SPECIES
Fe 55.845 Fe.UPF  // 需准备对应的赝势文件

LATTICE_CONSTANT
1.889726125  // 1.8897... Bohr = 1.0 Angstrom

LATTICE_VECTORS
1.43 0.00 0.00   // 2.86 Å 的一半，对应 BCC 原胞矢量
0.00 1.43 0.00
0.00 0.00 1.43

ATOMIC_POSITIONS
Direct
Fe
0.0 0.0 0.0 1.0 1.0 1.0 // 初始磁矩设为 1.0
```

**2. 输入参数 (INPUT)**
这是本章的重点。我们需要开启自旋极化（`nspin 2`）并激活 ELF 输出（`out_elf 1`）。

```text
INPUT_PARAMETERS
# 基础计算参数
calculation     scf
basis_type      pw          # 本例演示平面波基组，LCAO 同样适用
ecutwfc         60          # 铁通常需要较高的截断能
ks_solver       dav
smearing_method mp          # 金属体系推荐 Methfessel-Paxton
smearing_sigma  0.02

# 磁性设置 (关键)
nspin           2           # 开启自旋极化计算

# ELF 输出设置 (关键)
out_elf         1           # 1: 输出 ELF, 0: 不输出
```

> **参数详解：**
> *   **`nspin 2`**: 告诉 ABACUS 这是一个自旋极化计算。此时，电子密度 $\rho$ 被分解为 $\rho_\alpha$ (up) 和 $\rho_\beta$ (down)。
> *   **`out_elf`**: 该参数接受两个整数。第一个控制开关（1 为开启），第二个控制输出精度（默认为 3，即保留 3 位有效数字）。通常设置 `out_elf 1` 即可。

### 4.1.2 运行与输出
运行计算后，检查 `OUT` 文件夹。与非磁性体系不同，你会发现目录下生成了**三个** cube 文件：

1.  **`ELF.cube`**: 总电子局域函数（Total ELF）。
2.  **`ELF_SPIN1.cube`**: 自旋向上电子的 ELF（对应 Majority Spin）。
3.  **`ELF_SPIN2.cube`**: 自旋向下电子的 ELF（对应 Minority Spin）。

---

## 4.2 自旋通道的分辨与分析

在磁性材料中，总 ELF 往往掩盖了磁性电子的轨道特征。通过分别分析 `ELF_SPIN1` 和 `ELF_SPIN2`，我们可以更清晰地观察到未配对电子的空间分布。

### 4.2.1 理论背景：自旋极化下的 ELF 定义

为了正确理解这三个文件，我们需要回顾 ELF 在自旋极化下的定义差异。

**1. 自旋非极化（nspin=1）**
$$ELF = \frac{1}{1 + \chi^2}, \quad \chi = \frac{\tau_{KS} - \tau_{vW}}{\tau_{TF}}$$
其中 $\tau_{TF}$ 是 Thomas-Fermi 动能密度，与总密度 $\rho^{5/3}$ 成正比。

**2. 自旋极化（nspin=2）**
对于特定自旋通道 $\sigma$（$\alpha$ 或 $\beta$），ELF 定义为：
$$ELF_\sigma = \frac{1}{1 + \chi_\sigma^2}, \quad \chi_\sigma = \frac{\tau_{KS\sigma} - \tau_{vW\sigma}}{\tau_{TF\sigma}}$$

**关键差异（Critical）**：
此处需注意动能泛函的**自旋标度率（Spin Scaling）**。对于自旋极化体系，单自旋通道的 Thomas-Fermi 动能密度定义为：
$$\tau_{TF\sigma} = \frac{3}{10}(6\pi^2)^{2/3} \rho_\sigma^{5/3}$$
这与非极化公式中的系数不同（非极化为 $(3\pi^2)^{2/3}$）。ABACUS 内部已严格按照此公式处理，因此输出的 `ELF_SPIN1/2` 是物理意义明确的单自旋通道局域性度量。

**3. 总 ELF**
总 ELF 并非简单的两个自旋 ELF 相加，而是基于总动能密度构建：
$$\chi_{tot} = \frac{(\tau_{KS\alpha} + \tau_{KS\beta}) - (\tau_{vW\alpha} + \tau_{vW\beta})}{\tau_{TF\alpha} + \tau_{TF\beta}}$$
对应输出文件 `ELF.cube`。

### 4.2.2 结果分析示例

使用 VESTA 打开生成的 cube 文件，我们可以观察到以下特征：

*   **ELF.cube (Total)**: 展示了铁原子间的金属键特征以及原子核周围的壳层结构。在 BCC 结构中，电子云呈现出明显的各向异性。
*   **ELF_SPIN1.cube (Majority)**: 对于铁（d 轨道未满），多数自旋电子通常占据更多的 d 轨道态。你会发现其 ELF 分布更加饱满，反映了 d 电子的局域化特征。
*   **ELF_SPIN2.cube (Minority)**: 少数自旋电子数量较少，其 ELF 空间分布通常比 SPIN1 收缩或表现出不同的轨道形状（如 $t_{2g}$ 或 $e_g$ 的差异）。

**实战建议**：
在 VESTA 中，尝试将 `ELF_SPIN1` 和 `ELF_SPIN2` 设置为不同的颜色（例如红色和蓝色），并叠加显示，可以直观地看到磁矩主要来源于哪些空间区域的电子局域化差异。

---

## 附录：常见问题与进阶建议

### 1. LCAO 基组下的稳定性修正（重要）

如果你使用原子轨道基组（`basis_type lcao`）计算 ELF，可能会在远离原子的真空区域观察到 ELF 不为 0 的异常现象（例如出现数值为 0.5 左右的“鬼影”）。

*   **原因**: ELF 的核心变量 $\chi$ 是一个分式。在远离原子核的区域，分子（动能密度差）和分母（Thomas-Fermi 动能）理论上都趋于 0。然而，由于 LCAO 基组的高斯尾部衰减特性，分母 $\tau_{TF}$ 可能比分子趋于 0 的速度更快，导致 $\chi$ 数值发散或趋于有限值，进而使 ELF 出现非物理的数值。
*   **ABACUS 的解决方案**: 为了消除这种数值不稳定性，ABACUS 开发者在 $\chi$ 的分子上增加了一个极小的修正值 $\epsilon$（目前默认取值为 $10^{-5}$）。
    $$ \chi_{modified} = \frac{\tau_{KS} - \tau_{vW} + \epsilon}{\tau_{TF}} $$
    这确保了在真空区（$\tau_{TF} \to 0$），分子始终略大于分母的衰减速度，使得 $\chi \to \infty$，从而强制 ELF 正确地趋于 0。
*   **用户注意**: 这一修正对原子附近的化学成键区域影响极小，但在分析极低密度区域（如范德华吸附的远端）时，需知晓此数值处理的存在。

### 2. 真空区的数值噪音

即便有了上述修正，在平面波（PW）或 LCAO 计算中，真空区偶尔仍会出现 $10^{-4}$ 量级的微小数值波动。这属于数值计算的正常噪音（Numerical Noise），不代表该处有真实的电子局域化。在可视化时，通过调高 Isosurface level 可以轻松过滤掉这些背景噪音。

### 3. 可视化技巧：Isosurface Level 的选择

初学者常问：“ELF 的等值面应该设为多少？”
*   **默认值**: VESTA 打开时通常默认为 0.5（对应均匀电子气），但这往往看不清成键细节。
*   **推荐值**:
    *   **0.7 - 0.85**: 观察共价键（Covalent Bond）和孤对电子（Lone Pair）的最佳区间。
    *   **0.2 - 0.4**: 观察金属键或弱相互作用区域。
*   **操作**: 在 VESTA 的 `Properties` -> `Isosurfaces` 面板中，手动修改 `Isosurface level` 的值，直到物理图像清晰为止。

### 4. 文件覆盖风险

**警示**: 每次运行 ABACUS，`OUT` 文件夹下的 `ELF.cube` 等文件会被直接覆盖。如果你需要对比不同参数（如不同磁矩设置）下的 ELF，请务必在计算完成后手动重命名或备份这些文件。