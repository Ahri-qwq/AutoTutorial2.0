# 第四章：场景实战：不同体系的配置策略

这是一个非常关键的章节。在计算材料学中，电子占据数的展宽（Smearing）处理不仅影响自洽迭代（SCF）的收敛性，更直接决定了总能量的物理意义。

作为开发者，我经常看到用户（尤其是从 VASP 迁移过来的用户）因为单位混淆或策略错误导致计算结果偏差巨大。本章将手把手教你如何针对不同体系“对症下药”。

---

# 第四章：场景实战：不同体系的配置策略

在 DFT 计算中，我们通过数值积分来计算布里渊区内的电子密度。对于绝缘体，费米能级位于带隙中，积分相对简单；但对于金属，费米面穿过能带，导致被积函数发生阶跃，数值积分极难收敛。

为此，我们引入了 **Smearing（展宽）** 技术。本章将指导你如何在 ABACUS 中针对绝缘体、半导体和金属选择最正确的配置。

## 4.0 核心原则与单位陷阱（必读）

在开始具体场景前，请务必铭记以下两条“铁律”，它们是 ABACUS 用户最容易犯错的地方：

### 1. 单位是 Ry，不是 eV！
这是 ABACUS 与 VASP 等软件最大的区别之一。
- **参数**: `smearing_sigma`
- **单位**: **Rydberg (Ry)**
- **换算**: 1 Ry $\approx$ 13.605 eV

> **警示**：如果你习惯了 VASP 中的 `SIGMA = 0.05` (eV)，直接在 ABACUS 中设置 `smearing_sigma 0.05`，你实际上设置了 $0.05 \times 13.605 \approx 0.68$ eV 的巨大展宽！这会导致物理结果完全失真。
> **推荐**: ABACUS 中常用的 Sigma 范围通常在 **0.01 Ry (约 0.136 eV)** 到 **0.02 Ry (约 0.27 eV)** 之间。

### 2. 能量对比一致性原则
当你计算形成能、吸附能或反应势垒时，参与加减运算的所有体系（如：吸附物、基底、吸附体系）**必须使用完全相同的 Smearing 方法和 Sigma 值**。严禁为了让某个结构好收敛而单独增大它的 Sigma，这会导致能量基准面不一致，计算出的 $\Delta E$ 毫无意义。

---

## 4.1 绝缘体与半导体体系

对于有带隙的体系，我们的目标是尽量减少人为引入的电子分数占据，还原其绝缘/半导体本征特性。

### 场景 A：宽带隙绝缘体 (Wide Bandgap Insulators)
**典型体系**: SiO$_2$, NaCl, MgO
**物理特征**: 费米能级深深位于带隙之中，价带全满，导带全空。

**推荐配置**:
使用 `fixed` 方法。这对应于绝对零度下的阶跃函数，不引入任何展宽，物理上最严谨。

```bash
# INPUT file snippet
smearing_method  fixed  # 固定占据数
# smearing_sigma  # fixed 模式下该参数无效，可不写
```

### 场景 B：小带隙半导体 (Small Bandgap Semiconductors)
**典型体系**: Si, Ge, GaAs, TiO$_2$
**物理特征**: 虽然有带隙，但在 SCF 迭代初期，电荷密度波动可能导致能隙闭合或误判。

**推荐配置**:
推荐使用 `gauss` (高斯展宽) 配合极小的 Sigma。这样既能提供一定的数值稳定性，帮助 SCF 收敛，又因为 Sigma 很小，对总能量的修正（Entropy term）微乎其微。

```bash
# INPUT file snippet
smearing_method  gauss  # 高斯展宽
smearing_sigma   0.01   # 0.01 Ry ≈ 0.136 eV，足够小以保证精度
```

---

## 4.2 金属体系的标准配置

**典型体系**: Al, Cu, Fe, Pt, 石墨烯
**物理特征**: 费米面穿过能带，电子在费米面附近存在部分占据。

### 为什么不推荐 Gauss？
对于金属，如果使用 Gauss 方法，为了收敛通常需要较大的 Sigma。但 Gauss 方法引入的能量误差与 $\sigma^2$ 成正比，会导致总能量严重偏离基态。

### 推荐方案：Methfessel-Paxton (MP)
ABACUS 提供了 `mp` (一阶) 和 `mp2` (二阶) 方法。MP 方法通过构造特殊的数学函数（允许出现非物理的负占据数作为数学修正），使得总能量对 Sigma 的依赖关系变为 $\sigma^3$ 或更高阶。这意味着即使使用稍大的 Sigma 换取收敛性，计算出的能量依然非常接近真实的基态能量。

**推荐配置**:

```bash
# INPUT file snippet
smearing_method  mp     # 推荐用于金属弛豫
smearing_sigma   0.015  # 约 0.2 eV，兼顾收敛与精度
```

*   **常规金属**: `smearing_sigma` 建议在 0.01 ~ 0.02 Ry 之间。
*   **磁性金属**: 往往需要稍大的 Sigma (如 0.02 Ry) 来抑制磁矩翻转带来的电荷震荡。

---

## 4.3 难收敛体系的“两步走”策略

在处理复杂的磁性界面、大表面积的金属 Slab 或缺陷体系时，经常遇到 SCF 不收敛（Charge Sloshing）的情况。此时，“增大 Sigma”是有效的救急手段，但会牺牲精度。

为了平衡收敛性与精度，建议采用 **“两步走” (Two-Step)** 策略：

### 第一步：粗略收敛 (Coarse SCF)
使用较大的 Sigma 快速获得一个大致正确的电荷密度分布。

*   **INPUT 设置**:
    ```bash
    calculation     scf
    smearing_method mp
    smearing_sigma  0.03   # 较大展宽 (约 0.4 eV) 助收敛
    out_chg         1      # 必须输出电荷密度
    ```
*   **结果**: 运行结束后，目录下会生成 `SPIN1_CHG.cube` (或 `SPIN*_CHG`) 等文件。

### 第二步：高精度计算 (Fine SCF)
读取上一步的电荷密度作为初猜，减小 Sigma 进行最终计算。

*   **INPUT 设置**:
    ```bash
    calculation     scf
    init_chg        file   # 关键：从文件读取电荷密度
    smearing_method mp
    smearing_sigma  0.01   # 调回标准精度
    ```

> **专家提示**: 如果 Smearing 策略失效，请检查 `mixing_beta`。对于难收敛体系，通常需要同时降低 `mixing_beta` (例如从默认值降至 0.2 或 0.1)。

---

## 附录：常见问题与避坑指南

### 1. DOS 计算与 SCF 计算的区别
很多用户混淆了 SCF 过程中的展宽和画态密度（DOS）时的展宽。
*   **`smearing_sigma` (Ry)**: 用于 SCF 自洽迭代，改变电子占据数，**直接影响总能量和受力**。
*   **`dos_sigma` (eV)**: 仅用于 `calculation nscf` 或 `calculation dos` 后处理阶段。它只是为了让画出来的 DOS 曲线更平滑，**不改变系统的物理状态**。
*   **注意**: `dos_sigma` 的单位通常是 **eV**（请以具体版本文档为准，ABACUS 部分后处理工具使用 eV），而 INPUT 中的 `smearing_sigma` 永远是 **Ry**。

### 2. 结构弛豫中的“受力滞后”
当你使用较大的 `smearing_sigma` (如 > 0.02 Ry) 进行结构弛豫时，可能会发现能量已经收敛，但原子受力很难降到极低（如 0.01 eV/Ang）。
*   **原因**: 高温（大 Sigma）引入了电子熵效应，原子感受到的是“自由能面”上的力，而非单纯的势能面。
*   **对策**: 在结构弛豫的最后阶段，务必减小 Sigma 再次优化，以获得准确的基态结构。

### 3. 报错排查：Parameter Conflict
如果你在 INPUT 中同时设置了：
```bash
smearing_method  fixed
smearing_sigma   0.01
```
ABACUS 可能会在日志中给出警告或忽略 `smearing_sigma`。请养成良好的习惯：**只有在 `gauss`, `mp`, `mp2`, `fd` 模式下才设置 sigma 值**。

### 4. 关于 Fermi-Dirac (`fd`)
你可能会在文档中看到 `smearing_method fd`。这是费米-狄拉克分布。
*   **用途**: 主要用于**有限温度 DFT** (Finite Temperature DFT) 模拟，例如温稠密物质（Warm Dense Matter）。
*   **常规计算**: 对于常规的材料基态计算（模拟 0K 性质），**不推荐**使用 `fd`，因为它的能量拖尾效应比 `mp` 严重得多，需要极小的 Sigma 才能准确，但这又会导致收敛困难。金属计算请首选 `mp`。