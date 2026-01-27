# 第二章：参数详解：ABACUS 中的 Smearing 设置

在第一章完成了环境搭建与基本运行后，我们现在进入实战核心——`INPUT` 参数文件的深度解析。

在 DFT 计算中，**Smearing（展宽）** 是一个看似不起眼，实则对**收敛性（Convergence）**和**能量精度（Accuracy）**有着决定性影响的参数。特别是在处理金属体系或复杂磁性体系时，错误的 Smearing 设置往往是 SCF（自洽场）不收敛或结果荒谬的罪魁祸首。

本章将带你深入理解 ABACUS 中的展宽机制，并重点解决一个困扰了无数新手的“隐形杀手”——**单位问题**。

---

## 2.1 展宽方法的选择 (Method)

在 ABACUS 的 `INPUT` 文件中，`smearing_method` 参数决定了电子在费米能级附近的占据方式。简单来说，它决定了我们如何处理布里渊区积分中的不连续性。

### 核心参数：`smearing_method`

#### 1. `fixed`：绝缘体与半导体的首选
- **适用场景**：具有清晰带隙（Band Gap）的绝缘体、半导体、孤立分子。
- **原理**：使用固定的整数占据数。此时电子数必须严格匹配能带数，不进行任何展宽处理。
- **注意**：如果你计算的是金属，或者带隙非常小的体系，使用 `fixed` 会导致 SCF 极难收敛，甚至震荡发散。

#### 2. `gauss` (Gaussian)：通用且安全的“万金油”
- **适用场景**：任何体系，但主要用于**难以确定带隙的绝缘体**或**初步测试**。
- **原理**：使用高斯函数对费米能级附近的电子态密度进行平滑处理。
- **优缺点**：
    - **优点**：非常稳健，几乎不会导致负的占据数，收敛性较好。
    - **缺点**：为了获得精确的总能量，需要将展宽宽度（Sigma）设得非常小，或者进行 Sigma $\to$ 0 的能量外推。否则，计算出的能量会包含较大的“非物理”熵贡献。

#### 3. `mp` (Methfessel-Paxton)：金属体系的标准配置
- **适用场景**：**金属**、合金、表面吸附体系（涉及金属基底）。
- **原理**：在 Gaussian 基础上引入高阶厄米多项式修正。
- **核心优势**：它能比 Gaussian 方法更好地近似 T=0K 时的阶梯函数。这意味着，即使使用稍大的 Sigma（为了好收敛），算出来的总能量依然非常接近真实的基态能量。
- **警告**：MP 方法可能会产生非物理的**负占据数**或大于1的占据数。这在数学上是为了能量收敛，但在做态密度（DOS）分析时需留意。对于宽带隙绝缘体，**严禁**使用 `mp`，否则可能错误地预测出金属性。

#### 4. `fd` (Fermi-Dirac)：有限温度模拟
- **适用场景**：**高温分子动力学（AIMD）**、需要考虑电子温度效应的体系。
- **原理**：使用真实的费米-狄拉克分布函数。此时 Sigma 具有明确的物理意义——电子温度 ($k_B T$)。

---

## 2.2 展宽宽度的设定与单位陷阱 (Sigma)

这是本章最关键的部分。请务必打起十二分精神，因为这里有一个**极易踩中的陷阱**。

### 核心参数：`smearing_sigma`

该参数定义了展宽的能量范围（即高斯函数的宽度或电子温度）。

### ⚠️ 红色警报：单位是 Rydberg (Ry) 不是 eV！

许多从 VASP 或其他软件转到 ABACUS 的用户，习惯了 eV 单位。在 VASP 中，`SIGMA = 0.05` 意味着 0.05 eV。

**但是，在 ABACUS 中，`smearing_sigma` 的默认单位是 Ry (Rydberg)！**

让我们算一笔账：
- $1 \text{ Ry} \approx 13.605 \text{ eV}$
- 如果你在 ABACUS 中填入 `smearing_sigma 0.05`（误以为是 0.05 eV）：
    - 实际效果 $= 0.05 \times 13.605 \approx \mathbf{0.68 \text{ eV}}$
- **后果**：0.68 eV 的展宽对于大多数金属体系来说**太大了**！这相当于给电子加热到了几千度，会导致总能量严重偏离基态，原子受力计算错误，甚至导致结构优化得到错误的晶格常数。

#### 推荐设置值 (Reference Values)

| 体系类型 | 推荐方法 (`smearing_method`) | 推荐宽度 (`smearing_sigma`) | 对应 eV 值 | 说明 |
| :--- | :--- | :--- | :--- | :--- |
| **绝缘体/半导体** | `fixed` | (忽略) | - | 最准确 |
| **绝缘体 (难收敛)**| `gauss` | `0.001` ~ `0.005` | ~0.01 - 0.07 eV | 尽量小 |
| **金属 (常规)** | `mp` | `0.01` ~ `0.02` | ~0.13 - 0.27 eV | 兼顾收敛与精度 |
| **金属 (难收敛)** | `mp` | `0.02` ~ `0.03` | ~0.27 - 0.40 eV | 需测试能量对 Sigma 的收敛性 |
| **高温 MD** | `fd` | `0.002` ~ `0.01` | 取决于物理温度 | 对应真实温度 |

### 替代参数：`smearing_sigma_temp` (Kelvin)

如果你不想进行 Ry 和 eV 的换算，ABACUS 提供了一个更直观的参数：`smearing_sigma_temp`，单位是 **Kelvin (K)**。

*   **互斥规则**：`smearing_sigma` 和 `smearing_sigma_temp` **不能同时出现**在 INPUT 文件中。
*   **换算关系**：$1 \text{ Ry} \approx 157887 \text{ K}$。
    *   例如：设置 `smearing_sigma 0.01` (Ry) 等效于 `smearing_sigma_temp 1579` (K)。

---

## 2.3 实战策略：收敛性与精度的权衡

在实际计算中，我们经常面临一个 Trade-off（权衡）：
- **Sigma 越大**：电荷密度混合越平滑，SCF 迭代**越容易收敛**，但总能量**越不准确**。
- **Sigma 越小**：结果越接近真实的基态物理，但 SCF 容易在费米面附近发生**电荷震荡**（Charge Sloshing），导致不收敛。

### 策略一：分步收敛法 (Step-down Strategy)
对于极难收敛的磁性金属或团簇体系，建议采用“先粗后细”的策略：

1.  **第一步**：设置较大的 Sigma（如 `0.02` Ry 或 `0.03` Ry），先让 SCF 跑通，得到一个收敛的电荷密度（`SPIN1_CHG`）。
2.  **第二步**：读取上一步的电荷密度，减小 Sigma（如 `0.01` Ry），继续进行 SCF 计算。

### 策略二：能量对比的一致性原则
**这是科研红线**：当你计算形成能、吸附能或反应势垒时，涉及相减的两个体系（例如：`表面+分子` 与 `洁净表面`），**必须使用完全相同的 Smearing 设置**。

*   **错误示范**：
    *   洁净表面用 `smearing_sigma 0.01`。
    *   吸附体系因为难收敛，改用了 `smearing_sigma 0.02`。
    *   **结果**：计算出的吸附能包含了两者的“熵贡献差”，结果无效。

---

## 2.4 输入文件示例

### 场景 A：常规金属计算 (如 Cu, Al, Pt)
这是最常见的金属计算设置，兼顾效率与准确性。

```bash
# INPUT file snippet
calculation        scf
basis_type         lcao

# Smearing settings
smearing_method    mp      # Methfessel-Paxton for metals
smearing_sigma     0.01    # 0.01 Ry ~= 0.136 eV. DO NOT use 0.1!
```

### 场景 B：宽带隙半导体 (如 Si, TiO2)
如果结构合理，首选 `fixed`。如果遇到收敛困难，退化为 `gauss`。

```bash
# INPUT file snippet
calculation        scf

# Smearing settings
smearing_method    fixed   # Strict occupation
# smearing_sigma is ignored here
```

### 场景 C：高温分子动力学 (3000K)
模拟极端条件下的物质状态。

```bash
# INPUT file snippet
calculation        md

# Smearing settings
smearing_method    fd      # Fermi-Dirac
smearing_sigma_temp 3000   # Temperature in Kelvin
# Do NOT set smearing_sigma
```

---

## 本章总结

1.  **选对方法**：绝缘体用 `fixed`，金属用 `mp`，高温/MD 用 `fd`。
2.  **盯紧单位**：`smearing_sigma` 的单位是 **Ry**。**0.01 Ry** 是一个好的起点，千万别填成 0.1 或 0.05（除非你真的知道自己在做什么）。
3.  **保持一致**：对比能量时，Smearing 参数必须“神圣不可侵犯”，严禁随意更改。

下一章，我们将讨论决定计算速度与精度的另一个核心参数：**基组与截断能 (Basis Set & Ecut)**。