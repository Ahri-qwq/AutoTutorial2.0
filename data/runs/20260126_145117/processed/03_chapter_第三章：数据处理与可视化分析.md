# 第三章：数据处理与可视化分析

算出数据只是完成了一半。在上一章中，我们通过“两步法”（SCF + NSCF）成功完成了 DOS 计算。然而，ABACUS 输出的原始数据通常是一堆枯燥的数字，且能量标度往往未对齐。

本章将指导你如何像专家一样处理这些原始数据，特别是至关重要的**费米能级对齐**，并结合 Python 脚本生成可发表级别的图表。我们将以硅（Si）为例，不仅教你怎么画，更教你怎么看懂图中的物理含义。

---

## 3.1 输出文件结构解析

当你的 NSCF 计算顺利结束后（即 `calculation` 设置为 `nscf` 且 `out_dos` 设置为 `1`），ABACUS 会在 `OUT.${suffix}` 文件夹中生成一系列文件。

### 3.1.1 核心文件概览

进入输出目录（假设 `suffix` 为默认的 `ABACUS`）：
```bash
cd OUT.ABACUS
ls
```

你会看到以下关键文件：

1.  **`TDOS` (Total Density of States)**
    *   **描述**：系统的总态密度。
    *   **格式**：纯文本文件，包含三列数据。
        *   第 1 列：能量 (Energy)，单位通常为 **eV**。
        *   第 2 列：总态密度 (Total DOS)，单位为 States/eV。
        *   第 3 列：积分态密度 (Integrated DOS)，即总电子数的累积值。

2.  **`PDOS` (Projected Density of States)**
    *   **描述**：分波态密度（如果使用了 LCAO 基组或开启了投影功能）。
    *   **格式**：包含能量列以及投影到各个原子轨道（如 Si 的 3s, 3p）上的 DOS 分量。
    *   **用途**：用于分析能带主要由哪些轨道贡献（例如分析杂化情况）。

3.  **`running_nscf.log` (或类似日志)**
    *   **描述**：计算运行日志。
    *   **用途**：**非常重要**，我们需要从中提取费米能级（Fermi Energy）的数值。

> **专家提示 (Developer Note)**：
> 请务必区分标准 **KS-DFT** 和 **SDFT (Stochastic DFT)** 的输出。
> *   如果你使用的是标准的平面波或原子轨道基组（绝大多数初学者的场景），请关注上述的 `TDOS` 文件。
> *   如果你在做高温稠密物质或超大体系并使用了 SDFT（随机波函数方法），输出文件结构会有所不同。本教程主要针对标准 KS-DFT 流程。

---

## 3.2 数据后处理：费米能级对齐

这是新手最容易犯错的步骤。

### 3.2.1 为什么要对齐？
ABACUS（以及大多数 DFT 软件）输出的原始能量值是相对于计算原本的参考零点（通常是平均静电势或真空能级）。这意味着你的价带顶（VBM）或费米能级（$E_F$）可能显示为 -5.2 eV 或 +3.1 eV，而不是物理直觉上的 0 eV。

为了方便分析，我们需要将费米能级平移到 0 eV。

### 3.2.2 提取费米能级
费米能级通常记录在 SCF 计算的日志中（因为 NSCF 步骤通常不重新更新费米能级，或者使用 SCF 的电荷密度确定的费米能级）。

执行以下命令查找费米能级：
```bash
# 在 SCF 或 NSCF 的日志文件中查找
grep "Fermi" OUT.ABACUS/running_*.log
```

输出示例：
```text
Fermi energy = 6.2345 eV
```
记下这个数值，例如 $E_F = 6.2345$。

### 3.2.3 数据处理公式
在绘图时，我们需要对 `TDOS` 文件中的第一列（能量）进行如下操作：

$$ E_{\text{plot}} = E_{\text{raw}} - E_F $$

*   对于**金属**：$E_{\text{plot}} = 0$ 处即为费米面，DOS 值通常不为零。
*   对于**半导体/绝缘体**：$E_{\text{plot}} = 0$ 通常位于带隙中（或价带顶，取决于具体定义）。

---

## 3.3 绘图实战与案例分析：硅 (Si)

我们将使用 Python (Matplotlib) 编写一个通用的脚本来处理数据并绘图。

### 3.3.1 Python 绘图脚本 (`plot_dos.py`)

将以下代码保存为 `plot_dos.py`，放在 `OUT.ABACUS` 同级目录下。

```python
import numpy as np
import matplotlib.pyplot as plt

# ================= 参数设置 =================
fermi_energy = 6.2345  # <--- 请替换为你 grep 到的实际数值
filename = 'OUT.ABACUS/TDOS'
# ===========================================

# 读取数据
# ABACUS TDOS 文件通常有头信息，使用 comments='#' 跳过注释行
try:
    data = np.loadtxt(filename, comments=['#', 'S']) # 兼容不同版本的头文件注释
except Exception as e:
    print(f"Error reading file: {e}")
    exit()

energy_raw = data[:, 0]
dos_total = data[:, 1]

# 费米能级对齐
energy_aligned = energy_raw - fermi_energy

# 绘图
plt.figure(figsize=(8, 6))

# 绘制 DOS 曲线
plt.plot(energy_aligned, dos_total, color='k', linewidth=1.5, label='Total DOS')

# 填充价带部分 (E < 0)
plt.fill_between(energy_aligned, 0, dos_total, 
                 where=(energy_aligned < 0), 
                 facecolor='skyblue', alpha=0.3)

# 装饰图表
plt.axvline(x=0, color='r', linestyle='--', linewidth=1, label='$E_F$')
plt.xlabel(r'$E - E_F$ (eV)', fontsize=14)
plt.ylabel('DOS (states/eV)', fontsize=14)
plt.title('Total Density of States (Si)', fontsize=16)
plt.xlim(-10, 10)  # 根据需要调整显示范围
plt.ylim(0, max(dos_total) * 1.1)
plt.legend()
plt.grid(alpha=0.3)

# 保存图片
plt.savefig('Si_DOS.png', dpi=300)
plt.show()
```

### 3.3.2 物理图像分析

运行脚本后，你会得到一张 DOS 图。结合第一章的物理知识，我们来“看图说话”：

1.  **带隙 (Band Gap)**：
    *   观察 $E=0$ 附近。对于硅（半导体），你应该看到一段 DOS 为 0 的区域。
    *   **验证**：如果 $E=0$ 处 DOS 明显不为 0，说明要么计算的是金属，要么费米能级提取错误，或者 smearing（展宽）参数设置过大导致带隙被“涂抹”了。

2.  **范霍夫奇点 (Van Hove Singularities)**：
    *   你会注意到 DOS 曲线中有许多尖锐的峰。
    *   **物理含义**：这些峰对应于能带结构中 $E(k)$ 比较平坦的区域（即 $\nabla_k E \approx 0$）。在这些能量处，电子态密度非常高。
    *   **K点提示**：如果你的 DOS 曲线看起来非常粗糙或充满杂乱的毛刺，通常是因为 NSCF 计算时 **K 点网格不够密**。DOS 计算通常需要比 SCF 计算更密的 K 点（例如 SCF 用 8x8x8，DOS 建议用 12x12x12 或更高）。

3.  **轨道贡献 (PDOS)**：
    *   如果你绘制了 PDOS，你会发现硅的价带深处（低能区）主要由 s 轨道贡献，而靠近费米能级的价带顶则由 s 和 p 轨道强烈的杂化（sp3 杂化）贡献。

---

## 3.4 自动化工具链 (Bohrium & Atomkit)

对于高通量计算或不想手动编写脚本的用户，ABACUS 生态提供了便捷的工具。

### 3.4.1 Bohrium 平台
如果你在 Bohrium 平台上运行 ABACUS，可以使用平台集成的 `abacus_dos_run` 或相关 Notebook 模板。这些工具会自动读取 `INPUT` 中的参数和 `OUT` 文件夹中的数据，一键生成交互式的 DOS 图表。

### 3.4.2 Atomkit 与 ASE
高级用户可以使用 Python 库来处理数据。
*   **ASE (Atomic Simulation Environment)**：虽然 ASE 主要用于结构操作，但结合自定义的 ABACUS 解析器，可以批量处理多个结构的 DOS 数据。
*   **dpdata / Atomkit**：DeepModeling 社区开发的工具包，能够方便地将 ABACUS 数据转换为其他格式以便分析。

---

**下一章预告**：
掌握了电子结构的静态性质（DOS）后，我们将让电子“动”起来。下一章我们将进入**能带结构（Band Structure）**的计算，这是连接理论与角分辨光电子能谱（ARPES）实验的桥梁。我们将学习如何选取高对称路径（K-path），并绘制漂亮的能带图。