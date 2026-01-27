# 第三章：预分析与质量控制 (Quality Control)

在上一章中，我们成功运行了 ABACUS 的分子动力学（MD）模拟。看着输出文件夹里不断增加的 `STRU` 文件，你可能迫不及待地想计算径向分布函数（RDF）或均方位移（MSD）。

**请停下来。**

计算材料学界有一句名言：“**Garbage In, Garbage Out**”。如果你的轨迹数据包含未弛豫的非平衡态结构，或者时间步长换算错误，那么后续所有的分析结果都将是错误的。

本章将带你完成从“生产运行”到“数据清洗”的关键步骤。我们将学习如何判断体系是否达到热力学平衡，以及如何在分析工具 Candela 中正确截断“垃圾数据”。

---

## 3.0 生产阶段：确保有数据可分析 (Production Run)

在进入分析之前，我们必须回顾一下 ABACUS 的 `INPUT` 文件设置。很多初学者在分析阶段遇到困难，往往是因为生产阶段的参数设置不当。

### 关键参数：`md_dumpfreq`
在 ABACUS 进行 AIMD 生产运行时，`md_dumpfreq` 决定了轨迹文件的输出频率。

*   **如果设为 1 (默认)**：每一步都输出结构。对于长轨迹，这会产生巨大的 I/O 开销和硬盘占用。
*   **如果设太大 (如 1000)**：采样点太稀疏，导致 MSD 曲线不够平滑，统计误差大。

**推荐设置**：
```bash
# INPUT 文件片段
calculation     md
md_type         nvt      # 或 nve, npt
md_nstep        10000    # 总步数
md_dt           1.0      # 时间步长 (fs)
md_dumpfreq     10       # 每 10 步输出一帧结构
```

> **注意**：请记住这里的 `md_dt` (1.0 fs) 和 `md_dumpfreq` (10)，它们是后续设置分析参数的基石。

---

## 3.1 热力学平衡判定

拿到轨迹后的第一件事，**绝对不是**直接扔进分析软件，而是绘制能量和温度随时间的变化曲线。

### 为什么需要判定平衡？
MD 模拟通常从一个静态结构（如晶体或随机堆积）开始。在模拟的初期，体系会经历一个剧烈的**弛豫过程（Relaxation）**，此时的结构、压力和能量都在快速漂移，不代表体系在目标温度下的真实状态。这段时间的数据必须被剔除。

### 如何判定？
你需要提取 ABACUS 输出目录（通常为 `OUT.suffix`）中的能量信息。

1.  **查看日志文件**：
    通常信息位于 `OUT.suffix/running_md.log` (或 `running_scf.log`) 中。关注 `E_Fermi`, `E_KohnSham`, `Kinetic Energy`, `Temperature` 等关键词。

2.  **判定标准**：
    *   **温度 (Temperature)**：在设定温度（如 300 K）上下波动，且没有明显的长期上升或下降趋势。
    *   **势能 (Potential Energy)**：能量应在一个常数值附近波动（Fluctuation）。如果能量还在持续下降，说明体系仍在寻找稳态，尚未平衡。

### 实战操作
你可以使用简单的 `grep` 命令提取数据并绘图（使用 Python 或 Gnuplot）：

```bash
# 提取温度 (假设在 OUT.ABACUS/running_md.log 中)
grep "Temperature" OUT.ABACUS/running_md.log | awk '{print $2}' > temp.dat

# 提取总能
grep "total energy" OUT.ABACUS/running_md.log | awk '{print $NF}' > etot.dat
```

**观察结论**：
*   如果前 1000 步温度剧烈震荡，从第 1001 步开始稳定在 300K ± 20K，那么**前 1000 步就是必须被截断的“预平衡期”**。

---

## 3.2 轨迹截断与 Candela 设置

确认了平衡点后，我们进入分析阶段。我们将使用 **Candela** (一款高效的 MD 后处理工具) 来读取 ABACUS 的轨迹。

### 核心策略：平衡期截断 (Burn-in Truncation)

假设你的总步数 (`md_nstep`) 为 10000，`md_dumpfreq` 为 10，那么你总共拥有 1000 帧结构。
根据 3.1 节的分析，体系在第 2 ps（即第 2000 步，对应第 200 帧）达到平衡。

我们需要告诉 Candela：**忽略前 200 帧**。

### Candela 参数详解

在 Candela 的配置文件（通常是 `.toml` 或 Python 脚本）中，有几个极易混淆的参数用于控制读取范围。

#### 1. `geo_ignore`: 忽略初始帧（**最常用**）
这是实现平衡期截断的标准参数。它表示在读取轨迹时，跳过开头的多少帧。

*   **物理含义**：剔除未弛豫的非平衡态数据。
*   **设置示例**：如果前 200 帧未平衡，则设置 `geo_ignore = 200`。

#### 2. `geo_1` 和 `geo_2`: 文件索引范围（**风险提示**）
ABACUS 在长时间模拟时，可能会将轨迹拆分为多个文件（取决于版本和设置）。
*   `geo_1`: 起始文件索引。
*   `geo_2`: 结束文件索引。
*   **区别**：`geo_ignore` 是丢弃帧数（Frame），而 `geo_1/2` 是选择读取哪些文件（File）。通常情况下，我们读取所有文件，但丢弃第一个文件中的前几百帧，因此主要使用 `geo_ignore`。

---

## 3.3 锦囊妙计：MSD 时间步长的正确换算

这是新手在计算扩散系数（Diffusion Coefficient）时最容易犯错的地方，直接导致结果差 1000 倍甚至更多。

### 问题背景
ABACUS 的 `md_dt` 单位是 **飞秒 (fs)**。
Candela (以及大多数分析软件) 计算 MSD 时，需要知道两帧之间的时间间隔，通常单位要求为 **皮秒 (ps)**。

### 计算公式
你需要设置 Candela 中的 `msd_dt` 参数。

$$ \text{msd\_dt (ps)} = \frac{\text{md\_dt (fs)} \times \text{md\_dumpfreq}}{1000} $$

### 示例
*   **ABACUS 设置**:
    *   `md_dt = 2.0` (fs)
    *   `md_dumpfreq = 10` (每 10 步存一次)
*   **物理间隔**: 实际上每 $2.0 \times 10 = 20$ fs 存储一帧。
*   **Candela 设置**:
    *   `msd_dt = 0.02` (ps)

**代码示例 (Candela 配置)**：

```toml
[input]
# 假设 ABACUS 输出了 STRU_MD_0, STRU_MD_1...
fmt = "abacus"
path = "./OUT.ABACUS"

# 质量控制：截断前 200 帧不平衡轨迹
geo_ignore = 200  

[calc]
# 必须开启 MSD 计算
do_msd = true

# 关键：时间步长换算 (2.0 fs * 10 / 1000)
msd_dt = 0.02 

# 平滑处理
msd_smooth = true
```

---

## 本章总结

1.  **先检查，后分析**：拿到轨迹先画温度/能量曲线，确认体系已达到热力学平衡。
2.  **舍得丢弃**：未平衡的轨迹是噪声，必须通过 `geo_ignore` 剔除。
3.  **算对时间**：严格按照公式 $\text{md\_dt} \times \text{dumpfreq} / 1000$ 设置分析软件的时间步长。

做好了这些质量控制，你的 RDF 和 MSD 曲线才具有物理意义。下一章，我们将正式进入**结构特性分析**，解读径向分布函数背后的微观结构。