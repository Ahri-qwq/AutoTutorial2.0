# 第二章：sDFT 基础流程：单点能计算 (SCF)

在第一章中，我们已经熟悉了 ABACUS 的基本安装与运行。本章我们将进入随机波函数密度泛函理论（Stochastic DFT, sDFT）的核心实战环节。

sDFT 是 ABACUS 处理温稠密物质（Warm Dense Matter, WDM）和超大体系的“杀手锏”。与传统的 Kohn-Sham DFT (KSDFT) 不同，sDFT 不依赖于哈密顿矩阵的对角化，而是通过随机向量在希尔伯特空间中的演化来获取物理量。这种算法特性的改变，意味着我们需要一套全新的参数设置逻辑。

本章将以最简单的 **硅双原子（Si2）** 体系为例，带你跑通 sDFT 的标准流程，并掌握控制随机误差的核心技巧。

---

## 2.1 输入文件概览与核心开关

首先，我们需要明确告诉 ABACUS：“请放弃传统的对角化方法，启动随机求解器。”

我们以 `pw_Si2` 案例为基础。与 KSDFT 一样，你需要准备 `STRU`（结构文件）、`KPT`（K点文件）和赝势文件。唯一的剧变发生在 `INPUT` 文件中。

### 2.1.1 标准 sDFT INPUT 模板

以下是一个标准的 sDFT 计算输入文件示例：

```bash
INPUT_PARAMETERS
#Parameters (1. General)
suffix          Si2_sdft
calculation     scf             # 计算类型：自洽场计算
esolver_type    sdft            # 核心开关：指定求解器为 sdft
pseudo_dir      ../../PP_ORB    # 赝势路径
symmetry        1               # 开启对称性分析

#Parameters (2. Stochastic & Basis)
nbands          0               # KS轨道数：设为0表示纯sDFT
nbands_sto      64              # 随机轨道数：sDFT精度的关键
nche_sto        100             # 切比雪夫展开阶数
method_sto      2               # 算法模式：2为速度优先

#Parameters (3. Iteration & Accuracy)
ecutwfc         50              # 平面波截断能 (Ry)
scf_nmax        20              # 最大电子迭代步数
scf_thr         1e-6            # 收敛阈值

#Parameters (4. Smearing - CRITICAL)
smearing_method fd              # 必须设为 fd (Fermi-Dirac)
smearing_sigma  0.6             # 电子温度 (Ry)，0.6 Ry ≈ 8.16 eV
```

### 2.1.2 核心参数解析

1.  **`esolver_type` (必须为 sdft)**
    *   这是 sDFT 的总开关。默认值为 `ksdft`。一旦设置为 `sdft`，ABACUS 将不再进行波函数对角化，而是启用切比雪夫（Chebyshev）多项式展开求解器。

2.  **`calculation` (scf)**
    *   sDFT 同样支持 `scf`（自洽计算）和 `md`（分子动力学）。但在进行 MD 之前，必须先通过 SCF 确定合适的参数设置。

3.  **`smearing_method` (必须为 fd)**
    *   **原理**：sDFT 天然基于费米-狄拉克（Fermi-Dirac）分布算符的展开。
    *   **限制**：在 sDFT 模式下，`smearing_method` **必须** 设置为 `fd`。不支持 `gaussian` 或 `mp` 等其他展宽方法。
    *   **`smearing_sigma`**：这里的物理意义直接对应**电子温度**。注意单位是 Ry（1 Ry $\approx$ 13.6 eV）。例如温稠密物质模拟中常用的 10 eV 温度，此处应设为 $\approx 0.735$。

---

## 2.2 随机空间与切比雪夫展开设置

sDFT 的计算不再是“精确”的，而是伴随着**随机误差**。作为用户，你的核心任务就是在“计算精度”和“计算成本”之间通过以下两个参数找到平衡点。

### 2.2.1 随机轨道数量 (`nbands_sto`) 与误差控制

`nbands_sto` 定义了用于采样希尔伯特空间的随机向量个数。
*   **规律**：`nbands_sto` 越大，随机噪声越小，结果越精确，但计算耗时线性增加。
*   **经验值**：对于小体系，通常从 64 或 128 开始测试。

#### ⚠️ 关键流程：误差测试 (Error Test)
由于 sDFT 的结果随随机数种子（`seed_sto`）波动，你**不能**仅凭一次计算就断定结果。必须执行以下标准测试流程来确定合适的 `nbands_sto` 和 `ecutwfc`：

1.  **固定参数**：设定一个待测的 `nbands_sto`（如 64）。
2.  **多种子运行**：修改 `seed_sto` 参数（默认为 0，即随机），手动设置 10 个不同的整数（如 1001, 1002, ..., 1010），运行 10 次 SCF 计算。
3.  **统计分析**：收集这 10 次计算的总能量（Total Energy）。
    *   计算**平均值** $\bar{E}$。
    *   计算**标准差** $\sigma_E$。
4.  **判据**：如果标准差 $\sigma_E$ 小于你的目标精度（通常建议 $< 10^{-4}$ Ry/atom），则该 `nbands_sto` 足够；否则需增加轨道数重测。

> **提示**：`ecutwfc` 的收敛性测试也需遵循此流程，比较不同 `ecut` 下的**平均能量**差值。

### 2.2.2 切比雪夫展开阶数 (`nche_sto`)

sDFT 不计算本征值，而是直接计算密度矩阵。这一过程依赖于将费米-狄拉克函数展开为切比雪夫多项式，`nche_sto` 即为展开的阶数。

*   **物理依赖关系（至关重要）**：
    1.  **温度 ($T$)**：温度**越高**，费米分布函数越平滑，需要的 `nche_sto` **越小**（计算越快）。
    2.  **截断能 (`ecut`)**：截断能越大，哈密顿量的谱范围越宽，需要的 `nche_sto` **越大**。

*   **如何检查是否足够？**
    查看输出文件 `running_scf.log`，搜索 `Chebyshev Precision`。
    ```text
    Chebyshev Precision: 1.23e-09
    ```
    推荐该值小于 $10^{-8}$。如果精度不够，请增加 `nche_sto`。

### 2.2.3 进阶技巧：MDFT (混合 DFT)

**何时使用？**
当模拟温度较高，但体系中仍有部分深能级或半核心态（Semi-core states）被完全占据且对化学键有重要贡献时，纯 sDFT 需要极大的 `nbands_sto` 才能准确描述这些定域态。

**如何设置？**
此时应使用 **MDFT**。
*   **设置方法**：在 INPUT 中设置 `nbands > 0`（例如 `nbands 4`）。
*   **效果**：ABACUS 会先用 KSDFT 精确计算最低的 4 条能带（处理深能级），然后用 sDFT 处理剩余的高能级部分。这能显著降低对 `nbands_sto` 的需求，加速收敛。

---

## 2.3 算法选择与并行策略

sDFT 的计算量主要集中在哈密顿量对随机向量的作用上。合理的算法选择和并行设置能带来数倍的性能提升。

### 2.3.1 算法模式 (`method_sto`)

*   **`method_sto 1` (省内存模式)**：
    *   适用于内存受限的大体系。
    *   计算速度相对较慢。
*   **`method_sto 2` (速度优先模式)**：
    *   **默认推荐**。
    *   通过空间换时间，速度较快，但内存消耗较大。

### 2.3.2 并行策略

sDFT 支持多级并行，优先级如下：

1.  **K 点并行 (`kpar`)**：
    *   **优先级：最高**。
    *   如果体系包含多个 K 点，优先设置 `kpar` 等于 K 点总数（或其约数）。这是效率最高的并行方式。

2.  **能带并行 (`bndpar`)**：
    *   **优先级：次高**。
    *   将 `nbands_sto` 个随机轨道分配到不同的进程组中计算。
    *   **设置建议**：默认为 1。在 K 点并行用满后，若还有空闲核数，可尝试增加 `bndpar`。注意 `bndpar` 并非越大越好，需进行简单的 benchmark。

### 2.3.3 风险提示：应力计算 (Stress)

虽然 sDFT 在计算力（Force）和进行 MD 模拟方面已经非常成熟，但涉及晶胞尺寸优化（`calculation cell-relax`）时需要计算应力张量（Stress Tensor）。

*   **注意**：sDFT 的应力计算对随机误差非常敏感。若必须进行晶胞优化，建议：
    1.  务必使用更严格的 `nbands_sto`（通过误差测试确定）。
    2.  查阅最新的官方文档，确认当前版本 sDFT 对应力计算的支持状态及限制。

---

**本章小结**
通过本章，你已经掌握了 sDFT 的核心输入逻辑：
1.  开启 `esolver_type sdft` 和 `smearing_method fd`。
2.  利用“10种子测试法”确定 `nbands_sto`。
3.  根据温度和截断能调整 `nche_sto`。
4.  在深能级主导时启用 MDFT (`nbands > 0`)。

下一章，我们将利用这些知识，实战运行一个温稠密物质的分子动力学（MD）模拟。