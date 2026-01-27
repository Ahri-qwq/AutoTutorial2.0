# 第二章：确立“标准答案”——平面波 (PW) 基准计算

在进入 LCAO（原子轨道中心线性组合）的高效计算世界之前，我们必须先停下来，问自己一个关键问题：**我们怎么知道算得准不准？**

LCAO 方法的核心优势是效率，但这种效率是建立在对基组进行截断（Rcut）和近似基础上的。如果没有一个“标准答案”作为参照，盲目追求 LCAO 的速度是非常危险的。在 ABACUS 中，我们拥有一个独特的优势：**同一个软件框架既支持 LCAO 也支持平面波（Plane Wave, PW）基组**。

平面波基组在数学上具有完备性（Systematically Improvable），只要动能截断（Energy Cutoff）足够大，它就能给出在该泛函和赝势下的“精确解”。因此，本章的目标非常明确：**运行一个高精度的 PW 计算，获取总能量和电子结构的“标准答案”，作为后续评估 LCAO 轨道精度的标尺。**

---

## 2.1 准备 PW 计算环境

与 LCAO 模式相比，PW 模式最大的特点是**不需要轨道文件（.orb）**。它只需要赝势文件（.upf/.vwr）和足够密集的积分网格。

### 2.1.1 修改 INPUT 文件

要将 ABACUS 切换到平面波模式，我们需要在 `INPUT` 文件中进行关键设置。以下是一个用于基准测试的标准 `INPUT` 模板：

```bash
INPUT_PARAMETERS
# 1. 基础任务设置
suffix          benchmark_pw   # 任务后缀，区分输出文件
calculation     scf            # 自洽场计算
basis_type      pw             # 【核心参数】设置为平面波模式
symmetry        1              # 开启对称性分析以加速

# 2. 精度控制（基准测试核心）
ecutwfc         100            # 【关键】平面波动能截断，单位 Rydberg
scf_thr         1.0e-8         # 自洽收敛阈值，基准测试建议设严一点
scf_nmax        100            # 最大迭代步数

# 3. 涂抹与混合（根据体系调整）
smearing_method gaussian
smearing_sigma  0.01
mixing_type     pulay
mixing_beta     0.7

# 4. 路径设置
pseudo_dir      ./             # 赝势文件路径
# orbital_dir   ./             # PW 模式下此参数会被忽略，但保留也无妨
```

#### 关键参数详解：
*   **`basis_type pw`**: 这是告诉 ABACUS 放弃原子轨道，改用平面波基组展开波函数的开关。
*   **`ecutwfc` (Energy Cutoff for Wavefunctions)**:
    *   **含义**: 限制平面波基组的动能上限。
    *   **基准策略**: 在日常计算中，我们可能会取 50-60 Ry（取决于赝势），但在做 **Benchmark（基准）** 时，我们必须“杀鸡用牛刀”。建议设置为 **80-100 Ry** 甚至更高，以确保基组完备性误差降到最低，从而保证我们看到的误差纯粹来自 LCAO 基组本身，而不是 PW 没算准。

### 2.1.2 检查 STRU 文件

在 PW 模式下，`STRU` 文件中的 `NUMERICAL_ORBITAL` 部分虽然可以保留，但程序在运行时会**自动忽略**它。程序仅读取 `ATOMIC_SPECIES` 中的赝势信息。

**STRU 示例**：
```bash
ATOMIC_SPECIES
Si  28.086  Si_ONCV_PBE-1.0.upf  # 仅依赖赝势

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb       # 在 basis_type=pw 时，此行被忽略

LATTICE_CONSTANT
10.2                             # Bohr

LATTICE_VECTORS
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5

ATOMIC_POSITIONS
Direct
Si  0.00  0.00  0.00 1 1 1
Si  0.25  0.25  0.25 1 1 1
```

---

## 2.2 获取基准物理量与收敛标准

运行计算后（例如执行 `abacus`），我们需要从输出文件（如 `OUT.benchmark_pw/running_scf.log` 或标准输出）中提取关键数据。

### 2.2.1 提取总能量 (Total Energy)

这是最核心的指标。
```bash
grep "final etot" OUT.benchmark_pw/running_scf.log
```
输出示例：
```text
!FINAL_ETOT_IS -231.56781234 Ry
```
记录下这个数值，记为 $E_{PW}$。

### 2.2.2 确立“高精度”标准 (The Gold Standard)

在后续章节中，我们将运行 LCAO 计算并获得 $E_{LCAO}$。我们定义能量误差 Delta 为：
$$ \Delta E = \frac{|E_{LCAO} - E_{PW}|}{N_{atoms}} $$

根据计算材料学的经验法则（参考资料4），我们确立以下**收敛标准**：

1.  **高精度标准**: $\Delta E < 1 \text{ meV/atom}$
    *   这是发表级数据通常追求的精度，意味着 LCAO 基组的质量已经非常接近完备基。
2.  **常规计算标准**: $\Delta E \approx 10 \text{ meV/atom}$
    *   对于一些定性的趋势分析或大体系预筛选，这个精度有时是可以接受的，但在 ABACUS 中我们通常争取做到更好。

### 2.2.3 力的基准（可选但推荐）

如果你的后续任务涉及结构弛豫或分子动力学，仅看能量是不够的，还需要检查力（Force）的误差。
在 `OUT.benchmark_pw/` 目录下会生成力的输出文件（通常在日志或单独文件中）。记录下原子受力，作为 $F_{PW}$。

---

## 2.3 专家锦囊：避坑指南

在这一步，很多初学者容易混淆概念，导致后续测试失效。请务必阅读以下核心逻辑：

### 🔴 核心强调：文件名即参数
在 ABACUS 的 LCAO 模式中，**轨道截断半径（Rcut）是写死在轨道文件里的**。
*   **错误理解**: 在 `INPUT` 里找一个叫 `rcut` 的参数修改数值。
*   **正确操作**: 如果你想测试 6.0 au 和 7.0 au 的区别，你必须**更换轨道文件**。
    *   例如：从使用 `Si_gga_6au.orb` 改为在 STRU 中引用 `Si_gga_7au.orb`。
    *   *风险提示*: 目前 ABACUS 不支持在 INPUT 中直接对长轨道文件进行截断（除非使用特定开发版的高级功能，一般不推荐）。请务必准备好不同 Rcut 的轨道文件进行测试。

### 💡 概念辨析：ecutwfc vs. Rcut
*   **`ecutwfc` (Energy Cutoff)**: 控制积分网格的密度（Grid Mesh）。在 LCAO 中，它决定了实空间积分网格的精细度。
*   **`Rcut` (Cutoff Radius)**: 控制原子轨道的空间广度。
*   **常见误区**: 发现 LCAO 结果不准，拼命增加 `INPUT` 里的 `ecutwfc`。
*   **专家建议**: 如果 `ecutwfc` 已经达到 50-60 Ry，继续增加它对 LCAO 结果的改善微乎其微。此时**瓶颈通常在轨道基组本身**（即 Rcut 太小或 Zeta 数不足）。请优先检查轨道文件！

### ⚡ 效率与精度的折中建议
根据经验（参考资料7），不同体系对 Rcut 的敏感度不同：
*   **固体/体相材料 (Solids)**: 由于近邻原子的轨道重叠提供了额外的基组完备性，通常 **6.0 - 8.0 au** 的 Rcut 就能达到很好的精度。
*   **分子/团簇/表面 (Molecules/Surfaces)**: 或者是低维材料，原子周围有大量真空。为了描述真空中的波函数衰减，通常需要更大的 Rcut，建议 **9.0 - 10.0 au** 甚至更大，以减少基组重叠误差（BSSE）和描述长程效应。

---

**本章小结**
现在，你已经拿到了 $E_{PW}$ 这个“标准答案”。在下一章，我们将正式切换回 LCAO 模式，通过调整轨道文件，看看你的 LCAO 计算能多接近这个标准答案。