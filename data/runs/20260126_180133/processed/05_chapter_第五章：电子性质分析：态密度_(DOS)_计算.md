# 第五章：电子性质分析：态密度 (DOS) 计算

在完成了体系的能量收敛和结构确定后，深入理解材料物理性质的关键一步是分析电子结构。对于温稠密物质（WDM）或高温体系，态密度（Density of States, DOS）是揭示电子在不同能量区间分布情况的核心指标。

在 sDFT（随机密度泛函理论）框架下，DOS 的计算与传统 KSDFT 略有不同。由于我们使用随机轨道（Stochastic Orbitals）而非确定的本征态，DOS 的获取需要依赖切比雪夫（Chebyshev）展开技术进行重构。本章将基于 `186_PW_SDOS_10D10S` 案例，详细讲解如何在 ABACUS 中通过 sDFT/MDFT 方法计算 DOS，并着重讨论误差控制与大体系的内存优化。

---

## 5.1 DOS 计算流程与参数配置

sDFT 的 DOS 计算通常作为 SCF（自洽场）计算的一部分或紧随其后的后处理步骤进行。与传统方法不同，sDFT 不需要非自洽（NSCF）步骤来加密 K 点，因为随机波函数本身覆盖了整个希尔伯特空间。

### 5.1.1 核心输入文件配置

以案例 `186_PW_SDOS_10D10S` 为例，我们需要在 `INPUT` 文件中添加特定的 DOS 参数。以下是一个典型的用于计算 DOS 的输入文件示例：

```bash
INPUT_PARAMETERS
# 1. General Parameters
suffix         autotest
calculation    scf          # 进行自洽计算
esolver_type   sdft         # 核心：指定求解器为 sDFT
method_sto     2            # 推荐：速度较快且支持内存分块的方法
nbands         10           # 若 >0 则为 MDFT，若 =0 则为纯 sDFT
nbands_sto     10           # 随机轨道数量
nche_sto       120          # SCF 过程中的切比雪夫展开阶数
seed_sto       20000        # 固定随机数种子以复现结果

# 2. Basis & Accuracy
basis_type     pw
ecutwfc        20
scf_thr        1e-6
scf_nmax       20

# 3. Smearing (Critical for sDFT)
smearing_method fd          # 必须使用 Fermi-Dirac
smearing_sigma  0.6         # 电子温度 (Ry)

# 4. DOS Parameters (Key Section)
out_dos         1           # 开启 DOS 输出
dos_emin_ev     -20         # DOS 能量窗口下限 (eV)
dos_emax_ev     100         # DOS 能量窗口上限 (eV)
dos_edelta_ev   0.1         # DOS 能量间隔/分辨率 (eV)
dos_sigma       4           # DOS 重构的高斯展宽因子 (eV)
dos_nche        240         # DOS 计算专用的切比雪夫展开阶数
npart_sto       2           # 内存分块参数 (见 5.2 节)
```

### 5.1.2 关键参数详解

1.  **`out_dos`**:
    *   设置为 `1` 时，程序在 SCF 收敛后会自动调用随机切比雪夫方法计算 DOS。
    *   **输出文件**: 结果将保存在输出目录下的 `DOS1_smearing.dat` 文件中。

2.  **`dos_nche` vs `nche_sto`**:
    *   **`nche_sto`**: 用于 SCF 迭代中计算电荷密度。由于只需要积分性质，通常较低的阶数（如 100-120）即可满足能量收敛要求。
    *   **`dos_nche`**: 用于重构 DOS 曲线。为了获得平滑且细节丰富的 DOS 曲线，**通常需要比 SCF 过程更高的展开阶数**（如本例中的 240）。如果该值过小，输出的 DOS 曲线可能会出现非物理的震荡。

3.  **`dos_sigma`**:
    *   这是在重构 DOS 时引入的高斯展宽因子。注意这与电子温度 `smearing_sigma` 不同。`dos_sigma` 用于平滑随机噪声，通常取值在几 eV 量级（取决于 `dos_edelta_ev` 和切比雪夫展开的收敛性）。

### 5.1.3 误差控制流程（**Critical**）

sDFT 的本质决定了结果带有随机误差。**在进行任何物理分析之前，必须验证随机轨道数量（`nbands_sto`）和截断能（`ecutwfc`）的可靠性。**

**标准测试步骤**：
1.  **固定参数**：设定一组 `nbands_sto` 和 `ecutwfc`。
2.  **多种子测试**：使用约 **10 个不同的 `seed_sto`**（如 100, 200, ..., 1000）分别运行 SCF 计算。
3.  **统计分析**：
    *   收集所有计算的总能量。
    *   计算平均值 $\bar{E}$ 和标准差 $\sigma_E$。
4.  **判据**：
    *   如果 $\sigma_E / N_{atom} < 10^{-4}$ Ry（或根据具体精度需求），则认为 `nbands_sto` 足够。
    *   如果误差过大，需增加 `nbands_sto`。注意：误差随 $\frac{1}{\sqrt{N_{sto}}}$ 衰减。

---

## 5.2 大体系 DOS 计算的内存优化

在处理包含数百甚至上千原子的高温体系时，存储巨大的随机波函数矩阵可能会导致内存溢出（OOM）。ABACUS 提供了分块计算策略来解决这一问题。

### 5.2.1 内存分块参数

*   **`method_sto`**: 建议设置为 `2`。这是启用内存优化和更高效算法的前提。
*   **`npart_sto`**: 内存分块数（默认为 1）。
    *   **作用**: 当 `method_sto=2` 时，程序会将 DOS 计算所需的内存需求降为原来的 $1/\text{npart\_sto}$。
    *   **代价**: 计算时间会随着 `npart_sto` 的增加而线性增加。这是一个典型的“时间换空间”策略。

### 5.2.2 配置建议

如果您的计算节点内存为 256 GB，而标准计算预估需要 400 GB 内存：
1.  设置 `method_sto 2`。
2.  设置 `npart_sto 2`（理论内存需求降至 ~200 GB）。
3.  运行计算并监控内存峰值。

---

## 附录：常见问题与专家建议

### 1. MDFT 与 sDFT 的选择策略 (**Critical**)
很多用户困惑何时使用混合方法（MDFT，即 `nbands > 0`）。
*   **纯 sDFT (`nbands=0`)**: 适用于极高温度（如 > 50 eV），此时几乎所有电子都处于高度离域状态，费米狄拉克分布极其平缓。
*   **MDFT (`nbands > 0`)**: 适用于“温稠密”区间（如 1 eV - 20 eV）。在此温度下，虽然价电子被激发，但深层能级（Deep states）或部分价带底的电子仍然具有较强的束缚性质。
    *   **专家建议**: 将低能的束缚态用确定的 Kohn-Sham 轨道（`nbands`）处理，而将高能的连续态用随机轨道（`nbands_sto`）处理。这能显著降低随机误差，加速收敛。

### 2. 切比雪夫阶数 (`nche_sto`) 的依赖关系 (**Critical**)
`nche_sto` 的选取并非任意，它遵循以下物理规律：
*   **与温度成反比**: 温度越高，费米分布越平滑，需要的切比雪夫展开阶数越**小**（计算越快）。
*   **与 Ecut 成正比**: 截断能 `ecutwfc` 越大，哈密顿量的谱范围越宽，需要的展开阶数越**大**。
*   **监控指标**: 检查输出文件 `running_scf.log` 中的 `Chebyshev Precision`，目标是使其小于 `1e-8`。

### 3. 高温赝势的陷阱
*   **内壳层电离**: 在高温下（如 > 50 eV），常温下被视为“冻结”的内壳层电子可能会发生电离。
*   **警告**: 使用标准赝势库（如 SG15, ONCVPSP）进行极高温计算是危险的。务必检查赝势的截断半径和价电子设置，必要时需生成专门的高温赝势（包含更多内层电子作为价电子）。

### 4. Smearing 方法的限制
*   **严格限制**: sDFT 目前**仅支持** `smearing_method fd` (Fermi-Dirac)。
*   **禁止操作**: 切勿使用 `gauss` 或 `mp` (Methfessel-Paxton) 方法，这些方法在随机波函数框架下未被适配，会导致计算错误或崩溃。

### 5. 风险提示：应力与晶胞优化
*   虽然 sDFT 支持计算力（Force）并进行分子动力学（MD）模拟，但**应力（Stress）张量**的计算在 sDFT 中可能存在精度或实现上的限制。
*   **建议**: 若需进行变胞优化（Cell Relaxation, `calculation cell-relax`），请务必先查阅最新的 ABACUS 官方文档确认当前版本 sDFT 是否完全支持应力计算。如果不支持，请仅固定晶胞体积进行离子步弛豫。