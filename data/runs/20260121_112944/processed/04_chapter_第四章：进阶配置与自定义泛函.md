# 第四章：进阶配置与自定义泛函

欢迎来到《ABACUS 实战教程》的进阶篇章。作为一名资深开发者，我深知在实际科研中，标准的 HSE06 或 PBE0 往往只是起步。面对强关联体系、宽禁带半导体或复杂的异质结，我们经常需要“魔改”泛函参数。同时，杂化泛函计算的高昂成本和收敛困难也是困扰无数研究生的噩梦。

本章将带你深入 ABACUS 的内核配置，掌握自定义杂化泛函的“手术刀”，并学会如何驯服那些难以收敛的计算任务。

---

## 4.1 自定义杂化参数：打造你的专属泛函

在 ABACUS 中，虽然我们可以通过 `dft_functional` 关键词直接调用 `hse` 或 `pbe0`，但这些预设值并不总是万能的。例如，在计算某些宽禁带氧化物时，标准的 25% 混合比例可能低估了带隙。此时，我们需要手动调控 Hartree-Fock (HF) 交换项的混合比例（Alpha）和范围分离参数（Omega）。

### 4.1.1 核心参数详解

ABACUS 的 LCAO 模块基于 LibRI 和 Libxc 库，通过以下参数控制杂化泛函的物理行为：

1.  **`exx_hybrid_alpha` (混合比例)**
    *   **物理意义**：决定了在 DFT 交换关联泛函中掺入多少比例的精确 HF 交换能（Exact Exchange）。
    *   **默认值**：
        *   若 `dft_functional` 为 `hf`，默认为 `1.0`。
        *   若 `dft_functional` 为 `pbe0` 或 `hse`，默认为 `0.25`（即 25%）。
    *   **实战场景**：如果你需要使用 PBE0-10%（即 hybrid mixing 为 0.1），需显式设置此参数。

2.  **`exx_hse_omega` (范围分离参数)**
    *   **物理意义**：HSE 泛函特有的参数，用于控制库伦相互作用的衰减范围（Range-separation parameter, $\omega$）。它将库伦势分为短程（Short-range）和长程（Long-range）部分，HF 交换项仅作用于短程。
    *   **默认值**：`0.11` ($bohr^{-1}$)，对应标准的 HSE06 泛函。
    *   **实战场景**：调整 $\omega$ 可以改变屏蔽长度。$\omega \to 0$ 回退到 PBE0；$\omega \to \infty$ 回退到 PBE。

### 4.1.2 实战案例：自定义 PBE0 和 HSE

假设我们需要计算一个特殊的半导体材料，文献建议将 HF 混合比例提高到 35%，或者微调 HSE 的屏蔽参数。

#### 场景 A：自定义 PBE0 (Mixing = 0.35)

```bash
INPUT_PARAMETERS
# ... 基础参数 ...
basis_type      lcao
dft_functional  pbe0           # 基础泛函模板
exx_hybrid_alpha 0.35          # 覆盖默认的 0.25，设置为 35%
```

#### 场景 B：自定义 HSE (Omega = 0.2, Mixing = 0.3)

```bash
INPUT_PARAMETERS
# ... 基础参数 ...
basis_type      lcao
dft_functional  hse            # 基础泛函模板
exx_hybrid_alpha 0.30          # 设置混合比例为 30%
exx_hse_omega    0.20          # 修改范围分离参数
```

> **教授提示 (Professor's Note)**：
> 在修改 `exx_hse_omega` 时请务必小心。这个参数对带隙的影响非常敏感。除非你有明确的文献支持或拟合需求，否则建议保持默认值 `0.11`。此外，SCAN0 泛函在不同文献中定义不同（有的取 0.1，有的取 0.25），使用时请务必通过 `exx_hybrid_alpha` 明确指定你的物理定义。

---

## 4.2 收敛性问题处理：双层循环与单层循环

杂化泛函计算最令人头疼的问题之一就是 SCF 不收敛。ABACUS 针对 LCAO 基组下的杂化泛函提供了两种独特的自洽迭代策略：**双层循环**和**单层循环**。理解它们的区别是解决收敛问题的关键。

### 4.2.1 双层循环 (Double Loop)
这是 ABACUS 的默认或常用策略（取决于具体版本配置，建议显式指定）。

*   **逻辑**：
    *   **外层循环 (Blue Loop)**：更新构建精确交换势（EXX Potential）。
    *   **内层循环 (Red Loop)**：固定 EXX 势，仅更新 GGA 部分的电荷密度，直到达到 `scf_thr`。
*   **参数设置**：
    ```bash
    exx_separate_loop 1
    ```
*   **优点**：**内存占用较低**。因为内层循环不需要反复计算昂贵的 EXX 积分。
*   **缺点**：对于某些复杂体系，外层循环可能发生震荡，导致总能量无法收敛。

### 4.2.2 单层循环 (Single Loop)

*   **逻辑**：
    *   在每一步 SCF 迭代中，同时更新 GGA 密度和 EXX 势。
*   **参数设置**：
    ```bash
    exx_separate_loop 0
    ```
*   **优点**：**收敛稳定性更好**。适用于那些在双层循环中外层震荡剧烈的“钉子户”体系。
*   **缺点**：**内存消耗大**。需要额外存储历史步数的密度矩阵用于 Mixing，峰值内存可能比双层循环高出 10%~20%。

### 4.2.3 教授推荐的收敛策略 (Best Practice)

面对一个未知的杂化泛函计算任务，建议遵循以下流程：

1.  **起手式**：先用纯 GGA (如 PBE) 跑完弛豫或自洽，得到收敛的电荷密度文件 (`SPIN*_CHG.cube` 或密度矩阵)。
2.  **读入密度**：在杂化泛函计算中设置 `start_charge atomic` (如果从头算) 或 `file` (读取 PBE 结果，推荐)。
3.  **尝试双层循环**：设置 `exx_separate_loop 1`。这是效率优先的选择。
4.  **切换单层循环**：如果发现 SCF 能量在震荡（即 `dEtot` 不下降），果断将 `exx_separate_loop` 改为 `0`。

```bash
INPUT_PARAMETERS
# 针对难收敛体系的推荐配置
dft_functional   hse
exx_separate_loop 0          # 开启单层循环模式以增强稳定性
scf_thr          1e-6        # 杂化泛函通常比 GGA 难收敛，建议初期阈值不要设得过低（如 1e-8）
mixing_type      pulay       # 确保使用 Pulay 混合
mixing_beta      0.4         # 如果震荡依旧，适当降低混合因子
```

---

<!-- APPENDIX_START -->
## 附录：常见问题排查 (Troubleshooting)

在杂化泛函的征途上，报错是家常便饭。以下是 ABACUS 开发组总结的“急救清单”。

### A.1 常见报错解析

| 错误信息 (Error Message) | 原因分析 (Root Cause) | 解决方案 (Solution) |
| :--- | :--- | :--- |
| `Unrecognized exchange-correlation functional 'HSE'` | 编译时未链接 Libxc 库。 | 重新编译 ABACUS，确保 `LIBXC_DIR` 正确指向 Libxc 安装路径。 |
| `compile with libri to use hybrid functional in lcao basis` | 编译时未链接 LibRI 库。 | LCAO 杂化泛函依赖 LibRI。请检查 `make` 或 `cmake` 配置，确保包含了 `-DUSE_LIBRI=ON` (具体视 CMakeList 而定) 并链接了 LibRI。 |
| `segmentation fault` (在 EXX 计算阶段) | 内存溢出 (OOM)。 | 杂化泛函对内存极度饥渴。请参考下文的“性能清单”调整并行策略。 |

### A.2 性能与内存优化清单 (Checklist)

如果你的计算太慢或者节点内存爆红，请按顺序检查以下设置：

1.  **并行策略 (至关重要)**：
    *   **现象**：内存不足。
    *   **对策**：EXX 计算部分对 OpenMP 线程并行更友好，而对 MPI 进程并行内存开销大。
    *   **黄金法则**：**少 MPI，多 OpenMP**。
    *   **示例**：在一个 56 核的节点上，推荐 `mpirun -np 1 -env OMP_NUM_THREADS=56 abacus` (极端省内存) 或 `mpirun -np 4 -env OMP_NUM_THREADS=14 abacus` (平衡速度与内存)。千万不要使用 56 个 MPI 进程！

2.  **实数/复数优化**：
    *   **参数**：`exx_real_number`
    *   **检查**：如果体系具有中心反演对称性且 `gamma_only` 为 1，或者你能确定密度矩阵虚部为零，设置 `exx_real_number 1`。这能将 EXX 计算速度提升约 2-3 倍。

3.  **基组与截断半径**：
    *   **检查**：EXX 计算量与基组数目的四次方 ($N^4$) 成正比，与截断半径的九次方 ($R_{cut}^9$) 成正比。
    *   **对策**：不要盲目使用 `tzdp` 或过大的 `Rcut`。对于大多数杂化泛函计算，`dzp` 和默认的 `Rcut` (如 6-7 a.u.) 已经足够。

4.  **辅助格点加速**：
    *   **参数**：`exx_ccp_rmesh_times`
    *   **技巧**：该参数控制库伦势积分的径向格点密度。默认值（HSE 为 1.5）通常是安全的。如果极慢，可确认是否被误设为了高值（如 5.0）。

通过本章的配置，你应该能够驾驭绝大多数杂化泛函计算任务。下一章，我们将探讨如何利用 ABACUS 进行激发态性质的计算。