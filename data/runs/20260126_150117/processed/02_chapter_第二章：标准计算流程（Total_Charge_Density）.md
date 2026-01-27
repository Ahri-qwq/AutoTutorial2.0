# 第二章：标准计算流程（Total Charge Density）

在第一性原理计算中，**差分电荷密度（Charge Density Difference, CDD）** 是分析原子间成键、电荷转移以及界面相互作用最直观的工具之一。

然而，新手最容易陷入的误区是试图寻找一个直接输出“差分电荷”的开关。在此必须明确一个核心概念：**ABACUS 内核不直接输出差分电荷密度文件**。我们需要通过严格控制的“三步走”流程，分别计算整体与分体系的电荷密度，最后在后处理阶段进行网格点对点的相减。

本章将带你完成这一流程的核心部分：获取三个在空间上严格对齐的 `.cube` 电荷密度文件。

---

## 2.1 步骤一：整体体系（AB）的 SCF 计算

这是所有后续分析的基准（Reference）。我们需要对由 A 部分和 B 部分组成的完整体系（System AB）进行自洽场（SCF）计算。

### 2.1.1 核心参数设置

在 `INPUT` 文件中，除了常规的收敛参数（如 `ecutwfc`, `scf_thr`）外，必须特别关注以下三个控制项：

```bash
INPUT_PARAMETERS
# ... (常规参数如 ecutwfc, mixing_type 等) ...

calculation     scf         # 进行自洽计算
symmetry        -1          # [关键] 强制关闭对称性分析
out_chg         1           # 输出电荷密度文件 (通常为 density.cube 或 CHARGE)
basis_type      lcao        # 或 pw，根据实际需求选择
```

### 2.1.2 为什么必须设置 `symmetry -1`？
这是本教程最重要的经验总结之一。
*   **默认行为**：ABACUS 默认会分析晶体对称性以加速计算（`symmetry 1`）。这可能会导致程序自动旋转坐标系、移动原点或通过对称性操作减少存储的网格点。
*   **灾难后果**：如果 System AB 开启了对称性，而 System A（由于原子缺失）对称性降低，两者输出的电荷密度网格（FFT Grid）在实空间中可能发生**旋转或错位**。此时进行差分相减（$\rho_{AB} - \rho_A - \rho_B$），得到的将是毫无物理意义的噪声。
*   **解决方案**：设置 `symmetry -1` 强制程序使用原始输入的 `STRU` 坐标系，不进行任何对称性约化，确保空间网格的绝对锚定。

### 2.1.3 运行与记录
提交任务运行完成后，请务必查看输出日志（如 `OUT.AB/running_scf.log`），找到 FFT 网格的维度信息。这对于下一步至关重要。

> **日志示例（需记录的数值）：**
> ```text
> FFT grid for charge/potential:  108 * 108 * 144
> ```
> *注：这代表 x, y, z 方向的网格点数分别为 nx=108, ny=108, nz=144。*

---

## 2.2 步骤二：分体系（A 与 B）的独立计算

接下来，我们需要分别计算单独存在的 A 部分和 B 部分的电荷密度。此步骤的成败在于：**控制变量**。

### 2.2.1 修改 STRU 文件：移除而非移动
以计算 System A 为例。我们需要基于 System AB 的 `STRU` 文件进行修改：
1.  **保持晶胞不变**：`ATOMIC_POSITIONS` 上方的晶格矢量（Lattice Vectors）必须与 System AB **完全一致**。
2.  **移除原子**：直接删除（或注释掉）属于 B 部分的原子行。
3.  **禁止移动**：绝对不要改变剩余 A 原子的坐标数值。

**System AB 的 STRU (示例):**
```text
ATOMIC_POSITIONS
Cartesian
Fe 0.00 0.00 0.00  # A部分
O  1.50 0.00 0.00  # B部分
```

**System A 的 STRU (修改后):**
```text
ATOMIC_POSITIONS
Cartesian
Fe 0.00 0.00 0.00  # 保留 A
# O  1.50 0.00 0.00  <-- 删除 B
```

### 2.2.2 锁定 FFT 网格（Grid Consistency）
这是最容易被忽略的隐形陷阱。
虽然理论上只要晶胞大小（Cell）和截断能（`ecutwfc`）不变，FFT 网格应该是相同的。但在某些极端情况下（或不同版本的自动并行策略下），软件自动生成的网格可能会有微小差异（例如从 108 变为 110）。

为了万无一失，建议在 System A 和 System B 的 `INPUT` 文件中**手动固定 FFT 网格维度**，数值来源于 2.1.3 节中记录的 System AB 的日志。

**System A/B 的 INPUT 追加参数：**
```bash
# 手动固定 FFT 网格，确保与 System AB 完全一致
# 注意：参数名可能随版本略有不同，请查阅文档确认是 nx/ny/nz 或 nr1/nr2/nr3
nx              108 
ny              108
nz              144

symmetry        -1    # 同样必须保持关闭
out_chg         1     # 同样需要输出电荷
```

### 2.2.3 LCAO 基组的特殊考量（Ghost Atom）
如果你使用的是 `basis_type lcao`：
*   **直接移除原子**：意味着同时也移除了该原子所携带的原子轨道基组。这在物理上对应的是“将 B 移至无穷远”。这是计算成键电荷密度的标准做法。
*   **Ghost Atom（鬼原子）**：如果你需要极高精度的基组重叠误差（BSSE）校正，或者分析的是特定轨道杂化，可能需要将 B 原子设为“Ghost Atom”（保留基组但不含核电荷与电子）。
    *   *本教程采用标准做法（直接移除原子），适合绝大多数界面与成键分析。*

---

## 2.3 常见陷阱自检（Self-Check）

在进行最终的减法处理之前（将在第三章介绍），请务必执行以下“三查”：

1.  **查收敛性**：
    检查 `OUT.A/running_scf.log` 和 `OUT.B/running_scf.log`，确认 `TOTAL-FORCE` 和 `ETOT` 均已达到收敛标准。孤立体系有时比完整体系更难收敛（因为存在悬挂键），若不收敛，需调整 `mixing_beta`。

2.  **查网格维度（Grid Match）**：
    使用 `grep` 命令快速对比三个计算的网格是否严格一致。
    ```bash
    grep "FFT grid" OUT.*/running_scf.log
    ```
    **预期输出（必须完全相同）：**
    ```text
    OUT.AB/running_scf.log: FFT grid for charge/potential:  108 * 108 * 144
    OUT.A/running_scf.log:  FFT grid for charge/potential:  108 * 108 * 144
    OUT.B/running_scf.log:  FFT grid for charge/potential:  108 * 108 * 144
    ```

3.  **查输出文件**：
    确认目录下生成了对应的电荷密度文件（通常为 `SPIN1_CHG.cube` 或 `density.cube`，具体取决于版本和设置）。

---

**下一步预告**：
现在你已经手握三份在空间上完美对齐的“数字地图”（$\rho_{AB}, \rho_A, \rho_B$）。在第三章中，我们将使用 Python 脚本或后处理工具，对这些地图进行像素级的减法运算，揭示出隐藏在原子间隙中的电子转移秘密。