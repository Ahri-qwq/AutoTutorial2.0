# 第四章：任务类型切换与参数定制

## 第四章：任务类型切换与参数定制

在上一章中，我们介绍了如何快速生成默认的 SCF（自洽场）计算任务。然而在实际科研中，我们往往需要对晶体结构进行几何优化（Relaxation），或者针对特定体系（如磁性材料、强关联体系）调整计算参数。本章将深入介绍如何利用 `abacustest` 的 `--jtype`、`--input` 和 `--kpt` 等参数，灵活地切换任务类型并定制输入文件。

### 4.1 结构优化任务（Relax/Cell-Relax）

默认情况下，`abacustest` 生成的是静态计算（SCF）的输入文件。如果需要优化原子位置或晶胞参数，可以通过 `--jtype` 参数指定任务类型。常见的类型包括 `relax`（仅优化原子位置）和 `cell-relax`（同时优化原子位置和晶胞形状/体积）。

#### 案例：MgO 的变胞优化

假设我们需要对 MgO 的晶胞结构进行全优化。我们可以使用 `--jtype cell-relax` 参数，并指定一个新的文件夹命名规则以便区分：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax
```

执行该命令后，生成的 `INPUT` 文件会自动调整为适合变胞优化的配置。以下是生成的关键参数解析：

```python
INPUT_PARAMETERS
calculation     cell-relax    # 计算类型切换为变胞优化
# ... (省略基础参数) ...
cal_force       1             # 开启力计算
cal_stress      1             # 开启应力计算（变胞优化必须）
relax_method    cg            # 优化算法，默认使用共轭梯度法(CG)
relax_nmax      60            # 最大离子步数
force_thr_ev    0.01          # 力收敛标准，单位 eV/A
stress_thr      0.5           # 应力收敛标准，单位 kbar
fixed_axes      None          # 不固定任何轴，允许全自由度弛豫
```

**参数详解：**
*   **calculation**: 自动被设置为 `cell-relax`。如果指定 `--jtype relax`，此处则为 `relax`，且会自动移除 `cal_stress` 和 `stress_thr` 参数。
*   **relax_method**: 对于 `cell-relax`，工具默认推荐使用 `cg` (Conjugate Gradient) 算法，这在 ABACUS 的变胞优化中通常表现稳健。
*   **收敛标准**: 默认设置了较为合理的收敛阈值（力 0.01 eV/A，应力 0.5 kbar）。用户若需更高精度的结果，可参照下一节的方法进行覆盖。

### 4.2 自定义 INPUT 模板与 K 点设置

虽然 `abacustest` 提供了一套通用的默认参数，但在处理复杂体系时，我们经常需要手动干预。例如，对于磁性体系或金属体系，可能需要更精细的展宽（Smearing）设置或更密集的 K 点网格。

此时，我们可以通过准备一个 `INPUT_template` 文件和使用 `--kpt` 参数来实现定制化。

#### 案例：Fe2O3 的高精度磁性计算

以 Fe2O3 为例，假设我们需要进行如下定制：
1.  **收敛性调整**：将 `mixing_beta` 降低至 0.2 以防止电荷震荡。
2.  **电子温度**：将 `smearing_sigma` 设置为更小的 0.001 Ry。
3.  **K点设置**：不再使用默认的 `kspacing`，而是强制指定 5x5x5 的 Gamma-centered 网格。

**步骤 1：准备模板文件**

在当前目录下创建一个名为 `INPUT_template` 的文件，写入需要覆盖或新增的参数：

```python
INPUT_PARAMETERS
smearing_sigma     0.001
mixing_beta        0.2
```

**步骤 2：执行生成命令**

在命令中加入 `--input INPUT_template` 和 `--kpt 5 5 5`，同时结合磁性（`--nspin 2`）和 DFT+U 的设置：

```bash
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao \
    --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0 \
    --input INPUT_template --kpt 5 5 5
```

**步骤 3：结果验证**

检查生成的 `Fe2O3/INPUT` 文件，你会发现：
*   `smearing_sigma` 和 `mixing_beta` 已被更新为模板中的值。
*   **关键变化**：原有的 `kspacing` 参数被自动移除。这是因为我们显式指定了 K 点网格，ABACUS 将读取生成的 `KPT` 文件，而不是根据间距自动生成。

生成的 `KPT` 文件内容将如下所示：

```text
K_POINTS
0
Gamma
5 5 5 0 0 0
```

通过这种方式，我们既利用了自动化工具处理了繁琐的结构转换（CIF 到 STRU）和赝势配置，又保留了对关键物理参数的完全控制权。在下一章中，我们将探讨如何批量处理多个结构文件，以应对高通量计算的需求。