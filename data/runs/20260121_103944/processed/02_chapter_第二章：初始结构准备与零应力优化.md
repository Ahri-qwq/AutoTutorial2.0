# 第二章：初始结构准备与零应力优化

在计算材料的弹性性质之前，最关键的一步是获得一个处于“零应力”状态的基态结构。弹性常数（Elastic Constants）描述的是材料在平衡位置附近对应变的响应，如果初始结构本身含有残余应力（Residual Stress），根据广义胡克定律，这将直接导致计算出的应力-应变关系产生偏差，甚至得到非物理的结果。

本章将指导你如何按照 IEEE 标准准备晶体结构，并使用 `abacustest` 工具流自动化生成 ABACUS 任务，完成高精度的变胞优化（Cell Relaxation）。

---

## 2.1 晶体取向与惯用胞选择

在开始计算之前，我们必须明确一个概念：**弹性张量 $C_{ij}$ 的数值高度依赖于坐标系的选择。**

### 为什么选择惯用胞（Conventional Cell）？
对于硅（Si）、铜（Cu）等常见的立方晶系材料，其原胞（Primitive Cell）通常是非正交的（例如 Si 的面心立方原胞包含 2 个原子，晶格矢量互成 60 度）。然而，文献中汇报的弹性模量（如 $C_{11}, C_{12}, C_{44}$）通常是基于**立方惯用胞**定义的。

*   **原胞计算的后果**：如果你直接使用原胞计算，得到的弹性张量将是基于原胞基矢方向的，这会导致矩阵元极其复杂，难以与实验值直接对比。
*   **惯用胞的优势**：使用正方体形式的惯用胞（对于 Si 为 8 个原子），并将晶格矢量 $\vec{a}, \vec{b}, \vec{c}$ 分别与笛卡尔坐标系的 $x, y, z$ 轴对齐。这种取向符合 IEEE 标准（IEEE 176-1987），能确保计算出的 $C_{ij}$ 矩阵具有标准的对称性形式（如立方晶系仅有 3 个独立分量）。

> **专家提示**：在准备 `STRU` 或 `CIF` 文件时，请务必检查晶格矢量是否已旋转至标准笛卡尔方向。如果晶体被随意旋转，你计算出的将是旋转后的张量 $C'_{ijkl}$。

---

## 2.2 生成初始结构与优化输入文件

我们将使用 Python 的 `ase` 库生成标准的 Si 惯用胞，并利用 `abacustest` 自动生成 ABACUS 的输入文件。

### 2.2.1 第一步：生成标准结构文件
创建一个名为 `generate_structure.py` 的 Python 脚本，内容如下：

```python
from ase.build import bulk

# 生成 Si 的金刚石结构
# cubic=True 确保生成的是正方体惯用胞（8个原子），而非原胞
# 此时晶格矢量默认与 Cartesian 坐标轴对齐
si = bulk('Si', 'diamond', a=5.43, cubic=True)

# 保存为 CIF 格式
si.write('Si.cif')
print("Si.cif generated successfully with cubic orientation.")
```

运行该脚本后，当前目录下将生成 `Si.cif` 文件。

### 2.2.2 第二步：使用 `abacustest` 生成任务文件夹
`abacustest` 是 ABACUS 的官方辅助工具，能够根据结构文件自动匹配赝势（PP）和轨道（Orbital），并生成规范的 `INPUT` 文件。

请确保你的环境变量 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH` 已正确指向赝势库和轨道库目录。

在终端执行以下命令：

```bash
abacustest submit -f Si.cif --ftype cif --jtype cell-relax --lcao --folder-syntax Si
```

**参数详解**：
*   `-f Si.cif`: 指定输入的结构文件。
*   `--ftype cif`: 指定文件格式为 CIF。
*   `--jtype cell-relax`: **关键参数**。指定任务类型为“变胞优化”，即同时优化原子位置和晶胞尺寸，以消除应力。
*   `--lcao`: 指定使用 LCAO（线性组合原子轨道）基组，计算效率更高。
*   `--folder-syntax Si`: 指定生成的任务文件夹名称为 `Si`。

**警告**：
> **请注意，`abacustest` 的生成命令通常具有覆盖性。** 如果你重复执行上述命令，它可能会删除并重新生成 `Si/` 文件夹。请务必在操作前确认数据备份，防止误删已有的计算结果。

执行完毕后，你会看到一个名为 `Si/` 的文件夹，其中包含：
*   `INPUT`: 控制参数文件
*   `STRU`: 由 CIF 转换而来的 ABACUS 结构文件
*   `*.upf` / `*.orb`: 自动链接的赝势和轨道文件

---

## 2.3 执行结构优化（Cell-Relax）

为了获得准确的弹性常数，初始结构的残余应力必须收敛至忽略不计（通常建议小于 1 kBar）。

### 2.3.1 检查 INPUT 文件
进入 `Si/` 文件夹，检查自动生成的 `INPUT` 文件。重点关注以下参数（如果未生成，需手动调整）：

```plaintext
INPUT_PARAMETERS
# ... (省略系统参数)

calculation     cell-relax  # 必须：进行变胞优化
basis_type      lcao        # 基组类型

# 收敛标准
force_thr_ev    0.001       # 力收敛阈值 (eV/Angstrom)
stress_thr      1           # 应力收敛阈值 (kBar)，弹性计算要求该值尽可能小
relax_nmax      100         # 最大离子步数

# 输出设置
out_stru        1           # 输出优化后的结构文件
cal_stress      1           # 显式开启应力计算（cell-relax模式下通常默认开启）
```

### 2.3.2 提交计算与结果验证
使用提交脚本（如 `sbatch`）或直接运行 ABACUS 可执行程序：

```bash
# 示例运行命令（需根据实际环境调整）
mpirun -np 8 abacus | tee output.log
```

计算完成后，检查输出目录（通常为 `OUT.Si/`）下的 `running_cell-relax.log` 或主日志文件。

1.  **验证收敛性**：确认任务正常结束，且最后一步的 `Max Force` 和 `Max Stress` 均低于阈值。
2.  **获取基态结构**：优化后的结构保存在 `OUT.Si/STRU_ION_D`（最后一步的结构）。

### 2.3.3 为什么必须做 Cell-Relax？
在计算弹性常数时，我们实际上是在计算能量（或应力）对应变的二阶导数。
*   **非零应力的影响**：如果初始结构存在残余应力 $\sigma_0$，则应力-应变关系变为 $\sigma = \sigma_0 + C \cdot \epsilon$。虽然可以通过拟合扣除截距，但非平衡态下的 $C$ 值本身也会发生物理偏移。
*   **原子弛豫（Relaxed-ion）**：在后续施加应变计算弹性常数时，通常建议允许原子在变形后的晶胞内移动（即 `relax`），这对应于**松弛离子（Relaxed-ion）**弹性常数，比“固定离子（Clamped-ion）”结果更符合宏观物理实测。因此，从一个完全松弛的基态出发是至关重要的。

### 下一步准备
优化完成后，你需要使用优化后的结构（`STRU_ION_D`）作为后续弹性模量计算的起点。

**操作流：**
1.  将 `OUT.Si/STRU_ION_D` 复制出来，重命名为 `STRU_optimized`。
2.  确认该结构的晶格常数（对于 Si 应接近 5.43 Å）和原子位置正确。
3.  此结构将作为第三章批量生成形变结构（Deformation）的母本。

> **注意**：在本教程的案例流程中，我们将在下一章直接使用 `abacustest` 的弹性工作流，它会自动处理后续的形变生成，但前提是你必须提供这个优化好的结构。