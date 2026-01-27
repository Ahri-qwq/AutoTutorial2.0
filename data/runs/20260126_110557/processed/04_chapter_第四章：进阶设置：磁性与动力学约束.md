# 第四章：进阶设置：磁性与动力学约束

你好！我是你的 ABACUS 导师。在上一章中，我们已经成功构建了最基础的晶体结构。但在实际的科研战场上，我们面对的往往不是完美的、静止的、非磁性的理想晶体。

本章我们将深入 `STRU` 文件的“深水区”。我们将学习如何处理**磁性体系**（这是 ABACUS 与 VASP 逻辑差异最大的地方之一），以及如何在结构优化中施加**动力学约束**（例如固定表面模型的底层原子）。

---

# 第四章：进阶设置：磁性与动力学约束

## 4.1 初始磁矩的设置 (Magnetic Moments)

在处理过渡金属氧化物、铁磁/反铁磁材料或自旋极化体系时，正确设置初始磁矩至关重要。这不仅影响计算的收敛性，更直接决定了你得到的是铁磁态、反铁磁态还是非磁态。

### 4.1.1 核心差异：告别 INCAR，拥抱 STRU

如果你是 VASP 用户，请务必注意：**ABACUS 不在 `INPUT` 文件中设置每个原子的磁矩。**
你不需要像在 `INCAR` 中那样编写冗长的 `MAGMOM` 数组。在 ABACUS 中，磁矩是原子的内禀属性，直接定义在 `STRU` 文件中。

### 4.1.2 语法详解

在 `ATOMIC_POSITIONS` 块中，每一行原子坐标的定义实际上包含了更多信息。标准格式如下：

```text
<X坐标> <Y坐标> <Z坐标> <move_x> <move_y> <move_z> mag <初始磁矩值>
```

*   **`mag`**: 这是一个关键词，必须显式写出。
*   **`<初始磁矩值>`**: 紧跟在 `mag` 之后，单位通常为玻尔磁子 ($\mu_B$)。

### 4.1.3 实战案例：反铁磁氧化镍 (NiO)

假设我们要计算沿 [111] 方向反铁磁排列的 NiO。我们需要明确指定哪些 Ni 原子自旋向上，哪些向下。

**STRU 文件片段示例：**

```text
ATOMIC_POSITIONS
Direct

Ni
1.0   # 磁矩缩放因子（通常设为0，由mag关键字覆盖，或设为1.0）
2     # 原子数量
#  X    Y    Z    mx my mz  mag  value
0.00 0.00 0.00    1  1  1   mag  2.0   # Ni1: 自旋向上
0.50 0.50 0.50    1  1  1   mag -2.0   # Ni2: 自旋向下

O
0.0   # 氧原子通常无磁性
2
0.25 0.25 0.25    1  1  1   mag  0.0
0.75 0.75 0.75    1  1  1   mag  0.0
```

> **教授提示**：
> 1. **INPUT 配合**：仅在 `STRU` 中写 `mag` 是不够的，你必须在 `INPUT` 文件中设置 `nspin 2`（开启自旋极化计算）。否则，ABACUS 会忽略 `mag` 设置，强制进行非磁计算。
> 2. **收敛技巧**：初始磁矩建议设置得比预期值稍大一些（例如预期 1.5 $\mu_B$，初始设为 2.0 or 3.0），这有助于打破对称性，引导体系落入正确的磁性基态。

---

## 4.2 结构优化中的原子约束 (Constraints)

在进行表面吸附、界面弛豫或寻找过渡态（NEB）计算时，我们经常需要固定一部分原子，只允许另一部分原子移动。

### 4.2.1 坐标后的“0”与“1”

在 `ATOMIC_POSITIONS` 中，紧跟在 X, Y, Z 坐标后的三个整数，分别代表该原子在 **x, y, z 三个方向上的移动自由度**。

*   `1`: **允许移动** (True)
*   `0`: **固定不动** (False)

### 4.2.2 实战案例：金属表面吸附

假设我们有一个 4 层的 Al(111) 表面，我们希望固定最底下的两层原子以模拟体相环境，只弛豫上面两层和吸附分子。

**STRU 文件片段示例：**

```text
ATOMIC_POSITIONS
Direct

Al
0.0
4
# 底层原子（固定）
0.00 0.00 0.00    0  0  0   mag 0.0  # x,y,z 全固定
0.50 0.50 0.00    0  0  0   mag 0.0

# 表层原子（弛豫）
0.00 0.00 0.25    1  1  1   mag 0.0  # x,y,z 全自由
0.50 0.50 0.25    1  1  1   mag 0.0

# 仅允许 Z 方向弛豫的特殊情况（例如层间距优化）
# 0.50 0.50 0.50  0  0  1   mag 0.0 
```

---

## 4.3 高效工具：使用 ASE 生成 STRU

手动编辑几十上百个原子的 `0` 和 `1` 既痛苦又容易出错。作为开发者，我强烈推荐使用 Python 的 ASE 库来自动化这一过程。

### 4.3.1 关键警告：ASE 的安装

<div style="background-color: #ffe6e6; padding: 15px; border-left: 5px solid #ff0000; margin-bottom: 20px;">
<strong>🚨 红色警报：不要直接 pip install ase！</strong>
<br><br>
官方的 ASE 库（PyPI版）目前对 ABACUS 的支持非常有限且过时。你必须安装带有 ABACUS 接口的分支版本。
<br><br>
<strong>正确安装命令：</strong>
<pre><code>git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .</code></pre>
</div>

### 4.3.2 Python 脚本模板

以下脚本展示了如何读取 CIF 文件，设置磁矩，并生成包含赝势（PP）和轨道（Orbital）信息的完整 `STRU` 文件。

```python
from ase.io import read, write
import os

# 1. 读取结构 (例如从 Materials Project 下载的 CIF)
atoms = read("Fe3O4.cif")

# 2. 设置磁矩 (对应原子顺序)
# 假设前3个是Fe，后4个是O。这里演示简单的铁磁设置
# 实际操作中需根据原子索引精确设置
mag_list = [4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0]
atoms.set_initial_magnetic_moments(mag_list)

# 3. 设置约束 (固定 Z < 5.0 Angstrom 的原子)
# constraint 需使用 ASE 的 FixAtoms 类，但在输出 STRU 时，
# ase-abacus 接口会自动识别 atoms 对象的 constraint 属性吗？
# 实际上，ase-abacus 目前主要通过直接修改 atoms 对象的数组来处理，
# 但更稳妥的方式是直接在生成后检查 STRU，或者使用专门的 constraint 设置函数。
# 这里展示最通用的方法：直接生成，后续微调。
# (注：高级用户可直接操作 atoms.constraints)

# 4. 指定赝势和轨道 (关键步骤！)
# 必须提供字典，键为元素符号，值为文件名
pp_dict = {
    "Fe": "Fe_ONCV_PBE-1.0.upf",
    "O":  "O_ONCV_PBE-1.0.upf"
}

# 如果是 LCAO 计算，必须指定 basis；如果是 PW 计算，basis 设为 None 或空字典
orb_dict = {
    "Fe": "Fe_gga_8au_100Ry_4s2p2d.orb",
    "O":  "O_gga_7au_100Ry_2s2p1d.orb"
}

# 5. 输出 STRU
write("STRU", atoms, format="abacus", pp=pp_dict, basis=orb_dict)

print("STRU file generated successfully!")
```

### 4.3.3 懒人神器：ATOMKIT

如果你不想写 Python 代码，可以使用 **ATOMKIT**。这是一个强大的命令行工具，支持 VASP POSCAR 与 ABACUS STRU 的互转。

*   **转换命令**: 交互式运行 `atomkit`，选择 `113` (CIF) 或 `175` (VASP) 导入，然后选择 `101` (ABACUS STRU) 导出。
*   **优点**: 自动处理晶格格式，速度极快。

---

<!-- APPENDIX_START -->
## 附录：常见问题与排错指南 (Troubleshooting)

在配置 `STRU` 文件时，新手最容易在以下三个地方“翻车”。

### 1. 晶格常数的单位陷阱 (The 1.8897 Confusion)

这是 ABACUS 论坛提问率 Top 1 的问题。

*   **现象**: 结构看起来是对的，但计算出的体积巨大或极小，能量异常。
*   **原因**: 混淆了 Bohr 和 Angstrom。
*   **原理**: ABACUS 内部计算单位是 Bohr (1 Bohr $\approx$ 0.529177 Angstrom)。
    *   `1 Angstrom / 0.529177 ≈ 1.889726125`

**最佳实践规则**:

| 你的 LATTICE_VECTORS 单位 | LATTICE_CONSTANT 设置值 | 解释 |
| :--- | :--- | :--- |
| **Angstrom (埃)** (推荐) | `1.8897261258369282` | ABACUS 会把你的矢量乘以这个常数，转为内部的 Bohr。 |
| **Bohr (波尔)** | `1.0` | 你的矢量已经是 Bohr 了，不需要缩放。 |

> **风险提示**: 虽然存在 `Cartesian_angstrom` 选项，但由于部分后处理工具支持不佳，**不推荐**在生产环境中使用，建议坚持使用 `LATTICE_CONSTANT` + `Direct` 的组合。

### 2. LCAO 计算崩溃：缺少 NUMERICAL_ORBITAL

*   **报错**: 程序启动即崩溃，提示读取轨道错误。
*   **检查**:
    *   如果是平面波 (PW) 计算 (`basis_type pw`)，`STRU` 中**不需要** `NUMERICAL_ORBITAL` 块。
    *   如果是原子轨道 (LCAO) 计算 (`basis_type lcao`)，`STRU` 中**必须**包含 `NUMERICAL_ORBITAL` 块，且轨道文件的顺序必须与 `ATOMIC_SPECIES` 中的元素顺序严格一致。

### 3. 文件路径未找到

*   **报错**: `file not found: element.upf`
*   **解决**: 
    1. 绝对路径：在 `STRU` 中写死 `/home/user/pp/Fe.upf` (不推荐，迁移性差)。
    2. 环境变量 (推荐)：
       ```bash
       export ABACUS_PP_PATH=/path/to/your/pseudopotentials
       export ABACUS_ORBITAL_PATH=/path/to/your/orbitals
       ```
       设置好环境变量后，`STRU` 中只需写文件名 `Fe.upf` 即可。

---
**下一章预告**：
搞定了结构和磁性，我们将进入第五章：**电子结构计算核心参数**。我们将剖析 `INPUT` 文件中的 `ecutwfc`、`kspacing` 以及基组选择的艺术。