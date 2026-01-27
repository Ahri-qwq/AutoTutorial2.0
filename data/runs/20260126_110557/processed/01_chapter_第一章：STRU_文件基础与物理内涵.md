# 第一章：STRU 文件基础与物理内涵

欢迎来到 ABACUS 实战教程。作为一名深耕计算材料学多年的开发者，我深知从 VASP 或 Quantum Espresso 迁移到 ABACUS 时，最大的认知阻力往往来自于输入文件的逻辑差异。

在 ABACUS 中，`STRU` 文件不仅仅是几何结构的容器，它还承载了基组定义（LCAO 模式）和赝势映射。本章将带你深入理解 STRU 文件的物理内涵，避开新手最容易踩的“单位制”和“基组匹配”这两个深坑。

---

## ⚠️ 课前准备：避坑指南（关键！）

在开始编写文件之前，我们通常需要借助工具进行格式转换。这里有一个**极其关键**的注意事项，请务必仔细阅读：

### 关于 ASE (Atomic Simulation Environment) 的安装
很多教程会告诉你直接运行 `pip install ase`，但在使用 ABACUS 时，**请不要这样做**！官方源中的 ASE 尚未完全包含最新的 ABACUS 接口支持。

**正确做法**是安装带有 ABACUS 接口的分支版本：

```bash
# 卸载可能存在的官方版本
pip uninstall ase

# 安装支持 ABACUS 的分支版本
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```
*只有这样，你才能在 Python 中正确使用 `format='abacus'` 进行读写操作。*

---

## 1.1 STRU 文件在计算流程中的角色

在 DFT 计算的“三驾马车”（控制参数、结构、K点）中，`STRU` 文件的地位等同于 VASP 的 `POSCAR`，但它包含的信息量略大。

### 物理意义映射
| 概念 | VASP (POSCAR + POTCAR) | ABACUS (STRU) | 备注 |
| :--- | :--- | :--- | :--- |
| **晶格矢量** | POSCAR header | `LATTICE_VECTORS` | 需配合缩放系数 |
| **原子位置** | POSCAR body | `ATOMIC_POSITIONS` | 支持分数/笛卡尔坐标 |
| **元素种类** | POSCAR (元素行) | `ATOMIC_SPECIES` | 显式定义元素与赝势文件的对应关系 |
| **赝势文件** | POTCAR (独立文件) | `ATOMIC_SPECIES` (文件名) | STRU 中只写文件名，路径在环境变量或 INPUT 中指定 |
| **磁矩设置** | INCAR (MAGMOM) | `ATOMIC_POSITIONS` | **差异点**：ABACUS 是逐个原子设置 |

**核心逻辑**：STRU 负责描述**实空间**（Real Space）的所有物理属性。与 `INPUT`（控制计算流程的哈密顿量参数）和 `KPT`（倒空间采样）共同构成了计算的初始态。

---

## 1.2 核心语法与单位制详解：破解 1.8897 之谜

ABACUS 的内核是基于原子单位制（Atomic Units, a.u.）开发的，这意味着其底层长度单位是 **Bohr（波尔半径）**。然而，材料学界通用的单位是 **Angstrom（埃）**。

为了解决这个矛盾，`STRU` 文件设计了 `LATTICE_CONSTANT` 这一参数。

### 标准 STRU 文件结构解析

```abacus
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf  # 元素符号 质量 赝势文件名

LATTICE_CONSTANT
1.8897261258369282  # 关键缩放因子

LATTICE_VECTORS
2.715 2.715 0.000   # 矢量 v1
2.715 0.000 2.715   # 矢量 v2
0.000 2.715 2.715   # 矢量 v3

ATOMIC_POSITIONS
Direct              # 坐标类型：Direct(分数坐标) 或 Cartesian(笛卡尔坐标)
Si                  # 元素标签（必须与 ATOMIC_SPECIES 对应）
0.0                 # 初始磁矩（低优先级，通常设为0，在原子行具体设置）
2                   # 原子数量
0.00 0.00 0.00 1 1 1 mag 0.0  # x y z move_x move_y move_z mag moment
0.25 0.25 0.25 1 1 1 mag 0.0
```

### 深度解析：LATTICE_CONSTANT 的机制

这是新手最容易困惑的地方。ABACUS 读取晶格矢量的逻辑是：
$$ \vec{V}_{actual} (Bohr) = \text{LATTICE\_CONSTANT} \times \vec{V}_{input} $$

因此，我们有两种常见的设置策略：

1.  **推荐策略（埃单位输入）**：
    *   如果你习惯用 **Angstrom** 思考（例如从 CIF 文件读取的数值），你需要将 `LATTICE_CONSTANT` 设为 **1.889726...**（即 $1 \text{ Angstrom} \approx 1.8897 \text{ Bohr}$）。
    *   此时，`LATTICE_VECTORS` 中的数值填埃单位的数值。
    *   *物理结果*：程序内部会自动将其乘起来，得到正确的 Bohr 单位数值。

2.  **原子单位策略（Bohr 单位输入）**：
    *   如果你直接填入 **Bohr** 单位的数值，则将 `LATTICE_CONSTANT` 设为 **1.0**。

> **风险提示**：虽然存在 `Cartesian_angstrom` 这一坐标类型选项，允许直接使用埃单位而不依赖缩放因子，但考虑到该功能在部分后处理工具（如某些可视化软件）中的兼容性较差，我**强烈建议**坚持使用标准的 `LATTICE_CONSTANT` + `Direct/Cartesian` 组合。

---

## 1.3 基组依赖性：LCAO 与 PW 的区别

ABACUS 的独特之处在于同时支持平面波（Plane Wave, PW）和数值原子轨道（Numerical Atomic Orbitals, LCAO）两种基组。这直接导致了 `STRU` 写法的不同。

### 关键差异点：NUMERICAL_ORBITAL 模块

*   **平面波计算 (PW)**：
    *   不需要轨道文件。
    *   `STRU` 文件中**不要**包含 `NUMERICAL_ORBITAL` 模块。

*   **原子轨道计算 (LCAO)**：
    *   必须包含 `NUMERICAL_ORBITAL` 模块。
    *   **严格顺序规则**：轨道文件的排列顺序必须与 `ATOMIC_SPECIES` 中的元素顺序**完全一致**。

**LCAO 模式下的 STRU 片段示例**：

```abacus
ATOMIC_SPECIES
Ga 69.723 Ga_ONCV_PBE-1.0.upf
As 74.921 As_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Ga_gga_8au_100Ry_4s4p1d.orb  # 对应第一个元素 Ga
As_gga_8au_100Ry_4s4p1d.orb  # 对应第二个元素 As
```

---

## 1.4 进阶设置：磁矩与原子弛豫

在 VASP 中，我们习惯在 `INCAR` 中通过 `MAGMOM` 数组来设置磁矩。但在 ABACUS 中，这一物理量回归到了它所属的几何空间——`STRU` 文件。

### 原子位置行的详细定义
每一行原子位置包含以下信息：
```text
coord_x  coord_y  coord_z  m_x  m_y  m_z  mag  moment
```

*   `coord_x/y/z`: 坐标数值。
*   `m_x/y/z`: 移动限制（1=可动，0=固定）。这对应 VASP POSCAR 中的 Selective Dynamics。
*   `mag`: 关键词，标志着后面跟着的是磁矩。
*   `moment`: 初始磁矩值（玻尔磁子 $\mu_B$）。

**示例**：设置一个反铁磁构型（第一个原子自旋向上，第二个向下）：
```abacus
ATOMIC_POSITIONS
Direct
Fe
0.0
2
0.0 0.0 0.0 1 1 1 mag  2.0
0.5 0.5 0.5 1 1 1 mag -2.0
```

---

## 1.5 实战工具箱：如何快速生成 STRU

手动编写 STRU 文件容易出错，推荐使用以下两种“黑盒”方案。

### 方案 A：Python 脚本（最灵活）
利用我们开头安装的 ASE-ABACUS 接口。这是转换 CIF 到 STRU 的标准范式，**请务必同时指定 `pp` (赝势) 和 `basis` (轨道) 字典**，否则生成的 STRU 会缺失关键信息。

```python
from ase.io import read, write

# 1. 读取结构 (支持 cif, poscar, xsd 等)
atoms = read("structure.cif")

# 2. 定义赝势和轨道映射 (Key必须是元素符号)
# 注意：这里只写文件名，不需要写绝对路径
pp_dict = {
    "Si": "Si_ONCV_PBE-1.0.upf",
    "O":  "O_ONCV_PBE-1.0.upf"
}

basis_dict = {
    "Si": "Si_gga_8au_100Ry_4s4p1d.orb",
    "O":  "O_gga_8au_100Ry_2s2p1d.orb"
}

# 3. 输出 STRU 文件
# 如果是做平面波计算，去掉 basis=basis_dict 参数即可
write("STRU", atoms, format="abacus", pp=pp_dict, basis=basis_dict)

print("STRU file generated successfully!")
```

### 方案 B：ATOMKIT（最快速）
如果你不想写代码，可以使用 VASPKIT 团队开发的 `ATOMKIT`。这是一个交互式的命令行工具。

1.  运行 `atomkit`。
2.  选择 `101` (ABACUS STRU) 作为导出格式。
3.  输入导入文件名（如 `POSCAR` 或 `test.cif`）。
4.  工具会自动生成标准的 STRU 文件（默认使用 1.8897 缩放因子）。

---

**本章总结**：
理解 `STRU` 文件，关键在于掌握 **1.8897 的单位换算**以及 **LCAO 模式下的轨道对应关系**。当你能熟练地在 STRU 中定义磁矩和固定原子时，你就已经迈出了掌握 ABACUS 的第一步。下一章，我们将进入 `INPUT` 文件，探讨如何控制量子力学的演化过程。