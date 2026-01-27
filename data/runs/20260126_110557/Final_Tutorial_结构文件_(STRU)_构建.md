# ABACUS 结构文件 (STRU) 指南

## 前言

欢迎踏入 ABACUS 的第一性原理计算世界。作为国产开源电子结构计算软件的佼佼者，ABACUS 凭借其独特的数值原子轨道（LCAO）与平面波（PW）双基组优势，在处理大规模体系和高精度模拟中展现了卓越的性能。然而，对于习惯了 VASP 或 Quantum Espresso 逻辑的研究者而言，ABACUS 的 `STRU` 文件结构、单位制以及基组匹配往往是入门时的第一道“摩擦力”。

**学习路线图**
本书聚焦结构与前处理。我们将从**第一章**夯实基础，深入剖析 `STRU` 文件的物理内涵与单位制避坑指南；在**第二章**引入 ASE (Atomic Simulation Environment)，教你如何利用 Python 脚本实现从 POSCAR 转换到批量建模的自动化流程；**第三章**将为你装备 ATOMKIT 与 abacustest 两个高效工具箱，实现从“手动操作”到“工程化生产”的跨越；最后，在**第四章**攻克磁性设置与动力学约束等进阶课题，助你应对真实的科研战场。

**知识体系定位**
在整个 ABACUS/DFT 知识体系中，本书定位于“输入构建与预处理”这一核心支柱。它是连接原始材料数据与最终物理特性分析的桥梁。只有掌握了精准的结构描述与参数定义，后续的电子自洽迭代（SCF）、能带结构计算及分子动力学模拟才具有物理意义。

**前置知识**
在开始本书的学习前，我们建议你已具备以下基础：
1. 密度泛函理论（DFT）的基本物理概念。
2. 基础的 Linux 命令行操作经验。
3. 简单的 Python 语法基础（有助于理解 ASE 脚本）。

---

# 第一章：STRU 文件基础与物理内涵

从 VASP 或 Quantum Espresso 迁移到 ABACUS 时，最大的认知阻力往往来自于输入文件的逻辑差异。
在 ABACUS 中，`STRU` 文件不仅仅是几何结构的容器，它还承载了基组定义（LCAO 模式）和赝势映射。本章将带你深入理解 STRU 文件的物理内涵，避开新手最容易踩的“单位制”和“基组匹配”这两个深坑。

---

## 避坑指南（关键！）

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

> **风险提示**：虽然存在 `Cartesian_angstrom` 这一坐标类型选项，允许直接使用埃单位而不依赖缩放因子，但考虑到该功能在部分后处理工具（如某些可视化软件）中的兼容性较差，在某些旧工具链里可能不如 `LATTICE_CONSTANT` + `Direct/Cartesian` 组合常用。

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
理解 `STRU` 文件，关键在于掌握 **1.8897 的单位换算**以及 **LCAO 模式下的轨道对应关系**。当你能熟练地在 STRU 中定义磁矩和固定原子时，你就已经迈出了掌握 ABACUS 的第一步。

---

# 第二章：基于 ASE 的 Python 脚本化构建

在计算材料学的实际工作中，手动编写结构文件不仅效率低下，而且极易出错。推荐使用 **ASE (Atomic Simulation Environment)** 作为构建和操作原子结构的首选工具。

ASE 是 Python 生态中最强大的原子模拟接口库。本章将带你掌握如何利用 ASE 的灵活性，完成从零建模、格式转换（如 POSCAR 转 STRU）以及批量生成输入文件的全流程。


## 2.1 环境配置与 ASE-ABACUS 接口安装

### 2.1.1 关键警告：安装源的选择
这也是非常容易踩的一个坑。虽然 PyPI 官方源提供了 `ase` 包（即 `pip install ase`），但官方版本对 ABACUS 的支持非常有限且更新滞后。

为了获得完整的 STRU 文件读写支持（包括磁矩设置、多元素混合、LCAO 轨道映射等），**必须安装带有 ABACUS 接口的分支版本**。

**推荐安装步骤：**

```bash
# 1. 如果之前安装过官方 ase，建议先卸载以避免冲突
pip uninstall ase

# 2. 克隆包含 ABACUS 接口的 ASE 分支
git clone https://gitlab.com/1041176461/ase-abacus.git

# 3. 进入目录并安装
cd ase-abacus
pip install .
```

### 2.1.2 环境变量配置
ASE 在生成输入文件时，需要知道你的赝势（PP）和数值轨道（ORB）文件存放在哪里，以便进行路径检查或自动关联。请在你的 `~/.bashrc` 或 `~/.zshrc` 中添加以下环境变量：

```bash
# 假设你的赝势和轨道文件分别存放在 ~/pseudopotentials 和 ~/orbitals
export ABACUS_PP_PATH="$HOME/pseudopotentials"
export ABACUS_ORBITAL_PATH="$HOME/orbitals"
```

> **教授提示**：设置这两个变量后，ASE 在转换结构时，如果检测到字典中指定的文件名存在于这些路径下，可以协助进行校验，虽然不是强制的，但在批量处理时能避免很多“找不到文件”的低级错误。

---

## 2.2 从零构建：以 FCC 铝为例

让我们从最简单的例子开始：使用 Python 脚本创建一个面心立方（FCC）的铝晶体，并将其导出为 ABACUS 标准的 STRU 文件。

### 2.2.1 编写构建脚本

创建一个名为 `build_fcc_al.py` 的文件：

```python
from ase.build import bulk
from ase.io import write
import os

# 1. 使用 ASE 内置函数构建 FCC 铝
# a=4.05 是晶格常数 (Angstrom)，cubic=True 表示生成常规晶胞而非原胞
al_fcc = bulk('Al', 'fcc', a=4.05, cubic=True)

# 2. 定义计算所需的物理参数映射
# 注意：文件名必须与你实际使用的文件一致
pp_map = {'Al': 'Al_ONCV_PBE-1.0.upf'}
basis_map = {'Al': 'Al_gga_7au_100Ry_4s4p1d.orb'}

# 3. 导出为 ABACUS STRU 格式
# format='abacus' 是调用接口的关键
write('STRU', al_fcc, format='abacus', pp=pp_map, basis=basis_map)

print("STRU file generated successfully.")
```

### 2.2.2 深入解析生成的 STRU 文件

运行上述脚本后，打开生成的 `STRU` 文件，你会看到如下关键部分。这里有几个**核心概念**必须理解：

```text
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb

LATTICE_CONSTANT
1.8897261258

LATTICE_VECTORS
4.050000 0.000000 0.000000
0.000000 4.050000 0.000000
0.000000 0.000000 4.050000

ATOMIC_POSITIONS
Direct
Al
0.0
4
0.000000 0.000000 0.000000 1 1 1 mag 0.0
0.000000 0.500000 0.500000 1 1 1 mag 0.0
...
```

#### 关键点解析：

1.  **LATTICE_CONSTANT 为什么是 1.8897...？**
    *   ABACUS 内核默认使用 **Bohr**（原子单位）作为长度单位。
    *   1 Angstrom $\approx$ 1.8897261258 Bohr。
    *   ASE 中的建模单位是 Angstrom。为了适配 ABACUS，接口自动将 `LATTICE_CONSTANT` 设为 1.8897...，这意味着 `LATTICE_VECTORS` 中的数值（如 4.05）会被乘以这个常数，从而在内部转换为 Bohr。
    *   **切记**：如果你手动修改 STRU，若 Vector 填埃，Constant 设 1.8897；若 Vector 填 Bohr，Constant 设 1.0。

2.  **LCAO vs PW 计算模式**
    *   上述例子生成了 `NUMERICAL_ORBITAL` 块，这是 **LCAO（线性组合原子轨道）** 基组计算所必需的。
    *   如果你做的是 **PW（平面波）** 计算，STRU 文件中**不应该**包含 `NUMERICAL_ORBITAL` 块。在 Python 脚本中，只需不传递 `basis` 参数即可。

3.  **坐标类型推荐**
    *   虽然 ABACUS 支持 `Cartesian_angstrom`（直接写埃单位坐标），但为了最大化兼容性（如兼容旧版本或可视化软件），推荐始终使用标准的 `LATTICE_CONSTANT` 缩放配合 `Direct`（分数坐标）或 `Cartesian`（以晶格常数为缩放单位的笛卡尔坐标）。

---

## 2.3 格式转换：CIF/POSCAR 转 STRU

这是实战中最常见的场景：你从 Materials Project 下载了 CIF 文件，或者从 VASP 迁移项目（POSCAR），需要转换为 STRU。

### 2.3.1 转换的核心难点
CIF 或 POSCAR 文件只包含几何结构信息（晶格和坐标），而 ABACUS 的 STRU 文件需要额外的物理信息：**赝势（PP）**、**轨道（Basis）** 以及 **磁矩（Mag）**。

因此，转换脚本的核心在于**“几何信息 + 物理定义”的拼接**。

### 2.3.2 通用转换脚本模板

以下是一个能够处理磁矩、赝势映射的完整转换脚本。建议保存为 `convert_to_stru.py` 复用。

```python
from ase.io import read, write
import numpy as np

# 1. 读取输入文件 (支持 .cif, POSCAR, .xsf 等)
# 假设我们有一个反铁磁的氧化镍结构
input_file = "NiO.cif"
atoms = read(input_file)

# 2. 设置磁矩 (ABACUS 的磁矩是在 STRU 中逐原子设置的，不同于 VASP 的 INCAR)
# 假设前两个原子是 Ni，后两个是 O，我们给 Ni 设置反铁磁序
# 注意：这里设置的是初始磁矩猜测值
mag_moments = [2.0, -2.0, 0.0, 0.0] 
atoms.set_initial_magnetic_moments(mag_moments)

# 3. 定义物理文件映射 (字典 Key 必须匹配元素符号)
# 平面波计算(PW)不需要 basis_dict
pp_dict = {
    'Ni': 'Ni_ONCV_PBE-1.0.upf',
    'O':  'O_ONCV_PBE-1.0.upf'
}

basis_dict = {
    'Ni': 'Ni_gga_8au_100Ry_4s3p1d.orb',
    'O':  'O_gga_7au_100Ry_2s2p1d.orb'
}

# 4. 写入 STRU
write(
    'STRU', 
    atoms, 
    format='abacus',
    pp=pp_dict,       # 必需：指定赝势
    basis=basis_dict, # LCAO必需，PW计算请去掉此参数
    units='Angstrom'  # 显式指定单位，触发 1.8897 缩放逻辑
)

print(f"Converted {input_file} to STRU with magnetic moments.")
```

### 2.3.3 磁矩设置的差异 (VASP vs ABACUS)
*   **VASP**: 在 `INCAR` 中使用 `MAGMOM = 2.0 -2.0 ...` 统一设置。
*   **ABACUS**: 磁矩属于结构属性，直接写在 `STRU` 文件的原子坐标行末尾（例如 `mag 2.0`）。
*   **ASE 接口的作用**: 当你调用 `atoms.set_initial_magnetic_moments` 后，ASE-ABACUS 接口会自动将这些值写入 STRU 文件对应的原子行中，无需手动编辑。

### 2.3.4 替代方案：ATOMKIT
如果你不想编写 Python 脚本，或者只需要进行快速的格式转换（不涉及复杂的磁矩或批量操作），可以使用 **ATOMKIT** 工具。

*   **特点**: 交互式命令行工具，由 VASPKIT 团队开发。
*   **适用性**: 适合“黑盒”转换，但对赝势和轨道的自动匹配灵活性不如 ASE 脚本高。
*   **命令示例**:
    ```bash
    atomkit
    # 选择 113 (Import CIF) -> 输入文件名
    # 选择 101 (Export ABACUS STRU)
    ```

---

**本章小结**：
掌握了 ASE 脚本化构建，你就拥有了批量处理成百上千个结构的能力。无论是从零构建晶体，还是将现有的 CIF/POSCAR 库转化为 ABACUS 格式，Python 字典映射（`pp` 和 `basis`）都是连接几何结构与 DFT 计算参数的桥梁。在下一章，我们将进入具体的计算参数文件 `INPUT` 的配置详解。

# 第三章：高效工具箱：ATOMKIT 与 abacustest

在前两章中，我们深入了解了 ABACUS 的输入参数和物理原理。然而，在实际科研工作中，手动编写 `STRU` 文件不仅效率低下，而且极易出错——尤其是在处理复杂超胞或需要批量测试数十种结构时。

本章将为你打开 ABACUS 的“高效工具箱”。对于不想编写 Python 脚本的初学者，我们有交互式的 **ATOMKIT**；对于需要高通量生成任务的进阶用户，我们有工程化的 **abacustest**。

但在介绍这些“黑盒”工具之前，必须先了解几个关于 `STRU` 文件的核心概念。无论你使用什么工具生成结构，理解这些底层逻辑都是排查报错（Debug）的关键。

---

## 3.0 教授的课前预警：STRU 文件的四大核心坑点

在开始使用自动化工具前，请务必理解以下四个手动修改或检查 `STRU` 时最容易混淆的概念。

### 1. ASE 安装避坑指南 (Critical)
很多用户习惯直接使用 `pip install ase` 安装 ASE，**这是错误的**。官方 ASE 尚未完全合并 ABACUS 的所有接口功能。
要使用 Python 脚本处理 ABACUS 结构，必须安装带有 ABACUS 接口的分支版本：

```bash
# ❌ 不要这样做
# pip install ase

# ✅ 正确的安装方式 (推荐在虚拟环境中进行)
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```

### 2. 神秘的 `1.8897...` 与 LATTICE_CONSTANT
打开一个生成的 `STRU` 文件，你常会看到这样的设置：
```text
LATTICE_CONSTANT
1.8897261258369282
```
**这是什么？**
ABACUS 内部计算使用原子单位制（Bohr），而晶体结构通常使用埃（Angstrom）。
*   `1 Bohr ≈ 0.529177 Angstrom`
*   `1 Angstrom ≈ 1.889726 Bohr`

**逻辑如下：**
*   **场景 A (推荐)**：你的 `LATTICE_VECTORS` 是以**埃 (Angstrom)** 为单位写的。此时 `LATTICE_CONSTANT` 必须设为 `1.8897...`，ABACUS 会自动将矢量乘以这个常数，转换为内部的 Bohr。
*   **场景 B**：你的 `LATTICE_VECTORS` 是以**Bohr** 为单位写的。此时 `LATTICE_CONSTANT` 设为 `1.0`。

> **风险提示**：虽然 `ATOMIC_POSITIONS` 支持 `Cartesian_angstrom` 选项（强制使用埃），但为了保证与后处理工具的最大兼容性，建议始终使用标准的 `LATTICE_CONSTANT` (1.8897...) 配合 `Direct`（分数坐标）或 `Cartesian`（以晶格常数为缩放因子的笛卡尔坐标）。

### 3. LCAO 与 PW 的区别
*   **平面波 (PW) 计算**：`STRU` 文件中**不需要** `NUMERICAL_ORBITAL` 这一栏。
*   **原子轨道 (LCAO) 计算**：`STRU` 文件中**必须包含** `NUMERICAL_ORBITAL`，且轨道文件的顺序必须与 `ATOMIC_SPECIES` 中的元素顺序严格对应。

### 4. 磁矩在哪里设置？
VASP 用户习惯在 `INCAR` 中通过 `MAGMOM` 设置磁矩。但在 ABACUS 中，**初始磁矩是结构的一部分**，必须在 `STRU` 文件中逐个原子设置。

```text
ATOMIC_POSITIONS
Direct
Fe
1.0  # 初始磁矩 (mag) - 低优先级，通常被忽略或作为默认值
2    # 原子数量
0.0 0.0 0.0 1 1 1 mag 2.5  # <--- 高优先级：针对该特定原子的磁矩
0.5 0.5 0.5 1 1 1 mag -2.5 # <--- 反铁磁设置示例
```

---

## 3.1 ATOMKIT：交互式快速转换

如果你手头有一个 `CIF` 文件，不想写 Python 代码，只想“一键”转成 `STRU`，**ATOMKIT** 是目前最便捷的交互式工具。它由 VASPKIT 团队开发，操作逻辑非常符合直觉。

### 3.1.1 安装与启动
下载并解压 ATOMKIT 后，直接运行：
```bash
~/software/atomkit/atomkit
```

### 3.1.2 实战演练：CIF 转 STRU
假设你当前目录下有一个 `Si.cif` 文件。启动 ATOMKIT 后，你会看到一个数字菜单。请按照以下“密码”顺序操作：

1.  **输入 `1`**：选择 `1) Export Structures`（导出结构）。
2.  **输入 `101`**：选择目标格式 `101) ABACUS (*.stru)`。
3.  **输入 `113`**：选择源文件格式 `113) CIF Format (*.cif)`。
    *   *注：你也可以直接输入 `113 Si.cif` 跳过文件名询问。*
4.  **输入文件名**：`Si.cif`。

**屏幕输出示例**：
```text
-->> (01) Reading Structure from Si.cif File...
...
-->> (02) Written Si.STRU File!
```

此时目录下会生成 `Si.STRU`。打开检查，你会发现 ATOMKIT 已经自动处理了：
*   晶格常数转换为 `1.889726`。
*   坐标转换为 `Direct` 分数坐标。
*   生成了 `ATOMIC_SPECIES` 和 `NUMERICAL_ORBITAL` 的模板。

> **注意**：ATOMKIT 无法预知你具体使用哪个版本的赝势和轨道文件。它生成的 `STRU` 文件中，赝势文件名通常是占位符（如 `Si.upf`），你需要手动修改为实际的文件名（如 `Si_ONCV_PBE-1.0.upf`）。

---

## 3.2 abacustest：工程化输入准备

当你需要对一个文件夹里的 50 个结构进行批量计算，或者希望自动匹配赝势和轨道文件时，交互式的 ATOMKIT 就显得太慢了。这时我们需要使用 **abacustest**。

它是 ABACUS 官方提供的测试与前处理框架，具备强大的 `model inputs` 功能。

### 3.2.1 核心命令：model inputs
该命令可以扫描指定目录，自动匹配结构、赝势和轨道，生成完整的输入文件集。

**基本语法**：
```bash
abacustest model inputs -f [结构文件列表] --pp [赝势目录] --orb [轨道目录]
```

### 3.2.2 实战演练：批量准备任务
假设你的工作目录结构如下：
```text
/work
├── structures/
│   ├── Si.cif
│   └── GaAs.cif
├── pp_dir/          # 存放 .upf 赝势文件
└── orb_dir/         # 存放 .orb 轨道文件
```

**步骤 1：准备环境**
确保 `pp_dir` 和 `orb_dir` 中的文件名包含元素符号（如 `Si_ONCV.upf`），或者在目录下提供 `element.json` 映射文件，以便程序能自动识别。

**步骤 2：执行命令**
```bash
abacustest model inputs \
    -f structures/*.cif \
    --pp ./pp_dir \
    --orb ./orb_dir \
    --d ./jobs  # 指定输出目录
```

**步骤 3：检查结果**
程序会自动在 `jobs` 目录下生成 `Si` 和 `GaAs` 两个子文件夹，每个文件夹内包含：
*   `STRU`：已自动填好对应的赝势和轨道文件名。
*   `INPUT`：基于默认模板生成的输入参数。
*   `KPT`：自动生成的 K 点文件。

### 3.2.3 进阶技巧：自定义参数模板
如果你不希望使用默认的 `INPUT` 参数，可以先准备一个 `INPUT_template` 文件，然后通过 `--input` 参数传递给 abacustest：

```bash
abacustest model inputs -f *.cif --pp ./pp --orb ./orb --input INPUT_template
```
这样，所有生成的任务都会继承你模板中的参数（如 `ecutwfc`, `mixing_type` 等），同时针对每个结构自动调整 `STRU`。

---

## 3.3 附录：ASE 高级脚本模板 (Power User)

对于习惯使用 Python 的用户，这里提供一份**能够同时处理赝势和轨道映射**的 ASE 转换脚本。这是实现“完美转换”的关键。

```python
from ase.io import read, write
import os

# 1. 读取结构 (支持 cif, poscar, xsd 等)
atoms = read("Si.cif")

# 2. 定义赝势映射 (必须与实际文件名一致)
# 键是元素符号，值是文件名
pp_map = {
    "Si": "Si_ONCV_PBE-1.0.upf",
    "O":  "O_ONCV_PBE-1.0.upf"
}

# 3. 定义轨道映射 (仅 LCAO 需要)
orb_map = {
    "Si": "Si_gga_7au_100Ry_4s4p1d.orb",
    "O":  "O_gga_6au_100Ry_2s2p1d.orb"
}

# 4. 设置磁矩 (可选)
# 注意：ABACUS 接口会读取 atoms 对象的初始磁矩并写入 STRU
# atoms.set_initial_magnetic_moments([1.0, -1.0]) 

# 5. 输出 STRU
# format='abacus' 是关键
write("STRU", atoms, format='abacus', pp=pp_map, basis=orb_map)

print("STRU file generated successfully with PP and ORB mappings.")
```

**总结**：
*   **新手/偶尔转换**：使用 **ATOMKIT** (113 -> 101)。
*   **批量/工程化**：使用 **abacustest** (`model inputs`)。
*   **定制化/脚本控**：使用 **ASE-ABACUS** (注意安装正确的分支)。

掌握了这些工具，你已经跨过了 ABACUS 使用门槛中最高的一级台阶。接下来，我们将正式运行计算，并解读输出文件。

# 第四章：进阶设置：磁性与动力学约束

上一章中，我们已经成功构建了最基础的晶体结构。但在实际的科研战场上，我们面对的往往不是完美的、静止的、非磁性的理想晶体。
本章我们将深入 `STRU` 文件的“深水区”。我们将学习如何处理**磁性体系**（这是 ABACUS 与 VASP 逻辑差异最大的地方之一），以及如何在结构优化中施加**动力学约束**（例如固定表面模型的底层原子）。


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

手动编辑几十上百个原子的 `0` 和 `1` 既痛苦又容易出错，推荐使用 Python 的 ASE 库来自动化这一过程。

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

## 附录：进阶学习指南

恭喜你完成了本书的学习！掌握 `STRU` 文件的构建与进阶设置只是第一步，计算材料学的海洋广阔无垠，以下是我们为你推荐的进阶探索路径：

**1. 延伸阅读与后续学习 (Extended Reading)**
- **基组优化进阶**：本书讨论了 LCAO 基组的使用，但基组本身的生成与优化（如针对特定高压或激发态体系）是提升精度的关键，建议查阅 ABACUS 关于数值原子轨道生成的专项教程。
- **高通量计算流 (High-throughput)**：结合 `abacustest` 与 `dpgen` 等工具，你可以尝试构建自动化的主动学习工作流，这在开发机器学习势函数（Deep Potential）时至关重要。
- **性能调优**：深入理解 `INPUT` 文件中关于并行效率（如 `KPAR`、`NPROC`）的设置，优化在大规模超算集群上的运行表现。

**2. 通用调试建议 (Troubleshooting)**
- **单位制对齐**：当计算结果异常时，首要检查 `LATTICE_CONSTANT` 的单位（Bohr vs Angstrom）。
- **赝势与基组匹配**：确保 `ATOMIC_SPECIES` 中引用的赝势与 `NUMERICAL_ORBITAL` 中的基组文件在元素种类和价电子设置上保持逻辑一致。
- **磁性收敛**：对于复杂的磁性体系，若 SCF 难以收敛，可尝试从简单的铁磁态开始，或逐步调整 `mixing_type` 参数。

**3. 官方资源与社区支持**
- **官方在线文档**：ABACUS 官方文档提供了最权威的参数手册，是解决疑难杂症的首选。
- **GitHub 社区**：积极参与 ABACUS 的 GitHub Issue 讨论，不仅能解决问题，还能与全球开发者共同推动软件的进步。

科研之路，贵在实践。愿此书能成为你探索微观世界的一把利刃！
