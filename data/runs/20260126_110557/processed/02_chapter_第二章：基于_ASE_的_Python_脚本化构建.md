# 第二章：基于 ASE 的 Python 脚本化构建

在计算材料学的实际工作中，手动编写结构文件不仅效率低下，而且极易出错。作为 ABACUS 的开发者，我强烈推荐使用 **ASE (Atomic Simulation Environment)** 作为构建和操作原子结构的首选工具。

ASE 是 Python 生态中最强大的原子模拟接口库。本章将带你掌握如何利用 ASE 的灵活性，完成从零建模、格式转换（如 POSCAR 转 STRU）以及批量生成输入文件的全流程。

---

## 2.1 环境配置与 ASE-ABACUS 接口安装

### 2.1.1 关键警告：安装源的选择
这是新手最容易踩的第一个坑。虽然 PyPI 官方源提供了 `ase` 包（即 `pip install ase`），但官方版本对 ABACUS 的支持非常有限且更新滞后。

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