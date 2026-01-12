# 第三章：准备输入文件

在开始使用 ABACUS 进行计算之前，必须准备好所有必要的输入文件。这些文件包括 `INPUT`、`STRU` 和 `KPT` 文件，它们直接影响到最终结果的准确性。本章将详细介绍如何生成或获取合适的结构文件、设置 `INPUT` 文件中的关键参数以及合理设置 K 点采样和选择赝势文件。

## Section 3.1: 结构文件准备

### 晶胞优化的重要性
在施加应变前，必须进行晶胞优化，消除残余应力，确保初始结构处于势能面的局部极小值点。这一步骤对于获得准确的弹性常数至关重要。

### 生成或获取结构文件
结构文件通常包含原子种类、原子位置、晶格常数以及晶格向量等信息。可以通过以下几种方式生成或获取结构文件：

1. **从实验数据或数据库中获取**：例如，可以从Materials Project或Crystallography Open Database (COD)下载CIF文件，并将其转换为ABACUS所需的格式。
2. **使用材料建模软件**：如VASP、Quantum ESPRESSO等，可以生成POSCAR或QE格式的结构文件，然后使用外部工具（如`dpdata`）将其转换为ABACUS所需的`STRU`格式。

### STRU文件格式
`STRU`文件的格式如下：

```plaintext
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf
# 元素名称、相对原子质量、原子赝势的文件名

LATTICE_CONSTANT
1.8897259886  # 晶格常数的单位，默认单位为Bohr

LATTICE_VECTORS
a1 a2 a3
b1 b2 b3
c1 c2 c3
# 晶格三个边对应的向量

ATOMIC_POSITIONS
Direct
x1 y1 z1
x2 y2 z2
...
# 原子的位置坐标
```

### 使用 `dpdata` 转换结构文件
`dpdata` 是一个非常有用的库，可以方便地将不同格式的结构文件转换为ABACUS所需的`STRU`格式。例如，将POSCAR文件转换为`STRU`文件：

```python
from dpdata import LabeledSystem

# 读取POSCAR文件
system = LabeledSystem('POSCAR', fmt='vasp/poscar')
# 转换为STRU文件
system.to('abacus/STRU', 'STRU')
```

## Section 3.2: 设置 INPUT 文件

### 关键参数设置
`INPUT`文件是ABACUS的主要配置文件，其中包含了计算所需的各种参数。以下是几个关键参数及其设置方法：

#### cal_stress
需要在 `INPUT` 文件中开启应力计算功能（请查阅官方文档确认具体参数名）。该参数用于计算系统的应力张量，对于弹性常数的计算非常重要。

#### cal_force
同样，需要在 `INPUT` 文件中开启力的计算功能（请查阅官方文档确认具体参数名）。该参数用于计算每个原子上的受力，这对于结构优化和分子动力学模拟是必需的。

#### gamma_only
`gamma_only` 参数决定是否只在Gamma点进行计算。对于某些对称性较高的系统，可以设置为 `1` 以提高计算效率。默认值为 `0`。

#### smearing_method
`smearing_method` 参数指定电子态密度展宽的方法。常用的选项有 `gaussian`、`fermi_dirac` 等。推荐使用 `gaussian` 方法，因为它在大多数情况下表现良好。

#### smearing_sigma
`smearing_sigma` 参数指定展宽的宽度。典型的值为 `0.002` 到 `0.01`，单位为 Hartree。较小的值适用于金属体系，较大的值适用于绝缘体和半导体。

#### mixing_type
`mixing_type` 参数指定电荷密度混合的方法。常见的选项有 `pulay`、`anderson` 等。推荐使用 `pulay` 方法，它在大多数情况下表现稳定且收敛速度快。

### 示例 INPUT 文件
以下是一个示例 `INPUT` 文件，展示了上述参数的设置：

```plaintext
INPUT_PARAMETERS
suffix       Mn           # 输出后缀
ntype        1            # 元素种类
ecutwfc      20           # 展开截止能量
scf_thr      1e-7         # 电荷密度收敛阈值
basis_type   pw           # 基函数类型
calculation  scf          # 计算类型
cal_stress   [PARAMETER_MISSING]  # 应力计算
cal_force    [PARAMETER_MISSING]  # 力的计算
gamma_only   0            # 是否只在Gamma点计算
smearing_method  gaussian  # 展宽方法
smearing_sigma  0.002     # 展宽宽度
mixing_type  pulay        # 电荷密度混合方法
```

## Section 3.3: K点采样与赝势选择

### K点采样
K点采样是DFT计算中的一个重要步骤，它决定了布里渊区的采样密度。合理的K点采样可以提高计算精度，但也会增加计算成本。K点采样的设置取决于系统的对称性和计算需求。

### KPT文件格式
`KPT`文件的格式如下：

```plaintext
K_POINTS
MP
0 0 0 1
# K点的坐标和权重
```

### 赝势选择
赝势文件的选择对于计算结果的准确性至关重要。赝势文件天然会带着交换关联泛函，目前最常用的是PBE交换关联泛函。赝势文件可以从ABACUS官方网站下载：

- [ABACUS赝势下载](https://abacus.ustc.edu.cn/pseudo/list.htm)

### 示例 KPT 文件
以下是一个示例 `KPT` 文件，展示了K点的设置：

```plaintext
K_POINTS
MP
0 0 0 1
```

### 外部工具的使用
- **pymatgen**: 用于弹性常数计算。可以通过 `pymatgen` 库进行弹性常数的分析和可视化。
- **dpdata**: 用于结构文件的转换和处理。
- **monty**: 用于自动化脚本编写和任务管理。

### 风险提示
- 在设置应力和力的计算参数时，请务必查阅ABACUS官方文档，确认具体的参数名和推荐值。

通过以上步骤，您可以准备好所有必要的输入文件，从而顺利进行ABACUS计算。