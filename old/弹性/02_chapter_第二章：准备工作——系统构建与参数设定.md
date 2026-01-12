# 第二章：准备工作——系统构建与参数设定

# 第二章：准备工作——系统构建与参数设定

在开始任何模拟之前，正确地准备初始结构并合理设置计算参数是非常重要的步骤。本章将详细介绍如何使用pymatgen库生成应变后的晶体结构、讨论超胞尺寸选择的原则、配置ABACUS输入文件以及处理磁性材料的特别方法。

## Section 2.1: 结构文件创建

### 使用pymatgen库生成应变后的晶体结构

在进行力学性质研究时，我们经常需要对晶体施加应变来观察其响应。这里以pymatgen库为例，展示如何生成应变后的晶体结构。

```python
from pymatgen.core import Structure, Lattice
import numpy as np

# 加载原始结构
structure = Structure.from_file("POSCAR")

# 定义应变张量
strain_tensor = np.array([[1.05, 0, 0],
                          [0, 1.05, 0],
                          [0, 0, 1.05]])

# 应用应变
strained_lattice = Lattice(np.dot(strain_tensor, structure.lattice.matrix))
strained_structure = Structure(strained_lattice, structure.species, structure.frac_coords)

# 保存应变后的结构
strained_structure.to(filename="STRU_strained")
```

### 超胞尺寸选择的原则

选择合适的超胞尺寸对于计算结果的准确性至关重要。通常情况下，超胞尺寸的选择应满足以下原则：

- **避免周期性边界条件的影响**：确保超胞尺寸足够大，使得相邻单元之间的相互作用可以忽略不计。
- **考虑计算资源**：超胞越大，计算成本越高。因此，需要在精度和计算效率之间找到平衡点。

### 关键参数

- `cal_stress`：是否计算应力张量。默认值为 `.FALSE.`，建议在力学性质研究中设置为 `.TRUE.`。
- `cal_force`：是否计算原子受力。默认值为 `.FALSE.`，建议在结构优化或分子动力学模拟中设置为 `.TRUE.`。

## Section 2.2: ABACUS输入文件配置

### 主要的INPUT文件选项解析

#### 应力和力的计算开关

- `cal_stress`：控制是否计算应力张量。设置为 `.TRUE.` 可以获得应力信息。
- `cal_force`：控制是否计算原子受力。设置为 `.TRUE.` 可以获得受力信息。

#### 展宽方法的选择

- `smearing_method`：展宽方法，常见的有 `gaussian` 和 `methfessel-paxton`。
- `smearing_sigma`：展宽参数，单位为 eV。推荐值为 `0.01` 到 `0.1` 之间。

#### 混合类型

- `mixing_type`：混合类型，常用的有 `pulay` 和 `anderson`。`pulay` 是默认值，适用于大多数情况。

### 示例 INPUT 文件

```plaintext
INPUT_PARAMETERS
suffix       example
ntype        1
ecutwfc      20
scf_thr      1e-7
basis_type   pw
calculation  scf
cal_stress   .TRUE.
cal_force    .TRUE.
smearing_method  gaussian
smearing_sigma  0.05
mixing_type  pulay
```

## Section 2.3: 磁性材料特别处理

### 初始化磁矩方向和大小

对于具有磁性的体系，初始化磁矩的方向和大小非常重要。这可以通过在 `STRU` 文件中添加 `magmom` 标签来实现。

#### 示例 STRU 文件

```plaintext
ATOMIC_SPECIES
Fe 55.845 Fe.upf
NUMERICAL_ORBITAL
Fe.orb
LATTICE_CONSTANT
1.0
LATTICE_VECTORS
  5.4309  0.0000  0.0000
  0.0000  5.4309  0.0000
 -2.7155 -2.7155  4.6274
ATOMIC_POSITIONS
Fe  0.0000  0.0000  0.0000  1.0  0.0  0.0
Fe  0.5000  0.5000  0.5000  1.0  0.0  0.0
```

在这个例子中，`magmom` 标签用于指定每个原子的磁矩方向和大小。例如，`1.0 0.0 0.0` 表示磁矩沿 x 方向，大小为 1.0 μB。

### 强调晶胞优化的重要性

确保初始结构处于势能面的局部极小值点，消除残余应力。这可以通过设置 `calculation` 参数为 `relax` 来实现。

### 详细说明应变-应力计算流程

1. **生成应变后的结构文件**：使用pymatgen库生成应变后的结构文件。
2. **运行ABACUS计算**：使用生成的结构文件进行ABACUS计算。
3. **处理结果**：使用pymatgen和dpdata处理计算结果，提取应力张量等信息。

### 强调K点采样和基组截断的敏感性

对于应力张量的计算，K点采样和基组截断对结果有较大影响。建议选择合适的K点网格和基组截断能量，以确保结果的准确性。

通过以上步骤，您可以正确地准备初始结构并合理设置计算参数，从而进行高质量的ABACUS模拟。