# 第三章：基于应力-应变法的弹性常数计算

一旦完成了所有前期准备工作，就可以进入到实际的弹性常数计算环节了。本章将重点介绍如何利用应力-应变法来获取材料的弹性常数，并通过具体的例子展示整个过程。

## Section 3.1: 应力-应变法简介

### 基本原理
应力-应变法是通过施加微小的应变到晶体结构上，然后计算相应的应力响应，从而拟合得到弹性张量的方法。这种方法的基本假设是在线弹性范围内，应力和应变之间存在线性关系，即胡克定律：
\[ \sigma_{ij} = C_{ijkl} \epsilon_{kl} \]
其中，\(\sigma_{ij}\) 是应力张量，\(\epsilon_{kl}\) 是应变张量，\(C_{ijkl}\) 是弹性刚度张量。

### 优势
1. **高信息密度**：单次 DFT 计算即可输出包含 6 个独立分量的应力张量，相比能量法仅输出标量能量，大大减少了所需的形变结构数量。
2. **数值稳定性**：基于广义胡克定律的线性拟合对 DFT 计算中的微小数值噪声具有更强的容忍度，拟合结果更加稳定可靠。

### 关键参数
- `[PARAMETER_MISSING]`：需要在 INPUT 文件中开启应力计算功能（请查阅官方文档确认具体参数名）。

## Section 3.2: 使用 pymatgen 生成变形后的结构

### 生成变形结构
为了进行应力-应变法计算，我们需要生成不同应变状态下的晶体结构文件。这可以通过 `pymatgen` 库轻松实现。

```python
from pymatgen import Structure
from pymatgen.core.tensors import symmetry_reduce, SquareTensor
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.abacus import AbacusStructure

# 读取原始结构
structure = Structure.from_file("POSCAR")

# 定义应变
strain_tensor = SquareTensor([[0.01, 0, 0],
                               [0, 0, 0],
                               [0, 0, 0]])

# 应用应变
deformed_structure = structure.apply_strain(strain_tensor)

# 保存变形后的结构
deformed_structure.to(fmt="abacus", filename="STRU_deformed")
```

## Section 3.3: 执行 ABACUS 计算

### 设置 ABACUS 输入文件
在生成变形结构后，我们需要设置 ABACUS 的输入文件并提交计算任务。以下是详细的步骤指南：

#### 1. 晶胞优化
在施加应变前，我们首先进行晶胞优化，确保初始结构处于势能面的局部极小值点。

```plaintext
INPUT_PARAMETERS
calculation cell-relax
basis_type lcao
pseudo_dir ../
orbital_dir ../
ecutwfc 50
kpoints_mp 8 8 8
gamma_only .true.
smearing_method gaussian
smearing_sigma 0.01
[PARAMETER_MISSING] .true.
[PARAMETER_MISSING] .true.
END
```

#### 2. 应力-应变计算
晶胞优化完成后，我们需要通过一个 JSON 配置文件来定义计算参数、指定输入结构以及设置形变幅度。

##### 配置文件准备 (`config.json`)
```json
{
  "base_dir": "/path/to/your/project",
  "structures": [
    {
      "name": "original",
      "file": "STRU"
    },
    {
      "name": "deformed_1",
      "file": "STRU_deformed_1"
    }
  ],
  "abacus_params": {
    "calculation": "relax",
    "basis_type": "lcao",
    "pseudo_dir": "../",
    "orbital_dir": "../",
    "ecutwfc": 50,
    "kpoints_mp": [8, 8, 8],
    "gamma_only": true,
    "smearing_method": "gaussian",
    "smearing_sigma": 0.01,
    "cal_stress": true,
    "cal_force": true
  }
}
```

#### 3. 提交计算任务
使用脚本或手动方式提交计算任务。

```bash
#!/bin/bash
for struct in $(cat config.json | jq -r '.structures[].name'); do
  cp ${struct}.STRU ./STRU
  abacus
done
```

### 监控进度
可以通过查看输出文件 `OUT.ABACUS/running_scf.log` 来监控计算进度。

## Section 3.4: 数据处理与结果分析

### 处理 ABACUS 输出数据
使用 `dpdata` 和 `numpy` 等工具来处理 ABACUS 输出的数据，并基于应力-应变曲线拟合出弹性常数。

```python
import dpdata
import numpy as np
from scipy.optimize import least_squares

# 读取数据
data = dpdata.LabeledSystem('OUT.ABACUS', fmt='abacus')

# 提取应力和应变
stress = data['virial'] / data['cells'].volume
strain = data['cells'] / data['orig_cells'] - 1

# 拟合弹性常数
def residual(C, strain, stress):
    return (np.dot(C, strain) - stress).flatten()

C_init = np.eye(6)
result = least_squares(residual, C_init, args=(strain, stress))
C = result.x.reshape(6, 6)

print("Elastic Constants (in GPa):")
print(C)
```

### 常见数据分析技巧
- **收敛检查**：确保每一步骤都充分收敛，特别是晶胞优化和应力张量计算。
- **K 点采样和基组截断**：选择合适的 K 点网格和基组截断对于提高计算精度至关重要。

## 附录：常见问题与进阶建议

### 内存消耗问题
由于涉及多次 DFT 计算，可能会导致较大的内存需求，建议使用高性能计算资源。

### 收敛问题
如果晶胞优化不充分或应力张量计算不准确，可能导致弹性常数拟合失败。务必仔细检查每一步骤。

### K 点采样和基组截断的选择
这对提高计算精度至关重要，建议根据实际情况选择合适参数。

### 后续学习方向
探索更复杂的材料体系、考虑温度效应等非零K条件下弹性常数的变化等。

通过以上步骤，您可以成功地使用 ABACUS 和应力-应变法计算材料的弹性常数。希望本章内容对您有所帮助！