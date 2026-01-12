# 第一章：弹性常数的基础理论

在开始具体计算之前，首先需要建立对弹性常数物理意义的理解，这是后续所有计算工作的基础。本章将介绍弹性常数的概念、胡克定律与Voigt标记法、应力张量与应变张量等内容。

## 1.1 弹性常数的概念与重要性

### 1.1.1 定义与物理意义

弹性常数是表征材料在弹性极限内抵抗外力导致的可逆形变能力的物理量。它描述了晶体微观化学键的强度及各向异性特征，宏观上决定了材料的刚度（Stiffness）、硬度以及机械稳定性。通过弹性常数，我们可以评估材料在不同方向上的力学性能，从而为材料设计和应用提供重要的参考依据。

### 1.1.2 在材料力学性能评估中的作用

- **刚度**：弹性常数直接反映了材料的刚度，即材料在外力作用下抵抗变形的能力。
- **硬度**：虽然硬度通常由其他实验方法测量，但弹性常数可以间接反映材料的硬度特性。
- **机械稳定性**：弹性常数可以帮助我们判断材料在受力时是否会发生塑性变形或断裂，从而评估其机械稳定性。

## 1.2 胡克定律与Voigt标记法

### 1.2.1 广义胡克定律

在连续介质力学的线弹性近似下，晶体内部的应力张量 \(\sigma\) 与应变张量 \(\epsilon\) 呈线性关系。对于微小的形变，其本构关系由广义胡克定律描述：

\[
\sigma_{ij} = C_{ijkl} \epsilon_{kl}
\]

其中：
- \(\sigma_{ij}\) 为二阶应力张量（Stress Tensor）。
- \(\epsilon_{kl}\) 为二阶应变张量（Strain Tensor）。
- \(C_{ijkl}\) 为四阶弹性刚度张量（Elastic Stiffness Tensor），包含 81 个分量。

由于 \(\sigma_{ij}\) 和 \(\epsilon_{kl}\) 均为对称张量（即 \(\sigma_{ij} = \sigma_{ji}\) 和 \(\epsilon_{kl} = \epsilon_{lk}\)），利用张量对称性，可引入 Voigt 标记法将四阶张量 \(C_{ijkl}\) 降维映射为 6×6 的对称矩阵 \(C_{\alpha\beta}\)。

### 1.2.2 Voigt 标记法

Voigt 标记法将应力张量和应变张量的笛卡尔坐标分量映射到一个简化的 6 维向量中，具体规则如下：

\[
\begin{aligned}
\sigma_1 &= \sigma_{xx}, & \sigma_2 &= \sigma_{yy}, & \sigma_3 &= \sigma_{zz}, \\
\sigma_4 &= \sigma_{yz}, & \sigma_5 &= \sigma_{xz}, & \sigma_6 &= \sigma_{xy},
\end{aligned}
\]

\[
\begin{aligned}
\epsilon_1 &= \epsilon_{xx}, & \epsilon_2 &= \epsilon_{yy}, & \epsilon_3 &= \epsilon_{zz}, \\
\epsilon_4 &= \epsilon_{yz}, & \epsilon_5 &= \epsilon_{xz}, & \epsilon_6 &= \epsilon_{xy}.
\end{aligned}
\]

此时，胡克定律可简化为矩阵形式：

\[
\sigma_\alpha = C_{\alpha\beta} \epsilon_\beta
\]

### 1.2.3 独立弹性常数的数量

虽然 \(C_{\alpha\beta}\) 的 Voigt 矩阵包含 36 个分量，但因其是对称矩阵，所以对于最一般的晶体，其独立弹性常数并非 36 个，而是缩减为 21 个（6 个对角项 + 15 个非对角项）。此外，晶体点群对称性会进一步施加几何约束，使得独立分量数量继续减少：

- 三斜晶系 (Triclinic)：无额外对称性，保持 21 个独立分量。
- 单斜晶系 (Monoclinic)：13 个独立分量。
- 正交晶系 (Orthorhombic)：9 个独立分量。

## 1.3 应力张量与应变张量

### 1.3.1 定义与关系

- **应力张量 \(\sigma_{ij}\)**：描述了材料内部单位面积上的力，表示为：

\[
\sigma_{ij} = \frac{F_i}{A_j}
\]

- **应变张量 \(\epsilon_{kl}\)**：描述了材料在外部力作用下的相对位移，表示为：

\[
\epsilon_{kl} = \frac{1}{2} \left( \frac{\partial u_k}{\partial x_l} + \frac{\partial u_l}{\partial x_k} \right)
\]

其中，\(u_k\) 表示位移矢量。

### 1.3.2 重要性

应力张量和应变张量是描述材料内部状态的重要工具。它们之间的关系不仅揭示了材料在外力作用下的响应行为，还为我们提供了计算弹性常数的基础。通过这些张量，我们可以更深入地理解材料的力学性能，并为实际应用提供指导。

## 1.4 晶胞优化的重要性

在进行弹性常数计算之前，必须进行晶胞优化，消除残余应力，确保初始结构处于势能面的局部极小值点。这一步骤至关重要，因为未优化的晶胞可能会导致计算结果的偏差。在 ABACUS 中，可以通过设置 `relax` 参数来实现晶胞优化。

### 1.4.1 晶胞优化的步骤

1. **准备输入文件**：创建 `INPUT` 文件并设置必要的参数。
2. **运行优化**：使用 ABACUS 运行优化任务。
3. **检查结果**：确保优化后的晶胞没有残余应力。

### 1.4.2 示例代码

```plaintext
&GLOBAL
  run_type = 'relax'
/
&CELL_RELAX
  cell_relax_method = 'bfgs'
/
&ELECTRON
  pseudopot_dir = './pp_orb'
  pseudopot_file = 'C_ONCV_PBE-1.0.upf'
  orbital_file = 'C_gga_7au_100Ry_2s2p1d.orb'
  smearing_method = 'gaussian'
  mixing_type = 'pulay'
  mixing_beta = 0.7
  scf_thr = 1.0e-6
/
```

## 1.5 使用 `pymatgen` 进行弹性常数计算

### 1.5.1 外部工具的作用

- **`pymatgen`**：用于计算弹性常数的主要工具包。
- **`dpdata`**：用于处理 ABACUS 输出的数据。
- **`monty`**：用于处理脚本执行过程中的临时文件和日志。

### 1.5.2 安装依赖库

```bash
pip install monty numpy dpdata pymatgen
```

### 1.5.3 计算流程

1. **准备输入文件**：创建 `INPUT` 文件并设置必要的参数。
2. **运行 ABACUS 计算**：使用 `run_task.sh` 和 `sub.sh` 批量运行 ABACUS 计算。
3. **处理数据**：使用 `dpdata` 处理 ABACUS 输出的数据。
4. **计算弹性常数**：使用 `pymatgen` 计算弹性常数。

### 1.5.4 示例代码

#### 1.5.4.1 `gene_dfm.py` 脚本

```python
from pymatgen import Structure
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor
import dpdata

# 读取结构
structure = Structure.from_file("POSCAR")

# 生成应变
deformations = [Deformation([[1.01, 0, 0], [0, 1, 0], [0, 0, 1]]),
                Deformation([[1, 0, 0], [0, 1.01, 0], [0, 0, 1]]),
                Deformation([[1, 0, 0], [0, 1, 0], [0, 0, 1.01]])]

# 保存应变后的结构
for i, deformation in enumerate(deformations):
    deformed_structure = deformation.apply_to_structure(structure)
    deformed_structure.to(filename=f"POSCAR_{i+1}")
```

#### 1.5.4.2 `compute_dfm.py` 脚本

```python
from pymatgen import Structure
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor
import dpdata

# 读取应力数据
stress_data = dpdata.LabeledSystem("path/to/stress_data").get_stresses()

# 读取应变数据
strains = [Strain.from_deformation(Deformation.from_index(i)) for i in range(len(stress_data))]

# 计算弹性常数
elastic_tensor = ElasticTensor.from_stress_strain(stress_data, strains)
print(elastic_tensor)
```

### 1.5.5 参数设置的准确性

在进行弹性常数计算时，参数设置的准确性至关重要。特别是 `cal_stress` 和 `cal_force` 参数的设置，这两个参数分别控制应力和力的计算。请查阅 ABACUS 官方文档以确认具体的参数名和推荐值。

```plaintext
&FORCE_EVAL
  cal_stress = .TRUE.
  cal_force = .TRUE.
/
```

请注意，上述参数名可能需要根据 ABACUS 官方文档进行调整。如果不确定具体参数名，请查阅官方文档。

## 1.6 风险提示

- **关于应力计算参数的资料缺失**：撰写时请提醒读者查阅 ABACUS 官方文档，确认具体的参数名和推荐值。

通过以上内容，我们介绍了弹性常数的基础理论、胡克定律与Voigt标记法、应力张量与应变张量、晶胞优化的重要性以及如何使用 `pymatgen` 进行弹性常数计算。希望这些知识能够帮助读者更好地理解和应用弹性常数计算方法。