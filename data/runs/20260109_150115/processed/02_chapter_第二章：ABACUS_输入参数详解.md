# 第二章：ABACUS 输入参数详解

在掌握了基本理论之后，下一步是了解如何正确配置 ABACUS 的输入文件以进行弹性常数计算。本章将逐一解析主要的输入参数，包括它们的物理意义以及推荐值。

## Section 2.1: 应力与力的计算相关参数

### 内容
在材料科学中，应力和力的计算对于理解材料的力学性质至关重要。ABACUS 提供了专门的参数来控制这些计算。本节将详细讲解 `cal_stress` 和 `cal_force` 参数的作用，强调这两个选项对于获取准确应力张量数据的重要性。

### 关键参数
- **[PARAMETER_MISSING]**: 开启应力计算（请查阅官方文档确认具体参数名）。推荐值为 1。
- **[PARAMETER_MISSING]**: 开启力的计算（请查阅官方文档确认具体参数名）。推荐值为 1。

#### 物理意义
- **应力计算**：应力张量描述了材料内部各点上的内力分布情况。通过开启应力计算，可以得到每个原子上的应力分量，这对于研究材料的力学性能非常重要。
- **力的计算**：力的计算提供了每个原子上受到的外力信息，这对于结构优化和分子动力学模拟非常关键。

#### 推荐值
- **[PARAMETER_MISSING]**: 1（开启应力计算）
- **[PARAMETER_MISSING]**: 1（开启力的计算）

## Section 2.2: 电子结构计算相关参数

### 内容
电子结构计算是 DFT 计算的核心部分，直接影响到计算结果的精度和效率。本节将解释 `gamma_only`, `smearing_method`, `smearing_sigma` 和 `mixing_type` 等参数如何影响电子结构计算的精度和效率。

### 关键参数
- **`gamma_only`**: 是否只使用 \(\Gamma\) 点进行计算（推荐值: 0）。
- **`smearing_method`**: 能量展宽方法（推荐值: gaussian）。
- **`smearing_sigma`**: 能量展宽参数（推荐值: 0.002）。
- **`mixing_type`**: 电子密度混合类型（需查阅文档确认推荐设置）。

#### 物理意义
- **`gamma_only`**: 该参数决定是否只使用倒空间中的 \(\Gamma\) 点进行计算。通常情况下，设置为 0 表示使用更多的 k 点进行采样，从而提高计算精度。
- **`smearing_method`**: 为了处理费米能级附近的不连续性，DFT 计算中通常引入能量展宽方法。常用的展宽方法有高斯展宽（gaussian）、甲基展宽（methfessel-paxton）等。
- **`smearing_sigma`**: 能量展宽参数，用于控制展宽的宽度。较小的值会提高计算精度，但可能增加收敛难度。
- **`mixing_type`**: 电子密度混合类型，不同的混合方法会影响电子密度的收敛速度和稳定性。常见的混合方法有 Pulay 混合、Broyden 混合等。

#### 推荐值
- **`gamma_only`**: 0（使用多个 k 点进行采样）
- **`smearing_method`**: gaussian
- **`smearing_sigma`**: 0.002
- **`mixing_type`**: 请查阅 ABACUS 官方文档确认推荐设置

## 风险提示
- **`mixing_type` 参数的推荐设置资料不足**，撰写时请提醒读者查阅 ABACUS 官方文档。
- **关于开启应力计算的具体参数名资料不足**，撰写时请提醒读者查阅 ABACUS 官方文档。

## 文件准备

使用 ABACUS 进行计算最重要的是需要准备 `INPUT`, `KPT`, `STRU` 三个文件。以下是一个简单的示例：

### `INPUT` 文件
```plaintext
INPUT_PARAMETERS
suffix       Mn           # 输出后缀
ntype        1            # 元素种类
ecutwfc      20           # 展开截止能量 (Ry)
scf_thr      1e-7         # 电荷密度收敛阈值
basis_type   pw           # 基函数类型 (pw 或 lcao)
calculation  scf          # 计算类型 (自洽场计算)
[PARAMETER_MISSING] 1     # 开启应力计算
[PARAMETER_MISSING] 1     # 开启力的计算
gamma_only   0            # 是否只使用 \(\Gamma\) 点进行计算
smearing_method gaussian  # 能量展宽方法
smearing_sigma 0.002      # 能量展宽参数
mixing_type  [MIXING_TYPE]  # 电子密度混合类型 (请查阅官方文档确认推荐设置)
```

### `STRU` 文件
```plaintext
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf   # 元素名称、相对原子质量、原子赝势的文件名

LATTICE_CONSTANT    # 晶格常数的单位，默认单位为Bohr
1.8897259886  # 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS   # 晶格三个边对应
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0

ATOMIC_POSITIONS
Direct
0.0 0.0 0.0
```

### `KPT` 文件
```plaintext
K_POINTS
0
Automatic
4 4 4 0 0 0  # K 点网格及偏移
```

### 小技巧
- 在 ABACUS 的输入文件 `INPUT` 中可以设置变量 `stru_file`，该变量对应的原子构型文件为 `STRU-001`，则 ABACUS 可以直接读取该结构文件。

## 结论
通过正确配置 ABACUS 的输入参数，可以有效地进行弹性常数计算。特别需要注意的是，晶胞优化的重要性，确保初始结构无残余应力。此外，准备 JSON 配置文件时，需要定义计算参数、指定输入结构以及设置形变幅度。最后，提醒读者注意应力-应变法的数值稳定性，特别是在 K 点采样和基组截断方面。