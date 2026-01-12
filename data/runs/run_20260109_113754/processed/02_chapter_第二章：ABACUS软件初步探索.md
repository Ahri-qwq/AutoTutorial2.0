# 第二章：ABACUS软件初步探索

在掌握了基本理论后，接下来我们将重点放在如何实际使用ABACUS软件上。本章将详细介绍ABACUS的核心输入文件——INPUT文件中各个参数的意义及其对计算结果的影响，并教授如何根据目标材料制备正确的STRUCTURE文件，确保所有原子坐标、晶格矢量等信息准确无误。最后，我们将演示如何运行简单的SCF循环来获取材料的基态电子密度，并解释输出数据的意义。

## Section 2.1: INPUT文件详解

### 关键参数

#### `basis_type`
- **物理意义**: 指定计算中使用的基函数类型。
- **可选值**:
  - `pw`: 平面波基组
  - `lcao`: 原子轨道线性组合
- **默认值**: `pw`
- **推荐值**: 根据具体问题选择，平面波基组适用于大多数情况，而LCAO基组在某些情况下可以提供更高的效率。

#### `ecutwfc`
- **物理意义**: 平面波基组的能量截止值，用于控制平面波展开的精度。
- **单位**: Rydberg (Ry)
- **默认值**: 未指定（必须手动设置）
- **推荐值**: 通常为40-80 Ry，具体取决于材料和所需的精度。

#### `scf_nmax`
- **物理意义**: 自洽场迭代的最大步数。
- **默认值**: 50
- **推荐值**: 根据具体问题调整，一般建议从默认值开始，必要时增加以确保收敛。

#### `relax_nmax`
- **物理意义**: 结构优化过程中的最大步数。
- **默认值**: 50
- **推荐值**: 根据具体问题调整，一般建议从默认值开始，必要时增加以确保结构优化收敛。

### 完整的INPUT文件示例

```plaintext
INPUT_PARAMETERS
suffix       Mn           # 输出后缀
ntype        1            # 元素种类
ecutwfc      60           # 展开截止能量
scf_thr      1e-7         # 电荷密度收敛阈值
basis_type   pw           # 基函数类型
calculation  scf          # 计算类型
scf_nmax     100          # SCF迭代最大步数
relax_nmax   50           # 结构优化最大步数
```

## Section 2.2: 结构文件（STRU）准备

### 关键参数

- **`ATOMIC_SPECIES`**: 定义元素名称、相对原子质量以及赝势文件名。
- **`LATTICE_CONSTANT`**: 晶格常数，默认单位为Bohr。
- **`LATTICE_VECTORS`**: 晶格向量。
- **`ATOMIC_POSITIONS`**: 原子位置。

### 完整的STRU文件示例

```plaintext
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf

LATTICE_CONSTANT
1.8897259886  # 1.8897259886 Bohr = 1 Angstrom

LATTICE_VECTORS
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0

ATOMIC_POSITIONS
U 0.0 0.0 0.0
```

## Section 2.3: 执行基本的第一性原理计算

### 运行简单的SCF循环

#### 命令
```bash
mpirun -n 4 abacus
```

#### 解释输出数据

- **`TOTAL Time`**: 总计算时间。
- **`SEE INFORMATION IN`**: 输出文件路径。
- **`FINISH Time`**: 计算结束时间。

### 将ABACUS的结果传递给Wannier90

为了将ABACUS的结果传递给Wannier90，需要进行以下步骤：

1. **生成Wannier90所需的输入文件**:
   - 使用ABACUS提供的工具或脚本生成Wannier90所需的输入文件（如`wannier90.win`）。
   - 请参考最新的ABACUS和Wannier90官方文档获取最准确的信息。

2. **配置Wannier90输入文件**:
   - 在`wannier90.win`文件中，配置必要的参数，如`num_wann`、`dis_win_min`、`dis_win_max`等。
   - 请参考Wannier90官方文档获取详细的参数说明。

3. **运行Wannier90**:
   - 使用Wannier90执行计算，生成能带结构和其他相关结果。

### 风险提示

关于与Wannier90接口的具体参数配置资料较为缺乏，请读者务必参考最新的ABACUS和Wannier90官方文档获取最准确的信息。

通过以上步骤，您可以成功地使用ABACUS进行第一性原理计算，并将结果传递给Wannier90进行进一步分析。希望本章的内容能够帮助您更好地理解和使用ABACUS软件。