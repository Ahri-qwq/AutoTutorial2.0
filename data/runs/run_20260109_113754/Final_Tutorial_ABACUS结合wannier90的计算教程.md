# ABACUS 实战教程：ABACUS结合wannier90的计算教程

# 第一章：第一性原理计算与Wannier函数概览

## 本章逻辑
作为入门章节，首先介绍必要的背景知识，包括第一性原理计算的基本思想以及Wannier函数的概念，为后续具体操作打下坚实的理论基础。

### Section 1.1: 第一性原理计算简介

#### 内容
第一性原理计算是基于量子力学基本原理的一种计算方法，它通过求解薛定谔方程来预测材料的性质。这种方法不需要任何实验参数，只需要知道原子核和电子的基本物理常数即可。第一性原理计算可以用于研究材料的各种性质，如电子结构、光学性质、磁学性质等。

#### 关键参数
- `basis_type`：基组类型，例如平面波（`pw`）或数值原子轨道线性组合（`lcao`）。
- `ecutwfc`：波函数截断能，单位为Rydberg (Ry) 或电子伏特 (eV)，用于确定平面波基组的大小。
- `scf_thr`：自洽场（SCF）计算的能量收敛阈值，单位为Rydberg (Ry) 或电子伏特 (eV)。

### Section 1.2: Wannier函数的作用与优势

#### 内容
Wannier函数是一种局域化的电子态，它们在描述复杂系统的电子结构时具有显著的优势。Wannier函数可以通过将扩展的Bloch态投影到一组局域化的基函数上来构造。这种局域化使得Wannier函数非常适合于分析电子间的相互作用、化学键的形成以及电子传输特性。

#### 关键参数
- `num_wann`：Wannier函数的数量。
- `dis_win_min` 和 `dis_win_max`：定义了能量窗口，用于选择参与Wannier函数构造的能带。
- `dis_froz_min` 和 `dis_froz_max`：定义了冻结窗口，用于指定哪些能带被冻结以保持其原始形状。

### Section 1.3: 将ABACUS的结果传递给wannier90

#### 内容
为了将ABACUS的第一性原理计算结果传递给wannier90进行进一步分析，需要执行以下步骤：

1. **准备ABACUS输入文件**：
   - `INPUT` 文件：包含基本的计算参数，如 `basis_type`, `ecutwfc`, `scf_thr` 等。
   - `STRU` 文件：包含晶体结构信息，如原子种类、位置、晶格常数等。
   - `KPT` 文件：包含k点路径信息。

2. **运行ABACUS计算**：
   - 使用ABACUS进行自洽迭代计算，生成波函数和能带结构。

3. **后处理步骤**：
   - 从ABACUS输出中提取必要的数据，如波函数和能带结构。
   - 使用ABACUS提供的工具（如 `abacus_to_wannier90.py`）将这些数据转换为wannier90所需的格式。

4. **构建wannier90输入文件**：
   - 创建 `wannier90.win` 文件，包含Wannier函数的相关参数，如 `num_wann`, `dis_win_min`, `dis_win_max`, `dis_froz_min`, `dis_froz_max` 等。
   - 运行wannier90进行Wannier函数的构造。

#### 示例代码

**ABACUS `INPUT` 文件示例**
```plaintext
INPUT_PARAMETERS
ntype 1
pseudo_dir ../../PP_ORB
ecutwfc 100
scf_thr 1e-6
basis_type pw
calculation scf
```

**ABACUS `STRU` 文件示例**
```plaintext
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf
LATTICE_CONSTANT
1.8897259886
LATTICE_VECTORS
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
ATOMIC_POSITIONS
U 0.0 0.0 0.0
```

**ABACUS `KPT` 文件示例**
```plaintext
K_POINTS
Automatic
0 0 0 0 0 0
```

**wannier90 `wannier90.win` 文件示例**
```plaintext
num_bands = [NUM_BANDS]
num_wann = [NUM_WANN]
dis_win_min = [DIS_WIN_MIN]
dis_win_max = [DIS_WIN_MAX]
dis_froz_min = [DIS_FROZ_MIN]
dis_froz_max = [DIS_FROZ_MAX]

begin unit_cell_cart
[UNIT_CELL_CART]
end unit_cell_cart

begin atoms_cart
[ATOMS_CART]
end atoms_cart

begin kpoints
[KPOINTS]
end kpoints

mp_grid : [MP_GRID]
```

#### 风险提示
关于与wannier90接口的具体参数配置资料较为缺乏，请读者参考最新的ABACUS和wannier90官方文档获取最准确的信息。

通过以上步骤，您可以将ABACUS的第一性原理计算结果成功传递给wannier90，并进行进一步的Wannier函数分析。希望本章的内容能够帮助您更好地理解和应用这些强大的计算工具。

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

# 第三章：连接ABACUS与wannier90

## 第三章：连接ABACUS与wannier90

在本章中，我们将探讨如何有效地将ABACUS与wannier90结合起来工作，以生成并分析Wannier函数。通过这一过程，我们可以更深入地理解材料的电子结构和拓扑性质。

### Section 3.1: wannier90简介

**内容**: 简要回顾wannier90的功能特点及其在固体物理学研究中的重要地位。

wannier90是一款开源软件包，用于计算和分析Wannier函数。Wannier函数是一种局域化的单电子波函数，可以用来描述固体中的电子态。它们在凝聚态物理中有着广泛的应用，特别是在能带结构、电荷密度、磁性以及拓扑性质的研究中。wannier90能够从第一性原理计算（如DFT）得到的波函数出发，构造出最大局域化的Wannier函数，并进一步用于计算各种物理量。

**关键参数**:
- `num_wann`: Wannier函数的数量。
- `dis_win_min` 和 `dis_win_max`: 能带窗口的最小值和最大值。
- `dis_froz_min` 和 `dis_froz_max`: 内部冻结窗口的最小值和最大值。
- `write_hr`: 是否写入哈密顿量矩阵到文件中。
- `write_tb`: 是否写入紧束缚模型到文件中。

### Section 3.2: 设置wannier90环境

**内容**: 提供详细的步骤指导用户如何安装并配置好wannier90软件，确保能够顺利与ABACUS对接。

#### 安装wannier90

1. **下载源代码**:
   ```bash
   git clone https://github.com/wannier-developers/wannier90.git
   cd wannier90
   ```

2. **编译源代码**:
   ```bash
   make config
   make
   sudo make install
   ```

3. **设置环境变量**:
   将wannier90的可执行文件路径添加到系统PATH中：
   ```bash
   export PATH=$PATH:/path/to/wannier90/bin
   ```

#### 配置wannier90

- **检查依赖库**: 确保所有依赖库已正确安装。
- **测试安装**: 运行一个简单的测试案例来验证安装是否成功。

### Section 3.3: 从ABACUS到wannier90的数据传递

**内容**: 详细说明如何正确地将ABACUS产生的波函数等信息转换成适合wannier90处理的形式，并构造相应的输入文件。

#### ABACUS输出文件准备

1. **运行SCF计算**:
   使用ABACUS进行自洽场（SCF）计算，生成波函数文件。例如：
   ```bash
   abacus INPUT
   ```

2. **生成wannier90输入文件**:
   ABACUS提供了工具`abacus2wannier90`，可以将ABACUS的输出文件转换为wannier90所需的输入文件。具体步骤如下：
   ```bash
   abacus2wannier90 -i INPUT -o wannier90_input
   ```

#### 构造wannier90输入文件

1. **创建`wannier90.win`文件**:
   根据ABACUS的输出结果，手动编辑或使用脚本生成`wannier90.win`文件。以下是一个示例：
   ```plaintext
   num_bands = [PARAMETER_MISSING]  # 总带数
   num_wann = [PARAMETER_MISSING]  # Wannier函数数量
   dis_win_min = [PARAMETER_MISSING]  # 能带窗口最小值
   dis_win_max = [PARAMETER_MISSING]  # 能带窗口最大值
   dis_froz_min = [PARAMETER_MISSING]  # 内部冻结窗口最小值
   dis_froz_max = [PARAMETER_MISSING]  # 内部冻结窗口最大值
   write_hr = .true.
   write_tb = .true.
   ```

2. **运行wannier90**:
   在命令行中运行wannier90：
   ```bash
   wannier90.x -pp wannier90_input
   wannier90.x wannier90_input
   ```

### Section 3.4: 分析Wannier函数结果

**内容**: 介绍如何解读由wannier90生成的Wannier函数相关图表，从中提取有价值的信息。

#### 查看Wannier函数中心

- **Wannier函数中心分布图**:
  `wannier90`会生成一个`wannier90_centres.xyz`文件，其中包含了Wannier函数中心的位置。可以使用可视化软件（如VMD或XCrySDen）打开该文件，查看Wannier函数的局域化情况。

#### 查看能带结构

- **紧束缚模型能带结构**:
  `wannier90`会生成一个`wannier90_band.dat`文件，其中包含了紧束缚模型的能带结构。可以使用Matplotlib或其他绘图工具绘制能带图。

#### 查看哈密顿量矩阵

- **哈密顿量矩阵**:
  `wannier90`会生成一个`wannier90_hr.dat`文件，其中包含了哈密顿量矩阵。可以使用Python或其他编程语言读取并分析该文件。

### 附录：常见问题与进阶建议

- **内存消耗问题**:
  当`ecutwfc`值较大时，可能会遇到内存不足的情况。建议尝试降低该参数或增加可用RAM。

- **SCF不收敛解决办法**:
  如果发现自洽场迭代无法正常结束，考虑调整初始电荷密度猜测或者适当提高截断能量。

- **版本兼容性检查**:
  确保所使用的ABACUS和wannier90版本之间存在良好的互操作性；对于特定功能，请参考官方发布的最新文档。

- **进一步学习方向**:
  鼓励读者探索更复杂的模型如磁性体系、表面效应等，同时关注领域内最新的研究成果和技术进展。

通过以上步骤，您可以成功地将ABACUS与wannier90结合使用，生成并分析Wannier函数，从而深入理解材料的电子结构和拓扑性质。

