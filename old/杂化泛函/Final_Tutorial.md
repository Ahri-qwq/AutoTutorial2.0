# ABACUS 实战教程：ABACUS原子轨道基组杂化泛函计算

# 第一章：密度泛函理论与杂化泛函概述

# 第一章：密度泛函理论与杂化泛函概述

本章将为读者提供必要的理论背景，帮助理解为什么需要使用杂化泛函以及它们如何改进了传统DFT方法的局限性。我们将从密度泛函理论的基础概念开始，逐步深入到杂化泛函的概念及其优势，并简要介绍LCAO方法在执行杂化泛函计算时的特点。

## Section 1.1: 密度泛函理论简介

### 基础概念
密度泛函理论（Density Functional Theory, DFT）是目前最广泛使用的量子力学方法之一，用于研究物质的电子结构。DFT的核心思想是通过电子密度而非波函数来描述多电子系统，从而大大简化了问题的复杂度。DFT的理论基础主要包括Hohenberg-Kohn定理和Kohn-Sham方程。

- **Hohenberg-Kohn定理**：该定理表明，对于一个给定的外部势场，系统的基态能量仅依赖于电子密度，而与其他变量无关。这意味着，我们可以通过优化电子密度来找到系统的基态能量。
- **Kohn-Sham方程**：为了求解电子密度，引入了一组无相互作用的单电子轨道（Kohn-Sham轨道），这些轨道满足薛定谔方程形式的方程。通过求解这些方程，我们可以得到电子密度，进而计算出系统的总能量。

### 关键参数
虽然本节不涉及具体的`INPUT`文件参数，但理解这些基本概念对于后续章节中选择合适的`dft_functional`等参数至关重要。

## Section 1.2: 杂化泛函的概念与优势

### 杂化泛函的概念
杂化泛函（Hybrid Functionals）是在传统的交换关联泛函中引入一部分Hartree-Fock（HF）交换能，以改善对电子交换相互作用的描述。这种混合方式可以显著提高带隙预测、表面性质和吸附行为等方面的计算精度。

- **Hartree-Fock交换能**：HF交换能是一种精确的交换能形式，它考虑了电子之间的库仑排斥力。然而，纯HF方法忽略了相关能，导致其在某些情况下表现不佳。
- **杂化比例**：通过调整HF交换能的比例（通常用`exx_hybrid_alpha`参数表示），可以在保持HF优点的同时，引入局域或半局域泛函的相关能部分，从而达到更好的平衡。

### 关键参数
- **`dft_functional`**：指定所使用的杂化泛函类型，如`pbe0`, `hse`, `scan0`等。
- **`exx_hybrid_alpha`**：定义HF交换能的比例，默认值取决于具体泛函，例如PBE0为0.25。
- **`exx_hse_omega`**：HSE泛函中的范围分离参数，默认值为0.11。

### 编译注意事项
- **LibXC支持**：确保ABACUS编译时链接了LibXC库，这对于支持高级功能如HSE泛函至关重要。
- **实验性功能**：平面波下的杂化泛函是实验性功能，尚未经过全面测试，请参考官方文档获取最新信息。

## Section 1.3: LCAO方法简述

### LCAO方法
线性组合原子轨道（Linear Combination of Atomic Orbitals, LCAO）方法是一种常用的电子结构计算方法。它通过将原子轨道线性组合来构造分子轨道，从而简化了计算过程。相比于平面波方法，LCAO方法在处理大体系和复杂材料时具有更高的效率。

### 关键参数
- **`basis_type`**：设置为`lcao`，表示使用数值原子轨道基组。

### 计算策略
- **内存管理**：对于大规模计算任务，建议评估并适当调整内存分配策略，以避免因资源不足导致的任务失败。

### 示例代码
以下是一个简单的`INPUT`文件示例，展示了如何设置杂化泛函和LCAO方法：

```plaintext
&GLOBAL
  run_type = "scf"
  out_level = 1
  basis_type = "lcao"
/
&ELECTRON
  dft_functional = "pbe0"
  exx_hybrid_alpha = 0.25
  ks_solver = "genelpa"
  nband = 30
  smearing_method = "gaussian"
  smearing_sigma = 0.01
/
&FORCE_EVAL
  stress = .TRUE.
  do_force = .TRUE.
/
```

### 总结
本章介绍了密度泛函理论的基础概念，解释了杂化泛函的优势及其关键参数，并简要讨论了LCAO方法的特点。通过这些内容，读者将能够更好地理解和应用ABACUS软件中的相关功能，进行更准确的材料模拟和计算。

# 第二章：ABACUS软件环境搭建与基本配置

## 第二章：ABACUS软件环境搭建与基本配置

在掌握了必要的理论知识后，本章将指导用户完成ABACUS软件的安装，并设置好进行杂化泛函计算所需的基本环境。我们将详细说明如何下载源码并根据不同的操作系统配置编译选项，特别是关于是否链接LibXC库的重要性。此外，我们还将深入讲解输入文件的主要组成部分，以及如何构建或导入初始晶体结构模型。

### Section 2.1: 安装与编译

#### 下载源码
首先，你需要从ABACUS的官方GitHub仓库下载最新的源代码。你可以通过以下命令克隆仓库：

```bash
git clone https://github.com/abacusmodeling/abacus-develop.git
cd abacus-develop
```

#### 配置编译选项
在开始编译之前，你需要根据你的系统和需求配置编译选项。这里特别强调的是，是否链接LibXC库对于支持某些高级功能（如HSE泛函）至关重要。如果你需要使用这些高级功能，请确保在编译时指定了`--with-libxc`选项。

##### 编译步骤
1. **安装依赖**：确保你的系统已经安装了所有必要的依赖库。对于大多数Linux发行版，你可以使用包管理器来安装这些依赖。例如，在Ubuntu上，你可以运行以下命令：

   ```bash
   sudo apt-get update
   sudo apt-get install build-essential cmake libfftw3-dev libhdf5-dev libopenmpi-dev
   ```

2. **配置CMake**：使用CMake生成Makefile。你可以通过以下命令来配置CMake：

   ```bash
   mkdir build
   cd build
   cmake .. -DWITH_LIBXC=ON
   ```

   这里，`-DWITH_LIBXC=ON`选项指定了要链接LibXC库。如果你不需要LibXC库，可以省略这个选项。

3. **编译**：运行`make`命令来编译ABACUS：

   ```bash
   make -j $(nproc)
   ```

   这里的`-j $(nproc)`选项表示使用所有可用的核心来加速编译过程。

4. **安装**：编译完成后，你可以将生成的可执行文件安装到系统路径中：

   ```bash
   sudo make install
   ```

#### 检查安装
安装完成后，你可以通过以下命令来检查ABACUS是否正确安装：

```bash
which abacus
abacus --version
```

如果一切正常，你应该能看到ABACUS的版本信息。

### Section 2.2: 输入文件结构解析

#### INPUT文件
ABACUS的输入文件主要是`INPUT`文件，它包含了计算所需的大部分参数。以下是`INPUT`文件的一些关键部分及其物理意义：

```plaintext
INPUT_PARAMETERS
suffix       Mn           # 输出后缀
ntype        1            # 元素种类数
ecutwfc      20           # 平面波展开截止能量 (Ry)
scf_thr      1e-7         # 电荷密度收敛阈值
basis_type   pw           # 基函数类型 (pw 或 lcao)
calculation  scf          # 计算类型 (scf, nscf, relax, md 等)
ks_solver    cg           # Kohn-Sham 方程求解器 (cg, hf, etc.)
dft_functional PBE         # DFT 泛函 (PBE, HSE, etc.)
exx_ccp_rmesh_times 1.5   # EXX 截断半径倍数
```

- `ecutwfc`：平面波展开的截止能量，单位为Ry。这个参数决定了计算的精度，较高的值会带来更高的精度但也会增加计算成本。
- `scf_thr`：自洽场迭代的收敛阈值。较小的值意味着更严格的收敛标准。
- `basis_type`：基函数类型，可以选择`pw`（平面波）或`lcao`（数值原子轨道线性组合）。
- `ks_solver`：Kohn-Sham方程的求解器，`cg`表示共轭梯度法，`hf`表示哈特里-福克方法。
- `dft_functional`：DFT泛函，常用的有PBE、HSE等。
- `exx_ccp_rmesh_times`：EXX截断半径倍数，用于控制EXX计算的精度。

#### STRU文件
`STRU`文件包含了晶体结构的信息，包括原子种类、位置、晶格常数和晶格向量等。

```plaintext
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf
LATTICE_CONSTANT
1.8897259886
LATTICE_VECTORS
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
ATOM_POSITIONS
Direct
0.0 0.0 0.0
```

- `ATOMIC_SPECIES`：定义了原子种类、相对原子质量和赝势文件名。
- `LATTICE_CONSTANT`：晶格常数，默认单位为Bohr。
- `LATTICE_VECTORS`：晶格向量，定义了晶胞的形状。
- `ATOM_POSITIONS`：原子位置，可以是直接坐标或笛卡尔坐标。

### Section 2.3: 初始结构准备

#### 构建初始结构
你可以使用各种工具来构建或导入初始晶体结构模型。常见的工具包括ASE、VESTA等。以下是一个简单的例子，展示如何使用ASE来创建一个立方体晶胞：

```python
from ase.build import bulk
from ase.io import write

# 创建一个立方体晶胞
atoms = bulk('U', 'sc', a=5.468, cubic=True)

# 将结构写入STRU文件
write('STRU', atoms, format='stru')
```

#### 调整磁性设置
如果你需要考虑磁性，可以在`STRU`文件中添加磁矩信息。例如：

```plaintext
MAGNETIC_MOMENTS
1.0
```

这表示每个原子有一个磁矩为1.0 μB。

#### 创建超胞
如果你需要创建超胞，可以使用ASE的`repeat`方法：

```python
from ase.build import bulk
from ase.io import write

# 创建一个立方体晶胞
atoms = bulk('U', 'sc', a=5.468, cubic=True)

# 创建2x2x2的超胞
supercell = atoms * (2, 2, 2)

# 将超胞写入STRU文件
write('STRU', supercell, format='stru')
```

### 总结
通过本章的学习，你应该能够成功安装和配置ABACUS软件，并准备好进行杂化泛函计算所需的基本环境。请务必检查官方文档获取最新信息，并注意实验性功能可能存在未预见的问题。对于大规模计算任务，建议评估并适当调整内存分配策略以避免因资源不足导致的任务失败。

# 第三章：执行杂化泛函计算

# 第三章：执行杂化泛函计算

在本章中，我们将详细介绍如何使用 ABACUS 软件执行杂化泛函计算。通过前两章的学习，我们已经掌握了基本的输入文件设置和计算流程。现在，我们将进一步探讨如何优化 SCF 循环、提高计算效率，并分析输出结果。

## Section 3.1: SCF循环优化

### 内容
SCF（自洽场）迭代是 DFT 计算的核心步骤，尤其在杂化泛函计算中，SCF 的收敛性往往更加困难。本节将讨论如何通过调整各种参数来加速 SCF 迭代过程，避免收敛困难。

### 关键参数
- `scf_nmax`: 控制 SCF 迭代的最大次数。
- `smearing_method`: 用于控制电子占据状态的展宽方法。
- `mixing_beta`: 用于控制新旧密度矩阵混合的比例。

#### 示例 INPUT 文件片段
```plaintext
# Parameters (SCF)
scf_nmax            100
smearing_method     gaussian
smearing_sigma      0.05
mixing_beta         0.7
```

### 参数解释
- **`scf_nmax`**: 设置 SCF 迭代的最大次数。默认值通常是 100，但在复杂体系或难以收敛的情况下，可能需要增加该值。
- **`smearing_method`**: 选择展宽方法，常见的有 `gaussian` 和 `methfessel-paxton`。展宽方法可以平滑费米面附近的能级，有助于 SCF 收敛。
- **`mixing_beta`**: 混合比例因子，用于控制新旧密度矩阵的混合程度。较小的值会使得新的密度矩阵对旧的依赖更少，但过小可能导致不稳定。通常建议从 0.7 开始尝试。

### 实用技巧
- 如果 SCF 收敛困难，可以尝试减小 `mixing_beta` 值，或者增加 `scf_nmax`。
- 对于杂化泛函，特别是 HSE 和 PBE0，建议使用 `gaussian` 展宽方法，并适当调整 `smearing_sigma`。

## Section 3.2: 高级特性与性能调优

### 内容
本节将介绍一些更高级的功能设置，包括内存管理策略、并行计算优化等，以提高计算效率。

### 关键参数
- `exx_ccp_rmesh_times`: 控制库仑势计算所需的径向格点数。
- 环境变量配置: 如 `OMP_NUM_THREADS` 和 `MKL_NUM_THREADS`。

#### 示例 INPUT 文件片段
```plaintext
# Parameters (EXX)
exx_ccp_rmesh_times  1.5
```

### 参数解释
- **`exx_ccp_rmesh_times`**: 控制库仑势计算所需的径向格点数。对于 HSE 泛函，默认值为 1.5，对于 PBE0/HF 默认值为 5。减小该值可以提高计算速度，但可能影响精度。

### 性能调优
- **内存管理**: 杂化泛函计算通常需要大量的内存。可以通过调整 `exx_ccp_rmesh_times` 来减少内存消耗。
- **并行计算**: 使用 OpenMP 和 MPI 并行化可以显著提高计算效率。确保合理设置环境变量如 `OMP_NUM_THREADS` 和 `MKL_NUM_THREADS`。

### 实用技巧
- 对于大规模计算任务，建议评估并适当调整内存分配策略以避免因资源不足导致的任务失败。
- 使用多节点并行时，注意 MPI 和 OpenMP 的混合并行设置，以充分利用计算资源。

## Section 3.3: 结果分析与验证

### 内容
本节将学习如何解读输出结果，并对比实验数据或其他计算方法的结果来进行验证。

### 关键参数
- 输出文件格式说明: 主要关注 `OUT.ABACUS/running_scf.log` 和 `OUT.ABACUS/STRU_ION_D` 文件。

#### 示例 OUTPUT 文件片段
```plaintext
# OUT.ABACUS/running_scf.log
SCF iteration 1, energy = -1234.56789 eV
SCF iteration 2, energy = -1234.56790 eV
...
SCF converged in 10 iterations
```

### 结果分析
- **能量收敛**: 检查 `running_scf.log` 文件中的能量变化，确保 SCF 迭代已经收敛。
- **结构优化**: 如果进行了结构优化，检查 `STRU_ION_D` 文件中的原子坐标变化，确保结构已经稳定。

### 实用技巧
- 对比实验数据或其他计算方法的结果，验证计算的准确性。
- 使用可视化工具如 VESTA 或 ASE 来查看结构和电子密度分布。

## 附录：常见问题与进阶建议

- **内存消耗过高**:
  - 减小体系规模或降低精度要求。
  - 调整 `exx_ccp_rmesh_times` 参数以减少内存消耗。
- **难以收敛**:
  - 修改 SCF 相关参数如 `scf_nmax` 和 `mixing_beta`。
  - 尝试不同的初始猜测。
- **特定泛函的使用**:
  - 在使用特定泛函之前，请确保其已被充分测试，并参考官方文档获取最新信息。
- **大规模计算任务**:
  - 合理规划资源分配，确保有足够的 CPU 核心和内存。
  - 探索更多功能，如使用 Python 脚本自动化工作流程，或与其他工具（如 Phonopy, ASE）结合使用以扩展研究范围。

通过本章的学习，你将能够有效地执行和优化杂化泛函计算，提高计算效率，并准确分析和验证计算结果。

