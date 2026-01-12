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