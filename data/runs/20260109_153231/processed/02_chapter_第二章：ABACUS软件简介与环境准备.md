# 第二章：ABACUS软件简介与环境准备

## 2.1 ABACUS概览

ABACUS（Atomic-orbital Based Ab-initio Computation at UStc）是一款强大的第一性原理计算软件，主要用于执行电子结构自洽迭代计算、原子结构优化以及分子动力学模拟等任务。它支持多种类型的材料模拟，并且能够处理复杂的物理问题，如晶格对称性分析、布里渊区的k点对称性分析、电荷密度对称性分析以及力的对称性分析。

ABACUS采用模守恒赝势和周期性边界条件，使得用户能够灵活地选择不同的基组类型进行计算，包括平面波基组（PW）和平面波叠加局域轨道基组（LCAO）。此外，该软件还提供了丰富的功能选项，以满足不同研究需求。

## 2.2 必要软件安装

在开始使用ABACUS之前，我们需要确保已经正确安装了所有必要的软件包。这不仅包括ABACUS本身，还有用于辅助分析的一些Python库。以下是详细的安装步骤：

### 2.2.1 安装ABACUS

首先，从官方GitHub仓库克隆ABACUS源代码并编译安装：
```bash
git clone https://github.com/abacusmodeling/abacus-develop.git
cd abacus-develop
mkdir build && cd build
cmake ..
make -j4
```

完成上述命令后，请将构建好的可执行文件路径添加到系统环境变量中：
```bash
export PATH=$PWD:$PATH
```

### 2.2.2 安装辅助Python库

为了更好地处理输入输出文件及数据分析，建议安装以下Python库：

- **pymatgen**: 用于材料科学的数据分析工具。
- **dpdata**: Deep Potential数据处理库。
- **monty**: 提供一系列实用函数，简化脚本编写过程。
- **numpy**: Python中用于数值计算的基础库。

可以通过pip轻松安装这些库：
```bash
pip install pymatgen dpdata monty numpy
```

### 2.2.3 配置文件准备

使用ABACUS进行计算时，需要准备几个关键文件：`INPUT`, `KPT`, 和 `STRU`。这里提供一个简单的示例来说明每个文件的基本结构及其参数意义。

#### INPUT 文件示例
```plaintext
INPUT_PARAMETERS
suffix       Mn           # 输出文件名后缀
ntype        1            # 原子种类数量
ecutwfc      20           # 平面波截止能量 (单位: Ry)
scf_thr      1e-7         # 自洽场收敛阈值
basis_type   pw           # 基组类型 (pw 或 lcao)
calculation  scf          # 计算类型 (scf, nscf, md 等)
```

- `suffix`: 指定输出文件名的后缀部分。
- `ntype`: 系统中存在的不同元素种类数目。
- `ecutwfc`: 设置平面波展开的能量截止值。
- `scf_thr`: 控制自洽循环停止的标准。
- `basis_type`: 选择使用的基组形式。
- `calculation`: 指明本次运行所执行的具体计算任务。

#### STRU 文件示例
```plaintext
ATOMIC_SPECIES
Mn 54.938045 Mn.upf
LATTICE_CONSTANT
1.8897259886
LATTICE_VECTORS
  4.0  0.0  0.0
  0.0  4.0  0.0
  0.0  0.0  4.0
ATOMIC_POSITIONS
Direct
  0.0  0.0  0.0
```

- `ATOMIC_SPECIES`: 列出所有参与计算的元素信息。
- `LATTICE_CONSTANT`: 设定晶格常数，默认单位为Bohr。
- `LATTICE_VECTORS`: 定义晶体的三个基本向量。
- `ATOMIC_POSITIONS`: 给出原子在晶胞中的坐标位置。

请注意，在实际应用中可能还需要根据具体情况进行调整，例如增加额外的赝势或数值原子轨道文件路径等。

### 2.2.4 结果验证

最后但同样重要的是，务必检查计算结果的合理性。特别是对于应力张量这类敏感指标，其准确性直接影响后续分析结论。如果发现任何异常情况，请参考[ABACUS官方文档](https://abacus.deepmodeling.com/)以获取更详细的指导信息。

通过遵循以上步骤，您应该能够成功设置好ABACUS的工作环境，并准备好开展您的第一次计算实验了！