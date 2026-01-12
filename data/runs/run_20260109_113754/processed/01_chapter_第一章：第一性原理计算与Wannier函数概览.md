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