# PYATB 能带模块介绍

> 来源: https://mp.weixin.qq.com/s/auO0I5gjhPP-lxKY6evhRQ

在标准的计算流程中，ABACUS 负责完成第一性原理的自洽计算，并输出构建紧束缚模型所需的关键矩阵（如哈密顿量、重叠矩阵等）；PYATB 则基于这些结果生成紧束缚模型，进而开展各类物性计算与分析。ABACUS 中的一些常用的计算功能（如能带结构与态密度）可利用PYATB 更便捷地完成。尤其在处理 HSE（杂化泛函）计算结果时，PYATB 可跳过 ABACUS 中复杂的非自洽矩阵重构过程，直接进行高效的处理分析，显著提升整体计算效率。

相较于流行的 Wannier90、WannierTools 等基于最局域 Wannier 函数的紧束缚模型软件，PYATB 无需构造 Wannier 函数，可直接从 ABACUS 提取哈密顿量及相关矩阵，显著简化了计算流程，降低了使用门槛。

目前，PYATB包含多个功能模块，包括：

**能带模块**：计算和分析电子能带结构

**几何模块**：分析能带的贝里曲率及相关几何性质

**光学模块**：计算线性与非线性光学响应

本文作为PYATB系列教程的第一篇，重点介绍**能带模块**的使用方法，包括能带、Fat band、PDOS、以及Spin texture等的计算与可视化。后续文章将继续介绍几何与光学等模块的使用方法，欢迎持续关注。

常用链接以及参考文献：

**GitHub仓库**：https://github.com/pyatb/pyatb.git

**在线使用手册**：https://pyatb.github.io/pyatb/

如使用PYATB软件，请引用以下文章：

Gan Jin, Hongsheng Pang, Yuyang Ji, Zujian Dai, and Lixin He, PYATB: An efficient Python package for electronic structure calculations using ab initio tight-binding model, Comput. Phys. Commun. 291, 108844 (2023).

典型的计算流程概括如下：

1.准备ABACUS的输入文件，通过 ABACUS 进行自洽计算，生成构建紧束缚模型所需的关键矩阵，包括哈密顿量矩阵（HR）、重叠矩阵（SR）以及偶极矩阵（rR）。

2.配置 PYATB 的 Input 文件，并提供 ABACUS 输出的 HR、SR、rR 文件。部分功能（如 Fat band 和 PDOS）还需要 ABACUS 的结构文件（STRU）以及原子轨道基组文件（如 .orb 文件）。

3.运行 PYATB 程序，完成对应功能的计算。在输出文件夹中可以查看各功能函数生成的数据文件及可视化图像结果。

相关示例文件可在此链接获取：https://github.com/pyatb/pyatb/tree/main/examples

在 PYATB 的能带模块中，你可以使用如下功能：

能带结构：能带计算，支持高对称k点路径、Monkhorst-Pack 均匀网格以及自定义 k 点。

能带反折叠：将超胞计算得到的能带结构反折叠回原胞布里渊区，生成对应的谱函数图。用于分析缺陷、掺杂、界面等超胞结构中的电子能态，可直接与 ARPES结果对比。

投影能带（Fat band）：沿高对称k点路径下绘制能带图，并可视化各原子轨道分辨的投影权重。用于分析不同轨道（如 s、p、d）对能带的贡献，揭示能带的轨道成分特征与混合情况。

费米能和费米面：根据体系的价电子数计算费米能，并搜索所有满足费米能的 k 点，绘制费米面。

寻找费米能附近简并点：在指定能量范围内自动识别所有能带简并交叉点。用于寻找 Dirac 点、Weyl 点等拓扑特征，帮助发现新型拓扑材料。

投影态密度（PDOS）：计算不同轨道或原子类型在给定能量下的态密度分布。

自旋纹理（Spin texture）：计算能带中各态的自旋算符期望值，支持、、 分量输出。用于分析自旋轨道耦合效应、Rashba 效应及自旋输运性质，是研究自旋电子学的重要工具。

我们以二维 TMDC 材料 为例，介绍如何使用 ABACUS 进行自洽计算，并生成 PYATB 所需的三类关键矩阵：

哈密顿量矩阵（HR）

重叠矩阵（SR）

偶极矩阵（rR）

完整示例可参考 GitHub 仓库：

https://github.com/pyatb/pyatb/tree/main/examples/WS2

在示例的`abacus/`

文件夹中，已准备好运行所需的 ABACUS 输入文件，包括：

`INPUT`

`STRU`

`*.upf`

`*.orb`

`KPT`

`INPUT`

文件的核心参数如下，已对关键设置作出注释，其他细节可参考 ABACUS 官方文档：https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html`INPUT_PARAMETERS`


# System variables

suffix WS2

ntype 2

calculation scf # 指定计算类型为自洽计算

esolver_type ksdft

symmetry 1

init_chg atomic


# Plane wave related variables

ecutwfc 120


# Electronic structure

basis_type lcao # 使用原子轨道作为基组

ks_solver genelpa

nspin 1

smearing_method gauss

smearing_sigma 0.01

mixing_type pulay

mixing_beta 0.4

mixing_gg0 1.5

scf_nmax 200

scf_thr 1e-8


# Variables related to output information

out_chg 1 # 输出自洽后的电荷密度

out_mat_hs2 1 # 输出HR、SR矩阵

out_mat_r 1 # 输出rR矩阵



执行 ABACUS 自洽计算后，若设置了`out_mat_hs2 = 1`

和`out_mat_r = 1`

，将在`OUT.*`

目录中生成所需矩阵文件：

`data-HR-sparse_SPIN0.csr`

：实空间哈密顿量`data-SR-sparse_SPIN0.csr`

：重叠矩阵`data-rR-sparse.csr`

：偶极矩阵注：若设置

`nspin = 2`

（即自旋极化计算），则还会生成`data-HR-sparse_SPIN1.csr`

文件，用于描述第二个自旋通道。

PYATB 的主控输入文件固定命名为Input，用于设置紧束缚哈密顿量路径、晶格参数以及各计算功能的具体参数。该文件主要由三个区域构成，以 为例，并在注释中说明了关键参数含义。

更多详细说明可参考官方手册：**PYATB用户手册（https://github.com/pyatb/pyatb/blob/main/doc/pyatb.pdf）**

`INPUT_PARAMETERS # 控制参数区域`

{

nspin 1 # 与 ABACUS 的 nspin 保持一致

package ABACUS # 指定数据来源为 ABACUS

fermi_energy 0.3484302262859574

fermi_energy_unit eV

HR_route data-HR-sparse_SPIN0.csr # 哈密顿量矩阵路径

SR_route data-SR-sparse_SPIN0.csr # 重叠矩阵路径

rR_route data-rR-sparse.csr # 偶极矩阵路径（可选）

HR_unit Ry # HR 单位为 Rydberg

rR_unit Bohr # rR 单位为 Bohr

}


LATTICE # 晶格参数区域

{

lattice_constant 1.8897162

lattice_constant_unit Bohr

lattice_vector

3.183820900165 0.0 0.0

-1.591910450082 2.757269780643 0.0

0.0 0.0 20.086904001384

}


BAND_STRUCTURE # 计算功能区域：能带结构

{

kpoint_mode line # 设置 k 点采样模式为高对称路径

kpoint_num 4 # 设置路径中的高对称点数

kpoint_label G, M, K, G # 设置高对称点名称

high_symmetry_kpoint # 设置高对称点坐标及插值数

0.0000000000 0.0000000000 0.0000000000 200

0.5000000000 0.0000000000 0.0000000000 200

0.3333333333 0.3333333333 0.0000000000 200

0.0000000000 0.0000000000 0.0000000000 1

}



PYATB 的参数设置格式与 ABACUS 类似，每个参数块使用`KEYWORD { ... }`

的结构组织。例如：

`KEYWORD`

{

parameter value

}



其中，`INPUT_PARAMETERS`

和`LATTICE`

只能各出现一次；而功能模块（如 `BAND_STRUCTURE`

、`PDOS`

、`FAT_BAND`

等）可以多次定义，实现多功能联合计算。

`HR`

、`SR`

和`rR`

文件是 PYATB 所需的核心输入，由`Input`

文件中的`HR_route`

、`SR_route`

和`rR_route`

参数指定。这些文件通常通过 ABACUS 的自洽计算生成。

提示：对于能带模块计算，

`rR_route`

可省略；而使用几何模块和光学模块等功能时，rR 为必需输入。

ABACUS 自洽计算完成后，可使用 PYATB 自带的辅助程序`pyatb_input`

快速生成`Input`

文件及相关配置文件：

`pyatb_input -h`



该工具支持自动识别 ABACUS 输出结构，快速配置能带、PDOS、Fat band、自旋纹理等计算模块，是推荐的入门方法。

在本节中，我们将以二维材料 为 例，演示如何使用 PYATB 的能带模块完成以下计算任务：

功能 | 物理意义 |
|---|---|

通过这些分析，研究人员可更深入理解材料的电子结构及其物性来源，是材料筛选与机制研究的重要基础。

本节内容基于示例：https://github.com/pyatb/pyatb/tree/main/examples/WS2

`abacus/`

文件夹下运行 ABACUS，生成哈密顿量（HR）、重叠矩阵（SR）和偶极矩阵（rR）等文件。`abacus/`

中执行命令：`pyatb_input -i ./ --band --dim=2 --fatband --bandrange=10 --pdos --erange=5`



参数说明：

`--band`

：启用能带计算；`--fatband`

：启用 Fat band，并设置上下各 10 条能带；`--pdos`

：启用投影态密度计算；`--dim=2`

：指定二维体系；`--erange=5`

：设置能量窗口（±5 eV）。对于Fat band 与 PDOS 功能，还需要轨道文件，因此需要将`W.orb`

和`S.orb`

从`abacus/`

拷贝至`pyatb/。`


使用pyatb_input工具生成的Input文件展示如下：

`INPUT_PARAMETERS`

{

nspin 1

package ABACUS

fermi_energy 0.36725115766

fermi_energy_unit eV

HR_route data-HR-sparse_SPIN0.csr

SR_route data-SR-sparse_SPIN0.csr

rR_route data-rR-sparse.csr

HR_unit Ry

rR_unit Bohr

max_kpoint_num 4000 # 设置内部k点循环计算最大值，避免内存溢出

}


LATTICE

{

lattice_constant 1.0

lattice_constant_unit Angstrom

lattice_vector

3.1838041770607197 0.0 0.0

-1.5919020885298598 2.757255298009863 0.0

0.0 0.0 20.08679849438445

}


BAND_STRUCTURE

{

wf_collect 0 # 是否输出基组表象下的本征波函数

kpoint_mode line

kpoint_num 4

high_symmetry_kpoint

0.0000000000 0.0000000000 0.0000000000 113 # G

0.5000000000 0.0000000000 0.0000000000 65 # M

0.3333333333 0.3333333333 0.0000000000 131 # K

0.0000000000 0.0000000000 0.0000000000 1 # G

kpoint_label G,M,K,G

}


FAT_BAND

{

band_range 3 23 # 设置计算能带范围

stru_file STRU # 设置结构文件名，与ABACUS的结构文件保持一致

kpoint_mode line

kpoint_num 4

high_symmetry_kpoint

0.0000000000 0.0000000000 0.0000000000 113 # G

0.5000000000 0.0000000000 0.0000000000 65 # M

0.3333333333 0.3333333333 0.0000000000 131 # K

0.0000000000 0.0000000000 0.0000000000 1 # G

kpoint_label G,M,K,G

}


PDOS

{

stru_file STRU # 设置结构文件名，与ABACUS的结构文件保持一致

e_range -4.63274884234 5.36725115766 # 设置PDOS的能量范围

de 0.01 # 设置PDOS能量间隔

sigma 0.10 # Gauss smearing参数

kpoint_mode mp # 设置k点模式为MP，即均匀撒点

mp_grid 45 45 1 # 设置MP模式下撒点网格

}



`pyatb/`

，执行：`mpirun -np 8 pyatb`



所有输出将保存至 Out/ 目录，每个功能对应一个子文件夹。

输出文件的说明可以参考**PYATB用户手册。 （https://github.com/pyatb/pyatb/blob/main/doc/pyatb.pdf）**

`Out/Band_Structure/`

`核心文件：`

`band_info.dat`

：`费米能、带隙、CBM/VBM 等基本信息；`

`band.pdf`

：默认绘图输出。`plot_band.py`

：绘图脚本。 对该文件进行修改可以调整能带图。例如进行如下修改，来改变能带窗口范围：`y_min = -5 # eV`

y_max = 5 # eV



`Out/Fat_Band/`

`fatband.xml`

：能带 + 轨道投影信息；`plot_fatband.py`

：绘图脚本。我们需要手动设置需要的参数来进行调整。例如，我们需要绘制W元素的d轨道，可以在脚本中进行如下修改：`# 使用species来绘制fat band。`

species = {"W": [2]}


# 替换atom_index参数，使用species。输出绘图数据到文件中。

pband.write(species=species)


# 替换atom_index参数，species。

pband.plot_contributions(fig, ax, species=species, efermi=efermi, energy_range=energy_range, colors=[])



该 Fat band 图为示意图，如有其他绘图需求，可使用脚本或 Origin 基于原始数据重新绘制。在

`plot_fatband.py`

脚本中，我们执行了`pband.write(species=species)`

函数，因此计算得到的轨道权重数据保存在`Fat_Band/PBAND1_FILE/`

目录下，`band.dat`

文件则包含能带结构数据。通过这些原始文件，可以绘制符合具体要求的自定义图像。

另外，示意图中的X轴没有显示高对称点，可以修改`Fat_Band/high_symmetry_kpoint.dat`

，使用注释的方式增加高对称点的标识。这样就会在X轴显示高对称点。

`K_POINTS`

4

Line

0.000000 0.000000 0.000000 113 # G

0.500000 0.000000 0.000000 65 # M

0.333333 0.333333 0.000000 131 # K

0.000000 0.000000 0.000000 1 # G



`Out/PDOS/`

`PDOS.xml`

：轨道投影信息；`TDOS.dat`

：总态密度（Total DOS）；`plot_dos.py`

：绘图脚本，我们需要手动设置需要的参数来进行调整。例如，我们需要绘制W和S元素的PDOS，可以在脚本中进行如下修改：`# 使用species来绘制PDOS。`

species = ["W", "S"]


# 替换atom_index参数，使用species。输出绘图数据到文件中。

pdos.write(species=species)


# 替换atom_index参数，使用species。

dosplots = pdos.plot(fig, ax, species=species, efermi=efermi, shift=False, energy_range=energy_range, dos_range=dos_range)



与 Fat Band 功能类似，该图也是示意图。如果需要进一步自定义绘图，可以使用脚本或 Origin 基于原始数据进行绘制。在

`plot_dos.py`

脚本中，我们执行了`pdos.write(species=species)`

函数，因此相应的轨道权重数据已保存在`PDOS/PDOS_FILE/`

目录下。

**Spin Texture** 描述的是电子在动量空间中的自旋极化方向分布，尤其在存在自旋轨道耦合（SOC）的体系中具有重要意义。通过自旋纹理图，可以直观展现电子态的自旋结构，是理解拓扑表面态、自旋输运、自旋动量锁定等现象的关键工具。

下面我们以 为例，介绍如何使用 PYATB 的能带模块生成带有自旋纹理的能带图。

相关示例可在以下链接中获取：https://github.com/pyatb/pyatb/tree/main/examples/Bi2Se3

计算流程

`abacus/`

文件夹中执行 ABACUS 计算，生成所需的紧束缚哈密顿量（HR）、重叠矩阵（SR）和位置算符偶极矩阵（rR）。`abacus/`

文件夹下执行命令，自动创建`pyatb/`

文件夹并生成带有自旋纹理功能的`Input`

文件：`pyatb_input -i ./ --spintexture`



说明：

`--spintexture`

：启用`SPIN_TEXTURE`

功能，默认为沿高对称路径绘制能带图。`pyatb/`

文件夹，编辑`Input`

文件中的`SPIN_TEXTURE`

功能区域。例如，可设置一个简化的高对称路径如下：`SPIN_TEXTURE`

{

band_range 73 83 # 设置能带范围（从1开始）

kpoint_mode line

kpoint_num 5

kpoint_label G, Z, F, G, L

high_symmetry_kpoint

0.00000 0.00000 0.0000 100 # G

0.00000 0.00000 0.5000 100 # Z

0.50000 0.50000 0.0000 100 # F

0.00000 0.00000 0.0000 100 # G

0.50000 0.00000 0.0000 1 # L

}



4.执行 PYATB 计算

`mpirun -np 8 pyatb`



`Out/Spin_Texture/`

`spin_texture_x.dat`

： 数据;`spin_texture_y.dat`

： 数据;`spin_texture_z.dat`

： 数据;`plot_spintexture_line.py`

：绘图脚本，用于生成带有自旋期望值的能带图 。`band_with_Sx.png`

`band_with_Sy.png`

`band_with_Sz.png`

以下为 `band_with_Sz.png`

示意图：

图中颜色表示不同 k 点态的自旋极化强度。通过这些图可清晰观察能带中的自旋分布特征，适用于分析 SOC 主导的能带劈裂、自旋动量锁定等现象。