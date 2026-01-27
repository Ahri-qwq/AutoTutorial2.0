# ABACUS入门教程 - 结构文件STRU

在使用ABACUS开展计算之前，我们需要准备好ABACUS计算用的输入文件，这些输入文件主要包括：

STRU文件：ABACUS的结构文件格式，可类比VASP的POSCAR。

INPUT文件：开展ABACUS计算的输入参数集合，可类比VASP的INCAR。

KPT文件：ABACUS的K点文件格式，代表了开展第一性原理计算时的K空间采样情况，可类比VASP的KPOINTS。

其中，结构文件STRU的准备是ABACUS输入文件准备的重要一环。一般来说，用户的结构文件来源主要包括：

自行通过Materials Studio Visualizer等晶体结构建模可视化软件构建结构模型，并保存为这一建模软件所支持的格式。

从Materials Project等晶体结构数据库上下载的结构文件（多为CIF格式）。

其他来源，如求解晶体表征谱图得到的CIF文件，从文献中获取的结构文件等。

这些结构文件的格式大多为CIF，POSCAR，XSD等，不会直接是ABACUS的STRU文件。那么，我们该如何准备ABACUS的结构文件呢？

与此同时，在ABACUS计算完成之后，我们有没有将计算结果转化为可被常见的晶体结构建模可视化软件导入并可视化的常见结构文件格式（如CIF）的方法？

ASE（Atomic Simulation Environment）是开源的原子尺度模拟Python工具库，由丹麦技术大学（DTU）等机构开发，广泛应用于计算材料科学和分子模拟领域。该工具库支持多种第一性原理软件，ABACUS也有ASE-ABACUS接口，不过这一接口是独立于ASE主仓库以外的分支，需要单独下载安装。安装方式：

`# install ase-abacus`

git clone https://gitlab.com/1041176461/ase-abacus.git

cd ase-abacus

pip install .

# use tsinghua mirror in CN internet if needed

# pip install . -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple


安装好了之后就可以使用这个ASE来准备STRU文件了。有两种可行方式，一种是直接基于ASE构建结构文件并导出为STRU格式。比如以下Python脚本（或IPython命令行）：

`from ase.build import bulk`

from ase.io import write


structure = bulk("Al", "fcc")

write("STRU", structure, format='abacus')


如此即可得到FCC Al的一个最简单的STRU文件，它是FCC Al的原胞。

注：通过pip install ase可以下载到PyPi上官方的ASE，但这个ASE里没有接入ABACUS的部分支持。ASE-ABACUS是ASE的一个派生分支，专用于ASE的ABACUS支持，同时兼容ASE本体。


通过上述最简单的方式准备的STRU文件内容如下：

`ATOMIC_SPECIES`

Al 26.9815385 Al_ONCV_PBE-1.0.upf


NUMERICAL_ORBITAL

Al_gga_7au_100Ry_4s4p1d.orb


LATTICE_CONSTANT

1.8897261258369282


LATTICE_VECTORS

0.0000000000 2.0250000000 2.0250000000

2.0250000000 0.0000000000 2.0250000000

2.0250000000 2.0250000000 0.0000000000


ATOMIC_POSITIONS

Direct


Al

0.0000000000

1

0.0000000000 0.0000000000 0.0000000000 1 1 1 mag 0.0



元素类别信息。包括元素符号，摩尔质量，赝势文件名称。

原子轨道文件名称。在使用原子轨道基组计算时使用。

晶格常数单位。ABACUS的默认长度单位为Bohr，这一晶格常数将所有长度放大1.889726倍，将长度单位Bohr转化为更为人熟知的埃（Angstrom）。

晶格矢量，为九分量形式。

原子坐标类别：有多种选项，可以在docs: atomic-positions（https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html#atomic-positions）处查看，常见选项包括分数坐标Direct，笛卡尔坐标Cartesian。

注：还有个Cartesian_angstrom选项可以避免对晶格坐标单位LATTICE_CONSTANT的显式指定，直接以埃为单位，不过这一选项在各个接口软件中支持的不多，也不广为人知。


`Al #元素类别`

0.0000000000 # 该元素每个原子的磁矩（低优先级）

1 # 该元素包含的原子数目

0.0000000000 0.0000000000 0.0000000000 1 1 1 mag 0.0

# 坐标x,y,z; 原子在x,y,z方向是否可动（1代表可动），原子磁矩（高优先级）


除了基于ASE构建的结构准备STRU外，更常见的是，用户已经有自己想做计算的输入结构，只是要把它转化为STRU格式，这一点ASE-ABACUS当然也支持，以如下FCC Al的CIF文件为例：

`data_image0`

_chemical_formula_structural Al

_chemical_formula_sum "Al1"

_cell_length_a 2.8637824638055176

_cell_length_b 2.8637824638055176

_cell_length_c 2.8637824638055176

_cell_angle_alpha 60.00000000000001

_cell_angle_beta 60.00000000000001

_cell_angle_gamma 60.00000000000001


_space_group_name_H-M_alt "P 1"

_space_group_IT_number 1


loop_

_space_group_symop_operation_xyz

'x, y, z'


loop_

_atom_site_type_symbol

_atom_site_label

_atom_site_symmetry_multiplicity

_atom_site_fract_x

_atom_site_fract_y

_atom_site_fract_z

_atom_site_occupancy

Al Al1 1.0 0.0 0.0 0.0 1.0000



比较简单的Python操作方式如下：

`from ase.io import read, write`

# read existing cif file

Al_fcc = read("Al.cif", format="cif")

# write this structure to STRU

write("STRU", Al_fcc, format='abacus')


如果想自行指定输出的赝势文件和轨道文件，可以这么写你的Python脚本：

`from ase.io import read, write`

Al_fcc = read("Al.cif", format="cif")

# specify pseudopot from APNS-v1.0

pp = {"Al": "Al.upf"}

# specify orbital from APNS-v1.0

basis = {"Al": "Al_gga_8au_100Ry_2s2p1d.orb"}

write("STRU", Al_fcc, format='abacus', pp=pp, basis=basis)



`# Set up the environment for ABACUS`

PP=${HOME}/pseudopotentials # 赝势文件路径

ORB=${HOME}/orbitals # 轨道文件路径

export ABACUS_PP_PATH=${PP} # 设置环境变量

export ABACUS_ORBITAL_PATH=${ORB} # 设置环境变量


ATOMKIT由VASPKIT开发团队开发，提供了另一种选项：通过命令行交互方式，集成晶体结构操作、建模、格式转化及可视化，并可输出多种常见材料模拟软件结构文件格式。同时，ATOMKIT兼具一定的可视化功能。

其最新版可在此处下载：

下载好了之后，直接解压，然后安装即可，如：

`# under target directory`

tar -zxvf atomkit.0.9.0.linux.x64.tar.gz

cd atomkit.0.9.0.linux.x64

bash setup.sh

# reload your system environment

source ~/.bashrc



` ~> atomkit`

\\\///

/ _ _ \

(| (o)(o) |)

o-----.OOOo--()--oOOO.------------------------------------------o

| ATOMKIT Standard Edition 0.9.0 (26 Mar. 2024) |

| Lead Developer: Vei WANG (wangvei@icloud.com) |

| Contributors: Zhao-Ke WANG, Gang TANG, Ya-Chao LIU |

o-----.oooO-----------------------------------------------------o

( ) Oooo. ATOMKIT Made Simple

\ ( ( )

\_) ) /

(_/

===================== Structural Utilities ======================

1) Export Structures 2) Symmetry Analysis

3) Generate K-Mesh & K-Path 4) Edit Structure

5) Advanced Structure Models 6) 2D-Material Kit

8) Visualize Structure 0) Quit

------------>>



`1`

+---------------------------- Tip ------------------------------+

| Supplementary parameters from GLOBAL_PARAMETERS.in will be |

| automatically copied to the exported file, if available. |

+---------------------------------------------------------------+

===================== Export Format Options =====================

101) ABACUS (*.stru) 102) ABINIT (*.in)

105) ASE (*.py) 107) ATAT (*.in)

109) AtomEye (*.cfg) 112) CASTEP (*.cell)

113) CIF Format (*.cif) 114) CP2K (*.inp)

115) CPMD (*.inp) 123) Elk (*.in)

124) exciting (*.xml) 126) DFTB+ (*.hsd)

127) DMol3 (*.incoor) 129) FHI-aims (*.in)

132) FLEUR (*.inp) 134) Gaussian (*.gjf)

136) LAMMPS (*.lmp) 138) JDFTx (*.in)

142) Octopus (*.in) 144) OpenMX (*.in)

146) ORCA (*.in) 149) PDB Format (*.pdb)

156) Quantum-Espresso (*.in) 165) RSPt (*.inp)

168) Siesta (*.fdf) 175) VASP (*.vasp)

180) Wannier90 (*.win) 182) WIEN2K (*.struct)

185) XCrySDen (*.xsf) 189) XYZ Format (*.xyz)


0) Quit

9) Back

------------>>



`101`


===================== Import Format Options =====================

101) ABACUS (STRU) 112) CASTEP (*.cell)

113) CIF Format (*.cif) 123) Elk (*.in)

134) Gaussian (*.gjf) 142) Materials Studio (*.xsd)

144) MDL Molfile Format (*.mol) 149) PDB Format (*.pdb)

156) Quantum-Espresso (*.in) 168) Siesta (*.fdf)

175) VASP (POSCAR, *.vasp) 185) XCrySDen (*.xsf)

189) XYZ Format (*.xyz)

190) Build Crystal from Wyckoff Positions & Space Group Number


+---------------------------- Tip ------------------------------+

| Examples of Usage: |

| 113 Si.cif: Read structure from Si.cif file with cif format |

| 175 POSCAR: Read structure from POSCAR file with VASP format |

| 190 INPUT: Build crystal From Wyckoff positions in INPUT file|

+---------------------------------------------------------------+

Input the code-ID listed above and the filename to be read:



`113 Al.cif`

-->> (01) Reading Structure from Al.cif File...

+---------------------------------------------------------------+

| This utility is in experimental stage & needs to be improved. |

+---------------------------------------------------------------+

+---------------------------------------------------------------+

Species and Wyckoff Sites

[Al] 0.000 0.000 0.000

+---------------------------------------------------------------+

| Cell Parameters After Applying Symmetry Operations |

+---------------------------------------------------------------+

Import File Format: CIF

Lattice Constants: 2.864 2.864 2.864

Lattice Angles: 60.000 60.000 60.000

Total Atoms: 1

Number of Atoms: [Al] 1

+---------------------------------------------------------------+

-->> (02) Written Al.STRU File!

o---------------------------------------------------------------o

| WARNING: Your are using a beta version of ATOMKIT |

| Please Upgrade to the latest stable version if available |

| Any Suggestions for Improvement are Welcome and Appreciated |

o---------------------------------------------------------------o



可以看到ATOMKIT自动读取了Al.cif文件，并将其转化为了ABACUS的STRU格式文件Al.STRU，这一文件内容为：

`ATOMIC_SPECIES`

Al 26.982 Al_ONCV_PBE-1.0.upf


NUMERICAL_ORBITAL

Al_gga_7au_100Ry_4s4p1d.orb


LATTICE_CONSTANT

1.889726


LATTICE_VECTORS

2.863782463806 0.000000000000 0.000000000000

1.431891231903 2.480108364568 0.000000000000

1.431891231903 0.826702788189 2.338268590218


ATOMIC_POSITIONS

Direct


Al

0.000

1

0.000000000000 0.000000000000 0.000000000000 1 1 1 mag 0.0


读者可自行思考为什么LATTICE_VECTORS看起来和之前的STRU文件不同。这里能给出的结论是二者等价。

`from ase.io import read, write`

# read existing STRU file

structure = read("STRU", format="abacus")

# write this structure to CIF

write("my_stru.cif", stru, format='cif')



`1`

113

101 STRU



然后通过如下方式：

`atomkit < example.txt`



除去上述两种方式以外，还有一个工具abacustest也能用来准备ABACUS的STRU等输入文件。

abacustest（https://github.com/pxlxingliang/abacus-test）用于对ABACUS进行前后处理，以及高通量任务的计算。目前abacustest有两种安装方式：通过pip从PyPi源安装，或者是git拉取仓库后源码安装：

`# install via pip`

pip install abacustest

# install via git+pip

git clone https://github.com/pxlxingliang/abacus-test.git

cd abacus-test

pip install .



就快速准备单个或几个ABACUS任务的输入文件而言，可以使用abacustest model inputs功能。使用abacustest model inputs -h可以查看准备ABACUS输入文件的命令：

`>>> abacustest model inputs -h`


------------------------------------------------------------

+ ABACUSTEST

+ version: v0.4.20

+ GITHUB: https://github.com/pxlxingliang/abacus-test/tree/develop

+ BohriumApp: https://bohrium.dp.tech/apps/abacustest

------------------------------------------------------------


usage: abacustest model inputs [-h] [-f [FILE ...]] [--ftype FTYPE] [--jtype JTYPE] [--pp PP] [--orb ORB] [--input INPUT] [--kpt [KPT ...]] {prepare,post} ...


Prepare the ABACUS inputs file.


positional arguments:

{prepare,post}

prepare Prepare the model

post Post-process the model


options:

-h, --help show this help message and exit

-f [FILE ...], --file [FILE ...]

the structure files

--ftype FTYPE the structure type, should be cif or dpdata supportted type like: poscar, abacus/stru, ..

--jtype JTYPE the job type, should be one of ['scf', 'relax', 'cell-relax', 'md', 'band']

--pp PP the path of pseudopotential library, or read from enviroment variable ABACUS_PP_PATH

--orb ORB the path of orbital library, or read from enviroment variable ABACUS_ORB_PATH

--input INPUT the template of input file, if not specified, the default input will be generated

--kpt [KPT ...] the kpoint setting, should be one or three integers



1.结构 (-f)：支持CIF文件，POSCAR，或直接使用ABACUS的STRU。

2.赝势库 (-pp): ABACUS 计算所需的赝势库路径。

3.轨道库 (-orb): ABACUS 计算所需的轨道库路径。

注：赝势库和轨道库的文件名需以元素名开头，或者在其中准备一个element.json文件，文件内容为{元素名：赝势文件名}的dict。如此程序才能自动查找元素对应的赝势和轨道。如果在赝势库的目录下同时提供了文件ecutwfc.json（其内容为每种元素推荐的ecutwfc的值，格式与element.json文件相同），abacustest会自动检查体系的元素，并自动设置ecutwfc为体系所有元素的最大值。

之后，用户可以用如下命令：

`abacustest model inputs -f 1.cif 2.cif \`

--ftype cif \

--pp /path/to/your/pp_dir \

--orb /path/to/your/orb_dir



同时，abacustest也支持基于准备好的设置，批量准备输入任务，用户需要准备类似于下面的param.json文件：

`{`

"prepare":{

"strus": ["1.vasp", "2.vasp"],

"stru_format": "poscar",

"input_template": "INPUT",

"kpt_template": "KPT",

"pp_dict": {"H": "H.upf", "O": "O.upf"},

"orb_dict": {"H": "H.orb", "O": "O.orb"},

"pp_path": "path/to/pp",

"orb_path": "path/to/orb"

}

}



`abacustest prepare -p param.json -s abacustest`



# ABACUS 使用教程｜如何转换 STRU 文件

经常有同学问 CIF 文件怎么转 STRU, POSCAR 怎么转 STRU，以及反过来又怎么做？这个文档就来回答这个问题。

总的回答是：我们推荐使用Python ASE-ABACUS 接口来做。具体的操作步骤和例子如下：

## 1 安装 python ase-abacus 模块

（在本篇 Notebook 镜像中已安装，可跳过）

```
%%bash
# 安装 ase-abacus
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
python3 setup.py install
```

## 2 设置环境变量（可选）

[ABACUS](http://abacus.ustc.edu.cn/main.htm)支持两种基组PW, LCAO。赝势和轨道文件的存放路径可以通过环境变量设置，分别为：ABACUS_PP_PATH 和 ABACUS_ORBITAL_PATH,设置方法如下：

```
%%bash
# Set up the environment for ABACUS
PP=${HOME}/pseudopotentials  # 赝势文件路径
ORB=${HOME}/orbitals  # 轨道文件路径
export ABACUS_PP_PATH=${PP}  # 设置环境变量
export ABACUS_ORBITAL_PATH=${ORB}  # 设置环境变量
```


PW 赝势计算只需要设置 ABACUS_PP_PATH 。 LCAO 赝势需要两个都设置：ABACUS_PP_PATH 和 ABACUS_ORBITAL_PATH 。

## 3 加入 ABACUS 计算器的方法

在 Python 脚本中引入 ABACUS 接口的方式为：

```python
from ase.calculators.abacus import Abacus
```

## 4 STRU 格式转化例子

在示例数据集中，我们为你准备好了格式转换的示例文件。我们可以直接访问：


```python
! tree /bohr/abacus-stru-q28v/v2/
```

    [01;34m/bohr/abacus-stru-q28v/v2/[0m
    └── [01;34mABACUS_STRU[0m
        ├── [01;34mCIF2STRU[0m
        │   └── [01;32mAMS_DATA-92.cif[0m
        ├── [01;34mPOSCAR2STRU[0m
        │   └── [01;32mPOSCAR[0m
        ├── [01;34mSTRU2CIF[0m
        │   └── [01;32mSTRU[0m
        └── [01;34mSTRU2POSCAR[0m
            └── [01;32mSTRU[0m
    
    5 directories, 4 files
    

出于安全考虑，我们没有数据集所在文件夹的写入权限，因此我们将其复制到 `/data/` 目录下:


```python
! cp -nr /bohr/ /data/
```

### 4.1 CIF 转 STRU

*.cif 是常见的晶体结构文件，这里提供了一种 cif 文件转 STRU 文件的方法：


```python
from ase.io import read, write
import os

parent_path = os.path.join("/", "data", "bohr", "abacus-stru-q28v", "v2", "ABACUS_STRU")

cs_dir = os.path.join(parent_path, "CIF2STRU")

# CIF文件路径
cs_vasp = os.path.join(cs_dir, "AMS_DATA-92.cif")

# 读取CIF文件
cs_atoms = read(cs_vasp, format="cif")

# 输出的STRU文件路径
cs_stru = os.path.join(cs_dir, "STRU")

# 赝势文件路径
pp = {"Si": "Si.PD04.PBE.UPF", "O": "O.upf"}

# 基组文件路径
basis = {"Si": "Si_gga_10au_100Ry_3s3p2d.orb", "O": "O.orb"}  

# 输出STRU文件
write(cs_stru, cs_atoms, format="abacus", pp=pp, basis=basis)  
```

可以看到，我们的 CIF2STRU 文件夹下除了 cif 文件外，生成了相应的 STRU 文件：


```python
! tree /data/bohr/abacus-stru-q28v/v2/ABACUS_STRU
```

    [34;42m/data/bohr/abacus-stru-q28v/v2/ABACUS_STRU[0m
    ├── [34;42mCIF2STRU[0m
    │   ├── [01;32mAMS_DATA-92.cif[0m
    │   └── [01;32mSTRU[0m
    ├── [34;42mPOSCAR2STRU[0m
    │   └── [01;32mPOSCAR[0m
    ├── [34;42mSTRU2CIF[0m
    │   └── [01;32mSTRU[0m
    └── [34;42mSTRU2POSCAR[0m
        └── [01;32mSTRU[0m
    
    4 directories, 5 files
    

### 4.2 STRU 转 CIF


```python
from ase.io import read, write
import os

cs_dir = os.path.join(parent_path, "STRU2CIF")
cs_stru = os.path.join(cs_dir, 'STRU')
cs_atoms= read( cs_stru, format='abacus')
cs_vasp = os.path.join(cs_dir, 'STRU.cif')
write(cs_vasp, cs_atoms, format='cif')
```

### 4.3 POSCAR 转 STRU

POSCAR 是常见第一性原理计算软件 VASP 的结构文件，这里提供了一种 POSCAR 转 STRU 文件的方法：


```python
from ase.io import read, write
import os

cs_dir = os.path.join(parent_path, "POSCAR2STRU")
cs_vasp = os.path.join(cs_dir, "POSCAR")
cs_atoms = read(cs_vasp, format="vasp")
cs_stru = os.path.join(cs_dir, "STRU")

# # 需要设置原子磁矩时，可添加代码：
# cs_atoms.set_initial_magnetic_moments([1.0,1.0,1.0,2.0])  # 后面列表里列出每个原子的磁矩；

# # 如果需要设置每个原子每个方向的磁矩，则需要设置为二维数组：
# cs_atoms.set_initial_magnetic_moments([[1.0,1.0,1.0],
#                                        [1.0,1.0,2.0],
#                                        [1.0,1.0,3.0],
#                                        [1.0,1.0,4.0]])

pp = {"Al": "Al.PD04.PBE.UPF"}
basis = {"Al": "Al_gga_10au_100Ry_3s3p2d.orb"}
write(cs_stru, cs_atoms, format="abacus", pp=pp, basis=basis)

```

### 4.4 STRU 转 POSCAR


```python
from ase.io import read, write
from pathlib import Path
from ase.calculators.abacus import Abacus

cs_dir = os.path.join(parent_path, "STRU2POSCAR")
cs_stru = os.path.join(cs_dir, 'STRU')
cs_atoms= read(cs_stru, format='abacus')
cs_vasp = os.path.join(cs_dir, 'POSCAR')
write(cs_vasp, cs_atoms, format='vasp')
```
