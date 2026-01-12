# ABACUS 使用教程｜如何转换 STRU 文件


```python
'''@ 作者: 赵天琦'''
```

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
