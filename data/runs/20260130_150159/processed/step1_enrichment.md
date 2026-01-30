# 使用abacustest准备ABACUS输入文件

从结构文件出发开始做计算是一个常见的需求。如果希望使用ABACUS进行想要的计算，通常需要将常见格式的结构文件（如CIF、VASP POSCAR等）转化为ABACUS使用的STRU文件，按需配置赝势文件和数值原子轨道文件（使用LCAO基组时），并修改INPUT文件，然后才可以开始计算。这是一个稍显繁琐的过程，并且有较多结构需要做计算时不够方便。

abacustest提供了一个使用结构文件作为输入，根据要求快速准备完整的输入文件夹的方式，可以自动配置赝势、轨道和其它输入文件（INPUT、KPT等）的工作流，均使用下面的命令作为开头：

```Python
abacustest model inputs
```

下面将用案例说明如何使用。

## 下载推荐的APNS-pp-orb-v1赝势轨道库

如果还没有下载当前对一般计算推荐使用的APNS-pp-orb-v1赝势轨道库，可使用下面的命令下载到当前文件夹：

```Python
abacustest model inputs --download-pporb
```

命令执行完成后，将在当前目录下出现两个文件夹：

```Python
$ ls
apns-orbitals-efficiency-v1  apns-pseudopotentials-v1
```

这两个文件夹是对原始的APNS-pp-orb-v1版本的重新打包，分别包含从H到Bi的83种元素推荐使用的赝势和轨道。所有轨道均为DZP水平，并且对镧系元素，4f电子被视为核电子，且选用了截断半径为8 au的轨道。后面的输入文件准备需要准备赝势和轨道，可以在环境变量中设置ABACUS\_PP\_PATH和ABACUS\_ORB\_PATH，值为下载的赝势和轨道的路径：

```Python
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

也可以在准备ABACUS输入文件夹时通过--pp和--orb选项指定赝势和轨道的路径。

## 准备MgO做SCF计算的输入文件

如果你已经有一个MgO的CIF文件MgO.cif，你可以使用下面的命令准备ABACUS对MgO做SCF计算的输入文件夹：

```Python
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO
```

该命令的各项参数含义为：

* -f: 用于转化为STRU文件的结构
* --ftype：结构文件格式
* --lcao：使用LCAO基组。如不使用该选项，将使用PW基组
* --folder-syntax：指定生成的ABACUS输入文件夹的名称。如果不设置，将从000000开始计数。

该命令将从ABACUS\_PP\_PATH和ABACUS\_ORB\_PATH中选取该结构所需的赝势和轨道，并自动配置适用于SCF计算的INPUT文件。执行完毕后的文件目录结构为：

```Bash
MgO
├── INPUT
├── Mg_gga_10au_100Ry_2s1p.orb -> /home/abc/apns-orbitals-efficiency-v1/Mg_gga_10au_100Ry_2s1p.orb
├── Mg.PD04.PBE.UPF -> /home/abc/apns-pseudopotentials-v1/Mg.PD04.PBE.UPF
├── O_gga_6au_100Ry_2s2p1d.orb -> /home/abc/apns-orbitals-efficiency-v1/O_gga_6au_100Ry_2s2p1d.orb
├── O.upf -> /home/abc/apns-pseudopotentials-v1/O.upf
├── STRU
└── struinfo.txt
MgO.cif
run.sh
setting.json
struinfo.json
```

其中：

* run.sh和setting.json是使用dpdispatcher提交计算的文件，如有需要可自行编辑内容，不需要可直接删除
* struinfo.txt和struinfo.json中记录了生成STRU的原始结构路径，如不需要也可删除
* UPF格式的赝势文件和orb格式的轨道文件为软链接到原始目录中的文件，如果不需要软链接，可使用--copy-pp-orb选项复制文件

INPUT文件的内容如下，正确使用了LCAO基组，并默认提供了一套适合很多体系的基础参数，可自行修改：

```Python
calculation     scf
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
#cal_force     1
#cal_stress     1
kspacing     0.14 # unit in 1/bohr
#gamma_only     0
```

## 准备Fe2O3做SCF计算的输入文件

Fe2O3是磁性材料，Fe的磁矩约为4 μB，在做第一性原理计算时，需要开启共线自旋极化，设置初始磁矩，并按需设置DFT+U参数。使用下面的命令可满足该要求，并设置对Fe的d轨道使用3 eV的Ueff：

```Python
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

新出现的参数的含义为：

* --nspin：设置自旋极化，与INPUT文件中的nspin含义相同。设置为2表示使用共线自旋极化
* --init\_mag Fe 4.0：对所有Fe原子设置4 μB的初猜磁矩
* --dftu：使用DFT+U
* --dftu\_param：对不同元素设置Ueff。这里设置的Fe 3.0将默认为Fe使用3 eV的Ueff，施加的轨道为对Fe常使用DFT+U的d轨道。

--init\_mag和--dftu\_param都可以设置多个元素的初猜磁矩和DFT+U参数，依次列出即可，下面是一个例子：

```Python
abacustest model inputs -f Co2FeAl.cif --ftype cif --lcao --nspin 2 --init_mag Co 1.5 Fe 2.0 --dftu --dftu_param Co 1.0 Fe 3.0
```

执行完毕后，INPUT文件的内容为如下，启用了nspin=2，相应适当调低了mixing\_beta到0.4，并设置onsite\_radius和out\_mul 1，在ABACUS运行结束后会在OUT.ABACUS/mulliken.txt中记录原子磁矩以及各条轨道对磁矩的贡献。此外还设置了DFT+U的相关参数，并且启用了uramping，促进使用DFT+U时的SCF收敛。

```Python
INPUT_PARAMETERS
calculation     scf
symmetry     0
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.4
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
#cal_force     1
#cal_stress     1
kspacing     0.14 # unit in 1/bohr
#gamma_only     0
nspin     2
onsite_radius     3
out_mul     1
dft_plus_u     1
orbital_corr     2 -1
hubbard_u     3.0 0
uramping     3.0
mixing_restart     0.001
```

STRU文件的内容为：

```Python
ATOMIC_SPECIES
Fe 55.845000 Fe_ONCV_PBE-1.2.upf
O 15.999400 O.upf

NUMERICAL_ORBITAL
Fe_gga_7au_100Ry_4s2p2d1f.orb
O_gga_6au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
    5.09190205000     0.00000000000     0.00000000000
   -2.54595102500     4.40971652888     0.00000000000
    0.00000000000     0.00000000000    13.77429765000

ATOMIC_POSITIONS
Cartesian

Fe
0.000000
12
    0.00000000000     0.00000000000     1.99971355156 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000    11.77458409843 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000     4.88743527343 1 1 1 mag   4.00000000
    0.00000000000     0.00000000000     8.88686237657 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     6.59114610157 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     2.59171899844 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963     9.47886782343 1 1 1 mag   4.00000000
    2.54595102500     1.46990550963    13.47829492657 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925    11.18257865157 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     7.18315154844 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     0.29600272344 1 1 1 mag   4.00000000
    0.00000000000     2.93981101925     4.29542982656 1 1 1 mag   4.00000000

O
0.000000
18
    4.31436631561     1.34673139667    10.33072323750 1 1 1 mag   0.00000000
   -1.76841529061     3.06298513222     3.44357441250 1 1 1 mag   0.00000000
    1.76841529061     3.06298513222    10.33072323750 1 1 1 mag   0.00000000
    0.77753573439     1.34673139667     3.44357441250 1 1 1 mag   0.00000000
    1.55507146878     0.00000000000    10.33072323750 1 1 1 mag   0.00000000
    0.99087955622     4.40971652888     3.44357441250 1 1 1 mag   0.00000000
    1.76841529061     2.81663690629     1.14785813750 1 1 1 mag   0.00000000
    3.32348675939     0.12317411296     8.03500696250 1 1 1 mag   0.00000000
    1.76841529061     0.12317411296     1.14785813750 1 1 1 mag   0.00000000
    3.32348675939     2.81663690629     8.03500696250 1 1 1 mag   0.00000000
    4.10102249378     1.46990550963     1.14785813750 1 1 1 mag   0.00000000
    0.99087955622     1.46990550963     8.03500696250 1 1 1 mag   0.00000000
   -0.77753573439     4.28654241592     5.73929068750 1 1 1 mag   0.00000000
    0.77753573439     1.59307962259    12.62643951250 1 1 1 mag   0.00000000
   -0.77753573439     1.59307962259     5.73929068750 1 1 1 mag   0.00000000
    0.77753573439     4.28654241592    12.62643951250 1 1 1 mag   0.00000000
    1.55507146878     2.93981101925     5.73929068750 1 1 1 mag   0.00000000
   -1.55507146878     2.93981101925    12.62643951250 1 1 1 mag   0.00000000
```

可见为所有Fe原子设置了4 μB的初始磁矩。

## 选择不同的任务类型

可以为常见的任务类型（scf，relax，cell-relax，md等）设置默认INPUT参数，对很多情况都是合适的。例如需要优化MgO的晶胞结构，可使用的命令为

```Python
abacustest model inputs -f MgO.cif --ftype cif --lcao --jtype cell-relax --folder-syntax MgO-cellrelax
```

生成的INPUT文件如下，将calculation设置为cell-relax，并且设置了在ABACUS中适合做变胞优化的cg优化算法以及合理的力收敛限与应力收敛限。如果设置为relax，会将calculation设置为relax，并去除cal\_stress和stress\_thr参数。

```Python
INPUT_PARAMETERS
calculation     cell-relax
symmetry     1
ecutwfc     100
scf_thr     1e-07
scf_nmax     100
smearing_method     gauss
smearing_sigma     0.015
mixing_type     broyden
mixing_beta     0.8
basis_type     lcao
ks_solver     genelpa
precision     double  # or single
cal_force     1
cal_stress     1
kspacing     0.14 # unit in 1/bohr
relax_method     cg # or bfgs, bfgs_trad, cg_bfgs, sd, fire
relax_nmax     60
force_thr_ev     0.01  # unit in eV/A
stress_thr     0.5 # unit in kbar
fixed_axes     None # or volume, shape, a, b, c, ab, ac, bc; only valid for cell-relax calculation to fix some axes
#gamma_only
```

## 自定义INPUT和KPT

前面的设置中都使用了默认的INPUT参数和k点设置（通过INPUT中的kspacing参数）。如果需要自定义INPUT参数和k点，可通过--input和--kpt选项指定。例如如果希望在前面的Fe2O3的例子中将smearing\_sigma设置为0.001，将mixing\_beta降低到0.2，并使用5×5×5的k点，可以先准备下面的一个INPUT\_template的文件：

```Python
INPUT_PARAMETERS
smearing_sigma     0.001
mixing_beta     0.2
```

然后使用下面的命令：

```Python
abacustest model inputs -f Fe2O3.cif --ftype cif --lcao --nspin 2 --input INPUT_template --kpt 5 5 5 --init_mag Fe 4.0 --dftu --dftu_param Fe 3.0
```

执行完毕后，将使用INPUT\_template中的参数替换默认参数，并生成以Gamma点为中心，采样密度为5×5×5的KPT文件，同时去掉INPUT文件中的kspacing参数。

## 批量使用结构准备ABACUS输入文件

部分情况下需要为相同的一批结构配置相同的输入文件。例如，在做Pd(100)面的表面能随层数的收敛性测试时，需要为一批不同厚度的Pd(100)面的结构文件做相同的优化计算。假设结构文件的名称如下：

```Python
Pd100_1layer.vasp  Pd100_3layer.vasp  Pd100_5layer.vasp  Pd100_7layer.vasp
Pd100_2layer.vasp  Pd100_4layer.vasp  Pd100_6layer.vasp  Pd100_8layer.vasp
```

可使用下面的命令准备ABACUS对这些结构做优化的任务文件夹：

```Python
abacustest model inputs -f Pd100_*layer.vasp --ftype poscar --lcao --jtype relax --folder-syntax "x[:-5]"
```

这里的--folder-syntax命令可使用Python的字符串切片语法，为每个结构根据结构文件名设置目录名称。这里使用的"x[:-5]"的含义为生成的目录名称为文件名除去后5个字符，在这个例子中的含义为去除".vasp"。因此会为每个结构生成名称为Pd100\_\*layer的ABACUS任务文件夹。
