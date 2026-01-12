# 写在前面

本文记录了笔者使用ABACUS对92号元素铀的平面波基组计算、收敛性判断、晶格弛豫过程的经验与遇到的问题，希望该笔记能对读者遇到的问题起到帮助。

# **铀的基本信息**

铀的元素符号为U，原子序数92，原子量为238.03，是自然界中能够找到的最重原生元素。自然界中存在三种铀的同位素，铀-234、铀-235和铀-238，都具有放射性。铀存在三种晶格结构： α-U为斜方结构、β-U为正方结构、γ-U为体心立方结构。

铀的单质为银白色的致密金属，铀的新切面呈发亮的钢灰色，但在室温空气中逐渐生成黑色氧化膜。铀密度高，相对密度约18.95，熔点1135℃，沸点4134℃。

铀有+3、+4、+5、+6四种价态，以+4和+6为主要价态，铀是正电性很强的活泼元素，与几乎所有非金属元素（惰性气体除外）反应生成化合物，常以$U^{3+}、U^{4+}、UO_2^+、UO_2^{2+}$离子存在，此外$UF_6$是气态铀同位素分离的原料。

铀-234的半衰期为$2.47*10^5$年，铀-235的半衰期为$7.00*10^8$年，铀-238的半衰期为$4.51*10^9$年。1942年前铀主要用作玻璃和陶瓷的着色剂，但是随着链式核裂变反应的被发现，铀开始被用作核燃料，最先被制作成原子弹、氢弹，也可用于制作铀核反应堆。

# 密度泛函理论的基本原理
密度泛函理论(Density functional theory ，DFT)是一种研究多电子体系电子结构的方法，是目前凝聚态物理计算材料学和计算化学领域最常用的方法之一。

## 主要目的
研究原子结构/多原子结构也就是要解对应的定态薛定谔方程：
$$
\hat{H}\psi(\vec{r})=E\psi(\vec{r})
$$
但是对于一个存在着多个电子、原子核，甚至彼此之间还有相互作用的体系，其哈密顿算符非常复杂乃至于无法直接求解，因此需要通过近似方法来进行求解。





## Hohenberg-Kohn定理

Hohenberg-Kohn定理是密度泛函理论的基础，它由两部分组成

定理一：体系的基态能量仅仅是电子函数的泛函，即$E_0$由基态的电子概率密度$ρ_0(x,y,z)$唯一确定。

定理二：以基态密度为变量，体系的能量可以通过采用变分法取得最小值，且这一最小值就是基态能量。

Hohenberg-Kohn定理证明了存在着一个联系体系能量与其电子密度分布的普适性密度泛函，但是没有给出这个泛函是什么，事实上这一泛函的具体形式至今未知。

## Born-Oppenheimer近似
Born-Oppenheimer近似是一种非常重要的近似方法，由于原子核的质量比电子大得多，运动比电子慢得多，因此在一般情况下，可以把原子核看做不动而把原子核和电子的相对运动看做是电子围绕不动的原子核运动的问题。

对于M个原子核和N个电子组成的体系，Born-Oppenheimer近似将哈密顿量中的变量由由 $3(M+N)$ 降低为 $3N$。但是这对于计算一个复杂体系来说还远远不够。

## Kohn-Sham方程
Kohn-Sham方程在进一步近似的前提下给出了Hohenberg-Kohn定理所预言的泛函的表达式，它标志着密度泛函理论走向实用。

Kohn-Sham方程基于Kohn-Sham假设：有相互作用系统的基态电荷密度，可以被一个无相互作用系统的基态电荷密度表示出来。这一假设使得多电子系统被转化为了无相互作用的单电子系统，从而能够较为简便地计算。但是应当注意的是，Kohn-Sham假设目前仍然停留在假设的状态，未被证实也未被证伪。

得到的Kohn-Sham方程形式如下：
$$
\left\{-\frac{\nabla^2}{2}-\sum_q \frac{Z_q}{\left|r-R_q\right|}+\int \frac{\rho(r)}{\left|r-r^{\prime}\right|} d r^{\prime}+V_{X C}(r)\right\} \phi_i(r)=\varepsilon_i \phi_i(r)
$$

## 求解Kohn-Sham方程-电子自洽迭代
电子自洽迭代(Self-Consistent Field, SCF)方法的基本思想为：\
1.随机给出一个基态电荷密度\
2.由基态电荷密度计算有效势能，代入求解Kohn-Sham方程，计算出新的基态电荷密度\
3.判断新的基态电荷密度和原来的基态电荷密度相差多少，如果相差不多则输出，否则继续用新得到的基态电荷密度重复2.3.两个步骤进行迭代

解出Kohn-Sham方程也就意味着得到了真实体系对应的虚拟体系的电子密度分布与能量分布，能够一定程度上反映真实体系的状态。



# ABACUS软件简介
原子算筹(ABACUS)是国内自主研发的开源第一性原理材料计算软件，支持多种基于密度泛函理论（Density Functional Theory，简称DFT）的算法以及平面波与数值原子轨道两种基组。

目前ABACUS的主要功能包括电子结构自洽迭代计算、原子结构优化以及分子动力学计算等。采用模守恒赝势和周期性边界条件, 可以对晶格对称性、布里渊区的 k 点对称性、电荷密度对称性以及力的对称性进行分析。

# ABACUS使用方法入门
## 需要准备的文件

使用ABACUS进行计算最重要的是需要准备INPUT、KPT、STRU三个文件(无后缀名，Windows端可以选择使用记事本打开并编辑)。

此外，在计算过程中可能会用到原子的赝势和数值原子轨道文件(lcao，原子轨道线性组合计算需要用到轨道文件；pw，平面波基组计算不需要用到轨道文件)，赝势和数值原子轨道都可以在https://github.com/kirk0830/ABACUS-Pseudopot-Nao-Square/tree/main/download 中下载。赝势生成时,是解了个一维的kohn-sham方程,解方程时是需要指定交换关联泛函的,所以赝势文件天然会带着交换关联泛函，目前最常用的是PBE交换关联泛函。

接下来我们以U的自洽迭代计算为例介绍INPUT、KPT、STRU三个文件。

`STRU`：结构文件，包含了原子种类、原子位置、晶格常数以及晶格向量等信息
```
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf   
#元素名称、相对原子质量、原子赝势的文件名，这里使用平面波基组计算，因此只需要准备赝势而不需要准备原子轨道

# 如果准备了原子轨道按下方的格式输入
# NUMERICAL_ORBITAL
# <原子轨道的文件名>  

LATTICE_CONSTANT    #晶格常数的单位，默认单位为Bohr
1.8897259886 	# 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS   #晶格三个边对应的三个向量，这里采用了α-U的晶格结构，正交晶系，三个晶格常数如下
2.8541 0.00000 0.00000
0.00000 5.8692 0.00000
0.00000 0.00000 4.9563

ATOMIC_POSITIONS 
Direct #原子位置的输入方式，直接输入
U #元素名称
0.0 
1 #选取的晶胞中的原子数量
0.0 0.0 0.0 1 1 1 #前三位为原子位置的对应向量，后三位代表着原子运动的自由度

# 注：在笔者这次计算中选取的原子位置有误，正确的应当是
# ATOMIC_POSITIONS 
# Direct
# U
# 0.0 
# 2
# 0.0 0.105 0.25 1 1 1
# 0.0 0.895 0.75 1 1 1

```

更多的STRU文件的可选信息参考官网https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html

`KPT`：包含了布里渊区积分所需的k点信息

```
K_POINTS
0
Gamma
4 4 4 0 0 0   #4 4 4表示在三个方向上各取4个k点
```

`INPUT`：包含了计算过程中所需的各种参数，定义和控制计算任务
```
INPUT_PARAMETERS
suffix 			U    #物质名称
ntype 			1    #计算的物质中有几种元素
pseudo_dir 		../../PP_ORB  #赝势的存放位置
orbital_dir 		../../PP_ORB  #数值原子轨道的存放位置
ecutwfc 		100  #波函数截断能，用于确定平面波基组的大小。这里的值为100 Rydberg，表示波函数截断能为100 Ry。
scf_thr 		1e-6  #自洽场（SCF）计算的能量收敛阈值，这里设为1e-4 Rydberg，表示能量收敛到1e-4 Ry时计算停止。
basis_type 		pw   #表示基组类型，这里设置为"pw"，表示使用平面波基组。
calculation 		scf    #计算类型，这里设置为"scf"，表示进行自洽性迭代计算
```
更多的INPUT文件的可选参数参考官网https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html


本次准备的文件夹存放结构为
```
ABACUS
    ├── U_CELLRELAX
    │   └── SCF
    │       ├── INPUT
    │       ├── KPT
    │       └── STRU
    ├── U_PW
    │   └── SCF
    │       ├── INPUT
    │       ├── KPT
    │       └── STRU
    ├── UO2_CELLRELAX
    │    └── SCF
    │        ├── INPUT
    │        ├── KPT
    │        └── STRU
    └── PP_ORB
         ├── O_gga_7au_100Ry_2s2p1d.orb
         ├── O_ONCV_PBE-1.0.upf
         └── U-5spdf.upf
```

## 自洽迭代计算
进入INPUT、KPT、STRU三个文件所在的目录，执行以下命令即可开始计算。运行成功后，输出文件将在同目录下OUT.N文件夹中，其中的running_scf.log文件可以查找到晶体的总能量。


```
# 进入工作文件夹
cd ./ABACUS/U_PW/SCF
# OMP_NUM_THREADS=1 表示使用单线程，如果你的机器配置比较高，可以使用多线程，比如 4 线程，就可以写成 OMP_NUM_THREADS=4
# mpirun -n 后面的数字表示计算所使用的 CPU 核心数，这里使用 2 个核心，你可以根据你的机器配置进行修改。
OMP_NUM_THREADS=1 mpirun -n 2 abacus
```

## 收敛性计算

分别改变INPUT中ecut和KPT中k的取值，观察计算结果的能量与变量的关系，结果如下。调节ecut时k为6，调节k时ecut为100。
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/31613/710b32742df043ecadf7ff36110ef426/a3cce6a3-350d-4a59-8f21-357d64a55196.jpeg)
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/31613/710b32742df043ecadf7ff36110ef426/3699ec4d-aa05-4e81-b562-4b972d033080.jpeg)
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/31613/710b32742df043ecadf7ff36110ef426/278454fb-f670-4931-b906-e9973f99d239.jpeg)
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/31613/710b32742df043ecadf7ff36110ef426/bedc9003-e432-4728-9191-17c2636333f2.jpeg)
可以看到随着ecut和k的增大，运行时间整体上都是上升趋势，但是在ecut为90、k为6时基本上可以视为已经收敛。

# 晶格弛豫
## U的晶格弛豫

### 准备的文件
（晶格常数都加了0.1埃以观察优化效果，但是似乎这一改变量过大导致计算花费时间很多。此外，晶格中原子位置仍存在错误，本例子仅供运行测试参考，正确的设定见上文）

`INPUT`:
```
INPUT_PARAMETERS
suffix 			U
ntype 			1
pseudo_dir 		../../PP_ORB
orbital_dir 		../../PP_ORB
ecutwfc 		100 
scf_thr 		1e-4
basis_type 		pw
calculation 		cell-relax
force_thr_ev		0.1	
stress_thr		20
relax_nmax		100		
out_stru		1
cell_factor     2
relax_scale_force 0.4
```
`KPT`:
```
K_POINTS
0
Gamma
4 4 4 0 0 0
```
`STRU`:
```
ATOMIC_SPECIES
U 238.0508 U-5spdf.upf

LATTICE_CONSTANT
1.8897259886 	# 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
2.9541 0.00000 0.00000
0.00000 5.9692 0.00000
0.00000 0.00000 5.0563

ATOMIC_POSITIONS
Direct 
U
0.0 
4
0.0 0.0 0.0 0 0 0
0.5 0.0 0.0 0 0 0
0.0 0.5 0.0 0 0 0
0.0 0.0 0.5 0 0 0
```


同样是进入STRU等三个文件所在的目录，利用abacus命令运行，使用的运行参数如下


```
OMP_NUM_THREADS=2 mpirun -n 8 abacus
#8核2线程
```

由于时间过长，在运行4285s后提前终止程序，此时应力为-233.633023 KBAR，得到的结果为

Lattice vectors: (Cartesian coordinate: in unit of a_0)

            +3.816061           +0.000000           -0.000000
            +0.000000           +5.871062           +0.000000
            -0.000000           +0.000000           +5.025222
这一错误的初始结构（面心正交晶系）下结构优化的趋势是向更大、更接近立方晶系的方向靠近。这可能是由于**晶格弛豫算法会寻找离给定的初始构型最接近的一个势能极小值、而非势能的最小值**。因此如果初始构型离实际构型太远，那么优化得到的构型可能不是真实构型，而是另一个势能极小值位置的构型。

## 错误处理

在U的晶格弛豫过程中，起初没有修改cell_factor和relax_scale_force参数，运行时报错
```
Not enough space allocated for radial FFT.
 It is due to the rapid change of the size of cell:
 Try reseting a larger cell_factor parameter in INPUT
 Or try reseting a smaller relax_scale_force parameter in INPUT
```
查阅官方文档https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html ，得知cell-factor项应当超过优化过程中的最大线性收缩，relax_scale_force控制共轭梯度下降的步长、对于较大的系统一个较小的步长有助于规避崩溃。于是将cell_factor从默认的1.2改为了2，将relax_scale_force从默认的0.5调整至0.3。

如果遇到计算不收敛的问题，可以参考https://mcresearch.github.io/abacus-user-guide/abacus-conv.html


## UO$_2$的结构优化

UO$_2$为面心立方结构，铀组成面心立方，氧占据正四面体空隙，晶格常数a=5.471埃。准备的文件如下

`INPUT`:
```
INPUT_PARAMETERS
suffix 			UO2
ntype 			2
pseudo_dir 		../../PP_ORB
orbital_dir 		../../PP_ORB
ecutwfc 		100 
scf_thr 		1e-4
basis_type 		pw
calculation 		cell-relax
force_thr_ev		0.1	
stress_thr		20	
relax_nmax		10		
out_stru		1
cell_factor     2
relax_scale_force 0.4
```
`KPT`:
```
K_POINTS
0
Gamma
4 4 4 0 0 0
```
`STRU`:
```ATOMIC_SPECIES
U 238.0508 U-5spdf.upf
O  15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
O_gga_7au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897259886 	# 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
5.571 0.00000 0.00000
0.00000 5.571 0.00000
0.00000 0.00000 5.571

ATOMIC_POSITIONS
Direct 
U
0.0 
4
0.0 0.0 0.0 0 0 0
0.5 0.0 0.0 0 0 0
0.0 0.5 0.0 0 0 0
0.0 0.0 0.5 0 0 0

O
0.0
8
0.25 0.25 0.05 0 0 0
0.25 0.25 0.75 0 0 0
0.25 0.75 0.25 0 0 0
0.25 0.75 0.75 0 0 0
0.75 0.25 0.25 0 0 0
0.75 0.25 0.75 0 0 0
0.75 0.75 0.25 0 0 0
0.75 0.75 0.75 0 0 0
```

运行参数同样为


```
OMP_NUM_THREADS=2 mpirun -n 8 abacus
#8核2线程
```

设置初始参数时给晶格参数加0.1埃用于观察优化效果、力收敛阈值设置为0.1 eV/Angstrom，应力收敛阈值设置为200kPa。这样设置的条件过于宽松，第一次迭代后的应力就只有52kPa，直接达到了收敛。

将应力收敛阈值调整为20kPa。

最终结果为LATTICE VECTORS: (CARTESIAN COORDINATE: IN UNIT OF A0)

            +5.595584           -0.107596           -0.027923
            -0.107596           +5.595584           -0.027923
            -0.027984           -0.027984           +5.643928
可以看到晶格的三个边不再两两正交，但是接近正交，晶格的边长不减反增，但是夹角变大，
也是一种舒缓应力的方式，但是和预期有些不符。

