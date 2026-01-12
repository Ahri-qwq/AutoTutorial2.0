# ABACUS 使用教程｜电子自洽迭代（LCAO 基组与 PW 基组）

<a href="https://nb.bohrium.dp.tech/detail/7417640496" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

## 1. ABACUS 自洽场计算（LCAO）

ABACUS 是一款基于密度泛函理论（DFT）的第一性原理计算软件，自洽场计算（SCF, self-consistent field calculation）是其计算的核心过程。在某些教程中，这个过程也被称为「电子步」，即通过优化电荷密度来寻找势能最低点的过程。

在 DFT 计算中，首要任务是分析体系的势能面。通过 SCF 计算，可以得到一个体系的基态结构和基态能量。其具体的流程如下：

```mermaid
graph TD;
  A[随机给一个初始电荷密度] --> B[建立有效势能 Veff];
  B[建立有效势能 Veff] --> C[求解 Kohn-Sham (KS)  方程];
  C[求解 Kohn-Sham (KS)  方程] --> D[计算出一个新的电荷密度];
  D[计算出一个新的电荷密度] --> E[判断是否达到自洽];
  A[随机给一个初始电荷密度] --> B[建立有效势能 Veff];
  E[判断是否达到自洽] -- 成立 --> F[输出各种性质物理量];
  E[判断是否达到自洽] -- 不成立 --> B[建立有效势能 Veff];
```

有关更多的 DFT 理论背景介绍，推荐阅读：[了解第一性原理计算与密度泛函理论](https://nb.bohrium.dp.tech/detail/7317630945)

ABACUS 以支持 LCAO（Linear Combination of Atomic Orbital，原子轨道线性组合）基组在计算周期性凝聚态系统方面而闻名，因此从 LCAO 自洽场（SCF）计算示例开始是一个不错的选择。这里，选择了 FCC MgO 作为快速开始示例。

### 1.1 输入文件

在正式开始 ABACUS SCF 计算之前，请确保这些文件已经准备好并存储在工作目录中：

- `INPUT`：包含了计算过程中所需的各种参数，定义和控制计算任务；
- `STRU`：结构文件，包含了原子种类、原子位置、晶格常数以及晶格向量等信息；
- `KPT`：包含了布里渊区积分所需的k点信息；
- `*.upf`：包含了原子的赝势信息；
- `*.orb`：包含了原子轨道的数值表示；

在示例数据集中，我们为你准备好了格式转换的示例文件。我们可以直接访问：（你可以在左侧点击数据集查看相应文件）：




```python
! tree /bohr/
```

出于安全考虑，我们没有数据集所在文件夹的写入权限，因此我们将其复制到 `/data/` 目录下:


```python
! cp -nr /bohr/ /data/
```

我们在这里定义一些路径，并切换到工作路径，方便后续调用：


```python
import os

bohr_dataset_url = "/bohr/abacus01-9qud/v4/"  # url 可从左侧数据集复制
work_path = os.path.join("/data", bohr_dataset_url[1:])
os.chdir(work_path)
print(f"当前路径为：{os.getcwd()}")
```

    当前路径为：/data/bohr/abacus01-9qud/v4
    



#### 1.1.1 STRU 文件

##### 文件示例

STRU 文件包含晶格几何信息，赝势和数值轨道文件的名称和/或位置，以及关于系统的结构信息。我们可以通过 cat 命令进行查看：



```python
! cat ./ABACUS_SCF/MgO_LCAO/SCF/STRU
```

    ATOMIC_SPECIES
    Mg 24.305 Mg_ONCV_PBE-1.0.upf
    O  15.999 O_ONCV_PBE-1.0.upf
    
    NUMERICAL_ORBITAL
    Mg_gga_8au_100Ry_4s2p1d.orb
    O_gga_7au_100Ry_2s2p1d.orb
    
    LATTICE_CONSTANT
    1.8897259886 	# 1.8897259886 Bohr = 1.0 Angstrom
    
    LATTICE_VECTORS
    4.27957 0.00000 0.00000
    0.00000 4.27957 0.00000
    0.00000 0.00000 4.27957
    
    ATOMIC_POSITIONS
    Direct 
    Mg 
    0.0 
    4 
    0.0 0.0 0.0 0 0 0 
    0.0 0.5 0.5 0 0 0
    0.5 0.0 0.5 0 0 0
    0.5 0.5 0.0 0 0 0
    
    O
    0.0 
    4 
    0.5 0.0 0.0 0 0 0
    0.5 0.5 0.5 0 0 0 
    0.0 0.0 0.5 0 0 0 
    0.0 0.5 0.0 0 0 0
    




下面显示了 LCAO 计算中的 FCC MgO 的 STRU 文件：

```sh
#This is the atom file containing all the information
#about the lattice structure.

ATOMIC_SPECIES
Mg 24.305  Mg_ONCV_PBE-1.0.upf  # element name, atomic mass, pseudopotential file
O  15.999 O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Mg_gga_8au_100Ry_4s2p1d.orb
O_gga_8au_100Ry_2s2p1d.orb

LATTICE_CONSTANT
1.8897259886 		# 1.8897259886 Bohr =  1.0 Angstrom

LATTICE_VECTORS
4.25648 0.00000 0.00000  
0.00000 4.25648 0.00000
0.00000 0.00000 4.25648

ATOMIC_POSITIONS
Direct                  #Cartesian(Unit is LATTICE_CONSTANT)
Mg                      #Name of element        
0.0                     #Magnetic for this element.(Be careful: value 1.0 refers to 1.0 bohr mag, but not fully spin up !!!)
4                       #Number of atoms
0.0  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
O                       #Name of element        
0.0                     #Magnetic for this element.
4                       #Number of atoms
0.5  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
```

##### 参数解释

- ATOMIC_SPECIES

从左到右依次为原子种类、相对原子质量、赝势文件

- NUMERICAL_ORBITAL

原子轨道文件

* 赝势文件 `Mg_ONCV_PBE-1.0.upf` 和 `O_ONCV_PBE-1.0.upf` 应放在 `pseudo_dir` 目录下，
* 轨道文件 `Mg_gga_8au_100Ry_4s2p1d.orb` 和 `O_gga_8au_100Ry_2s2p1d.orb` 应放在 `orbital_dir` 目录下。

`pseudo_dir` 和 `orbital_dir` 文件路径在 `INPUT` 文件中进行设置，赝势和轨道文件可以从 [ABACUS 网站](http://abacus.ustc.edu.cn/pseudo/list.htm) 下载。

- LATTICE_CONSTANT

晶格常量，晶胞长度的缩放系数。

- LATTICE_VECTORS

晶格向量，即晶胞的 x,y,z 轴长度（实际长度需乘以晶格常量）

- ATOMIC_POSITIONS

原子位置信息

有关更多 STRU 文件的参数信息，请参阅：[The STRU file](https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html)

 我们可以使用 Python ASE 库来可视化查看一下 STRU 文件对应的晶体结构：

如果你的电脑还没有安装 ase，我们可以运行以下命令安装 ase-abacus：

```bash
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
python3 setup.py install
```

ASE 模块还有许多妙用，例如快速生成 STRU 文件，请查看：[2.1.5 使用 ase 快速生成 STRU 文件]()

#### 1.1.2 INPUT 文件

接下来，需要 INPUT 文件，该文件设置了所有关键参数，以指导 ABACUS 如何进行计算以及输出什么内容：

##### 示例文件

```sh
# INPUT_PARAMETERS
suffix 			MgO
ntype 			2
pseudo_dir 		./PP_ORB
orbital_dir 	        ./PP_ORB
ecutwfc 		100     # Rydberg
scf_thr 		1e-6    # Rydberg
basis_type 		lcao
calculation 	        scf 
```

##### 参数解释

* suffix：系统名称，默认为ABACUS
* ntype：单位晶胞中元素的种类数目
* pseudo_dir：提供赝势文件的目录
* orbital_dir：提供轨道文件的目录
* ecutwfc：波函数展开的平面波能量截止（单位：Rydberg）
* scf_thr：电荷密度收敛阈值（单位：Rydberg）
* basis_type：用于展开电子波函数的基函数类型
* calculation：ABACUS要执行的计算类型

参数列表始终以关键字 `INPUT_PARAMETERS` 开头。在 `INPUT_PARAMETERS` 之前的任何内容都将被忽略。以 `#` 或 `/` 开头的任何行也将被忽略。

**注意：INPUT 参数设置的第一原则是：简洁明了。** 参数不是越多越好，对于大多数场景，默认值已不失为一个很好的选择。请确保你设置的参数你都明确清楚其意义。

**再注意：用户无法将文件名“INPUT”更改为其他名称。** 布尔参数（如 out_chg）可以通过使用 True 和 False、1 和 0 或 T 和 F 来设置。大小写不敏感，因此也支持使用 true 和 false、TRUE 和 FALSE 以及 t 和 f 设置布尔值的其他偏好。

有关 INPUT 文件的所有可能参数信息，请参阅：[Full List of INPUT Keywords](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html)

#### 1.1.3 KPT 文件

最后一个必需的输入文件称为 KPT，此文件包含关于布里渊区采样的k点网格设置的信息。以下是一个示例：

```
K_POINTS
0           # k点的总数，`0'表示自动生成
Gamma       # Monkhorst-Pack方法的类型，`Gamma'或`MP'
4 4 4 0 0 0 # 前三个数字：沿着倒数向量的细分
            # 后三个数字：网格的偏移
```

#### 1.1.4 赝势和原子轨道文件

##### 官网信息

**以下是 ABACUS 官网的赝势信息，如果初次接触 DFT 计算，并不清楚如何选择赝势和原子轨道文件，较通用的推荐是 [SG15-V1.0_Pseudopotential.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0_Pseudopotential.zip) 赝势文件与 [SG15-V1.0__StandardOrbitals-V2.0.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__StandardOrbitals-Version2p0.zip) 原子轨道文件。**

在[ABACUS 网站](http://abacus.ustc.edu.cn/pseudo/list.htm)上，你可以找到赝势文件以及与ABACUS相对应的优化原子基组。

尽管本页面上呈现的基组经过了仔细测试，但我们并不提供任何保证。用户在将其应用于实际应用之前，应进行充分测试。

目前，ABACUS支持UPF格式的守恒赝势，这是Quantum ESPRESSO的标准赝势格式。特别是，ABACUS支持M. Schlipf和F. Gygi开发的SG15赝势。

SG15优化的Norm-Conserving Vanderbilt (ONCV)多投影仪赝势，在此网站正式发布：
http://www.quantum-simulation.org/potentials/sg15_oncv/ , 
包括周期表中的大部分元素。有关赝势的构造和测试的详细信息，请参阅：[M. Schlipf, F. Gygi, Comp. Phys. Comm. 196, 36 (2015) . ](http://dx.doi.org/10.1016/j.cpc.2015.05.011)此赝势在Creative Commons Attribution-ShareAlike 4.0 International License.下获得许可。

你可以通过点击以下链接下载赝势和基组的压缩包：

[Dojo-NC-FR](https://github.com/abacusmodeling/ABACUS-orbitals/tree/main/Dojo-NC-FR) : Dojo守恒全相对论赝势和优化DZP原子轨道。

[SG15-V1.0_Pseudopotential.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0_Pseudopotential.zip) : 适用于PBE GGA功能的ONCV类型多投影仪赝势。

[SG15-V1.0__StandardOrbitals-V2.0.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__StandardOrbitals-Version2p0.zip) : 使用“PyTorch-Gradient”最小化方法和结合参考波函数梯度[1]的SG15-V1.0赝势的分层优化标准原子轨道。这些轨道更精确，适用于DFT+U计算。

[SG15-V1.0__AllOrbitals-V2.0.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__AllOrbitals-Version2p0.zip) : 使用“PyTorch-Gradient”最小化方法和结合参考波函数梯度[1]的SG15-V1.0赝势的分层优化所有原子轨道。这些轨道更精确，适用于DFT+U计算。

[SG15-V1.0__Orbitals-V1.0.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__Orbitals-Version1p0.zip) : 使用“模拟退火”方法[2]为SG15-V1.0赝势生成的旧原子轨道。在上述版本2.0原子轨道基组发布之前，这些轨道得到了广泛使用。

“重要提示：请注意，原子基组是针对特定赝势进行优化的，这意味着当赝势发生变化时，必须重新生成原子基组。当使用原子基组时，应选择参考能量最接近计算中所用能量的基组。”

我们还为一些旧的赝势优化了原子基组，只有单投影。可用元素如下所示。你可以从以下地址下载赝势和相应的原子基组：[Traditional_PP_Orb.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Traditional_PP_Orb.zip)



##### 从 0 开始设置赝势信息

从官网下载 [SG15-V1.0_Pseudopotential.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0_Pseudopotential.zip) 赝势文件并解压：


```bash
%%bash
# 下载
wget http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0_Pseudopotential.zip
# 解压
unzip SG15-Version1p0_Pseudopotential.zip

# 将下载的赝势文件移动到 ABACUS 的PP_ORB文件夹中，我们已经准备好了相关文件，因此这一步可以省略
# mv ./SG15-Version1p0_Pseudopotential/SG15_ONCV_v1.0_upf/ ./ABACUS/PP_ORB
```

从官网下载 [SG15-V1.0__StandardOrbitals-V2.0.zip](http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__StandardOrbitals-Version2p0.zip) 原子轨道文件并解压：


```bash
%%bash
# 下载
wget http://abacus.ustc.edu.cn/_upload/tpl/0c/d8/3288/template3288/download/Libs/SG15-Version1p0__StandardOrbitals-Version2p0.zip
# 解压
unzip SG15-Version1p0__StandardOrbitals-Version2p0.zip

# 将下载的原子轨道文件移动到 ABACUS 的PP_ORB文件夹中，我们已经准备好了相关文件，因此这一步可以省略
# mv ./SG15-Version1p0__StandardOrbitals-Version2p0/ ./ABACUS/PP_ORB/
```

**注意再次确认 INPUT 文件中赝势文件与原子轨道文件的路径设置正确。**

### 1.2 运行计算

**进入已经准备好输入文件的文件夹**

执行以下命令即可快速运行 ABACUS 计算：


```bash
%%bash
# 进入工作文件夹
cd ./ABACUS_SCF/MgO_LCAO/SCF/
# OMP_NUM_THREADS=1 表示使用单线程，如果你的机器配置比较高，可以使用多线程，比如 4 线程，就可以写成 OMP_NUM_THREADS=4
# mpirun -n 后面的数字表示计算所使用的 CPU 核心数，这里使用 2 个核心，你可以根据你的机器配置进行修改。
OMP_NUM_THREADS=1 mpirun -n 2 abacus
```

                                                                                         
                                  ABACUS v3.3.2
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: e39b50efe (Fri Aug 18 16:14:25 2023 +0800)
    
     Wed Aug 23 13:57:41 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 2 for Mg: [Ne] 3s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 54 * 54 * 54
     UNIFORM GRID DIM(BIG): 18 * 18 * 18
     DONE(0.356784   SEC) : SETUP UNITCELL
     DONE(0.460117   SEC) : SYMMETRY
     DONE(0.608863   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       10              2           112         
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Mg      4s2p1d-8au      15          4           
     O       2s2p1d-7au      13          4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(3.9617     SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.654009e+03  0.000000e+00   1.169e-01  6.376e+00  
     GE2    -7.660990e+03  -6.981374e+00  9.731e-02  5.088e+00  
     GE3    -7.664311e+03  -3.321079e+00  1.647e-02  5.130e+00  
     GE4    -7.664324e+03  -1.312409e-02  2.085e-03  5.029e+00  
     GE5    -7.664326e+03  -1.498631e-03  2.882e-04  5.013e+00  
     GE6    -7.664326e+03  -1.031777e-05  2.253e-05  5.206e+00  
     GE7    -7.664326e+03  -5.778648e-08  5.863e-07  5.036e+00  
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               40.886         9         4.5       1e+02     %
       Driver              driver_line         40.844         1         41        1e+02     %
       ORB_control         set_orb_tables      2.8855         1         2.9       7.1       %
       ORB_gen_tables      gen_tables          2.8855         1         2.9       7.1       %
       ORB_table_phi       init_Table          1.9757         1         2         4.8       %
       ORB_table_phi       cal_ST_Phi12_R      1.961          278       0.0071    4.8       %
       ORB_table_beta      init_Table_Beta     0.68329        1         0.68      1.7       %
       ORB_table_beta      VNL_PhiBeta_R       0.67974        120       0.0057    1.7       %
       Ions                opt_ions            37.068         1         37        91        %
       ESolver_KS_LCAO     Run                 37.06          1         37        91        %
       ESolver_KS_LCAO     beforescf           0.16367        1         0.16      0.4       %
       Potential           update_from_charge  0.45088        8         0.056     1.1       %
       Potential           cal_v_eff           0.4482         8         0.056     1.1       %
       PW_Basis            real2recip          0.12167        62        0.002     0.3       %
       PotXC               cal_v_eff           0.39044        8         0.049     0.95      %
       XC_Functional       v_xc                0.38922        8         0.049     0.95      %
       HSolverLCAO         solve               36.288         7         5.2       89        %
       HamiltLCAO          updateHk            19.255         70        0.28      47        %
       OperatorLCAO        init                17.339         140       0.12      42        %
       Veff                contributeHk        17.338         70        0.25      42        %
       Gint_interface      cal_gint            32.166         14        2.3       79        %
       Gint_interface      cal_gint_vlocal     17.009         7         2.4       42        %
       Gint_Tools          cal_psir_ylm        4.601          40824     0.00011   11        %
       Gint_k              folding_vl_k        0.32851        70        0.0047    0.8       %
       Gint_k              Distri              0.24767        70        0.0035    0.61      %
       Nonlocal<LCAO>      contributeHR        0.95627        1         0.96      2.3       %
       LCAO_gen_fixedH     b_NL_mu_new         0.95583        1         0.96      2.3       %
       OperatorLCAO        folding_fixed       0.76567        70        0.011     1.9       %
       LCAO_nnr            folding_fixedH      0.76326        70        0.011     1.9       %
       HSolverLCAO         hamiltSolvePsiK     0.88287        70        0.013     2.2       %
       DiagoElpa           elpa_solve          0.76798        70        0.011     1.9       %
       ElecStateLCAO       psiToRho            16.15          7         2.3       40        %
       LCAO_Charge         cal_dk_k            0.68474        7         0.098     1.7       %
       Gint_interface      cal_gint_rho        15.157         7         2.2       37        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Wed Aug 23 13:57:41 2023
     FINISH Time  : Wed Aug 23 13:58:22 2023
     TOTAL  Time  : 41
     SEE INFORMATION IN : OUT.MgO/
    

运行过程中的输出日志如下：

```log
                                                                                     
                              ABACUS v3.3.2

               Atomic-orbital Based Ab-initio Computation at UStc                    

                     Website: http://abacus.ustc.edu.cn/                             
               Documentation: https://abacus.deepmodeling.com/                       
                  Repository: https://github.com/abacusmodeling/abacus-develop       
                              https://github.com/deepmodeling/abacus-develop         
                      Commit: 0d06d5b (Tue Aug 15 22:43:59 2023 +0800)

 Fri Aug 18 17:01:21 2023
 MAKE THE DIR         : OUT.MgO/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Warning: the number of valence electrons in pseudopotential > 2 for Mg: [Ne] 3s2
 Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
 If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 UNIFORM GRID DIM     : 54 * 54 * 54
 UNIFORM GRID DIM(BIG): 18 * 18 * 18
 DONE(0.253323   SEC) : SETUP UNITCELL
 DONE(0.655541   SEC) : SYMMETRY
 DONE(0.823616   SEC) : INIT K-POINTS
 ---------------------------------------------------------
 Self-consistent calculations for electrons
 ---------------------------------------------------------
 SPIN    KPOINTS         PROCESSORS  NBASE       
 1       10              12          112         
 ---------------------------------------------------------
 Use Systematically Improvable Atomic bases
 ---------------------------------------------------------
 ELEMENT ORBITALS        NBASE       NATOM       XC          
 Mg      4s2p1d-8au      15          4           
 O       2s2p1d-8au      13          4           
 ---------------------------------------------------------
 Initial plane wave basis and FFT box
 ---------------------------------------------------------
 -------------------------------------------
 SELF-CONSISTENT : 
 -------------------------------------------
 START CHARGE      : atomic
 DONE(6.59195    SEC) : INIT SCF
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
 GE1    -7.653810e+03  0.000000e+00   1.172e-01  3.174e+00  
 GE2    -7.660698e+03  -6.887227e+00  9.635e-02  2.919e+00  
 GE3    -7.663893e+03  -3.195615e+00  1.673e-02  2.894e+00  
 GE4    -7.663908e+03  -1.521795e-02  2.058e-03  2.909e+00  
 GE5    -7.663910e+03  -1.464449e-03  2.765e-04  2.900e+00  
 GE6    -7.663910e+03  -1.052007e-05  2.365e-05  2.884e+00  
 GE7    -7.663910e+03  -6.828454e-08  4.643e-07  2.947e+00  

  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                       total               27.253         9         3         1e+02     %
   Driver              driver_line         27.239         1         27        1e+02     %
   PW_Basis            setuptransform      0.18666        1         0.19      0.68      %
   ORB_control         read_orb_first      0.13776        1         0.14      0.51      %
   LCAO_Orbitals       Read_Orbitals       0.13772        1         0.14      0.51      %
   ORB_control         set_orb_tables      5.2963         1         5.3       19        %
   ORB_gen_tables      gen_tables          5.2963         1         5.3       19        %
   ORB_table_phi       init_Table          3.8418         1         3.8       14        %
   ORB_table_phi       cal_ST_Phi12_R      3.8253         278       0.014     14        %
   ORB_table_beta      init_Table_Beta     1.105          1         1.1       4.1       %
   ORB_table_beta      VNL_PhiBeta_R       1.1004         120       0.0092    4         %
   Ions                opt_ions            20.816         1         21        76        %
   ESolver_KS_LCAO     Run                 20.814         1         21        76        %
   ESolver_KS_LCAO     beforescf           0.15939        1         0.16      0.58      %
   Potential           update_from_charge  0.1562         8         0.02      0.57      %
   Potential           cal_v_eff           0.15501        8         0.019     0.57      %
   PotXC               cal_v_eff           0.13665        8         0.017     0.5       %
   XC_Functional       v_xc                0.13574        8         0.017     0.5       %
   HSolverLCAO         solve               20.341         7         2.9       75        %
   HamiltLCAO          updateHk            10.338         70        0.15      38        %
   OperatorLCAO        init                9.1455         140       0.065     34        %
   Veff                contributeHk        9.144          70        0.13      34        %
   Gint_interface      cal_gint            15.224         14        1.1       56        %
   Gint_interface      cal_gint_vlocal     8.0903         7         1.2       30        %
   Gint_Tools          cal_psir_ylm        1.6025         9072      0.00018   5.9       %
   Gint_k              folding_vl_k        1.0535         70        0.015     3.9       %
   Gint_k              Distri              0.92132        70        0.013     3.4       %
   Nonlocal<LCAO>      contributeHR        0.20092        1         0.2       0.74      %
   LCAO_gen_fixedH     b_NL_mu_new         0.20061        1         0.2       0.74      %
   OperatorLCAO        folding_fixed       0.92411        70        0.013     3.4       %
   LCAO_nnr            folding_fixedH      0.92311        70        0.013     3.4       %
   HSolverLCAO         hamiltSolvePsiK     1.1197         70        0.016     4.1       %
   DiagoElpa           elpa_solve          1.0458         70        0.015     3.8       %
   ElecStateLCAO       psiToRho            8.883          7         1.3       33        %
   LCAO_Charge         cal_dk_k            0.95741        7         0.14      3.5       %
   Gint_interface      cal_gint_rho        7.1334         7         1         26        %
 ----------------------------------------------------------------------------------------

 START  Time  : Fri Aug 18 17:01:21 2023
 FINISH Time  : Fri Aug 18 17:01:48 2023
 TOTAL  Time  : 27
 SEE INFORMATION IN : OUT.MgO/
```

让我们来分析一下运行日志，这个 ABACUS 运行输出日志可以分为以下几个部分：

1. 软件信息和时间戳
2. 输入参数和警告信息
3. 计算设置信息
4. 自洽场（SCF）迭代过程
5. 各模块运行时间分析

下面是对各部分的详细分析：

1. 软件信息和时间戳

这部分显示了 ABACUS 的版本信息、官方网站、文档网站、代码仓库以及当前使用的代码提交版本。同时，还给出了计算开始的时间。

2. 输出文件夹和警告信息

这部分首先给出了输出文件夹的名称（OUT.MgO/），然后显示了一个关于镁（Mg）原子赝势中价电子数大于 2 的警告。这说明用户选择了一个具有额外电子的赝势，可能会导致更准确的结果，但计算效率可能较低。如果用户对所选赝势有信心，可以忽略此警告。

3. 计算设置信息

这部分给出了计算所需的各种参数和设置。首先是均匀网格的维度（54 * 54 * 54 和 18 * 18 * 18），

接着是晶胞、对称性和 K 点的初始化时间。

然后列出了自洽场计算的设置，包括自旋、 K 点、处理器数量和原子基组数量。

接下来是原子和交换相关泛函的信息，

最后是平面波基组和 FFT 网格的初始化信息。

4. 自洽场（SCF）迭代过程

这部分展示了自洽场计算的迭代过程。首先给出了初始电荷密度的来源（原子），然后列出了每一步迭代的总能量（ETOT）、能量差（EDIFF）、电荷密度差（DRHO）和计算时间。最后一次迭代的能量差和电荷密度差已经非常小，说明自洽场计算已经收敛。

5. 各模块运行时间分析

    在这部分输出日志中，主要展示了程序中各个子模块的性能分析，包括每个子模块的运行时间、调用次数等信息。这有助于了解整个程序的运行状况和性能瓶颈所在。日志按照不同的类别（CLASS_NAME）、名称（NAME）、时间（TIME(Sec)）、调用次数（CALLS）、平均耗时（AVG）和占总时间百分比（PER%）进行了汇总。以下是一些关键子模块的分析：

    1. Driver: 是整个程序的驱动模块，负责执行整个流程。它的运行时间为27.239秒，占总时间的100%。

    2. PW_Basis: 模块负责设置平面波基组和FFT（快速傅里叶变换）格点。它的运行时间为0.19秒，占总时间的0.68%。

    3. ORB_control: 模块负责读取和设置轨道表格。它的运行时间为5.3秒，占总时间的19%。

    4. Ions: 模块负责优化离子。它的运行时间为21秒，占总时间的76%。

    5. ESolver_KS_LCAO: 模块负责Kohn-Sham方程的数值求解。它的运行时间为21秒，占总时间的76%。

    6. HSolverLCAO: 模块负责求解哈密顿矩阵。它的运行时间为20.341秒，占总时间的75%。

    7. Gint_interface: 模块负责计算G积分。它的运行时间为15.224秒，占总时间的56%。

    从上述分析可以看出，Ions、ESolver_KS_LCAO、HSolverLCAO 和 Gint_interface 这几个子模块的运行时间较长，占据了程序总运行时间的大部分。这些模块可能是程序的性能瓶颈所在，可以针对这些模块进行优化以提高程序的运行效率。

6. 开始时间、结束时间、整体耗时与输出文件夹

### 1.3 计算结果与输出文件

输出文件储存在当前文件夹下的 `OUT.MgO` 文件夹中，输出文件夹的名称为 `OUT.<suffix>`，其中 `<suffix>` 由 INPUT 文件中的 suffix 值指定。

计算结果主要储存在 `OUT.MgO/running_scf.log` 文件中，它的开始部分如下所示：

```log
                                                                                     
                              ABACUS v3.3.2

               Atomic-orbital Based Ab-initio Computation at UStc                    

                     Website: http://abacus.ustc.edu.cn/                             
               Documentation: https://abacus.deepmodeling.com/                       
                  Repository: https://github.com/abacusmodeling/abacus-develop       
                              https://github.com/deepmodeling/abacus-develop         
                      Commit: 0d06d5b (Tue Aug 15 22:43:59 2023 +0800)

    Start Time is Fri Aug 18 17:01:21 2023
                                                                                     
 ------------------------------------------------------------------------------------

 READING GENERAL INFORMATION
                           global_out_dir = OUT.MgO/
                           global_in_card = INPUT
                               pseudo_dir = 
                              orbital_dir = 
                                    DRANK = 1
                                    DSIZE = 12
                                   DCOLOR = 1
                                    GRANK = 1
                                    GSIZE = 1
 The esolver type has been set to : ksdft_lcao




 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 |                                                                    |
 | Reading atom information in unitcell:                              |
 | From the input file and the structure file we know the number of   |
 | different elments in this unitcell, then we list the detail        |
 | information for each element, especially the zeta and polar atomic |
 | orbital number for each element. The total atom number is counted. |
 | We calculate the nearest atom distance for each atom and show the  |
 | Cartesian and Direct coordinates for each atom. We list the file   |
 | address for atomic orbitals. The volume and the lattice vectors    |
 | in real and reciprocal space is also shown.                        |
 |                                                                    |
 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 ......
 
```

如果 ABACUS 自洽场计算成功完成，能量会输出在 `OUT.MgO/running_scf.log` 文件中:

```
 --------------------------------------------
 !FINAL_ETOT_IS -7663.897267807250 eV
 --------------------------------------------
```

## 2. ABACUS 自洽场计算（PW）

除了使用 LCAO 基组进行计算之外，我们也可以使用平面波基组进行 SCF 计算。

为了使用PW（平面波）基组进行SCF计算，只需**在 INPUT 文件中将 basis_type 标签从 lcao 更改为 pw，而无需在 STRU 文件中的NUMERICAL_ORBITAL 下提供轨道文件。**

### 2.1 输入文件

#### 2.1.1 INPUT 文件

```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
pseudo_dir              ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              pw              # changes the type of basis set
calculation             scf		# this is the key parameter telling abacus to do a scf
```

#### 2.1.2 STRU 文件

```
#This is the atom file containing all the information
#about the lattice structure.

ATOMIC_SPECIES
Mg 24.305  Mg_ONCV_PBE-1.0.upf  # element name, atomic mass, pseudopotential file
O  15.999 O_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.8897259886 		# 1.8897259886 Bohr =  1.0 Angstrom

LATTICE_VECTORS
4.25648 0.00000 0.00000  
0.00000 4.25648 0.00000
0.00000 0.00000 4.25648

ATOMIC_POSITIONS
Direct                  #Cartesian(Unit is LATTICE_CONSTANT)
Mg                      #Name of element        
0.0                     #Magnetic for this element.
4                       #Number of atoms
0.0  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
O                       #Name of element        
0.0                     #Magnetic for this element.
4                       #Number of atoms
0.5  0.0  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
0.5  0.5  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.0  0.5  0 0 0    #x,y,z, move_x, move_y, move_z
0.0  0.5  0.0  0 0 0    #x,y,z, move_x, move_y, move_z
```

#### 2.1.3 KPT 文件

和第 1 节中的 KPT 文件相同

```
K_POINTS
0           # k点的总数，`0'表示自动生成
Gamma       # Monkhorst-Pack方法的类型，`Gamma'或`MP'
4 4 4 0 0 0 # 前三个数字：沿着倒数向量的细分
            # 后三个数字：网格的偏移
```


### 2.2 运行


**进入已经准备好输入文件的文件夹**

执行以下命令运行 ABACUS 计算：


```bash
%%bash
# 进入工作文件夹
cd ./ABACUS_SCF/MgO_PW/SCF/
# OMP_NUM_THREADS=1 表示使用单线程，如果你的机器配置比较高，可以使用多线程，比如 4 线程，就可以写成 OMP_NUM_THREADS=4
# mpirun -n 后面的数字表示计算所使用的 CPU 核心数，这里使用 2 个核心，你可以根据你的机器配置进行修改。
OMP_NUM_THREADS=1 mpirun -n 2 abacus
```


### 2.3 结果

结果分析同第 1 小节，最终获得的能量为：

```
 --------------------------------------------
 !FINAL_ETOT_IS -7665.688319476949 eV
 --------------------------------------------
```

## 3. 练习

请根据以上示例，计算 Ti 金属单质的基态能量。

### 3.1 准备输入文件

已知 Ti 在常温下为 hcp 结构，Ti 的相对原子质量为：47.867。原子位置信息如下：

```POSCAR
Ti
1.0  # 晶胞常数
   2.2836874152833575   -3.9554626318764212    0.0000000000000000  # 晶胞向量
   2.2836874152833575    3.9554626318764212    0.0000000000000000  # 晶胞向量
   0.0000000000000000    0.0000000000000000    2.8262442700000001  # 晶胞向量
Ti  # 元素种类
3  # 原子个数
direct  # 坐标系
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Ti  # 原子 1 位置
   0.3333333333333333    0.6666666666666666    0.5000000000000000 Ti  # 原子 2 位置
   0.6666666666666667    0.3333333333333334    0.5000000000000000 Ti  # 原子 3 位置
```

有关直接从 POSCAR 文件转换为 STRU 文件，可参考：[ABACUS 使用教程｜如何转换 STRU 文件](https://nb.bohrium.dp.tech/detail/1311634750)


```python
import os

# 创建 STRU 文件
stru = """
/* 删掉本句话，然后在这里输入你的 STRU 文件内容 */
"""
with open('./ABACUS_SCF/Practice/STRU', 'w') as f:
    f.write(stru)
    
# 创建 INPUT 文件
os.mkdir('./ABACUS_SCF/Practice', is_exist_ok=True)
input = """
/* 删掉本句话，然后在这里输入你的 INPUT 文件内容 */
"""
with open('./ABACUS_SCF/Practice/INPUT', 'w') as f:
    f.write(input)
    
# 创建 KPT 文件
kpt = """
/* 删掉本句话，然后在这里输入你的 KPT 文件内容 */
"""
with open('./ABACUS_SCF/Practice/KPT', 'w') as f:
    f.write(kpt)
```

### 3.2 运行计算


```bash
%%bash
# 在这里输入你的运行脚本
```

### 3.3 结果分析

找到输出结果文件，获取你得到的能量信息。


```bash
%%bash
# 修改 OUT.<your_system_name>
cat ./ABACUS_SCF/Practice/OUT.<your_system_name>/running_scf.log|grep "FINAL_ETOT_IS"
```

## 4. 本章小节

到这里，你已经基本了解了使用 ABACUS 的方法并进行了基础的 SCF 计算。接下来，请带着你的好奇心继续阅读下一小节，进一步了解 DFT 计算的另一大应用场景——[晶格弛豫｜结构优化｜几何优化｜以上三个都是同一种意思]()。


```python
# Author Information
import ipywidgets as widgets
import webbrowser

author ='''Haohui Que<br/>
'''
reference = '''Reference: <br/>
[1] https://abacus.deepmodeling.com/en/latest/
'''

license = '''License: Apache License 2.0 <br/>
'''

preamble = widgets.Accordion(
    children=[
        widgets.HTML(value=author),
        widgets.HTML(value=reference),
        widgets.HTML(value=license),
        ], 
    selected_index=0
    )
preamble.set_title(0, '作者')
preamble.set_title(1, 'Reference')
preamble.set_title(2, 'License')

display(preamble)
```
