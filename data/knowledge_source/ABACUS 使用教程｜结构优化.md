# ABACUS 使用教程｜结构优化/晶格弛豫/几何优化

在上一章节中，我们为大家带来了 [ABACUS 系列｜电子自洽迭代](https://nb.bohrium.dp.tech/detail/7417640496)，介绍了基础的 ABACUS 的输入输出文件与基础操作。在本章节中，我们将为大家带来 ABACUS 的另一大应用场景——结构优化（又称为晶格弛豫/几何优化等）。

在第一性原理计算中，晶格弛豫（Lattice relaxation）指的是一种过程，即通过调整晶格内原子的位置和晶格参数（如晶格常数），使晶体结构达到最低能量状态。

晶格弛豫的主要目的是找到晶体结构的稳定构型，这对于预测和理解材料的性能至关重要。在这个过程中，计算机模拟会根据所采用的理论框架来迭代地优化原子位置和晶格参数，最终得到一个能量最低的结构。这个结构通常被认为是晶体在实际条件下所处的状态。

在材料科学和凝聚态物理研究中，晶格弛豫已经成为了研究新材料和预测其性质的重要工具。通过这些计算方法，研究人员可以在实验之前获得对材料性能的深入了解，从而为实验设计提供有价值的指导。

让我们来看一段视频，更直观地了解电子自洽迭代与结构优化：

> 视频出处：[4分钟了解什么是 AI for Science | 深势科技](https://www.bilibili.com/video/BV1c44y1c7hJ/?spm_id_from=333.337.search-card.all.click&vd_source=4de0e1da1c75a9593bb4be49f9c341c8)


```python
# Video of Relaxation. Source: DP Technology
from IPython.display import Video

video_path = "/bohr/abacus02-ecia/v4/ABACUS_Relax/DFT.mp4"
video = Video(video_path, width=480)
display(video)
```


<video src="/bohr/abacus02-ecia/v4/ABACUS_Relax/DFT.mp4" controls  width="480" >
      Your browser does not support the <code>video</code> element.
    </video>


> *若视频无法在线播放，可点击右侧数据集下载后播放。

直观来说，上一章节中的电子自洽迭代我们解决了“某一帧”中原子之间相互作用的问题，那么结构优化/晶格弛豫就是在计算得到的相互作用下，进一步计算原子的移动，再计算新结构下的相互作用…… 循环往复并最终得到稳定结构的过程。**这一步又被称为“离子步”。**

**本篇 Notebook 将先给出一个结构优化操作示例。并以 $H_2$ 分子键长计算作为练习。**

## 结构优化计算示例

### 输入文件

与上一章节的[电子自洽迭代](https://nb.bohrium.dp.tech/detail/7417640496)中类似，我们需要准备以下五个类型地输入文件：

- `INPUT`：包含了计算过程中所需的各种参数，定义和控制计算任务；
- `STRU`：结构文件，包含了原子种类、原子位置、晶格常数以及晶格向量等信息；
- `KPT`：包含了布里渊区积分所需的k点信息；
- `*.upf`：包含了原子的赝势信息；
- `*.orb`：包含了原子轨道的数值表示；

其中，我们**只需要调整 INPUT 文件中的设置项**，其它文件与上一章节中的输入文件相同。

在示例数据集中，我们为你准备好了结构优化的示例文件。我们可以直接访问：（你可以在左侧点击数据集查看相应文件）：


```python
! tree /bohr/
```

    [01;34m/bohr/[0m
    └── [01;34mabacus02-ecia[0m
        └── [01;34mv4[0m
            └── [01;34mABACUS_Relax[0m
                ├── [01;32mDFT.mp4[0m
                ├── [01;34mMgO_LCAO[0m
                │   └── [01;34moptimization[0m
                │       ├── [01;32mINPUT[0m
                │       ├── [01;32mKPT[0m
                │       └── [01;32mSTRU[0m
                ├── [01;34mMgO_PW[0m
                │   └── [01;34moptimization[0m
                │       ├── [01;32mINPUT[0m
                │       ├── [01;32mKPT[0m
                │       └── [01;32mSTRU[0m
                ├── [01;34mPP_ORB[0m
                │   ├── [01;32mMg_ONCV_PBE-1.0.upf[0m
                │   ├── [01;32mMg_gga_8au_100Ry_4s2p1d.orb[0m
                │   ├── [01;32mO_ONCV_PBE-1.0.upf[0m
                │   └── [01;32mO_gga_7au_100Ry_2s2p1d.orb[0m
                └── [01;34mPractice[0m
                    └── [01;34mPP_ORB[0m
                        ├── [01;32mH_ONCV_PBE-1.0.upf[0m
                        └── [01;32mH_gga_6au_100Ry_2s1p.orb[0m
    
    10 directories, 13 files
    

出于安全考虑，我们没有数据集所在文件夹的写入权限，因此我们将其复制到 `/data/` 目录下:


```python
! cp -nr /bohr/ /personal/
```

我们在这里定义一些路径，并切换到工作路径，方便后续调用：


```python
import os

bohr_dataset_url = "/bohr/abacus02-ecia/v4/"  # url 可从左侧数据集复制
work_path = os.path.join("/personal", bohr_dataset_url[1:])
os.chdir(work_path)
print(f"当前路径为：{os.getcwd()}")
```

    当前路径为：/personal/bohr/abacus02-ecia/v4
    

#### INPUT 文件

为了在ABACUS中进行完整的几何优化，需要在INPUT中将 `calculation`标签设置为`cell-relax`。

此外，原子力和晶胞应力的收敛标准分别可以通过标签`force_thr_ev`和`stress_thr`来设置。

离子步骤的最大数量由`relax_nmax`控制。

同样的，结构优化也可以选择使用 LCAO 基组或 PW 基组。

* INPUT (LCAO 基组)

以下是一个使用 LCAO 基组的 INPUT 文件示例：

```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
orbital_dir             ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              lcao 
calculation             cell-relax	# this is the key parameter telling abacus to do a optimization calculation
force_thr_ev		0.01		# the threshold of the force convergence, in unit of eV/Angstrom
stress_thr		2		# the threshold of the stress convergence, in unit of kBar
relax_nmax		100		# the maximal number of ionic iteration steps
out_stru		1
```

1. suffix: 用于区分同一文件夹下不同计算的后缀，通常与计算对象有关。

2. ntype: 表示体系中原子种类的数量，这里的值为2，表示有两种原子（例如Mg和O）。

3. nelec: 表示体系中电子的数量，这里的值为0.0，表示电子数量由原子种类和价态自动计算得出。

4. pseudo_dir: 赝势文件的路径，这里设置为"./"，表示赝势文件位于当前文件夹。

5. orbital_dir: 轨道文件的路径，这里设置为"./"，表示轨道文件位于当前文件夹。

6. ecutwfc: 波函数截断能，用于确定平面波基组的大小。这里的值为100 Rydberg，表示波函数截断能为100 Ry。

7. scf_thr: 自洽场（SCF）计算的能量收敛阈值，这里设为1e-4 Rydberg，表示能量收敛到1e-4 Ry时计算停止。

8. basis_type: 表示基组类型，这里设置为"lcao"，表示使用 LCAO 基组。

9. calculation: 计算类型，这里设置为"cell-relax"，表示进行晶胞优化计算。

10. force_thr_ev: 力收敛阈值，用于判断优化过程是否收敛。这里设置为0.01 eV/Angstrom，表示当所有原子上的力小于0.01 eV/Angstrom时，认为已经收敛。默认值：0.0257112 eV/Angstrom。(你也可以使用 force_thr 参数（默认值：0.001 Ry/Bohr），两者只是单位不同。)

11. stress_thr: 应力收敛阈值，用于判断晶胞优化过程中应力是否收敛。这里设置为 2 kBar，表示当晶胞内的应力小于 2 kBar 时，认为已经收敛。默认值：0.5 kBar

12. relax_nmax: 优化过程中的最大迭代次数，这里设置为100，表示最多进行100次迭代。

13. out_stru: 表示是否输出优化后的结构信息。这里设置为1，表示在优化结束后将优化后的结构信息输出到文件。

---

* INPUT (PW 基组)

以下是一个使用 PW 基组的 INPUT 文件示例：

```
INPUT_PARAMETERS
suffix                  MgO
ntype                   2
nelec                   0.0
pseudo_dir              ./
ecutwfc                 100             # Rydberg
scf_thr                 1e-4		# Rydberg
basis_type              pw
calculation             cell-relax	# this is the key parameter telling abacus to do a optimization calculation
force_thr_ev		0.01		# the threshold of the force convergence, in unit of eV/Angstrom
stress_thr		2		# the threshold of the stress convergence, in unit of kBar
relax_nmax		100		# the maximal number of ionic iteration steps
out_stru		1
```



#### STRU 文件

```
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

##### 参数解释

- ATOMIC_SPECIES

从左到右依次为原子种类、相对原子质量、赝势文件

- NUMERICAL_ORBITAL

原子轨道文件

* 赝势文件 `Mg_ONCV_PBE-1.0.upf` 和 `O_ONCV_PBE-1.0.upf` 应放在 `pseudo_dir` 目录下，
* 轨道文件 `Mg_gga_8au_100Ry_4s2p1d.orb` 和 `O_gga_8au_100Ry_2s2p1d.orb` 应放在 `orbital_dir` 目录下。

`pseudo_dir` 和 `orbital_dir` 文件路径在 `INPUT` 文件中进行设置，赝势和轨道文件可以从 [ABACUS 网站](http://abacus.ustc.edu.cn/pseudo/list.htm) 下载。

- LATTICE_CONSTANT

晶格常量。

- LATTICE_VECTORS

晶格向量，即晶胞的 x,y,z 轴长度。

- ATOMIC_POSITIONS

原子位置信息

有关更多 STRU 文件的参数信息，请参阅：[The STRU file](https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html)

---

使用 PW 基组时，删去 STRU 文件中的 `NUMERICAL_ORBITAL` 部分，其它与以上示例一致。

#### KPT 文件

```
K_POINTS
0 
Gamma
4 4 4 0 0 0
```

#### 赝势和原子轨道文件

* 赝势文件 `Mg_ONCV_PBE-1.0.upf` 和 `O_ONCV_PBE-1.0.upf` 应放在 `pseudo_dir` 目录下，
* 轨道文件 `Mg_gga_8au_100Ry_4s2p1d.orb` 和 `O_gga_7au_100Ry_2s2p1d.orb` 应放在 `orbital_dir` 目录下。

`pseudo_dir` 和 `orbital_dir` 文件路径在 `INPUT` 文件中进行设置，赝势和轨道文件可以从 [ABACUS 网站](http://abacus.ustc.edu.cn/pseudo/list.htm) 下载。

### 运行
**进入已经准备好输入文件的文件夹**，执行以下命令即可快速运行 ABACUS 计算：


```bash
%%bash
# 进入工作文件夹
cd ./ABACUS_Relax/MgO_LCAO/optimization
# OMP_NUM_THREADS=1 表示使用单线程，如果你的机器配置比较高，可以使用多线程，比如 4 线程，就可以写成 OMP_NUM_THREADS=4
# mpirun -n 后面的数字表示计算所使用的 CPU 核心数，这里使用 2 个核心，你可以根据你的机器配置进行修改。
OMP_NUM_THREADS=2 mpirun -n 2 abacus
```

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WARNING: Total thread number(4) is larger than hardware availability(2).
    %% WARNING: The results may be INCORRECT. Please be sure what you are doing.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                                                                                         
                                  ABACUS v3.3.2
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: e39b50efe (Fri Aug 18 16:14:25 2023 +0800)
    
     Thu Feb 29 12:05:18 2024
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 2 for Mg: [Ne] 3s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 54 * 54 * 54
     UNIFORM GRID DIM(BIG): 18 * 18 * 18
     DONE(0.174909   SEC) : SETUP UNITCELL
     DONE(0.291938   SEC) : SYMMETRY
     DONE(0.456703   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Cell relaxation calculations
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
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(5.33552    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.654024e+03  0.000000e+00   1.168e-01  8.884e+00  
     GE2    -7.661008e+03  -6.984143e+00  9.703e-02  7.191e+00  
     GE3    -7.664311e+03  -3.302513e+00  1.644e-02  7.122e+00  
     GE4    -7.664324e+03  -1.310531e-02  2.066e-03  7.428e+00  
     GE5    -7.664326e+03  -1.473753e-03  2.860e-04  7.452e+00  
     GE6    -7.664326e+03  -1.013484e-05  2.229e-05  7.405e+00  
     GE7    -7.664326e+03  -5.586382e-08  5.773e-07  7.198e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     4.652e+00      5.816e-16      -2.181e-15     
     5.816e-16      4.652e+00      2.908e-16      
     -2.181e-15     2.908e-16      4.652e+00      
     TOTAL-PRESSURE: 4.652e+00 KBAR
     ETOT DIFF (eV)       : 0.000e+00
     LARGEST GRAD (eV/A)  : 0.000e+00
     DONE(8.428e+01  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 2
     -------------------------------------------
     DONE(8.433e+01  SEC) : LOCAL POTENTIAL
     DONE(8.454e+01  SEC) : SYMMETRY
     DONE(8.454e+01  SEC) : INIT K-POINTS
     DONE(8.494e+01  SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.664326e+03  0.000000e+00   9.340e-05  8.746e+00  
     GE2    -7.664326e+03  -1.863489e-06  2.615e-05  7.222e+00  
     GE3    -7.664326e+03  -9.361320e-08  1.239e-05  7.270e+00  
     GE4    -7.664326e+03  -1.025212e-07  3.036e-06  7.039e+00  
     GE5    -7.664326e+03  3.088937e-09   1.846e-07  7.219e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -2.070e+00     7.332e-19      -1.018e-15     
     7.332e-19      -2.070e+00     -2.896e-16     
     -1.018e-15     -2.896e-16     -2.070e+00     
     TOTAL-PRESSURE: -2.070e+00 KBAR
     ETOT DIFF (eV)       : -3.241e-04
     LARGEST GRAD (eV/A)  : 0.000e+00
     DONE(1.491e+02  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 3
     -------------------------------------------
     DONE(1.492e+02  SEC) : LOCAL POTENTIAL
     DONE(1.494e+02  SEC) : SYMMETRY
     DONE(1.494e+02  SEC) : INIT K-POINTS
     DONE(1.498e+02  SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.664326e+03  0.000000e+00   2.437e-05  8.736e+00  
     GE2    -7.664326e+03  -1.383603e-08  6.737e-06  7.177e+00  
     GE3    -7.664326e+03  3.146168e-09   2.605e-06  7.482e+00  
     GE4    -7.664326e+03  1.314770e-09   8.809e-07  7.293e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -4.350e-01     5.805e-16      5.800e-16      
     5.805e-16      -4.350e-01     2.903e-16      
     5.800e-16      2.903e-16      -4.350e-01     
     TOTAL-PRESSURE: -4.350e-01 KBAR
     ETOT DIFF (eV)       : -3.477e-05
     LARGEST GRAD (eV/A)  : 0.000e+00
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               206.77         29        7.1       1e+02     %
       Driver              driver_line         206.7          1         2.1e+02   1e+02     %
       Symmetry            analy_sys           84.16          2         42        41        %
       ORB_control         read_orb_first      0.30106        1         0.3       0.15      %
       LCAO_Orbitals       Read_Orbitals       0.30102        1         0.3       0.15      %
       NOrbital_Lm         extra_uniform       0.1333         12        0.011     0.064     %
       Mathzone_Add1       Uni_Deriv_Phi       0.11465        12        0.0096    0.055     %
       ORB_control         set_orb_tables      4.1285         1         4.1       2         %
       ORB_gen_tables      gen_tables          4.1285         1         4.1       2         %
       ORB_table_phi       init_Table          2.767          1         2.8       1.3       %
       ORB_table_phi       cal_ST_Phi12_R      2.7318         278       0.0098    1.3       %
       ORB_table_beta      init_Table_Beta     0.98341        1         0.98      0.48      %
       ORB_table_beta      VNL_PhiBeta_R       0.9777         120       0.0081    0.47      %
       Ions                opt_ions            201.74         1         2e+02     98        %
       ESolver_KS_LCAO     Run                 123.4          3         41        60        %
       ESolver_KS_LCAO     beforescf           1.563          3         0.52      0.76      %
       ESolver_KS_LCAO     beforesolver        0.33075        3         0.11      0.16      %
       ESolver_KS_LCAO     set_matrix_grid     0.31241        3         0.1       0.15      %
       Grid_Technique      init                0.17764        3         0.059     0.086     %
       Record_adj          for_2d              0.11388        3         0.038     0.055     %
       Charge              atomic_rho          0.13536        6         0.023     0.065     %
       PW_Basis            recip2real          0.45549        113       0.004     0.22      %
       PW_Basis            gathers_scatterp    0.19563        113       0.0017    0.095     %
       Potential           init_pot            0.28476        3         0.095     0.14      %
       Potential           update_from_charge  1.6927         19        0.089     0.82      %
       Potential           cal_v_eff           1.6682         19        0.088     0.81      %
       H_Hartree_pw        v_hartree           0.24776        19        0.013     0.12      %
       PW_Basis            real2recip          0.78129        164       0.0048    0.38      %
       PW_Basis            gatherp_scatters    0.40169        164       0.0024    0.19      %
       PotXC               cal_v_eff           1.4114         19        0.074     0.68      %
       XC_Functional       v_xc                1.4076         19        0.074     0.68      %
       Symmetry            rho_symmetry        0.27148        19        0.014     0.13      %
       HSolverLCAO         solve               118.59         16        7.4       57        %
       HamiltLCAO          updateHk            61.83          160       0.39      30        %
       OperatorLCAO        init                55.008         320       0.17      27        %
       Veff                contributeHk        54.767         160       0.34      26        %
       Gint_interface      cal_gint            124.02         35        3.5       60        %
       Gint_interface      cal_gint_vlocal     49.175         16        3.1       24        %
       Gint_Tools          cal_psir_ylm        12.821         46656     0.00027   6.2       %
    

       Gint_k              folding_vl_k        5.5915         160       0.035     2.7       %
       Gint_k              Distri              5.2444         160       0.033     2.5       %
       Overlap             contributeHR        0.38136        3         0.13      0.18      %
       LCAO_gen_fixedH     calculate_S_no      0.38135        3         0.13      0.18      %
       Ekin<LCAO>          contributeHR        0.38151        3         0.13      0.18      %
       Nonlocal<LCAO>      contributeHR        3.8553         3         1.3       1.9       %
       LCAO_gen_fixedH     b_NL_mu_new         14.302         6         2.4       6.9       %
       ORB_gen_tables      snap_psibeta_half   0.32913        29614     1.1e-05   0.16      %
       OperatorLCAO        folding_fixed       2.203          160       0.014     1.1       %
       LCAO_nnr            folding_fixedH      2.1881         160       0.014     1.1       %
       HSolverLCAO         hamiltSolvePsiK     7.7372         160       0.048     3.7       %
       DiagoElpa           elpa_solve          7.5961         160       0.047     3.7       %
       ElecStateLCAO       psiToRho            49.019         16        3.1       24        %
       elecstate           cal_dm              0.68966        19        0.036     0.33      %
       psiMulPsiMpi        pdgemm              0.4961         190       0.0026    0.24      %
       LCAO_Charge         cal_dk_k            2.2343         16        0.14      1.1       %
       Gint_interface      cal_gint_rho        41.695         16        2.6       20        %
       Charge              mix_rho             0.16164        13        0.012     0.078     %
       Charge              Pulay_mixing        0.15542        13        0.012     0.075     %
       Force_Stress_LCAO   getForceStress      78.183         3         26        38        %
       Stress_Func         stress_gga          0.12853        3         0.043     0.062     %
       Force_LCAO_k        ftable_k            77.684         3         26        38        %
       Force_LCAO_k        allocate_k          14.601         3         4.9       7.1       %
       Force_LCAO_k        cal_foverlap_k      0.30319        3         0.1       0.15      %
       Force_LCAO_k        cal_edm_2d          0.21397        3         0.071     0.1       %
       Local_Orbital_Chargecal_dm_R            0.23968        6         0.04      0.12      %
       Force_LCAO_k        cal_fvl_dphi_k      33.147         3         11        16        %
       Gint_interface      cal_gint_force      33.147         3         11        16        %
       Gint_Tools          cal_dpsir_ylm       15.983         4374      0.0037    7.7       %
       Gint_Tools          cal_dpsirr_ylm      1.9517         4374      0.00045   0.94      %
       Force_LCAO_k        cal_fvnl_dbeta_k_new28.503         3         9.5       14        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Thu Feb 29 12:05:18 2024
     FINISH Time  : Thu Feb 29 12:08:45 2024
     TOTAL  Time  : 207
     SEE INFORMATION IN : OUT.MgO/
    

### 计算结果与输出文件

输出文件储存在当前文件夹下的 `OUT.MgO` 文件夹中，输出文件夹的名称为 `OUT.<suffix>`，其中 `<suffix>` 由 INPUT 文件中的 suffix 值指定。

最终优化的结构可以在 `STRU_NOW.cif` 和 `OUT.MgO/running_cell-relax.log` 中找到。

**结构优化的收敛要求 一般是原子受力很小、体系的压强小 即平衡。**

让我们来查看一下结构弛豫的标准输出文件（即运行日志文件），来查看结构优化的过程：

其中，每一离子步的内容如下，让我们来查看其中两步（第 2～3 步）。
```
 -------------------------------------------
 STEP OF RELAXATION : 2
 -------------------------------------------
 DONE(3.406e+01  SEC) : LOCAL POTENTIAL
 DONE(3.411e+01  SEC) : SYMMETRY
 DONE(3.412e+01  SEC) : INIT K-POINTS
 DONE(3.432e+01  SEC) : INIT SCF
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
 GE1    -7.664326e+03  0.000000e+00   9.340e-05  3.690e+00  
 GE2    -7.664326e+03  -1.863560e-06  2.615e-05  3.270e+00  
 GE3    -7.664326e+03  -9.364259e-08  1.239e-05  3.270e+00  
 GE4    -7.664326e+03  -1.024454e-07  3.036e-06  3.266e+00  
 GE5    -7.664326e+03  3.054908e-09   1.846e-07  3.270e+00  
 ><><><><><><><><><><><><><><><><><><><><><><
 TOTAL-STRESS (KBAR):
 ><><><><><><><><><><><><><><><><><><><><><><
 -2.070e+00     -7.239e-16     2.031e-15      
 -7.239e-16     -2.070e+00     5.800e-16      
 2.031e-15      5.800e-16      -2.070e+00     
 TOTAL-PRESSURE: -2.070e+00 KBAR
 ETOT DIFF (eV)       : -3.241e-04
 LARGEST GRAD (eV/A)  : 0.000e+00
 DONE(5.932e+01  SEC) : SETUP UNITCELL
 -------------------------------------------
 STEP OF RELAXATION : 3
 -------------------------------------------
 DONE(5.933e+01  SEC) : LOCAL POTENTIAL
 DONE(5.938e+01  SEC) : SYMMETRY
 DONE(5.938e+01  SEC) : INIT K-POINTS
 DONE(5.959e+01  SEC) : INIT SCF
 ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
 GE1    -7.664326e+03  0.000000e+00   2.437e-05  3.711e+00  
 GE2    -7.664326e+03  -1.390563e-08  6.737e-06  3.283e+00  
 GE3    -7.664326e+03  3.350344e-09   2.605e-06  3.298e+00  
 GE4    -7.664326e+03  1.101314e-09   8.809e-07  3.278e+00  
 ><><><><><><><><><><><><><><><><><><><><><><
 TOTAL-STRESS (KBAR):
 ><><><><><><><><><><><><><><><><><><><><><><
 -4.350e-01     2.904e-16      2.757e-15      
 2.904e-16      -4.350e-01     5.804e-16      
 2.757e-15      5.804e-16      -4.350e-01     
 TOTAL-PRESSURE: -4.350e-01 KBAR
 ETOT DIFF (eV)       : -3.477e-05
 LARGEST GRAD (eV/A)  : 0.000e+00
```

每一步的输出内容格式是相同的。

其中 `TOTAL-PRESSURE` 参数为晶胞内应力，对应 INPUT 文件中的 `stress_thr` 参数，

`LARGEST GRAD (eV/A)` 参数为最大原子受力，对应 INPUT 文件中的 `force_thr_ev` 参数，

当 `TOTAL-PRESSURE` 小于 `stress_thr` 参数中设置的阈值时，且 `LARGEST GRAD (eV/A)` 小于 `force_thr_ev` 参数中设置的阈值时，认为结构优化收敛，优化结束。

在本例中，我们设置的 $stress_the = 0.5 KBAR$，在第 2 步中，$TOTAL-PRESSURE = -2.070e+00 KBAR > 0.5 KBAR$，尚未收敛，继续进行下一离子步；在第 3 步中，$TOTAL-PRESSURE: -4.350e-01 KBAR < 0.5 KBAR$，结构优化收敛，优化结束。


### 优化算法

ABACUS中可变单元格弛豫的当前实现遵循嵌套过程：首先进行固定单元格结构弛豫，然后更新单元格参数，重复该过程直到收敛。

在上述嵌套过程中，我们使用了共轭梯度（CG）方法进行单元格弛豫，同时为固定单元格结构弛豫提供了四种不同的算法：BFGS，SD（最陡下降），CG（共轭梯度）以及混合CG-BFGS方法。可以使用关键字 [relax_method](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#relax-method) 选择优化算法。我们还提供了一个用于控制弛豫过程的[关键字列表](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#geometry-relaxation)。 

* BFGS方法 

[BFGS方法](https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm) 是一种求解非线性优化问题的拟牛顿方法。它属于在优化过程中对海森矩阵进行近似的拟牛顿方法。如果初始点离极值点不远，BFGS方法往往比基于梯度的方法效果更好。在ABACUS中，我们实现了用于固定单元格结构弛豫的BFGS方法。 

* SD方法

[SD（最陡下降）方法](https://en.wikipedia.org/wiki/Gradient_descent) 是最简单的一阶优化方法之一，在每一步中，运动沿着梯度方向进行，函数在该方向上下降最快。在实践中，SD方法可能需要很多次迭代才能收敛，通常不使用。 

* CG方法

[CG（共轭梯度）方法](https://en.wikipedia.org/wiki/Conjugate_gradient_method) 是求解优化问题的最广泛使用的方法之一。在ABACUS中，我们实现了用于固定单元格结构弛豫以及单元格参数优化的CG方法。

## 带约束的结构优化
除了允许所有自由度移动的传统优化之外，我们还在ABACUS中提供了受限优化的选项。

### 固定原子位置
用户可能注意到，在上述示例中，STRU文件中的原子位置是与三个整数一起给出的：

```
Mg                      # Name of element        
0.0                     # Magnetic for this element.
4                       # Number of atoms
0.0  0.0  0.0  0 0 0    # x,y,z, move_x, move_y, move_z
0.0  0.5  0.5  0 0 0    # x,y,z, move_x, move_y, move_z
0.5  0.0  0.5  0 0 0    # x,y,z, move_x, move_y, move_z
0.5  0.5  0.0  0 0 0    # x,y,z, move_x, move_y, move_z
```

对于弛豫计算，这三个整数表示是否允许相应的自由度移动。例如，如果我们用以下STRU文件替换：

```
Mg                      # Name of element        
0.0                     # Magnetic for this element.
4                       # Number of atoms
0.0  0.0  0.0  0 0 1    # x,y,z, move_x, move_y, move_z
0.0  0.5  0.5  0 0 0    # x,y,z, move_x, move_y, move_z
0.5  0.0  0.5  0 0 0    # x,y,z, move_x, move_y, move_z
0.5  0.5  0.0  0 0 0    # x,y,z, move_x, move_y, move_z
```

那么第一个 Mg 原子将允许在 z 方向上移动。
在分子/团簇的弛豫过程中，固定原子位置有时有助于防止系统在空间中漂移。

### 固定晶胞参数
有时我们希望在固定一些晶胞自由度的情况下进行可变晶胞弛豫。这可以通过诸如[fixed_axes](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#fixed-axes)、[fixed_ibrav](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#fixed-ibrav)和[fixed_atoms](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#fixed-atoms)等关键字实现。

具体来说，如果用户熟悉VASP中的ISIF选项，那么我们提供如下对应关系：
* ISIF = 0 : calculation = “relax”
* ISIF = 1, 2 : calculation = “relax”, cal_stress = 1
* ISIF = 3 : calculation = “cell-relax”
* ISIF = 4 : calculation = “cell-relax”, fixed_axes = “volume”
* ISIF = 5 : calculation = “cell-relax”, fixed_axes = “volume”, fixed_atoms = True
* ISIF = 6 : calculation = “cell-relax”, fixed_atoms = True
* ISIF = 7 : calculation = “cell-realx”, fixed_axes = “shape”, fixed_atoms = True

[ISIF 参数表](https://www.vasp.at/wiki/index.php/ISIF)
| ISIF | calculate |   calculate   | degrees-of-freedom | degrees-of-freedom | degrees-of-freedom |
| :--: | :-------: | :-----------: | :----------------: | :----------------: | :----------------: |
|      |  forces   | Stress tensor | positions          |     cell shape     |    cell volume     |
|  0   |    yes    |      no       | yes                |         no         |         no         |
|  1   |    yes    |  trace only   | yes                |         no         |         no         |
|  2   |    yes    |      yes      | yes                |         no         |         no         |
|  3   |    yes    |      yes      | yes                |        yes         |        yes         |
|  4   |    yes    |      yes      | yes                |        yes         |         no         |
|  5   |    yes    |      yes      | no                 |        yes         |         no         |
|  6   |    yes    |      yes      | no                 |        yes         |        yes         |
|  7   |    yes    |      yes      | no                 |         no         |        yes         |

## 【练习】$H_2$ 分子键长优化（LCAO 基组）

在这一小节中，我们将练习如何优化氢气（$H_2$）分子键长。键长是分子中两个原子之间的平均距离，对于氢分子来说，它是两个氢原子之间的距离。在实验室中，可以通过各种光谱学方法测量键长，而在理论计算领域，我们可以通过量子化学方法精确计算它。在本教程中，我们将展示如何使用量子化学方法计算并优化氢气分子的键长，以获得最低能量状态。这对于了解化学反应和分子性质具有重要意义。

$H_2$ 分子的键长在实验中测量为大约 0.74 Å（即 0.74 x 10⁻¹⁰米，或 74 pm）。

### STRU 文件

我们先建立一个 STRU 文件以模拟氢原子。在这个例子中，晶格是一个 15x15x15 的立方体。这样设置的目的是为了避免周期性边界条件下相邻原子的相互作用。对于 H2 分子，它是一个非常小的分子，其实际尺寸远小于 15 Å。通过设置一个较大的晶格，可以确保分子处于足够大的空间中，从而避免模拟中的周期性影响。这是一种常见的方法，用于处理孤立分子或原子的计算。


```python
#  新建练习文件夹
import os

practice_path = os.path.join("ABACUS_Relax", "Practice")
if not os.path.exists(practice_path):
    os.mkdir(practice_path)

# 新建一个 STRU 文件
stru = """# H2 Relax STRU
ATOMIC_SPECIES
H 1.008  H_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
H_gga_6au_100Ry_2s1p.orb

LATTICE_CONSTANT
1.8897259886 				# 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
10.0000 0.00000 0.00000  
0.00000 10.0000 0.00000
0.00000 0.00000 10.0000

ATOMIC_POSITIONS
Cartesian                   # Cartesian(Unit is LATTICE_CONSTANT)
H                       	# Name of element        
1.0                         # Magnetic for this element.
2                       	# Number of atoms
0.00  0.00  0.00  0 0 1 	# x,y,z, move_x, move_y, move_z
0.00  0.00  0.74  0 0 1 	# x,y,z, move_x, move_y, move_z
"""
stru_path = os.path.join(practice_path, "STRU")
with open(stru_path, "w") as f:
    f.write(stru)
```

接下来，准备 H2 分子结构优化的 INPUT 文件与 K 点文件：(出于演示目的，这里仅使用较低设置的参数进行练习)


```python
# 新建一个 INPUT 文件
input_content = """INPUT_PARAMETERS
suffix                  H2                  # Suffix for output files
pseudo_dir              ./PP_ORB            # Directory path for pseudopotential files
orbital_dir             ./PP_ORB            # Directory path for orbital files
ecutwfc                 100                 # Plane-wave energy cutoff in Rydberg
scf_thr                 1e-7                # Self-consistency threshold in Rydberg
basis_type              lcao                # Basis set type (e.g. lcao for linear combination of atomic orbitals)
calculation             relax               # Type of calculation (e.g. relax for geometry optimization)
force_thr_ev            0.1                 # Force convergence threshold in eV/Angstrom
stress_thr              2.0                 # Stress convergence threshold in kBar
relax_nmax              100                 # Max ionic iteration steps for relaxation calculation
out_stru                True                # Output final optimized structure if set to True
cal_stress              True                # Calculate stress tensor if set to True
nspin                   2                   # Number of spin channels (e.g. 2 for spin-polarized calculations)
"""

input_path = os.path.join(practice_path, "INPUT")
with open(input_path, "w") as f:
    f.write(input_content)
    
# 新建一个 KPT 文件
kpt_content = """K_POINTS
0 
Gamma
1 1 1 0 0 0
"""

kpt_path = os.path.join(practice_path, "KPT")
with open(kpt_path, "w") as f:
    f.write(kpt_content)
```


```bash
%%bash
# 运行氢气分子结构优化
cd ./ABACUS_Relax/Practice
OMP_NUM_THREADS=2 mpirun -n 2 abacus
```

                                                                                         
                                  ABACUS v3.3.2
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: e39b50efe (Fri Aug 18 16:14:25 2023 +0800)
    
     Thu Feb 29 12:22:48 2024
     MAKE THE DIR         : OUT.H2/
     UNIFORM GRID DIM     : 125 * 125 * 125
     UNIFORM GRID DIM(BIG): 25 * 25 * 25
     DONE(0.289009   SEC) : SETUP UNITCELL
     DONE(0.334426   SEC) : SYMMETRY
     DONE(0.438936   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     2       2               2           10          
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     H       2s1p-6au        5           2           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.70998    SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -3.24e-05 8.26e-01  6.925057e-01   0.000000e+00   1.002e+00  2.196e+00  
     GE2    -6.08e-06 2.49e-01  -3.348568e+01  -3.417819e+01  4.101e-01  1.963e+00  
     GE3    -2.02e-06 7.37e-02  -3.333818e+01  1.475023e-01   2.345e-01  1.961e+00  
     GE4    2.69e-07  1.11e-02  -3.063830e+01  2.699877e+00   7.943e-02  1.963e+00  
     GE5    -1.00e-08 4.90e-04  -3.152012e+01  -8.818140e-01  1.176e-02  1.997e+00  
     GE6    -4.25e-08 1.54e-03  -3.169529e+01  -1.751700e-01  4.847e-03  1.950e+00  
     GE7    8.53e-09  3.41e-04  -3.162243e+01  7.286222e-02   1.095e-03  1.993e+00  
     GE8    9.32e-10  3.71e-05  -3.163447e+01  -1.204232e-02  1.193e-04  1.973e+00  
     GE9    -1.37e-11 2.40e-07  -3.163598e+01  -1.514257e-03  1.514e-06  1.956e+00  
     GE10   2.90e-13  1.51e-08  -3.163599e+01  -3.326612e-06  5.458e-08  1.807e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -4.274e+01     0.000e+00      -7.083e-32     
     0.000e+00      -4.274e+01     0.000e+00      
     -7.083e-32     0.000e+00      -1.761e+01     
     TOTAL-PRESSURE: -3.437e+01 KBAR
     ETOT DIFF (eV)       : 0.000e+00
     LARGEST GRAD (eV/A)  : 1.659e+01
     -------------------------------------------
     STEP OF RELAXATION : 2
     -------------------------------------------
     DONE(3.011e+01  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -3.77e-07 1.96e-01  -3.482052e+01  0.000000e+00   1.104e+00  1.849e+00  
     GE2    -9.84e-08 7.99e-02  -3.678444e+01  -1.963921e+00  3.661e-01  1.955e+00  
     GE3    -1.35e-08 9.52e-03  -3.687538e+01  -9.093861e-02  2.883e-02  1.990e+00  
     GE4    -3.19e-09 2.03e-03  -3.685591e+01  1.946647e-02   8.231e-03  1.977e+00  
     GE5    3.42e-11  2.74e-05  -3.685408e+01  1.831917e-03   1.908e-04  1.975e+00  
     GE6    7.44e-12  1.16e-05  -3.685498e+01  -9.054716e-04  5.038e-05  2.029e+00  
     GE7    -4.92e-13 8.15e-07  -3.685533e+01  -3.487319e-04  1.664e-05  1.973e+00  
     GE8    1.60e-15  2.36e-07  -3.685522e+01  1.105357e-04   6.553e-07  1.967e+00  
     GE9    2.88e-15  2.17e-09  -3.685522e+01  -1.882580e-06  2.553e-08  1.807e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.340e+01     0.000e+00      -6.419e-32     
     0.000e+00      -3.340e+01     0.000e+00      
     -6.419e-32     0.000e+00      -3.132e+01     
     TOTAL-PRESSURE: -3.271e+01 KBAR
     ETOT DIFF (eV)       : -5.219e+00
     LARGEST GRAD (eV/A)  : 3.606e-01
     -------------------------------------------
     STEP OF RELAXATION : 3
     -------------------------------------------
     DONE(5.539e+01  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -2.85e-10 1.04e-02  -3.686309e+01  0.000000e+00   2.746e-02  1.838e+00  
     GE2    -5.98e-11 4.35e-03  -3.687071e+01  -7.618362e-03  1.123e-02  1.937e+00  
     GE3    3.97e-11  9.48e-05  -3.687228e+01  -1.565757e-03  4.632e-04  1.956e+00  
     GE4    7.79e-12  5.09e-06  -3.687190e+01  3.785484e-04   1.920e-05  1.950e+00  
     GE5    2.32e-12  1.31e-06  -3.687192e+01  -2.146043e-05  3.874e-06  1.938e+00  
     GE6    2.47e-14  1.50e-08  -3.687193e+01  -1.007637e-05  3.094e-07  1.967e+00  
     GE7    -2.18e-15 5.79e-10  -3.687193e+01  2.620974e-06   2.547e-08  1.732e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.337e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.337e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.145e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.670e-02
     LARGEST GRAD (eV/A)  : 3.300e-01
     -------------------------------------------
     STEP OF RELAXATION : 4
     -------------------------------------------
     DONE(7.641e+01  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    9.49e-13  1.85e-04  -3.687212e+01  0.000000e+00   4.910e-04  1.892e+00  
     GE2    1.37e-12  7.79e-05  -3.687219e+01  -7.240341e-05  2.009e-04  1.974e+00  
     GE3    7.07e-13  1.63e-06  -3.687222e+01  -2.720243e-05  8.240e-06  1.996e+00  
     GE4    1.13e-13  6.45e-08  -3.687221e+01  7.242570e-06   1.884e-07  1.975e+00  
     GE5    1.46e-14  8.27e-09  -3.687221e+01  -8.042551e-07  2.477e-08  1.799e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.337e+01     0.000e+00      -2.213e-33     
     0.000e+00      -3.337e+01     0.000e+00      
     -2.213e-33     0.000e+00      -3.145e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -2.874e-04
     LARGEST GRAD (eV/A)  : 3.295e-01
     -------------------------------------------
     STEP OF RELAXATION : 5
     -------------------------------------------
     DONE(9.361e+01  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    6.94e-12  5.56e-04  -3.687279e+01  0.000000e+00   1.473e-03  1.839e+00  
     GE2    5.84e-12  2.34e-04  -3.687301e+01  -2.236015e-04  6.025e-04  1.991e+00  
     GE3    2.17e-12  4.90e-06  -3.687309e+01  -8.152375e-05  2.472e-05  2.137e+00  
     GE4    3.42e-13  1.94e-07  -3.687307e+01  2.163208e-05   5.698e-07  2.101e+00  
     GE5    5.33e-14  2.97e-08  -3.687307e+01  -2.312848e-06  8.874e-08  1.777e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -5.312e-32     
     0.000e+00      -3.336e+01     0.000e+00      
     -5.312e-32     0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -8.590e-04
     LARGEST GRAD (eV/A)  : 3.279e-01
     -------------------------------------------
     STEP OF RELAXATION : 6
     -------------------------------------------
     DONE(1.113e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    9.25e-14  5.94e-06  -3.687308e+01  0.000000e+00   1.573e-05  1.875e+00  
     GE2    7.01e-14  2.50e-06  -3.687308e+01  -2.288324e-06  6.435e-06  1.974e+00  
     GE3    2.32e-14  5.22e-08  -3.687308e+01  -8.726475e-07  2.641e-07  2.002e+00  
     GE4    4.97e-15  2.04e-09  -3.687308e+01  2.327435e-07   5.960e-09  1.778e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -8.987e-06
     LARGEST GRAD (eV/A)  : 3.279e-01
     -------------------------------------------
     STEP OF RELAXATION : 7
     -------------------------------------------
     DONE(1.265e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.58e-13  1.78e-05  -3.687310e+01  0.000000e+00   4.719e-05  1.972e+00  
     GE2    2.01e-13  7.49e-06  -3.687311e+01  -6.845099e-06  1.931e-05  1.982e+00  
     GE3    7.01e-14  1.57e-07  -3.687311e+01  -2.614423e-06  7.922e-07  2.013e+00  
     GE4    1.04e-14  6.15e-09  -3.687311e+01  6.935118e-07   1.794e-08  1.818e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -2.741e-05
     LARGEST GRAD (eV/A)  : 3.279e-01
     -------------------------------------------
     STEP OF RELAXATION : 8
     -------------------------------------------
     DONE(1.419e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.41e-13  9.50e-06  -3.687312e+01  0.000000e+00   2.516e-05  1.859e+00  
     GE2    1.09e-13  3.99e-06  -3.687312e+01  -3.638551e-06  1.029e-05  1.986e+00  
     GE3    3.67e-14  8.35e-08  -3.687312e+01  -1.386847e-06  4.224e-07  1.967e+00  
     GE4    6.03e-15  3.28e-09  -3.687312e+01  3.699939e-07   9.554e-09  1.773e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.336e+01     0.000e+00      
     8.854e-33      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.466e-05
     LARGEST GRAD (eV/A)  : 3.278e-01
     -------------------------------------------
     STEP OF RELAXATION : 9
     -------------------------------------------
     DONE(1.574e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    4.25e-13  2.85e-05  -3.687315e+01  0.000000e+00   7.548e-05  1.857e+00  
     GE2    3.28e-13  1.20e-05  -3.687316e+01  -1.095211e-05  3.088e-05  1.982e+00  
     GE3    1.10e-13  2.50e-07  -3.687317e+01  -4.167716e-06  1.267e-06  1.972e+00  
     GE4    1.66e-14  9.83e-09  -3.687317e+01  1.111410e-06   2.868e-08  1.766e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.336e+01     0.000e+00      
     8.854e-33      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -4.385e-05
     LARGEST GRAD (eV/A)  : 3.278e-01
     -------------------------------------------
     STEP OF RELAXATION : 10
     -------------------------------------------
     DONE(1.729e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.32e-13  1.52e-05  -3.687318e+01  0.000000e+00   4.024e-05  1.853e+00  
     GE2    1.77e-13  6.39e-06  -3.687319e+01  -5.822519e-06  1.646e-05  1.965e+00  
     GE3    5.85e-14  1.34e-07  -3.687319e+01  -2.218398e-06  6.755e-07  1.972e+00  
     GE4    9.67e-15  5.24e-09  -3.687319e+01  5.917352e-07   1.528e-08  1.794e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      4.427e-33      
     0.000e+00      -3.336e+01     0.000e+00      
     4.427e-33      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -2.345e-05
     LARGEST GRAD (eV/A)  : 3.277e-01
     -------------------------------------------
     STEP OF RELAXATION : 11
     -------------------------------------------
     DONE(1.884e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    7.06e-13  4.56e-05  -3.687324e+01  0.000000e+00   1.207e-04  1.867e+00  
     GE2    5.36e-13  1.92e-05  -3.687326e+01  -1.753182e-05  4.939e-05  2.578e+00  
     GE3    1.79e-13  4.01e-07  -3.687326e+01  -6.663448e-06  2.027e-06  1.972e+00  
     GE4    2.82e-14  1.57e-08  -3.687326e+01  1.776526e-06   4.589e-08  1.808e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      3.099e-32      
     0.000e+00      -3.336e+01     0.000e+00      
     3.099e-32      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -7.010e-05
     LARGEST GRAD (eV/A)  : 3.276e-01
     -------------------------------------------
     STEP OF RELAXATION : 12
     -------------------------------------------
     DONE(2.043e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    3.87e-13  2.43e-05  -3.687329e+01  0.000000e+00   6.433e-05  1.875e+00  
     GE2    2.90e-13  1.02e-05  -3.687330e+01  -9.313979e-06  2.632e-05  1.971e+00  
     GE3    9.46e-14  2.13e-07  -3.687330e+01  -3.546279e-06  1.080e-06  1.968e+00  
     GE4    1.41e-14  8.37e-09  -3.687330e+01  9.457588e-07   2.442e-08  1.773e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     -3.541e-32     
     0.000e+00      -3.541e-32     -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -3.750e-05
     LARGEST GRAD (eV/A)  : 3.275e-01
     -------------------------------------------
     STEP OF RELAXATION : 13
     -------------------------------------------
     DONE(2.195e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.19e-12  7.29e-05  -3.687338e+01  0.000000e+00   1.930e-04  1.868e+00  
     GE2    8.83e-13  3.06e-05  -3.687340e+01  -2.807405e-05  7.896e-05  1.976e+00  
     GE3    2.86e-13  6.41e-07  -3.687341e+01  -1.064443e-05  3.240e-06  1.961e+00  
     GE4    4.37e-14  2.51e-08  -3.687341e+01  2.838254e-06   7.338e-08  1.827e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.120e-04
     LARGEST GRAD (eV/A)  : 3.273e-01
     -------------------------------------------
     STEP OF RELAXATION : 14
     -------------------------------------------
     DONE(2.346e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    6.67e-13  3.88e-05  -3.687345e+01  0.000000e+00   1.028e-04  1.852e+00  
     GE2    4.82e-13  1.63e-05  -3.687347e+01  -1.488979e-05  4.206e-05  1.969e+00  
     GE3    1.53e-13  3.41e-07  -3.687347e+01  -5.665483e-06  1.726e-06  1.944e+00  
     GE4    2.27e-14  1.34e-08  -3.687347e+01  1.509253e-06   3.903e-08  1.793e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     1.771e-32      
     0.000e+00      1.771e-32      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -5.988e-05
     LARGEST GRAD (eV/A)  : 3.272e-01
     -------------------------------------------
     STEP OF RELAXATION : 15
     -------------------------------------------
     DONE(2.498e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.18e-12  1.16e-04  -3.687359e+01  0.000000e+00   3.084e-04  1.873e+00  
     GE2    1.52e-12  4.90e-05  -3.687364e+01  -4.500853e-05  1.262e-04  2.009e+00  
     GE3    4.59e-13  1.02e-06  -3.687365e+01  -1.699971e-05  5.178e-06  1.971e+00  
     GE4    7.32e-14  4.02e-08  -3.687365e+01  4.532034e-06   1.173e-07  1.992e+00  
     GE5    8.75e-15  4.95e-09  -3.687365e+01  -5.076063e-07  1.486e-08  1.783e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      2.656e-32      
     0.000e+00      -3.336e+01     0.000e+00      
     2.656e-32      0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.794e-04
     LARGEST GRAD (eV/A)  : 3.269e-01
     -------------------------------------------
     STEP OF RELAXATION : 16
     -------------------------------------------
     DONE(2.673e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.21e-12  6.20e-05  -3.687371e+01  0.000000e+00   1.642e-04  1.841e+00  
     GE2    8.33e-13  2.61e-05  -3.687374e+01  -2.385476e-05  6.717e-05  1.973e+00  
     GE3    2.44e-13  5.45e-07  -3.687375e+01  -9.051603e-06  2.757e-06  1.979e+00  
     GE4    3.72e-14  2.14e-08  -3.687375e+01  2.414179e-06   6.230e-08  1.761e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -6.640e-33     
     0.000e+00      -3.336e+01     0.000e+00      
     -6.640e-33     0.000e+00      -3.146e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -9.501e-05
     LARGEST GRAD (eV/A)  : 3.267e-01
     -------------------------------------------
     STEP OF RELAXATION : 17
     -------------------------------------------
     DONE(2.824e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    3.87e-12  1.86e-04  -3.687394e+01  0.000000e+00   4.925e-04  1.840e+00  
     GE2    2.60e-12  7.82e-05  -3.687401e+01  -7.225212e-05  2.015e-04  1.960e+00  
     GE3    7.36e-13  1.64e-06  -3.687404e+01  -2.712801e-05  8.271e-06  1.976e+00  
     GE4    1.13e-13  6.41e-08  -3.687403e+01  7.228285e-06   1.873e-07  1.952e+00  
     GE5    1.53e-14  8.23e-09  -3.687403e+01  -8.037960e-07  2.467e-08  1.744e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      1.771e-32      
     0.000e+00      -3.336e+01     0.000e+00      
     1.771e-32      0.000e+00      -3.147e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -2.859e-04
     LARGEST GRAD (eV/A)  : 3.262e-01
     -------------------------------------------
     STEP OF RELAXATION : 18
     -------------------------------------------
     DONE(2.996e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.23e-12  9.89e-05  -3.687413e+01  0.000000e+00   2.619e-04  1.821e+00  
     GE2    1.46e-12  4.16e-05  -3.687417e+01  -3.811212e-05  1.072e-04  1.955e+00  
     GE3    3.93e-13  8.69e-07  -3.687419e+01  -1.441511e-05  4.398e-06  1.929e+00  
     GE4    5.97e-14  3.40e-08  -3.687418e+01  3.844212e-06   9.928e-08  1.740e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -1.771e-32     
     0.000e+00      -3.336e+01     0.000e+00      
     -1.771e-32     0.000e+00      -3.147e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.513e-04
     LARGEST GRAD (eV/A)  : 3.259e-01
     -------------------------------------------
     STEP OF RELAXATION : 19
     -------------------------------------------
     DONE(3.146e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    7.18e-12  2.97e-04  -3.687449e+01  0.000000e+00   7.857e-04  1.851e+00  
     GE2    4.58e-12  1.25e-04  -3.687460e+01  -1.161166e-04  3.214e-04  2.450e+00  
     GE3    1.18e-12  2.61e-06  -3.687465e+01  -4.319740e-05  1.320e-05  1.960e+00  
     GE4    1.82e-13  1.02e-07  -3.687464e+01  1.150487e-05   2.991e-07  1.966e+00  
     GE5    2.49e-14  1.40e-08  -3.687464e+01  -1.261983e-06  4.179e-08  1.754e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.336e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.147e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -4.550e-04
     LARGEST GRAD (eV/A)  : 3.251e-01
     -------------------------------------------
     STEP OF RELAXATION : 20
     -------------------------------------------
     DONE(3.325e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    4.12e-12  1.58e-04  -3.687480e+01  0.000000e+00   4.171e-04  1.846e+00  
     GE2    2.56e-12  6.62e-05  -3.687486e+01  -6.090269e-05  1.706e-04  1.945e+00  
     GE3    6.29e-13  1.38e-06  -3.687488e+01  -2.291766e-05  7.005e-06  1.971e+00  
     GE4    9.66e-14  5.40e-08  -3.687488e+01  6.109882e-06   1.579e-07  1.946e+00  
     GE5    1.24e-14  6.86e-09  -3.687488e+01  -6.813236e-07  2.055e-08  1.734e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -8.854e-33     
     0.000e+00      -3.336e+01     0.000e+00      
     -8.854e-33     0.000e+00      -3.147e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -2.408e-04
     LARGEST GRAD (eV/A)  : 3.246e-01
     -------------------------------------------
     STEP OF RELAXATION : 21
     -------------------------------------------
     DONE(3.494e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.37e-11  4.73e-04  -3.687536e+01  0.000000e+00   1.251e-03  1.833e+00  
     GE2    8.22e-12  1.99e-04  -3.687555e+01  -1.871920e-04  5.118e-04  1.962e+00  
     GE3    1.90e-12  4.16e-06  -3.687562e+01  -6.863467e-05  2.102e-05  1.952e+00  
     GE4    2.89e-13  1.63e-07  -3.687560e+01  1.826736e-05   4.767e-07  1.944e+00  
     GE5    4.38e-14  2.42e-08  -3.687560e+01  -1.967347e-06  7.233e-08  1.743e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -1.771e-32     
     0.000e+00      -3.336e+01     7.083e-32      
     -1.771e-32     7.083e-32      -3.148e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -7.209e-04
     LARGEST GRAD (eV/A)  : 3.233e-01
     -------------------------------------------
     STEP OF RELAXATION : 22
     -------------------------------------------
     DONE(3.666e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    7.92e-12  2.50e-04  -3.687585e+01  0.000000e+00   6.622e-04  1.822e+00  
     GE2    4.64e-12  1.05e-04  -3.687595e+01  -9.714434e-05  2.709e-04  1.941e+00  
     GE3    1.01e-12  2.20e-06  -3.687599e+01  -3.626887e-05  1.113e-05  1.950e+00  
     GE4    1.52e-13  8.55e-08  -3.687598e+01  9.660643e-06   2.500e-07  2.031e+00  
     GE5    2.03e-14  1.14e-08  -3.687598e+01  -1.070362e-06  3.420e-08  1.791e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.336e+01     0.000e+00      
     8.854e-33      0.000e+00      -3.148e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -3.806e-04
     LARGEST GRAD (eV/A)  : 3.226e-01
     -------------------------------------------
     STEP OF RELAXATION : 23
     -------------------------------------------
     DONE(3.840e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.90e-11  7.50e-04  -3.687673e+01  0.000000e+00   1.986e-03  1.864e+00  
     GE2    1.62e-11  3.16e-04  -3.687703e+01  -3.031927e-04  8.126e-04  1.948e+00  
     GE3    3.09e-12  6.60e-06  -3.687714e+01  -1.086491e-04  3.338e-05  1.948e+00  
     GE4    4.58e-13  2.58e-07  -3.687711e+01  2.884301e-05   7.589e-07  1.948e+00  
     GE5    7.61e-14  4.28e-08  -3.687712e+01  -3.026725e-06  1.276e-07  1.942e+00  
     GE6    7.93e-16  3.75e-10  -3.687712e+01  6.016736e-08   1.073e-08  1.736e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      -7.083e-32     
     0.000e+00      -3.336e+01     0.000e+00      
     -7.083e-32     0.000e+00      -3.149e+01     
     TOTAL-PRESSURE: -3.273e+01 KBAR
     ETOT DIFF (eV)       : -1.137e-03
     LARGEST GRAD (eV/A)  : 3.205e-01
     -------------------------------------------
     STEP OF RELAXATION : 24
     -------------------------------------------
     DONE(4.033e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.78e-11  3.96e-04  -3.687751e+01  0.000000e+00   1.047e-03  1.856e+00  
     GE2    9.57e-12  1.66e-04  -3.687767e+01  -1.549441e-04  4.282e-04  1.952e+00  
     GE3    1.65e-12  3.48e-06  -3.687773e+01  -5.713872e-05  1.759e-05  1.944e+00  
     GE4    2.41e-13  1.35e-07  -3.687771e+01  1.522105e-05   3.938e-07  1.929e+00  
     GE5    3.61e-14  1.94e-08  -3.687771e+01  -1.655344e-06  5.799e-08  1.731e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.336e+01     0.000e+00      3.541e-32      
     0.000e+00      -3.336e+01     0.000e+00      
     3.541e-32      0.000e+00      -3.149e+01     
     TOTAL-PRESSURE: -3.274e+01 KBAR
     ETOT DIFF (eV)       : -5.964e-04
     LARGEST GRAD (eV/A)  : 3.194e-01
     -------------------------------------------
     STEP OF RELAXATION : 25
     -------------------------------------------
     DONE(4.206e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    5.94e-11  1.19e-03  -3.687887e+01  0.000000e+00   3.140e-03  1.832e+00  
     GE2    3.13e-11  4.99e-04  -3.687936e+01  -4.940370e-04  1.284e-03  1.944e+00  
     GE3    5.05e-12  1.04e-05  -3.687953e+01  -1.708631e-04  5.279e-05  1.934e+00  
     GE4    7.32e-13  4.07e-07  -3.687949e+01  4.521380e-05   1.206e-06  2.062e+00  
     GE5    1.41e-13  7.75e-08  -3.687949e+01  -4.553539e-06  2.302e-07  1.928e+00  
     GE6    1.73e-15  7.16e-10  -3.687949e+01  -1.065360e-07  2.116e-08  1.771e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.335e+01     0.000e+00      1.771e-32      
     0.000e+00      -3.335e+01     0.000e+00      
     1.771e-32      0.000e+00      -3.151e+01     
     TOTAL-PRESSURE: -3.274e+01 KBAR
     ETOT DIFF (eV)       : -1.779e-03
     LARGEST GRAD (eV/A)  : 3.161e-01
     -------------------------------------------
     STEP OF RELAXATION : 26
     -------------------------------------------
     DONE(4.396e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    3.35e-11  6.21e-04  -3.688010e+01  0.000000e+00   1.643e-03  1.825e+00  
     GE2    1.74e-11  2.61e-04  -3.688035e+01  -2.458226e-04  6.722e-04  1.947e+00  
     GE3    2.65e-12  5.46e-06  -3.688044e+01  -8.891916e-05  2.763e-05  1.964e+00  
     GE4    3.77e-13  2.09e-07  -3.688041e+01  2.366646e-05   6.146e-07  1.942e+00  
     GE5    6.12e-14  3.35e-08  -3.688042e+01  -2.514660e-06  9.979e-08  1.787e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.335e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.335e+01     0.000e+00      
     8.854e-33      0.000e+00      -3.152e+01     
     TOTAL-PRESSURE: -3.274e+01 KBAR
     ETOT DIFF (eV)       : -9.249e-04
     LARGEST GRAD (eV/A)  : 3.144e-01
     -------------------------------------------
     STEP OF RELAXATION : 27
     -------------------------------------------
     DONE(4.569e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.31e-10  1.86e-03  -3.688215e+01  0.000000e+00   4.928e-03  1.920e+00  
     GE2    6.52e-11  7.84e-04  -3.688296e+01  -8.099909e-04  2.016e-03  1.977e+00  
     GE3    8.25e-12  1.64e-05  -3.688323e+01  -2.653723e-04  8.293e-05  1.959e+00  
     GE4    1.15e-12  6.38e-07  -3.688316e+01  6.989851e-05   1.914e-06  1.980e+00  
     GE5    2.55e-13  1.40e-07  -3.688316e+01  -6.656164e-06  4.157e-07  1.985e+00  
     GE6    3.23e-15  1.32e-09  -3.688317e+01  -5.401389e-07  3.920e-08  1.784e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.335e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.335e+01     -3.541e-32     
     8.854e-33      -3.541e-32     -3.154e+01     
     TOTAL-PRESSURE: -3.274e+01 KBAR
     ETOT DIFF (eV)       : -2.749e-03
     LARGEST GRAD (eV/A)  : 3.094e-01
     -------------------------------------------
     STEP OF RELAXATION : 28
     -------------------------------------------
     DONE(4.766e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    7.69e-11  9.66e-04  -3.688408e+01  0.000000e+00   2.552e-03  1.845e+00  
     GE2    3.75e-11  4.06e-04  -3.688447e+01  -3.887082e-04  1.044e-03  1.939e+00  
     GE3    4.29e-12  8.48e-06  -3.688461e+01  -1.365226e-04  4.295e-05  1.989e+00  
     GE4    5.81e-13  3.21e-07  -3.688457e+01  3.628014e-05   9.465e-07  1.976e+00  
     GE5    1.05e-13  5.84e-08  -3.688457e+01  -3.727222e-06  1.733e-07  1.955e+00  
     GE6    9.40e-16  4.96e-10  -3.688457e+01  -2.223662e-08  1.548e-08  1.748e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.334e+01     0.000e+00      -5.533e-33     
     0.000e+00      -3.334e+01     0.000e+00      
     -5.533e-33     0.000e+00      -3.155e+01     
     TOTAL-PRESSURE: -3.274e+01 KBAR
     ETOT DIFF (eV)       : -1.407e-03
     LARGEST GRAD (eV/A)  : 3.068e-01
     -------------------------------------------
     STEP OF RELAXATION : 29
     -------------------------------------------
     DONE(4.964e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.18e-10  2.90e-03  -3.688708e+01  0.000000e+00   7.652e-03  1.885e+00  
     GE2    1.07e-10  1.22e-03  -3.688842e+01  -1.340740e-03  3.132e-03  1.945e+00  
     GE3    1.28e-11  2.56e-05  -3.688883e+01  -4.061253e-04  1.289e-04  1.949e+00  
     GE4    1.78e-12  9.91e-07  -3.688872e+01  1.061390e-04   3.052e-06  1.952e+00  
     GE5    4.53e-13  2.49e-07  -3.688873e+01  -9.388299e-06  7.334e-07  2.328e+00  
     GE6    3.71e-15  2.26e-09  -3.688873e+01  -1.439993e-06  6.808e-08  1.731e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.334e+01     0.000e+00      8.854e-33      
     0.000e+00      -3.334e+01     0.000e+00      
     8.854e-33      0.000e+00      -3.158e+01     
     TOTAL-PRESSURE: -3.275e+01 KBAR
     ETOT DIFF (eV)       : -4.161e-03
     LARGEST GRAD (eV/A)  : 2.991e-01
     -------------------------------------------
     STEP OF RELAXATION : 30
     -------------------------------------------
     DONE(5.157e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.22e-10  1.48e-03  -3.689005e+01  0.000000e+00   3.899e-03  1.857e+00  
     GE2    5.89e-11  6.21e-04  -3.689066e+01  -6.085491e-04  1.596e-03  1.951e+00  
     GE3    6.44e-12  1.30e-05  -3.689087e+01  -2.045479e-04  6.568e-05  1.943e+00  
     GE4    8.65e-13  4.80e-07  -3.689081e+01  5.420652e-05   1.430e-06  1.953e+00  
     GE5    1.83e-13  1.01e-07  -3.689082e+01  -5.323244e-06  2.978e-07  1.935e+00  
     GE6    2.26e-15  8.77e-10  -3.689082e+01  -2.981132e-07  2.803e-08  1.799e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.333e+01     0.000e+00      1.771e-32      
     0.000e+00      -3.333e+01     0.000e+00      
     1.771e-32      0.000e+00      -3.160e+01     
     TOTAL-PRESSURE: -3.275e+01 KBAR
     ETOT DIFF (eV)       : -2.085e-03
     LARGEST GRAD (eV/A)  : 2.952e-01
     -------------------------------------------
     STEP OF RELAXATION : 31
     -------------------------------------------
     DONE(5.347e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.29e-10  4.44e-03  -3.689423e+01  0.000000e+00   1.170e-02  1.840e+00  
     GE2    1.19e-10  1.87e-03  -3.689646e+01  -2.234379e-03  4.789e-03  1.955e+00  
     GE3    1.77e-11  3.92e-05  -3.689707e+01  -6.062386e-04  1.973e-04  1.928e+00  
     GE4    2.62e-12  1.52e-06  -3.689691e+01  1.564350e-04   4.939e-06  1.962e+00  
     GE5    7.42e-13  4.19e-07  -3.689693e+01  -1.266252e-05  1.230e-06  2.012e+00  
     GE6    6.85e-15  3.63e-09  -3.689693e+01  -2.993824e-06  1.107e-07  2.100e+00  
     GE7    -3.12e-16 1.87e-10  -3.689693e+01  9.812776e-07   1.085e-08  1.745e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.332e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.332e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.165e+01     
     TOTAL-PRESSURE: -3.276e+01 KBAR
     ETOT DIFF (eV)       : -6.111e-03
     LARGEST GRAD (eV/A)  : 2.839e-01
     -------------------------------------------
     STEP OF RELAXATION : 32
     -------------------------------------------
     DONE(5.557e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.10e-10  2.21e-03  -3.689873e+01  0.000000e+00   5.825e-03  1.855e+00  
     GE2    5.67e-11  9.30e-04  -3.689967e+01  -9.385306e-04  2.385e-03  1.981e+00  
     GE3    8.31e-12  1.93e-05  -3.689996e+01  -2.964169e-04  9.812e-05  1.930e+00  
     GE4    1.20e-12  6.95e-07  -3.689988e+01  7.825694e-05   2.105e-06  1.937e+00  
     GE5    2.93e-13  1.68e-07  -3.689989e+01  -7.258368e-06  4.923e-07  1.946e+00  
     GE6    2.53e-15  1.38e-09  -3.689989e+01  -8.572351e-07  4.673e-08  1.766e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.331e+01     0.000e+00      7.083e-32      
     0.000e+00      -3.331e+01     0.000e+00      
     7.083e-32      0.000e+00      -3.167e+01     
     TOTAL-PRESSURE: -3.277e+01 KBAR
     ETOT DIFF (eV)       : -2.964e-03
     LARGEST GRAD (eV/A)  : 2.785e-01
     -------------------------------------------
     STEP OF RELAXATION : 33
     -------------------------------------------
     DONE(5.751e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.62e-10 6.64e-03  -3.690409e+01  0.000000e+00   1.749e-02  1.844e+00  
     GE2    -4.08e-11 2.80e-03  -3.690781e+01  -3.723930e-03  7.165e-03  1.961e+00  
     GE3    1.86e-11  5.88e-05  -3.690868e+01  -8.720926e-04  2.952e-04  1.952e+00  
     GE4    3.54e-12  2.32e-06  -3.690846e+01  2.205591e-04   8.220e-06  1.961e+00  
     GE5    1.09e-12  6.59e-07  -3.690848e+01  -1.601398e-05  1.926e-06  1.964e+00  
     GE6    9.14e-15  5.45e-09  -3.690848e+01  -5.269195e-06  1.679e-07  1.970e+00  
     GE7    6.66e-16  2.84e-10  -3.690848e+01  1.491741e-06   1.521e-08  1.749e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.330e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.330e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.174e+01     
     TOTAL-PRESSURE: -3.278e+01 KBAR
     ETOT DIFF (eV)       : -8.589e-03
     LARGEST GRAD (eV/A)  : 2.627e-01
     -------------------------------------------
     STEP OF RELAXATION : 34
     -------------------------------------------
     DONE(5.965e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.15e-10 3.21e-03  -3.691074e+01  0.000000e+00   8.439e-03  1.834e+00  
     GE2    -3.59e-11 1.35e-03  -3.691215e+01  -1.410327e-03  3.457e-03  1.944e+00  
     GE3    8.02e-12  2.80e-05  -3.691256e+01  -4.093916e-04  1.420e-04  1.961e+00  
     GE4    1.53e-12  9.59e-07  -3.691245e+01  1.074134e-04   3.006e-06  1.969e+00  
     GE5    4.27e-13  2.59e-07  -3.691246e+01  -9.345263e-06  7.572e-07  1.958e+00  
     GE6    3.87e-15  1.95e-09  -3.691246e+01  -1.754612e-06  7.117e-08  1.759e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.329e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.329e+01     -1.771e-32     
     0.000e+00      -1.771e-32     -3.177e+01     
     TOTAL-PRESSURE: -3.278e+01 KBAR
     ETOT DIFF (eV)       : -3.982e-03
     LARGEST GRAD (eV/A)  : 2.554e-01
     -------------------------------------------
     STEP OF RELAXATION : 35
     -------------------------------------------
     DONE(6.154e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -3.22e-10 9.62e-03  -3.691677e+01  0.000000e+00   2.534e-02  1.888e+00  
     GE2    -1.00e-10 4.05e-03  -3.692290e+01  -6.128768e-03  1.038e-02  1.974e+00  
     GE3    2.36e-11  8.56e-05  -3.692409e+01  -1.192113e-03  4.283e-04  1.944e+00  
     GE4    4.85e-12  3.57e-06  -3.692380e+01  2.921512e-04   1.421e-05  1.952e+00  
     GE5    1.55e-12  9.55e-07  -3.692382e+01  -1.839156e-05  2.782e-06  2.042e+00  
     GE6    1.11e-14  7.61e-09  -3.692383e+01  -8.093536e-06  2.356e-07  2.429e+00  
     GE7    -8.24e-16 4.08e-10  -3.692383e+01  2.100644e-06   1.943e-08  1.800e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.327e+01     0.000e+00      2.656e-32      
     0.000e+00      -3.327e+01     -1.771e-32     
     2.656e-32      -1.771e-32     -3.186e+01     
     TOTAL-PRESSURE: -3.280e+01 KBAR
     ETOT DIFF (eV)       : -1.136e-02
     LARGEST GRAD (eV/A)  : 2.347e-01
     -------------------------------------------
     STEP OF RELAXATION : 36
     -------------------------------------------
     DONE(6.375e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.46e-10  4.45e-03  -3.692635e+01  0.000000e+00   1.171e-02  1.851e+00  
     GE2    7.78e-11  1.87e-03  -3.692839e+01  -2.037618e-03  4.797e-03  1.965e+00  
     GE3    1.31e-11  3.88e-05  -3.692892e+01  -5.303631e-04  1.974e-04  1.958e+00  
     GE4    1.98e-12  1.25e-06  -3.692878e+01  1.377966e-04   4.140e-06  1.980e+00  
     GE5    5.99e-13  3.64e-07  -3.692879e+01  -1.123596e-05  1.056e-06  2.004e+00  
     GE6    3.06e-15  2.40e-09  -3.692880e+01  -2.903344e-06  9.830e-08  1.811e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.326e+01     0.000e+00      -1.771e-32     
     0.000e+00      -3.326e+01     0.000e+00      
     -1.771e-32     0.000e+00      -3.190e+01     
     TOTAL-PRESSURE: -3.281e+01 KBAR
     ETOT DIFF (eV)       : -4.972e-03
     LARGEST GRAD (eV/A)  : 2.257e-01
     -------------------------------------------
     STEP OF RELAXATION : 37
     -------------------------------------------
     DONE(6.566e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    1.52e-09  1.33e-02  -3.693174e+01  0.000000e+00   3.504e-02  1.839e+00  
     GE2    6.96e-10  5.61e-03  -3.694150e+01  -9.764590e-03  1.435e-02  1.942e+00  
     GE3    5.47e-11  1.20e-04  -3.694302e+01  -1.520187e-03  5.990e-04  1.952e+00  
     GE4    7.60e-12  5.60e-06  -3.694267e+01  3.547937e-04   2.512e-05  1.943e+00  
     GE5    2.17e-12  1.27e-06  -3.694269e+01  -1.829509e-05  3.691e-06  1.949e+00  
     GE6    1.32e-14  1.00e-08  -3.694270e+01  -1.095675e-05  3.012e-07  1.954e+00  
     GE7    -6.24e-16 4.86e-10  -3.694269e+01  2.680791e-06   2.000e-08  1.766e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.323e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.323e+01     0.000e+00      
     0.000e+00      0.000e+00      -3.202e+01     
     TOTAL-PRESSURE: -3.283e+01 KBAR
     ETOT DIFF (eV)       : -1.390e-02
     LARGEST GRAD (eV/A)  : 2.003e-01
     -------------------------------------------
     STEP OF RELAXATION : 38
     -------------------------------------------
     DONE(6.776e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    6.24e-10  5.84e-03  -3.694507e+01  0.000000e+00   1.534e-02  1.882e+00  
     GE2    2.82e-10  2.46e-03  -3.694784e+01  -2.771375e-03  6.289e-03  1.952e+00  
     GE3    1.98e-11  5.10e-05  -3.694848e+01  -6.313987e-04  2.603e-04  1.947e+00  
     GE4    2.40e-12  1.51e-06  -3.694831e+01  1.616446e-04   5.446e-06  1.997e+00  
     GE5    7.53e-13  4.52e-07  -3.694833e+01  -1.241235e-05  1.308e-06  1.969e+00  
     GE6    4.48e-15  2.72e-09  -3.694833e+01  -4.013831e-06  1.214e-07  1.980e+00  
     GE7    -1.36e-16 5.90e-11  -3.694833e+01  1.056198e-06   3.946e-09  1.763e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.322e+01     0.000e+00      2.822e-32      
     0.000e+00      -3.322e+01     0.000e+00      
     2.822e-32      0.000e+00      -3.207e+01     
     TOTAL-PRESSURE: -3.284e+01 KBAR
     ETOT DIFF (eV)       : -5.635e-03
     LARGEST GRAD (eV/A)  : 1.900e-01
     -------------------------------------------
     STEP OF RELAXATION : 39
     -------------------------------------------
     DONE(6.987e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.48e-09 1.76e-02  -3.694769e+01  0.000000e+00   4.612e-02  1.839e+00  
     GE2    -5.94e-10 7.41e-03  -3.696231e+01  -1.462750e-02  1.892e-02  1.938e+00  
     GE3    1.86e-11  1.59e-04  -3.696408e+01  -1.767977e-03  7.944e-04  1.968e+00  
     GE4    6.83e-12  8.79e-06  -3.696370e+01  3.829886e-04   4.278e-05  1.959e+00  
     GE5    2.35e-12  1.54e-06  -3.696371e+01  -1.491732e-05  4.468e-06  1.952e+00  
     GE6    1.30e-14  1.25e-08  -3.696373e+01  -1.292834e-05  3.423e-07  1.962e+00  
     GE7    6.38e-16  1.83e-10  -3.696372e+01  2.945808e-06   6.490e-09  1.804e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.319e+01     0.000e+00      -8.854e-33     
     0.000e+00      -3.319e+01     0.000e+00      
     -8.854e-33     0.000e+00      -3.220e+01     
     TOTAL-PRESSURE: -3.286e+01 KBAR
     ETOT DIFF (eV)       : -1.539e-02
     LARGEST GRAD (eV/A)  : 1.623e-01
     -------------------------------------------
     STEP OF RELAXATION : 40
     -------------------------------------------
     DONE(7.196e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -2.65e-10 7.33e-03  -3.696545e+01  0.000000e+00   1.923e-02  1.835e+00  
     GE2    -9.55e-11 3.09e-03  -3.696899e+01  -3.534998e-03  7.886e-03  1.948e+00  
     GE3    1.06e-11  6.36e-05  -3.696967e+01  -6.862130e-04  3.264e-04  1.964e+00  
     GE4    2.28e-12  1.71e-06  -3.696950e+01  1.715880e-04   6.854e-06  1.969e+00  
     GE5    7.79e-13  5.01e-07  -3.696951e+01  -1.251507e-05  1.446e-06  1.958e+00  
     GE6    4.44e-15  3.28e-09  -3.696952e+01  -4.754705e-06  1.340e-07  1.961e+00  
     GE7    -4.47e-16 4.71e-11  -3.696952e+01  1.160544e-06   2.702e-09  1.945e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.318e+01     0.000e+00      -8.854e-33     
     0.000e+00      -3.318e+01     0.000e+00      
     -8.854e-33     0.000e+00      -3.225e+01     
     TOTAL-PRESSURE: -3.287e+01 KBAR
     ETOT DIFF (eV)       : -5.793e-03
     LARGEST GRAD (eV/A)  : 1.519e-01
     -------------------------------------------
     STEP OF RELAXATION : 41
     -------------------------------------------
     DONE(7.410e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    2.71e-10  2.19e-02  -3.696289e+01  0.000000e+00   5.749e-02  1.841e+00  
     GE2    1.51e-10  9.24e-03  -3.698337e+01  -2.047838e-02  2.356e-02  1.983e+00  
     GE3    5.03e-11  2.00e-04  -3.698524e+01  -1.863425e-03  1.006e-03  1.965e+00  
     GE4    9.87e-12  1.33e-05  -3.698488e+01  3.587230e-04   6.759e-05  1.962e+00  
     GE5    2.74e-12  1.72e-06  -3.698489e+01  -9.533539e-06  4.999e-06  1.967e+00  
     GE6    1.74e-14  1.50e-08  -3.698490e+01  -1.319539e-05  3.438e-07  1.964e+00  
     GE7    -2.00e-15 6.59e-10  -3.698490e+01  3.050176e-06   1.883e-08  1.767e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.316e+01     0.000e+00      0.000e+00      
     0.000e+00      -3.316e+01     -1.771e-32     
     0.000e+00      -1.771e-32     -3.239e+01     
     TOTAL-PRESSURE: -3.290e+01 KBAR
     ETOT DIFF (eV)       : -1.538e-02
     LARGEST GRAD (eV/A)  : 1.227e-01
     -------------------------------------------
     STEP OF RELAXATION : 42
     -------------------------------------------
     DONE(7.619e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.85e-10 8.46e-03  -3.698552e+01  0.000000e+00   2.210e-02  1.843e+00  
     GE2    -6.54e-11 3.56e-03  -3.698954e+01  -4.019700e-03  9.070e-03  1.965e+00  
     GE3    9.74e-12  7.32e-05  -3.699020e+01  -6.578104e-04  3.776e-04  1.969e+00  
     GE4    2.05e-12  1.76e-06  -3.699004e+01  1.588879e-04   7.749e-06  1.979e+00  
     GE5    7.20e-13  4.79e-07  -3.699005e+01  -1.120695e-05  1.389e-06  1.956e+00  
     GE6    7.51e-15  4.92e-09  -3.699005e+01  -4.700079e-06  1.280e-07  1.961e+00  
     GE7    -6.03e-16 1.46e-10  -3.699005e+01  1.130235e-06   7.514e-09  1.761e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.315e+01     0.000e+00      -8.854e-33     
     0.000e+00      -3.315e+01     0.000e+00      
     -8.854e-33     0.000e+00      -3.244e+01     
     TOTAL-PRESSURE: -3.291e+01 KBAR
     ETOT DIFF (eV)       : -5.156e-03
     LARGEST GRAD (eV/A)  : 1.128e-01
     -------------------------------------------
     STEP OF RELAXATION : 43
     -------------------------------------------
     DONE(7.830e+02  SEC) : INIT SCF
     ITER   TMAG      AMAG      ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    3.04e-10  2.53e-02  -3.697672e+01  0.000000e+00   6.642e-02  1.831e+00  
     GE2    1.47e-10  1.07e-02  -3.700195e+01  -2.522468e-02  2.724e-02  1.950e+00  
     GE3    5.37e-11  2.31e-04  -3.700368e+01  -1.731240e-03  1.172e-03  1.948e+00  
     GE4    1.05e-11  1.76e-05  -3.700340e+01  2.783769e-04   9.065e-05  1.954e+00  
     GE5    2.74e-12  1.73e-06  -3.700341e+01  -5.575368e-06  5.046e-06  1.948e+00  
     GE6    2.51e-14  1.90e-08  -3.700342e+01  -1.119620e-05  3.057e-07  2.233e+00  
     GE7    -3.74e-16 7.63e-10  -3.700341e+01  2.538882e-06   1.893e-08  1.761e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -3.312e+01     0.000e+00      1.771e-32      
     0.000e+00      -3.312e+01     0.000e+00      
     1.771e-32      0.000e+00      -3.256e+01     
     TOTAL-PRESSURE: -3.293e+01 KBAR
     ETOT DIFF (eV)       : -1.336e-02
     LARGEST GRAD (eV/A)  : 8.789e-02
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               800.86         93        8.6       1e+02     %
       Driver              driver_line         800.83         1         8e+02     1e+02     %
       PW_Basis            setup_struc_factor  1.1803         43        0.027     0.15      %
       ORB_control         set_orb_tables      0.30925        1         0.31      0.039     %
       ORB_gen_tables      gen_tables          0.30925        1         0.31      0.039     %
       Ions                opt_ions            799.82         1         8e+02     1e+02     %
       ESolver_KS_LCAO     Run                 740.79         43        17        92        %
       ESolver_KS_LCAO     beforescf           155.76         43        3.6       19        %
       ESolver_KS_LCAO     beforesolver        0.27254        43        0.0063    0.034     %
       ESolver_KS_LCAO     set_matrix_grid     0.27101        43        0.0063    0.034     %
       Grid_Technique      init                0.2624         43        0.0061    0.033     %
       Grid_BigCell        grid_expansion_index0.17392        86        0.002     0.022     %
       Charge              atomic_rho          9.2667         86        0.11      1.2       %
       PW_Basis            recip2real          96.358         3630      0.027     12        %
       PW_Basis            gathers_scatterp    32.333         3630      0.0089    4         %
       Potential           init_pot            45.661         43        1.1       5.7       %
       Potential           update_from_charge  307.32         287       1.1       38        %
       Potential           cal_fixed_v         1.4098         43        0.033     0.18      %
       PotLocal            cal_fixed_v         1.3405         43        0.031     0.17      %
       Potential           cal_v_eff           305.91         287       1.1       38        %
       H_Hartree_pw        v_hartree           33.821         287       0.12      4.2       %
       PW_Basis            real2recip          120.53         4721      0.026     15        %
       PW_Basis            gatherp_scatters    43.307         4721      0.0092    5.4       %
       PotXC               cal_v_eff           269.3          287       0.94      34        %
       XC_Functional       v_xc                268.28         287       0.93      33        %
       Symmetry            rhog_symmetry       75.775         574       0.13      9.5       %
       H_Ewald_pw          compute_ewald       1.1308         43        0.026     0.14      %
       HSolverLCAO         solve               37.607         244       0.15      4.7       %
       HamiltLCAO          updateHk            19.694         488       0.04      2.5       %
       OperatorLCAO        init                19.68          976       0.02      2.5       %
       Veff                contributeHk        19.668         488       0.04      2.5       %
       Gint_interface      cal_gint            59.04          818       0.072     7.4       %
       Gint_interface      cal_gint_vlocal     19.577         488       0.04      2.4       %
       Gint_Tools          cal_psir_ylm        9.8252         1062852   9.2e-06   1.2       %
       HSolverLCAO         hamiltSolvePsiK     0.65061        488       0.0013    0.081     %
       DiagoElpa           elpa_solve          0.43684        488       0.0009    0.055     %
       ElecStateLCAO       psiToRho            17.259         244       0.071     2.2       %
       Gint_interface      cal_gint_rho        14.766         244       0.061     1.8       %
       Charge_Mixing       rhog_dot_product    0.811          244       0.0033    0.1       %
       Charge              mix_rho             28.429         201       0.14      3.5       %
       Charge              Pulay_mixing        27.104         201       0.13      3.4       %
       Charge              plain_mixing        0.18694        43        0.0043    0.023     %
       Force_Stress_LCAO   getForceStress      57.585         43        1.3       7.2       %
       Forces              cal_force_loc       2.3613         43        0.055     0.29      %
       Forces              cal_force_ew        1.1267         43        0.026     0.14      %
       Forces              cal_force_scc       3.1375         43        0.073     0.39      %
       Stress_Func         stress_loc          3.3395         43        0.078     0.42      %
       Stress_Func         stress_har          1.8808         43        0.044     0.23      %
       Stress_Func         stress_ewa          1.6309         43        0.038     0.2       %
       Stress_Func         stress_gga          19.282         43        0.45      2.4       %
       Force_LCAO_k        ftable_k            24.726         43        0.58      3.1       %
       Force_LCAO_k        cal_fvl_dphi_k      24.697         43        0.57      3.1       %
       Gint_interface      cal_gint_force      24.697         86        0.29      3.1       %
       Gint_Tools          cal_dpsir_ylm       12.516         124554    0.0001    1.6       %
       Gint_Tools          cal_dpsirr_ylm      0.90486        124554    7.3e-06   0.11      %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Thu Feb 29 12:22:48 2024
    

     FINISH Time  : Thu Feb 29 12:36:09 2024
     TOTAL  Time  : 801
     SEE INFORMATION IN : OUT.H2/
    

---

让我们来查看一下结构优化得到的 H2 键长：


```python
! cat ./ABACUS_Relax/Practice/OUT.H2/STRU_NOW.cif
```

    data_none
    
    _audit_creation_method generated by ABACUS
    
    _cell_length_a 10
    _cell_length_b 10
    _cell_length_c 10
    _cell_angle_alpha 90
    _cell_angle_beta 90
    _cell_angle_gamma 90
    
    loop_
    _atom_site_label
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    H 0 0 0.88994
    H 0 0 0.18406
    

可以看到，计算的键长为 0.706 AA。

**修改 STRU 文件中的键长与 INPUT 文件中参数，重新进行计算，分析你得到的不同结果。**
