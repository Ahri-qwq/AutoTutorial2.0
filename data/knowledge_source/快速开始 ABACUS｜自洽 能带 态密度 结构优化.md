# 快速开始 ABACUS｜自洽 能带 态密度 结构优化

<a href="https://bohrium.dp.tech/notebook/4e39364ea5c5455e94bdbdfd2fb9b037" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

<div style="color:black; background-color:#FFF3E9; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:1rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2023 @ Authors</i></b><br/>
        <i>作者：<a href="mailto:zhaotq@bjaisi.com"><b>赵天琦 📨 </b></a></i><br/>
        <i>日期：2023-04-01</i><br/>
        <i>共享协议：</a>本作品采用<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议</a>进行许可。</i>
    </p>
</div>

🎯 <b style='color:purple'>本教程旨在<u>快速掌握</u>使用 ABACUS 进行第一性原理材料模拟计算。</b>

* 一键运行，你可以快速在实践中检验你的想法。

* 丰富完善的注释，对于入门者友好。

**在 [Bohrium Notebook](https://bohrium-doc.dp.tech/docs/userguide/Notebook) 界面，你可以点击界面上方蓝色按钮 `开始连接`，选择 `abacus:3.2.0-toolkit-notebook` 镜像及任何一款节点配置，稍等片刻即可运行。**

<div style="width:auto; height:2px; background:linear-gradient(244deg,rgba(0,0,0,0) 0%,rgba(0,0,0,0.5) 50%,rgba(0,0,0,1) 100%)"></div>

## 目标

> **上手使用ABACUS做第一性原理材料模拟计算**

在学习本教程后，你将能够：

- 使用ABACUS做自洽计算
- 计算得到电子能带结构图
- 计算得到电子态密度图
- 使用ABACUS做结构优化计算

**阅读该教程【最多】约需 15 分钟，让我们开始吧！**

## 目录

* [背景](#background)
* [自洽计算](#1)
    * [准备输入文件](#1-1)
    * [运行自洽场计算](#1-2)
    * [得到总能量](#1-3)
* [电子能带结构](#2)
    * [自洽计算](#2-1)
    * [非自洽计算](#2-2)
    * [画能带图](#2-3)
* [态密度计算](#3)
    * [自洽计算](#3-1)
    * [非自洽计算](#3-2)
    * [画总态密度图](#3-3)
    * [画投影态密度图](#3-4)
* [结构优化](#4)
* [进一步阅读](#furtherreading)

## 背景 <a id ='background'></a>

密度泛函理论算法是用于第一性原理材料性质模拟的重要基础工具。目前，密度泛函理论算法已被实现在许多软件中，例如VASP，Quantum Espresso，ABIINT，CP2K等。

ABACUS是一款基于密度泛函理论发展起来的国内开源软件，支持平面波和数值原子轨道两种基矢量，主要采用模守恒赝势，功能较齐全，可适用于从小体系到上千原子大体系的电子结构优化、原子结构弛豫、分子动力学模拟等计算。

本教程将以面心立方氧化镁为例，介绍ABACUS的基本用法。这里我们采用ABACUS（原子算筹）软件的3.2.0版本。

## 1. 自洽场计算<a id ='1'></a>

自洽，是自洽场计算(SCF, self-consistent field calculation)的简称，它是密度泛函理论（DFT）算法的基础，也是DFT被称为第一性原理方法的原因。

从自洽场计算中我们可以得到一个体系的基态结构和基态能量。

基态结构使得我们可以得到体系的结构参数，例如晶胞参数，化学键键长键角等，而基态能量使得我们可以进一步得到原子受力和晶体材料的应力-应变性质等。

因此自洽场计算是开始做DFT计算的基础。

### 1.1 准备输入文件<a id ='1-1'></a>

这里我们在gitee上为大家准备好了快速上手的案例，大家可以直接使用git命令来下载。


```python
! git clone https://gitee.com/deepmodeling/colombo-ABACUS.git
```

    Cloning into 'colombo-ABACUS'...
    remote: Enumerating objects: 300, done.[K
    remote: Counting objects: 100% (118/118), done.[K
    remote: Compressing objects: 100% (81/81), done.[K
    remote: Total 300 (delta 29), reused 98 (delta 19), pack-reused 182[K
    Receiving objects: 100% (300/300), 3.79 MiB | 821.00 KiB/s, done.
    Resolving deltas: 100% (77/77), done.
    Updating files: 100% (160/160), done.
    

注意：该临时文件将在本案例末执行删除操作。

### 1.2 运行自洽场计算<a id ='1-2'></a>

运行自洽场计算需在**INPUT**中将calculation设置为**scf**(详细的参数说明请参考[文档说明](http://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html))。

这里的STRU结构为MgO的单胞(conventional cell)。

STRU文件的详细说明请参考[链接](http://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html)。

需注意：
- 在此镜像中，由于abacus是基于intel oneapi的环境安装的，我们运行abacus之前需要使用命令`. /opt/intel/oneapi/setvars.sh`装载相应的环境。
- 设置`OMP_NUM_THREADS=1`使用单线程进行计算。
- `abacus`为ABACUS可执行程序的命令。


```python
cd colombo-ABACUS/ABACUS/MgO_LCAO/SCF
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/SCF
    


```python
cat INPUT
```

    INPUT_PARAMETERS
    suffix 			MgO
    ntype 			2
    pseudo_dir 		../../PP_ORB
    orbital_dir 		../../PP_ORB
    ecutwfc 		100 
    scf_thr 		1e-6 
    basis_type 		lcao
    calculation 		scf 
    


```python
cat KPT
```

    K_POINTS
    0
    Gamma
    4 4 4 0 0 0
    


```python
cat STRU
```

    ATOMIC_SPECIES
    Mg 24.305 Mg_ONCV_PBE-1.0.upf
    O  15.999 O_ONCV_PBE-1.0.upf
    
    NUMERICAL_ORBITAL
    Mg_gga_8au_100Ry_4s2p1d.orb
    O_gga_8au_100Ry_2s2p1d.orb
    
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
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:12:21 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 54 * 54 * 54
     UNIFORM GRID DIM(BIG): 18 * 18 * 18
     DONE(0.083602   SEC) : SETUP UNITCELL
     DONE(0.0944471  SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       36              1           112         
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
     DONE(2.59466    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.653810e+03  0.000000e+00   1.172e-01  9.182e+00  
     GE2    -7.660698e+03  -6.887227e+00  9.635e-02  8.233e+00  
     GE3    -7.663893e+03  -3.195615e+00  1.673e-02  8.213e+00  
     GE4    -7.663908e+03  -1.521795e-02  2.058e-03  8.213e+00  
     GE5    -7.663910e+03  -1.464449e-03  2.765e-04  8.192e+00  
     GE6    -7.663910e+03  -1.052038e-05  2.365e-05  8.226e+00  
     GE7    -7.663910e+03  -6.776327e-08  4.643e-07  8.209e+00  
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               61.091         7         8.7       1e+02     %
       Driver              driver_line         61.056         1         61        1e+02     %
       ORB_control         set_orb_tables      2.2134         1         2.2       3.6       %
       ORB_gen_tables      gen_tables          2.2134         1         2.2       3.6       %
       ORB_table_phi       init_Table          1.528          1         1.5       2.5       %
       ORB_table_phi       cal_ST_Phi12_R      1.5131         278       0.0054    2.5       %
       ORB_table_beta      init_Table_Beta     0.48746        1         0.49      0.8       %
       ORB_table_beta      VNL_PhiBeta_R       0.48358        120       0.004     0.79      %
       ESolver_KS_LCAO     Run                 58.629         1         59        96        %
       ESolver_KS_LCAO     beforescf           0.15308        1         0.15      0.25      %
       Potential           update_from_charge  0.61232        8         0.077     1         %
       Potential           cal_v_eff           0.60904        8         0.076     1         %
       PW_Basis            real2recip          0.1081         62        0.0017    0.18      %
       PotXC               cal_v_eff           0.55798        8         0.07      0.91      %
       XC_Functional       v_xc                0.55683        8         0.07      0.91      %
       HSolverLCAO         solve               57.812         7         8.3       95        %
       HamiltLCAO          updateHk            30.389         252       0.12      50        %
       OperatorLCAO        init                27.198         504       0.054     45        %
       Veff                contributeHk        27.193         252       0.11      45        %
       Gint_interface      cal_gint            49.622         14        3.5       81        %
       Gint_interface      cal_gint_vlocal     26.86          7         3.8       44        %
       Gint_Tools          cal_psir_ylm        6.6024         81648     8.1e-05   11        %
       Gint_k              folding_vl_k        0.33298        252       0.0013    0.55      %
       Gint_k              Distri              0.10942        252       0.00043   0.18      %
       Overlap             contributeHR        0.10519        1         0.11      0.17      %
       LCAO_gen_fixedH     calculate_S_no      0.10518        1         0.11      0.17      %
       Ekin<LCAO>          contributeHR        0.10553        1         0.11      0.17      %
       Nonlocal<LCAO>      contributeHR        0.65638        1         0.66      1.1       %
       LCAO_gen_fixedH     b_NL_mu_new         0.65617        1         0.66      1.1       %
       OperatorLCAO        folding_fixed       2.3242         252       0.0092    3.8       %
       LCAO_nnr            folding_fixedH      2.3144         252       0.0092    3.8       %
       HSolverLCAO         hamiltSolvePsiK     1.512          252       0.006     2.5       %
       DiagoElpa           elpa_solve          1.4488         252       0.0057    2.4       %
       ElecStateLCAO       psiToRho            25.909         7         3.7       42        %
       LCAO_Charge         cal_dk_k            3.0993         7         0.44      5.1       %
       Gint_interface      cal_gint_rho        22.761         7         3.3       37        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:12:21 2023
     FINISH Time  : Mon Apr  3 15:13:22 2023
     TOTAL  Time  : 61
     SEE INFORMATION IN : OUT.MgO/
    

### 1.3 得到总能量<a id ='1-3'></a>


```python
! grep ETOT OUT.MgO/running_scf.log
```

     !FINAL_ETOT_IS -7663.909829526411 eV
    

## 2. 电子能带结构<a id ='2'></a>

电子的能带结构(Electronic band structure)是材料的一个重要性质，能带结构可以通过密度泛函理论算法来计算。

电子能带结构本质是电子的本征能量。得到电子的本征能量后，还可以计算得到电子有效质量，材料的光学，电学和磁性性质，因此电子能带结构也是分析材料性质的重要工具。

做电子能带结构计算首先需要通过**自洽计算**得到体系的电子密度，然后基于自洽计算的电子密度做一步**非自洽计算**，以得到布里渊区里高对称k点路径上电子的本征值，即能带结构。

### 2.1 自洽计算<a id ='2-1'></a>

这里的STRU结构为MgO的原胞(primitive cell)。


```python
cd ../../MgO_LCAO/band_structure
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/band_structure
    


```python
cat STRU
```

    ATOMIC_SPECIES
    Mg 24.305 Mg_ONCV_PBE-1.0.upf upf201
    O 15.999 O_ONCV_PBE-1.0.upf upf201
    
    NUMERICAL_ORBITAL
    Mg_gga_8au_100Ry_4s2p1d.orb
    O_gga_8au_100Ry_2s2p1d.orb
    
    LATTICE_CONSTANT
    8.087214649032902
    
    ATOMIC_POSITIONS
    Direct
    
    Mg #label
    0 #magnetism
    1 #number of atoms
    0  0  0  m  0  0  0
    
    O #label
    0 #magnetism
    1 #number of atoms
    0.5  0.5  0.5  m  0  0  0
    


```python
cp INPUT_scf INPUT && cp KPT_scf KPT
```


```python
cat INPUT
```

    INPUT_PARAMETERS
    suffix 			MgO
    ntype 			2
    pseudo_dir 		../../PP_ORB
    orbital_dir 		../../PP_ORB
    ecutwfc 		100 
    scf_thr 		1e-6 
    basis_type 		lcao
    calculation 		scf 
    out_chg			1
    symmetry		1
    latname			fcc
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:13:25 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 40 * 40 * 40
     UNIFORM GRID DIM(BIG): 10 * 10 * 10
     DONE(0.0639761  SEC) : SETUP UNITCELL
     DONE(0.0997776  SEC) : SYMMETRY
     DONE(0.173368   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       29              1           28          
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Mg      4s2p1d-8au      15          1           
     O       2s2p1d-8au      13          1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.52608    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.913453e+03  0.000000e+00   1.172e-01  3.067e+00  
     GE2    -1.915174e+03  -1.721656e+00  9.638e-02  2.854e+00  
     GE3    -1.915974e+03  -7.993561e-01  1.673e-02  2.857e+00  
     GE4    -1.915977e+03  -3.794312e-03  2.061e-03  2.856e+00  
     GE5    -1.915978e+03  -3.672082e-04  2.767e-04  2.848e+00  
     GE6    -1.915978e+03  -2.633162e-06  2.365e-05  2.853e+00  
     GE7    -1.915978e+03  -1.705373e-08  4.561e-07  2.848e+00  
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               22.823         9         2.5       1e+02     %
       Driver              driver_line         22.788         1         23        1e+02     %
       ORB_control         set_orb_tables      2.1871         1         2.2       9.6       %
       ORB_gen_tables      gen_tables          2.1871         1         2.2       9.6       %
       ORB_table_phi       init_Table          1.5112         1         1.5       6.6       %
       ORB_table_phi       cal_ST_Phi12_R      1.4965         278       0.0054    6.6       %
       ORB_table_beta      init_Table_Beta     0.48066        1         0.48      2.1       %
       ORB_table_beta      VNL_PhiBeta_R       0.47688        120       0.004     2.1       %
       ESolver_KS_LCAO     Run                 20.345         1         20        89        %
       Potential           update_from_charge  0.22907        8         0.029     1         %
       Potential           cal_v_eff           0.22801        8         0.029     1         %
       PotXC               cal_v_eff           0.21081        8         0.026     0.92      %
       XC_Functional       v_xc                0.21048        8         0.026     0.92      %
       HSolverLCAO         solve               19.402         7         2.8       85        %
       HamiltLCAO          updateHk            11.516         203       0.057     50        %
       OperatorLCAO        init                10.811         406       0.027     47        %
       Veff                contributeHk        10.81          203       0.053     47        %
       Gint_interface      cal_gint            18.225         14        1.3       80        %
       Gint_interface      cal_gint_vlocal     10.745         7         1.5       47        %
       Gint_Tools          cal_psir_ylm        2.4918         14000     0.00018   11        %
       Nonlocal<LCAO>      contributeHR        0.16289        1         0.16      0.71      %
       LCAO_gen_fixedH     b_NL_mu_new         0.16283        1         0.16      0.71      %
       OperatorLCAO        folding_fixed       0.48838        203       0.0024    2.1       %
       LCAO_nnr            folding_fixedH      0.48693        203       0.0024    2.1       %
       HSolverLCAO         hamiltSolvePsiK     0.17813        203       0.00088   0.78      %
       DiagoElpa           elpa_solve          0.12734        203       0.00063   0.56      %
       ElecStateLCAO       psiToRho            7.7075         7         1.1       34        %
       LCAO_Charge         cal_dk_k            0.22205        7         0.032     0.97      %
       Gint_interface      cal_gint_rho        7.4794         7         1.1       33        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:13:25 2023
     FINISH Time  : Mon Apr  3 15:13:48 2023
     TOTAL  Time  : 23
     SEE INFORMATION IN : OUT.MgO/
    

### 2.2 非自洽计算<a id ='2-2'></a>

运行非自洽计算需在**INPUT**中将calculation设置为**nscf**，为了输出能带结构文件，还需设置**out_band**为**true**或者**1**。高对称k点路径在**KPT**文件中设置，选用**Line**模式。


```python
cp INPUT_nscf INPUT && cp KPT_nscf KPT
```


```python
cat INPUT
```

    INPUT_PARAMETERS
    suffix 			MgO
    ntype 			2
    pseudo_dir 		../../PP_ORB
    orbital_dir 		../../PP_ORB
    ecutwfc 		100 
    scf_thr 		1e-6 
    basis_type 		lcao
    calculation 		nscf 
    init_chg		file
    out_band		1
    out_dos			1
    latname			fcc
    symmetry		0
    


```python
cat KPT
```

    K_POINTS			
    7
    Line		
    0.375  0.375  0.750  20		# K
    0.000  0.000  0.000  20		# G
    0.500  0.000  0.500  20		# X
    0.625  0.250  0.625  20         # U
    0.500  0.500  0.500  20		# L
    0.000  0.000  0.000  20		# G
    0.500  0.250  0.750  1		# W
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:13:51 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 40 * 40 * 40
     UNIFORM GRID DIM(BIG): 10 * 10 * 10
     DONE(0.0518589  SEC) : SETUP UNITCELL
     DONE(0.0620326  SEC) : INIT K-POINTS
     ---------------------------------------------------------
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       121             1           28          
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Mg      4s2p1d-8au      15          1           
     O       2s2p1d-8au      13          1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     NONSELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : file
     NON-SELF CONSISTENT CALCULATIONS
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               8.4843         5         1.7       1e+02     %
       Driver              driver_line         8.4595         1         8.5       1e+02     %
       ORB_control         set_orb_tables      2.178          1         2.2       26        %
       ORB_gen_tables      gen_tables          2.178          1         2.2       26        %
       ORB_table_phi       init_Table          1.5037         1         1.5       18        %
       ORB_table_phi       cal_ST_Phi12_R      1.4886         278       0.0054    18        %
       ORB_table_beta      init_Table_Beta     0.47611        1         0.48      5.6       %
       ORB_table_beta      VNL_PhiBeta_R       0.47223        120       0.0039    5.6       %
       ESolver_KS_LCAO     othercalculation    2.2761         1         2.3       27        %
       HSolverLCAO         solve               2.1905         1         2.2       26        %
       HamiltLCAO          updateHk            2.0828         121       0.017     25        %
       OperatorLCAO        init                1.5669         242       0.0065    18        %
       Veff                contributeHk        1.5661         121       0.013     18        %
       Gint_interface      cal_gint            1.5275         1         1.5       18        %
       Gint_interface      cal_gint_vlocal     1.5275         1         1.5       18        %
       Gint_Tools          cal_psir_ylm        0.1782         1000      0.00018   2.1       %
       Nonlocal<LCAO>      contributeHR        0.16534        1         0.17      1.9       %
       LCAO_gen_fixedH     b_NL_mu_new         0.16529        1         0.17      1.9       %
       OperatorLCAO        folding_fixed       0.29656        121       0.0025    3.5       %
       LCAO_nnr            folding_fixedH      0.59353        242       0.0025    7         %
       HSolverLCAO         hamiltSolvePsiK     0.10727        121       0.00089   1.3       %
       ModuleIO            nscf_band           0.51157        1         0.51      6         %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:13:51 2023
     FINISH Time  : Mon Apr  3 15:14:00 2023
     TOTAL  Time  : 9
     SEE INFORMATION IN : OUT.MgO/
    

### 2.3 画能带图<a id ='2-3'></a>


```python
cd OUT.MgO
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/band_structure/OUT.MgO
    


```python
cp ../KPT ../config.json .
```

在running_scf.log里面搜索E_Fermi，找到费米能级的值


```python
cat config.json
```

    {
    	"bandfile" : "BANDS_1.dat",
    	"efermi" : 8.27,
    	"energy_range"	: [-4, 9],
    	"kptfile" : "KPT"
    }
    
    
    
    


```python
!/opt/miniconda3/bin/abacus-plot -b
```


```python
import matplotlib.pyplot as plt
from PIL import Image
plt.axis('off')
plt.imshow(Image.open('band.png'))
plt.show()
```



## 3. 态密度计算<a id ='3'></a>

态密度(Density of states)是另外一个用于分析材料性质的重要性质。

它本质上描述的是在不同能量上电子态的密度，利用态密度可以得到材料的光学性质等。

通过将波函数投影到每个原子的不同轨道，我们还可以分析不同原子轨道对电子态密度的贡献，用于从原子和电子角度分析和设计材料性质。

### 3.1 自洽计算<a id ='3-1'></a>

与能带计算的自洽计算相同。


```python
cd ../../../MgO_LCAO/dos
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/dos
    


```python
cp INPUT_scf INPUT && cp KPT_scf KPT
```


```python
cat KPT
```

    K_POINTS
    0
    Gamma
    8 8 8 0 0 0
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:14:08 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 40 * 40 * 40
     UNIFORM GRID DIM(BIG): 10 * 10 * 10
     DONE(0.0685985  SEC) : SETUP UNITCELL
     DONE(0.104701   SEC) : SYMMETRY
     DONE(0.178032   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       29              1           28          
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Mg      4s2p1d-8au      15          1           
     O       2s2p1d-8au      13          1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.51577    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.913453e+03  0.000000e+00   1.172e-01  3.064e+00  
     GE2    -1.915174e+03  -1.721656e+00  9.638e-02  2.841e+00  
     GE3    -1.915974e+03  -7.993561e-01  1.673e-02  2.842e+00  
     GE4    -1.915977e+03  -3.794312e-03  2.061e-03  2.842e+00  
     GE5    -1.915978e+03  -3.672082e-04  2.767e-04  2.829e+00  
     GE6    -1.915978e+03  -2.633162e-06  2.365e-05  2.848e+00  
     GE7    -1.915978e+03  -1.705373e-08  4.561e-07  2.845e+00  
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               22.745         9         2.5       1e+02     %
       Driver              driver_line         22.709         1         23        1e+02     %
       ORB_control         set_orb_tables      2.1734         1         2.2       9.6       %
       ORB_gen_tables      gen_tables          2.1734         1         2.2       9.6       %
       ORB_table_phi       init_Table          1.506          1         1.5       6.6       %
       ORB_table_phi       cal_ST_Phi12_R      1.4913         278       0.0054    6.6       %
       ORB_table_beta      init_Table_Beta     0.47219        1         0.47      2.1       %
       ORB_table_beta      VNL_PhiBeta_R       0.46835        120       0.0039    2.1       %
       ESolver_KS_LCAO     Run                 20.271         1         20        89        %
       Potential           update_from_charge  0.22668        8         0.028     1         %
       Potential           cal_v_eff           0.22561        8         0.028     0.99      %
       PotXC               cal_v_eff           0.20855        8         0.026     0.92      %
       XC_Functional       v_xc                0.20823        8         0.026     0.92      %
       HSolverLCAO         solve               19.347         7         2.8       85        %
       HamiltLCAO          updateHk            11.549         203       0.057     51        %
       OperatorLCAO        init                10.84          406       0.027     48        %
       Veff                contributeHk        10.838         203       0.053     48        %
       Gint_interface      cal_gint            18.179         14        1.3       80        %
       Gint_interface      cal_gint_vlocal     10.774         7         1.5       47        %
       Gint_Tools          cal_psir_ylm        2.5668         14000     0.00018   11        %
       Nonlocal<LCAO>      contributeHR        0.16355        1         0.16      0.72      %
       LCAO_gen_fixedH     b_NL_mu_new         0.16349        1         0.16      0.72      %
       OperatorLCAO        folding_fixed       0.49287        203       0.0024    2.2       %
       LCAO_nnr            folding_fixedH      0.4914         203       0.0024    2.2       %
       HSolverLCAO         hamiltSolvePsiK     0.16938        203       0.00083   0.74      %
       DiagoElpa           elpa_solve          0.12265        203       0.0006    0.54      %
       ElecStateLCAO       psiToRho            7.6272         7         1.1       34        %
       LCAO_Charge         cal_dk_k            0.2167         7         0.031     0.95      %
       Gint_interface      cal_gint_rho        7.4048         7         1.1       33        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:14:08 2023
     FINISH Time  : Mon Apr  3 15:14:31 2023
     TOTAL  Time  : 23
     SEE INFORMATION IN : OUT.MgO/
    

### 3.2 非自洽计算<a id ='3-2'></a>

做态密度计算中的非自洽计算需设置更加密的k点网格。为了输出能带结构文件，还需设置**out_dos**为**true**或者**1**。


```python
cp INPUT_nscf INPUT && cp KPT_nscf KPT
```


```python
cat KPT
```

    K_POINTS
    0
    Gamma
    10 10 10 0 0 0
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:14:34 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 40 * 40 * 40
     UNIFORM GRID DIM(BIG): 10 * 10 * 10
     DONE(0.0614797  SEC) : SETUP UNITCELL
     DONE(0.0851821  SEC) : INIT K-POINTS
     ---------------------------------------------------------
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       504             1           28          
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Mg      4s2p1d-8au      15          1           
     O       2s2p1d-8au      13          1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     NONSELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : file
     NON-SELF CONSISTENT CALCULATIONS
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               12.769         5         2.6       1e+02     %
       Driver              driver_line         12.741         1         13        1e+02     %
       ORB_control         set_orb_tables      2.253          1         2.3       18        %
       ORB_gen_tables      gen_tables          2.253          1         2.3       18        %
       ORB_table_phi       init_Table          1.5323         1         1.5       12        %
       ORB_table_phi       cal_ST_Phi12_R      1.5173         278       0.0055    12        %
       ORB_table_beta      init_Table_Beta     0.50964        1         0.51      4         %
       ORB_table_beta      VNL_PhiBeta_R       0.50553        120       0.0042    4         %
       ESolver_KS_LCAO     othercalculation    3.6888         1         3.7       29        %
       HSolverLCAO         solve               3.582          1         3.6       28        %
       HamiltLCAO          updateHk            3.1582         504       0.0063    25        %
       OperatorLCAO        init                1.7253         1008      0.0017    14        %
       Veff                contributeHk        1.722          504       0.0034    13        %
       Gint_interface      cal_gint            1.5624         1         1.6       12        %
       Gint_interface      cal_gint_vlocal     1.5624         1         1.6       12        %
       Gint_Tools          cal_psir_ylm        0.17806        1000      0.00018   1.4       %
       Gint_k              folding_vl_k        0.1592         504       0.00032   1.2       %
       Nonlocal<LCAO>      contributeHR        0.16364        1         0.16      1.3       %
       LCAO_gen_fixedH     b_NL_mu_new         0.16358        1         0.16      1.3       %
       OperatorLCAO        folding_fixed       1.2155         504       0.0024    9.5       %
       LCAO_nnr            folding_fixedH      2.4493         1008      0.0024    19        %
       HSolverLCAO         hamiltSolvePsiK     0.42219        504       0.00084   3.3       %
       DiagoElpa           elpa_solve          0.31119        504       0.00062   2.4       %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:14:34 2023
     FINISH Time  : Mon Apr  3 15:14:46 2023
     TOTAL  Time  : 12
     SEE INFORMATION IN : OUT.MgO/
    

### 3.3 画总态密度图<a id ='3-3'></a>

总态密度(Total density of states)是指在一个能量范围内电子在每个能量点上总的态密度值。

这里的“总”是相对下面的投影态密度(projected density of states)讲的，而投影态密度指的是在一个能量范围内将电子波函数投影到每个原子的每个轨道得到的态密度值。


```python
cd OUT.MgO
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/dos/OUT.MgO
    


```python
cp ../config_tdos.json ./config.json
```


```python
cat config.json
```

    {
        "tdosfile": "TDOS",
        "efermi": 8.27771540465,
        "energy_range": [-5,7],
        "dos_range": [0,5],
        "figsize":[14,10],
        "tdosfig": "tdos.png"
    }
    


```python
!/opt/miniconda3/bin/abacus-plot -d
```


```python
import matplotlib.pyplot as plt
from PIL import Image
plt.axis('off')
plt.imshow(Image.open('tdos.png'))
plt.show()
```



### 3.4 画投影态密度图<a id ='3-4'></a>


```python
cp ../config_pdos.json ./config.json
```


```python
cat config.json
```

    {
        "pdosfile": "PDOS",
        "efermi": 8.27771540465,
        "energy_range": [-5,7],
        "dos_range": [0,5],
        "figsize":[14,10],
        "species":{"Mg":[0,1,2],"O":[0,1,2]},
        "pdosfig": "pdos.png"
    }
    


```python
!/opt/miniconda3/bin/abacus-plot -d -p -o
```


```python
import matplotlib.pyplot as plt
from PIL import Image
plt.axis('off')
plt.imshow(Image.open('pdos.png'))
plt.show()
```



## 4. 结构优化<a id ='4'></a>

结构优化是从一个一般的结构出发，利用优化算法得到基态结构的过程。对于周期性体系，优化计算可以根据需要针对体系的不同自由度来进行优化。

在ABACUS中，如果只优化原子位置，而晶胞固定不变，则需在INPUT中将calculation设置为**relax**;

如果需要同时做原子位置优化和晶胞优化，则需要将calculation设置为**cell-relax**；

此外，需要注意的参数还有原子受力的收敛标准**force_thr_ev**和应力收敛标准**stress_thr**的设置；

更多的优化设置，请阅读[ABACUS Documentation](http://abacus.deepmodeling.com/en/latest/)。关于结构优化，可参考[文档](http://abacus.deepmodeling.com/en/latest/quick_start/hands_on.html#running-geometry-optimization)中对这部分的详细说明。


```python
cd ../../../MgO_LCAO/optimization
```

    /data/bohrium-notebook/Proem/colombo-ABACUS/ABACUS/MgO_LCAO/optimization
    


```python
cat INPUT
```

    INPUT_PARAMETERS
    suffix 			MgO
    ntype 			2
    pseudo_dir 		../../PP_ORB
    orbital_dir 		../../PP_ORB
    ecutwfc 		100 
    scf_thr 		1e-6 
    basis_type 		lcao
    calculation 		cell-relax
    force_thr_ev		0.01
    stress_thr		5
    relax_nmax		100
    out_stru		1
    
    
    
    


```python
!. /opt/intel/oneapi/setvars.sh && OMP_NUM_THREADS=1 abacus
```

     
    :: initializing oneAPI environment ...
       dash: SH_VERSION = unknown
       args: Using "$@" for setvars.sh arguments: 
    :: compiler -- latest
    :: debugger -- latest
    :: dev-utilities -- latest
    :: mkl -- latest
    :: mpi -- latest
    :: tbb -- latest
    :: vtune -- latest
    :: oneAPI environment initialized ::
     
    WARNING: Total thread number on this node mismatches with hardware availability. This may cause poor performance.
    Info: Local MPI proc number: 1,OpenMP thread number: 1,Total thread number: 1,Local thread limit: 52
                                                                                         
                                  ABACUS v3.1                                            
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
    
     Mon Apr  3 15:15:01 2023
     MAKE THE DIR         : OUT.MgO/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: number valence electrons > 2 for Mg: [Ne] 3s2
     Please make sure the pseudopotential file is what you need
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 54 * 54 * 54
     UNIFORM GRID DIM(BIG): 18 * 18 * 18
     DONE(0.0830077  SEC) : SETUP UNITCELL
     DONE(0.0940724  SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Cell relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       36              1           112         
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
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.56792    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.653810e+03  0.000000e+00   1.172e-01  9.129e+00  
     GE2    -7.660698e+03  -6.887227e+00  9.635e-02  8.247e+00  
     GE3    -7.663893e+03  -3.195615e+00  1.673e-02  8.224e+00  
     GE4    -7.663908e+03  -1.521795e-02  2.058e-03  8.142e+00  
     GE5    -7.663910e+03  -1.464449e-03  2.765e-04  8.164e+00  
     GE6    -7.663910e+03  -1.052038e-05  2.365e-05  8.149e+00  
     GE7    -7.663910e+03  -6.776327e-08  4.643e-07  8.150e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     1.202e+00      4.377e-12      1.201e-12      
     4.377e-12      1.202e+00      1.541e-12      
     1.201e-12      1.541e-12      1.202e+00      
     TOTAL-PRESSURE: 1.202e+00 KBAR
     LARGEST GRAD (eV/A)  : 2.106e-11
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               79.81          7         11        1e+02     %
       Driver              driver_line         79.772         1         80        1e+02     %
       ORB_control         set_orb_tables      2.1614         1         2.2       2.7       %
       ORB_gen_tables      gen_tables          2.1614         1         2.2       2.7       %
       ORB_table_phi       init_Table          1.4899         1         1.5       1.9       %
       ORB_table_phi       cal_ST_Phi12_R      1.4749         278       0.0053    1.8       %
       ORB_table_beta      init_Table_Beta     0.476          1         0.48      0.6       %
       ORB_table_beta      VNL_PhiBeta_R       0.47219        120       0.0039    0.59      %
       ESolver_KS_LCAO     Run                 58.369         1         58        73        %
       ESolver_KS_LCAO     beforescf           0.15329        1         0.15      0.19      %
       Potential           update_from_charge  0.61303        8         0.077     0.77      %
       Potential           cal_v_eff           0.6098         8         0.076     0.76      %
       PW_Basis            real2recip          0.11709        68        0.0017    0.15      %
       PotXC               cal_v_eff           0.55894        8         0.07      0.7       %
       XC_Functional       v_xc                0.55778        8         0.07      0.7       %
       HSolverLCAO         solve               57.551         7         8.2       72        %
       HamiltLCAO          updateHk            30.066         252       0.12      38        %
       OperatorLCAO        init                26.874         504       0.053     34        %
       Veff                contributeHk        26.87          252       0.11      34        %
       Gint_interface      cal_gint            61.396         15        4.1       77        %
       Gint_interface      cal_gint_vlocal     26.539         7         3.8       33        %
       Gint_Tools          cal_psir_ylm        6.5462         81648     8e-05     8.2       %
       Gint_k              folding_vl_k        0.33013        252       0.0013    0.41      %
       Gint_k              Distri              0.10859        252       0.00043   0.14      %
       Overlap             contributeHR        0.10385        1         0.1       0.13      %
       LCAO_gen_fixedH     calculate_S_no      0.10385        1         0.1       0.13      %
       Ekin<LCAO>          contributeHR        0.105          1         0.11      0.13      %
       Nonlocal<LCAO>      contributeHR        0.66882        1         0.67      0.84      %
       LCAO_gen_fixedH     b_NL_mu_new         2.1443         2         1.1       2.7       %
       OperatorLCAO        folding_fixed       2.3136         252       0.0092    2.9       %
       LCAO_nnr            folding_fixedH      2.3033         252       0.0091    2.9       %
       HSolverLCAO         hamiltSolvePsiK     1.521          252       0.006     1.9       %
       DiagoElpa           elpa_solve          1.4572         252       0.0058    1.8       %
       ElecStateLCAO       psiToRho            25.963         7         3.7       33        %
       LCAO_Charge         cal_dk_k            3.1115         7         0.44      3.9       %
       Gint_interface      cal_gint_rho        22.803         7         3.3       29        %
       Force_Stress_LCAO   getForceStress      19.022         1         19        24        %
       Force_LCAO_k        ftable_k            18.865         1         19        24        %
       Force_LCAO_k        allocate_k          2.3928         1         2.4       3         %
       Force_LCAO_k        cal_foverlap_k      0.11431        1         0.11      0.14      %
       Force_LCAO_k        cal_edm_2d          0.10043        1         0.1       0.13      %
       Local_Orbital_Chargecal_dm_R            0.18659        2         0.093     0.23      %
       Force_LCAO_k        cal_fvl_dphi_k      12.054         1         12        15        %
       Gint_interface      cal_gint_force      12.054         1         12        15        %
       Gint_Tools          cal_dpsir_ylm       4.9742         5832      0.00085   6.2       %
       Gint_Tools          cal_dpsirr_ylm      0.69358        5832      0.00012   0.87      %
       Force_LCAO_k        cal_fvnl_dbeta_k_new4.1995         1         4.2       5.3       %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Mon Apr  3 15:15:01 2023
     FINISH Time  : Mon Apr  3 15:16:21 2023
     TOTAL  Time  : 80
     SEE INFORMATION IN : OUT.MgO/
    


```python
# 在案例完成后，删除该案例中下载的 colombo-ABACUS 文件夹
import os

# 获取当前工作目录的绝对路径
current_dir = os.getcwd()

if current_dir.endswith("colombo-ABACUS/ABACUS/MgO_LCAO/optimization"):
    # 获取 colombo-ABACUS 上级目录的绝对路径
    parent_dir = os.path.abspath(os.path.join(current_dir, '../../../..'))
    os.chdir(parent_dir)

! if [ -d "colombo-ABACUS" ]; then rm -r colombo-ABACUS; echo "删除成功"; fi
```

    删除成功
    

## 进一步阅读 <a id='furtherreading'></a>

如果您希望更深入学习ABACUS，本节提供有关该主题的更多资源。

**ABACUS项目**

- [ABACUS Homepage](http://abacus.ustc.edu.cn/main.htm)
- [ABACUS Documentation](http://abacus.deepmodeling.com/en/latest/)
- [ABACUS on Github](https://github.com/deepmodeling/abacus-develop)
- [ABACUS on Gitee](https://gitee.com/deepmodeling/abacus-develop)



<a href="https://bohrium.dp.tech/notebook/4e39364ea5c5455e94bdbdfd2fb9b037" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>
