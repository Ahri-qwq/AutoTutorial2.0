# ABACUS+pyatb 能带反折叠计算

<div style="color:black; background-color:#FFF3E9; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:1rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2023 @ Authors</i></b><br/>
        <i>作者：<b><a href="19820221153869@stu.xmu.edu.cn">王璟航 📨 </a></b></i><br/>
        <i>日期：2023-09-01</i><br/>
        <i>共享协议：</a>本作品采用<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议</a>进行许可。</i>
    </p>
</div>

<a href="" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

<font color="grey">


推荐镜像：`abacus:3.2.3` （注意，本教程涉及环境搭建，选择不同的镜像很可能导致报错！）\
推荐计算资源：`CPU` \
预计运行用时：30min\
内容：本教程主要介绍如何用ABACUS+pyatb做能带反折叠计算。 \
使用方式：您可在 Bohrium Notebook上直接运行。您可以点击界面上方蓝色按钮 `开始连接`，选择 `abacus:3.2.3` 镜像及`c16_m32_cpu`款节点配置，稍等片刻即可运行。如您遇到任何问题，请联系 [bohrium@dp.tech](mailto:bohrium@dp.tech) 。 \
共享协议：本作品采用[知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议](https://creativecommons.org/licenses/by-nc-sa/4.0/)进行许可。

</font>

## 一、介绍

PYATB(Python ab initio tight binding simuation package),基于从头算紧束缚哈密顿量计算电子结构和相关性质的python包。\
github库地址：https://github.com/pyatb/pyatb ；安装及参数教程：https://pyatb.github.io/pyatb/ \
本文主要使用其中的两个模块：BandStructure和Bandunfolding。\
本文主要内容包括：**，abacus自洽计算，pyatb软件包安装及所需环境安装，能带计算及unfolding**。



#### 1.Unfolding原理：

首先，之所以需要对能带进行unfolding（反折叠），是因为当我们想要得到含有杂质、缺陷或低组分合金材料的性质时，通常的做法是将仅含少量原子的PC（primitive cell）扩胞为LC（Large cell）后进行替位或间隙式掺杂。然而此做法会导致原胞的空间平移对称性被破坏，在原胞的相差倒格矢的k点间引入耦合。反折叠就是将K（超胞k点）用原胞的一组完备的布洛赫波函数进行展开，通过谱权重（与展开系数直接相关）便可以得到不同k点对K点的贡献。详细原理请见文献：https://doi.org/10.1016/j.commatsci.2022.111656

 <center>  $$ |\Psi_{K N}\rangle = \sum\limits_{n,k_p}|\psi_{k_p~n}\rangle\langle\psi_{k_p~n}|\Psi_{KN}|$$
    $$ |\Psi_{K N}\rangle = \frac{1}{\sqrt{\mathcal{N}}}\sum\limits_{R}\sum\limits_{\mu,i}C_{\mu,i}(K)e^{i~K\cdot R}\phi_\mu (r-\tau_{\alpha i}-R) $$

#### 2.输入文件介绍：

本文中的案例文件放在https://gitee.com/Luc1anooo/abacus_wjh.git， 更多案例请见 https://github.com/pyatb/pyatb/tree/main/examples


```python
cd ~
```

    /root
    


```python
! git clone https://gitee.com/Luc1anooo/abacus_wjh.git
```

    Cloning into 'abacus_wjh'...
    remote: Enumerating objects: 133, done.[K
    remote: Counting objects: 100% (133/133), done.[K
    remote: Compressing objects: 100% (122/122), done.[K
    remote: Total 133 (delta 40), reused 0 (delta 0), pack-reused 0[K
    Receiving objects: 100% (133/133), 2.76 MiB | 1.88 MiB/s, done.
    Resolving deltas: 100% (40/40), done.
    

**使用pyatb软件包进行计算时需要准备两个文件夹，一个用于做自洽计算，命名为abacus；一个用于使用pyatb软件包的模块，命名为pyatb**

注：1.pyatb文件夹中需要准备赝势和轨道文件 2.pyatb文件夹中的输入参数文件命名为Input，以免和abacus的输入文件“INPUT”混淆 3.两个文件夹中的STRU文件相同


```python
! tree abacus_wjh -L 3
```

    [01;34mabacus_wjh[0m
    ├── [01;34mGeC[0m
    │   ├── [01;34mabacus[0m
    │   │   ├── INPUT
    │   │   ├── KPT
    │   │   └── STRU
    │   └── [01;34mpyatb[0m
    │       ├── C_ONCV_PBE-1.0.upf
    │       ├── C_gga_8au_100Ry_2s2p1d.orb
    │       ├── Ge_ONCV_PBE-1.0.upf
    │       ├── Ge_gga_8au_100Ry_2s2p2d.orb
    │       ├── Input
    │       └── STRU
    ├── [01;34mPP_ORB[0m
    │   ├── C_ONCV_PBE-1.0.upf
    │   ├── C_gga_8au_100Ry_2s2p1d.orb
    │   ├── Ge_ONCV_PBE-1.0.upf
    │   └── Ge_gga_8au_100Ry_2s2p2d.orb
    ├── README.en.md
    ├── README.md
    └── [01;34meigen3[0m
        └── [01;31meigen-3.4.0.tar.gz[0m
    
    5 directories, 16 files
    

## 二、自洽计算

自洽计算教程见[快速开始 ABACUS｜自洽 能带 态密度 结构优化](https://nb.bohrium.dp.tech/detail/4641406377#2)


```python
cd ~/abacus_wjh/GeC/abacus/
```

    /root/abacus_wjh/GeC/abacus
    


```python
cat INPUT
```

    INPUT_PARAMETERS
    
    # System variables
    suffix                GeC
    pseudo_dir 		../../PP_ORB
    orbital_dir 	../../PP_ORB
    ntype                 2
    calculation           scf
    symmetry              1
    init_chg              atomic
    
    # Plane wave related variables
    ecutwfc               50
    
    # Electronic structure
    basis_type            lcao
    ks_solver             genelpa
    nspin                 1
    gamma_only            0
    scf_nmax              300
    scf_thr               1e-7
    nbands                500
    
    # Variables related to output information
    out_chg               1
    out_level             ie
    out_mat_hs2           1
    out_mat_r             1
    
    


```python
cat KPT
```

    K_POINTS
    0
    Gamma
    1 1 1 0 0 0
    


```python
cat STRU #为了节约时间，离子弛豫已经达到
```

    ATOMIC_SPECIES
    Ge 72.63 Ge_ONCV_PBE-1.0.upf upf201
    C 12.011 C_ONCV_PBE-1.0.upf upf201
    
    NUMERICAL_ORBITAL
    Ge_gga_8au_100Ry_2s2p2d.orb
    C_gga_8au_100Ry_2s2p1d.orb
    
    LATTICE_CONSTANT
    1.8897259886
    
    LATTICE_VECTORS
    11.3149995804 0 0 #latvec1
    0 11.3149995804 0 #latvec2
    0 0 11.3149995804 #latvec3
    
    ATOMIC_POSITIONS
    Direct
    
    Ge #label
    0 #magnetism
    63 #number of atoms
    1.06541996772e-20  2.04995561528e-20  8.85139430939e-21  m  1  1  1
    0.000161097984705  0.250005684723  0.250005684723  m  1  1  1
    0.250005684723  0.000161097984705  0.250005684723  m  1  1  1
    0.250005684723  0.250005684723  0.000161097984705  m  1  1  1
    0.374879555379  0.1254052439  0.374879555379  m  1  1  1
    0.124695128451  0.124695128451  0.124695128451  m  1  1  1
    0.1254052439  0.374879555379  0.374879555379  m  1  1  1
    0.374879555379  0.374879555379  0.1254052439  m  1  1  1
    0.5  9.56096752097e-22  1.41608590038e-20  m  1  1  1
    0.50006231713  0.255431956156  0.255431956156  m  1  1  1
    0.749994315277  0.999838902015  0.250005684723  m  1  1  1
    0.749994315277  0.250005684723  0.999838902015  m  1  1  1
    0.872413789244  0.127586210756  0.374894019293  m  1  1  1
    0.625105980707  0.127586210756  0.127586210756  m  1  1  1
    0.60810021808  0.39189978192  0.39189978192  m  1  1  1
    0.872413789244  0.374894019293  0.127586210756  m  1  1  1
    0  0.5  0  m  1  1  1
    0.999838902015  0.749994315277  0.250005684723  m  1  1  1
    0.255431956156  0.50006231713  0.255431956156  m  1  1  1
    0.250005684723  0.749994315277  0.999838902015  m  1  1  1
    0.39189978192  0.60810021808  0.39189978192  m  1  1  1
    0.127586210756  0.625105980707  0.127586210756  m  1  1  1
    0.127586210756  0.872413789244  0.374894019293  m  1  1  1
    0.374894019293  0.872413789244  0.127586210756  m  1  1  1
    0.5  0.5  0  m  1  1  1
    0.49993768287  0.744568043844  0.255431956156  m  1  1  1
    0.744568043844  0.49993768287  0.255431956156  m  1  1  1
    0.749994315277  0.749994315277  0.000161097984705  m  1  1  1
    0.8745947561  0.625120444621  0.374879555379  m  1  1  1
    0.625120444621  0.625120444621  0.1254052439  m  1  1  1
    0.625120444621  0.8745947561  0.374879555379  m  1  1  1
    0.875304871549  0.875304871549  0.124695128451  m  1  1  1
    0  9.94619022158e-21  0.5  m  1  1  1
    0.999838902015  0.250005684723  0.749994315277  m  1  1  1
    0.250005684723  0.999838902015  0.749994315277  m  1  1  1
    0.255431956156  0.255431956156  0.50006231713  m  1  1  1
    0.374894019293  0.127586210756  0.872413789244  m  1  1  1
    0.127586210756  0.127586210756  0.625105980707  m  1  1  1
    0.127586210756  0.374894019293  0.872413789244  m  1  1  1
    0.39189978192  0.39189978192  0.60810021808  m  1  1  1
    0.5  0  0.5  m  1  1  1
    0.49993768287  0.255431956156  0.744568043844  m  1  1  1
    0.749994315277  0.000161097984705  0.749994315277  m  1  1  1
    0.744568043844  0.255431956156  0.49993768287  m  1  1  1
    0.875304871549  0.124695128451  0.875304871549  m  1  1  1
    0.625120444621  0.1254052439  0.625120444621  m  1  1  1
    0.625120444621  0.374879555379  0.8745947561  m  1  1  1
    0.8745947561  0.374879555379  0.625120444621  m  1  1  1
    0  0.5  0.5  m  1  1  1
    0.000161097984705  0.749994315277  0.749994315277  m  1  1  1
    0.255431956156  0.49993768287  0.744568043844  m  1  1  1
    0.255431956156  0.744568043844  0.49993768287  m  1  1  1
    0.374879555379  0.625120444621  0.8745947561  m  1  1  1
    0.1254052439  0.625120444621  0.625120444621  m  1  1  1
    0.124695128451  0.875304871549  0.875304871549  m  1  1  1
    0.374879555379  0.8745947561  0.625120444621  m  1  1  1
    0.50006231713  0.744568043844  0.744568043844  m  1  1  1
    0.744568043844  0.50006231713  0.744568043844  m  1  1  1
    0.744568043844  0.744568043844  0.50006231713  m  1  1  1
    0.872413789244  0.625105980707  0.872413789244  m  1  1  1
    0.60810021808  0.60810021808  0.60810021808  m  1  1  1
    0.625105980707  0.872413789244  0.872413789244  m  1  1  1
    0.872413789244  0.872413789244  0.625105980707  m  1  1  1
    
    C #label
    0 #magnetism
    1 #number of atoms
    0.5  0.5  0.5  m  1  1  1


```python
!mpirun -n 16 abacus
```

                                                                                         
                                  ABACUS v3.2.3
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 9df02abab (Fri May 19 05:10:31 2023)
    
     Fri Sep  1 02:08:44 2023
     MAKE THE DIR         : OUT.GeC/
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 4 for Ge: [Ar] 3d10 4s2 4p2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM     : 100 * 100 * 100
     UNIFORM GRID DIM(BIG): 25 * 25 * 25
     DONE(0.0692687  SEC) : SETUP UNITCELL
     DONE(0.120118   SEC) : SYMMETRY
     DONE(0.248293   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       1               16          1147        
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Ge      2s2p2d-8au      18          63          
     C       2s2p1d-8au      13          1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(3.67992    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.214879e+05  0.000000e+00   4.154e-02  4.879e+00  
     GE2    -1.214874e+05  4.984676e-01   3.170e-02  3.594e+00  
     GE3    -1.214870e+05  4.653295e-01   1.868e-02  3.549e+00  
     GE4    -1.214872e+05  -1.908288e-01  4.611e-03  3.786e+00  
     GE5    -1.214874e+05  -1.819687e-01  2.323e-03  3.586e+00  
     GE6    -1.214874e+05  -7.376873e-03  1.001e-03  3.583e+00  
     GE7    -1.214874e+05  -4.828837e-04  3.002e-04  3.663e+00  
     GE8    -1.214874e+05  -2.935446e-04  9.217e-05  3.743e+00  
     GE9    -1.214874e+05  -1.882474e-05  3.697e-05  3.754e+00  
     GE10   -1.214874e+05  -3.990765e-06  1.578e-05  3.561e+00  
     GE11   -1.214874e+05  -1.667488e-06  5.803e-06  3.768e+00  
     GE12   -1.214874e+05  -6.828144e-07  2.018e-06  3.562e+00  
     GE13   -1.214874e+05  -1.463386e-07  7.830e-07  3.572e+00  
     GE14   -1.214874e+05  -1.781901e-08  4.262e-07  3.542e+00  
     GE15   -1.214874e+05  1.145861e-08   2.336e-07  3.745e+00  
     GE16   -1.214874e+05  8.761012e-09   1.334e-07  3.542e+00  
     GE17   -1.214874e+05  6.162407e-09   4.935e-08  3.550e+00  
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               98.597         9         11        1e+02     %
       Driver              driver_line         98.58          1         99        1e+02     %
       NOrbital_Lm         extra_uniform       0.2136         221       0.00097   0.22      %
       Mathzone_Add1       Cubic_Spline_Interpolation0.10593        221       0.00048   0.11      %
       ORB_control         set_orb_tables      2.881          1         2.9       2.9       %
       ORB_gen_tables      gen_tables          2.881          1         2.9       2.9       %
       ORB_table_phi       init_Table          1.9237         1         1.9       2         %
       ORB_table_phi       cal_ST_Phi12_R      14.48          766       0.019     15        %
       ORB_table_beta      init_Table_Beta     0.80341        1         0.8       0.81      %
       ORB_table_beta      VNL_PhiBeta_R       0.79846        158       0.0051    0.81      %
       ORB_gaunt_table     Get_Gaunt_SH        71.461         7285      0.0098    72        %
       Ions                opt_ions            95.209         1         95        97        %
       ESolver_KS_LCAO     Run                 95.209         1         95        97        %
       ESolver_KS_LCAO     beforescf           0.30069        1         0.3       0.3       %
       PW_Basis            recip2real          0.26739        92        0.0029    0.27      %
       PW_Basis            gathers_scatterp    0.13292        92        0.0014    0.13      %
       Potential           update_from_charge  1.3031         18        0.072     1.3       %
       Potential           cal_v_eff           1.2992         18        0.072     1.3       %
       H_Hartree_pw        v_hartree           0.14918        18        0.0083    0.15      %
       PW_Basis            real2recip          0.47889        142       0.0034    0.49      %
       PW_Basis            gatherp_scatters    0.28188        142       0.002     0.29      %
       PotXC               cal_v_eff           1.1429         18        0.063     1.2       %
       XC_Functional       v_xc                1.1383         18        0.063     1.2       %
       Symmetry            rho_symmetry        1.0659         18        0.059     1.1       %
       HSolverLCAO         solve               60.127         17        3.5       61        %
       HamiltLCAO          updateHk            19.562         17        1.2       20        %
       OperatorLCAO        init                18.047         34        0.53      18        %
       Veff                contributeHk        18.04          17        1.1       18        %
       Gint_interface      cal_gint            28.153         34        0.83      29        %
       Gint_interface      cal_gint_vlocal     15.148         17        0.89      15        %
       Gint_Tools          cal_psir_ylm        6.9915         42500     0.00016   7.1       %
       Gint_k              folding_vl_k        2.8919         17        0.17      2.9       %
       Gint_k              Distri              2.6139         17        0.15      2.7       %
       Overlap             contributeHR        0.13056        1         0.13      0.13      %
       LCAO_gen_fixedH     calculate_S_no      0.13056        1         0.13      0.13      %
       Ekin<LCAO>          contributeHR        0.13096        1         0.13      0.13      %
       Nonlocal<LCAO>      contributeHR        0.82751        1         0.83      0.84      %
       LCAO_gen_fixedH     b_NL_mu_new         0.82659        1         0.83      0.84      %
       OperatorLCAO        folding_fixed       0.42597        17        0.025     0.43      %
       LCAO_nnr            folding_fixedH      0.40381        17        0.024     0.41      %
       HSolverLCAO         hamiltSolvePsiK     6.6462         17        0.39      6.7       %
       OperatorLCAO        get_hs_pointers     31.666         9         3.5       32        %
       DiagoElpa           elpa_solve          6.5898         17        0.39      6.7       %
       ElecStateLCAO       psiToRho            33.918         17        2         34        %
       elecstate           cal_dm              0.35598        17        0.021     0.36      %
       psiMulPsiMpi        pdgemm              0.34283        17        0.02      0.35      %
        Local_Orbital_wfc  wfc_2d_to_grid      0.20742        18        0.012     0.21      %
       LCAO_Charge         cal_dk_k            16.272         17        0.96      17        %
       Gint_interface      cal_gint_rho        13.005         17        0.77      13        %
       Charge              mix_rho             0.20219        16        0.013     0.21      %
       Charge              Pulay_mixing        0.19376        16        0.012     0.2       %
       ModuleIO            output_HS_R         6.7545         1         6.8       6.9       %
       ModuleIO            save_HSR_sparse     5.416          1         5.4       5.5       %
       cal_r_overlap_R     init                14.227         1         14        14        %
       cal_r_overlap_R     out_rR_other        10.284         1         10        10        %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Fri Sep  1 02:08:44 2023
     FINISH Time  : Fri Sep  1 02:10:22 2023
     TOTAL  Time  : 98
     SEE INFORMATION IN : OUT.GeC/
    


```python
!grep -i efermi OUT.GeC/running_scf.log  #需要手动得到自洽计算后的材料的费米能级，用于后续赋值给pyatb文件夹的Input文件中参数fermi_energy
```

     EFERMI = 11.08291763890926 eV
    


```python
cd ./OUT.GeC
```

    /root/abacus_wjh/GeC/abacus/OUT.GeC
    


```python
cp data-HR-sparse_SPIN0.csr data-SR-sparse_SPIN0.csr data-rR-sparse.csr ../../pyatb
```

## 三、pyatb软件包安装与环境搭建

根据安装教程，我们需要：\
1.安装pybind11和mpi4py模块 \
2.安装eigen库 \
3.克隆pyatb软件包 \
4.修改siteconfig.py中的路径与本地匹配


```python
cd ~
```

    /root
    


```python
! git clone https://github.com/pyatb/pyatb.git
```

    Cloning into 'pyatb'...
    remote: Enumerating objects: 965, done.[K
    remote: Counting objects: 100% (117/117), done.[K
    remote: Compressing objects: 100% (77/77), done.[K
    remote: Total 965 (delta 47), reused 106 (delta 40), pack-reused 848[K
    Receiving objects: 100% (965/965), 57.09 MiB | 1.22 MiB/s, done.
    Resolving deltas: 100% (229/229), done.
    

安装前可先使用pip list命令检查自己的环境中是否已经安装这两个模块。对于abacus:3.2.3镜像，需要自己安装。


```python
pip install pybind11
```

    Looking in indexes: https://pypi.tuna.tsinghua.edu.cn/simple
    Collecting pybind11
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/06/55/9f73c32dda93fa4f539fafa268f9504e83c489f460c380371d94296126cd/pybind11-2.11.1-py3-none-any.whl (227 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m227.7/227.7 kB[0m [31m5.2 MB/s[0m eta [36m0:00:00[0ma [36m0:00:01[0m
    [?25hInstalling collected packages: pybind11
    Successfully installed pybind11-2.11.1
    [33mWARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.pypa.io/warnings/venv[0m[33m
    [0mNote: you may need to restart the kernel to use updated packages.
    


```python
pip list | grep pybind11
```

    pybind11                 2.11.1
    Note: you may need to restart the kernel to use updated packages.
    


```python
pip install mpi4py
```

    Looking in indexes: https://pypi.tuna.tsinghua.edu.cn/simple
    Collecting mpi4py
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/bc/f2/749af7fd0e7703ddca6cea525ab40f26c3ca6cbe6c23658441c6f9705860/mpi4py-3.1.4.tar.gz (2.5 MB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m2.5/2.5 MB[0m [31m19.5 MB/s[0m eta [36m0:00:00[0ma [36m0:00:01[0m
    [?25h  Installing build dependencies ... [?25ldone
    [?25h  Getting requirements to build wheel ... [?25ldone
    [?25h  Preparing metadata (pyproject.toml) ... [?25ldone
    [?25hBuilding wheels for collected packages: mpi4py
      Building wheel for mpi4py (pyproject.toml) ... [?25ldone
    [?25h  Created wheel for mpi4py: filename=mpi4py-3.1.4-cp310-cp310-linux_x86_64.whl size=632036 sha256=65ff0a590b5025f4ef723003219fb867299360406088c20c7e3ec3a0c4cba076
      Stored in directory: /root/.cache/pip/wheels/58/40/60/82a1db2d10eb254532456c6dd135c0124267e56798110ddc64
    Successfully built mpi4py
    Installing collected packages: mpi4py
    Successfully installed mpi4py-3.1.4
    [33mWARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.pypa.io/warnings/venv[0m[33m
    [0mNote: you may need to restart the kernel to use updated packages.
    


```python
pip list | grep mpi4py
```

    mpi4py                   3.1.4
    Note: you may need to restart the kernel to use updated packages.
    


```python
cd /root/abacus_wjh/eigen3
```

    /root/abacus_wjh/eigen3
    


```python
!tar -xzf eigen-3.4.0.tar.gz
```

将文件复制到siteconfig.py的默认路径中：（当然也可以修改默认路径为自定义路径）


```python
cp -r eigen-3.4.0 /usr/local/include/eigen3
```


```python
cd ~/pyatb
```

    /root/pyatb
    


```python
pip install ./
```

    Looking in indexes: https://pypi.tuna.tsinghua.edu.cn/simple
    Processing /root/pyatb
      Preparing metadata (setup.py) ... [?25ldone
    [?25hRequirement already satisfied: numpy in /opt/mamba/lib/python3.10/site-packages (from pyatb==1.0.0) (1.24.1)
    Requirement already satisfied: scipy in /opt/mamba/lib/python3.10/site-packages (from pyatb==1.0.0) (1.9.3)
    Requirement already satisfied: matplotlib in /opt/mamba/lib/python3.10/site-packages/matplotlib-3.7.1-py3.10-linux-x86_64.egg (from pyatb==1.0.0) (3.7.1)
    Requirement already satisfied: contourpy>=1.0.1 in /opt/mamba/lib/python3.10/site-packages/contourpy-1.0.7-py3.10-linux-x86_64.egg (from matplotlib->pyatb==1.0.0) (1.0.7)
    Requirement already satisfied: cycler>=0.10 in /opt/mamba/lib/python3.10/site-packages/cycler-0.11.0-py3.10.egg (from matplotlib->pyatb==1.0.0) (0.11.0)
    Requirement already satisfied: fonttools>=4.22.0 in /opt/mamba/lib/python3.10/site-packages/fonttools-4.39.4-py3.10.egg (from matplotlib->pyatb==1.0.0) (4.39.4)
    Requirement already satisfied: kiwisolver>=1.0.1 in /opt/mamba/lib/python3.10/site-packages/kiwisolver-1.4.4-py3.10-linux-x86_64.egg (from matplotlib->pyatb==1.0.0) (1.4.4)
    Requirement already satisfied: packaging>=20.0 in /opt/mamba/lib/python3.10/site-packages (from matplotlib->pyatb==1.0.0) (22.0)
    Requirement already satisfied: pillow>=6.2.0 in /opt/mamba/lib/python3.10/site-packages/Pillow-9.5.0-py3.10-linux-x86_64.egg (from matplotlib->pyatb==1.0.0) (9.5.0)
    Requirement already satisfied: pyparsing>=2.3.1 in /opt/mamba/lib/python3.10/site-packages/pyparsing-3.1.0b2-py3.10.egg (from matplotlib->pyatb==1.0.0) (3.1.0b2)
    Requirement already satisfied: python-dateutil>=2.7 in /opt/mamba/lib/python3.10/site-packages (from matplotlib->pyatb==1.0.0) (2.8.2)
    Requirement already satisfied: six>=1.5 in /opt/mamba/lib/python3.10/site-packages (from python-dateutil>=2.7->matplotlib->pyatb==1.0.0) (1.16.0)
    Building wheels for collected packages: pyatb
      Building wheel for pyatb (setup.py) ... [?25ldone
    [?25h  Created wheel for pyatb: filename=pyatb-1.0.0-cp310-cp310-linux_x86_64.whl size=1071176 sha256=fb9bdac8bc2e3e6f282f4b8fbaf43f939c3c3f7b2b3177f856c2cff6575b8606
      Stored in directory: /tmp/pip-ephem-wheel-cache-hkge1dh2/wheels/2a/a9/13/33af4c71b5d32f01441404fadc8d21bfcd3f0b93d200240770
    Successfully built pyatb
    Installing collected packages: pyatb
    Successfully installed pyatb-1.0.0
    [33mWARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.pypa.io/warnings/venv[0m[33m
    [0mNote: you may need to restart the kernel to use updated packages.
    

## 四、能带计算及unfolding

**首先，有三个操作需要手动完成：\
1.从自洽计算的输出文件中读取费米能级，赋值给pyatb文件夹中的Input文件中的参数fermi_energy\
2.将自洽计算中得到的包含紧束缚哈密顿量信息的三个文件data-HR-sparse_SPIN0.csr, data-SR-sparse_ SPIN0.csr, data-rR-sparse.csr复制到pyatb文件夹中\
3.在运行完pyatb命令后，根据能带图，在Out/Bandunfolding文件中将参数energy_range修改为合适的值**\
以上操作在https://pyatb.github.io/pyatb/ 中都有详细介绍。


```python
cd ~/abacus_wjh/GeC/pyatb/
```

    /root/abacus_wjh/GeC/pyatb
    


```python
ls #在第二步（自洽计算）的最后已经将三个文件复制完成
```

    C_ONCV_PBE-1.0.upf           STRU
    C_gga_8au_100Ry_2s2p1d.orb   data-HR-sparse_SPIN0.csr
    Ge_ONCV_PBE-1.0.upf          data-SR-sparse_SPIN0.csr
    Ge_gga_8au_100Ry_2s2p2d.orb  data-rR-sparse.csr
    Input
    


```python
cat Input #已经将上面所得的费米能级赋值给参数fermi_energy
```

    INPUT_PARAMETERS
    {
    
        nspin                           1   
        package                         ABACUS
        fermi_energy                    11.08291763890926
        fermi_energy_unit               eV  
        HR_route                        data-HR-sparse_SPIN0.csr
        SR_route                        data-SR-sparse_SPIN0.csr
        rR_route                        data-rR-sparse.csr
        HR_unit                         Ry  
        rR_unit                         Bohr
    
      
    }
    
    LATTICE
    {
        lattice_constant                1.8897162
        lattice_constant_unit           Bohr
        lattice_vector
    11.3149995804 0 0 #latvec1
    0 11.3149995804 0 #latvec2
    0 0 11.3149995804 #latvec3
    }
    
    BAND_STRUCTURE
    {
        wf_collect                     0
        kpoint_mode                    line
        kpoint_num                     5
        high_symmetry_kpoint		
        0.50000  0.50000 0.5000 50  # L
        0.00000  0.00000 0.0000 50  # G
        0.50000  0.00000 0.5000 50  # X
        0.37500 -0.37500 0.0000 50  # K
        0.00000  0.00000 0.0000 1    # G
    }
    
    BANDUNFOLDING
    {
        stru_file                       STRU
        ecut                            50
        band_range                      1 500
        m_matrix                        -2    2    2   2  -2  2  2  2  -2
    
        kpoint_mode                     line
        kpoint_num                      5
        high_symmetry_kpoint
        0.50000  0.50000 0.5000 50  # L
        0.00000  0.00000 0.0000 50  # G
        0.50000  0.00000 0.5000 50  # X
        0.37500 -0.37500 0.0000 50  # K
        0.00000  0.00000 0.0000 1    # G
    }
    


```python
import os

# Set the environment variable
os.environ['LD_PRELOAD'] = '/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so'

# Verify the environment variable is set
print(os.environ['LD_PRELOAD'])
```

    /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so
    

这段代码的作用在安装教程的最后有提到。如果不运行这段代码，就会“顺利”地碰到报错：undefined symbol: mkl_sparse_optimize_bsr_trsm_i8。


```python
cd ~/abacus_wjh/GeC/pyatb/
```

    /root/abacus_wjh/GeC/pyatb
    


```python
!mpirun -n 16 pyatb
```

对于能带结构和unfolding计算，默认是没有实时计算过程输出在屏幕上的。对于此案例，计算时间约为1min，耐心等待即可。

计算完成后得到的Out文件如下：


```python
! tree Out -L 2
```

    [01;34mOut[0m
    ├── [01;34mBand_Structure[0m
    │   ├── band.dat
    │   ├── band.pdf
    │   ├── high_symmetry_kpoint.dat
    │   ├── kpt.dat
    │   └── plot_band.py
    ├── [01;34mBandunfolding[0m
    │   ├── kpt.dat
    │   ├── plot_unfold.py
    │   └── spectral_weight.dat
    ├── input.json
    └── running.log
    
    2 directories, 10 files
    

对于unfolding文件夹，需要手动设置plot_unfold.py中的参数energy_range和能带图中的能量范围匹配：


```python
!sed -i 's/energy_range = \[-4, 6\]/energy_range = [-14, 12]/g' plot_unfold.py
```


```python
!python plot_unfold.py
```


```python
ls
```

    kpt.dat  plot_unfold.py  spectral_weight.dat  unfold.pdf
    

完成！
