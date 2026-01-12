# ABACUS+ShengBTE 计算晶格热导率

<a href="" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

<font color="grey">


推荐镜像：`abacus-user-guide:3.3.2` \
推荐计算资源：`CPU` \
内容：本教程主要介绍如何用ABACUS+ShengBTE计算晶格热导率。 \
使用方式：您可在 Bohrium Notebook上直接运行。您可以点击界面上方蓝色按钮 `开始连接`，选择 `abacus-user-guide:3.3.2` 镜像及`c4_m8_cpu`款节点配置，稍等片刻即可运行。如您遇到任何问题，请联系 [bohrium@dp.tech](mailto:bohrium@dp.tech) 。 \
共享协议：本作品采用[知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议](https://creativecommons.org/licenses/by-nc-sa/4.0/)进行许可。

</font>

本notebook改编自abacus使用指南，更多信息详见[这里](https://github.com/MCresearch/abacus-user-guide)。

**请注意：运行本notebook需要选择** `Kernel`**为Bash**

## 一、介绍

本教程旨在介绍采用 ABACUS做密度泛函理论计算，并且结合ShengBTE软件计算晶格的热导率的流程。其中，整个计算过程中还用到了：1）采用Phonopy程序来计算二阶力常数，2）采用 ASE 程序进行原子结构的转换，3）采用ShengBTE的thirdorder程序计算三阶力常数，4）最后使用ShengBTE来计算材料的晶格热导率。

上述提到了一些需要结合的外部软件，这里推荐大家阅读这些软件的相关文档和说明：

ShengBTE：https://bitbucket.org/sousaw/shengbte/src/master/

phonopy：http://abacus.deepmodeling.com/en/latest/advanced/interface/phonopy.html

ASE：http://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html

thirdorder: https://bitbucket.org/sousaw/thirdorder/src/master/

## 二、准备

ABACUS的软件包中提供了一个 ABACUS+ShengBTE 计算晶格热导率的算例，可以从Gitee上[下载](https://gitee.com/mcresearch/abacus-user-guide/tree/master/examples/interface_ShengBTE)，或者在 linux 终端执行如下命令`git clone https://gitee.com/mcresearch/abacus-user-guide.git`得到算例，**在镜像中已经提前下载了本例程，大家可以直接使用。**

算例中包含采用数值原子轨道的LCAO（Linear Combination of Atomic Orbitals）和采用平面波基矢量的 PW（Plane Wave，平面波）两个文件夹。每个文件夹下分别又包含了`2nd`、`3rd`和`shengbte`这三个文件夹，分别保存了使用phonopy计算二阶力常数（`2nd`）、thirdorder计算三阶力常数（`3rd`）和`ShengBTE`计算晶格热导率（shengbte）的相关文件。

## 三、流程

以`LCAO`文件夹为例，我们这里提供的测试案例是包含2个原子的金刚石结构Si结构，采用的模守恒赝势是`Si_ONCV_PBE-1.0.upf`，以及原子轨道文件采用的是 `Si_gga_7au_100Ry_2s2p1d.orb`（GGA泛函，7 au截断半径，100 Ry能量截断，以及包含2s2p1d的DZP轨道）。

### 1. 计算二阶力常数

要计算二阶力常数，除了 ABACUS 之外，还需要结合 Phonopy 和 ASE（镜像中已提前安装）。

在镜像中默认是 python 3.10.6 的环境，采用conda安装Phonopy


```bash
echo -e "y\n" | conda install -c conda-forge phonopy
```

#### 1.1 结构优化

做晶格热导率计算之前要先对模拟的材料体系的进行原子构型的优化。下面是采用 ABACUS 做结构优化（relax）后得到的原子构型文件`STRU`。在这个例子里，为了简单起见，结构优化过程采用的是 2\*2\*2 的布里渊区 k 点采样，平面波的能量截断值 ecut（LCAO 里面也用到了平面波基矢量）为 100 Ry，注意实际计算中应该要采用更收敛的k点采样。

注意：第一行 Si 的质量 28.0855 在计算中不起作用。

#### 1.2 计算二阶力常数

首先，进入`2nd`文件夹， 调用 Phonopy 软件产生需要计算的超胞及相应微扰的多个原子构型，命令如下：

其中 setting.conf 文件的内容为：


```bash
cd /abacus-user-guide/examples/interface_ShengBTE/LCAO/2nd && cat setting.conf
```

    DIM = 2 2 2
    ATOM_NAME = Si
    


```bash
phonopy setting.conf --abacus -d
```

            _
      _ __ | |__   ___  _ __   ___   _ __  _   _
     | '_ \| '_ \ / _ \| '_ \ / _ \ | '_ \| | | |
     | |_) | | | | (_) | | | | (_) || |_) | |_| |
     | .__/|_| |_|\___/|_| |_|\___(_) .__/ \__, |
     |_|                            |_|    |___/
                                          2.20.0
    
    Compiled with OpenMP support (max 4 threads).
    Python version 3.10.6
    Spglib version 2.0.2
    
    "setting.conf" was read as phonopy configuration file.
    Calculator interface: abacus
    Crystal structure was read from "STRU".
    Unit of length: au
    Displacements creation mode
    Settings:
      Supercell: [2 2 2]
    Number of symmetry operations in supercell: 384
    Use -v option to watch primitive cell, unit cell, and supercell structures.
    
    "phonopy_disp.yaml" and supercells have been created.
    
    Summary of calculation was written in "phonopy_disp.yaml".
                     _
       ___ _ __   __| |
      / _ \ '_ \ / _` |
     |  __/ | | | (_| |
      \___|_| |_|\__,_|
    
    

这里我们采用的Si的例子只需要产生1个微扰构型 STRU-001 即可，对所有微扰构型（这里Si的例子只有1个）进行SCF计算（SCF代表 Self-consistent field，这里代表进行密度泛函理论的电子迭代自洽计算）获得原子受力，算完之后用以下命令产生`FORCE_SET`文件：


```bash
mpirun -n 4 abacus
```

                                                                                         
                                  ABACUS v3.3.2
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: e39b50efe (Fri Aug 18 16:14:25 2023 +0800)
    
     Tue Aug 22 11:16:42 2023
     MAKE THE DIR         : OUT.DIA-50/
     UNIFORM GRID DIM     : 96 * 96 * 96
     UNIFORM GRID DIM(BIG): 24 * 24 * 24
     DONE(0.104537   SEC) : SETUP UNITCELL
     DONE(0.155984   SEC) : SYMMETRY
     DONE(0.284879   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       6               4           208         
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Si      2s2p1d-7au      13          16          
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.83309    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -1.712832e+03  0.000000e+00   1.468e-01  3.016e+00  
     GE2    -1.712947e+03  -1.153429e-01  3.635e-02  2.481e+00  
     GE3    -1.712949e+03  -2.368474e-03  2.783e-03  2.643e+00  
     GE4    -1.712949e+03  6.023439e-05   1.179e-03  2.567e+00  
     GE5    -1.712950e+03  -3.517155e-04  8.774e-05  2.487e+00  
     GE6    -1.712950e+03  -1.863484e-06  3.634e-05  2.683e+00  
     GE7    -1.712950e+03  -3.357846e-07  8.614e-06  2.647e+00  
     GE8    -1.712950e+03  8.730656e-09   2.682e-06  2.542e+00  
     GE9    -1.712950e+03  7.298716e-09   1.205e-06  2.471e+00  
     GE10   -1.712950e+03  5.380891e-09   2.496e-07  2.542e+00  
     ><><><><><><><><><><><><><><><><><><><><><><
     TOTAL-STRESS (KBAR):
     ><><><><><><><><><><><><><><><><><><><><><><
     -5.542e+01     1.991e-01      1.991e-01      
     1.991e-01      -5.542e+01     6.105e-03      
     1.991e-01      6.105e-03      -5.542e+01     
     TOTAL-PRESSURE: -5.542e+01 KBAR
    
      |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------
                           total               36.241         9         4         1e+02     %
       Driver              driver_line         36.212         1         36        1e+02     %
       ORB_control         read_orb_first      0.12401        1         0.12      0.34      %
       LCAO_Orbitals       Read_Orbitals       0.10835        1         0.11      0.3       %
       ORB_control         set_orb_tables      0.95583        1         0.96      2.6       %
       ORB_gen_tables      gen_tables          0.95583        1         0.96      2.6       %
       ORB_table_phi       init_Table          0.49321        1         0.49      1.4       %
       ORB_table_phi       cal_ST_Phi12_R      0.49008        70        0.007     1.4       %
       ORB_table_beta      init_Table_Beta     0.18565        1         0.19      0.51      %
       ORB_table_beta      VNL_PhiBeta_R       0.18476        26        0.0071    0.51      %
       Ions                opt_ions            34.733         1         35        96        %
       ESolver_KS_LCAO     Run                 26.419         1         26        73        %
       ESolver_KS_LCAO     beforescf           0.33547        1         0.34      0.93      %
       PW_Basis            recip2real          0.4075         60        0.0068    1.1       %
       PW_Basis            gathers_scatterp    0.18545        60        0.0031    0.51      %
       Potential           init_pot            0.19096        1         0.19      0.53      %
       Potential           update_from_charge  2.2187         11        0.2       6.1       %
       Potential           cal_v_eff           2.2104         11        0.2       6.1       %
       H_Hartree_pw        v_hartree           0.22705        11        0.021     0.63      %
       PW_Basis            real2recip          0.80024        92        0.0087    2.2       %
       PW_Basis            gatherp_scatters    0.4748         92        0.0052    1.3       %
       PotXC               cal_v_eff           1.9713         11        0.18      5.4       %
       XC_Functional       v_xc                1.9655         11        0.18      5.4       %
       Symmetry            rho_symmetry        0.47357        11        0.043     1.3       %
       HSolverLCAO         solve               22.996         10        2.3       63        %
       HamiltLCAO          updateHk            10.166         60        0.17      28        %
       OperatorLCAO        init                9.5304         120       0.079     26        %
       Veff                contributeHk        9.5288         60        0.16      26        %
       Gint_interface      cal_gint            23.762         21        1.1       66        %
       Gint_interface      cal_gint_vlocal     9.0094         10        0.9       25        %
       Gint_Tools          cal_psir_ylm        5.5046         69120     8e-05     15        %
       Gint_k              folding_vl_k        0.51927        60        0.0087    1.4       %
       Gint_k              Distri              0.4549         60        0.0076    1.3       %
       Nonlocal<LCAO>      contributeHR        0.18389        1         0.18      0.51      %
       LCAO_gen_fixedH     b_NL_mu_new         0.61683        2         0.31      1.7       %
       OperatorLCAO        folding_fixed       0.37848        60        0.0063    1         %
       LCAO_nnr            folding_fixedH      0.37374        60        0.0062    1         %
       HSolverLCAO         hamiltSolvePsiK     3.152          60        0.053     8.7       %
       DiagoElpa           elpa_solve          3.0949         60        0.052     8.5       %
       ElecStateLCAO       psiToRho            9.6773         10        0.97      27        %
       elecstate           cal_dm              0.1036         11        0.0094    0.29      %
       LCAO_Charge         cal_dk_k            0.45646        10        0.046     1.3       %
       Gint_interface      cal_gint_rho        8.6161         10        0.86      24        %
       Charge              mix_rho             0.23411        9         0.026     0.65      %
       Charge              Pulay_mixing        0.22376        9         0.025     0.62      %
       Force_Stress_LCAO   getForceStress      8.3136         1         8.3       23        %
       Force_LCAO_k        ftable_k            7.9764         1         8         22        %
       Force_LCAO_k        allocate_k          0.715          1         0.71      2         %
       Force_LCAO_k        cal_fvl_dphi_k      6.137          1         6.1       17        %
       Gint_interface      cal_gint_force      6.137          1         6.1       17        %
       Gint_Tools          cal_dpsir_ylm       3.7506         3456      0.0011    10        %
       Gint_Tools          cal_dpsirr_ylm      0.53103        3456      0.00015   1.5       %
       Force_LCAO_k        cal_fvnl_dbeta_k_new1.0662         1         1.1       2.9       %
     ----------------------------------------------------------------------------------------
    
     START  Time  : Tue Aug 22 11:16:42 2023
     FINISH Time  : Tue Aug 22 11:17:19 2023
     TOTAL  Time  : 37
     SEE INFORMATION IN : OUT.DIA-50/
    


```bash
phonopy -f OUT.DIA-50/running_scf.log
```

            _
      _ __ | |__   ___  _ __   ___   _ __  _   _
     | '_ \| '_ \ / _ \| '_ \ / _ \ | '_ \| | | |
     | |_) | | | | (_) | | | | (_) || |_) | |_| |
     | .__/|_| |_|\___/|_| |_|\___(_) .__/ \__, |
     |_|                            |_|    |___/
                                          2.20.0
    
    Compiled with OpenMP support (max 4 threads).
    Python version 3.10.6
    Spglib version 2.0.2
    
    Calculator interface: abacus
    Displacements were read from "phonopy_disp.yaml".
    1. Drift force of "OUT.DIA-50/running_scf.log" to be subtracted
      0.00000000   0.00000000   0.00000000
    "FORCE_SETS" has been created.
                     _
       ___ _ __   __| |
      / _ \ '_ \ / _` |
     |  __/ | | | (_| |
      \___|_| |_|\__,_|
    
    

小技巧：在ABACUS的输入文件`INPUT`中可以设置变量stru_file，该变量对应的原子构型文件为`STRU-001`则ABACUS可以直接读取该结构文件。

下一步，设置`band.conf`文件,其内容如下（具体参数含义可以参考[notebook链接](https://nb.bohrium.dp.tech/detail/8741867512)）



```bash
cat band.conf
```

    ATOM_NAME = Si
    DIM = 2 2 2
    MESH = 8 8 8
    PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
    BAND = 0.0 0.0 0.0  0.5 0.0 0.5  0.625  0.25  0.625, 0.375 0.375 0.75  00 0.0 0.0  0.5 0.5 0.5
    BAND_POINTS = 101
    BAND_CONNECTION = .TRUE.
    FORCE_CONSTANTS = WRITE
    FULL_FORCE_CONSTANTS = .TRUE.
    

计算得到声子谱以及二阶力常数：


```bash
phonopy -p band.conf --abacus
```

            _
      _ __ | |__   ___  _ __   ___   _ __  _   _
     | '_ \| '_ \ / _ \| '_ \ / _ \ | '_ \| | | |
     | |_) | | | | (_) | | | | (_) || |_) | |_| |
     | .__/|_| |_|\___/|_| |_|\___(_) .__/ \__, |
     |_|                            |_|    |___/
                                          2.20.0
    
    Compiled with OpenMP support (max 4 threads).
    Python version 3.10.6
    Spglib version 2.0.2
    
    "band.conf" was read as phonopy configuration file.
    Calculator interface: abacus
    Crystal structure was read from "STRU".
    Unit of length: au
    Band structure and mesh sampling mode
    Settings:
      Sampling mesh: [8 8 8]
      Supercell: [2 2 2]
      Primitive matrix:
        [1. 0. 0.]
        [0. 1. 0.]
        [0. 0. 1.]
    Number of symmetry operations in supercell: 384
    Use -v option to watch primitive cell, unit cell, and supercell structures.
    
    Forces and displacements were read from "FORCE_SETS".
    Computing force constants...
    Max drift of force constants: -0.000000 (zz) -0.000000 (zz)
    Force constants are written into "FORCE_CONSTANTS".
      Array shape of force constants is (16, 16, 3, 3).
    
    Reciprocal space paths in reduced coordinates:
    [ 0.000  0.000  0.000] --> [ 0.500  0.000  0.500]
    [ 0.500  0.000  0.500] --> [ 0.625  0.250  0.625]
    [ 0.375  0.375  0.750] --> [ 0.000  0.000  0.000]
    [ 0.000  0.000  0.000] --> [ 0.500  0.500  0.500]
    Mesh numbers: [8 8 8]
    Number of irreducible q-points on sampling mesh: 60/512
    Calculating phonons on sampling mesh...
    Calculating DOS...
    
    Summary of calculation was written in "phonopy.yaml".
                     _
       ___ _ __   __| |
      / _ \ '_ \ / _` |
     |  __/ | | | (_| |
      \___|_| |_|\__,_|
    
    

这一步结束之后，Phonopy软件会产生`band.yaml`（用于绘制声子谱）和`FORCE_CONSTANTS`文件。其中，`FORCE_CONSTANTS`文件包含的数据即为二阶力常数，注意这里务必设置FULL_FORCE_CONSTANTS = .TRUE.，输出全部的二阶力常数，否则ShengBTE读取数据会报错。

此外，可以使用如下命令输出gnuplot格式的声子谱，用于绘制声子谱：


```bash
phonopy-bandplot --gnuplot > pho.dat
```

    FORCE_CONSTANTS  STRU                         au2si.py   phonopy.yaml
    FORCE_SETS       STRU-001                     band.conf  phonopy_disp.yaml
    INPUT            STRU.in                      band.yaml  setting.conf
    KPT              Si_ONCV_PBE-1.0.upf          mesh.yaml  total_dos.dat
    [0m[01;34mOUT.DIA-50[0m       Si_gga_7au_100Ry_2s2p1d.orb  pho.dat
    

#### 1.3 后处理

注意ShengBTE软件要求`FORCE_CONSTANTS_2ND`文件里数据的单位为 eV/Å^2，但是ABACUS结合phonopy计算的FORCE_CONSTANTS单位为 $eV/(Å*au)$，其中$au$是原子单位制，$1 au = 0.52918 Å$。可以使用`2nd`目录下提供的au2si.py脚本进行单位转换，生成`FORCE_CONSTANTS_2ND`文件，命令如下：


```bash
python au2si.py
```

在 `shengbte` 文件夹中提供了`FORCE_CONSTANTS_2ND` 文件供参考计算结果。

### 2. 计算三阶力常数

要计算三阶力常数，需要结合 thirdorder 程序，计算后输出三阶力常数文件`FORCE_CONSTANTS_3RD`。但是，thirdorder目前只支持读取VASP和QE的输入输出文件。因此，这里我们是通过将ABACUS的结构文件和输出受力分别转换为`POSCAR` 和 `vasprun.xml`来使用 thirdorder，请先进入`3rd`文件夹，具体步骤将在以下叙述。

#### 2.1 获得微扰构型

首先将 ABACUS 软件进行结构优化（relax）后的`STRU`文件转化为`POSCAR`（目录下已给出转化过的`POSCAR`，或者需要自己动手进行这个转换）。

之后，运行thirdorder_vasp程序，产生微扰过后的一系列原子构型文件`3RD.POSCAR.*`，例如这个例子一共产生了40个构型：

备注：thirdorder 程序在镜像`py37`的环境中提前编译好，编译教程可参考[链接](https://www.shengbte.org/documentation)


```bash
cd /abacus-user-guide/examples/interface_ShengBTE/LCAO/3rd 
/opt/mamba/envs/py37/bin/python /thirdorder/thirdorder_vasp.py sow 2 2 2 -2
```

    Reading POSCAR
    Analyzing the symmetries
    - Symmetry group b'Fd-3m' detected
    - 48 symmetry operations
    Creating the supercell
    Computing all distances in the supercell
    - Automatic cutoff: 0.43260363256848766 nm
    Looking for an irreducible set of third-order IFCs
    - 5 triplet equivalence classes found
    - 40 DFT runs are needed
    
    .d88888b   .88888.  dP   dP   dP
    88.    "' d8'   `8b 88   88   88
    `Y88888b. 88     88 88  .8P  .8P
          `8b 88     88 88  d8'  d8'
    d8'   .8P Y8.   .8P 88.d8P8.d8P
     Y88888P   `8888P'  8888' Y88'
    ooooooooooooooooooooooooooooooooo
    
    Writing undisplaced coordinates to 3RD.SPOSCAR
    Writing displaced coordinates to 3RD.POSCAR.*
    
    888888ba   .88888.  888888ba   88888888b
    88    `8b d8'   `8b 88    `8b  88
    88     88 88     88 88     88 a88aaaa
    88     88 88     88 88     88  88
    88    .8P Y8.   .8P 88     88  88
    8888888P   `8888P'  dP     dP  88888888P
    ooooooooooooooooooooooooooooooooooooooooo
    
    

运行 pos2stru.py，将上述`POSCAR`转化为`STRU`文件，注意该脚本里调用了ASE软件包的函数（需提前安装好ASE）：


```bash
python pos2stru.py
```

    3RD.POSCAR.32 32
    STRU_32
    3RD.POSCAR.37 37
    STRU_37
    3RD.POSCAR.22 22
    STRU_22
    3RD.POSCAR.20 20
    STRU_20
    3RD.POSCAR.33 33
    STRU_33
    3RD.POSCAR.27 27
    STRU_27
    3RD.POSCAR.38 38
    STRU_38
    3RD.POSCAR.26 26
    STRU_26
    3RD.POSCAR.03 03
    STRU_03
    3RD.POSCAR.08 08
    STRU_08
    3RD.POSCAR.09 09
    STRU_09
    3RD.POSCAR.29 29
    STRU_29
    3RD.POSCAR.14 14
    STRU_14
    3RD.POSCAR.15 15
    STRU_15
    3RD.POSCAR.25 25
    STRU_25
    3RD.POSCAR.23 23
    STRU_23
    3RD.POSCAR.07 07
    STRU_07
    3RD.POSCAR.35 35
    STRU_35
    3RD.POSCAR.11 11
    STRU_11
    3RD.POSCAR.24 24
    STRU_24
    3RD.POSCAR.39 39
    STRU_39
    3RD.POSCAR.21 21
    STRU_21
    3RD.POSCAR.19 19
    STRU_19
    3RD.POSCAR.18 18
    STRU_18
    3RD.POSCAR.30 30
    STRU_30
    3RD.POSCAR.13 13
    STRU_13
    3RD.POSCAR.02 02
    STRU_02
    3RD.POSCAR.31 31
    STRU_31
    3RD.POSCAR.05 05
    STRU_05
    3RD.POSCAR.06 06
    STRU_06
    3RD.POSCAR.01 01
    STRU_01
    3RD.POSCAR.16 16
    STRU_16
    3RD.POSCAR.36 36
    STRU_36
    3RD.POSCAR.10 10
    STRU_10
    3RD.POSCAR.34 34
    STRU_34
    3RD.POSCAR.40 40
    STRU_40
    3RD.POSCAR.12 12
    STRU_12
    3RD.POSCAR.04 04
    STRU_04
    3RD.POSCAR.28 28
    STRU_28
    3RD.POSCAR.17 17
    STRU_17
    

注意：这里**不能调用 dpdata 软件进行转化**。因为 dpdata 会强制将晶格改为下三角矩阵，相当于旋转了晶格，会导致原子间受力方向也相应旋转，从而发生错误。

#### 2.2 计算微扰构型的原子受力

这里需要采用ABACUS对40个原子构型分别进行SCF计算，会有些耗时。建议每个SCF单独在`SCF-*`文件夹内运行，这里的`INPUT`中的scf_thr 需要至少小到1e-8才能得到收敛的结果。

这里为了更好演示，这里我们采用**notebook的新功能数据集**，提前将40个SCF的计算结果文件和提交脚本打包在`/bohr/abacus-shengbte-v11u/v1`下，供大家参考和学习。



```bash
cp -r /bohr/abacus-shengbte-v11u/v1/  .  &&  cd v1 && tar zxvf scf.tgz  &&  cp batch.sh ..
```

    SCF-01/
    SCF-01/.completed
    SCF-01/.lbg-11490-8462671_base.sh
    SCF-01/.lbg-compute-completed
    SCF-01/INPUT
    SCF-01/KPT
    SCF-01/OUT.DIA-50-01/
    SCF-01/OUT.DIA-50-01/INPUT
    SCF-01/OUT.DIA-50-01/STRU_READIN_ADJUST.cif
    SCF-01/OUT.DIA-50-01/STRU_SIMPLE.cif
    SCF-01/OUT.DIA-50-01/istate.info
    SCF-01/OUT.DIA-50-01/kpoints
    SCF-01/OUT.DIA-50-01/running_scf.log
    SCF-01/OUT.DIA-50-01/warning.log
    SCF-01/STDOUTERR
    SCF-01/STRU_01
    SCF-01/Si_ONCV_PBE-1.0.upf
    SCF-01/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-01/job.json
    SCF-01/lbg-11490-8462671.sh
    SCF-01/log
    SCF-02/
    SCF-02/.completed
    SCF-02/.lbg-11490-8462672_base.sh
    SCF-02/.lbg-compute-completed
    SCF-02/INPUT
    SCF-02/KPT
    SCF-02/OUT.DIA-50-02/
    SCF-02/OUT.DIA-50-02/INPUT
    SCF-02/OUT.DIA-50-02/STRU_READIN_ADJUST.cif
    SCF-02/OUT.DIA-50-02/STRU_SIMPLE.cif
    SCF-02/OUT.DIA-50-02/istate.info
    SCF-02/OUT.DIA-50-02/kpoints
    SCF-02/OUT.DIA-50-02/running_scf.log
    SCF-02/OUT.DIA-50-02/warning.log
    SCF-02/STDOUTERR
    SCF-02/STRU_02
    SCF-02/Si_ONCV_PBE-1.0.upf
    SCF-02/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-02/job.json
    SCF-02/lbg-11490-8462672.sh
    SCF-02/log
    SCF-03/
    SCF-03/.completed
    SCF-03/.lbg-11490-8462673_base.sh
    SCF-03/.lbg-compute-completed
    SCF-03/INPUT
    SCF-03/KPT
    SCF-03/OUT.DIA-50-03/
    SCF-03/OUT.DIA-50-03/INPUT
    SCF-03/OUT.DIA-50-03/STRU_READIN_ADJUST.cif
    SCF-03/OUT.DIA-50-03/STRU_SIMPLE.cif
    SCF-03/OUT.DIA-50-03/istate.info
    SCF-03/OUT.DIA-50-03/kpoints
    SCF-03/OUT.DIA-50-03/running_scf.log
    SCF-03/OUT.DIA-50-03/warning.log
    SCF-03/STDOUTERR
    SCF-03/STRU_03
    SCF-03/Si_ONCV_PBE-1.0.upf
    SCF-03/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-03/job.json
    SCF-03/lbg-11490-8462673.sh
    SCF-03/log
    SCF-04/
    SCF-04/.completed
    SCF-04/.lbg-11490-8462674_base.sh
    SCF-04/.lbg-compute-completed
    SCF-04/INPUT
    SCF-04/KPT
    SCF-04/OUT.DIA-50-04/
    SCF-04/OUT.DIA-50-04/INPUT
    SCF-04/OUT.DIA-50-04/STRU_READIN_ADJUST.cif
    SCF-04/OUT.DIA-50-04/STRU_SIMPLE.cif
    SCF-04/OUT.DIA-50-04/istate.info
    SCF-04/OUT.DIA-50-04/kpoints
    SCF-04/OUT.DIA-50-04/running_scf.log
    SCF-04/OUT.DIA-50-04/warning.log
    SCF-04/STDOUTERR
    SCF-04/STRU_04
    SCF-04/Si_ONCV_PBE-1.0.upf
    SCF-04/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-04/job.json
    SCF-04/lbg-11490-8462674.sh
    SCF-04/log
    SCF-05/
    SCF-05/.completed
    SCF-05/.lbg-11490-8462675_base.sh
    SCF-05/.lbg-compute-completed
    SCF-05/INPUT
    SCF-05/KPT
    SCF-05/OUT.DIA-50-05/
    SCF-05/OUT.DIA-50-05/INPUT
    SCF-05/OUT.DIA-50-05/STRU_READIN_ADJUST.cif
    SCF-05/OUT.DIA-50-05/STRU_SIMPLE.cif
    SCF-05/OUT.DIA-50-05/istate.info
    SCF-05/OUT.DIA-50-05/kpoints
    SCF-05/OUT.DIA-50-05/running_scf.log
    SCF-05/OUT.DIA-50-05/warning.log
    SCF-05/STDOUTERR
    SCF-05/STRU_05
    SCF-05/Si_ONCV_PBE-1.0.upf
    SCF-05/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-05/job.json
    SCF-05/lbg-11490-8462675.sh
    SCF-05/log
    SCF-06/
    SCF-06/.completed
    SCF-06/.lbg-11490-8462676_base.sh
    SCF-06/.lbg-compute-completed
    SCF-06/INPUT
    SCF-06/KPT
    SCF-06/OUT.DIA-50-06/
    SCF-06/OUT.DIA-50-06/INPUT
    SCF-06/OUT.DIA-50-06/STRU_READIN_ADJUST.cif
    SCF-06/OUT.DIA-50-06/STRU_SIMPLE.cif
    SCF-06/OUT.DIA-50-06/istate.info
    SCF-06/OUT.DIA-50-06/kpoints
    SCF-06/OUT.DIA-50-06/running_scf.log
    SCF-06/OUT.DIA-50-06/warning.log
    SCF-06/STDOUTERR
    SCF-06/STRU_06
    SCF-06/Si_ONCV_PBE-1.0.upf
    SCF-06/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-06/job.json
    SCF-06/lbg-11490-8462676.sh
    SCF-06/log
    SCF-07/
    SCF-07/.completed
    SCF-07/.lbg-11490-8462677_base.sh
    SCF-07/.lbg-compute-completed
    SCF-07/INPUT
    SCF-07/KPT
    SCF-07/OUT.DIA-50-07/
    SCF-07/OUT.DIA-50-07/INPUT
    SCF-07/OUT.DIA-50-07/STRU_READIN_ADJUST.cif
    SCF-07/OUT.DIA-50-07/STRU_SIMPLE.cif
    SCF-07/OUT.DIA-50-07/istate.info
    SCF-07/OUT.DIA-50-07/kpoints
    SCF-07/OUT.DIA-50-07/running_scf.log
    SCF-07/OUT.DIA-50-07/warning.log
    SCF-07/STDOUTERR
    SCF-07/STRU_07
    SCF-07/Si_ONCV_PBE-1.0.upf
    SCF-07/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-07/job.json
    SCF-07/lbg-11490-8462677.sh
    SCF-07/log
    SCF-08/
    SCF-08/.completed
    SCF-08/.lbg-11490-8462678_base.sh
    SCF-08/.lbg-compute-completed
    SCF-08/INPUT
    SCF-08/KPT
    SCF-08/OUT.DIA-50-08/
    SCF-08/OUT.DIA-50-08/INPUT
    SCF-08/OUT.DIA-50-08/STRU_READIN_ADJUST.cif
    SCF-08/OUT.DIA-50-08/STRU_SIMPLE.cif
    SCF-08/OUT.DIA-50-08/istate.info
    SCF-08/OUT.DIA-50-08/kpoints
    SCF-08/OUT.DIA-50-08/running_scf.log
    SCF-08/OUT.DIA-50-08/warning.log
    SCF-08/STDOUTERR
    SCF-08/STRU_08
    SCF-08/Si_ONCV_PBE-1.0.upf
    SCF-08/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-08/job.json
    SCF-08/lbg-11490-8462678.sh
    SCF-08/log
    SCF-09/
    SCF-09/.completed
    SCF-09/.lbg-11490-8462679_base.sh
    SCF-09/.lbg-compute-completed
    SCF-09/INPUT
    SCF-09/KPT
    SCF-09/OUT.DIA-50-09/
    SCF-09/OUT.DIA-50-09/INPUT
    SCF-09/OUT.DIA-50-09/STRU_READIN_ADJUST.cif
    SCF-09/OUT.DIA-50-09/STRU_SIMPLE.cif
    SCF-09/OUT.DIA-50-09/istate.info
    SCF-09/OUT.DIA-50-09/kpoints
    SCF-09/OUT.DIA-50-09/running_scf.log
    SCF-09/OUT.DIA-50-09/warning.log
    SCF-09/STDOUTERR
    SCF-09/STRU_09
    SCF-09/Si_ONCV_PBE-1.0.upf
    SCF-09/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-09/job.json
    SCF-09/lbg-11490-8462679.sh
    SCF-09/log
    SCF-10/
    SCF-10/.completed
    SCF-10/.lbg-11490-8462680_base.sh
    SCF-10/.lbg-compute-completed
    SCF-10/INPUT
    SCF-10/KPT
    SCF-10/OUT.DIA-50-10/
    SCF-10/OUT.DIA-50-10/INPUT
    SCF-10/OUT.DIA-50-10/STRU_READIN_ADJUST.cif
    SCF-10/OUT.DIA-50-10/STRU_SIMPLE.cif
    SCF-10/OUT.DIA-50-10/istate.info
    SCF-10/OUT.DIA-50-10/kpoints
    SCF-10/OUT.DIA-50-10/running_scf.log
    SCF-10/OUT.DIA-50-10/warning.log
    SCF-10/STDOUTERR
    SCF-10/STRU_10
    SCF-10/Si_ONCV_PBE-1.0.upf
    SCF-10/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-10/job.json
    SCF-10/lbg-11490-8462680.sh
    SCF-10/log
    SCF-11/
    SCF-11/.completed
    SCF-11/.lbg-11490-8462681_base.sh
    SCF-11/.lbg-compute-completed
    SCF-11/INPUT
    SCF-11/KPT
    SCF-11/OUT.DIA-50-11/
    SCF-11/OUT.DIA-50-11/INPUT
    SCF-11/OUT.DIA-50-11/STRU_READIN_ADJUST.cif
    SCF-11/OUT.DIA-50-11/STRU_SIMPLE.cif
    SCF-11/OUT.DIA-50-11/istate.info
    SCF-11/OUT.DIA-50-11/kpoints
    SCF-11/OUT.DIA-50-11/running_scf.log
    SCF-11/OUT.DIA-50-11/warning.log
    SCF-11/STDOUTERR
    SCF-11/STRU_11
    SCF-11/Si_ONCV_PBE-1.0.upf
    SCF-11/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-11/job.json
    SCF-11/lbg-11490-8462681.sh
    SCF-11/log
    SCF-12/
    SCF-12/.completed
    SCF-12/.lbg-11490-8462682_base.sh
    SCF-12/.lbg-compute-completed
    SCF-12/INPUT
    SCF-12/KPT
    SCF-12/OUT.DIA-50-12/
    SCF-12/OUT.DIA-50-12/INPUT
    SCF-12/OUT.DIA-50-12/STRU_READIN_ADJUST.cif
    SCF-12/OUT.DIA-50-12/STRU_SIMPLE.cif
    SCF-12/OUT.DIA-50-12/istate.info
    SCF-12/OUT.DIA-50-12/kpoints
    SCF-12/OUT.DIA-50-12/running_scf.log
    SCF-12/OUT.DIA-50-12/warning.log
    SCF-12/STDOUTERR
    SCF-12/STRU_12
    SCF-12/Si_ONCV_PBE-1.0.upf
    SCF-12/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-12/job.json
    SCF-12/lbg-11490-8462682.sh
    SCF-12/log
    SCF-13/
    SCF-13/.completed
    SCF-13/.lbg-11490-8462683_base.sh
    SCF-13/.lbg-compute-completed
    SCF-13/INPUT
    SCF-13/KPT
    SCF-13/OUT.DIA-50-13/
    SCF-13/OUT.DIA-50-13/INPUT
    SCF-13/OUT.DIA-50-13/STRU_READIN_ADJUST.cif
    SCF-13/OUT.DIA-50-13/STRU_SIMPLE.cif
    SCF-13/OUT.DIA-50-13/istate.info
    SCF-13/OUT.DIA-50-13/kpoints
    SCF-13/OUT.DIA-50-13/running_scf.log
    SCF-13/OUT.DIA-50-13/warning.log
    SCF-13/STDOUTERR
    SCF-13/STRU_13
    SCF-13/Si_ONCV_PBE-1.0.upf
    SCF-13/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-13/job.json
    SCF-13/lbg-11490-8462683.sh
    SCF-13/log
    SCF-14/
    SCF-14/.completed
    SCF-14/.lbg-11490-8462684_base.sh
    SCF-14/.lbg-compute-completed
    SCF-14/INPUT
    SCF-14/KPT
    SCF-14/OUT.DIA-50-14/
    SCF-14/OUT.DIA-50-14/INPUT
    SCF-14/OUT.DIA-50-14/STRU_READIN_ADJUST.cif
    SCF-14/OUT.DIA-50-14/STRU_SIMPLE.cif
    SCF-14/OUT.DIA-50-14/istate.info
    SCF-14/OUT.DIA-50-14/kpoints
    SCF-14/OUT.DIA-50-14/running_scf.log
    SCF-14/OUT.DIA-50-14/warning.log
    SCF-14/STDOUTERR
    SCF-14/STRU_14
    SCF-14/Si_ONCV_PBE-1.0.upf
    SCF-14/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-14/job.json
    SCF-14/lbg-11490-8462684.sh
    SCF-14/log
    SCF-15/
    SCF-15/.completed
    SCF-15/.lbg-11490-8462685_base.sh
    SCF-15/.lbg-compute-completed
    SCF-15/INPUT
    SCF-15/KPT
    SCF-15/OUT.DIA-50-15/
    SCF-15/OUT.DIA-50-15/INPUT
    SCF-15/OUT.DIA-50-15/STRU_READIN_ADJUST.cif
    SCF-15/OUT.DIA-50-15/STRU_SIMPLE.cif
    SCF-15/OUT.DIA-50-15/istate.info
    SCF-15/OUT.DIA-50-15/kpoints
    SCF-15/OUT.DIA-50-15/running_scf.log
    SCF-15/OUT.DIA-50-15/warning.log
    SCF-15/STDOUTERR
    SCF-15/STRU_15
    SCF-15/Si_ONCV_PBE-1.0.upf
    SCF-15/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-15/job.json
    SCF-15/lbg-11490-8462685.sh
    SCF-15/log
    SCF-16/
    SCF-16/.completed
    SCF-16/.lbg-11490-8462686_base.sh
    SCF-16/.lbg-compute-completed
    SCF-16/INPUT
    SCF-16/KPT
    SCF-16/OUT.DIA-50-16/
    SCF-16/OUT.DIA-50-16/INPUT
    SCF-16/OUT.DIA-50-16/STRU_READIN_ADJUST.cif
    SCF-16/OUT.DIA-50-16/STRU_SIMPLE.cif
    SCF-16/OUT.DIA-50-16/istate.info
    SCF-16/OUT.DIA-50-16/kpoints
    SCF-16/OUT.DIA-50-16/running_scf.log
    SCF-16/OUT.DIA-50-16/warning.log
    SCF-16/STDOUTERR
    SCF-16/STRU_16
    SCF-16/Si_ONCV_PBE-1.0.upf
    SCF-16/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-16/job.json
    SCF-16/lbg-11490-8462686.sh
    SCF-16/log
    SCF-17/
    SCF-17/.completed
    SCF-17/.lbg-11490-8462687_base.sh
    SCF-17/.lbg-compute-completed
    SCF-17/INPUT
    SCF-17/KPT
    SCF-17/OUT.DIA-50-17/
    SCF-17/OUT.DIA-50-17/INPUT
    SCF-17/OUT.DIA-50-17/STRU_READIN_ADJUST.cif
    SCF-17/OUT.DIA-50-17/STRU_SIMPLE.cif
    SCF-17/OUT.DIA-50-17/istate.info
    SCF-17/OUT.DIA-50-17/kpoints
    SCF-17/OUT.DIA-50-17/running_scf.log
    SCF-17/OUT.DIA-50-17/warning.log
    SCF-17/STDOUTERR
    SCF-17/STRU_17
    SCF-17/Si_ONCV_PBE-1.0.upf
    SCF-17/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-17/job.json
    SCF-17/lbg-11490-8462687.sh
    SCF-17/log
    SCF-18/
    SCF-18/.completed
    SCF-18/.lbg-11490-8462688_base.sh
    SCF-18/.lbg-compute-completed
    SCF-18/INPUT
    SCF-18/KPT
    SCF-18/OUT.DIA-50-18/
    SCF-18/OUT.DIA-50-18/INPUT
    SCF-18/OUT.DIA-50-18/STRU_READIN_ADJUST.cif
    SCF-18/OUT.DIA-50-18/STRU_SIMPLE.cif
    SCF-18/OUT.DIA-50-18/istate.info
    SCF-18/OUT.DIA-50-18/kpoints
    SCF-18/OUT.DIA-50-18/running_scf.log
    SCF-18/OUT.DIA-50-18/warning.log
    SCF-18/STDOUTERR
    SCF-18/STRU_18
    SCF-18/Si_ONCV_PBE-1.0.upf
    SCF-18/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-18/job.json
    SCF-18/lbg-11490-8462688.sh
    SCF-18/log
    SCF-19/
    SCF-19/.completed
    SCF-19/.lbg-11490-8462689_base.sh
    SCF-19/.lbg-compute-completed
    SCF-19/INPUT
    SCF-19/KPT
    SCF-19/OUT.DIA-50-19/
    SCF-19/OUT.DIA-50-19/INPUT
    SCF-19/OUT.DIA-50-19/STRU_READIN_ADJUST.cif
    SCF-19/OUT.DIA-50-19/STRU_SIMPLE.cif
    SCF-19/OUT.DIA-50-19/istate.info
    SCF-19/OUT.DIA-50-19/kpoints
    SCF-19/OUT.DIA-50-19/running_scf.log
    SCF-19/OUT.DIA-50-19/warning.log
    SCF-19/STDOUTERR
    SCF-19/STRU_19
    SCF-19/Si_ONCV_PBE-1.0.upf
    SCF-19/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-19/job.json
    SCF-19/lbg-11490-8462689.sh
    SCF-19/log
    SCF-20/
    SCF-20/.completed
    SCF-20/.lbg-11490-8462690_base.sh
    SCF-20/.lbg-compute-completed
    SCF-20/INPUT
    SCF-20/KPT
    SCF-20/OUT.DIA-50-20/
    SCF-20/OUT.DIA-50-20/INPUT
    SCF-20/OUT.DIA-50-20/STRU_READIN_ADJUST.cif
    SCF-20/OUT.DIA-50-20/STRU_SIMPLE.cif
    SCF-20/OUT.DIA-50-20/istate.info
    SCF-20/OUT.DIA-50-20/kpoints
    SCF-20/OUT.DIA-50-20/running_scf.log
    SCF-20/OUT.DIA-50-20/warning.log
    SCF-20/STDOUTERR
    SCF-20/STRU_20
    SCF-20/Si_ONCV_PBE-1.0.upf
    SCF-20/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-20/job.json
    SCF-20/lbg-11490-8462690.sh
    SCF-20/log
    SCF-21/
    SCF-21/.completed
    SCF-21/.lbg-11490-8462691_base.sh
    SCF-21/.lbg-compute-completed
    SCF-21/INPUT
    SCF-21/KPT
    SCF-21/OUT.DIA-50-21/
    SCF-21/OUT.DIA-50-21/INPUT
    SCF-21/OUT.DIA-50-21/STRU_READIN_ADJUST.cif
    SCF-21/OUT.DIA-50-21/STRU_SIMPLE.cif
    SCF-21/OUT.DIA-50-21/istate.info
    SCF-21/OUT.DIA-50-21/kpoints
    SCF-21/OUT.DIA-50-21/running_scf.log
    SCF-21/OUT.DIA-50-21/warning.log
    SCF-21/STDOUTERR
    SCF-21/STRU_21
    SCF-21/Si_ONCV_PBE-1.0.upf
    SCF-21/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-21/job.json
    SCF-21/lbg-11490-8462691.sh
    SCF-21/log
    SCF-22/
    SCF-22/.completed
    SCF-22/.lbg-11490-8462692_base.sh
    SCF-22/.lbg-compute-completed
    SCF-22/INPUT
    SCF-22/KPT
    SCF-22/OUT.DIA-50-22/
    SCF-22/OUT.DIA-50-22/INPUT
    SCF-22/OUT.DIA-50-22/STRU_READIN_ADJUST.cif
    SCF-22/OUT.DIA-50-22/STRU_SIMPLE.cif
    SCF-22/OUT.DIA-50-22/istate.info
    SCF-22/OUT.DIA-50-22/kpoints
    SCF-22/OUT.DIA-50-22/running_scf.log
    SCF-22/OUT.DIA-50-22/warning.log
    SCF-22/STDOUTERR
    SCF-22/STRU_22
    SCF-22/Si_ONCV_PBE-1.0.upf
    SCF-22/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-22/job.json
    SCF-22/lbg-11490-8462692.sh
    SCF-22/log
    SCF-23/
    SCF-23/.completed
    SCF-23/.lbg-11490-8462693_base.sh
    SCF-23/.lbg-compute-completed
    SCF-23/INPUT
    SCF-23/KPT
    SCF-23/OUT.DIA-50-23/
    SCF-23/OUT.DIA-50-23/INPUT
    SCF-23/OUT.DIA-50-23/STRU_READIN_ADJUST.cif
    SCF-23/OUT.DIA-50-23/STRU_SIMPLE.cif
    SCF-23/OUT.DIA-50-23/istate.info
    SCF-23/OUT.DIA-50-23/kpoints
    SCF-23/OUT.DIA-50-23/running_scf.log
    SCF-23/OUT.DIA-50-23/warning.log
    SCF-23/STDOUTERR
    SCF-23/STRU_23
    SCF-23/Si_ONCV_PBE-1.0.upf
    SCF-23/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-23/job.json
    SCF-23/lbg-11490-8462693.sh
    SCF-23/log
    SCF-24/
    SCF-24/.completed
    SCF-24/.lbg-11490-8462694_base.sh
    SCF-24/.lbg-compute-completed
    SCF-24/INPUT
    SCF-24/KPT
    SCF-24/OUT.DIA-50-24/
    SCF-24/OUT.DIA-50-24/INPUT
    SCF-24/OUT.DIA-50-24/STRU_READIN_ADJUST.cif
    SCF-24/OUT.DIA-50-24/STRU_SIMPLE.cif
    SCF-24/OUT.DIA-50-24/istate.info
    SCF-24/OUT.DIA-50-24/kpoints
    SCF-24/OUT.DIA-50-24/running_scf.log
    SCF-24/OUT.DIA-50-24/warning.log
    SCF-24/STDOUTERR
    SCF-24/STRU_24
    SCF-24/Si_ONCV_PBE-1.0.upf
    SCF-24/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-24/job.json
    SCF-24/lbg-11490-8462694.sh
    SCF-24/log
    SCF-25/
    SCF-25/.completed
    SCF-25/.lbg-11490-8462695_base.sh
    SCF-25/.lbg-compute-completed
    SCF-25/INPUT
    SCF-25/KPT
    SCF-25/OUT.DIA-50-25/
    SCF-25/OUT.DIA-50-25/INPUT
    SCF-25/OUT.DIA-50-25/STRU_READIN_ADJUST.cif
    SCF-25/OUT.DIA-50-25/STRU_SIMPLE.cif
    SCF-25/OUT.DIA-50-25/istate.info
    SCF-25/OUT.DIA-50-25/kpoints
    SCF-25/OUT.DIA-50-25/running_scf.log
    SCF-25/OUT.DIA-50-25/warning.log
    SCF-25/STDOUTERR
    SCF-25/STRU_25
    SCF-25/Si_ONCV_PBE-1.0.upf
    SCF-25/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-25/job.json
    SCF-25/lbg-11490-8462695.sh
    SCF-25/log
    SCF-26/
    SCF-26/.completed
    SCF-26/.lbg-11490-8462696_base.sh
    SCF-26/.lbg-compute-completed
    SCF-26/INPUT
    SCF-26/KPT
    SCF-26/OUT.DIA-50-26/
    SCF-26/OUT.DIA-50-26/INPUT
    SCF-26/OUT.DIA-50-26/STRU_READIN_ADJUST.cif
    SCF-26/OUT.DIA-50-26/STRU_SIMPLE.cif
    SCF-26/OUT.DIA-50-26/istate.info
    SCF-26/OUT.DIA-50-26/kpoints
    SCF-26/OUT.DIA-50-26/running_scf.log
    SCF-26/OUT.DIA-50-26/warning.log
    SCF-26/STDOUTERR
    SCF-26/STRU_26
    SCF-26/Si_ONCV_PBE-1.0.upf
    SCF-26/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-26/job.json
    SCF-26/lbg-11490-8462696.sh
    SCF-26/log
    SCF-27/
    SCF-27/.completed
    SCF-27/.lbg-11490-8462697_base.sh
    SCF-27/.lbg-compute-completed
    SCF-27/INPUT
    SCF-27/KPT
    SCF-27/OUT.DIA-50-27/
    SCF-27/OUT.DIA-50-27/INPUT
    SCF-27/OUT.DIA-50-27/STRU_READIN_ADJUST.cif
    SCF-27/OUT.DIA-50-27/STRU_SIMPLE.cif
    SCF-27/OUT.DIA-50-27/istate.info
    SCF-27/OUT.DIA-50-27/kpoints
    SCF-27/OUT.DIA-50-27/running_scf.log
    SCF-27/OUT.DIA-50-27/warning.log
    SCF-27/STDOUTERR
    SCF-27/STRU_27
    SCF-27/Si_ONCV_PBE-1.0.upf
    SCF-27/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-27/job.json
    SCF-27/lbg-11490-8462697.sh
    SCF-27/log
    SCF-28/
    SCF-28/.completed
    SCF-28/.lbg-11490-8462698_base.sh
    SCF-28/.lbg-compute-completed
    SCF-28/INPUT
    SCF-28/KPT
    SCF-28/OUT.DIA-50-28/
    SCF-28/OUT.DIA-50-28/INPUT
    SCF-28/OUT.DIA-50-28/STRU_READIN_ADJUST.cif
    SCF-28/OUT.DIA-50-28/STRU_SIMPLE.cif
    SCF-28/OUT.DIA-50-28/istate.info
    SCF-28/OUT.DIA-50-28/kpoints
    SCF-28/OUT.DIA-50-28/running_scf.log
    SCF-28/OUT.DIA-50-28/warning.log
    SCF-28/STDOUTERR
    SCF-28/STRU_28
    SCF-28/Si_ONCV_PBE-1.0.upf
    SCF-28/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-28/job.json
    SCF-28/lbg-11490-8462698.sh
    SCF-28/log
    SCF-29/
    SCF-29/.completed
    SCF-29/.lbg-11490-8462699_base.sh
    SCF-29/.lbg-compute-completed
    SCF-29/INPUT
    SCF-29/KPT
    SCF-29/OUT.DIA-50-29/
    SCF-29/OUT.DIA-50-29/INPUT
    SCF-29/OUT.DIA-50-29/STRU_READIN_ADJUST.cif
    SCF-29/OUT.DIA-50-29/STRU_SIMPLE.cif
    SCF-29/OUT.DIA-50-29/istate.info
    SCF-29/OUT.DIA-50-29/kpoints
    SCF-29/OUT.DIA-50-29/running_scf.log
    SCF-29/OUT.DIA-50-29/warning.log
    SCF-29/STDOUTERR
    SCF-29/STRU_29
    SCF-29/Si_ONCV_PBE-1.0.upf
    SCF-29/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-29/job.json
    SCF-29/lbg-11490-8462699.sh
    SCF-29/log
    SCF-30/
    SCF-30/.completed
    SCF-30/.lbg-11490-8462701_base.sh
    SCF-30/.lbg-compute-completed
    SCF-30/INPUT
    SCF-30/KPT
    SCF-30/OUT.DIA-50-30/
    SCF-30/OUT.DIA-50-30/INPUT
    SCF-30/OUT.DIA-50-30/STRU_READIN_ADJUST.cif
    SCF-30/OUT.DIA-50-30/STRU_SIMPLE.cif
    SCF-30/OUT.DIA-50-30/istate.info
    SCF-30/OUT.DIA-50-30/kpoints
    SCF-30/OUT.DIA-50-30/running_scf.log
    SCF-30/OUT.DIA-50-30/warning.log
    SCF-30/STDOUTERR
    SCF-30/STRU_30
    SCF-30/Si_ONCV_PBE-1.0.upf
    SCF-30/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-30/job.json
    SCF-30/lbg-11490-8462701.sh
    SCF-30/log
    SCF-31/
    SCF-31/.completed
    SCF-31/.lbg-11490-8462702_base.sh
    SCF-31/.lbg-compute-completed
    SCF-31/INPUT
    SCF-31/KPT
    SCF-31/OUT.DIA-50-31/
    SCF-31/OUT.DIA-50-31/INPUT
    SCF-31/OUT.DIA-50-31/STRU_READIN_ADJUST.cif
    SCF-31/OUT.DIA-50-31/STRU_SIMPLE.cif
    SCF-31/OUT.DIA-50-31/istate.info
    SCF-31/OUT.DIA-50-31/kpoints
    SCF-31/OUT.DIA-50-31/running_scf.log
    SCF-31/OUT.DIA-50-31/warning.log
    SCF-31/STDOUTERR
    SCF-31/STRU_31
    SCF-31/Si_ONCV_PBE-1.0.upf
    SCF-31/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-31/job.json
    SCF-31/lbg-11490-8462702.sh
    SCF-31/log
    SCF-32/
    SCF-32/.completed
    SCF-32/.lbg-11490-8462704_base.sh
    SCF-32/.lbg-compute-completed
    SCF-32/INPUT
    SCF-32/KPT
    SCF-32/OUT.DIA-50-32/
    SCF-32/OUT.DIA-50-32/INPUT
    SCF-32/OUT.DIA-50-32/STRU_READIN_ADJUST.cif
    SCF-32/OUT.DIA-50-32/STRU_SIMPLE.cif
    SCF-32/OUT.DIA-50-32/istate.info
    SCF-32/OUT.DIA-50-32/kpoints
    SCF-32/OUT.DIA-50-32/running_scf.log
    SCF-32/OUT.DIA-50-32/warning.log
    SCF-32/STDOUTERR
    SCF-32/STRU_32
    SCF-32/Si_ONCV_PBE-1.0.upf
    SCF-32/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-32/job.json
    SCF-32/lbg-11490-8462704.sh
    SCF-32/log
    SCF-33/
    SCF-33/.completed
    SCF-33/.lbg-11490-8462705_base.sh
    SCF-33/.lbg-compute-completed
    SCF-33/INPUT
    SCF-33/KPT
    SCF-33/OUT.DIA-50-33/
    SCF-33/OUT.DIA-50-33/INPUT
    SCF-33/OUT.DIA-50-33/STRU_READIN_ADJUST.cif
    SCF-33/OUT.DIA-50-33/STRU_SIMPLE.cif
    SCF-33/OUT.DIA-50-33/istate.info
    SCF-33/OUT.DIA-50-33/kpoints
    SCF-33/OUT.DIA-50-33/running_scf.log
    SCF-33/OUT.DIA-50-33/warning.log
    SCF-33/STDOUTERR
    SCF-33/STRU_33
    SCF-33/Si_ONCV_PBE-1.0.upf
    SCF-33/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-33/job.json
    SCF-33/lbg-11490-8462705.sh
    SCF-33/log
    SCF-34/
    SCF-34/.completed
    SCF-34/.lbg-11490-8462706_base.sh
    SCF-34/.lbg-compute-completed
    SCF-34/INPUT
    SCF-34/KPT
    SCF-34/OUT.DIA-50-34/
    SCF-34/OUT.DIA-50-34/INPUT
    SCF-34/OUT.DIA-50-34/STRU_READIN_ADJUST.cif
    SCF-34/OUT.DIA-50-34/STRU_SIMPLE.cif
    SCF-34/OUT.DIA-50-34/istate.info
    SCF-34/OUT.DIA-50-34/kpoints
    SCF-34/OUT.DIA-50-34/running_scf.log
    SCF-34/OUT.DIA-50-34/warning.log
    SCF-34/STDOUTERR
    SCF-34/STRU_34
    SCF-34/Si_ONCV_PBE-1.0.upf
    SCF-34/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-34/job.json
    SCF-34/lbg-11490-8462706.sh
    SCF-34/log
    SCF-35/
    SCF-35/.completed
    SCF-35/.lbg-11490-8462707_base.sh
    SCF-35/.lbg-compute-completed
    SCF-35/INPUT
    SCF-35/KPT
    SCF-35/OUT.DIA-50-35/
    SCF-35/OUT.DIA-50-35/INPUT
    SCF-35/OUT.DIA-50-35/STRU_READIN_ADJUST.cif
    SCF-35/OUT.DIA-50-35/STRU_SIMPLE.cif
    SCF-35/OUT.DIA-50-35/istate.info
    SCF-35/OUT.DIA-50-35/kpoints
    SCF-35/OUT.DIA-50-35/running_scf.log
    SCF-35/OUT.DIA-50-35/warning.log
    SCF-35/STDOUTERR
    SCF-35/STRU_35
    SCF-35/Si_ONCV_PBE-1.0.upf
    SCF-35/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-35/job.json
    SCF-35/lbg-11490-8462707.sh
    SCF-35/log
    SCF-36/
    SCF-36/.completed
    SCF-36/.lbg-11490-8462708_base.sh
    SCF-36/.lbg-compute-completed
    SCF-36/INPUT
    SCF-36/KPT
    SCF-36/OUT.DIA-50-36/
    SCF-36/OUT.DIA-50-36/INPUT
    SCF-36/OUT.DIA-50-36/STRU_READIN_ADJUST.cif
    SCF-36/OUT.DIA-50-36/STRU_SIMPLE.cif
    SCF-36/OUT.DIA-50-36/istate.info
    SCF-36/OUT.DIA-50-36/kpoints
    SCF-36/OUT.DIA-50-36/running_scf.log
    SCF-36/OUT.DIA-50-36/warning.log
    SCF-36/STDOUTERR
    SCF-36/STRU_36
    SCF-36/Si_ONCV_PBE-1.0.upf
    SCF-36/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-36/job.json
    SCF-36/lbg-11490-8462708.sh
    SCF-36/log
    SCF-37/
    SCF-37/.completed
    SCF-37/.lbg-11490-8462709_base.sh
    SCF-37/.lbg-compute-completed
    SCF-37/INPUT
    SCF-37/KPT
    SCF-37/OUT.DIA-50-37/
    SCF-37/OUT.DIA-50-37/INPUT
    SCF-37/OUT.DIA-50-37/STRU_READIN_ADJUST.cif
    SCF-37/OUT.DIA-50-37/STRU_SIMPLE.cif
    SCF-37/OUT.DIA-50-37/istate.info
    SCF-37/OUT.DIA-50-37/kpoints
    SCF-37/OUT.DIA-50-37/running_scf.log
    SCF-37/OUT.DIA-50-37/warning.log
    SCF-37/STDOUTERR
    SCF-37/STRU_37
    SCF-37/Si_ONCV_PBE-1.0.upf
    SCF-37/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-37/job.json
    SCF-37/lbg-11490-8462709.sh
    SCF-37/log
    SCF-38/
    SCF-38/.completed
    SCF-38/.lbg-11490-8462710_base.sh
    SCF-38/.lbg-compute-completed
    SCF-38/INPUT
    SCF-38/KPT
    SCF-38/OUT.DIA-50-38/
    SCF-38/OUT.DIA-50-38/INPUT
    SCF-38/OUT.DIA-50-38/STRU_READIN_ADJUST.cif
    SCF-38/OUT.DIA-50-38/STRU_SIMPLE.cif
    SCF-38/OUT.DIA-50-38/istate.info
    SCF-38/OUT.DIA-50-38/kpoints
    SCF-38/OUT.DIA-50-38/running_scf.log
    SCF-38/OUT.DIA-50-38/warning.log
    SCF-38/STDOUTERR
    SCF-38/STRU_38
    SCF-38/Si_ONCV_PBE-1.0.upf
    SCF-38/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-38/job.json
    SCF-38/lbg-11490-8462710.sh
    SCF-38/log
    SCF-39/
    SCF-39/.completed
    SCF-39/.lbg-11490-8462711_base.sh
    SCF-39/.lbg-compute-completed
    SCF-39/INPUT
    SCF-39/KPT
    SCF-39/OUT.DIA-50-39/
    SCF-39/OUT.DIA-50-39/INPUT
    SCF-39/OUT.DIA-50-39/STRU_READIN_ADJUST.cif
    SCF-39/OUT.DIA-50-39/STRU_SIMPLE.cif
    SCF-39/OUT.DIA-50-39/istate.info
    SCF-39/OUT.DIA-50-39/kpoints
    SCF-39/OUT.DIA-50-39/running_scf.log
    SCF-39/OUT.DIA-50-39/warning.log
    SCF-39/STDOUTERR
    SCF-39/STRU_39
    SCF-39/Si_ONCV_PBE-1.0.upf
    SCF-39/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-39/job.json
    SCF-39/lbg-11490-8462711.sh
    SCF-39/log
    SCF-40/
    SCF-40/.completed
    SCF-40/.lbg-11490-8462712_base.sh
    SCF-40/.lbg-compute-completed
    SCF-40/INPUT
    SCF-40/KPT
    SCF-40/OUT.DIA-50-40/
    SCF-40/OUT.DIA-50-40/INPUT
    SCF-40/OUT.DIA-50-40/STRU_READIN_ADJUST.cif
    SCF-40/OUT.DIA-50-40/STRU_SIMPLE.cif
    SCF-40/OUT.DIA-50-40/istate.info
    SCF-40/OUT.DIA-50-40/kpoints
    SCF-40/OUT.DIA-50-40/running_scf.log
    SCF-40/OUT.DIA-50-40/warning.log
    SCF-40/STDOUTERR
    SCF-40/STRU_40
    SCF-40/Si_ONCV_PBE-1.0.upf
    SCF-40/Si_gga_7au_100Ry_2s2p1d.orb
    SCF-40/job.json
    SCF-40/lbg-11490-8462712.sh
    SCF-40/log
    batch.sh
    

使用脚本这里可以利用脚本`batch.sh`批量提交40个SCF的任务，批量产生`SCF-*`文件夹并提交计算，**注意需要修改脚本id号** (将第54行`project_id`后面的16280改为用户自己的id)，其内容为：


```bash
cd .. && pwd && cat batch.sh
```

    /abacus-user-guide/examples/interface_ShengBTE/LCAO/3rd
    #!/bin/bash
    
    
    for i in `ls STRU*`
    do
    stru=$(echo $i|cut -d"_" -f2)
    echo "stru=$stru"
    mkdir SCF-$stru
    cd SCF-$stru
    pwd
    cat > INPUT <<EOF
    INPUT_PARAMETERS
    #Parameters     (General)
    suffix          DIA-50-$stru
    calculation     scf
    esolver_type    ksdft
    pseudo_dir      ./
    orbital_dir     ./
    nbands          45
    symmetry        1
    cal_force       1
    cal_stress      1
    
    #Parameters (Accuracy)
    ecutwfc         100
    scf_thr         1e-8
    scf_nmax        100
    basis_type      lcao
    ks_solver       genelpa
    gamma_only      0
    smearing_method gauss
    smearing_sigma  0.01
    mixing_type     pulay
    mixing_beta     0.7
    
    stru_file       STRU_$stru
    EOF
    cat > KPT <<EOF
    K_POINTS
    0
    Gamma
    2 2 2 0 0 0
    EOF
    cp ../STRU_$stru .
    cp ../Si_ONCV_PBE-1.0.upf .
    cp ../Si_gga_7au_100Ry_2s2p1d.orb .
    
    cat > job.json <<EOF
    {
        "job_name": "ABACUS scf_$stru",
        "command": "OMP_NUM_THREADS=1 mpirun -np 8 abacus > log",
        "log_file": "log",
        "backward_files": ["*"],
        "project_id": 16280,
        "platform": "ali",
        "job_type": "container",
        "machine_type": "c8_m8_cpu",
        "image_address": "registry.dp.tech/dptech/abacus:3.2.1"
    }
    EOF
    
    lbg job submit -i ./job.json -p ./ -r /data/SCF-$stru
    #mpirun -n 96 ABACUS.mpi
    # sbatch ../sub.sh
    cd ../
    done
    
    
    
    

为了避免大家计算出错，无法走完整个流程。这里提供两种方式供大家选择：

方式1：批量提交任务，执行`bash batch.sh`, 计算完成后，所有文件存储在`/data/SCF-*`， 大家可以把它复制到当前路径下，采用指令`cp -r /data-SCF-* .`

方式2：这里大家也可以直接用数据集中计算好的的SCF文件来提取三阶力常数。

接下来，以方式2进行演示：


```bash
ls v1
```

    [0m[01;34mSCF-01[0m  [01;34mSCF-06[0m  [01;34mSCF-11[0m  [01;34mSCF-16[0m  [01;34mSCF-21[0m  [01;34mSCF-26[0m  [01;34mSCF-31[0m  [01;34mSCF-36[0m  batch.sh
    [01;34mSCF-02[0m  [01;34mSCF-07[0m  [01;34mSCF-12[0m  [01;34mSCF-17[0m  [01;34mSCF-22[0m  [01;34mSCF-27[0m  [01;34mSCF-32[0m  [01;34mSCF-37[0m  [01;32mscf.tgz[0m
    [01;34mSCF-03[0m  [01;34mSCF-08[0m  [01;34mSCF-13[0m  [01;34mSCF-18[0m  [01;34mSCF-23[0m  [01;34mSCF-28[0m  [01;34mSCF-33[0m  [01;34mSCF-38[0m  [01;34mv1[0m
    [01;34mSCF-04[0m  [01;34mSCF-09[0m  [01;34mSCF-14[0m  [01;34mSCF-19[0m  [01;34mSCF-24[0m  [01;34mSCF-29[0m  [01;34mSCF-34[0m  [01;34mSCF-39[0m
    [01;34mSCF-05[0m  [01;34mSCF-10[0m  [01;34mSCF-15[0m  [01;34mSCF-20[0m  [01;34mSCF-25[0m  [01;34mSCF-30[0m  [01;34mSCF-35[0m  [01;34mSCF-40[0m
    

我们将计算完成的`SCF-*`文件替换掉当前路径下的计算文件


```bash
cp -r v1/SCF-* .
```

运行 aba2vasp.py，将ABACUS计算的原子受力包装成`vasprun.xml` 格式，放置在每个`SCF-*`文件夹中，命令如下：


```bash
python aba2vasp.py
```

    SCF-13
    SCF-30
    SCF-39
    SCF-03
    SCF-26
    SCF-10
    SCF-07
    SCF-37
    SCF-29
    SCF-24
    SCF-15
    SCF-20
    SCF-22
    SCF-34
    SCF-40
    SCF-19
    SCF-05
    SCF-06
    SCF-18
    SCF-11
    SCF-33
    SCF-01
    SCF-28
    SCF-35
    SCF-27
    SCF-16
    SCF-04
    SCF-21
    SCF-12
    SCF-32
    SCF-23
    SCF-02
    SCF-08
    SCF-36
    SCF-25
    SCF-31
    SCF-09
    SCF-38
    SCF-17
    SCF-14
    

`vasprun.xml` 格式示意：


```bash
cat SCF-01/vasprun.xml
```

    <modeling>
    	<calculation>
    		<varray name="forces">
    			<v>0.0024483157 -0.11539483 0.11539483</v>
    			<v>3.091558e-06 -0.0014941451 0.0014941451</v>
    			<v>-0.0048437967 -0.0091443637 -0.003390013</v>
    			<v>-0.0048437967 0.003390013 0.0091443637</v>
    			<v>0.0049111965 0.0033763345 0.0091727893</v>
    			<v>0.0049111965 -0.0091727893 -0.0033763345</v>
    			<v>-6.6667389e-05 0.0082588156 -0.0082588156</v>
    			<v>-9.2824796e-06 0.0013161694 -0.0013161694</v>
    			<v>0.0 -0.0026203459 0.0026145983</v>
    			<v>0.0 -0.0026145983 0.0026203459</v>
    			<v>-0.0028072403 0.00018617914 -0.00018617914</v>
    			<v>0.029352752 0.045301167 -0.045301167</v>
    			<v>0.0027926385 0.00017918469 -0.00017918469</v>
    			<v>-0.031522491 0.046993534 -0.046993534</v>
    			<v>-0.00016280154 0.015791904 -0.015647767</v>
    			<v>-0.00016280154 0.015647767 -0.015791904</v>
    		</varray>
    	</calculation>
    </modeling>
    

最后执行如下命令：


```bash
find SCF-* -name vasprun.xml|sort -n|/opt/mamba/envs/py37/bin/python /thirdorder/thirdorder_vasp.py reap 2 2 2 -2
```

    Reading POSCAR
    Analyzing the symmetries
    - Symmetry group b'Fd-3m' detected
    - 48 symmetry operations
    Creating the supercell
    Computing all distances in the supercell
    - Automatic cutoff: 0.43260363256848766 nm
    Looking for an irreducible set of third-order IFCs
    - 5 triplet equivalence classes found
    - 40 DFT runs are needed
    
     888888ba   88888888b  .d888888   888888ba
     88    `8b  88        d8'    88   88    `8b
    a88aaaa8P' a88aaaa    88aaaaa88a a88aaaa8P'
     88   `8b.  88        88     88   88
     88     88  88        88     88   88
     dP     dP  88888888P 88     88   dP
    oooooooooooooooooooooooooooooooooooooooooooo
    
    XML ElementTree implementation: cElementTree
    Waiting for a list of vasprun.xml files on stdin
    - 40 filenames read
    Reading the forces
    /thirdorder/thirdorder_vasp.py:158: DeprecationWarning: This method will be removed in future versions.  Use 'list(elem)' or iteration over elem instead.
      for i in a.getchildren():
    - SCF-01/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 1.95693375e-08 -2.48125001e-10  2.48125000e-10] eV/(A * atom)
    - SCF-02/vasprun.xml read successfully
    - 	 Average force:
    - 	 [0. 0. 0.] eV/(A * atom)
    - SCF-03/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.48125002e-10  1.95693812e-08  2.48125001e-10] eV/(A * atom)
    - SCF-04/vasprun.xml read successfully
    - 	 Average force:
    - 	 [0. 0. 0.] eV/(A * atom)
    - SCF-05/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.20625001e-10 -2.49999760e-12 -6.44131252e-11] eV/(A * atom)
    - SCF-06/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-3.81937500e-11 -8.74999986e-11  8.74999988e-11] eV/(A * atom)
    - SCF-07/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.20625001e-10 -6.43700000e-11 -2.49999934e-12] eV/(A * atom)
    - SCF-08/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-3.82812501e-11  8.74999992e-11 -8.74999994e-11] eV/(A * atom)
    - SCF-09/vasprun.xml read successfully
    - 	 Average force:
    - 	 [1.1925e-10 0.0000e+00 0.0000e+00] eV/(A * atom)
    - SCF-10/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-1.48675000e-10 -2.16840434e-19 -2.16840434e-19] eV/(A * atom)
    - SCF-11/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.13546375e-08  1.86875000e-10  1.86875001e-10] eV/(A * atom)
    - SCF-12/vasprun.xml read successfully
    - 	 Average force:
    - 	 [2.12500002e-10 0.00000000e+00 4.33680869e-19] eV/(A * atom)
    - SCF-13/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 1.93124999e-10 -2.13546188e-08  1.93125000e-10] eV/(A * atom)
    - SCF-14/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 0.00000000e+00 -4.33680869e-19  2.12500003e-10] eV/(A * atom)
    - SCF-15/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.20625001e-10  2.49999782e-12  6.43743748e-11] eV/(A * atom)
    - SCF-16/vasprun.xml read successfully
    - 	 Average force:
    - 	 [6.00999999e-11 1.95000001e-10 1.95000000e-10] eV/(A * atom)
    - SCF-17/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.14375000e-10  6.43099999e-11  2.49999869e-12] eV/(A * atom)
    - SCF-18/vasprun.xml read successfully
    - 	 Average force:
    - 	 [6.01749998e-11 1.95000001e-10 1.94999999e-10] eV/(A * atom)
    - SCF-19/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-1.35187499e-11 -3.11249997e-11 -2.12500002e-10] eV/(A * atom)
    - SCF-20/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 6.12500001e-10 -1.08420217e-19  1.08420217e-19] eV/(A * atom)
    - SCF-21/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.13546500e-08 -1.93124999e-10 -1.93124999e-10] eV/(A * atom)
    - SCF-22/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.12500004e-10  0.00000000e+00  0.00000000e+00] eV/(A * atom)
    - SCF-23/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-1.93124999e-10 -2.13546437e-08 -1.93125000e-10] eV/(A * atom)
    - SCF-24/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-2.16840434e-19 -2.16840434e-19 -2.12500004e-10] eV/(A * atom)
    - SCF-25/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 3.74999901e-12 -1.22500001e-10  8.15674998e-11] eV/(A * atom)
    - SCF-26/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 6.02000000e-11 -1.95000001e-10 -1.94999999e-10] eV/(A * atom)
    - SCF-27/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 3.74999814e-12  8.14968748e-11 -1.22500001e-10] eV/(A * atom)
    - SCF-28/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 6.02125000e-11 -1.95000000e-10 -1.95000001e-10] eV/(A * atom)
    - SCF-29/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-1.34000003e-11  3.09999999e-11  2.12500002e-10] eV/(A * atom)
    - SCF-30/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-4.99999999e-10 -1.08420217e-19  1.08420217e-19] eV/(A * atom)
    - SCF-31/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 1.95692875e-08  2.48125001e-10 -2.48125001e-10] eV/(A * atom)
    - SCF-32/vasprun.xml read successfully
    - 	 Average force:
    - 	 [0. 0. 0.] eV/(A * atom)
    - SCF-33/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 2.48750000e-10  1.95694438e-08 -2.48750001e-10] eV/(A * atom)
    - SCF-34/vasprun.xml read successfully
    - 	 Average force:
    - 	 [0. 0. 0.] eV/(A * atom)
    - SCF-35/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 3.74999836e-12  1.22500000e-10 -8.15475001e-11] eV/(A * atom)
    - SCF-36/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-3.82812503e-11  8.74999994e-11 -8.75000001e-11] eV/(A * atom)
    - SCF-37/vasprun.xml read successfully
    - 	 Average force:
    - 	 [ 3.74999836e-12 -8.14943753e-11  1.22500000e-10] eV/(A * atom)
    - SCF-38/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-3.82500000e-11 -8.75000005e-11  8.74999986e-11] eV/(A * atom)
    - SCF-39/vasprun.xml read successfully
    - 	 Average force:
    - 	 [1.1925e-10 0.0000e+00 0.0000e+00] eV/(A * atom)
    - SCF-40/vasprun.xml read successfully
    - 	 Average force:
    - 	 [-1.48612499e-10 -2.16840434e-19  0.00000000e+00] eV/(A * atom)
    Computing an irreducible set of anharmonic force constants
    Reconstructing the full array
    - Storing the coefficients in a dense matrix
    Writing the constants to FORCE_CONSTANTS_3RD
    
    888888ba   .88888.  888888ba   88888888b
    88    `8b d8'   `8b 88    `8b  88
    88     88 88     88 88     88 a88aaaa
    88     88 88     88 88     88  88
    88    .8P Y8.   .8P 88     88  88
    8888888P   `8888P'  dP     dP  88888888P
    ooooooooooooooooooooooooooooooooooooooooo
    
    

即可得到三阶力常数文件`FORCE_CONSTANTS_3RD`。在`shengbte`文件夹中提供了`FORCE_CONSTANTS_3rd`文件供参考计算结果。

### 3. 运行 ShengBTE 得到晶格热导率

进入`shengbte`文件夹，里面已经准备好`CONTROL`（ShengBTE 的参数文件）、`FORCE_CONSTANTS_2ND`（二阶力常数文件）、`FORCE_CONSTANTS_3RD`（三阶力常数文件）这三个文件，使用如下命令运行 ShengBTE 即可得到晶格热导率，其中 `Ref` 文件夹中给出了计算结果供参考：


```bash
cd /abacus-user-guide/examples/interface_ShengBTE/LCAO/shengbte
```

`CONTROL`中主要是输入计算体系相关的信息，比如原子数、原子名、盒子大小和坐标信息，计算热导率的温度范围等


```bash
cat CONTROL
```

    &allocations
        nelements=1   
        natoms=2    
        ngrid(:)=10 10 10  
    &end
    &crystal
        lfactor=0.100000 
        lattvec(:,1)=0 2.81594778072 2.81594778072
        lattvec(:,2)=2.81594778072 0 2.81594778072
        lattvec(:,3)=2.81594778072 2.81594778072 0
        elements="Si"   
        types=1 1
        positions(:,1)=0.8750000000000000  0.8750000000000000  0.8750000000000000 
        positions(:,2)=0.1250000000000000  0.1250000000000000  0.1250000000000000
        scell(:)=2 2 2
    &end
    &parameters
            !T=300,
            T_min=200
            T_max=500
            T_step=50
            scalebroad=1.0
    &end
    &flags
            !espresso=.true.
            nonanalytic=.true.,
            isotopes=.true.
    &end
    


```bash
cp ../2nd/FORCE_CONSTANTS_2ND .
cp ../3rd/FORCE_CONSTANTS_3RD .
```


```bash
mpirun -n 4 ShengBTE
```

     Info: symmetry group Fd-3m      detected
     Info:           48  symmetry operations
     Info:          48 duplicated rotations will be discarded
     Info: This calculation is running on 4 MPI process(es)
     Info: Ntot =        1000
     Info: Nlist =          47
     Info: about to obtain the spectrum
     Info: expecting Phonopy 2nd-order format
     Info: about to set the acoustic frequencies at Gamma to zero
     Info: original values:
     Info: omega(1,1) = -4.790893013093741E-003 rad/ps
     Info: omega(1,2) = -4.790892587327632E-003 rad/ps
     Info: omega(1,3) = -4.790892249633965E-003 rad/ps
     Info: spectrum calculation finished in .138 seconds
     Info: start calculating specific heat and kappa in the small-grain limit 
     Info: Temperature =   200.000000000000     
     Info: Temperature =   250.000000000000     
     Info: Temperature =   300.000000000000     
     Info: Temperature =   350.000000000000     
     Info: Temperature =   400.000000000000     
     Info: Temperature =   450.000000000000     
     Info: Temperature =   500.000000000000     
     Info: Ntotal_plus =                394947
     Info: Ntotal_minus =                460218
     Info: max(N_plus), max(N_minus)        3632        2952
     Info: calculating Vp_plus and Vp_minus at .461 seconds
     Info:   0% don  1% don  3% don  4% don  6% don  7% don  8% don 10% don 11% don 13% don 14% don 15% don 17% don 18% don 20% don 21% don 23% don 24% don 25% don 27% don 28% don 30% don 31% don 32% don 34% don 35% don 37% don 38% don 39% don 41% don 42% don 44% don 45% don 46% don 48% don 49% don 51% don 52% don 54% don 55% don 56% don 58% don 59% don 61% don 62% don 63% don 65% don 66% don 68% don 69% don 70% don 72% don 73% don 75% don 76% don 77% don 79% don 80% don 82% don 83% don 85% don 86% don 87% don 89% don 90% don 92% don 93% don 94% don 96% don 97% don 99% don100% don100% done.
     Info: start calculating kappa at 2.883 seconds
     Info: Temperature=   200.000000000000     
     Info: Iteration : 1          Timer : 3.493 seconds
     Info: Relative change =  8.477415246079950E-002
     Info: Iteration : 2          Timer : 3.638 seconds
     Info: Relative change =  3.385986549366961E-002
     Info: Iteration : 3          Timer : 3.883 seconds
     Info: Relative change =  4.623778121298880E-003
     Info: Iteration : 4          Timer : 4.030 seconds
     Info: Relative change =  2.334048481057918E-003
     Info: Iteration : 5          Timer : 4.177 seconds
     Info: Relative change =  2.305418577536500E-004
     Info: Iteration : 6          Timer : 4.322 seconds
     Info: Relative change =  2.059825752046488E-004
     Info: Iteration : 7          Timer : 4.471 seconds
     Info: Relative change =  2.729413333871115E-006
     Info: Temperature=   250.000000000000     
     Info: Iteration : 1          Timer : 4.784 seconds
     Info: Relative change =  8.574845115047762E-002
     Info: Iteration : 2          Timer : 4.929 seconds
     Info: Relative change =  3.444083846868982E-002
     Info: Iteration : 3          Timer : 5.075 seconds
     Info: Relative change =  4.448964600534328E-003
     Info: Iteration : 4          Timer : 5.221 seconds
     Info: Relative change =  2.393367648766057E-003
     Info: Iteration : 5          Timer : 5.366 seconds
     Info: Relative change =  1.731631397906174E-004
     Info: Iteration : 6          Timer : 5.513 seconds
     Info: Relative change =  2.193543930882834E-004
     Info: Iteration : 7          Timer : 5.658 seconds
     Info: Relative change =  9.664034408258441E-006
     Info: Temperature=   300.000000000000     
     Info: Iteration : 1          Timer : 5.968 seconds
     Info: Relative change =  8.707926291454768E-002
     Info: Iteration : 2          Timer : 6.114 seconds
     Info: Relative change =  3.492683468633746E-002
     Info: Iteration : 3          Timer : 6.360 seconds
     Info: Relative change =  4.476358669278945E-003
     Info: Iteration : 4          Timer : 6.507 seconds
     Info: Relative change =  2.450407378745750E-003
     Info: Iteration : 5          Timer : 6.654 seconds
     Info: Relative change =  1.545509839335893E-004
     Info: Iteration : 6          Timer : 6.799 seconds
     Info: Relative change =  2.309266500521971E-004
     Info: Iteration : 7          Timer : 6.946 seconds
     Info: Relative change =  1.604756813141203E-005
     Info: Iteration : 8          Timer : 7.093 seconds
     Info: Relative change =  2.785838281895868E-005
     Info: Iteration : 9          Timer : 7.238 seconds
     Info: Relative change =  6.094532140089185E-006
     Info: Temperature=   350.000000000000     
     Info: Iteration : 1          Timer : 7.702 seconds
     Info: Relative change =  8.825093332029750E-002
     Info: Iteration : 2          Timer : 7.847 seconds
     Info: Relative change =  3.532871054638942E-002
     Info: Iteration : 3          Timer : 7.997 seconds
     Info: Relative change =  4.558904213786180E-003
     Info: Iteration : 4          Timer : 8.142 seconds
     Info: Relative change =  2.499226891382388E-003
     Info: Iteration : 5          Timer : 8.287 seconds
     Info: Relative change =  1.509406908234823E-004
     Info: Iteration : 6          Timer : 8.432 seconds
     Info: Relative change =  2.400448268802154E-004
     Info: Iteration : 7          Timer : 8.578 seconds
     Info: Relative change =  1.948439490063963E-005
     Info: Iteration : 8          Timer : 8.722 seconds
     Info: Relative change =  3.014418412656969E-005
     Info: Iteration : 9          Timer : 8.867 seconds
     Info: Relative change =  7.202792929955250E-006
     Info: Temperature=   400.000000000000     
     Info: Iteration : 1          Timer : 9.176 seconds
     Info: Relative change =  8.921130913722577E-002
     Info: Iteration : 2          Timer : 9.321 seconds
     Info: Relative change =  3.566672475858301E-002
     Info: Iteration : 3          Timer : 9.479 seconds
     Info: Relative change =  4.651981399847585E-003
     Info: Iteration : 4          Timer : 9.624 seconds
     Info: Relative change =  2.541134732108795E-003
     Info: Iteration : 5          Timer : 9.769 seconds
     Info: Relative change =  1.534282536407042E-004
     Info: Iteration : 6          Timer : 9.914 seconds
     Info: Relative change =  2.473911714796271E-004
     Info: Iteration : 7          Timer : 10.060 seconds
     Info: Relative change =  2.139124169039141E-005
     Info: Iteration : 8          Timer : 10.258 seconds
     Info: Relative change =  3.192333259480310E-005
     Info: Iteration : 9          Timer : 10.403 seconds
     Info: Relative change =  7.979542718355900E-006
     Info: Temperature=   450.000000000000     
     Info: Iteration : 1          Timer : 10.712 seconds
     Info: Relative change =  8.999158840301218E-002
     Info: Iteration : 2          Timer : 10.958 seconds
     Info: Relative change =  3.595529785049464E-002
     Info: Iteration : 3          Timer : 11.105 seconds
     Info: Relative change =  4.741486740580683E-003
     Info: Iteration : 4          Timer : 11.252 seconds
     Info: Relative change =  2.577595330046745E-003
     Info: Iteration : 5          Timer : 11.458 seconds
     Info: Relative change =  1.583822052349165E-004
     Info: Iteration : 6          Timer : 11.603 seconds
     Info: Relative change =  2.535072665682756E-004
     Info: Iteration : 7          Timer : 11.748 seconds
     Info: Relative change =  2.245899472208437E-005
     Info: Iteration : 8          Timer : 11.894 seconds
     Info: Relative change =  3.335478069935806E-005
     Info: Iteration : 9          Timer : 12.158 seconds
     Info: Relative change =  8.544499734196537E-006
     Info: Temperature=   500.000000000000     
     Info: Iteration : 1          Timer : 12.566 seconds
     Info: Relative change =  9.063047144963478E-002
     Info: Iteration : 2          Timer : 12.711 seconds
     Info: Relative change =  3.620455528701171E-002
     Info: Iteration : 3          Timer : 12.856 seconds
     Info: Relative change =  4.823294043678694E-003
     Info: Iteration : 4          Timer : 13.004 seconds
     Info: Relative change =  2.609673797576762E-003
     Info: Iteration : 5          Timer : 13.149 seconds
     Info: Relative change =  1.642228226383847E-004
     Info: Iteration : 6          Timer : 13.296 seconds
     Info: Relative change =  2.587358855828025E-004
     Info: Iteration : 7          Timer : 13.443 seconds
     Info: Relative change =  2.304443208249280E-005
     Info: Iteration : 8          Timer : 13.589 seconds
     Info: Relative change =  3.454185013654348E-005
     Info: Iteration : 9          Timer : 13.734 seconds
     Info: Relative change =  8.969125652120961E-006
     Info: normal exit after 13.749 seconds
    


```bash
cat BTE.KappaTensorVsT_CONV
```

      200.0   0.15308E+03   0.15529E-19  -0.55339E-19  -0.92609E-19   0.15308E+03   0.67763E-20  -0.46869E-19  -0.67763E-19   0.15308E+03     7
      250.0   0.12031E+03  -0.18070E-19  -0.49693E-19  -0.42916E-19   0.12031E+03   0.47646E-20  -0.27839E-18  -0.40658E-19   0.12031E+03     7
      300.0   0.99540E+02  -0.67763E-20  -0.22588E-19  -0.22588E-20   0.99540E+02  -0.25411E-20  -0.81880E-20  -0.10729E-19   0.99540E+02     9
      350.0   0.85108E+02   0.16037E-18  -0.17844E-18  -0.16263E-18   0.85108E+02   0.16517E-18  -0.12056E-18  -0.99950E-19   0.85108E+02     9
      400.0   0.74448E+02   0.21458E-18  -0.17788E-18  -0.25072E-18   0.74448E+02   0.16037E-18   0.16235E-18  -0.26653E-18   0.74448E+02     9
      450.0   0.66226E+02   0.25072E-18   0.17505E-19   0.18070E-19   0.66226E+02   0.19199E-18   0.26696E-18   0.22446E-19   0.66226E+02     9
      500.0   0.59677E+02  -0.58728E-19   0.70021E-19  -0.13411E-19   0.59677E+02  -0.58728E-19  -0.13256E-18   0.84703E-21   0.59677E+02     9
    

可以看到：300 K下热导率为99.54 W/mK

## 四、结尾

对于 ABACUS 中使用平面波（PW）来做 ShengBTE 的计算也是采用以上类似的流程，但要注意使用平面波时，计算三阶力常数的 `INPUT` 中scf_thr 需要至少小到1e-12。通过计算结果可以发现，PW 和 LCAO 基组计算出的 Si 的晶格热导率是接近的，300 K 下均在 100 W/(mK) 左右，而实验中 Si 在 300 K 的热导率在 150 W/(m K) 附近。这是因为作为教学例子，这里使用的是 2\*2\*2 的扩胞以及 2*2\*2 的 K 点，导致计算结果偏小，实际科研中需要测试扩胞的大小以及 K 点的采样方案来达到收敛的结果。以上就是 ABACUS+ShengBTE 计算晶格热导率的全部流程，如果有什么问题，欢迎通过邮件联系。
