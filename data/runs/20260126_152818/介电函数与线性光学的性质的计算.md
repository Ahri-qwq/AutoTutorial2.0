<div style="color:black; background-color:#E4EAF2; border: 5px solid #335A87; border-radius: 12px; margin-bottom:1rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2023 @ Authors</i></b><br/>
        <i>共享协议：</a>本作品采用<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议</a>进行许可。</i><br/>  
        <i>镜像：<span style='color:rgb(85,91,228); font-weight:bold'> registry.dp.tech/dptech/prod-19853/abacus-pyatb-open:v0.0.1 </span></i><br/>
        <i>推荐配置：<span style='color:rgb(85,91,228); font-weight:bold'> `c16_m32_cpu` </span></i><br/>
        <i>点击界面上方蓝色按钮 <span style="background-color:rgb(85, 91, 228); color:white; padding: 3px; border-radius: 5px;box-shadow: 2px 2px 3px rgba(0, 0, 0, 0.3); font-size:0.75rem;"> 开始连接 </span> ，选择上述镜像和推荐配置后，稍等片刻即可运行。</i><br/>
    </p>
</div>

# 背景
在当今迅速发展的材料科学领域中，理解和预测材料的电子性质变得越来越重要。介电函数是表征材料对电磁波响应能力的重要参数，直接关联到材料的光学、电子及光谱特性。随着计算方法和计算机技术的进步，第一性原理计算已经成为研究物性的有力工具。
常见的线性光学性质如下：
- 光导率（Optical Conductivity）：描述材料对电磁波（特别是可见光）的传导能力。

- 介电函数：介电函数是表征介质对电磁场作用响应的复数函数，它由实部和虚部组成：它分为实部（$\epsilon_1$）和虚部（$\epsilon_2$），其中实部与材料的折射特性有关，虚部则与材料的吸收特性相关，也因此与材料的带隙有关。

- 折射率（$n$）：描述光从一个介质（例如真空或空气）进入另一个介质（例如玻璃或水）时速度减慢的程度。折射率越高，光在材料中传播的速度差异越大。

- 吸收系数（$\alpha$）：描述光在材料中传播时衰减的程度，与材料吸收光的能力有关。吸收系数越高，光在材料内被吸收得越多，传播得越短。

- 消光系数（$\kappa$）：用于表示在介电质材料中光波传播时由于散射和吸收引起的能量衰减的指标。

本notebook将以晶态二氧化硅为例子，演示使用ABACUS结合pyatb计算材料的介电函数的方法，并且通过后处理脚本展示了通过介电函数计算各线性光学性质的流程。


# 计算原理
## 原子轨道基组与紧束缚模型的构建

在周期性系统中，给定 $k$点的 Kohn-Sham 方程可以写为：
$$
H|\Psi_{n\mathbf{k}}\rangle = E_{nk} |\Psi_{n\mathbf{k}}\rangle.
$$
这里 $\Psi_{n\mathbf{k}}$ 是第$n$个带的布洛赫波函数，可以用数值原子轨道（numeric atomic orbital, NAO）表达为：
$$
|\Psi_{n\mathbf{k}}\rangle = \frac{1}{\sqrt{N}} \sum_{\mu} C_{\mu n}(\mathbf{k}) \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}|\mathbf{R}\mu\rangle,
$$
其中 $|\mu\mathbf{R}\rangle$ 是第 $\mu$ 个原子轨道，在$R^\mathbf{th}$ 单元胞中，而$\mathbf{r}_{\mu}$表示这个原子轨道的中心位置。复合指标$\mu = (a, i, l, c, m)$，其中 a 是元素类型，$i$是每种元素类型的原子的索引，$c$ 是角动量$l$的径向函数的多重性，$m$是磁量子数。NAO 的系数由$C_{\mu n}(\mathbf{k})$ 给出。

在 NAO 基下，Kohn-Sham 方程成为一个本征值问题，
$$
H(\mathbf{k})C_n(\mathbf{k}) = E_{n\mathbf{k}}S(\mathbf{k})C_n .
$$
其中$H(\mathbf{k})$、$S(\mathbf{k})$和 $C_n(\mathbf{k})$ 分别是哈密顿矩阵、重叠矩阵和第 n 个带的本征矢量。

为了获得$H(\mathbf{k})$和$S(\mathbf{k})$，我们首先使用基于 NAO 的第一性原理软件计算实空间中的紧束缚哈密顿量，
$$
H_{\nu\mu}(\mathbf{R}) = \langle \nu | H | \mathbf{R}\mu  \rangle,
$$

$$
S_{\nu\mu}(\mathbf{R}) = \langle \nu | \mathbf{R} \mu  \rangle.
$$

一旦我们有了$H_{\mu\nu}(\mathbf{R})$ 和$S_{\mu\nu}(\mathbf{R})$，我们就可以使用下列关系获得任意 k 点的哈密顿矩阵和重叠矩阵，
$$
H_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}H_{\mu\nu}(\mathbf{R}),
$$

$$
S_{\nu\mu}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}S_{\mu\nu}(\mathbf{R}),
$$

我们还需要 NAO 之间的偶极矩阵，即：
$$
r_{\nu\mu,a}(\mathbf{R}) = \langle \mu | r_a | \nu \mathbf{R} \rangle \quad a = x, y, z.
$$
然后可以通过傅立叶变换获得偶极矩阵 $A^{R}_{\nu\mu,a}(\mathbf{k})$，
$$
A^{R}_{\nu\mu,a}(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}r_{\mu\nu,a}(\mathbf{R}).
$$
在获得紧束缚参数$H_{\nu\mu}(\mathbf{R})$、$S_{\nu\mu}(\mathbf{R})$和 $r_{\nu\mu,a}(\mathbf{R})$之后，可以使用 PYATB 计算电子结构和相关物理属性：

1. **Hamiltonian 和 Overlap 矩阵**：$H_{\mu\nu}(\mathbf{R})$和$S_{\mu\nu}(\mathbf{R})$分别是用于构建$\mathbf{k}$点上的 Hamiltonian 矩阵$H_{\mu\nu}(\mathbf{k})$和重叠矩阵 $S_{\mu\nu}(\mathbf{k})$ 。通过对所有的$\mathbf{R}$点求和得到对于每个$\mathbf{k}$点的矩阵。
2. **计算本征值**：通过解决$H_{\mu\nu}(\mathbf{k})$和 $S_{\mu\nu}(\mathbf{k})$ 形成的本征值问题，可以得到电子的能带结构，即各个能带的能量$E_{nk}$以及波函数的系数。
3. **偶极矩阵**：使用之前得到的波函数系数和 $A^{R}_{\mu\nu,a}(\mathbf{k})$，可以计算出速度矩阵元：

$$
\langle n\mathbf{k}|v{\alpha}|m\mathbf{k}\rangle = \frac{1}{\hbar} \sum_{\mu,\nu} C_{n\mu}(k) C_{m\nu}(k) A^\alpha_{\mu,\nu}(k) (E_{nk} - E_{mk})
$$

这里的 $C_{n\mu}(k)$ 和$C_{m\nu}(k)$是哈密顿量和重叠矩阵的本征向量在原子轨道基态下的系数。



## 利用紧束缚模型计算线性光学性质

现在我们有了所有必要的矩阵参数，可以计算Kubo-Greenwood方程，含频的的光电导率（optical conductivity）由Kubo-Greenwood公式表达为：
$$
\sigma_{\alpha\beta}(\hbar\omega) = -\frac{ie^2\hbar}{NV_\text{cell}} \sum_{\mathbf{k}} \sum_{n,m} (\frac{f_{n\mathbf{k}} - f_{m\mathbf{k}}}{E_{n\mathbf{k}} - E_{m\mathbf{k}}}) \frac{\langle n\mathbf{k}|v_{\alpha}|m\mathbf{k}\rangle \langle m\mathbf{k}|v_{\beta}|n\mathbf{k}\rangle}{\hbar\omega + E_{n\mathbf{k}} - E_{m\mathbf{k}} + i\eta}
$$
其中$f_{nk}$和$f_{mk}$ 是Fermi-Dirac占据函数，其值取决于电子的能量和温度。速度矩阵元是通过波函数的系数从偶极矩阵元计算出来的。之后，在给定的频率 ( $\omega$ )下，对所有的 k 点和能带指数 ( $n$ )和 ( $m$ )求和，从而得到光电导率 ($\sigma_{\alpha\beta}(\hbar\omega)$)
光导率与介电函数的关系是
$$
\epsilon(\omega)=\epsilon_0+ i\frac{\sigma(\omega)}{\omega}
$$
而在pyatb程序中，介电函数的虚部由下式单独计算以避免在$\omega=0$处的额外处理： 
$$
\epsilon^{\alpha\beta}_2(\omega) = \frac{e^2\pi}{\epsilon_0\hbar} \int \frac{d\mathbf{k}}{(2\pi)^3} \sum_{n,m} f_{nm} r^{\alpha}_{nm}r^{\beta}_{mn}\delta(\omega_{nm} - \omega)
$$

介电函数的实部通过Kramer-Kronig变换得到
$$
\epsilon^{\alpha\beta}_1(\omega) = \delta_{\alpha\beta} + \frac{2}{\pi} \mathbf{P} \int_0^\infty \frac{\omega' \epsilon^{\alpha\beta}_i(\omega')}{\omega'^2 - \omega^2} d\omega'
$$
线性光谱可以通过介电函数计算
例如折射率$n(\omega)$
消光系数$\kappa(\omega)$
吸光系数$\alpha(\omega)$
能量损失函数$L(\omega)$
以及反射率$R(\omega)$
$$
n(\omega) =  \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} + \epsilon_1}{2}}
$$
$$
\kappa(\omega) =  \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1}{2}}
$$
$$
\alpha(\omega) = \sqrt{\frac{2\omega^2}{c^2} \left( \sqrt{\epsilon_1^2 + \epsilon_2^2} + \epsilon_1 \right)}
$$

$$
L(\omega) = \text{Im} \left( \frac{-1}{\epsilon(\omega)} \right) = \frac{\epsilon_2}{\epsilon_1^2 + \epsilon_2^2}
$$

$$
R(\omega) =  \frac{(n - 1)^2 + k^2}{(n + 1)^2 + k^2}
$$

## References
- pyatb Manual: https://pyatb.github.io/pyatb/functions/optical_conductivity.html
- https://en.wikipedia.org/wiki/Optical_conductivity
- https://www.openmx-square.org/tech_notes/Dielectric_Function_YTL.pdf


# 对晶态$\mathbf{SiO_2}$单胞的计算
## ABACUS自洽计算



```
cd ~/abacus-pyatb_tutorial/silica_PrimaryCell
```

    /root/abacus-pyatb_tutorial/silica_PrimaryCell
    /usr/local/lib/python3.10/dist-packages/IPython/core/magics/osm.py:393: UserWarning: using bookmarks requires you to install the `pickleshare` library.
      bkms = self.shell.db.get('bookmarks', {})
    /usr/local/lib/python3.10/dist-packages/IPython/core/magics/osm.py:417: UserWarning: using dhist requires you to install the `pickleshare` library.
      self.shell.db['dhist'] = compress_dhist(dhist)[-100:]
    

* 查看ABACUS的输入文件。
1. STRU文件可以基于POSCAR文件利用dpdata库转换得到
2. INPUT文件中应当注意out_chg、out_mat_hs2、out_mat_r均应当设置为1
3. KPT中记录了k点的选取方式和密度


```
ls && cat INPUT && cat STRU && cat KPT
```

    INPUT               O_gga_10au_100Ry_3s3p2d.orb  Si_gga_10au_100Ry_3s3p2d.orb
    KPT                 O_gga_7au_100Ry_2s2p1d.orb   Si_gga_7au_100Ry_2s2p1d.orb
    [0m[01;34mOUT.silica[0m/         STRU                         [01;34mpyatb_OpticConductivity[0m/
    O_ONCV_PBE-1.0.upf  Si_ONCV_PBE-1.0.upf          time.json
    INPUT_PARAMETERS
    
    suffix                  silica
    calculation             scf	
    esolver_type     ksdft
    symmetry          0
    init_chg             atomic
    
    
    pseudo_dir              ./
    orbital_dir		./
    #kspacing 0.25
    #gamma_only 1
    
    basis_type              lcao  
    ks_solver               genelpa
    #nspin
    smearing_method    gaussian
    smearing_sigma       0.01
    mixing_type             broyden
    mixing_beta             0.1
    #scf_max                200
    ecutwfc                 100             # Rydberg
    scf_thr                 1e-8		# Rydberg
    mixing_gg0              1.5  
    mixing_ndim              20
    
    
    out_chg              1
    out_mat_hs2      1
    out_mat_r        1
    ATOMIC_SPECIES
    Si 28.086  Si_ONCV_PBE-1.0.upf
    O 15.999   O_ONCV_PBE-1.0.upf
    
    NUMERICAL_ORBITAL
    Si_gga_7au_100Ry_2s2p1d.orb
    O_gga_7au_100Ry_2s2p1d.orb
    
    LATTICE_CONSTANT
    1.8897261246257702
    
    LATTICE_VECTORS
    7.1199998856 0.0 0.0 
    0.0 7.1199998856 0.0 
    0.0 0.0 7.1199998856 
    
    ATOMIC_POSITIONS
    Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)
    Si
    0.0
    8
    0.000000000000 0.000000000000 0.000000000000 1 1 1
    0.000000000000 3.559999943000 3.559999943000 1 1 1
    3.559999943000 3.559999943000 0.000000000000 1 1 1
    3.559999943000 0.000000000000 3.559999943000 1 1 1
    5.339999914000 1.779999971000 5.339999914000 1 1 1
    1.779999971000 1.779999971000 1.779999971000 1 1 1
    1.779999971000 5.339999914000 5.339999914000 1 1 1
    5.339999914000 5.339999914000 1.779999971000 1 1 1
    O
    0.0
    16
    0.889999986000 0.889999986000 0.889999986000 1 1 1
    6.229999900000 2.669999957000 4.449999928000 1 1 1
    2.669999957000 4.449999928000 6.229999900000 1 1 1
    4.449999928000 6.229999900000 2.669999957000 1 1 1
    0.889999986000 4.449999928000 4.449999928000 1 1 1
    6.229999900000 6.229999900000 0.889999986000 1 1 1
    2.669999957000 0.889999986000 2.669999957000 1 1 1
    4.449999928000 2.669999957000 6.229999900000 1 1 1
    4.449999928000 0.889999986000 4.449999928000 1 1 1
    2.669999957000 2.669999957000 0.889999986000 1 1 1
    6.229999900000 4.449999928000 2.669999957000 1 1 1
    0.889999986000 6.229999900000 6.229999900000 1 1 1
    4.449999928000 4.449999928000 0.889999986000 1 1 1
    2.669999957000 6.229999900000 4.449999928000 1 1 1
    6.229999900000 0.889999986000 6.229999900000 1 1 1
    0.889999986000 2.669999957000 2.669999957000 1 1 1
    K_POINTS
    0 
    Gamma
    6 6 6 0 0 0
    


```
! export OMP_NUM_THREADS=1 && mpirun -np 16 abacus
```

                                                                                         
                                  ABACUS v3.4.3
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 67412ea (Tue Nov 28 09:18:45 2023 +0800)
    
     Thu Jan 11 19:37:12 2024
     MAKE THE DIR         : OUT.silica/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
     UNIFORM GRID DIM        : 90 * 90 * 90
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0530675  SEC) : SETUP UNITCELL
     DONE(0.0576298  SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  NBASE       
     1       112             16          312         
     ---------------------------------------------------------
     Use Systematically Improvable Atomic bases
     ---------------------------------------------------------
     ELEMENT ORBITALS        NBASE       NATOM       XC          
     Si      2s2p1d-7au      13          8           
     O       2s2p1d-7au      13          16          
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.134316   SEC) : INIT PLANEWAVE
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.20528    SEC) : INIT SCF
     ITER   ETOT(eV)       EDIFF(eV)      DRHO       TIME(s)    
     GE1    -7.851283e+03  0.000000e+00   2.125e-01  1.011e+01  
     GE2    -7.826465e+03  2.481833e+01   1.691e-01  9.909e+00  
     GE3    -7.835437e+03  -8.971779e+00  4.789e-02  1.018e+01  
     GE4    -7.834963e+03  4.732887e-01   3.151e-02  9.615e+00  
     GE5    -7.835196e+03  -2.321661e-01  2.469e-03  9.415e+00  
     GE6    -7.835165e+03  3.109536e-02   6.826e-03  9.413e+00  
     GE7    -7.835176e+03  -1.125703e-02  5.259e-04  1.001e+01  
     GE8    -7.835176e+03  -7.333078e-05  6.149e-05  9.734e+00  
     GE9    -7.835176e+03  -7.536110e-07  3.825e-06  1.109e+01  
     GE10   -7.835176e+03  -1.180200e-09  5.955e-07  9.737e+00  
     GE11   -7.835176e+03  -3.557614e-11  1.918e-07  9.821e+00  
     GE12   -7.835176e+03  3.093578e-11   1.843e-07  1.017e+01  
     GE13   -7.835176e+03  -1.392110e-11  7.527e-10  9.462e+00  
    TIME STATISTICS
    ------------------------------------------------------------------------------------
         CLASS_NAME                 NAME            TIME(Sec)  CALLS   AVG(Sec) PER(%)
    ------------------------------------------------------------------------------------
                         total                      155.45           9  17.27   100.00
    Driver               reading                      0.02           1   0.02     0.01
    Input                Init                         0.00           1   0.00     0.00
    Input_Conv           Convert                      0.00           1   0.00     0.00
    Driver               driver_line                155.43           1 155.43    99.99
    UnitCell             check_tau                    0.00           1   0.00     0.00
    PW_Basis_Sup         setuptransform               0.01           1   0.01     0.01
    PW_Basis_Sup         distributeg                  0.01           1   0.01     0.01
    mymath               heapsort                     0.00           3   0.00     0.00
    PW_Basis_K           setuptransform               0.03           1   0.03     0.02
    PW_Basis_K           distributeg                  0.01           1   0.01     0.00
    PW_Basis             setup_struc_factor           0.03           1   0.03     0.02
    ORB_control          read_orb_first               0.12           1   0.12     0.08
    LCAO_Orbitals        Read_Orbitals                0.12           1   0.12     0.08
    NOrbital_Lm          extra_uniform                0.18         181   0.00     0.12
    Mathzone_Add1        SplineD2                     0.00         181   0.00     0.00
    Mathzone_Add1        Cubic_Spline_Interpolation   0.08         181   0.00     0.05
    Sphbes               Spherical_Bessel             0.06        8040   0.00     0.04
    ppcell_vl            init_vloc                    0.04           1   0.04     0.02
    Ions                 opt_ions                   154.91           1 154.91    99.65
    ESolver_KS_LCAO      Run                        154.91           1 154.91    99.65
    ESolver_KS_LCAO      beforescf                    0.70           1   0.70     0.45
    ESolver_KS_LCAO      beforesolver                 0.57           1   0.57     0.37
    ESolver_KS_LCAO      set_matrix_grid              0.37           1   0.37     0.24
    atom_arrange         search                       0.00           1   0.00     0.00
    Grid_Technique       init                         0.29           1   0.29     0.18
    Grid_BigCell         grid_expansion_index         0.02           2   0.01     0.01
    Record_adj           for_2d                       0.08           1   0.08     0.05
    Grid_Driver          Find_atom                    0.01         552   0.00     0.01
    LCAO_Hamilt          grid_prepare                 0.00           1   0.00     0.00
    Veff                 initialize_HR                0.00           1   0.00     0.00
    OverlapNew           initialize_SR                0.00           1   0.00     0.00
    EkineticNew          initialize_HR                0.00           1   0.00     0.00
    NonlocalNew          initialize_HR                0.01           1   0.01     0.01
    Charge               set_rho_core                 0.00           1   0.00     0.00
    Charge               atomic_rho                   0.05           1   0.05     0.03
    PW_Basis_Sup         recip2real                   0.17          72   0.00     0.11
    PW_Basis_Sup         gathers_scatterp             0.08          72   0.00     0.05
    Potential            init_pot                     0.07           1   0.07     0.05
    Potential            update_from_charge           0.98          14   0.07     0.63
    Potential            cal_fixed_v                  0.00           1   0.00     0.00
    PotLocal             cal_fixed_v                  0.00           1   0.00     0.00
    Potential            cal_v_eff                    0.98          14   0.07     0.63
    H_Hartree_pw         v_hartree                    0.09          14   0.01     0.06
    PW_Basis_Sup         real2recip                   0.38          84   0.00     0.25
    PW_Basis_Sup         gatherp_scatters             0.29          84   0.00     0.18
    PotXC                cal_v_eff                    0.88          14   0.06     0.57
    XC_Functional        v_xc                         0.88          14   0.06     0.57
    Potential            interpolate_vrs              0.00          14   0.00     0.00
    H_Ewald_pw           compute_ewald                0.00           1   0.00     0.00
    HSolverLCAO          solve                      127.54          13   9.81    82.05
    HamiltLCAO           updateHk                    15.29        1456   0.01     9.84
    OperatorLCAO         init                        13.68        4368   0.00     8.80
    Veff                 contributeHR                 8.67          13   0.67     5.58
    Gint_interface       cal_gint                    16.22          26   0.62    10.43
    Gint_interface       cal_gint_vlocal              8.02          13   0.62     5.16
    Gint_Tools           cal_psir_ylm                 4.22       46800   0.00     2.72
    Gint_k               transfer_pvpR                0.65          13   0.05     0.42
    OverlapNew           calculate_SR                 0.11           1   0.11     0.07
    OverlapNew           contributeHk                 1.50        1456   0.00     0.96
    EkineticNew          contributeHR                 0.11          13   0.01     0.07
    EkineticNew          calculate_HR                 0.11           1   0.11     0.07
    NonlocalNew          contributeHR                 0.12          13   0.01     0.08
    NonlocalNew          calculate_HR                 0.11           1   0.11     0.07
    OperatorLCAO         contributeHk                 4.53        1456   0.00     2.91
    HSolverLCAO          hamiltSolvePsiK             96.18        1456   0.07    61.87
    DiagoElpa            elpa_solve                  94.49        1456   0.06    60.78
    ElecStateLCAO        psiToRho                    16.07          13   1.24    10.33
    elecstate            cal_dm                       2.12          13   0.16     1.36
    psiMulPsiMpi         pdgemm                       2.10        1456   0.00     1.35
    DensityMatrix        cal_DMR                      3.96          13   0.30     2.55
     Local_Orbital_wfc   wfc_2d_to_grid               1.18        1568   0.00     0.76
    Gint                 transfer_DMR                 0.60          13   0.05     0.39
    Gint_interface       cal_gint_rho                 8.20          13   0.63     5.28
    Charge_Mixing        get_drho                     0.00          13   0.00     0.00
    Charge               mix_rho                      0.07          12   0.01     0.04
    Charge               Broyden_mixing               0.04          12   0.00     0.02
    ModuleIO             output_HS_R                  1.40           1   1.40     0.90
    ModuleIO             save_HSR_sparse              1.37           1   1.37     0.88
    cal_r_overlap_R      init                        21.21           1  21.21    13.64
    ORB_gaunt_table      init_Gaunt_CH                0.00           1   0.00     0.00
    ORB_gaunt_table      Calc_Gaunt_CH                0.00         601   0.00     0.00
    ORB_gaunt_table      init_Gaunt                   0.01           1   0.01     0.01
    ORB_gaunt_table      Get_Gaunt_SH                 0.01        6272   0.00     0.00
    ORB_table_phi        cal_ST_Phi12_R              16.73         384   0.04    10.76
    cal_r_overlap_R      out_rR_other                 2.39           1   2.39     1.54
    ModuleIO             write_istate_info            0.03           1   0.03     0.02
    ------------------------------------------------------------------------------------
    
     ----------------------------------------------------------
    
     START  Time  : Thu Jan 11 19:37:12 2024
     FINISH Time  : Thu Jan 11 19:39:47 2024
     TOTAL  Time  : 155
     SEE INFORMATION IN : OUT.silica/
    


```
#可以按需查看输出文件的内容
#! cat OUT*/running_scf.log
! ls *
! grep 'occupied bands' ./OUT*/running_scf.log && grep 'E_Fermi' ./OUT*/running_scf.log
```

    INPUT			     STRU
    KPT			     Si_ONCV_PBE-1.0.upf
    O_ONCV_PBE-1.0.upf	     Si_gga_10au_100Ry_3s3p2d.orb
    O_gga_10au_100Ry_3s3p2d.orb  Si_gga_7au_100Ry_2s2p1d.orb
    O_gga_7au_100Ry_2s2p1d.orb   time.json
    
    OUT.silica:
    INPUT			data-HR-sparse_SPIN0.csr  kpoints
    SPIN1_CHG.cube		data-SR-sparse_SPIN0.csr  running_scf.log
    STRU_READIN_ADJUST.cif	data-rR-sparse.csr	  warning.log
    STRU_SIMPLE.cif		istate.info
    
    pyatb_OpticConductivity:
    Input
                               occupied bands = 64
    E_Fermi            0.2511464169         3.4170223014
    E_Fermi            0.2808148069         3.8206814565
    E_Fermi            0.3718619760         5.0594417433
    E_Fermi            0.3822231610         5.2004128978
    E_Fermi            0.4081010939         5.5525002372
    E_Fermi            0.4067718295         5.5344146672
    E_Fermi            0.4069691655         5.5370995607
    E_Fermi            0.4070841441         5.5386639258
    E_Fermi            0.4070695480         5.5384653351
    E_Fermi            0.4070746772         5.5385351216
    E_Fermi            0.4070749118         5.5385383140
    E_Fermi            0.4070748996         5.5385381470
    E_Fermi                0.4070749075         5.5385382545
    

## 将产生HR, SR, rR文件复制到pyatb的目录下


```
! cp OUT*/data* ./pyatb_OpticConductivity
```

## 线性光学的计算
### 介电函数的计算


```
cd pyatb_OpticConductivity 

```

    /root/abacus-pyatb_tutorial/silica_PrimaryCell/pyatb_OpticConductivity
    

**检查Fermi energy、占据能级数与SCF输出文件的一致性，将输入文件的晶格信息按以下格式填入Input文件中**


```
!cat Input 
```

    INPUT_PARAMETERS
    {
    nspin 1
    package ABACUS
    fermi_energy 5.5385382545
    fermi_energy_unit eV
    HR_route data-HR-sparse_SPIN0.csr
    SR_route data-SR-sparse_SPIN0.csr
    rR_route data-rR-sparse.csr
    HR_unit Ry
    rR_unit Bohr
    }
    LATTICE
    {
    lattice_constant 1.8897261246257702
    lattice_constant_unit Bohr
    lattice_vector
    7.1199998856 0.0 0.0 
    0.0 7.1199998856 0.0 
    0.0 0.0 7.1199998856 
    
    }
    OPTICAL_CONDUCTIVITY
    {
     occ_band 64
     omega 0 30
     domega 0.01
     eta 0.1
     grid 20 20 20
    }
    
    

接下来，运行pyatb程序进行计算，运行过程大约持续数分钟且不会再屏幕上输出信息，请耐心等待至运行完毕即可。输出文件为./Out目录下的running.log文件，用命令行提交时可以追踪此文件的输出。


```
! export OMP_NUM_THREADS=1 && mpirun -np 16 pyatb 
```


```
cd Out/Optical_Conductivity
```

    /root/abacus-pyatb_tutorial/silica_PrimaryCell/pyatb_OpticConductivity/Out/Optical_Conductivity
    


```
!head -n 5 diele*
```

    ==> dielectric_function_imag_part.dat <==
    # omega(eV)            xx             xy             xz             yx             yy             yz             zx             zy             zz
       0.00000   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
       0.01000   9.564624e-06  -6.324391e-15  -8.598452e-15  -6.324414e-15   9.564624e-06  -4.617589e-15  -8.598467e-15  -4.617598e-15   9.564624e-06
       0.02000   1.912933e-05  -1.264891e-14  -1.719705e-14  -1.264892e-14   1.912933e-05  -9.235269e-15  -1.719707e-14  -9.235286e-15   1.912933e-05
       0.03000   2.869419e-05  -1.897362e-14  -2.579590e-14  -1.897364e-14   2.869419e-05  -1.385316e-14  -2.579593e-14  -1.385317e-14   2.869419e-05
    
    ==> dielectric_function_real_part.dat <==
    # omega(eV)            xx             xy             xz             yx             yy             yz             zx             zy             zz
       0.00000   1.896644e+00   1.012614e-10  -1.632027e-10   1.012614e-10   1.896644e+00   1.692814e-10  -1.632027e-10   1.692814e-10   1.896644e+00
       0.01000   1.896644e+00   1.012611e-10  -1.632031e-10   1.012611e-10   1.896644e+00   1.692812e-10  -1.632031e-10   1.692812e-10   1.896644e+00
       0.02000   1.896646e+00   1.012601e-10  -1.632044e-10   1.012601e-10   1.896646e+00   1.692805e-10  -1.632044e-10   1.692805e-10   1.896646e+00
       0.03000   1.896648e+00   1.012585e-10  -1.632066e-10   1.012585e-10   1.896648e+00   1.692793e-10  -1.632066e-10   1.692793e-10   1.896648e+00
    


```
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RealPartData = '/root/abacus-pyatb_tutorial/silica_PrimaryCell/pyatb_OpticConductivity/Out/Optical_Conductivity/dielectric_function_real_part.dat'
ImagPartData = '/root/abacus-pyatb_tutorial/silica_PrimaryCell/pyatb_OpticConductivity/Out/Optical_Conductivity/dielectric_function_imag_part.dat'

# 读取虚部数据和标签
with open(ImagPartData, 'r') as file:
    labels_imag = file.readline().strip().split()[1:]
data_imag = np.loadtxt(ImagPartData, skiprows=1)
omega = data_imag[:, 0]
imag_parts = data_imag[:, 1:]

# 读取实部数据和标签
with open(RealPartData, 'r') as file:
    labels_real = file.readline().strip().split()[1:]
data_real = np.loadtxt(RealPartData, skiprows=1)
real_parts = data_real[:, 1:]

# 定义不同的标记形状
markers = ['o', '*', '^', 'd', 's', 'p', '+', 'x', 'h']

# 设置图形布局
fig = plt.figure(constrained_layout=True)
gs = GridSpec(1, 2, figure=fig)

# 绘制虚部图表
ax_imag = fig.add_subplot(gs[0, 0])
for i, (imag_part, marker) in enumerate(zip(imag_parts.T, markers), start=1):
    ax_imag.plot(omega, imag_part, linestyle='-', marker=marker, markersize=5, fillstyle="full", label=labels_imag[i])
ax_imag.legend()
ax_imag.set_title('Dielectric Function Imaginary Part')
ax_imag.set_xlabel('Energy (eV)')
ax_imag.set_ylabel('Imaginary Part')

# 绘制实部图表
ax_real = fig.add_subplot(gs[0, 1])
for i, (real_part, marker) in enumerate(zip(real_parts.T, markers), start=1):
    ax_real.plot(omega, real_part, linestyle='-', marker=marker, markersize=5, fillstyle="none",mew=0.4, label=labels_real[i])
ax_real.legend()
ax_real.set_title('Dielectric Function Real Part')
ax_real.set_xlabel('Energy (eV)')
ax_real.set_ylabel('Real Part')

# 显示图表
plt.show()
```


    
![png](output_18_0.png)
    



```
import matplotlib.pyplot as plt
import numpy as np


# 定义一个函数用于读取文件并计算结果
def read_and_process(filename):
    # 读取文件
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 去掉标签行，只处理数据行
    data_lines = lines[1:]

    # 初始化两个数组用于存储计算结果和能量
    result = []
    energy = []

    # 遍历每行数据
    for line in data_lines:
        # 以空格为分隔符分割数据
        split_line = line.split()

        # 判断是否是一个完整的数据行
        if len(split_line) < 10:
            continue
        
        # 将第一列的能量值存储到能量数组中
        energy.append(float(split_line[0]))

        # 解析xx, yy, zz三列的值
        xx = float(split_line[1])
        yy = float(split_line[5])
        zz = float(split_line[9])
        xy = float(split_line[2])
        yx = float(split_line[4])
        yz = float(split_line[6])
        zy = float(split_line[8])
        xz = float(split_line[5])
        zx = float(split_line[3])      
        
        value = (xx + yy + zz)/3
      
        result.append(value)
    
    return np.array(energy), np.array(result)

# 分别计算并存储结果
energy_imag, imag = read_and_process(ImagPartData)
energy_real, real = read_and_process(RealPartData)
tot_eps=np.sqrt(imag**2 + real**2) 

# 创建一个图形框，在其中创建两个子图
fig, ax1 = plt.subplots()

# 第一个子图显示imag（虚部）
ax1.plot(energy_imag, imag, label='Imaginary Part', color='g')
ax1.plot(energy_real, real, label='Real Part', color='r')
ax1.plot(energy_real, tot_eps, label='Total epsilon', linestyle='-',color='black')
#ax1.set_title('Dielectric Function')
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel('Dielectric Function')
ax1.legend()
ax1.grid(True)
plt.show()

```


    
![png](output_19_0.png)
    


### 后处理：其它线性光学性质的计算

单位`absorption_coefficient`的单位是`cm^-1`，单位转换如下：

首先，由于`E = ħω`。并且`1 eV = 1.602×10^-19 J`，`ħ = h / (2π) = 6.626×10^-34 J·s / (2π) = 1.055×10^-34 J·s`。

转换单位：我们需要将`ħ`从`J·s`转换为`eV·s`，那么`ħ`约为`4.135667696×10^-15 eV·s`。

然后我们将`E=ħω`中的`ω`换算成以`eV`为单位的能量，得到`ω = E / ħ`，带入计算得到一个`ω`值，单位为`s^-1`。

光速`c`在SI单位中是`m/s`，为了与`α(ω)`单位中的`cm^-1`匹配，需要将光速`c`从`m/s`转换为`cm/s`。由于`1 m = 100 cm`，`c`大约为`3×10^8 m/s`，转换后`c = 3×10^10 cm/s`。

将以上转换后的`ω`与`c`值带入给定的`α(ω)`表达式中，就可以得到以`cm^-1`为单位的吸收系数。这里的`ω`需要转换为能量单位`eV`，`ε1`和`ε2`则是介电函数的实部和虚部。

此处，能量`E`给定为`1 eV`，那么`ω = E / ħ = (1 eV) / (4.135667696×10^-15 eV·s)` 约等于`2.418×10^14 s^-1`。然后使用`α(ω)`表达式，并将`ω`与`c`（以`cm/s`为单位）的值带入，则可以计算出吸收系数`α(ω)`的值。


```
refraction_index = np.sqrt((np.sqrt(imag**2 + real**2) + real) / 2)
extinction_coefficient = np.sqrt((np.sqrt(imag**2 + real**2) - real) / 2)
absorption_coefficient = np.sqrt(2)*energy_real* np.sqrt((np.sqrt(imag**2 + real**2) - real)) / (3.0*10**8) * (2.418*10**14) 

energy_loss = imag/(real**2 + imag**2)

# 创建一个2x2的图表布局
fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))

# 在每个子图上绘制对应的数据
axs[0, 0].plot(energy_real, refraction_index, 'b-')
axs[0, 0].set_ylabel('Refraction Index')
axs[0, 0].set_xlabel('Energy (eV)')
axs[0, 0].grid(True)

axs[0, 1].plot(energy_real, extinction_coefficient, 'm-')
axs[0, 1].set_ylabel('Extinction Coefficient')
axs[0, 1].set_xlabel('Energy (eV)')
axs[0, 1].grid(True)

axs[1, 0].plot(energy_real, absorption_coefficient, 'y-')
axs[1, 0].set_ylabel(r'Absorption Coefficient (cm$^{-1}$)')
axs[1, 0].set_xlabel('Energy (eV)')
axs[1, 0].set_yscale('log')
axs[1, 0].grid(True)

axs[1, 1].plot(energy_real, energy_loss, 'r-')
axs[1, 1].set_ylabel('Energy Loss')
axs[1, 1].set_xlabel('Energy (eV)')
axs[1, 1].grid(True)

# 为整个图表设定标题
plt.suptitle('Optical Properties vs Energy')
# 调整布局避免子图之间的重叠
plt.tight_layout()

# 显示图示
plt.show()
```


    
![png](output_21_0.png)
    

