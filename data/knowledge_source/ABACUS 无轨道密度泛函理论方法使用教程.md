# ABACUS 无轨道密度泛函理论方法使用教程

> 来源: https://mcresearch.github.io/abacus-user-guide/abacus-ofdft.html

# __ABACUS 无轨道密度泛函理论方法使用教程

**作者：孙亮，邮箱：l.sun@pku.edu.cn**

**审核：陈默涵，邮箱：mohanchen@pku.edu.cn**

**最后更新时间：2023/07/19**

**在Bohrium Notebook上快速学习：**[![](https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg)](<https://nb.bohrium.dp.tech/detail/6416644691>)

# __一、无轨道密度泛函理论背景知识

无轨道密度泛函理论(Orbital free density functional theory, OFDFT)是一种第一性原理模拟方法，相比于 Kohn Sham DFT (KSDFT)，它的优势之一在于O(NlnN)O(Nln⁡N)的算法复杂度，这使得 OFDFT 可以用于上万原子甚至更大体系的电子基态计算，或者大体系、长时间的第一性原理分子动力学等。

目前，OFDFT 已被应用于简单金属、合金、半导体、小分子、温稠密物质等体系。

## __1\. 无轨道密度泛函理论

在 OFDFT 的框架下，体系的总能量泛函可以写为

EOF[ρ]=Ts[ρ]+Eext[ρ]+EH[ρ]+Exc[ρ]+EII.EOF[ρ]=Ts[ρ]+Eext[ρ]+EH[ρ]+Exc[ρ]+EII.

依次为无相互作用动能Ts[ρ]Ts[ρ]，电子-离子相互作用能Eext[ρ]Eext[ρ]，电子-电子相互作用能EH[ρ]EH[ρ]，交换关联能Exc[ρ]Exc[ρ]，离子-离子相互作用能EIIEII，其中ρρ为电荷密度。

为了在粒子数守恒的条件下求解其极小值，利用拉格朗日乘子法，定义 

L=EOF[ρ]−μ(∫ρ(r)dr−N),L=EOF[ρ]−μ(∫ρ(r)dr−N),

可以证明，这里的乘子μμ就是无相互作用体系中最高占据态的能量，即化学势。为了保证密度处处为正，一般对ρ−−√ρ进行优化，因此求LL对ρ−−√ρ的导数，为了方便，定义ϕ(r)=ρ(r)−−−−√ϕ(r)=ρ(r)，则有

δLδϕ=δEOF[ρ]δϕ−2μϕ=0→V(r)ϕ(r)=μϕ(r). (∗)δLδϕ=δEOF[ρ]δϕ−2μϕ=0→V(r)ϕ(r)=μϕ(r). (∗)

其中V(r)=δEOF[ρ]δρ(r)V(r)=δEOF[ρ]δρ(r)为势能。这就是 OFDFT 求解的方程，一般用共轭梯度（CG）法、截断牛顿（TN）法或 L-BFGS 等优化算法求解。目前 ABACUS 中实现了 TN 法和两种 CG 法（Polak-Ribire 形式和 Hager-Zhang 形式），默认采用 TN 法。

ABACUS 基于平面波基矢量，实现了上述流程，可以进行基于 OFDFT 的自洽计算，分子动力学计算，以及结构弛豫。

## __2\. 动能泛函

OFDFT 的精度高度依赖于动能泛函（kinetic energy density functional, 简称 KEDF）的精度，目前 ABACUS 中实现了 Thomas-Fermi (TF) [1], von Weizsäcker (vW) [2], TFλvW [3], Wang-Teter (WT) [4], Luo-Karasiev-Trickey (LKT) [5]共五种动能泛函。

下面我们对这些泛函做简单介绍，并且介绍在 ABACUS 的 `INPUT` 文件中如何设置相关的参数。

### __2.1 Thomas-Fermi KEDF

设置 `INPUT` 文件中的 `of_kinetic tf` 参数

TTF=310(3π2)2/3∫ρ5/3(r)dr,TTF=310(3π2)2/3∫ρ5/3(r)dr,

对均匀电子气精确成立，可用于极高温体系，比如处于温稠密状态的金属。

可通过 `of_tf_weight` 调整其权重，默认为 1。

### __2.2 von Weizsäcker KEDF

设置 `of_kinetic vw`

TvW=−12∫ρ(r)−−−−√∇2ρ(r)−−−−√dr,TvW=−12∫ρ(r)∇2ρ(r)dr,

对单电子、双电子体系（只有一条轨道）严格成立，一般不单独使用。

可通过 `of_vw_weight` 调整其权重，默认为 1。

### __2.3 TFλvW KEDF

设置 `of_kinetic tf+`

T=TTF+λTvW,T=TTF+λTvW,

当λ=19λ=19时就是TF KEDF的二阶梯度展开，一般λ=15λ=15时表现最好。

参数λλ可通过`of_vw_weight`设置，默认为1。

### __2.4 Wang-Teter KEDF

设置`of_kinetic wt`

TWT=310(3π2)2/3∬ρα(r)W(r−r′)ρβ(r′)drdr′+TvW+TTF,TWT=310(3π2)2/3∬ρα(r)W(r−r′)ρβ(r′)drdr′+TvW+TTF,

基于 Lindhard 响应函数推导，在简单金属 Li、Mg、Al 中有着不错的表现。**是 ABACUS 默认采用的 KEDF。**

参数α,βα,β可通过 `of_wt_alpha` 和 `of_wt_beta` 设置，默认值均为5656。

### __2.5 Luo-Karasiev-Trickey KEDF

设置 `of_kinetic lkt`

TLKT=∫τTF1coshasdr+TvW,s=12(3π2)1/3|∇ρ|ρ4/3.TLKT=∫τTF1cosh⁡asdr+TvW,s=12(3π2)1/3|∇ρ|ρ4/3.

可用于简单金属和半导体，计算效率较高，但在简单金属中精度低于 WT KEDF。

参数 aa 可通过 `of_lkt_a` 设置，默认值为 1.3。

## __3\. 局域赝势

由于 OFDFT 中舍弃了单电子轨道，无法采用常用的非局域赝势，如模守恒赝势，而必须采用局域赝势。

目前 ABACUS 支持 BLPS (bulk-derived local pseudopotential)。

> 下载地址： <https://github.com/EACcodes/local-pseudopotentials>  
> 实空间赝势：ABINIT, ABACUS； 倒空间赝势：PROFESS。  
> 赝势生成(需要和 ABINIT 7.0.5 结合)：<https://github.com/EACcodes/BLPSGenerator>  
> 覆盖 Li, Mg, Al, Si, P, Ga, As, In, Sb 九种元素

使用 BLPS 时，需要在 ABACUS 里调整的参数有：

INPUT 中：`pseudo_rcut 16`

STRU 中：赝势种类设置为 `blps`，比如 `Al 26.98 al.lda.lps blps`

# __二、ABACUS 中进行 OFDFT 计算的具体流程

## __1\. 自洽计算

### __1.1 示例

下面是输入文件的示例：

`INPUT` 文件记录 OFDFT 计算所需主要参数
    
    
    INPUT_PARAMETERS
    #Parameters (1.General)
    suffix      example
    calculation scf
    esolver_type    ofdft
    
    symmetry    1   
    pseudo_dir  ../../PP_ORB/
    pseudo_rcut 16
    nspin        1
    
    #Parameters (2.Iteration)
    ecutwfc     60
    scf_nmax    50
    
    #Parameters (3.Basis)
    basis_type  pw
    
    #OFDFT
    of_kinetic  wt
    of_method   tn __

`STRU` 文件记录元素种类、质量、赝势，晶格矢量，原子坐标等信息
    
    
    ATOMIC_SPECIES
    Al 26.98 al.lda.lps blps
    
    LATTICE_CONSTANT
    7.50241114482312  // add lattice constant
    
    LATTICE_VECTORS
    0.000000000000    0.500000000000    0.500000000000    
    0.500000000000    0.000000000000    0.500000000000    
    0.500000000000    0.500000000000    0.000000000000    
    
    ATOMIC_POSITIONS
    Direct
    
    Al
    0
    1
        0.000000000000    0.000000000000    0.000000000000 1 1 1 __

`KPT` 文件（因为 OFDFT 没有电子波函数，所以不需要布里渊区的多个 k 点，Gamma 点就可以）
    
    
    K_POINTS
    0
    Gamma
    1 1 1 0 0 0 __

如上所示，与 KSDFT 的自洽计算相比，OFDFT 自洽计算的输入文件有以下几个区别：

  * INPUT
    * 不需要设置 smearing 和 charge mixing 相关参数，如果设置了也没有关系，这些参数不会影响 OFDFT 计算；
    * 将 `esolver_type` 设置为 `ofdft`；
    * 将 `pseudo_rcut` 设置为 16，以适配 BLPS 赝势。
  * STRU
    * 将赝势种类设置为 `blps`。



做完以上调整后，即可使用默认参数进行 OFDFT 的自洽计算。

下面列举一些其它的重要参数：

  * `of_kinetic`：用于选择动能泛函，可选项有 `tf, vw, tf+, wt, lkt`，默认值为 `wt`，具体介绍见 1.2 节；
  * `of_method`：用于选择优化方法，可选项有 `tn, cg1, cg2`，分别对应截断牛顿法和两种 CG 法（Polak-Ribire 形式和 Hager-Zhang 形式），默认为 `tn`。一般而言，效率上 `tn > cg2 > cg1`；
  * `of_full_pw`：做快速傅里叶变换（FFT）时，是否使用全部的平面波，默认为 `True`。建议打开，可以保证计算的稳定性和精度；
  * `of_full_pw_dim`：控制 FFT 维数的奇偶性，可选项有 `0, 1, 2`，分别表示可奇可偶，保证为奇数，保证为偶数，默认为 `0`。FFT 维数为偶数时，可能导致微小的误差，但一般来说可以忽略。需要注意的是，如果打开了 `nbspline`，则需要设置 `of_full_pw_dim 1`，否则会导致计算不稳定。



### __1.2 注意事项

目前 ABACUS 的 OFDFT 模块并不是十分完善，使用时请注意以下几个注意事项：

  * 目前 OFDFT 不支持 gamma only，因此使用 OFDFT 功能时请关闭 `gamma_only`；
  * 目前 OFDFT 只支持自旋简并，即 `nspin 1` 的计算；
  * 如果使用 PBE 泛函，建议用 `dft_functional XC_GGA_X_PBE+XC_GGA_C_PBE` 调用 Libxc 中的 PBE，否则可能导致计算不稳定。



## __2\. 分子动力学与结构弛豫

ABACUS 中支持使用 OFDFT 作为能量、力和应力的求解器，进行分子动力学模拟与结构弛豫。

与使用 KSDFT 进行分子动力学或结构弛豫相比，使用 OFDFT 时，**不需要对 MD，relax，或 cell-relax 相关参数进行修改，只需要按照 2.1 中的方式，将能量、力和应力的求解器替换为 OFDFT** 。

下面是几个实际的 INPUT 例子：

### __2.1 分子动力学（MD）
    
    
    INPUT_PARAMETERS
    #Parameters (1.General)
    suffix          test
    calculation     md
    esolver_type    ofdft
    pseudo_dir      ../../PP_ORB
    pseudo_rcut     16
    
    #Parameters (2.Iteration)
    ecutwfc     60
    scf_nmax    100
    
    #OFDFT
    of_kinetic  wt
    of_method   tn
    of_full_pw_dim   1
    
    #Parameters (3.Basis)
    basis_type  pw
    
    md_restart  0
    md_type     nvt
    md_nstep    2
    md_dt       0.25
    md_tfirst   58022.52706
    md_dumpfreq 10
    md_tfreq    1.08
    md_tchain   1
    nbspline    10 __

### __2.2 原子结构弛豫（relax）
    
    
    INPUT_PARAMETERS
    #Parameters (1.General)
    suffix          test
    calculation     relax
    esolver_type    ofdft
    pseudo_dir      ../../PP_ORB
    pseudo_rcut     16
    
    #Parameters (2.Iteration)
    ecutwfc     60
    scf_nmax    100
    
    #OFDFT
    of_kinetic  wt
    of_method   tn
    
    #Parameters (3.Basis)
    basis_type  pw
    
    relax_nmax  50 __

### __2.3 晶格弛豫（cell-relax）
    
    
    INPUT_PARAMETERS
    #Parameters (1.General)
    suffix          test
    calculation     cell-relax
    esolver_type    ofdft
    pseudo_dir      ../../PP_ORB
    pseudo_rcut     16
    
    #Parameters (2.Iteration)
    ecutwfc     60
    scf_nmax    100
    
    #OFDFT
    of_kinetic  wt
    of_method   tn
    
    #Parameters (3.Basis)
    basis_type  pw
    
    relax_nmax  50 __

# __三、参考文献

[1] Fermi E. Statistical method to determine some properties of atoms[J]. Rend. Accad. Naz. Lincei, 1927, 6(602-607): 5.

[2] Weizsäcker C F. Zur theorie der kernmassen[J]. Zeitschrift für Physik, 1935, 96(7-8): 431-458.

[3] Berk A. Lower-bound energy functionals and their application to diatomic systems[J]. Physical Review A, 1983, 28(4): 1908.

[4] Wang L W, Teter M P. Kinetic-energy functional of the electron density[J]. Physical Review B, 1992, 45(23): 13196.

[5] Luo K, Karasiev V V, Trickey S B. A simple generalized gradient approximation for the noninteracting kinetic energy density functional[J]. Physical Review B, 2018, 98(4): 041111.

Copyright © mcresearch.gitee.io 2023 all right reserved，powered by Gitbook该文章修订时间： 2025-07-11 15:17:43 
