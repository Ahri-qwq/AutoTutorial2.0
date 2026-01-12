# ASE-ABACUS | 第二章：NEB过渡态计算与ATST-Tools

<p style="color:purple; font-weight:bold">
    【ATST-Tools】，全名【Advanced ASE Transition State Tools for ABACUS and Deep Potential】, 是基于ASE调用ABACUS和DP势函数进行过渡态搜索计算的工作流脚本集，涵盖ASE中包含的NEB，AutoNEB和Dimer方法，已接入高效单端鞍点优化器Sella，并预期未来基于ASE以及ASE-ABACUS接口进一步接入或开发其他的过渡态搜索方法. 
</p>。

在完成这个教程之后，你将能够：
* 学习过渡态理论与过渡态搜索中的基础知识
* 掌握在ATST-Tools工作流组织下，使用ASE调用ABACUS进行基于NEB的双端过渡态搜索计算的前处理，具体计算和后处理工作的方法

欢迎大家提出改进意见！



<div style="color:black; background-color:#FFF3E9; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:0rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2024 @ Authors</i></b><br/>
        <i>作者： 
            <b>
            量子御坂 （Quantum Misaka）
            <a href="mailto:quanmisaka@stu.pku.edu.cn">*** 📨 </a>
            </b>
        </i>
        <br/>
        <i>日期：2024-01-03 </i><br/>
        <i>共享协议：</a>本作品采用<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议</a>进行许可。</i><br/>
        <i>快速开始：点击上方的</i> <span style="background-color:rgb(85, 91, 228); color:white; padding: 3px; border-radius: 5px;box-shadow: 2px 2px 3px rgba(0, 0, 0, 0.3); font-size:0.75rem;">开始连接</span> <i>按钮，选择 <b><u>atst-tools:1.5.0 </u> 镜像</b> 
        和任意配置机型即可开始。
注：由于具体计算的案例都需要比较长的计算时间，本案例不会要求完成具体计算，只涉及预处理和结果分析，故使用c2_m4_cpu即可。如果你想进行具体计算，可使用c16_m32_cpu, c32_m64_cpu等具有较高核数的机型
    </p>
</div>



<a href="https://nb.bohrium.dp.tech/detail/39369325971" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

## 1. 过渡态与过渡态搜索

这部分的内容来自笔者的一些粗浅认知。一些参考资料如下：
- 侯文华等译. Atkins' Physical Chemistry 第十一版中文版 Page. 772. 
- 单斌. 计算材料学. Page. 224

### 1.1 简介

化学是一门研究物质组成，结构及其变化的学科，其主要的研究对象都是物质之间的相互反应，这些反应过程往往由一系列基元反应组成，而这一系列基元反应的发生都离不开反应的过渡态。

**过渡态**由过渡态理论引出，它认为：在A -> B的反应过程中，反应物A会经历一个能量较高的临界状态，并经此亚稳定态转化为产物B，这种临界状态即为过渡态。它是反应能垒的来源，并与反应速率直接相关。

<div>
    <img src="https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/20204/cc00542f53e6459e90e866f6c009e25e/2voIGX7BDbS6pZAOIqRohw.png", alt="reaction_diffusion" width="350" title="reaction_diffusion">
    <p style='font-size:0.8rem; font-weight:bold'> 化学反应图示。黑线为能量曲线，「反应物, 生成物」对应极点，「过渡态」对应鞍点</p>
</div>

在（基于BO近似的）势能面上，这种过渡态位置表现为马鞍点，其能量对坐标的一阶梯度为0，且二阶梯度仅在反应坐标方向为负，其他位置为正，并且也仅在反应坐标方向上存在一个虚频。反应经由这样一个**最小能量路径（Minimum Energy Path， MEP）**进行，从反应物经过渡态生成产物。

<div>
    <img src="https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/2a5aaf89c66d4e3ea3d0b390813beeab/3eed9e80-083d-4f64-adef-7973f5c38bec.jpeg", alt="PES" width="350" title="PES">
    <p style='font-size:0.8rem; font-weight:bold'> 势能面上的极小点（稳定结构），过渡态（马鞍点结构）与最小能量路径 </p>
</div>


基于统计热力学的方法可以推导得到Eyring过渡态理论，该理论指出，基元反应的活化自由能，即过渡态自由能与反应物自由能之差，可以直接与反应速率常数联系起来，即：
$$
k = \kappa \frac{k_B T}{h} e^{\frac{-\Delta G^{\ddagger}}{R T}}
$$
从而基元反应的反应速率可以通过反应物与过渡态的信息计算得到。

过渡态是一种典型的稀有事件，无论是扩散过程、解离过程、还是其他各种均相/多相的化学反应过程，其微观机制均涉及对过渡态的分析，从而过渡态的理论研究与实验检测都是化学研究的重要组成部分。同时，过渡态搜索则是通过理论模拟方法寻找过渡态的关键。过渡态搜索的方法有很多，其中最常用的是NEB和Dimer方法。当然，新的过渡态搜索方法也在不断涌现，尤其是近年来我们很惊喜地看到基于机器学习的生成扩散模型在过渡态搜索与反应生成方面起到了非常惊艳的效果（[戳这里](https://nb.bohrium.dp.tech/detail/1221842760)）

### 1.2 NEB方法与Dimer方法

NEB方法和Dimer方法代表着搜索过渡态的两类不同方法，其一是基于已知反应物和产物的信息去搜索过渡态，即双端（Double-End）过渡态方法，其二是已知过渡态的一个初始猜测去搜索过渡态，即单端（Single-End）过渡态方法。

#### 1.2.1 NEB方法

NEB方法，即Nudge Elastic Band方法，是chain-of-states大类过渡态方法中最著名的一个，其基本原理为：
1. 在初态（第0个态）和末态（第P个态）之间生成P-1个映像（images），这些映像的编号从1到P-1，每个映像代表势能面上连接初末态的一个点，它们同时出现在势能面连接初末态的反应路径上，并在彼此之间通过弹簧力耦合，组成一条映像链（chain）。
2. 对整条链进行优化，使这条映像链拟合到MEP上，其能量最高点即为反应过渡态，而整条映像链则代表反应路径，
3. 在整条链的优化过程中，第i个映像的受力并不是完整的势能梯度+弹簧力，而是二者的投影。如图所示：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/8896a6fc65f94af892baea55f95d7bea/5f0c3af5-799a-4129-8ca4-e09a6e74de84.png)







在目前最为常用的IT-NEB（improved-tangent NEB）方法中，映像在链上的切线方向一般定义为指向能量更高的邻近映像的方向（当然了这个切线方向最终要归一化），即：
$$
\boldsymbol{\tau}_i=\left\{\begin{array}{lll}
\boldsymbol{\tau}_i^+ & \text {if } & E_{i+1}>E_i>E_{i-1} \\
\boldsymbol{\tau}_i^- & \text {if } & E_{i+1}<E_i<E_{i-1}
\end{array}\right.
$$
$$
\boldsymbol{\tau}_i^+ = \boldsymbol{R}_{i+1} - \boldsymbol{R}_{i}, \boldsymbol{\tau}_i^- = \boldsymbol{R}_{i+1} - \boldsymbol{R}_{i}, \hat{\boldsymbol{\tau}_i} = \frac{\boldsymbol{\tau}_i}{\vert \boldsymbol{\tau}_i \vert }
$$
并在能量极值点做平均化处理，使各个点处的切线具有较好的连续性：
$$
\boldsymbol{\tau}_i=\left\{\begin{array}{lll}
\boldsymbol{\tau}_i^{+} \Delta E_i^{\max }+\boldsymbol{\tau}_i^{-} \Delta E_i^{\min } & \text { if } & E_{i+1}>E_{i-1} \\
\boldsymbol{\tau}_i^{+} \Delta E_i^{\min }+\boldsymbol{\tau}_i^{-} \Delta E_i^{\max } & \text { if } & E_{i+1}<E_{i-1}
\end{array}\right.
$$
$$
\Delta E_i^{\max }=\max \left(\left|E_{i+1}-E_i\right|,\left|E_{i-1}-E_i\right|\right)
$$
$$
\Delta E_i^{\min }=\min \left(\left|E_{i+1}-E_i\right|,\left|E_{i-1}-E_i\right|\right) 
$$

与此同时，为了使NEB映像链能更好地定位过渡态，在进行了几步NEB迭代计算之后，将NEB链上能量最高点的映像的受力进行调整，使之不受弹簧力约束，但沿切线方向上的受力反向，即：
$$
\mathbf{F}_{i_{\max }}=-\nabla E\left(\mathbf{R}_{i_{\max }}\right)+\left.2 \nabla E\left(\mathbf{R}_{i_{\max }}\right)\right|_{\|}
$$
即是Climbing-Image NEB方法，简写为CI-NEB方法。这一受力调整不影响计算量，且能明显提升过渡态的解析度。在实际的过渡态搜索中，IT-NEB结合CI-NEB是最常用的方式，能较为高效地定位过渡态与反应路径

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/24407d62f44b463cad88495071f21744/b8c5485f-ca98-46b4-a825-677f2cbe25dd.png)

#### 1.2.2 AutoNEB方法

IT-NEB+CI-NEB已经是非常好的NEB搜索方法了，诚然NEB方法存在一些细节上要注意的问题，比如初末态间原子需要匹配（这是个大难点），纯粹的线性插值容易使初始映像链的结构过于不合理（可以用IDPP插值方法解决）等，但这些都不是非常关键的问题。

一个较为关键的问题在于，NEB计算时映像链的数目都是从最初就确定不变的，虽然这使得每个映像的能量/受力可以单独计算从而明显提升并行度，但这种并行计算的过程中，越靠近过渡态的点其电子结构自洽收敛越困难，也即NEB映像链上各个映像并非同时计算结束，会导致计算资源分布不平衡。同时，传统NEB计算哪怕是采用IDPP等改进的插值方法生成初始NEB链，也会导致各个映像分布过于均衡，对过渡态位置缺乏解析度。

Hammer课题组提出的AutoNEB方法可以缓解上述问题，并明显提升NEB的效率，减小NEB的计算量。其算法示意图如下图所示：



![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/24407d62f44b463cad88495071f21744/ada0347f-ddea-499d-a918-67b71f7fe6a5.png)

AutoNEB的核心是NEB链上的映像会在计算过程中动态地加入，其算法逻辑可以简述为
1. 从初末态+少量初猜构成的NEB链出发，完成（较粗精度的）NEB计算。
2. 定位该NEB链上几何结构/能量差别最为明显的两个映像，在这两个映像间通过插值方法加入一个新的映像
3. 以该新的映像为中心，继续进行（较粗精度的）NEB计算（该过程中有些原有映像是不会参与计算的），得到一条新的连接初末态的NEB链。
4. 重复上述两步，直至NEB链上映像个数达到设定阈值
5. 以能量最高点为中心，最后进行一步（目标精度的）CI-NEB计算，得到最终的NEB链



AutoNEB与传统NEB的对比可从下图看出：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/24407d62f44b463cad88495071f21744/7a65d809-25f5-44de-aed2-9280522569f7.png)

通过AutoNEB这样的动态NEB过程，新加入的映像将对过渡态和反应路径有更好的解析度，并且显著缓解NEB的计算资源不平衡问题，以及NEB的计算量。

这两种方法都内置在了ASE中，并可在ATST-Tools中方便地使用

#### 1.2.3 Dimer等单端过渡态方法

该部分将在ASE-ABACUS教程第三章详细讲解：
https://bohrium.dp.tech/notebooks/29581597682

## 2. ATST-Tools及其使用

此前，我们介绍过ASE-ABACUS接口的一些基础使用（[回顾点这个](https://nb.bohrium.dp.tech/detail/6516485694)）。本期我们补上过渡态搜索的坑。同时，本教程还会涉及基于ASE进行过渡态计算信息收集，针对过渡态和反应初末态的振动校正等操作。

### 2.1 ASE中的过渡态搜索

ASE中内置了三种过渡态搜索方法
1. 约束优化法，在一维反应坐标上进行约束与爬升，放开其他自由度优化，逐步逼近过渡态。这一方法简单且直觉性强，但难以处理复杂过渡态。该方法可以通过`ase.optimize.climbfixinternals.BFGSClimbFixIntgernals`模块调用
2. NEB方法及其变体，可通过`ase.mep.neb`调用。其中，AutoNEB方法需要通过`ase.mep.autoneb`调用
3. Dimer方法，可通过`ase.mep.dimer`调用









### 2.2 过渡态搜索

#### 2.2.1 基础

在熟练了ASE和ASE-ABACUS接口的代码逻辑之后，直接写一个python脚本，调用ABACUS等计算组件，基于ASE进行过渡态搜索计算，并不是一件非常困难的事情。但我们或许有更好的办法。

在ATST-Tools中，与过渡态搜索相关的前处理，具体计算与后处理工作被合理模块化与耦合起来，我们只需要进行简单的编辑设置，就可以在服务器等命令行环境下快速开展过渡态搜索工作。

ATST-Tools可以通过如下链接访问：
https://github.com/QuantumMisaka/ATST-Tools 。其README中有详细的环境配置和使用方法介绍。

在本notebook的推荐镜像中已经装有ATST-Tools及其依赖环境。我们可以进入相关目录查看


```
cd /opt/ATST-Tools
```

    /opt/ATST-Tools
    


```
ls
```

    README.md  [0m[01;34mdimer[0m/     [01;34mimg[0m/  [01;34mrelax[0m/  [01;34msource[0m/
    [01;34mase-dp[0m/    [01;34mexamples[0m/  [01;34mneb[0m/  [01;34msella[0m/  [01;34mvibration[0m/
    

在ATST-Tools的母目录中我们可以看到这些件夹：
- `neb`：使用ATST-Tools进行NEB和其他改进版NEB计算时所需的脚本。
- `dimer`：使用ATST-Tools进行Dimer计算所需的脚本。
- `sella`：使用ATST-Tools进行Sella计算所需的脚本。
- `vibration`：进行振动计算和自由能校正时所需的程序，此处使用ASE进行振动计算。
- `relax`：使用ASE-ABACUS进行结构优化计算的脚本。
- `source`：ATST-Tools整合的Python库，包括用ase-abacus接口进行NEB, Dimer等运算时的工作流配置，以及自行微调/添加的过渡态计算方法等。
- `examples`： ATST-Tools进行计算的案例。这些案例可以在已经配好的服务器/Bohrium镜像下，在简单调整之后直接运行
- `ase-dp`：ATST-Tools基于DP势作为计算后端进行过渡态计算的工作流，均是独立脚本，且处于快速迭代阶段，新特性往往现在这里面更新，然后再同步到其他以ABACUS为计算后端的脚本中。

在实际使用ATST-Tools之前，我们需要让Python能识别到`source`目录下的Python库，在终端运行：
```shell
export PYTHONPATH=/opt/ATST-Tools/source:$PYTHONPATH
```
或是将上一行加入到 *~/.bashrc* 完成这个过程，你的镜像应当已经配置好了。

在README中有对ATST-Tools非常详细的介绍，我们可以看一看。这个README非常的长，我们设置为了滚动模式。


```
cat README.md
```

    # ATST-Tools
    Advanced ASE Transition State Tools for ABACUS and Deep-Potential, including:
    - NEB, including CI-NEB, IT-NEB and others.
    - Serial Dynamic NEB (DyNEB) calculation.
    - AutoNEB: an automatic NEB workflow.
    - Single-End TS search: Sella, Dimer.
    - IRC analysis by Sella.
    - Double-to-single (D2S) TS workflow: neb2dimer, neb2sella.
    - Vibration analysis and ideal gas thermochemistry analysis.
    
    Version 1.5.0
    
    Copyright @ QuantumMisaka from TMC-PKU & AISI
    
    ## Update log from 1.4.0 to 1.5.0
    - Add internal reaction coordination (IRC) calculation by using Sella packages.
    - Add *neb2dimer_dp.py*, *neb2sella_dp.py*, *relax_dp.py* for deepmd usage, including tf and torch version, and all D2S method are using pymatgen IDPP.
    - Change default relaxation optimizer from BFGS to QuasiNewton (BFGSLineSearch).
    - change the way to use `AbacusProfile` to march the newest version of ase-abacus interface, which use different command format.
    - Since constraints read-in problem in ASE-ABACUS interface is fixed [#28](https://gitlab.com/1041176461/ase-abacus/-/issues/28), the constraints is the two scripts above is modified, make Dimer and Sella useful.
    - Single-End TS search scripts update using Dimer method and Sella packages.
    - Default *neb_make.py* change to pymatgen version below.
    - *neb_make_pymatgen.py* scripts by @MoseyQAQ, which use IDPP of pymatgen to do NEB initial guess, avoiding IDPP problem from ASE including edge-crossing problem in NEB guess generation.
    - Support DeepPotential and DPA-2 Potential usage scripts in `ase-dp` directory.
    - New function: `neb2vib` method to use NEB chain information to do partial vibrational analysis.
    - Thermochemistry analysis for ideal gas system.
    - Examples Update
    
    ## Dependencies:
    - [ASE](https://wiki.fysik.dtu.dk/ase/about.html), but you should install ASE by [ASE-ABACUS interface](https://gitlab.com/1041176461/ase-abacus) as a separate ASE package for ABACUS usage.
    - [ABACUS](https://github.com/deepmodeling/abacus-develop), one can install ABACUS by ABACUS toolchain [ABACUS](https://github.com/deepmodeling/abacus-develop/tree/develop/toolchain), or refer to [ABACUS-docs](https://abacus.deepmodeling.com/en/latest/)
    - [pymatgen](https://pymatgen.org/) and [pymatgen-analysis-diffusion](https://github.com/materialsvirtuallab/pymatgen-analysis-diffusion) in the usage of new `neb_make.py` script and D2S TS method, which can be installed by `pip install pymatgen pymatgen-analysis-diffusion`
    - [Sella](https://github.com/zadorlab/sella) if one wants to use Sella method for Single-End TS search. which can be installed by `pip install sella`
    - [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) if one wants to run NEB images relaxation in parallel. The installation of GPAW need some efferts, and one can refer to [GPAW installation](https://wiki.fysik.dtu.dk/gpaw/install.html) for more details.
    - [deepmd-kit](https://github.com/deepmodeling/deepmd-kit) or [DPA-2](https://zenodo.org/records/10428497) if one wants to use Deep-Potential or DPA-2 potential, one can refer to [deepmd-docs](https://docs.deepmodeling.com/projects/deepmd/en/master/)
    
    
    Notice: GPAW and ABACUS should be dependent on same MPI and libraries environments. 
    > For instance, if your ABACUS is installed by Intel-OneAPI toolchain, your GPAW should NOT be dependent on gcc-toolchain like OpenMPI and OpenBLAS.
    
    ATST-Tools is Under actively development, please let me know if any problem occurs.
    
    ## Tutorials
    - One can get ASE-ABACUS usage from this [Bohrium Notebook 1](https://bohrium.dp.tech/notebooks/6516485694)
    - One can get NEB tutorial from this [Bohrium Notebook 2](https://nb.bohrium.dp.tech/detail/39369325971)
    - One can get Single-End TS search tutorial from this [Bohrium Notebook 3](https://bohrium.dp.tech/notebooks/29581597682)
    - D2S TS exploration tutorial is coming soon.
    
    ## Workflow
    
    ![ATST-NEB-workflow](./img/ATST-workflow.png)
    
    ## Workflow libraries and setting
    All workflow library files and re-constructed ASE libraries will be put in `./source` directory. including:
    ```bash
    source
    ├── abacus_autoneb.py
    ├── abacus_dimer.py
    ├── abacus_neb.py
    ├── my_autoneb.py
    └── neb2vib.py
    ```
    Before use running scripts, you should add these libraries into your PYTHONPATH:
    ```bash
    export PYTHONPATH=/path/to/source:$PYTHONPATH
    ```
    
    There are also other workflow in ATST-Tools
    - AutoNEB TS exploration
    - Double-to-Single (D2S) TS exploration
    
    ## Developing
    - [ ] More fiexible options for NEB, Dimer and AutoNEB, like full properties in trajectory file, and fiexibly utilize SCF wavefunction/charge output files from previous calculation.
    - [ ] Move workflow parts of Deep-Potential and DPA-2 potential to `source` directory
    - [ ] Parallel NEB calculation for D2S method and DP usage
    - [ ] [dflow](https://github.com/dptech-corp/dflow) version of ATST-Tools workflow
    - [ ] Bond soft scanning and automatic reaction locating (proposed by [CLAM workflow](https://github.com/lalaheihaihei/catalyticLAM))
    - [x] Test for Sella usage to get best performance
    - [x] Checkout the problem in Dimer-ABACUS calculation
    - [x] Inplement `neb2vib` method in  `neb2dimer` workflow for DPA-2 (and ABACUS)
    - [x] Other TS method usage, like [Sella](https://github.com/zadorlab/sella)
    
    
    ## NEB workflow
    
    ### Method
    - For serial NEB calculation, DyNEB, namely dynamic NEB method `ase.mep.neb.DyNEB` is for default used.
    - For parallel NEB calculation, `ase.mep.neb.NEB` traditional method is for default used.
    - The Improved Tangent NEB method `IT-NEB` and Climbing Image NEB method `CI-NEB` in ASE are also default used in this workflow, which is highly recommended by Sobervea. In `AutoNEB`, `eb` method is used for default, but Improved Tangent method is also recommended.
    - Users can change lots of parameter for different NEB setting. one can refer to [ASE NEB calculator](https://wiki.fysik.dtu.dk/ase/ase/neb.html#module-ase.neb) for more details: 
    - The workflow also support use `AutoNEB` method in ASE。You can view AutoNEB method in paper below. Also. One can refer to [AutoNEB](https://wiki.fysik.dtu.dk/ase/ase/neb.html#autoneb) to view it. 
    > E. L. Kolsbjerg, M. N. Groves, and B. Hammer, J. Chem. Phys, 145, 094107, 2016. (doi: 10.1063/1.4961868)
    
    The AutoNEB method in ASE lies in `ase.mep.autoneb.AutoNEB` object, which will do NEB calculation in following steps:
    1. Define a set of images and name them sequentially. Must at least have a relaxed starting and ending image. User can supply intermediate guesses which do not need to have previously determined energies (probably from another
    NEB calculation with a lower level of theory)
    1. AutoNEB will first evaluate the user provided intermediate images
    2. AutoNEB will then add additional images dynamically until n_max is reached
    3. A climbing image will attempt to locate the saddle point
    4. All the images between the highest point and the starting point are further relaxed to smooth the path
    5. All the images between the highest point and the ending point are further relaxed to smooth the path
    
    Step 4 and 5-6 are optional steps. Note that one can specify different `fmax` for CI-NEB in step 4 compared with other NEB calculation, which can be set as `fmax=[fmax1, fmax2]` in `AutoNEB` object.
    
    > Notice: in surface calculation and hexangonal system, the vaccum and c-axis should be set along y-direction but not z-direction, which is much more efficient for ABACUS calculation.
    
    
    ### Usage
    #### Basic NEB
    The NEB workflow is based on 3 main python scripts and 1 workflow submit script. Namely:
    
    - `neb_make.py` will make initial guess for NEB calculation, which is based on ABACUS (and other calculator) output files of initial and final state. This script will generate `init_neb_chain.traj` for neb calculation. Also, You can do continuation calculation by using this script. You can get more usage by `python neb_make.py`. 
    - `neb_run.py` is the key running script of NEB, which will run NEB calculation based on `init_neb_chain.traj` generated by `neb_make.py`. This script will generate `neb.traj` for neb calculation. Users should edit this file to set parameters for NEB calculation. sereal running can be performed by `python neb_run.py`, while parallel running can be performed by `mpirun gpaw python neb_run.py`.
    When running, the NEB trajectory will be output to `neb.traj`, and NEB images calculation will be doing in `NEB-rank{i}` directory for each rank which do calculation of each image. 
    - `neb_post.py` will post-process the NEB calculation result, which will based on `neb.traj` from neb calculation. This script will generate nebplots.pdf to view neb calculation result, and also print out the energy barrier and reaction energy. You can get more usage by `python neb_post.py`. Meanwhile, users can also view result by `ase -T gui neb.traj` or `ase -T gui neb.traj@-{n_images}:` by using ASE-GUI
    - `neb_submit.sh` will do all NEB process in one workflow scripts and running NEB calculation in parallel. Users should edit this file to set parameters for NEB calculation. Also this submit script can be used as a template for job submission in HPC. the Default setting is for `slurm` job submission system.
    
    #### AutoNEB method 
    In ATST-Tools, the AutoNEB method can be easily used by the following scripts
    - `autoneb_run.py` is the key running script for `AutoNEB` method, which is like `neb_run.py` but the NEB workflow in `AutoNEB` is enhanced and the I/O logic have some difference. Users can use it with `mpirun gpaw python autoneb_run.py` by existing `init_neb_chain.traj` which can only contain initial and final state or contain some initial-guess.
    - `autoneb_submit.sh` will do all AutoNEB process in one workflow and running AutoNEB calculation in parallel. Users should edit this file to set parameters for AutoNEB calculation. Also this submit script can be used as a template for job submission in HPC. the Default setting is for `slurm` job submission system.
    - `neb_make.py` and `neb_post.py` can be used for `AutoNEB` method, but the workflow have slight difference. 
    
    #### Running
    Users can run NEB each step respectively: 
    1. `python neb_make.py [INIT/result] [FINAL/result] [n_max]` to create initial guess of neb chain
       1. Also You can use `python neb_make.py -i [input_traj_file] [n_max]` to create initial guess from existing traj file, which can be used for continuation calculation.
    2. `python neb_run.py` or `mpirun -np [nprocs] gpaw python neb_run.py` to run NEB calculation
    3. `python neb_post.py neb.traj [n_max]` to post-process NEB calculation result
    
    Users can run AutoNEB each step respectively:
    1. `python neb_make.py -i [INIT/result] [FINAL/result] -n [nprocs]` to create initial guess of neb chain
    2. `mpirun -np [nprocs] gpaw python autoneb_run.py` to run AutoNEB calculation
    3. `python neb_post.py --autoneb run_autoneb???.traj` to post-process NEB calculation result
    
    Also, user can run each step in one script `neb_submit.sh` by `bash neb_submit.sh` or `sbatch neb_submit.sh`. AutoNEB scripts usage is like that. 
    
    > Notice: Before you start neb calculation process, make sure that you have check the nodes and cpus setting and other setting like n_max, constraints and initial magnetic moments in `*neb_submit.sh` and `*neb_run.py` to make sure that you can reach the highest performance and reach the simulation result you want !!!   
    
    #### Visualize
    You can simply use ASE-GUI to have a view of NEB or AutoNEB trajectory.
    
    For NEB, you can view all running trajectory by 
    ```bash
    ase -T gui neb.traj 
    ```
    You can also view the last 10 images by
    ```bash
    ase -T gui neb.traj@-10:
    ```
    For AutoNEB, the most recent NEB path can always be monitored by:
    ```bash
    ase -T gui -n -1 run_autoneb???.traj
    ```
    
    #### Continuation calculation for NEB
    If NEB or AutoNEB is break down somehow, you can do continuation calculation based on saved trajectory files and ATST-Tools scripts.
    
    For NEB, you can simply:
    ```bash
    python neb_make.py -i neb.traj -n [n_max] [fix and mag information]
    ```
    to generate `init_neb_chain.traj` for continuation calculation. You can also `python neb_post.py neb.traj` to generate the latest neb band `neb_latest.traj` and do continuation calculation by `python neb_make.py -i neb_latest.traj [n_max]`. note that `n_max = n_image - 2`
    
    For AutoNEB, you need to get `neb_latest.traj` in a more compicated way, especially for the first NEB step in AutoNEB running:
    ```bash
    python traj_collect.py ./AutoNEB_iter/run_autoneb???iter[index].traj
    ```
    to generate `collection.traj` from certain index (like 006) stage of AutoNEB calculation. Another way to do the same thing is:
    ```bash
    python traj_collect.py ./run_autoneb???.traj
    ```
    or
    ```bash
    python traj_collect.py ./AutoNEB_iter/run_autoneb???iter00[i].traj
    ```
    or
    ```bash
    python traj_collect.py ./AutoNEBrun_rank?/STRU
    ```
    When `collection.traj` is gotten, one can do
    ```bash
    python neb_make.py -i collection.traj [n_max] [fix and mag information]
    ```
    to generate `init_neb_chain.traj` for continuation calculation.
    
    If one want to continue calculation with the interrupted autoneb, one can:
    ```bash
    python traj_collect.py --no-calc ./run_autoneb???.traj 
    ```
    and `--no-calc` noted can be anywhere.
    
    > Note: Linux shell will automatically detect and sort number of index, so you will not be worried about using format like `run_autoneb???iter005.traj`, the consequence will be right, for example:
    ```bash
     test> ll run_autone*.traj
    -rw-r--r-- 1 james james 6.0K Nov 24 20:35 run_autoneb000.traj
    -rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb001.traj
    -rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb002.traj
    -rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb003.traj
    -rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb004.traj
    -rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb005.traj
    -rw-r--r-- 1 james james 531K Nov 24 20:35 run_autoneb006.traj
    -rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb007.traj
    -rw-r--r-- 1 james james 6.4K Nov 24 20:35 run_autoneb008.traj
    -rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb009.traj
    -rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb010.traj
    -rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb011.traj
    -rw-r--r-- 1 james james 6.5K Nov 24 20:35 run_autoneb012.traj
    -rw-r--r-- 1 james james 531K Nov 24 20:37 run_autoneb025.traj
    ```
    
    
    #### Other scripts
    Because ATST is originally based on ASE, the trajectory file can be directly read, view and analysis by `ase gui` and other ASE tools. Abide by `neb_make.py` and `neb_post.py`, We also offer some scripts to help you:
    - `neb_dist.py`: This script will give distance between initial and final state, which is good for you to check whether the atoms in two image is correspondent, and is also a reference for setting number of n_max
    - `traj_transform.py`: This script can transfer traj files into other format like `extxyz`, `abacus`(STRU), `cif` and so on (coming soon). Also if user specify `--neb` option, this script will automatically detect and cut the NEB trajectory when doing format transform. This script will be helpful for analysis and visualization of NEB trajectory.
    - `traj_collect.py`: This script can collect structure files into a trajectory file, which is specifically used for NEB continuation calculation.
    
    ## Single-End TS search
    
    Contrary to NEB, the single-end TS search is based on the exploration of saddle point by using the gradient and Hessian (or only Hessian eigenmode) information of points in PES, which have good efficiency if one have an approximate TS information. ATST-Tools support Dimer method and Sella method.
    
    ### Sella workflow
    The Sella workflow is based on 1 main python scripts and 1 submit script, namely:
    - `sella_run.py`  is the key running script of Dimer calculation, which will run Sella calculation based on supplying `STRU` files for initial state of Sella calculation. This script will generate `run_sella.traj` for Sella calculation trajectory. Users should edit this file to set parameters for Sella calculation, and run Sella calculation by `python sella_run.py`. When running, any Sella images calculation will be doing in `ABACUS` directory.
    - `sella_submit.sh` will do Sella workflow in one scripts. The Default setting is for `slurm` job submission system.
    
    Also, Sella package can be used for IRC calculation, one can refer to `sella_IRC.py` for more details.
    
    For more detail on Sella principle and usage, one can refer to [Sella](https://github.com/zadorlab/sella) and [Sella-wiki](https://github.com/zadorlab/sella/wiki)
    
    
    ### Dimer workflow
    
    The Dimer workflow is based on 2 main python scripts and 1 workflow submit script, namely:
    - `neb2dimer.py` can be used by `python neb2dimer [neb.traj] ([n_max])`, which will transform NEB trajetory `neb.traj` or NEB result trajectory `neb_result.traj` to Dimer input files,  including:
    - - `dimer_init.traj` for initial state of Dimer calculation, which is the highest energy image, namely, TS state. 
    - - `STRU` files, representing the initial state of Dimer calculation, stored as `STRU` files for ABACUS calculation.
    - - `displacement_vector.npy` for displacement vector of Dimer calculation, which will be generated from position minus of the nearest image before and after TS point, and be normalized to 0.01 Angstrom. 
    - `dimer_run.py` is the key running script of Dimer calculation, which will run Dimer calculation based on `dimer_init.traj` and `displacement_vector.npy` generated by `neb2dimer.py` or based on other setting. This script will generate `run_dimer.traj` for Dimer calculation trajectory. Users should edit this file to set parameters for Dimer calculation, and run Dimer calculation by `python dimer_run.py`. When running, any Dimer images calculation will be doing in `ABACUS` directory.
    - `dimer_submit.sh` will do Dimer workflow in one scripts. The Default setting is for `slurm` job submission system.
    
    
    ## Double-to-single (D2S) TS exploration
    
    ATST-Tools also offer a double-to-single (D2S) TS exploration method, which is named as `neb2dimer_abacus.py` and `neb2sella_abacus.py` in `dimer` and `sella` directories respectively, which run a rough NEB first and use the maximum information for initial guess of running single-ended method like dimer or Sella, obtaining a high-efficiency TS exploration. 
    
    In D2S workflow, Dynamic NEB acceleration is used for NEB calculation for serial usage for faster rough NEB calculation.
    
    This workflow is only for serial NEB calculation, but is efficient enough for TS exploration. The parallel NEB acceleration is in considering.
    
    ## Vibration Analysis 
    The vibration analysis is based on `ase.vibrations.Vibrations` object, which can be used by `python vib_analysis.py` to do vibration analysis by finite displacement method for initial, final and transition state. The result will be printed out  and saved in `running_vib.out` file. All force matrix for displaced and normal mode will also be saved and printed.
    
    Also, thermodynamic analysis will be performed based on `ase.thermochemistry.HarmonicThermo` object based on vibration analysis result and specified temperature.
    
    ## Relaxation
    Relaxation method offer by ASE can be used by scripts from `relax` directory, which use ABACUS as SCF calculator. 
    > Notes: QuasiNewton method in ASE is BFGSLineSearch, which is default usage in ATST-Tools and can have robust structural relaxation ability.
    
    ## ase-dp
    There are also scripts for Deep-Potential and DPA-2 potential usage in `ase-dp` directory for TS exploration and structural relaxation, some newest examples can be used by `neb2dimer_dp.py`, `neb2sella_dp.py` and `relax_dp.py` for Deep-Potential and DPA-2 potential usage. 
    
    Parallel NEB version for ase-dp is also in developing.
    
    ## Examples
    - Li-diffu-Si: Li diffusion in Si, very easy example for serial and parallel NEB calculation about diffusion system. Note that ATS-Tools will tune the parameters towards TS exploration, for those who want to do diffusion calculation, one should check the effectiveness of the parameters by yourself.
    - H2-Au111: H2 dissociation on Au(111) surface. which will have NEB, serial DyNEB, AutoNEB, Dimer, Sella and D2S example. The barrier is around 1.1 eV consistent with existing paper and calculation result. IRC results shows the reaction path is correct.
    - CO-Pt111 : CO dissociation on Pt(111) surface. which have high barrier as 1.5 eV and including diffusion part. AutoNEB will always fail due to the low starting image. By performing Dimer / Sella calculation, or use the D2S method, the TS can be explored in an efficient and accurate way. By using Sella package, one can also get the IRC of C-O dissociaition reaction process.
    - Cy-Pt_graphene: Cyclohexane dehydrogenation on Pt-doped graphene surface. The barrier is around 1.3 eV. Noted that the `IT-NEB` result is wrong, but which is consistent to the result in VTST-Tools when using 4 image to do IT-NEB calculation. Dimer calculation also get wrong result in this case, while Sella perform good result.
    
    More examples is welcomed from users. 
    
    ## Notices
    ### Property Loss in Trajectory 
    Some property should be get via specific way from trajectory files, and some will be lost in trajetory files, 
    - Stress property will not be stored in trajetory file
    - In NEB calculation, the Force property for fixed atoms and Stress property will NOT be stored in trajectory file, one should get it by `get_force(apply_constraint=False)`.
    - in Dimer calculation, the Energy, Forces and Stress property will be stored in trajetory file after specified which properties need to be stored. (Sella calculation should be likely) (The most easiest lost information is the stress information)
    - in AutoNEB calculation, all property in processing trajectory will be stored in AutoNEB_iter directory, but in the result `run_autoneb???.traj`, the forces and stress information will be lost.
    
    
    

`examples`文件夹中提供了几个简单的例子的输入与输出。这些例子在具体使用ATST-Tools并调用ABACUS进行计算时会耗费较长的运算时间。一般来说，在16个物理核上，`Li-diffu-Si`的NEB例子会花费1-2小时左右的时间，`H2-Au111`以及`Cy-Pt@graphene`会耗费12-24小时左右的时间，感兴趣的小伙伴可以自行拉取ATST-Tools源码并运行
```
git clone https://github.com/QuantumMisaka/ATST-Tools
```
这里我们直接看看



```
cd /opt/ATST-Tools/examples/
```

    /opt/ATST-Tools/examples
    


```
ls
```

     [0m[01;34mCH4-HOAuHOPd-ZSM5[0m/   [01;34mCO-Pt111[0m/  [01;34m'Cy-Pt@graphene'[0m/   [01;34mH2-Au111[0m/   [01;34mLi-diffu-Si[0m/
    

#### 2.2.2 案例一：扩散过程的串并行NEB计算

我们先看一个简单的例子，`Li-diffu-Si`跑的是Li在硅晶胞中的扩散。进入到文件夹中，其中有两个算例文件夹`dyneb`, `neb`和数据文件夹`data`。`dyneb`是ASE中实现了的一种针对串行NEB进行了优化的NEB算法，`neb`中则是传统的并行NEB方法。由于该过程仅设计Li原子的扩散，因而在初末态之间插入3个映像，即可很好地得到扩散过渡态结构与扩散路径。


```
cd /opt/ATST-Tools/examples/Li-diffu-Si
```

    /opt/ATST-Tools/examples/Li-diffu-Si
    


```
ls
```

    [0m[01;34mdata[0m/  [01;34mdyneb[0m/  [01;34mneb[0m/
    

这个例子如果采用3个映像插值，在`c32_m64_cpu`机器来运行，并采用串行优化的NEB方法，可以在1-2个小时内运行完毕。具体设置方法如下。


```
mkdir dyneb-example
```

    mkdir: cannot create directory ‘dyneb-example’: File exists
    


```
cp -r data/* dyneb-example
```


```
cd dyneb-example
```

    /opt/ATST-Tools/examples/Li-diffu-Si/dyneb-example
    


```
ls
```

    [0m[01;34mFINAL[0m/  Li_ONCV_PBE-1.2.upf        Si_ONCV_PBE-1.2.upf
    [01;34mINIT[0m/   Li_gga_8au_100Ry_4s1p.orb  Si_gga_8au_100Ry_2s2p1d.orb
    

通过上述方法准备好数据文件，包括：
- 初末态的计算结果
- 计算所用的赝势和轨道文件
这一准备过程一般需要用户自己完成。注意NEB计算的初末态是需要提前结构优化好的。

之后，通过ATST-Tools中的`neb_make.py`脚本进行插值，在Bohrium上我们需要通过`python3`来调用python。这一脚本的使用方法可以通过直接运行查看。不难发现，该脚本支持在加入初猜的时候显式地加入基底约束，磁矩初猜，过渡态初猜等信息，并允许加载已有的初猜NEB链，灵活性非常强。

该脚本经历过社区协作者的一次更新迭代，采用了Pymatgen的IDPP插值，能处理更为复杂的插值情况。

我们可以先通过-h设置来看一看该脚本的输入需求


```
! python3 /opt/ATST-Tools/neb/neb_make.py -h
```

    usage: neb_make.py [-h] -n N [-f FORMAT] -i INPUT INPUT [-m {IDPP,linear}]
                       [-o O] [-sort_tol SORT_TOL] [--fix FIX] [--mag MAG]
    
    Make input files for NEB calculation
    
    options:
      -h, --help            show this help message and exit
      -n N                  Number of images
      -f FORMAT, --format FORMAT
                            Format of the input files, default is abacus-out
      -i INPUT INPUT, --input INPUT INPUT
                            IS and FS file
      -m {IDPP,linear}, --method {IDPP,linear}
                            Method to generate images
      -o O                  Output file
      -sort_tol SORT_TOL    Sort tolerance for matching the initial and final
                            structures, default is 1.0
      --fix FIX             [height]:[direction] : fix atom below height
                            (fractional) in direction (0,1,2 for x,y,z), default
                            None
      --mag MAG             [element1]:[magmom1],[element2]:[magmom2],... : set
                            initial magmom for atoms of element, default None
    

然后，基于已经准备好了的初末态数据进行插值


```
! python3 /opt/ATST-Tools/neb/neb_make.py -i ./INIT/OUT.ABACUS/running_scf.log ./FINAL/OUT.ABACUS/running_scf.log -n 3
```

    Reading files: ./INIT/OUT.ABACUS/running_scf.log and ./FINAL/OUT.ABACUS/running_scf.log
    Generating path, number of images: 3, sort_tol: 1.0
    Optimizing path using IDPP method
    Writing path: init_neb_chain.traj,Number of images: 5
    

通过上述方法我们能得到`init_neb_chain.traj`。这一二进制轨迹格式是ASE官方指定的轨迹存储格式，并能直接被ATST-Tools利用。


```
ls
```

    [0m[01;34mFINAL[0m/               Li_gga_8au_100Ry_4s1p.orb    init_neb_chain.traj
    [01;34mINIT[0m/                Si_ONCV_PBE-1.2.upf
    Li_ONCV_PBE-1.2.upf  Si_gga_8au_100Ry_2s2p1d.orb
    

于是我们完成了NEB的插值过程，接下来可以利用ATST-Tools开展NEB计算。这里分别展示串行DyNEB计算和并行NEB计算的过程，这些计算默认并不直接进行。

一般来说，在准备好初猜之后，将`ATST-Tools/neb`目录中的`neb_run.py`拷贝到目录下，微调其设置，然后即可运行该脚本来跑NEB计算。


```
cp /opt/ATST-Tools/neb/neb_run.py .
```


```
cat neb_run.py
```

    # JamesMisaka in 2023-11-27
    # Run NEB calculation on NEB images by ASE-ABACUS
    # part of ATST-Tools scripts
    
    from ase.optimize import FIRE, BFGS
    from ase.io import read, write
    from abacus_neb import AbacusNEB
    
    # setting
    mpi = 8
    omp = 4
    fmax = 0.05  # eV / Ang
    neb_optimizer = FIRE # suited for CI-NEB
    neb_directory = "NEBrun"
    algorism = "improvedtangent" # IT-NEB is recommended
    climb = True
    dyneb = False  
    parallel = True
    k = 0.10 # eV/Ang^2, spring constant
    init_chain = "init_neb_chain.traj"
    abacus = "abacus"
    #lib_dir = "/lustre/home/2201110432/example/abacus"
    lib_dir = ""
    pseudo_dir = f"{lib_dir}/PP"
    basis_dir = f"{lib_dir}/ORB"
    # default pp and basis is supported by ase-abacus interface, need to check usage
    pp = {
          'C':'C_ONCV_PBE-1.0.upf',
          'H':'H_ONCV_PBE-1.0.upf',
          'Pt':'Pt_ONCV_PBE-1.0.upf',
          }
    basis = {
             'C': 'C_gga_7au_100Ry_2s2p1d.orb',
             'H': 'H_gga_6au_100Ry_2s1p.orb',
             'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb',
             }
    kpts = [2, 1, 2]
    parameters = {
        'calculation': 'scf',
        'nspin': 2,
        'xc': 'pbe',
        'ecutwfc': 100,
        'dft_functional': 'pbe',
        'ks_solver': 'genelpa',
        'symmetry': 0,
        'vdw_method': 'd3_bj',
        'smearing_method': 'gaussian',
        'smearing_sigma': 0.001,
        'basis_type': 'lcao',
        'mixing_type': 'broyden',
        'scf_thr': 1e-6,
        'scf_nmax': 100,
        'kpts': kpts,
        'pp': pp,
        'basis': basis,
        'pseudo_dir': pseudo_dir,
        'basis_dir': basis_dir,
        'init_wfc': 'atomic',
        'init_chg': 'atomic',
        'cal_force': 1,
        'cal_stress': 1,
        'out_stru': 1,
        'out_chg': 0,
        'out_mul': 0,
        'out_wfc_lcao': 0,
        'out_bandgap': 0,
        'efield_flag': 1,
        'dip_cor_flag': 1,
        'efield_dir': 1,
    }
    
    
    if __name__ == "__main__":
    # running process
    # read initial guessed neb chain
        init_chain = read(init_chain, index=':')
        neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                        directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                        algorism=algorism, k=k, dyneb=dyneb)
        neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
    
    
    

需要编辑的地方主要有；
- 赝势与轨道文件所在目录`pseudo_dir`, `basis_dir`
- 目标体系的赝势与轨道信息
- INPUT参数设置和K设置
- 其他关于NEB运行相关的设置。

如果直接在notebook中运行，默认跑的是串行版本的NEB计算，并行计算需要利用mpirun调用gpaw的python，需要在命令行进行。

以下是串行运行NEB的脚本，采用了Dynamic的串行优化方法。脚本第一行将脚本写入到`dyneb_run.py`文件中，并可在命令行调用。由于该程序不对NEB本身做并行，用户可以将该行去掉，并在c32_m64_cpu机器中直接运行该代码。


```
%%writefile dyneb_run.py

# erase the row above and you can run it in notebook

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from abacus_neb import AbacusNEB

# setting
mpi = 16
omp = 1
fmax = 0.05  # eV / Ang
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
climb = True
dyneb = True
parallel = False
k = 0.10 # eV/Ang^2, spring constant
init_chain = "init_neb_chain.traj"
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}"
basis_dir = f"{lib_dir}"
pp = {
      'Li':'C_ONCV_PBE-1.2.upf',
      'Si':'H_ONCV_PBE-1.2.upf',
      }
basis = {
         'Li': 'Li_gga_8au_100Ry_4s1p.orb',
         'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
         ,}
kpts = [2, 2, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 1,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': -1,
}


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)

```

    Writing dyneb_run.py
    

在终端内通过vim等方法，按如上方式编辑完成之后，即可调用python运行该程序，以实现ASE调用ABACUS进行过渡态计算，该过程耗时约2h，感兴趣的用户可以自行计算。去除下面的#号注释即可运行


```
# ！python3 dyneb_run.py
```

我们直接来看计算输出，ATST-Tools的计算输出实际就是ASE在计算NEB过程中的输出，只是通过ATST-Tools，我们方便地调用了ABACUS作为ASE的计算单元


```
cat ../dyneb/running_neb.out
```

    Notice: Parallel calculation is not being used
    Set NEB method to DyNEB automatically
    ----- Running Dynamic NEB -----
    ----- improvedtangent method is being used -----
    ----- Default scale_fmax = 1.0 -----
          Step     Time          Energy          fmax
    FIRE:    0 18:28:24    -7051.845345         0.733359
    FIRE:    1 18:34:09    -7051.879062         0.594320
    FIRE:    2 18:39:55    -7051.922971         0.350617
    FIRE:    3 18:45:23    -7051.949616         0.195919
    FIRE:    4 18:50:57    -7051.951295         0.224455
    FIRE:    5 18:56:39    -7051.952800         0.217455
    FIRE:    6 19:04:29    -7051.955533         0.203700
    FIRE:    7 19:12:16    -7051.958973         0.183747
    FIRE:    8 19:19:13    -7051.962508         0.159008
    FIRE:    9 19:24:39    -7051.965466         0.134910
    FIRE:   10 19:30:23    -7051.967701         0.103314
    FIRE:   11 19:35:58    -7051.968988         0.112794
    FIRE:   12 19:41:36    -7051.969789         0.129614
    FIRE:   13 19:47:16    -7051.970598         0.138944
    FIRE:   14 19:52:55    -7051.972043         0.136817
    FIRE:   15 19:58:38    -7051.974469         0.120972
    FIRE:   16 20:04:24    -7051.977592         0.091176
    FIRE:   17 20:09:53    -7051.980373         0.075211
    FIRE:   18 20:11:43    -7051.981981         0.052248
    FIRE:   19 20:13:33    -7051.982934         0.087704
    FIRE:   20 20:17:14    -7051.984384         0.084422
    FIRE:   21 20:22:44    -7051.985587         0.112840
    FIRE:   22 20:26:29    -7051.985587         0.086535
    FIRE:   23 20:28:27    -7051.985587         0.036730
    
    === n_max set to 0, automatically detect the images of chain by NEBTools ===
    Number of images per band guessed to be 5.
    num: 0; Energy: -7052.6037618 (eV)
    num: 1; Energy: -7052.2767038 (eV)
    num: 2; Energy: -7051.9855868 (eV)
    num: 3; Energy: -7052.2797724 (eV)
    num: 4; Energy: -7052.6035978 (eV)
    Reaction Barrier and Energy Difference: (0.61817499999961, 0.00016400000004068715) (eV)
    Number of images per band guessed to be 5.
    Processing band         23 /         24
    Appears to be only one band in the images.
    Processing band          0 /          1

直接运行`python3 neb_run.py`只会输出运行过程中的一系列时间-能量和力的结果，后续的NEB后处理通过调用`ATST-Tools/neb`目录下的`neb_post.py`完成


```
# ! python3 /opt/ATST-Tools/neb/neb_post.py neb.traj
```

`neb.traj`，存储了NEB计算过程中的所有计算轨迹。可以通过`ase -T gui neb.traj`可视化。它经过`neb_post.py`后处理之后，输出结果除绝对能量信息外，会包括：
- `nebplots_chain.pdf`，绘制了计算收敛时的NEB链能量-反应坐标曲线，并且该曲线上有各个映像的受力在反应坐标方向的投影，所做的能量-反应坐标曲线也是基于这个投影结果进行了拟合得到的。
- `neb_latest.traj`，存储了NEB计算收敛时的NEB轨迹。可以通过`ase -T gui neb_latest.traj`可视化

<div>
    <img src="https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/938e9557a9654e9981ba1d92b678a37d/5c3134c1-df0d-446a-b778-b3f7fa43cf60.png", alt="reaction_diffusion" width="450" title="reaction_diffusion">
    <p style='font-size:0.8rem; font-weight:bold'> Li-diffu-Si 案例的计算结果 </p>
</div>

直接终端可视化在Bohrium里面不能完成，可以采用`ase.visualizer`模块的`view`函数，注意此时需要调用ngl可视化器，如果没有的话需要`pip install nglview`安装


```
%pip install nglview
```


```
from ase.io import read
from ase.visualize import view

atoms = read("/opt/ATST-Tools/examples/Li-diffu-Si/dyneb/neb_latest.traj", ":")
view(atoms, viewer='ngl')
```




    HBox(children=(NGLWidget(max_frame=4), VBox(children=(Dropdown(description='Show', options=('All', 'Si', 'Li')…



并行运行需要对脚本进行一些修改，开启其并行特征。但最关键的还是在调用脚本的时候启动并行，并调整并行核数与NEB链的映像个数一致。


```
cd ..
```

    /opt/ATST-Tools/examples/Li-diffu-Si
    


```
mkdir neb_para_example
```


```
cd neb_para_example
```

    /opt/ATST-Tools/examples/Li-diffu-Si/neb_para_example
    

我们此前使用的`init_neb_chain.traj`是可以直接复用的。你没有想错，这也是利用ATST-Tools在计算NEB时做续算的关键。


```
cp ../dyneb-example/init_neb_chain.traj .
```

简单修改一下`neb_run.py`，改改核数与并行NEB相关参数即可


```
%%writefile neb_run.py

# parallel NEB cannot be directly run in notebook

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from abacus_neb import AbacusNEB

# setting
mpi = 5
omp = 1
fmax = 0.05  # eV / Ang
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "NEBrun"
algorism = "improvedtangent" # IT-NEB is recommended
climb = True
dyneb = False
parallel = True
k = 0.10 # eV/Ang^2, spring constant
init_chain = "init_neb_chain.traj"
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}"
basis_dir = f"{lib_dir}"
pp = {
      'Li':'C_ONCV_PBE-1.2.upf',
      'Si':'H_ONCV_PBE-1.2.upf',
      }
basis = {
         'Li': 'Li_gga_8au_100Ry_4s1p.orb',
         'Si': 'Si_gga_8au_100Ry_2s2p1d.orb',
         ,}
kpts = [2, 2, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 1,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'none',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
}


if __name__ == "__main__":
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory=neb_directory, mpi=mpi, omp=omp, abacus=abacus, 
                    algorism=algorism, k=k, dyneb=dyneb)
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
```

    Writing neb_run.py
    

在终端准备好脚本之后，即可通过`mpirun gpaw python`进行并行NEB计算。当然了，其中的ABACUS也是用mpirun调用的，这两个mpirun需要对应同一个mpirun程序，你的镜像应该帮你准备好了。如果你有兴趣的话，可以去掉下面的注释，直接跑跑——当然这需要至少16核的计算环境


```
# ! mpirun -np 3 gpaw python neb_run.py
```

我们直接来看计算结果


```
cat ../neb/running_neb.out
```

    ----- Running ASE-NEB -----
    ----- improvedtangent method is being used -----
    ----- Parallel calculation is being used -----
          Step     Time          Energy          fmax
    FIRE:    0 16:50:48    -7051.845345         0.733359
    FIRE:    1 16:53:44    -7051.879062         0.594320
    FIRE:    2 16:56:39    -7051.922971         0.350617
    FIRE:    3 16:59:18    -7051.949616         0.195919
    FIRE:    4 17:02:01    -7051.951295         0.224455
    FIRE:    5 17:04:47    -7051.952800         0.217455
    FIRE:    6 17:07:33    -7051.955533         0.203700
    FIRE:    7 17:10:25    -7051.958973         0.183747
    FIRE:    8 17:13:04    -7051.962508         0.159008
    FIRE:    9 17:15:42    -7051.965466         0.134910
    FIRE:   10 17:18:34    -7051.967701         0.103314
    FIRE:   11 17:21:25    -7051.968988         0.112794
    FIRE:   12 17:24:13    -7051.969789         0.129614
    FIRE:   13 17:27:03    -7051.970598         0.138944
    FIRE:   14 17:29:55    -7051.972043         0.136817
    FIRE:   15 17:32:49    -7051.974469         0.120972
    FIRE:   16 17:35:43    -7051.977592         0.091176
    FIRE:   17 17:38:31    -7051.980373         0.075211
    FIRE:   18 17:41:20    -7051.981840         0.058992
    FIRE:   19 17:44:06    -7051.982561         0.090015
    FIRE:   20 17:46:57    -7051.983863         0.090542
    FIRE:   21 17:49:49    -7051.985435         0.047680
    
    === n_max set to 0, automatically detect the images of chain by NEBTools ===
    Number of images per band guessed to be 5.
    num: 0; Energy: -7052.6037618 (eV)
    num: 1; Energy: -7052.277576 (eV)
    num: 2; Energy: -7051.985435 (eV)
    num: 3; Energy: -7052.2790721 (eV)
    num: 4; Energy: -7052.6035978 (eV)
    Reaction Barrier and Energy Difference: (0.6183268000004314, 0.00016400000004068715) (eV)
    Number of images per band guessed to be 5.
    Processing band         21 /         22
    Appears to be only one band in the images.
    Processing band          0 /          1

上述过程基本上就是调用ATST-Tools进行NEB计算的全过程，这一过程在Slurm队列管理的服务器上可以更容易地运行，用户只需要使用`ATST-Tools/neb`下的`neb_submit.sh`脚本即可。同理，ATST-Tools中还提供了`autoneb_submit.sh`、`dimer_submit.sh`等适用于Slurm队列的shell脚本，帮助用户在服务器上快速运行ASE-ABACUS过渡态计算。

#### 2.2.3 案例二：单原子催化的AutoNEB

我们来看一个更贴近催化计算实际的例子，`Cy-Pt@graphene`是环己烷在单原子Pt负载石墨烯上的解离

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/653e1e52311e491e87fb27dadbc674f7/e8a8e628-3a98-45ed-aa68-cef0806faedd.png)

进入到该文件夹中：


```
cd /opt/ATST-Tools/examples/Cy-Pt@graphene
```

    /opt/ATST-Tools/examples/Cy-Pt@graphene
    


```
ls
```

    [0m[01;34mautoneb[0m/  [01;34mdata[0m/  [01;34meb-neb[0m/  [01;34mit-neb[0m/  [01;34mneb2sella[0m/  [01;34msella[0m/
    

该案例的NEB采用8个image描述反应过程，并且由于环己烷分子的非极性性、负载石墨烯的柔性等因素，这一NEB计算是存在挑战的。如果直接去跑CI-NEB计算，一方面需要过大的计算量，另一方面也容易因为对过渡态的解析度不够，搜索到错误的反应路径和过渡态。本例中的`it-neb`采用的即是最常用的IT-NEB+CI-NEB的方法，但它所用NEB步数非常长，且收敛到了错误的结果

<div>
    <img src="https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/24c55852f0b146c4a38fe03bfd4b255c/dc6c9dd4-eb1d-4439-91a4-79e48635db1e.png", alt="reaction_diffusion" width="450" title="reaction_diffusion">
    <p style='font-size:0.8rem; font-weight:bold'> Cy-Pt@graphene 案例的错误计算结果 </p>
</div>

至于如何判断过渡态本身是正确的还是错误的，一个简单的方法是直接看图。ASE的NEB后处理模块`NEBTools`中的作图模块是会在图中以切线的形式，画出结构在NEB链切线方向的受力的投影的，这个力又称切线力，它代表了该点处受力在反应坐标下的投影。

在能量最高点，如果它是一个正确的过渡态，这个投影切线会是基本水平的。当然，最终还是要通过振动计算得到反应坐标方向的单一虚频来进行验证。这个振动计算可以通过phonopy或者ASE完成，在ATST-Tools当前版本（1.5.0）中，以独立脚本形式集成了ASE的振动计算和热力学较正方法，后续会有介绍。




于是，我们可以采用AutoNEB方法，一方面减少计算量，另一方面提升过渡态附近的解析度，从而在更高的计算效率下计算出正确的过渡态。

进行这一AutoNEB计算可以分为如下几个步骤
1. 准备所需的初末态计算结果，以及计算所需的赝势和轨道文件
2. 通过`neb_make.py`生成初猜结构
2. 通过`autoneb_run.py`设置计算参数并开展计算

让我们开始吧


```
mkdir autoneb_test
```


```
cd autoneb_test
```

    /opt/ATST-Tools/examples/Cy-Pt@graphene/autoneb_test
    


```
cp -r ../data/* .
```


```
ls
```

    C_ONCV_PBE-1.0.upf          H_gga_6au_100Ry_2s1p.orb
    C_gga_7au_100Ry_2s2p1d.orb  [0m[01;34mIS[0m/
    [01;34mFS[0m/                         Pt_ONCV_PBE-1.0.upf
    H_ONCV_PBE-1.0.upf          Pt_gga_7au_100Ry_4s2p2d1f.orb
    


```
# 生成初猜，采用4个image作为AutoNEB的初始image数与并行进程数
! python3 /opt/ATST-Tools/neb/neb_make.py -i IS/OUT.ABACUS/running_relax.log FS/OUT.ABACUS/running_relax.log -n 4
```

    Reading files: IS/OUT.ABACUS/running_relax.log and FS/OUT.ABACUS/running_relax.log
    Generating path, number of images: 4, sort_tol: 1.0
    Optimizing path using IDPP method
    /opt/mamba/lib/python3.10/site-packages/pymatgen/analysis/diffusion/neb/pathfinder.py:222: UserWarning: Auto sorting is turned off because it is unable to match the end-point structures!
      warnings.warn(
    Writing path: init_neb_chain.traj,Number of images: 6
    


```
# 复制已有的autoneb_run.py并查看
! cp /opt/ATST-Tools/neb/autoneb_run.py .
```


```
cat autoneb_run.py
```

    # JamesMisaka in 2023-11-27
    # Run AutoNEB calculation by ASE-ABACUS
    # part of ATST-Tools scripts
    
    
    from ase.optimize import FIRE, BFGS
    from ase.io import read, write
    from ase.parallel import world, parprint, paropen
    #from pathlib import Path
    from abacus_autoneb import AbacusAutoNEB
    
    # setting
    mpi = 16
    omp = 4
    neb_optimizer = FIRE # suited for CI-NEB
    neb_directory = "AutoNEBrun"
    # algorism = 'eb' # default for AutoNEB
    algorism = "improvedtangent" # IT-NEB
    init_chain = "init_neb_chain.traj"
    climb = True
    fmax = [0.20, 0.05]  # eV / Ang, 2 fmax for others and last CI-NEB in all-images
    n_simul = world.size # only for autoneb, number of simultaneous calculation
    n_images = 10 # only for autoneb, max number of all image, which should be reached
    smooth_curve = False # True to do more neb step to smooth the curve.
    k = 0.10 # eV/Ang^2, force constant of spring, 0.05 is from VTST-Tools
    abacus = "abacus"
    #lib_dir = "/lustre/home/2201110432/example/abacus"
    lib_dir = "/lustre/home/2201110432/example/abacus"
    pseudo_dir = f"{lib_dir}/PP"
    basis_dir = f"{lib_dir}/ORB"
    # default pp and basis is supported by ase-abacus interface, need to check usage
    pp = {
            "H":"H_ONCV_PBE-1.0.upf",
            "C":"C_ONCV_PBE-1.0.upf",
            "O":"O_ONCV_PBE-1.0.upf",
            "Fe":"Fe_ONCV_PBE-1.0.upf",
          }
    basis = {
            "H":"H_gga_6au_100Ry_2s1p.orb",
            "C":"C_gga_7au_100Ry_2s2p1d.orb",
            "O":"O_gga_7au_100Ry_2s2p1d.orb",
            "Fe":"Fe_gga_8au_100Ry_4s2p2d1f.orb",
            }
    kpts = [3, 1, 2]
    parameters = {
        'calculation': 'scf',
        'nspin': 2,
        'xc': 'pbe',
        'ecutwfc': 100,
        'ks_solver': 'genelpa',
        'symmetry': 0,
        'vdw_method': 'd3_bj',
        'smearing_method': 'mp',
        'smearing_sigma': 0.002,
        'basis_type': 'lcao',
        'mixing_type': 'broyden',
        'mixing_beta': 0.4,
        'mixing_gg0': 1.0,
        'mixing_ndim': 8,
        'scf_thr': 1e-7,
        'scf_nmax': 300,
        'kpts': kpts,
        'pp': pp,
        'basis': basis,
        'pseudo_dir': pseudo_dir,
        'basis_dir': basis_dir,
        'cal_force': 1,
        'cal_stress': 0,
        'init_wfc': 'atomic',
        'init_chg': 'atomic',
        'out_stru': 1,
        'out_chg': -1,
        'out_mul': 1,
        'out_wfc_lcao': 0,
        'out_bandgap': 1,
        'efield_flag': 1,
        'dip_cor_flag': 1,
        'efield_dir': 1,
    }
    
    
    
    if __name__ == "__main__": 
    # running process
    # read initial guessed neb chain
        init_chain = read(init_chain, index=':')
        neb = AbacusAutoNEB(init_chain, parameters, algorism=algorism, 
                            directory=neb_directory, k=k,
                            n_simul=n_simul, n_max=n_images, 
                            abacus=abacus,  mpi=mpi, omp=omp, )
        neb.run(optimizer=neb_optimizer, climb=climb, 
                    fmax=fmax, smooth_curve=smooth_curve)

同样的，需要对该脚本内的赝势、轨道和ABACUS计算设置进行一些自定义


```
%%writefile autoneb_run.py

# parallel AutoNEB cannot be directly run in notebook

from ase.optimize import FIRE, BFGS
from ase.io import read, write
from ase.parallel import world, parprint, paropen
#from pathlib import Path
from abacus_autoneb import AbacusAutoNEB

# setting
mpi = 16
omp = 4
neb_optimizer = FIRE # suited for CI-NEB
neb_directory = "AutoNEBrun"
# algorism = 'eb' # default for AutoNEB
algorism = "improvedtangent" # IT-NEB
init_chain = "init_neb_chain.traj"
climb = True
fmax = [0.20, 0.05]  # eV / Ang, 2 fmax for others and last CI-NEB in all-images
n_simul = world.size # only for autoneb, number of simultaneous calculation
n_images = 10 # only for autoneb, max number of all image, which should be reached
smooth_curve = False # True to do more neb step to smooth the curve.
k = 0.10 # eV/Ang^2, force constant of spring, 0.05 is from VTST-Tools
abacus = "abacus"
#lib_dir = "/lustre/home/2201110432/example/abacus"
lib_dir = ""
pseudo_dir = f"{lib_dir}/"
basis_dir = f"{lib_dir}/"
# default pp and basis is supported by ase-abacus interface, need to check usage
pp = {
      'C':'C_ONCV_PBE-1.0.upf',
      'H':'H_ONCV_PBE-1.0.upf',
      'Pt':'Pt_ONCV_PBE-1.0.upf',
      }
basis = {
         'C': 'C_gga_7au_100Ry_2s2p1d.orb',
         'H': 'H_gga_6au_100Ry_2s1p.orb',
         'Pt': 'Pt_gga_7au_100Ry_4s2p2d1f.orb',
         }
kpts = [2, 1, 2]
parameters = {
    'calculation': 'scf',
    'nspin': 2,
    'xc': 'pbe',
    'ecutwfc': 100,
    'dft_functional': 'pbe',
    'ks_solver': 'genelpa',
    'symmetry': 0,
    'vdw_method': 'd3_bj',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'basis_type': 'lcao',
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    'scf_nmax': 100,
    'kpts': kpts,
    'pp': pp,
    'basis': basis,
    'pseudo_dir': pseudo_dir,
    'basis_dir': basis_dir,
    'init_wfc': 'atomic',
    'init_chg': 'atomic',
    'cal_force': 1,
    'cal_stress': 1,
    'out_stru': 1,
    'out_chg': 0,
    'out_mul': 0,
    'out_wfc_lcao': 0,
    'out_bandgap': 0,
    'efield_flag': 1,
    'dip_cor_flag': 1,
    'efield_dir': 1,
}



if __name__ == "__main__": 
# running process
# read initial guessed neb chain
    init_chain = read(init_chain, index=':')
    neb = AbacusAutoNEB(init_chain, parameters, algorism=algorism, 
                        directory=neb_directory, k=k,
                        n_simul=n_simul, n_max=n_images, 
                        abacus=abacus,  mpi=mpi, omp=omp, )
    neb.run(optimizer=neb_optimizer, climb=climb, 
                fmax=fmax, smooth_curve=smooth_curve)
```

    Overwriting autoneb_run.py
    

这一脚本的内容已经基本满足了对该实例进行AutoNEB的需要，用户可以通过上面提到的方法，简单微调以下改脚本的相关设置即可进行AutoNEB计算。运行方式也是一样的

```
mpirun -np 4 gpaw python autoneb_run.py
```
我们直接来看结果。


```
cd /opt/ATST-Tools/examples/Cy-Pt@graphene/autoneb
```

    /opt/ATST-Tools/examples/Cy-Pt@graphene/autoneb
    


```
# 显示运行结果文件
! ls
```

    AutoNEB_iter	     run_autoneb000.traj  run_autoneb006.traj
    autoneb_run.py	     run_autoneb001.traj  run_autoneb007.traj
    autoneb_submit.sh    run_autoneb002.traj  run_autoneb008.traj
    init_neb_chain.traj  run_autoneb003.traj  run_autoneb009.traj
    neb_latest.traj      run_autoneb004.traj  running_autoneb.out
    nebplots_all.pdf     run_autoneb005.traj  vib_analysis_TS
    


```
# 查看整体运算结果，注意该运算结果实际是用autoneb_submit.sh得到的
! cat running_autoneb.out
```

    Loading mkl version 2023.0.0
    Loading tbb version 2021.8.0
    Loading compiler-rt version 2023.0.0
    Loading mpi version 2021.8.0
    Loading compiler version 2023.0.0
    Loading oclfpga version 2023.0.0
      Load "debugger" to debug DPC++ applications with the gdb-oneapi debugger.
      Load "dpl" for additional DPC++ APIs: https://github.com/oneapi-src/oneDPL
    
    Loading abacus/3.4.2-icx
      Loading requirement: tbb/latest compiler-rt/latest mkl/2023.0.0 mpi/latest
        oclfpga/latest compiler/latest
    NSIMUL is 4
    ===== AutoNEB Job Starting =====
     ===== Make Initial NEB Guess =====
    ---- Fix Atoms below 0.0 in direction y ----
    ---- Set initial magmom for [] to [] ----
    ---- Fix Atoms below 0.0 in direction y ----
    ---- Set initial magmom for [] to [] ----
    ---- Fix Atoms below 0.0 in direction y ----
    ---- Set initial magmom for [] to [] ----
    ---- Fix Atoms below 0.0 in direction y ----
    ---- Set initial magmom for [] to [] ----
    --- Successfully make guessed image chain by idpp method ! ---
    ===== Running AutoNEB =====
    Notice: AutoNEB method is set
    You manually set n_simul = 4, n_max = 10
    ----- Running AutoNEB -----
    ----- improvedtangent method is being used -----
    The NEB initially has 6 images  (including the end-points)
    Start of evaluation of the initial images
    Now starting iteration 1 on  [0, 1, 2, 3, 4, 5]
    Finished initialisation phase.
    ****Now adding another image until n_max is reached (6/10)****
    Adding image between 2 and 3. New image point is selected on the basis of the biggest spring length!
    Now starting iteration 2 on  [1, 2, 3, 4, 5, 6]
    ****Now adding another image until n_max is reached (7/10)****
    Adding image between 1 and 2. New image point is selected on the basis of the biggest spring length!
    Now starting iteration 3 on  [1, 2, 3, 4, 5, 6]
    ****Now adding another image until n_max is reached (8/10)****
    Adding image between 5 and 6. New image point is selected on the basis of the biggest spring length!
    Now starting iteration 4 on  [2, 3, 4, 5, 6, 7]
    ****Now adding another image until n_max is reached (9/10)****
    Adding image between 7 and 8. New image point is selected on the basis of the biggest spring length!
    Now starting iteration 5 on  [4, 5, 6, 7, 8, 9]
    n_max images has been reached
    ****Now doing the CI-NEB calculation****
    Now starting iteration 6 on  [2, 3, 4, 5, 6, 7]
    ----- AutoNEB calculation finished -----
    ===== AutoNEB Process Done ! =====
    ===== Running Post-Processing =====
    === n_max set to 0, automatically detect the images of chain by NEBTools ===
    Appears to be only one band in the images.
    num: 0; Energy: -11866.824886 (eV)
    num: 1; Energy: -11866.82722 (eV)
    num: 2; Energy: -11866.801706 (eV)
    num: 3; Energy: -11866.78262 (eV)
    num: 4; Energy: -11866.587954 (eV)
    num: 5; Energy: -11865.497 (eV)
    num: 6; Energy: -11866.330196 (eV)
    num: 7; Energy: -11866.396832 (eV)
    num: 8; Energy: -11866.429487 (eV)
    num: 9; Energy: -11866.435419 (eV)
    Reaction Barrier and Energy Difference: (1.3278860000009587, 0.3894670000008773) (eV)
    Appears to be only one band in the images.
    Processing band          0 /          1
    ===== Done ! =====
    


```
# 查看AutoNEB中第一步NEB的计算结果
! cat AutoNEB_iter/run_autoneb_log_iter001.log
```

          Step     Time          Energy          fmax
    FIRE:    0 00:35:21   -11864.578532         3.642360
    FIRE:    1 00:36:10   -11864.782954         2.870872
    FIRE:    2 00:36:56   -11865.029191         2.407401
    FIRE:    3 00:37:40   -11865.183575         3.420662
    FIRE:    4 00:38:26   -11865.149576         4.842759
    FIRE:    5 00:39:11   -11865.185455         3.786428
    FIRE:    6 00:40:01   -11865.238631         2.154321
    FIRE:    7 00:40:47   -11865.284643         1.979722
    FIRE:    8 00:41:33   -11865.307675         1.889517
    FIRE:    9 00:42:20   -11865.309583         1.742232
    FIRE:   10 00:43:06   -11865.304529         1.941274
    FIRE:   11 00:43:51   -11865.305786         1.824236
    FIRE:   12 00:44:39   -11865.321340         1.566756
    FIRE:   13 00:45:26   -11865.351452         1.232224
    FIRE:   14 00:46:12   -11865.385412         1.163525
    FIRE:   15 00:47:00   -11865.405953         1.863983
    FIRE:   16 00:47:46   -11865.411933         1.196430
    FIRE:   17 00:48:34   -11865.426725         1.069129
    FIRE:   18 00:49:19   -11865.431487         0.961849
    FIRE:   19 00:50:06   -11865.439472         0.760762
    FIRE:   20 00:50:51   -11865.448191         0.517175
    FIRE:   21 00:51:39   -11865.455199         0.459462
    FIRE:   22 00:52:25   -11865.459289         0.472398
    FIRE:   23 00:53:13   -11865.460873         0.431848
    FIRE:   24 00:54:00   -11865.461725         0.539257
    FIRE:   25 00:54:45   -11865.463809         0.619627
    FIRE:   26 00:55:31   -11865.468630         0.598042
    FIRE:   27 00:56:20   -11865.476511         0.468127
    FIRE:   28 00:57:07   -11865.485968         0.365569
    FIRE:   29 00:57:53   -11865.494549         0.412251
    FIRE:   30 00:58:43   -11865.501988         0.385664
    FIRE:   31 00:59:30   -11865.511484         0.422284
    FIRE:   32 01:00:16   -11865.525009         0.366845
    FIRE:   33 01:01:02   -11865.539757         0.231356
    FIRE:   34 01:01:50   -11865.554350         0.309826
    FIRE:   35 01:02:39   -11865.570950         0.294388
    FIRE:   36 01:03:24   -11865.590207         0.275131
    FIRE:   37 01:04:08   -11865.611798         0.247073
    FIRE:   38 01:04:53   -11865.636943         0.299314
    FIRE:   39 01:05:39   -11865.672444         0.295211
    FIRE:   40 01:06:26   -11865.724616         0.398352
    FIRE:   41 01:07:17   -11865.809261         0.672083
    FIRE:   42 01:08:02   -11865.836765         0.677076
    FIRE:   43 01:08:49   -11865.840971         0.681209
    FIRE:   44 01:09:39   -11865.848278         0.688700
    FIRE:   45 01:10:27   -11865.858156         0.697759
    FIRE:   46 01:11:13   -11865.871159         0.704303
    FIRE:   47 01:11:59   -11865.887878         0.703259
    FIRE:   48 01:12:47   -11865.907860         0.688451
    FIRE:   49 01:13:34   -11865.929898         0.651860
    FIRE:   50 01:14:24   -11865.955356         0.574260
    FIRE:   51 01:15:10   -11865.982482         0.436765
    FIRE:   52 01:15:55   -11866.006105         0.314752
    FIRE:   53 01:16:44   -11866.020034         0.413740
    FIRE:   54 01:17:31   -11866.025234         0.697364
    FIRE:   55 01:18:18   -11866.027275         0.653814
    FIRE:   56 01:19:01   -11866.030891         0.572655
    FIRE:   57 01:19:46   -11866.035293         0.464302
    FIRE:   58 01:20:32   -11866.039733         0.341400
    FIRE:   59 01:21:23   -11866.043688         0.219902
    FIRE:   60 01:22:07   -11866.046996         0.200859
    FIRE:   61 01:22:55   -11866.049775         0.236517
    FIRE:   62 01:23:41   -11866.052505         0.255194
    FIRE:   63 01:24:30   -11866.055451         0.245242
    FIRE:   64 01:25:18   -11866.058968         0.220082
    FIRE:   65 01:26:06   -11866.063303         0.193788
    


```
# 查看AutoNEB中最后一步CI-NEB的计算结果
! cat AutoNEB_iter/run_autoneb_log_iter006.log
```

          Step     Time          Energy          fmax
    FIRE:    0 01:57:54   -11866.084572         0.652487
    FIRE:    1 01:58:40   -11866.076705         0.676456
    FIRE:    2 01:59:28   -11866.059007         0.715741
    FIRE:    3 02:00:13   -11866.028561         0.745702
    FIRE:    4 02:00:59   -11865.982395         0.745921
    FIRE:    5 02:01:46   -11865.918309         0.713379
    FIRE:    6 02:02:36   -11865.835733         0.659343
    FIRE:    7 02:03:21   -11865.736940         0.602895
    FIRE:    8 02:04:09   -11865.616163         0.570285
    FIRE:    9 02:04:56   -11865.486295         0.576872
    FIRE:   10 02:05:41   -11865.366292         0.553076
    FIRE:   11 02:06:26   -11865.290706         0.491710
    FIRE:   12 02:07:10   -11865.268212         0.382165
    FIRE:   13 02:07:56   -11865.295327         0.349660
    FIRE:   14 02:08:42   -11865.381779         0.505350
    FIRE:   15 02:09:29   -11865.389766         0.478347
    FIRE:   16 02:10:14   -11865.402261         0.431838
    FIRE:   17 02:11:01   -11865.414176         0.381902
    FIRE:   18 02:11:46   -11865.421392         0.382459
    FIRE:   19 02:12:33   -11865.422792         0.386433
    FIRE:   20 02:13:18   -11865.419163         0.378440
    FIRE:   21 02:14:05   -11865.411475         0.372510
    FIRE:   22 02:14:52   -11865.399699         0.363692
    FIRE:   23 02:15:36   -11865.382875         0.311966
    FIRE:   24 02:16:22   -11865.360650         0.348324
    FIRE:   25 02:17:08   -11865.335546         0.494877
    FIRE:   26 02:17:55   -11865.316262         0.735960
    FIRE:   27 02:18:45   -11865.307194         1.100647
    FIRE:   28 02:19:30   -11865.309266         0.196186
    FIRE:   29 02:20:19   -11865.313339         0.909086
    FIRE:   30 02:21:08   -11865.313751         0.707882
    FIRE:   31 02:21:54   -11865.314565         0.350729
    FIRE:   32 02:22:42   -11865.315812         0.178845
    FIRE:   33 02:23:30   -11865.315906         0.178576
    FIRE:   34 02:24:16   -11865.316091         0.178013
    FIRE:   35 02:25:07   -11865.316373         0.177024
    FIRE:   36 02:25:57   -11865.316741         0.175838
    FIRE:   37 02:26:45   -11865.317208         0.174313
    FIRE:   38 02:27:32   -11865.317762         0.172504
    FIRE:   39 02:28:19   -11865.318410         0.170232
    FIRE:   40 02:29:08   -11865.319231         0.167451
    FIRE:   41 02:29:54   -11865.320258         0.163917
    FIRE:   42 02:30:40   -11865.321550         0.165854
    FIRE:   43 02:31:28   -11865.323145         0.171715
    FIRE:   44 02:32:18   -11865.325136         0.179361
    FIRE:   45 02:33:07   -11865.327600         0.188798
    FIRE:   46 02:33:54   -11865.330675         0.200416
    FIRE:   47 02:34:43   -11865.334504         0.213338
    FIRE:   48 02:35:32   -11865.339264         0.226164
    FIRE:   49 02:36:21   -11865.345094         0.234130
    FIRE:   50 02:37:09   -11865.352065         0.232614
    FIRE:   51 02:37:54   -11865.360183         0.216983
    FIRE:   52 02:38:43   -11865.369363         0.185817
    FIRE:   53 02:39:29   -11865.379371         0.148843
    FIRE:   54 02:40:14   -11865.389614         0.138991
    FIRE:   55 02:40:59   -11865.399380         0.137180
    FIRE:   56 02:41:45   -11865.408222         0.165171
    FIRE:   57 02:42:32   -11865.416086         0.203523
    FIRE:   58 02:43:19   -11865.423394         0.231605
    FIRE:   59 02:44:06   -11865.430772         0.247774
    FIRE:   60 02:44:53   -11865.438325         0.254307
    FIRE:   61 02:45:41   -11865.445080         0.267498
    FIRE:   62 02:46:29   -11865.448719         0.881362
    FIRE:   63 02:47:15   -11865.451780         0.231961
    FIRE:   64 02:48:02   -11865.452007         0.231221
    FIRE:   65 02:48:48   -11865.452438         0.229862
    FIRE:   66 02:49:34   -11865.453071         0.227879
    FIRE:   67 02:50:21   -11865.453876         0.225257
    FIRE:   68 02:51:08   -11865.454816         0.222156
    FIRE:   69 02:52:00   -11865.455865         0.218612
    FIRE:   70 02:52:44   -11865.457001         0.214634
    FIRE:   71 02:53:30   -11865.458341         0.209595
    FIRE:   72 02:54:15   -11865.459892         0.203370
    FIRE:   73 02:55:01   -11865.461651         0.195673
    FIRE:   74 02:55:49   -11865.463637         0.186376
    FIRE:   75 02:56:36   -11865.465861         0.175254
    FIRE:   76 02:57:22   -11865.468343         0.162491
    FIRE:   77 02:58:09   -11865.471156         0.147879
    FIRE:   78 02:58:58   -11865.474326         0.131731
    FIRE:   79 02:59:44   -11865.477956         0.114483
    FIRE:   80 03:00:31   -11865.482135         0.096221
    FIRE:   81 03:01:17   -11865.486787         0.084464
    FIRE:   82 03:02:05   -11865.491684         0.134461
    FIRE:   83 03:02:50   -11865.496380         0.273588
    FIRE:   84 03:03:36   -11865.496522         0.051535
    FIRE:   85 03:04:21   -11865.496806         0.203176
    FIRE:   86 03:05:08   -11865.496838         0.161525
    FIRE:   87 03:05:53   -11865.496901         0.086825
    FIRE:   88 03:06:42   -11865.497000         0.040039
    

这样计算出来的就是正确的过渡态，我们可以在命令行中使用`ase -T gui`来进行可视化: （前提是Bohrium上装有对应的可视化模块，至少在命令行里面这似乎是办不到的）


```
! ase -T gui ./run_autoneb005.traj
```

    Traceback (most recent call last):
      File "/opt/mamba/bin/ase", line 8, in <module>
        sys.exit(main())
      File "/opt/mamba/lib/python3.10/site-packages/ase/cli/main.py", line 102, in main
        f(args)
      File "/opt/mamba/lib/python3.10/site-packages/ase/gui/ag.py", line 106, in run
        gui = GUI(images, args.rotations, args.bonds, args.graph)
      File "/opt/mamba/lib/python3.10/site-packages/ase/gui/gui.py", line 45, in __init__
        self.window = ui.ASEGUIWindow(close=self.exit, menu=menu,
      File "/opt/mamba/lib/python3.10/site-packages/ase/gui/ui.py", line 597, in __init__
        MainWindow.__init__(self, 'ASE-GUI', close, menu)
      File "/opt/mamba/lib/python3.10/site-packages/ase/gui/ui.py", line 498, in __init__
        self.win = tk.Tk()
      File "/opt/mamba/lib/python3.10/tkinter/__init__.py", line 2299, in __init__
        self.tk = _tkinter.create(screenName, baseName, className, interactive, wantobjects, useTk, sync, use)
    _tkinter.TclError: no display name and no $DISPLAY environment variable
    

命令行可视化不可用的情况下可以通过python代码完成。在Bohrium上，可以采用ngl viewer来可视化ASE结果


```
%pip install nglview
```

    Looking in indexes: https://pypi.tuna.tsinghua.edu.cn/simple
    Collecting nglview
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/db/6e/f26f3350b49c1d813106a497e3e9089cc2b62e845f67d95bf82ed27d0539/nglview-3.1.2.tar.gz (5.5 MB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m5.5/5.5 MB[0m [31m23.9 MB/s[0m eta [36m0:00:00[0m00:01[0m00:01[0m
    [?25h  Installing build dependencies ... [?25ldone
    [?25h  Getting requirements to build wheel ... [?25ldone
    [?25h  Preparing metadata (pyproject.toml) ... [?25ldone
    [?25hCollecting ipywidgets>=8
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/22/2d/9c0b76f2f9cc0ebede1b9371b6f317243028ed60b90705863d493bae622e/ipywidgets-8.1.5-py3-none-any.whl (139 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m139.8/139.8 kB[0m [31m38.6 MB/s[0m eta [36m0:00:00[0m
    [?25hCollecting jupyterlab-widgets
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/a9/93/858e87edc634d628e5d752ba944c2833133a28fa87bb093e6832ced36a3e/jupyterlab_widgets-3.0.13-py3-none-any.whl (214 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m214.4/214.4 kB[0m [31m35.3 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: numpy in /opt/mamba/lib/python3.10/site-packages (from nglview) (1.26.4)
    Collecting notebook>=7
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/46/77/53732fbf48196af9e51c2a61833471021c1d77d335d57b96ee3588c0c53d/notebook-7.2.2-py3-none-any.whl (5.0 MB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m5.0/5.0 MB[0m [31m53.5 MB/s[0m eta [36m0:00:00[0m00:01[0m00:01[0m
    [?25hRequirement already satisfied: ipython>=6.1.0 in /opt/mamba/lib/python3.10/site-packages (from ipywidgets>=8->nglview) (8.11.0)
    Collecting comm>=0.1.3
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/e6/75/49e5bfe642f71f272236b5b2d2691cf915a7283cc0ceda56357b61daa538/comm-0.2.2-py3-none-any.whl (7.2 kB)
    Requirement already satisfied: traitlets>=4.3.1 in /opt/mamba/lib/python3.10/site-packages (from ipywidgets>=8->nglview) (5.9.0)
    Collecting widgetsnbextension~=4.0.12
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/21/02/88b65cc394961a60c43c70517066b6b679738caf78506a5da7b88ffcb643/widgetsnbextension-4.0.13-py3-none-any.whl (2.3 MB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m2.3/2.3 MB[0m [31m69.1 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: notebook-shim<0.3,>=0.2 in /opt/mamba/lib/python3.10/site-packages (from notebook>=7->nglview) (0.2.2)
    Requirement already satisfied: tornado>=6.2.0 in /opt/mamba/lib/python3.10/site-packages (from notebook>=7->nglview) (6.2)
    Collecting jupyterlab<4.3,>=4.2.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/fd/3f/24a0f0ce60959cfd9756a3291cd3a5581e51cbd6f7b4aa121f5bba5320e3/jupyterlab-4.2.5-py3-none-any.whl (11.6 MB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m11.6/11.6 MB[0m [31m87.1 MB/s[0m eta [36m0:00:00[0m00:01[0m00:01[0m
    [?25hCollecting jupyter-server<3,>=2.4.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/57/e1/085edea6187a127ca8ea053eb01f4e1792d778b4d192c74d32eb6730fed6/jupyter_server-2.14.2-py3-none-any.whl (383 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m383.6/383.6 kB[0m [31m56.0 MB/s[0m eta [36m0:00:00[0m
    [?25hCollecting jupyterlab-server<3,>=2.27.1
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/54/09/2032e7d15c544a0e3cd831c51d77a8ca57f7555b2e1b2922142eddb02a84/jupyterlab_server-2.27.3-py3-none-any.whl (59 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m59.7/59.7 kB[0m [31m17.6 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: stack-data in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (0.6.2)
    Requirement already satisfied: pygments>=2.4.0 in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (2.14.0)
    Requirement already satisfied: pickleshare in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (0.7.5)
    Requirement already satisfied: jedi>=0.16 in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (0.18.2)
    Requirement already satisfied: matplotlib-inline in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (0.1.6)
    Requirement already satisfied: backcall in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (0.2.0)
    Requirement already satisfied: pexpect>4.3 in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (4.8.0)
    Requirement already satisfied: prompt-toolkit!=3.0.37,<3.1.0,>=3.0.30 in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (3.0.38)
    Requirement already satisfied: decorator in /opt/mamba/lib/python3.10/site-packages (from ipython>=6.1.0->ipywidgets>=8->nglview) (5.1.1)
    Requirement already satisfied: argon2-cffi>=21.1 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (21.3.0)
    Collecting jupyter-events>=0.9.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/a5/94/059180ea70a9a326e1815176b2370da56376da347a796f8c4f0b830208ef/jupyter_events-0.10.0-py3-none-any.whl (18 kB)
    Requirement already satisfied: anyio>=3.1.0 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (3.6.2)
    Requirement already satisfied: pyzmq>=24 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (25.0.0)
    Requirement already satisfied: jupyter-client>=7.4.4 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (8.0.3)
    Requirement already satisfied: packaging>=22.0 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (23.0)
    Requirement already satisfied: prometheus-client>=0.9 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.16.0)
    Requirement already satisfied: nbformat>=5.3.0 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (5.7.3)
    Requirement already satisfied: nbconvert>=6.4.4 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (7.2.9)
    Collecting overrides>=5.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/2c/ab/fc8290c6a4c722e5514d80f62b2dc4c4df1a68a41d1364e625c35990fcf3/overrides-7.7.0-py3-none-any.whl (17 kB)
    Requirement already satisfied: jinja2>=3.0.3 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (3.1.2)
    Collecting send2trash>=1.8.2
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/40/b0/4562db6223154aa4e22f939003cb92514c79f3d4dccca3444253fd17f902/Send2Trash-1.8.3-py3-none-any.whl (18 kB)
    Requirement already satisfied: terminado>=0.8.3 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.17.1)
    Requirement already satisfied: jupyter-server-terminals>=0.4.4 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.4.4)
    Collecting websocket-client>=1.7
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/5a/84/44687a29792a70e111c5c477230a72c4b957d88d16141199bf9acb7537a3/websocket_client-1.8.0-py3-none-any.whl (58 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m58.8/58.8 kB[0m [31m18.5 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: jupyter-core!=5.0.*,>=4.12 in /opt/mamba/lib/python3.10/site-packages (from jupyter-server<3,>=2.4.0->notebook>=7->nglview) (5.2.0)
    Requirement already satisfied: setuptools>=40.1.0 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (65.5.0)
    Collecting jupyter-lsp>=2.0.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/07/e0/7bd7cff65594fd9936e2f9385701e44574fc7d721331ff676ce440b14100/jupyter_lsp-2.2.5-py3-none-any.whl (69 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m69.1/69.1 kB[0m [31m16.0 MB/s[0m eta [36m0:00:00[0m
    [?25hCollecting async-lru>=1.0.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/fa/9f/3c3503693386c4b0f245eaf5ca6198e3b28879ca0a40bde6b0e319793453/async_lru-2.0.4-py3-none-any.whl (6.1 kB)
    Requirement already satisfied: tomli>=1.2.2 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (2.0.1)
    Collecting httpx>=0.25.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/56/95/9377bcb415797e44274b51d46e3249eba641711cf3348050f76ee7b15ffc/httpx-0.27.2-py3-none-any.whl (76 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m76.4/76.4 kB[0m [31m23.8 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: ipykernel>=6.5.0 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (6.21.2)
    Requirement already satisfied: json5>=0.9.0 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (0.9.11)
    Requirement already satisfied: babel>=2.10 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (2.12.1)
    Collecting jsonschema>=4.18.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/69/4a/4f9dbeb84e8850557c02365a0eee0649abe5eb1d84af92a25731c6c0f922/jsonschema-4.23.0-py3-none-any.whl (88 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m88.5/88.5 kB[0m [31m27.4 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: requests>=2.31 in /opt/mamba/lib/python3.10/site-packages (from jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (2.32.3)
    Requirement already satisfied: idna>=2.8 in /opt/mamba/lib/python3.10/site-packages (from anyio>=3.1.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (3.4)
    Requirement already satisfied: sniffio>=1.1 in /opt/mamba/lib/python3.10/site-packages (from anyio>=3.1.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (1.3.0)
    Requirement already satisfied: argon2-cffi-bindings in /opt/mamba/lib/python3.10/site-packages (from argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (21.2.0)
    Requirement already satisfied: typing-extensions>=4.0.0 in /opt/mamba/lib/python3.10/site-packages (from async-lru>=1.0.0->jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (4.12.2)
    Collecting httpcore==1.*
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/06/89/b161908e2f51be56568184aeb4a880fd287178d176fd1c860d2217f41106/httpcore-1.0.6-py3-none-any.whl (78 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m78.0/78.0 kB[0m [31m17.5 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: certifi in /opt/mamba/lib/python3.10/site-packages (from httpx>=0.25.0->jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (2022.9.24)
    Collecting h11<0.15,>=0.13
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/95/04/ff642e65ad6b90db43e668d70ffb6736436c7ce41fcc549f4e9472234127/h11-0.14.0-py3-none-any.whl (58 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m58.3/58.3 kB[0m [31m18.0 MB/s[0m eta [36m0:00:00[0m
    [?25hRequirement already satisfied: nest-asyncio in /opt/mamba/lib/python3.10/site-packages (from ipykernel>=6.5.0->jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (1.5.6)
    Requirement already satisfied: psutil in /opt/mamba/lib/python3.10/site-packages (from ipykernel>=6.5.0->jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (5.9.4)
    Requirement already satisfied: debugpy>=1.6.5 in /opt/mamba/lib/python3.10/site-packages (from ipykernel>=6.5.0->jupyterlab<4.3,>=4.2.0->notebook>=7->nglview) (1.6.6)
    Requirement already satisfied: parso<0.9.0,>=0.8.0 in /opt/mamba/lib/python3.10/site-packages (from jedi>=0.16->ipython>=6.1.0->ipywidgets>=8->nglview) (0.8.3)
    Requirement already satisfied: MarkupSafe>=2.0 in /opt/mamba/lib/python3.10/site-packages (from jinja2>=3.0.3->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.1.2)
    Collecting referencing>=0.28.4
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/b7/59/2056f61236782a2c86b33906c025d4f4a0b17be0161b63b70fd9e8775d36/referencing-0.35.1-py3-none-any.whl (26 kB)
    Collecting rpds-py>=0.7.1
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/4a/6d/1166a157b227f2333f8e8ae320b6b7ea2a6a38fbe7a3563ad76dffc8608d/rpds_py-0.20.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (354 kB)
    [2K     [90m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[0m [32m354.8/354.8 kB[0m [31m67.4 MB/s[0m eta [36m0:00:00[0m
    [?25hCollecting jsonschema-specifications>=2023.03.6
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/d1/0f/8910b19ac0670a0f80ce1008e5e751c4a57e14d2c4c13a482aa6079fa9d6/jsonschema_specifications-2024.10.1-py3-none-any.whl (18 kB)
    Requirement already satisfied: attrs>=22.2.0 in /opt/mamba/lib/python3.10/site-packages (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (22.2.0)
    Requirement already satisfied: python-dateutil>=2.8.2 in /opt/mamba/lib/python3.10/site-packages (from jupyter-client>=7.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.8.2)
    Requirement already satisfied: platformdirs>=2.5 in /opt/mamba/lib/python3.10/site-packages (from jupyter-core!=5.0.*,>=4.12->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (3.0.0)
    Requirement already satisfied: pyyaml>=5.3 in /opt/mamba/lib/python3.10/site-packages (from jupyter-events>=0.9.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (6.0)
    Requirement already satisfied: rfc3339-validator in /opt/mamba/lib/python3.10/site-packages (from jupyter-events>=0.9.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.1.4)
    Requirement already satisfied: python-json-logger>=2.0.4 in /opt/mamba/lib/python3.10/site-packages (from jupyter-events>=0.9.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.0.7)
    Requirement already satisfied: rfc3986-validator>=0.1.1 in /opt/mamba/lib/python3.10/site-packages (from jupyter-events>=0.9.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.1.1)
    Requirement already satisfied: mistune<3,>=2.0.3 in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.0.5)
    Requirement already satisfied: defusedxml in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.7.1)
    Requirement already satisfied: nbclient>=0.5.0 in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.7.2)
    Requirement already satisfied: tinycss2 in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (1.2.1)
    Requirement already satisfied: jupyterlab-pygments in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.2.2)
    Requirement already satisfied: bleach in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (6.0.0)
    Requirement already satisfied: pandocfilters>=1.4.1 in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (1.5.0)
    Requirement already satisfied: beautifulsoup4 in /opt/mamba/lib/python3.10/site-packages (from nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (4.11.2)
    Requirement already satisfied: fastjsonschema in /opt/mamba/lib/python3.10/site-packages (from nbformat>=5.3.0->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.16.3)
    Requirement already satisfied: ptyprocess>=0.5 in /opt/mamba/lib/python3.10/site-packages (from pexpect>4.3->ipython>=6.1.0->ipywidgets>=8->nglview) (0.7.0)
    Requirement already satisfied: wcwidth in /opt/mamba/lib/python3.10/site-packages (from prompt-toolkit!=3.0.37,<3.1.0,>=3.0.30->ipython>=6.1.0->ipywidgets>=8->nglview) (0.2.6)
    Requirement already satisfied: urllib3<3,>=1.21.1 in /opt/mamba/lib/python3.10/site-packages (from requests>=2.31->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (1.26.11)
    Requirement already satisfied: charset-normalizer<4,>=2 in /opt/mamba/lib/python3.10/site-packages (from requests>=2.31->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (2.1.1)
    Requirement already satisfied: executing>=1.2.0 in /opt/mamba/lib/python3.10/site-packages (from stack-data->ipython>=6.1.0->ipywidgets>=8->nglview) (1.2.0)
    Requirement already satisfied: pure-eval in /opt/mamba/lib/python3.10/site-packages (from stack-data->ipython>=6.1.0->ipywidgets>=8->nglview) (0.2.2)
    Requirement already satisfied: asttokens>=2.1.0 in /opt/mamba/lib/python3.10/site-packages (from stack-data->ipython>=6.1.0->ipywidgets>=8->nglview) (2.2.1)
    Requirement already satisfied: six in /opt/mamba/lib/python3.10/site-packages (from asttokens>=2.1.0->stack-data->ipython>=6.1.0->ipywidgets>=8->nglview) (1.16.0)
    Requirement already satisfied: jsonpointer>1.13 in /opt/mamba/lib/python3.10/site-packages (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (2.3)
    Requirement already satisfied: isoduration in /opt/mamba/lib/python3.10/site-packages (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (20.11.0)
    Requirement already satisfied: fqdn in /opt/mamba/lib/python3.10/site-packages (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (1.5.1)
    Collecting webcolors>=24.6.0
      Downloading https://pypi.tuna.tsinghua.edu.cn/packages/f0/33/12020ba99beaff91682b28dc0bbf0345bbc3244a4afbae7644e4fa348f23/webcolors-24.8.0-py3-none-any.whl (15 kB)
    Requirement already satisfied: uri-template in /opt/mamba/lib/python3.10/site-packages (from jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (1.2.0)
    Requirement already satisfied: cffi>=1.0.1 in /opt/mamba/lib/python3.10/site-packages (from argon2-cffi-bindings->argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (1.15.1)
    Requirement already satisfied: soupsieve>1.2 in /opt/mamba/lib/python3.10/site-packages (from beautifulsoup4->nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.4)
    Requirement already satisfied: webencodings in /opt/mamba/lib/python3.10/site-packages (from bleach->nbconvert>=6.4.4->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (0.5.1)
    Requirement already satisfied: pycparser in /opt/mamba/lib/python3.10/site-packages (from cffi>=1.0.1->argon2-cffi-bindings->argon2-cffi>=21.1->jupyter-server<3,>=2.4.0->notebook>=7->nglview) (2.21)
    Requirement already satisfied: arrow>=0.15.0 in /opt/mamba/lib/python3.10/site-packages (from isoduration->jsonschema>=4.18.0->jupyterlab-server<3,>=2.27.1->notebook>=7->nglview) (1.2.3)
    Building wheels for collected packages: nglview
      Building wheel for nglview (pyproject.toml) ... [?25ldone
    [?25h  Created wheel for nglview: filename=nglview-3.1.2-py3-none-any.whl size=7493364 sha256=625fa45918fd765e65622d859b4a095eaf46f796fd0abbfdbd744853d3c5b76f
      Stored in directory: /root/.cache/pip/wheels/1c/d6/5e/fad246784fb1e0d89794da5c9ab85dab034cb01e1250c12204
    Successfully built nglview
    Installing collected packages: widgetsnbextension, websocket-client, webcolors, send2trash, rpds-py, overrides, jupyterlab-widgets, h11, comm, async-lru, referencing, httpcore, jsonschema-specifications, httpx, jsonschema, ipywidgets, jupyter-events, jupyter-server, jupyterlab-server, jupyter-lsp, jupyterlab, notebook, nglview
      Attempting uninstall: websocket-client
        Found existing installation: websocket-client 1.5.1
        Uninstalling websocket-client-1.5.1:
          Successfully uninstalled websocket-client-1.5.1
      Attempting uninstall: webcolors
        Found existing installation: webcolors 1.12
        Uninstalling webcolors-1.12:
          Successfully uninstalled webcolors-1.12
      Attempting uninstall: send2trash
        Found existing installation: Send2Trash 1.8.0
        Uninstalling Send2Trash-1.8.0:
          Successfully uninstalled Send2Trash-1.8.0
      Attempting uninstall: comm
        Found existing installation: comm 0.1.2
        Uninstalling comm-0.1.2:
          Successfully uninstalled comm-0.1.2
      Attempting uninstall: jsonschema
        Found existing installation: jsonschema 4.17.3
        Uninstalling jsonschema-4.17.3:
          Successfully uninstalled jsonschema-4.17.3
      Attempting uninstall: jupyter-events
        Found existing installation: jupyter-events 0.6.3
        Uninstalling jupyter-events-0.6.3:
          Successfully uninstalled jupyter-events-0.6.3
      Attempting uninstall: jupyter-server
        Found existing installation: jupyter_server 2.3.0
        Uninstalling jupyter_server-2.3.0:
          Successfully uninstalled jupyter_server-2.3.0
      Attempting uninstall: jupyterlab-server
        Found existing installation: jupyterlab_server 2.19.0
        Uninstalling jupyterlab_server-2.19.0:
          Successfully uninstalled jupyterlab_server-2.19.0
      Attempting uninstall: jupyterlab
        Found existing installation: jupyterlab 3.6.1
        Uninstalling jupyterlab-3.6.1:
          Successfully uninstalled jupyterlab-3.6.1
      Attempting uninstall: notebook
        Found existing installation: notebook 6.5.2
        Uninstalling notebook-6.5.2:
          Successfully uninstalled notebook-6.5.2
    Successfully installed async-lru-2.0.4 comm-0.2.2 h11-0.14.0 httpcore-1.0.6 httpx-0.27.2 ipywidgets-8.1.5 jsonschema-4.23.0 jsonschema-specifications-2024.10.1 jupyter-events-0.10.0 jupyter-lsp-2.2.5 jupyter-server-2.14.2 jupyterlab-4.2.5 jupyterlab-server-2.27.3 jupyterlab-widgets-3.0.13 nglview-3.1.2 notebook-7.2.2 overrides-7.7.0 referencing-0.35.1 rpds-py-0.20.0 send2trash-1.8.3 webcolors-24.8.0 websocket-client-1.8.0 widgetsnbextension-4.0.13
    [33mWARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.pypa.io/warnings/venv[0m[33m
    [0mNote: you may need to restart the kernel to use updated packages.
    


```
from ase.io import read
from ase.visualize import view

atoms = read("./run_autoneb005.traj")
view(atoms, viewer='ngl')
```




    HBox(children=(NGLWidget(), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'H', 'Pt'), value…



与此同时，在`neb_post.py`处理之后 


```
# ! python3 ~/opt/ATST-Tools/neb/neb_post.py --autoneb run_autoneb???.traj
```


```
ls
```

    [0m[01;34mAutoNEB_iter[0m/        run_autoneb000.traj  run_autoneb006.traj
    autoneb_run.py       run_autoneb001.traj  run_autoneb007.traj
    autoneb_submit.sh    run_autoneb002.traj  run_autoneb008.traj
    init_neb_chain.traj  run_autoneb003.traj  run_autoneb009.traj
    neb_latest.traj      run_autoneb004.traj  running_autoneb.out
    nebplots_all.pdf     run_autoneb005.traj  [01;34mvib_analysis_TS[0m/
    

所得到的*neb_latest.traj*内有收敛后的NEB链的所有信息，可以通过ASE可视化。


```
from ase.io import read
from ase.visualize import view

atoms = read("./neb_latest.traj", ":")
view(atoms, viewer='ngl')
```




    HBox(children=(NGLWidget(max_frame=9), VBox(children=(Dropdown(description='Show', options=('All', 'C', 'H', '…



且所生成的`nebplots_all.pdf`很好地给出了当前的NEB反应路径对应的势能-反应坐标曲线，从此时的曲线中不难看出如此计算得到的过渡态已然正确，并且过渡态邻近的点更为靠近过渡态。

<div>
    <img src="https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/13504/8259af60e2fb45afa29876c786449254/4d6d1414-f21e-486b-bfa5-15bf616832e5.png", alt="reaction_diffusion" width="450" title="reaction_diffusion">
    <p style='font-size:0.8rem; font-weight:bold'> Cy-Pt@graphene 案例经AutoNEB方法得到的正确计算结果 </p>
</div>


振动计算（频率分析）表明该过渡态确实只在反应坐标方向存在单一虚频，具体结过将在后续展示。

经过以上这些例子，相信你已经完全掌握了基于ATST-Tools的组织方法，采用ASE-ABACUS接口完成NEB以及AutoNEB过渡态搜索计算了。使用dimer等其他方法也是类似的道理。

### 2.3 振动计算和自由能较正

在计算均相或多相的化学反应时，我们在通过结构优化和过渡态搜索方法，得到反应的初末态和过渡态的能量和电子结构信息的基础之上，还会关心：
- 我们得到的初末态优化结构是否是稳定结构？即该优化结构是否具有虚频
- 我们得到的过渡态是否是正确的过渡态结构？即该过渡态是否具有沿反应路径的单一虚频
- 我们能否考虑除零温电子能量之外的其他统计热力学信息，比如零点能，以及适当的振动自由能校正。

这些事情需要通过振动计算（又称声子计算）来完成。其方法通常有两种，一是线性响应法，二是有限差分法。其中，最常用的是有限差分法，它通过生成一系列的微扰结构，基于结构间各原子的坐标差和受力差，通过有限差分的方法求得对应原子的Hessian贡献，进而求得振动模式。

ATST-Tools中集成了使用ASE-ABACUS基于有限差分进行振动计算的方法。该工作流实际上只使用了一个脚本，我们可以直接查看


```
cat /opt/ATST-Tools/vibration/vib_analysis.py
```

    # JamesMisaka in 2023-11-30
    # Vibrational analysis from finite displacement by using abacus
    # part of ATST-Tools scripts
    
    import os
    import numpy as np
    from ase.vibrations import Vibrations
    from ase.thermochemistry import HarmonicThermo
    from ase.io import read, write
    from ase.calculators.abacus import Abacus, AbacusProfile
    from ase.parallel import world, parprint
    from ase.mep.neb import NEBTools
    from neb2vib import neb2vib
    
    # reading usage can be changed by user
    # usage by neb
    neb_traj = read('neb_latest.traj', index=':')
    atoms, vib_indices = neb2vib(neb_traj)
    
    # traditional
    # stru = "STRU"
    # atoms = read(stru)
    # vib_indices = [0, 1, 37]
    
    # indices setting for which atoms to be displaced
    # vib_indices = [atom.index for atom in atoms if atom.symbol == 'H']
    
    T = 523.15 # K
    
    abacus = "abacus"
    mpi = 16
    omp = 4
    lib_dir = "/lustre/home/2201110432/example/abacus"
    pseudo_dir = f"{lib_dir}/PP"
    basis_dir = f"{lib_dir}/ORB"
    # default pp and basis is supported by ase-abacus interface
    pp = {
          'H': 'H_ONCV_PBE-1.0.upf',
          'C': 'C_ONCV_PBE-1.0.upf',
          'O': 'O_ONCV_PBE-1.0.upf',
          'Fe': 'Fe_ONCV_PBE-1.0.upf',
          }
    basis = {
             'H': 'H_gga_6au_100Ry_2s1p.orb',
             'C': 'C_gga_7au_100Ry_2s2p1d.orb',
             'O': 'O_gga_7au_100Ry_2s2p1d.orb',
             'Fe': 'Fe_gga_8au_100Ry_4s2p2d1f.orb',
             }
    kpts = [3, 1, 2]
    parameters = {
        'calculation': 'scf',
        'nspin': 2,
        'xc': 'pbe',
        'ecutwfc': 100,
        'ks_solver': 'genelpa',
        'symmetry': 0,
        'vdw_method': 'none',
        'smearing_method': 'mp',
        'smearing_sigma': 0.002,
        'basis_type': 'lcao',
        'mixing_type': 'broyden',
        'mixing_ndim': 8,
        'scf_thr': 1e-7,
        'scf_nmax': 300,
        'kpts': kpts,
        'pp': pp,
        'basis': basis,
        'pseudo_dir': pseudo_dir,
        'basis_dir': basis_dir,
        'init_chg': 'file',
        'init_wfc': 'atomic',
        'cal_force': 1,
        'cal_stress': 1,
        'out_stru': 1,
        'out_chg': 1,
        'out_mul': 0,
        'out_bandgap': 0,
        'out_wfc_lcao': 0,
        'efield_flag': 1,
        'dip_cor_flag': 1,
        'efield_dir': 1,
    }
    
    # developer only
    vib_name = 'vib'
    delta = 0.01
    nfree = 2
    
    def set_calculator(abacus, parameters, mpi=1, omp=1) -> Abacus:
        """Set Abacus calculators"""
        os.environ['OMP_NUM_THREADS'] = f'{omp}'
        profile = AbacusProfile(f"mpirun -np {mpi} {abacus}")
        out_directory = f"SCF-rank{world.rank}"
        calc = Abacus(profile=profile, directory=out_directory,
                    **parameters)
        return calc
    
    
    if __name__ == "__main__":
        print("==> Starting Vibrational Analysis <==")
    
        atoms.calc = set_calculator(abacus, parameters, mpi=mpi, omp=omp)
    
        vib = Vibrations(atoms, indices=vib_indices, 
                        name=vib_name, delta=delta, nfree=nfree)
    
        print("==> Running Vibrational Analysis <==")
        vib.run()
        # post-processing
        print("==> Done !!! <==")
        print(f"==> All force cache will be in {vib_name} directory <==")
        print("==> Vibrational Analysis Summary <==")
        vib.summary()
        print("==> Writing All Mode Trajectory <==")
        vib.write_mode()
        # thermochemistry
        print("==> Doing Harmonmic Thermodynamic Analysis <==")
        vib_energies = vib.get_energies()
        #real_vib_energies = np.array([energy for energy in vib.get_energies() if energy.imag == 0 and energy.real > 0], dtype=float)
        thermo = HarmonicThermo(vib_energies, ignore_imag_modes=True,)
        entropy = thermo.get_entropy(T)
        free_energy = thermo.get_helmholtz_energy(T)
        print(f"==> Entropy: {entropy:.6e} eV/K <==")
        print(f"==> Free Energy: {free_energy:.6f} eV <==")
        
    
    

最关键的一点在于选择需要进行有限差分和计算振动模的原子。实际上，严格的振动模计算肯定是需要对晶胞内所有原子进行有限差分的，但由于针对表面催化体系，我们实际上只关心那些直接参与反应的部分原子的振动模，从而我们只需要选择一部分原子即可。

在选定需要进行有限差分计算的原子序号，以及使用ASE进行振动校正与自由能计算，和调用ABACUS计算的相关参数之后，可以直接通过
```
python3 vib_analysis.py
```
调用ABACUS进行振动计算。

我们可以直接看看针对上一解离反应过渡态的计算结果


```
cd /opt/ATST-Tools/examples/Cy-Pt@graphene/autoneb/vib_analysis_TS
```

    /opt/ATST-Tools/examples/Cy-Pt@graphene/autoneb/vib_analysis_TS
    


```
ls
```

    JobRun.state      vib.14.traj  vib.26.traj  vib.38.traj  vib.5.traj
    ase_sort.dat      vib.15.traj  vib.27.traj  vib.39.traj  vib.50.traj
    neb_latest.traj   vib.16.traj  vib.28.traj  vib.4.traj   vib.51.traj
    [0m[01;34mneb_latest_CIFs[0m/  vib.17.traj  vib.29.traj  vib.40.traj  vib.52.traj
    running_vib.err   vib.18.traj  vib.3.traj   vib.41.traj  vib.53.traj
    running_vib.out   vib.19.traj  vib.30.traj  vib.42.traj  vib.6.traj
    [01;34mvib[0m/              vib.2.traj   vib.31.traj  vib.43.traj  vib.7.traj
    vib.0.traj        vib.20.traj  vib.32.traj  vib.44.traj  vib.8.traj
    vib.1.traj        vib.21.traj  vib.33.traj  vib.45.traj  vib.9.traj
    vib.10.traj       vib.22.traj  vib.34.traj  vib.46.traj  vib_analysis.py
    vib.11.traj       vib.23.traj  vib.35.traj  vib.47.traj
    vib.12.traj       vib.24.traj  vib.36.traj  vib.48.traj
    vib.13.traj       vib.25.traj  vib.37.traj  vib.49.traj
    

计算结果包括：
- `running_vib.out` 即标准输出，输出了频率计算过程与结果
- `vib.*.traj` 以轨迹形式输出了计算所得的各个振动模式的动画
- `vib`文件夹，里面以JSON格式存储了计算所得的各个不同结构的能量与原子受力，用于程序复用，在下一次输出振动分析结果，或更换温度等条件进行自由能校正计算时不需要重新调用ABACUS进行计算

至于`neb_latest.traj`则是上一步AutoNEB计算所得的NEB轨迹，可以通过ATST-Tools中的如下脚本将其拆分为CIF等格式的结构文件



```
# 先remove一下
! rm -rf neb_latest_CIFs/
```


```
! python3 /opt/ATST-Tools/neb/traj_transform.py neb_latest.traj cif
```

    ===> Successfully Transform NEB Traj to cif files ! <===
    

查看`running_vib.out`，我们即可获取关于此结构进行振动计算和自由能校正的所得信息


```
cat running_vib.out
```

    Loading mkl version 2023.0.0
    Loading tbb version 2021.8.0
    Loading compiler-rt version 2023.0.0
    Loading mpi version 2021.8.0
    Loading compiler version 2023.0.0
    Loading oclfpga version 2023.0.0
      Load "debugger" to debug DPC++ applications with the gdb-oneapi debugger.
      Load "dpl" for additional DPC++ APIs: https://github.com/oneapi-src/oneDPL
    
    Loading abacus/3.4.3-icx-dev
      Loading requirement: tbb/latest compiler-rt/latest mkl/latest mpi/latest
        oclfpga/latest compiler/latest
    Vibrational Calculation Start !
    ==> Running Vibrational Analysis <==
    ==> Get Energy and Frequency <==
    ==> Get ZPE <==
    ==> Writing Trajectory <==
    ==> Summary <==
    ---------------------
      #    meV     cm^-1
    ---------------------
      0   87.4i    705.0i
      1    2.4      19.3
      2    3.5      28.3
      3    6.1      48.9
      4   14.5     117.2
      5   17.4     140.7
      6   20.9     168.5
      7   30.4     244.9
      8   38.7     312.0
      9   51.0     411.6
     10   53.1     428.5
     11   57.8     466.1
     12   65.3     526.4
     13   79.5     641.3
     14   96.8     781.1
     15   98.2     791.9
     16  102.3     825.1
     17  107.3     865.6
     18  108.3     873.4
     19  111.7     900.6
     20  120.9     975.1
     21  123.8     998.4
     22  129.7    1045.8
     23  131.5    1060.4
     24  132.7    1070.3
     25  137.3    1107.4
     26  144.5    1165.2
     27  151.3    1220.2
     28  153.9    1241.2
     29  155.6    1255.2
     30  158.1    1275.4
     31  159.4    1285.7
     32  163.2    1316.1
     33  163.7    1320.7
     34  165.8    1337.5
     35  166.6    1343.8
     36  167.1    1348.0
     37  178.3    1438.0
     38  178.9    1442.9
     39  179.2    1445.6
     40  179.8    1450.5
     41  181.3    1462.1
     42  263.0    2121.5
     43  362.0    2920.0
     44  366.9    2959.6
     45  367.8    2966.6
     46  368.1    2968.8
     47  369.0    2976.1
     48  371.0    2992.1
     49  373.5    3012.7
     50  374.4    3019.4
     51  374.7    3022.4
     52  375.8    3031.3
     53  377.5    3044.9
    ---------------------
    Zero-point energy: 4.416 eV
    ===== Done at Wed Nov 29 23:13:33 CST 2023! =====
    

可以看到该结构确实只具有单一虚频。查看`vib.0.traj`轨迹文件，可证实该振动模式沿反应坐标方向。于是我们便完成了振动计算和自由能校正。

## 3. 总结

通过本Notebook，你学到了：
1. 过渡态理论与过渡态搜索的基础知识
2. 通过ATST-Tools，基于ASE和ASE-ABACUS接口完成以NEB和AutoNEB为例的过渡态搜索计算，包括预处理，具体计算与后处理过程。
3. 对计算所得过渡态进行振动分析和自由能较正，确保过渡态正确，且得到其包含有振动自由度项的自由能较正结果。

如果你对过渡态与过渡态搜索产生了更多兴趣，或者有更多见解，欢迎一起讨论!

