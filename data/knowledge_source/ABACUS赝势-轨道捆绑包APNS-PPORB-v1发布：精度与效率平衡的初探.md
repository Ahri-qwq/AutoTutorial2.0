# ABACUS赝势-轨道捆绑包APNS-PPORB-v1发布：精度与效率平衡的初探

> 来源: https://mp.weixin.qq.com/s/vUG-7uRjbMZJP2fs8Ep_8Q

ABACUS赝势轨道库项目（ABACUS Pseudopotential-NAO Square, APNS）已经开展数月，在对常见模守恒赝势，以及ABACUS新实现的对超软赝势进行大规模、高通量的效率和精度测试之外，也实现了轨道生成算法的部分改进。在经过兼顾效率与精度的甄选后，ABACUS-AISI开发团队推出了APNS-PPORB-v1版本赝势-轨道捆绑包，方便用户在仅需进行少量测试的情况下，更加放心地使用ABACUS对感兴趣体系进行计算，得到尽可能贴近全电子精度的结果。

本次捆绑包支持了从H到Bi共69个元素（除La外不含其他镧系、锕系元素）的赝势和对应ABACUS数值原子轨道。轨道分为efficiency和precision两个合集。efficiency系列的推荐使用场景为：

结构优化

过渡态搜索

声子谱计算

弹性模量计算

势能面扫描

分子动力学模拟

等仅需要占据态决定的性质。precision系列的推荐使用场景为：

高精度能量计算

能带计算（费米面以上能带精度有较高要求时）

TD-DFT和激发态相关计算

等需要轨道具有较高完备性，或同时需要占据态和非占据态决定的性质计算。

使用该捆绑包的学术引用方式和具体使用方法请参考readme.txt文件。

**注：**当前版本赝势-轨道捆绑包不支持旋轨耦合效应（Spin-orbit coupling, SOC）的计算，后续将单独提供面向该功能的赝势-轨道捆绑包版本。

目前赝势-轨道捆绑包在AIS-Square开放下载：https://aissquare.com/datasets/detail?pageType=datasets&name=ABACUS-APNS-PPORBs-v1%253Apre-release&id=326

该界面将随时有试验版本赝势-轨道捆绑包动态更新。

在DFT的Hamiltonian中，赝势代替了−Z/∣**r**−**R**∣项，因此元素的性质将全部由赝势所决定。赝势的质量和使用的正确性将直接决定计算结果的意义大小。

参数ecutwfc对应于展开波函数的平面波最大动能，单位是Ry，同时在模守恒赝势的情况下，常有ecutwfc * 4 = ecutrho，ecutrho为实空间电荷密度及其他实空间格点数据的格点密度控制参数。ecutwfc的正确设置对得到正确的计算结果十分重要（实际上，赝势很难同时实现“ecutwfc足够低”和“精度足够高”两个理想目标，因此对于部分精度较高的赝势，其ecutwfc很可能超过了常规取值，需要格外注意）。

APNS项目开展早期已经开放了使用

1. 原子平均Kohn-Sham能量ΔEKS

2. 晶格压强ΔP

3. 占据态能带相似程度η00

三个指标，对Materials Project在线开源材料数据库中热力学最稳定单质进行了赝势的ecutwfc收敛性测试。三指标的收敛标准分别为：

1. ΔEKS ≤ 1 meV/atom

2. ΔP ≤ 0.1 kbar

3. η 00 ≤ 10 meV

数据展示在APNS专属Github Pages：https://kirk0830.github.io/ABACUS-Pseudopot-Nao-Square/pseudopotential/pseudopotential.html

数据开源在AIS-Square科学智能广场：https://aissquare.com/datasets/detail?pageType=datasets&name=ABACUS-ncpp-efficiency-precision-testdata&id=288

注：自ABACUS v3.10 LTS发布以来，由于对称性模块的BUG修复和对角化算法的迭代，部分元素的收敛ecutwfc数值有所下降：

https://aissquare.com/datasets/detail?pageType=datasets&name=ABACUS-ncpp-efficiency-stress-update-20250311&id=314

赝势作为全电子的近似方法，其精度应当使用全电子数据进行对比评估。基于上述测得各赝势的收敛ecutwfc，我们使用了Bosoni等人[1]DFT软件实现的精度标定工作规范中开源的全电子数据，对所有元素均构建其体心立方、面心立方和金刚石单质，以及所有元素的六种氧化物晶体（记元素为X，则六种氧化物分别为：X2O, XO, X2O3, XO2, X2O5, XO3），进行从体积94%到106%范围的扫描，得到EOS曲线和全电子结果进行对比（计算delta值）。delta值用于衡量两EOS曲线间差异，亦即两EOS曲线数据点的均方根：

数据开源在AIS-Square科学智能广场：https://aissquare.com/datasets/detail?pageType=datasets&name=ABACUS-stable-psp-orb-test-data&id=318

如果赝势具有越高的ecutwfc收敛值，则其具有越低的“效率”（更高的计算成本）。相反如果赝势（相对于全电子计算结果）具有越低的EOS delta值，则其具有越高的精度。正如前文提到，对于所有元素，在众多赝势中，我们挑选可以使得ecutwfc按照标准收敛在≤100 Ry，但delta值尽可能小于1 meV/atom。具体筛选算法如下：

得到筛选后结果如下（数字为ecutwfc收敛值，颜色对应九种体系delta值平均）：

注：Hf和Ta元素赝势在氧化物体系中表现出精度稍差，体现出较弱的可迁移性，但在低氧化态和单质中表现仍比较可观（~1 meV/atom）：

和高斯型轨道（Gaussian Type Orbital, GTO）不同的是，数值原子轨道（Numerical Atomic Orbital, NAO）具有人为定义的截断半径，其取值同样面临精度和效率的抉择。为避免用户在实际使用过程中在截断半径的选择上有所疑虑，我们对H-Bi（除La外不含镧系元素）共69个元素进行了DZP（设赝势中价电子的组态为2s1p，则DZP为4s2p1d，即将各角动量的zeta函数数量乘2后增加一个角动量更高的zeta函数）基组，截断半径取值从6到10 au的EOS测试，分别相对于PW和全电子得到如下结果（数字为九种体系delta值平均，括号内为标准差）：

注1：元素Te, I, Xe, Cs呈现出较大的delta值与方差，体现DZP轨道的迁移性有限。这四个元素的具体数据为（相对于全电子）：

基于轨道的具体筛选算法如下：

最终得到轨道筛选结果如下（第一个数字为相对于PW的delta值，第二个数字为相对于全电子的delta值）：

efficiency系列轨道delta值测定结果与截断半径选取情况

注2：我们为Cs元素提供了10au和12au的轨道，供用户根据实际体系情况自主选择。

考虑到部分用户关心表面的分子吸附能量计算精度或距离费米面较远的非占据能带，这实际上对基组的完备性提出了更高的要求。前者相关的计算结果误差称为基组重叠误差（Basis-Set-Superposition-Error, BSSE）。使用较大的截断半径，能够降低该方面造成的误差。因此我们在轨道的precision集合中，统一将截断半径设置为10 au。更加全面地降低BSSE可以参考官网中counterpoise方法校正BSSE的步骤说明。以下展示precision系列轨道在能带计算任务中的误差测试结果，其中所使用指标 η 的定义为

，其中δ 为费米面平移量，δ 为smearing展宽，g(σ) 为展宽函数， 为两方法（A、B）计算能带（

和

）的占据数几何平均数，


。ω 用于最大化对齐两能带，以避免不同赝势之间的能量平移。

以下直接展示部分元素单质的能带计算结果（TZDP和PW对比）：

注：为保证能带计算精度，部分元素（B, Na, Mg, P, S, K, Se, Br, Cs, Sb, Te, I, Rb）的轨道使用了尚未发表的轨道生成算法生成，其各角动量的zeta函数数量构建形式暂时参考Dunning等人的“相关一致性”基组（cc-pVTZ(-f), cc-pVTZ和cc-pVQZ(-g)），下图分别展示了这三种基组对B和Se体系的能带精度的影响。

注2：通过修改对应于cc-pVTZ level轨道文件中的文件头，可以实现cc-pVTZ(-f)的使用：

B_gga_10au_100Ry_3s3p2d0f.orb

`B_gga_10au_100Ry_3s3p2d0f.orb`

`---------------------------------------------------------------------------`

`Element B`

`Energy Cutoff(Ry) 100`

`Radius Cutoff(a.u.) 10`

`Lmax 3`

`Number of Sorbital--> 3`

`Number of Porbital--> 3`

`Number of Dorbital--> 2`

`Number of Forbital--> 0`

`---------------------------------------------------------------------------`

`SUMMARY END`


注3：如果需要其他截断半径的轨道，请在ABACUS Github仓库以issue（https://github.com/deepmodeling/abacus-develop/issues）的形式提交需求。

AI for Science时代的到来，使得数据精度保证和计算效率提升的重要性与日俱增。ABACUS秉承算法开源、代码开源、测试数据开源的原则，在不断优化代码算法、代码结构的同时，着眼于帮助用户解决实际模拟问题。ABACUS开发团队将以此次赝势-轨道捆绑包v1版本的发布为起点，未来持续迭代，充分发挥ABACUS数值原子轨道基组的特色优势，使得用户使用ABACUS进行感兴趣体系相关性质计算时的效率-精度平衡能够登上新的台阶。

ABACUS赝势轨道库项目（一期、二期）得到北京科学智能研究院（AISI）院内培育项目资金支持，能带测试结果由DeePTB项目组友情提供。

[1] Bosoni E, Beal L, Bercx M, et al. How to verify the precision of density-functional-theory implementations via reproducible and universal workflows[J]. Nature Reviews Physics, 2024, 6(1): 45-58.