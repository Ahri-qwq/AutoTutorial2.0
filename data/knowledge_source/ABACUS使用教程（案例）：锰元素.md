# 前言
本notebook记录了笔者使用ABACUS完成对锰元素进行平面波基组计算、收敛性测试和结构弛豫的相关流程，希望能对读者有所帮助


# 元素介绍
## 锰元素
### 基本属性
- 原子序数：25  
- 原子量：54.938  
- 电子排布：[Ar]  3d^5 4s^2  
- 周期与族：第四周期，VIIB族（锰族）  
- 电负性：1.55  
- 常见氧化态：+2（Mn^(2+)盐）、+4（MnO_2）、+6（锰酸盐）、+7（高锰酸盐）。  

### 发现历史
1774年：瑞典矿物学家约翰·甘恩通过还原软锰矿（MnO_2）与碳首次制得金属锰。  
命名：源自拉丁语“Magnes”（磁石），因早期与磁铁矿混淆。

###物理性质
- 状态：银白色金属，硬而脆。  
- 熔点：1246°C  
- 沸点：2061°C  
- 密度：7.21 g/cm³  
- 磁性：顺磁性，某些合金（如Heusler合金）具铁磁性。  

### 化学性质
- 多价态特性：从-3到+7，以+2、+4、+7最常见。  
- 强氧化性：高锰酸钾（KMnO_4，+7价）在酸性条件下生成Mn^(2+)，中性生成MnO_2，碱性生成MnO_4^(2-)。  

### 同位素
唯一稳定同位素：Mn^55（100%）。  
放射性同位素：Mn^53（半衰期370万年，用于陨石年代测定）。  

### 应用领域
- 钢铁工业：脱氧剂（减少FeO）、合金元素（锰钢含12-14% Mn，耐磨，用于铁轨、挖掘机铲斗）。  
- 电池：二氧化锰作碱性电池和锂锰电池正极材料。  
- 化工：高锰酸钾用于水处理、消毒；硫酸锰作为肥料和动物饲料添加剂。  
- 陶瓷与玻璃：二氧化锰用于着色（紫色）和脱色（氧化铁杂质）。  

### 生物作用
- 必需微量元素：参与超氧化物歧化酶（SOD）催化超氧自由基分解，保护细胞。  
- 酶辅助因子：精氨酸酶、丙酮酸羧化酶等依赖Mn^(2+)。  
- 毒性：过量吸入锰粉尘导致“锰中毒”，症状类似帕金森病。  

### 环境与健康
- 污染来源：采矿、冶金排放，锰化合物通过空气、水传播。  
- 饮用水标准：WHO建议限值0.4 mg/L，防止神经毒性。  

### 晶体结构特性  

#### 金属锰的同素异形体  
- α-Mn（复杂立方，<727°C）**  
   空间群：I-43m  
   晶格常数：a = 8.894 Å  
   单位晶胞原子数：58个  
   配位环境：4种不同原子位置，配位数7-12。  

- β-Mn（立方晶系，727–1095°C）  
   空间群：P4_1 32  
   晶格常数：a = 6.315 Å  
   单位晶胞原子数：20个  
   结构特点：螺旋链状排列，较α-Mn简单。  

- γ-Mn（面心立方，1095–1133°C）  
   空间群：Fm-3m  
   晶格常数：a = 3.86 Å  
   原子密堆积：配位数12，高温下延展性增强。  

-δ-Mn（体心立方，1133–1244°C）
   空间群：Im-3m  
   晶格常数：a = 3.08 Å  
   松散排列：导热性优异，接近熔点。

#### 锰化合物的晶体结构参数  
- 二氧化锰（MnO_2，金红石型）
   空间群：P4_2/mnm  
   晶格常数：  
     a = 4.40 Å  
     c = 2.87 Å  
   键长：Mn-O键长1.89 Å（八面体配位）。

- 氧化锰（MnO，岩盐结构） 
   空间群：Fm-3m  
   晶格常数：a = 4.444 Å  
   配位：Mn^(2+)与O^(2-)形成八面体配位。

- 四氧化三锰（Mn_3 O_4，尖晶石结构）
   空间群：l4_1 amd  
   晶格常数：  
     a = 5.76 Å  
     c = 9.46 Å  
   混合价态：Mn^(2+)（四面体位）和 Mn^(3+)（八面体位）。  

- 高锰酸钾（KMnO_4，正交晶系）  
   空间群：Pnma  
   晶格常数：  
     a = 9.04 Å  
     b = 5.72 Å  
     c = 7.39 Å  
   结构：MnO_4^-四面体与K^+通过离子键连接。  

- 硫酸锰（MnSO_4，单斜晶系）  
   空间群：P2_1/c  
   晶格常数：  
     a = 7.92 Å  
     b = 8.62 Å  
     c = 5.21 Å  
     β = 98.5° 

# 密度泛函理论简介
密度泛函理论(Density functional theory ，DFT)是一种研究多电子体系电子结构的方法，是目前凝聚态物理计算材料学和计算化学领域最常用的方法之一。

## 主要目的
研究原子结构/多原子结构也就是要解对应的定态薛定谔方程：
$$
\hat{H}\psi(\vec{r})=E\psi(\vec{r})
$$
但是对于一个存在着多个电子、原子核，甚至彼此之间还有相互作用的体系，其哈密顿算符非常复杂乃至于无法直接求解，因此需要通过近似方法来进行求解。

## Hohenberg-Kohn定理

Hohenberg-Kohn定理是密度泛函理论的基础，它由两部分组成

定理一：对任意相互作用的电子系统，处在外部势场中，则该外部势场可由基态电子密度唯一决定（除了相差一个常数）

定理二：对任意给定的外势，可定义光宇电子密度的普适性泛函，其中，基态（非简并）电子密度使得此能量泛函取极小值

Hohenberg-Kohn定理证明了存在着一个联系体系能量与其电子密度分布的普适性密度泛函，但是没有给出这个泛函是什么，事实上这一泛函的具体形式至今未知。

## Kohn-Sham方程
Kohn-Sham方程在进一步近似的前提下给出了Hohenberg-Kohn定理所预言的泛函的表达式，它标志着密度泛函理论走向实用。
Kohn-Sham方程基于Kohn-Sham假设：有相互作用系统的基态电荷密度，可以被一个无相互作用系统的基态电荷密度表示出来。这一假设使得多电子系统被转化为了无相互作用的单电子系统，从而能够较为简便地计算。但是应当注意的是，Kohn-Sham假设目前仍然停留在假设的状态，未被证实也未被证伪。

Kohn-Sham方程的具体形式为：$$H_{KS}\psi_i(\textbf{r})=\varepsilon_i \psi_i(\textbf{r})$$

其中哈密顿量为
$$H_{KS}=-\frac{1}{2}\nabla^2+V_{Hartree}(\textbf{r})+V_{ext}(\textbf{r})+V_{xc}(\textbf{r})$$
其中$V_{ext}$为外势，$V_{Hartree}=\int \frac{\rho(r)}{\left|r-r^{\prime}\right|} d r^{\prime}$为电子库伦排斥势，$V_{xc}$为交换关联势，可用交换关联泛函来近似。

需要注意的是，Kohn-Sham方程的本征值并不严格等于体系在基态的能量。实际的基态能量需要先从Kohn-Sham方程的能量本征值减去重复计算的$V_{Hartree}$项，再加上真正的Hartree能量。
## 求解KS方程的电子自洽迭代方法
电子自洽迭代(Self-Consistent Field, SCF)方法的基本思想为：
1. 随机给出一个基态电荷密度
2. 由基态电荷密度计算有效势能，代入求解Kohn-Sham方程，计算出新的基态电荷密度
3. 判断新的基态电荷密度和原来的基态电荷密度相差多少，如果相差小于阈值则输出结果，否则继续用新得到的基态电荷密度重复2&3两个步骤进行迭代

解出Kohn-Sham方程也就意味着得到了真实体系对应的虚拟体系的电子密度分布与能量分布，能够一定程度上反映真实体系的状态。

# ABACUS安装
参照[ABACUS中文文档](https://mcresearch.github.io/abacus-user-guide/)进行安装，也可以在Bohrium中加载装有ABACUS镜像的节点进行使用。

前置：已有WSL Ubuntu 22.04系统、Git、Nvidia显卡和cuda
- 安装同时支持两种基组并支持CPU&GPU版本
1. 安装cmake
```
sudo apt update
sudo apt install cmake
```
2. 克隆仓库
```
git clone https://github.com/deepmodeling/abacus-develop.git
cd abacus-develop/
```
3.	安装PW基组依赖库
```
sudo apt update 
sudo apt install -y libopenblas-openmp-dev
sudo apt install -y liblapack-dev 
sudo apt install -y libfftw3-dev
```
4.	安装MPI library
```
sudo apt install -y libopenmpi-dev
```
5.	安装LCAO基组依赖
```
sudo apt install -y libscalapack-mpi-dev
sudo apt install -y libcereal-dev
```
6.	安装ELPA
```
sudo apt install -y libelpa-dev
```
7.	编译安装
```
cmake -B build
cd build && make -j`nproc`
```
8. GPU版本编译
```
cmake -B build -DUSE_CUDA=1
cd build && make -j`nproc`
```
可能遇到的问题：
- 退出软件配置/须知：可尝试按Tab键 切换焦点到确认按钮后再按Enter
- 无法识别abacus命令：export PATH=/abacus-develop/build/:$PATH


# 文件准备

使用ABACUS进行计算最重要的是需要准备INPUT、KPT、STRU三个文件，详细参数和说明可以参考
<https://abacus.deepmodeling.com/en/latest/advanced/input_files/index.html>

此外，在计算过程中可能会用到原子的赝势和数值原子轨道文件(lcao，原子轨道线性组合计算需要用到轨道文件；pw，平面波基组计算不需要用到轨道文件)，赝势和数值原子轨道都可以在<https://abacus.ustc.edu.cn/pseudo/list.htm> 中下载。赝势生成时,是解了个一维的kohn-sham方程,解方程时是需要指定交换关联泛函的,所以赝势文件天然会带着交换关联泛函，目前最常用的是PBE交换关联泛函。

三个文件的一组例子：

```INPUT```文件：
```
INPUT_PARAMETERS
suffix       Mn           #输出后缀
ntype        1                                    #元素种类
ecutwfc      20                                  #展开截止能量
scf_thr      1e-7                                 #电荷密度收敛阈值
basis_type   pw                                 #基函数类型  lcao/pw
calculation  scf                           #计算类型
relax_nmax   100
pseudo_dir   ../../       #赝势文件
orbital_dir  ../../       #轨道文件
out_stru     1
```
```KPT```文件
```
K_POINTS    
0           
Gamma     
2 2 2 0 0 0
```
```STRU```文件
```
ATOMIC_SPECIES
Mn 54.9380 Mn_ONCV_PBE-1.2.upf        # 名称; 相对原子质量; 赝势文件名

NUMERICAL_ORBITAL
Mn_gga_10au_100Ry_2s1p1d.orb         # 轨道文件名

LATTICE_CONSTANT
1.8897259886                         # 晶格常数的单位 (Bohr), 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
3.860 0.000 0.000
0.000 3.860 0.000
0.000 0.000 3.860

ATOMIC_POSITIONS
Direct                             # 笛卡尔坐标或直接坐标
Mn                                 # 元素名称
0.0                                # 磁矩
4                                  # 原子个数
0.000 0.000 0.000 0 0 0            # 原子位置和自由度
0.000 0.500 0.500 0 0 0
0.500 0.000 0.500 0 0 0
0.500 0.500 0.000 0 0 0
```

# 

# 收敛性测试

## Ecut收敛性测试
将INPUT文件中的ecutwfc参数从20逐步增大到100，步长为5，观察总能量的变化情况，利用python脚本自动提交和运行脚本。kpoints设为6
```python
import os
import re
import matplotlib.pyplot as plt

# 原始 INPUT 文件路径
input_file_path = "./Mn/INPUT"
abacus_command = "mpirun -np 2 abacus"  #2核心运行
working_directory = "./Mn"  # 指定工作目录

# 定义路径和参数
output_dir = "./Mn"  # 输出目录
ecutwfc_range = range(20, 101, 5)  # ecutwfc 的范围
energy_data = []  # 用于存储 ecutwfc 和对应的总能量

# 读取原始 INPUT 文件内容
with open(input_file_path, "r") as file:
    input_content = file.readlines()

# 修改 ecutwfc 和 suffix，运行 abacus，并提取总能量
for ecutwfc in ecutwfc_range:
    # 修改 INPUT 文件中的 ecutwfc 和 suffix
    new_content = []
    for line in input_content:
        if line.startswith("ecutwfc"):
            new_content.append(f"ecutwfc      {ecutwfc}                                  # 展开截止能量\n")
        elif line.startswith("suffix"):
            new_content.append(f"suffix       Mn_ecutwfc_{ecutwfc}                        # 输出后缀\n")
        else:
            new_content.append(line)
    
    # 写入修改后的内容到原始 INPUT 文件
    with open(input_file_path, "w") as file:
        file.writelines(new_content)
    
    # 在指定目录运行 abacus 脚本
    os.system(f"cd {working_directory} && {abacus_command}")
    
    # 提取 running_scf.log 中的 "final etot" 总能量
    log_file = os.path.join(output_dir, f"OUT.Mn_ecutwfc_{ecutwfc}/running_scf.log")
    if not os.path.exists(log_file):
        print(f"警告: 找不到文件 {log_file}")
        continue

    # 读取 running_scf.log 文件
    with open(log_file, "r") as file:
        lines = file.readlines()

    # 提取 "final etot" 的总能量值
    final_etot = None
    for line in lines:
        match = re.search(r"final etot is ([\-\d\.]+) eV", line)
        if match:
            final_etot = float(match.group(1))
            break

    if final_etot is not None:
        energy_data.append((ecutwfc, final_etot))
    else:
        print(f"警告: 在文件 {log_file} 中未找到 'final etot' 数据")

# 保存统计数据为 CSV 文件
output_csv = os.path.join(output_dir, "final_etot_vs_ecutwfc.csv")
with open(output_csv, "w") as file:
    file.write("ecutwfc,final_etot (eV)\n")
    for ecutwfc, final_etot in energy_data:
        file.write(f"{ecutwfc},{final_etot:.6f}\n")

print(f"\n统计结果已保存到 {output_csv}")

# 绘制折线图
ecutwfc_values, final_etot_values = zip(*energy_data)
plt.figure(figsize=(8, 6))
plt.plot(ecutwfc_values, final_etot_values, marker='o', linestyle='-', color='b', label='Final Total Energy')
plt.xlabel('ecutwfc')
plt.ylabel('Final Total Energy (eV)')
plt.title('Final Total Energy vs ecutwfc')
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(output_dir, "final_etot_vs_ecutwfc.png"))
plt.show()
```


```
#运行脚本
#%cd personal/convergence_test/
!python ecut.py
```

                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:07:23 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_20/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 24 * 24 * 24
     UNIFORM GRID DIM(BIG)   : 24 * 24 * 24
     DONE(0.0606174  SEC) : SETUP UNITCELL
     DONE(0.201874   SEC) : SYMMETRY
     DONE(0.495221   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.49668    SEC) : INIT PLANEWAVE
     DONE(0.49864    SEC) : LOCAL POTENTIAL
     DONE(0.523958   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 3.6377
     DONE(0.524031   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.770851   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.05023466e+04   0.00000000e+00   1.3582e+01   2.63
     CG2     -1.05544192e+04  -5.20726482e+01   1.7102e+00   1.20
     CG3     -1.05578312e+04  -3.41197387e+00   7.5377e-01   0.86
     CG4     -1.05590756e+04  -1.24434221e+00   2.9390e+00   0.84
     CG5     -1.05612454e+04  -2.16980814e+00   2.1147e-01   0.69
     CG6     -1.05614108e+04  -1.65448765e-01   3.8187e-02   0.73
     CG7     -1.05614604e+04  -4.95572966e-02   9.4653e-03   0.82
     CG8     -1.05614545e+04   5.84834542e-03   4.0272e-03   0.83
     CG9     -1.05614660e+04  -1.14637521e-02   3.2802e-05   0.84
     CG10    -1.05614667e+04  -7.15660527e-04   1.0858e-04   1.49
     CG11    -1.05614670e+04  -2.76978455e-04   4.3701e-07   0.99
     CG12    -1.05614670e+04  -2.85083072e-06   1.8108e-06   1.27
     CG13    -1.05614670e+04  -9.63329225e-07   8.6578e-08   0.86
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       14.88  17       0.88   100.00 
     Driver           reading                     0.04   1        0.04   0.24   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 14.84  1        14.84  99.76  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.95   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.01   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.01   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.02   1        0.02   0.16   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    14.34  1        14.34  96.34  
     ESolver_KS_PW    runner                      14.33  1        14.33  96.33  
     ESolver_KS_PW    before_scf                  0.25   1        0.25   1.66   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.00   2        0.00   0.02   
     PW_Basis_Sup     recip2real                  0.02   99       0.00   0.12   
     PW_Basis_Sup     gathers_scatterp            0.01   99       0.00   0.05   
     Potential        init_pot                    0.01   1        0.01   0.04   
     Potential        update_from_charge          0.08   14       0.01   0.52   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.08   14       0.01   0.52   
     H_Hartree_pw     v_hartree                   0.01   14       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.02   125      0.00   0.14   
     PW_Basis_Sup     gatherp_scatters            0.01   125      0.00   0.06   
     PotXC            cal_v_eff                   0.07   14       0.01   0.47   
     XC_Functional    v_xc                        0.07   14       0.01   0.47   
     Potential        interpolate_vrs             0.00   14       0.00   0.00   
     Symmetry         rhog_symmetry               0.03   14       0.00   0.17   
     Symmetry         group fft grids             0.01   14       0.00   0.04   
     PSIInit          initialize_psi              0.24   1        0.24   1.58   
     Nonlocal         getvnl                      0.11   280      0.00   0.74   
     pp_cell_vnl      getvnl                      0.11   280      0.00   0.74   
     Structure_Factor get_sk                      0.01   280      0.00   0.08   
     DiagoIterAssist  diagH_subspace              2.09   260      0.01   14.07  
     Operator         hPsi                        10.60  39368    0.00   71.21  
     Operator         EkineticPW                  0.08   39368    0.00   0.56   
     Operator         VeffPW                      7.51   39368    0.00   50.45  
     PW_Basis_K       recip2real                  4.44   59908    0.00   29.83  
     PW_Basis_K       gathers_scatterp            2.27   59908    0.00   15.25  
     PW_Basis_K       real2recip                  2.94   49508    0.00   19.76  
     PW_Basis_K       gatherp_scatters            1.29   49508    0.00   8.70   
     Operator         NonlocalPW                  2.94   39368    0.00   19.73  
     Nonlocal         add_nonlocal_pp             1.23   39368    0.00   8.28   
     DiagoIterAssist  diagH_LAPACK                0.19   260      0.00   1.28   
     ESolver_KS_PW    hamilt2density_single       13.95  13       1.07   93.77  
     HSolverPW        solve                       13.92  13       1.07   93.54  
     DiagoCG          diag_once                   11.00  260      0.04   73.95  
     DiagoCG_New      spsi_func                   0.12   78216    0.00   0.82   
     DiagoCG_New      hpsi_func                   8.92   39108    0.00   59.92  
     ElecStatePW      psiToRho                    0.87   13       0.07   5.87   
     Charge_Mixing    get_drho                    0.01   13       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   13       0.00   0.00   
     Charge           mix_rho                     0.01   12       0.00   0.04   
     Charge           Broyden_mixing              0.00   12       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.00   120      0.00   0.01   
     ESolver_KS_PW    after_scf                   0.05   1        0.05   0.32   
     ModuleIO         write_rhog                  0.04   1        0.04   0.29   
     ModuleIO         write_istate_info           0.01   1        0.01   0.08   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:07:23 2025
     FINISH Time  : Sat Apr 26 12:07:38 2025
     TOTAL  Time  : 15
     SEE INFORMATION IN : OUT.Mn_ecutwfc_20/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:07:39 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_25/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 24 * 24 * 24
     UNIFORM GRID DIM(BIG)   : 24 * 24 * 24
     DONE(0.448112   SEC) : SETUP UNITCELL
     DONE(0.591758   SEC) : SYMMETRY
     DONE(0.88171    SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.883616   SEC) : INIT PLANEWAVE
     DONE(0.88598    SEC) : LOCAL POTENTIAL
     DONE(0.914735   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 5.13916
     DONE(0.914814   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.13929    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.06594364e+04   0.00000000e+00   6.8617e+00   3.30
     CG2     -1.06839817e+04  -2.45452233e+01   3.8404e+00   1.40
     CG3     -1.06903008e+04  -6.31909987e+00   1.7204e+01   0.93
     CG4     -1.06983980e+04  -8.09718980e+00   3.2933e+00   0.79
     CG5     -1.07012656e+04  -2.86764746e+00   1.3330e+00   0.78
     CG6     -1.07011198e+04   1.45832546e-01   2.4669e+00   0.80
     CG7     -1.07026765e+04  -1.55672460e+00   3.3756e-02   0.78
     CG8     -1.07030399e+04  -3.63384332e-01   3.1785e-02   1.33
     CG9     -1.07029623e+04   7.76033094e-02   1.8381e-02   0.77
     CG10    -1.07030065e+04  -4.42515787e-02   1.1540e-03   0.91
     CG11    -1.07030172e+04  -1.06355346e-02   5.3827e-03   1.46
     CG12    -1.07030148e+04   2.40869780e-03   1.0173e-03   0.99
     CG13    -1.07030177e+04  -2.90008509e-03   6.9165e-06   0.96
     CG14    -1.07030181e+04  -4.75414137e-04   1.4211e-05   1.56
     CG15    -1.07030181e+04   3.69056224e-06   1.0091e-05   0.77
     CG16    -1.07030181e+04   2.25142918e-05   7.2973e-06   0.77
     CG17    -1.07030181e+04   1.60011835e-05   3.2143e-06   0.79
     CG18    -1.07030181e+04  -1.97905994e-05   4.4873e-07   1.09
     CG19    -1.07030181e+04   1.64988689e-06   4.8524e-07   1.00
     CG20    -1.07030181e+04   1.51008348e-07   6.6700e-08   0.83
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       23.24  17       1.37   100.00 
     Driver           reading                     0.42   1        0.42   1.82   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 22.82  1        22.82  98.18  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.62   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.01   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.03   1        0.03   0.12   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    22.30  1        22.30  95.98  
     ESolver_KS_PW    runner                      22.30  1        22.30  95.98  
     ESolver_KS_PW    before_scf                  0.22   1        0.22   0.97   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.00   2        0.00   0.02   
     PW_Basis_Sup     recip2real                  0.05   148      0.00   0.23   
     PW_Basis_Sup     gathers_scatterp            0.03   148      0.00   0.14   
     Potential        init_pot                    0.01   1        0.01   0.03   
     Potential        update_from_charge          0.15   21       0.01   0.63   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.14   21       0.01   0.62   
     H_Hartree_pw     v_hartree                   0.03   21       0.00   0.13   
     PW_Basis_Sup     real2recip                  0.03   188      0.00   0.12   
     PW_Basis_Sup     gatherp_scatters            0.01   188      0.00   0.05   
     PotXC            cal_v_eff                   0.11   21       0.01   0.49   
     XC_Functional    v_xc                        0.11   21       0.01   0.48   
     Potential        interpolate_vrs             0.00   21       0.00   0.00   
     Symmetry         rhog_symmetry               0.05   21       0.00   0.21   
     Symmetry         group fft grids             0.01   21       0.00   0.05   
     PSIInit          initialize_psi              0.21   1        0.21   0.92   
     Nonlocal         getvnl                      0.24   420      0.00   1.04   
     pp_cell_vnl      getvnl                      0.24   420      0.00   1.03   
     Structure_Factor get_sk                      0.03   420      0.00   0.11   
     DiagoIterAssist  diagH_subspace              3.39   400      0.01   14.58  
     Operator         hPsi                        16.12  51298    0.00   69.38  
     Operator         EkineticPW                  0.13   51298    0.00   0.58   
     Operator         VeffPW                      10.56  51298    0.00   45.44  
     PW_Basis_K       recip2real                  6.54   82898    0.00   28.14  
     PW_Basis_K       gathers_scatterp            3.36   82898    0.00   14.48  
     PW_Basis_K       real2recip                  4.08   66898    0.00   17.56  
     PW_Basis_K       gatherp_scatters            1.81   66898    0.00   7.79   
     Operator         NonlocalPW                  5.34   51298    0.00   22.97  
     Nonlocal         add_nonlocal_pp             2.45   51298    0.00   10.56  
     DiagoIterAssist  diagH_LAPACK                0.27   400      0.00   1.15   
     ESolver_KS_PW    hamilt2density_single       21.86  20       1.09   94.08  
     HSolverPW        solve                       21.80  20       1.09   93.81  
     DiagoCG          diag_once                   16.79  400      0.04   72.27  
     DiagoCG_New      spsi_func                   0.20   101796   0.00   0.88   
     DiagoCG_New      hpsi_func                   13.40  50898    0.00   57.66  
     ElecStatePW      psiToRho                    1.44   20       0.07   6.22   
     Charge_Mixing    get_drho                    0.01   20       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   20       0.00   0.00   
     Charge           mix_rho                     0.01   19       0.00   0.06   
     Charge           Broyden_mixing              0.01   19       0.00   0.03   
     Charge_Mixing    inner_product_recip_hartree 0.01   232      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.05   1        0.05   0.22   
     ModuleIO         write_rhog                  0.05   1        0.05   0.20   
     ModuleIO         write_istate_info           0.01   1        0.01   0.05   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:07:39 2025
     FINISH Time  : Sat Apr 26 12:08:03 2025
     TOTAL  Time  : 24
     SEE INFORMATION IN : OUT.Mn_ecutwfc_25/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:08:04 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_30/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 25 * 25 * 25
     UNIFORM GRID DIM(BIG)   : 25 * 25 * 25
     DONE(0.0648677  SEC) : SETUP UNITCELL
     DONE(0.20572    SEC) : SYMMETRY
     DONE(0.521883   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.524191   SEC) : INIT PLANEWAVE
     DONE(0.527623   SEC) : LOCAL POTENTIAL
     DONE(0.56017    SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 6.92139
     DONE(0.560234   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.857947   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.07502020e+04   0.00000000e+00   3.1671e+00   4.68
     CG2     -1.07674214e+04  -1.72194420e+01   3.3153e+00   1.89
     CG3     -1.07694165e+04  -1.99505037e+00   1.6351e+01   1.17
     CG4     -1.07782624e+04  -8.84596360e+00   2.8406e-01   1.02
     CG5     -1.07793179e+04  -1.05542369e+00   1.5375e+00   1.55
     CG6     -1.07801203e+04  -8.02399874e-01   1.8817e-02   1.04
     CG7     -1.07800515e+04   6.87644146e-02   3.4600e-01   1.69
     CG8     -1.07801726e+04  -1.21122271e-01   9.4139e-03   1.12
     CG9     -1.07801881e+04  -1.55184073e-02   7.8499e-04   1.16
     CG10    -1.07801998e+04  -1.16940677e-02   8.1571e-04   1.88
     CG11    -1.07801999e+04  -9.24367302e-05   1.0876e-04   1.04
     CG12    -1.07802010e+04  -1.01748334e-03   1.8674e-05   1.54
     CG13    -1.07802009e+04   2.72463131e-06   1.9868e-05   1.27
     CG14    -1.07802009e+04   4.60917612e-05   1.2006e-05   1.08
     CG15    -1.07802009e+04  -1.75193288e-05   2.5182e-06   1.10
     CG16    -1.07802009e+04  -4.60546924e-06   1.3873e-06   1.32
     CG17    -1.07802009e+04  -3.52829327e-06   8.5503e-08   1.22
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       26.71  17       1.57   100.00 
     Driver           reading                     0.04   1        0.04   0.13   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 26.67  1        26.67  99.87  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.01   1        0.01   0.02   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.53   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.01   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.01   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.03   1        0.03   0.12   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    26.13  1        26.13  97.82  
     ESolver_KS_PW    runner                      26.13  1        26.13  97.82  
     ESolver_KS_PW    before_scf                  0.30   1        0.30   1.11   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.03   
     PW_Basis_Sup     recip2real                  0.04   127      0.00   0.15   
     PW_Basis_Sup     gathers_scatterp            0.01   127      0.00   0.05   
     Potential        init_pot                    0.01   1        0.01   0.03   
     Potential        update_from_charge          0.13   18       0.01   0.48   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.13   18       0.01   0.48   
     H_Hartree_pw     v_hartree                   0.01   18       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.03   161      0.00   0.12   
     PW_Basis_Sup     gatherp_scatters            0.01   161      0.00   0.04   
     PotXC            cal_v_eff                   0.11   18       0.01   0.43   
     XC_Functional    v_xc                        0.11   18       0.01   0.43   
     Potential        interpolate_vrs             0.00   18       0.00   0.00   
     Symmetry         rhog_symmetry               0.05   18       0.00   0.19   
     Symmetry         group fft grids             0.01   18       0.00   0.04   
     PSIInit          initialize_psi              0.28   1        0.28   1.05   
     Nonlocal         getvnl                      0.26   360      0.00   0.96   
     pp_cell_vnl      getvnl                      0.26   360      0.00   0.96   
     Structure_Factor get_sk                      0.03   360      0.00   0.11   
     DiagoIterAssist  diagH_subspace              3.78   340      0.01   14.15  
     Operator         hPsi                        19.35  45762    0.00   72.43  
     Operator         EkineticPW                  0.14   45762    0.00   0.54   
     Operator         VeffPW                      12.52  45762    0.00   46.86  
     PW_Basis_K       recip2real                  7.72   72622    0.00   28.91  
     PW_Basis_K       gathers_scatterp            3.45   72622    0.00   12.91  
     PW_Basis_K       real2recip                  4.99   59022    0.00   18.69  
     PW_Basis_K       gatherp_scatters            1.86   59022    0.00   6.96   
     Operator         NonlocalPW                  6.60   45762    0.00   24.72  
     Nonlocal         add_nonlocal_pp             3.26   45762    0.00   12.19  
     DiagoIterAssist  diagH_LAPACK                0.22   340      0.00   0.83   
     ESolver_KS_PW    hamilt2density_single       25.63  17       1.51   95.97  
     HSolverPW        solve                       25.56  17       1.50   95.71  
     DiagoCG          diag_once                   20.06  340      0.06   75.12  
     DiagoCG_New      spsi_func                   0.21   90844    0.00   0.80   
     DiagoCG_New      hpsi_func                   16.23  45422    0.00   60.77  
     ElecStatePW      psiToRho                    1.58   17       0.09   5.92   
     Charge_Mixing    get_drho                    0.01   17       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   17       0.00   0.01   
     Charge           mix_rho                     0.02   16       0.00   0.06   
     Charge           Broyden_mixing              0.01   16       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.01   184      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.05   1        0.05   0.18   
     ModuleIO         write_rhog                  0.04   1        0.04   0.15   
     ModuleIO         write_istate_info           0.01   1        0.01   0.05   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:08:04 2025
     FINISH Time  : Sat Apr 26 12:08:30 2025
     TOTAL  Time  : 26
     SEE INFORMATION IN : OUT.Mn_ecutwfc_30/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:08:31 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_35/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 27 * 27 * 27
     UNIFORM GRID DIM(BIG)   : 27 * 27 * 27
     DONE(0.0655738  SEC) : SETUP UNITCELL
     DONE(0.201316   SEC) : SYMMETRY
     DONE(0.486474   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.489193   SEC) : INIT PLANEWAVE
     DONE(0.4925     SEC) : LOCAL POTENTIAL
     DONE(0.525526   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 8.52051
     DONE(0.525602   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.891573   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.07892336e+04   0.00000000e+00   1.6438e+00   5.78
     CG2     -1.07941935e+04  -4.95991412e+00   4.5557e+00   2.43
     CG3     -1.07995082e+04  -5.31471055e+00   2.1096e+01   1.53
     CG4     -1.08086703e+04  -9.16205400e+00   1.2008e+00   1.45
     CG5     -1.08089485e+04  -2.78247471e-01   1.2666e+00   1.26
     CG6     -1.08097708e+04  -8.22305775e-01   2.2955e-01   1.24
     CG7     -1.08102203e+04  -4.49482280e-01   4.1872e-01   1.55
     CG8     -1.08103755e+04  -1.55191025e-01   1.1426e-01   1.27
     CG9     -1.08103971e+04  -2.16267944e-02   1.8406e-02   1.24
     CG10    -1.08104089e+04  -1.17410370e-02   7.6158e-03   1.55
     CG11    -1.08104093e+04  -3.87997466e-04   2.1712e-03   1.37
     CG12    -1.08104181e+04  -8.81500933e-03   2.6205e-03   1.69
     CG13    -1.08104184e+04  -3.22656840e-04   2.9589e-04   1.22
     CG14    -1.08104190e+04  -6.12969968e-04   4.5771e-06   1.61
     CG15    -1.08104192e+04  -1.84269702e-04   1.0928e-05   2.41
     CG16    -1.08104192e+04  -1.13067754e-05   9.7409e-06   1.37
     CG17    -1.08104192e+04   6.12147692e-06   5.8015e-06   1.25
     CG18    -1.08104192e+04   3.48741164e-06   1.9846e-06   1.30
     CG19    -1.08104192e+04  -1.85347523e-06   5.3687e-08   1.26
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       33.76  17       1.99   100.00 
     Driver           reading                     0.04   1        0.04   0.12   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 33.72  1        33.72  99.88  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.40   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.03   1        0.03   0.09   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    33.22  1        33.22  98.38  
     ESolver_KS_PW    runner                      33.22  1        33.22  98.38  
     ESolver_KS_PW    before_scf                  0.37   1        0.37   1.08   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.02   
     PW_Basis_Sup     recip2real                  0.04   141      0.00   0.13   
     PW_Basis_Sup     gathers_scatterp            0.01   141      0.00   0.04   
     Potential        init_pot                    0.01   1        0.01   0.03   
     Potential        update_from_charge          0.17   20       0.01   0.51   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.17   20       0.01   0.50   
     H_Hartree_pw     v_hartree                   0.01   20       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.04   179      0.00   0.13   
     PW_Basis_Sup     gatherp_scatters            0.01   179      0.00   0.04   
     PotXC            cal_v_eff                   0.15   20       0.01   0.46   
     XC_Functional    v_xc                        0.15   20       0.01   0.45   
     Potential        interpolate_vrs             0.00   20       0.00   0.00   
     Symmetry         rhog_symmetry               0.07   20       0.00   0.21   
     Symmetry         group fft grids             0.02   20       0.00   0.05   
     PSIInit          initialize_psi              0.34   1        0.34   1.02   
     Nonlocal         getvnl                      0.36   400      0.00   1.06   
     pp_cell_vnl      getvnl                      0.36   400      0.00   1.06   
     Structure_Factor get_sk                      0.04   400      0.00   0.13   
     DiagoIterAssist  diagH_subspace              5.02   380      0.01   14.87  
     Operator         hPsi                        24.59  47272    0.00   72.84  
     Operator         EkineticPW                  0.18   47272    0.00   0.53   
     Operator         VeffPW                      15.85  47272    0.00   46.95  
     PW_Basis_K       recip2real                  9.87   77292    0.00   29.25  
     PW_Basis_K       gathers_scatterp            4.20   77292    0.00   12.43  
     PW_Basis_K       real2recip                  6.32   62092    0.00   18.71  
     PW_Basis_K       gatherp_scatters            2.29   62092    0.00   6.77   
     Operator         NonlocalPW                  8.47   47272    0.00   25.10  
     Nonlocal         add_nonlocal_pp             4.25   47272    0.00   12.58  
     DiagoIterAssist  diagH_LAPACK                0.25   380      0.00   0.74   
     ESolver_KS_PW    hamilt2density_single       32.59  19       1.72   96.53  
     HSolverPW        solve                       32.50  19       1.71   96.27  
     DiagoCG          diag_once                   25.08  380      0.07   74.28  
     DiagoCG_New      spsi_func                   0.26   93784    0.00   0.76   
     DiagoCG_New      hpsi_func                   20.41  46892    0.00   60.45  
     ElecStatePW      psiToRho                    2.17   19       0.11   6.43   
     Charge_Mixing    get_drho                    0.01   19       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   19       0.00   0.00   
     Charge           mix_rho                     0.02   18       0.00   0.06   
     Charge           Broyden_mixing              0.01   18       0.00   0.03   
     Charge_Mixing    inner_product_recip_hartree 0.01   216      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.06   1        0.06   0.18   
     ModuleIO         write_rhog                  0.05   1        0.05   0.16   
     ModuleIO         write_istate_info           0.01   1        0.01   0.04   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:08:31 2025
     FINISH Time  : Sat Apr 26 12:09:05 2025
     TOTAL  Time  : 34
     SEE INFORMATION IN : OUT.Mn_ecutwfc_35/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:09:06 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_40/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.136637   SEC) : SETUP UNITCELL
     DONE(0.333503   SEC) : SYMMETRY
     DONE(0.631188   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.63464    SEC) : INIT PLANEWAVE
     DONE(0.638754   SEC) : LOCAL POTENTIAL
     DONE(0.674418   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 10.4248
     DONE(0.674486   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.13691    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08068065e+04   0.00000000e+00   9.7374e-01   7.59
     CG2     -1.08097254e+04  -2.91897719e+00   4.6280e+00   3.54
     CG3     -1.08122443e+04  -2.51881836e+00   2.7777e+01   2.18
     CG4     -1.08234439e+04  -1.11996121e+01   1.1958e+00   2.05
     CG5     -1.08238148e+04  -3.70890234e-01   3.5184e-01   1.68
     CG6     -1.08239322e+04  -1.17458734e-01   1.8937e-01   1.64
     CG7     -1.08240611e+04  -1.28848335e-01   3.3918e-02   1.73
     CG8     -1.08242257e+04  -1.64619485e-01   1.6621e-03   2.22
     CG9     -1.08242378e+04  -1.21270852e-02   2.4657e-03   2.80
     CG10    -1.08242396e+04  -1.74852078e-03   3.1368e-03   1.67
     CG11    -1.08242402e+04  -5.86413208e-04   3.8480e-04   1.62
     CG12    -1.08242412e+04  -1.06344939e-03   1.8293e-04   2.05
     CG13    -1.08242413e+04  -1.30878919e-04   2.7394e-05   1.72
     CG14    -1.08242416e+04  -2.58969370e-04   2.7429e-05   2.11
     CG15    -1.08242413e+04   2.65389703e-04   2.1714e-04   1.68
     CG16    -1.08242414e+04  -1.03730007e-04   2.7426e-05   1.67
     CG17    -1.08242413e+04   9.80462710e-05   4.8389e-05   1.68
     CG18    -1.08242413e+04   6.76726573e-05   1.1867e-05   1.70
     CG19    -1.08242414e+04  -1.14403046e-04   1.1212e-05   2.49
     CG20    -1.08242414e+04   3.10985612e-05   9.0838e-06   1.62
     CG21    -1.08242414e+04   1.06050530e-05   3.4627e-06   1.67
     CG22    -1.08242414e+04  -5.48585192e-06   6.5919e-07   1.82
     CG23    -1.08242414e+04  -6.72677714e-06   3.0860e-07   2.62
     CG24    -1.08242414e+04   4.61297271e-07   1.6659e-07   1.68
     CG25    -1.08242414e+04   5.59320370e-07   6.2492e-08   1.87
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       56.47  17       3.32   100.00 
     Driver           reading                     0.04   1        0.04   0.07   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 56.43  1        56.43  99.93  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.20   1        0.20   0.35   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.03   1        0.03   0.06   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    55.77  1        55.77  98.76  
     ESolver_KS_PW    runner                      55.77  1        55.77  98.76  
     ESolver_KS_PW    before_scf                  0.46   1        0.46   0.82   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.02   
     PW_Basis_Sup     recip2real                  0.07   183      0.00   0.13   
     PW_Basis_Sup     gathers_scatterp            0.02   183      0.00   0.04   
     Potential        init_pot                    0.01   1        0.01   0.02   
     Potential        update_from_charge          0.31   26       0.01   0.54   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.30   26       0.01   0.54   
     H_Hartree_pw     v_hartree                   0.02   26       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.09   233      0.00   0.15   
     PW_Basis_Sup     gatherp_scatters            0.03   233      0.00   0.06   
     PotXC            cal_v_eff                   0.28   26       0.01   0.49   
     XC_Functional    v_xc                        0.28   26       0.01   0.49   
     Potential        interpolate_vrs             0.00   26       0.00   0.00   
     Symmetry         rhog_symmetry               0.11   26       0.00   0.19   
     Symmetry         group fft grids             0.02   26       0.00   0.04   
     PSIInit          initialize_psi              0.44   1        0.44   0.78   
     Nonlocal         getvnl                      0.52   520      0.00   0.92   
     pp_cell_vnl      getvnl                      0.52   520      0.00   0.92   
     Structure_Factor get_sk                      0.06   520      0.00   0.11   
     DiagoIterAssist  diagH_subspace              8.89   500      0.02   15.74  
     Operator         hPsi                        41.68  60259    0.00   73.82  
     Operator         EkineticPW                  0.28   60259    0.00   0.49   
     Operator         VeffPW                      28.28  60259    0.00   50.08  
     PW_Basis_K       recip2real                  17.30  99759    0.00   30.64  
     PW_Basis_K       gathers_scatterp            6.77   99759    0.00   12.00  
     PW_Basis_K       real2recip                  11.87  79759    0.00   21.03  
     PW_Basis_K       gatherp_scatters            3.91   79759    0.00   6.93   
     Operator         NonlocalPW                  13.01  60259    0.00   23.05  
     Nonlocal         add_nonlocal_pp             6.51   60259    0.00   11.53  
     DiagoIterAssist  diagH_LAPACK                0.33   500      0.00   0.59   
     ESolver_KS_PW    hamilt2density_single       54.74  25       2.19   96.94  
     HSolverPW        solve                       54.60  25       2.18   96.69  
     DiagoCG          diag_once                   41.20  500      0.08   72.97  
     DiagoCG_New      spsi_func                   0.35   119518   0.00   0.62   
     DiagoCG_New      hpsi_func                   34.02  59759    0.00   60.25  
     ElecStatePW      psiToRho                    4.08   25       0.16   7.22   
     Charge_Mixing    get_drho                    0.02   25       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   25       0.00   0.00   
     Charge           mix_rho                     0.03   24       0.00   0.06   
     Charge           Broyden_mixing              0.01   24       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.01   312      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.21   1        0.21   0.37   
     ModuleIO         write_rhog                  0.20   1        0.20   0.36   
     ModuleIO         write_istate_info           0.01   1        0.01   0.02   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:09:06 2025
     FINISH Time  : Sat Apr 26 12:10:03 2025
     TOTAL  Time  : 57
     SEE INFORMATION IN : OUT.Mn_ecutwfc_40/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:10:04 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_45/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 32 * 32 * 32
     UNIFORM GRID DIM(BIG)   : 32 * 32 * 32
     DONE(0.0672756  SEC) : SETUP UNITCELL
     DONE(0.207891   SEC) : SYMMETRY
     DONE(0.509588   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.513316   SEC) : INIT PLANEWAVE
     DONE(0.517921   SEC) : LOCAL POTENTIAL
     DONE(0.555761   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 12.2925
     DONE(0.555832   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.04247    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08151820e+04   0.00000000e+00   6.6876e-01   8.11
     CG2     -1.08228612e+04  -7.67916848e+00   2.2659e+00   3.75
     CG3     -1.08240118e+04  -1.15066495e+00   1.4230e+01   2.37
     CG4     -1.08289752e+04  -4.96332008e+00   2.3038e+00   1.85
     CG5     -1.08297008e+04  -7.25651078e-01   9.5993e-02   1.68
     CG6     -1.08298644e+04  -1.63579273e-01   1.7444e-02   2.13
     CG7     -1.08298824e+04  -1.80621956e-02   1.1206e-01   2.26
     CG8     -1.08299174e+04  -3.49966142e-02   1.0273e-02   1.69
     CG9     -1.08299179e+04  -4.17274762e-04   4.0103e-03   1.74
     CG10    -1.08299266e+04  -8.78743375e-03   5.6724e-04   1.94
     CG11    -1.08299299e+04  -3.26452544e-03   9.2006e-04   2.64
     CG12    -1.08299288e+04   1.10351663e-03   8.9530e-04   1.70
     CG13    -1.08299288e+04  -9.24040481e-06   1.4698e-04   1.66
     CG14    -1.08299290e+04  -2.30861113e-04   1.3323e-05   1.98
     CG15    -1.08299289e+04   1.54083255e-04   5.3854e-04   3.26
     CG16    -1.08299292e+04  -2.64373019e-04   9.9404e-06   2.71
     CG17    -1.08299291e+04   1.55399388e-05   7.1433e-05   1.68
     CG18    -1.08299292e+04  -1.78534182e-05   3.9029e-06   1.74
     CG19    -1.08299291e+04   1.69433232e-05   2.9750e-06   1.69
     CG20    -1.08299292e+04  -7.35199380e-06   4.2808e-07   1.93
     CG21    -1.08299292e+04  -2.01978287e-06   5.5175e-07   2.40
     CG22    -1.08299292e+04  -4.28723446e-08   2.1520e-08   1.69
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       53.77  17       3.16   100.00 
     Driver           reading                     0.04   1        0.04   0.07   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 53.73  1        53.73  99.93  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.26   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.07   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    53.19  1        53.19  98.93  
     ESolver_KS_PW    runner                      53.19  1        53.19  98.93  
     ESolver_KS_PW    before_scf                  0.49   1        0.49   0.90   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.02   
     PW_Basis_Sup     recip2real                  0.07   162      0.00   0.12   
     PW_Basis_Sup     gathers_scatterp            0.03   162      0.00   0.05   
     Potential        init_pot                    0.02   1        0.02   0.03   
     Potential        update_from_charge          0.30   23       0.01   0.56   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.30   23       0.01   0.55   
     H_Hartree_pw     v_hartree                   0.03   23       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.07   206      0.00   0.13   
     PW_Basis_Sup     gatherp_scatters            0.03   206      0.00   0.05   
     PotXC            cal_v_eff                   0.27   23       0.01   0.50   
     XC_Functional    v_xc                        0.27   23       0.01   0.50   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.11   23       0.00   0.21   
     Symmetry         group fft grids             0.02   23       0.00   0.04   
     PSIInit          initialize_psi              0.46   1        0.46   0.85   
     Nonlocal         getvnl                      0.51   460      0.00   0.94   
     pp_cell_vnl      getvnl                      0.51   460      0.00   0.94   
     Structure_Factor get_sk                      0.05   460      0.00   0.10   
     DiagoIterAssist  diagH_subspace              7.76   440      0.02   14.43  
     Operator         hPsi                        39.01  56056    0.00   72.54  
     Operator         EkineticPW                  0.30   56056    0.00   0.55   
     Operator         VeffPW                      24.35  56056    0.00   45.29  
     PW_Basis_K       recip2real                  14.73  90816    0.00   27.39  
     PW_Basis_K       gathers_scatterp            7.01   90816    0.00   13.03  
     PW_Basis_K       real2recip                  9.64   73216    0.00   17.93  
     PW_Basis_K       gatherp_scatters            3.65   73216    0.00   6.79   
     Operator         NonlocalPW                  14.25  56056    0.00   26.49  
     Nonlocal         add_nonlocal_pp             7.10   56056    0.00   13.20  
     DiagoIterAssist  diagH_LAPACK                0.29   440      0.00   0.55   
     ESolver_KS_PW    hamilt2density_single       52.29  22       2.38   97.25  
     HSolverPW        solve                       52.15  22       2.37   96.99  
     DiagoCG          diag_once                   40.41  440      0.09   75.16  
     DiagoCG_New      spsi_func                   0.37   111232   0.00   0.68   
     DiagoCG_New      hpsi_func                   32.47  55616    0.00   60.39  
     ElecStatePW      psiToRho                    3.57   22       0.16   6.64   
     Charge_Mixing    get_drho                    0.02   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.03   21       0.00   0.06   
     Charge           Broyden_mixing              0.01   21       0.00   0.03   
     Charge_Mixing    inner_product_recip_hartree 0.01   264      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.12   
     ModuleIO         write_rhog                  0.06   1        0.06   0.11   
     ModuleIO         write_istate_info           0.01   1        0.01   0.02   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:10:04 2025
     FINISH Time  : Sat Apr 26 12:10:58 2025
     TOTAL  Time  : 54
     SEE INFORMATION IN : OUT.Mn_ecutwfc_45/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:10:59 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_50/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0632191  SEC) : SETUP UNITCELL
     DONE(0.20267    SEC) : SYMMETRY
     DONE(0.486115   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.49049    SEC) : INIT PLANEWAVE
     DONE(0.496542   SEC) : LOCAL POTENTIAL
     DONE(0.536446   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 14.4165
     DONE(0.536541   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.27589    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08157476e+04   0.00000000e+00   6.7885e-01  11.23
     CG2     -1.08149586e+04   7.88909452e-01   7.5619e+00   5.22
     CG3     -1.08218266e+04  -6.86793715e+00   2.0346e+01   3.46
     CG4     -1.08289218e+04  -7.09519741e+00   7.9937e+00   3.20
     CG5     -1.08311010e+04  -2.17922452e+00   3.8796e-01   2.43
     CG6     -1.08311742e+04  -7.31837846e-02   1.4672e-01   2.45
     CG7     -1.08311779e+04  -3.75774277e-03   1.2032e-01   2.53
     CG8     -1.08311998e+04  -2.18247422e-02   1.8490e-02   2.40
     CG9     -1.08313150e+04  -1.15226507e-01   1.8871e-03   3.53
     CG10    -1.08313150e+04  -4.65669579e-05   2.7241e-03   3.36
     CG11    -1.08313156e+04  -5.35370978e-04   2.7847e-03   2.57
     CG12    -1.08313170e+04  -1.46934151e-03   1.3868e-04   2.41
     CG13    -1.08313183e+04  -1.25245937e-03   6.8606e-06   3.68
     CG14    -1.08313184e+04  -7.12976496e-05   6.2364e-05   4.28
     CG15    -1.08313184e+04   1.66361449e-05   1.1601e-05   2.40
     CG16    -1.08313183e+04   4.41833719e-06   4.5760e-06   2.43
     CG17    -1.08313183e+04   2.22096626e-05   2.8416e-05   2.46
     CG18    -1.08313184e+04  -2.97027467e-05   1.1313e-06   3.28
     CG19    -1.08313184e+04   1.51305022e-06   6.6884e-07   2.72
     CG20    -1.08313184e+04  -2.02346268e-07   8.2446e-08   2.74
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       70.14  17       4.13   100.00 
     Driver           reading                     0.04   1        0.04   0.05   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 70.10  1        70.10  99.95  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.20   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.05   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    69.58  1        69.58  99.21  
     ESolver_KS_PW    runner                      69.58  1        69.58  99.21  
     ESolver_KS_PW    before_scf                  0.74   1        0.74   1.05   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.01   
     PW_Basis_Sup     recip2real                  0.08   148      0.00   0.11   
     PW_Basis_Sup     gathers_scatterp            0.03   148      0.00   0.04   
     Potential        init_pot                    0.02   1        0.02   0.03   
     Potential        update_from_charge          0.38   21       0.02   0.54   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.38   21       0.02   0.54   
     H_Hartree_pw     v_hartree                   0.03   21       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.09   188      0.00   0.14   
     PW_Basis_Sup     gatherp_scatters            0.03   188      0.00   0.05   
     PotXC            cal_v_eff                   0.34   21       0.02   0.49   
     XC_Functional    v_xc                        0.34   21       0.02   0.49   
     Potential        interpolate_vrs             0.00   21       0.00   0.00   
     Symmetry         rhog_symmetry               0.12   21       0.01   0.17   
     Symmetry         group fft grids             0.02   21       0.00   0.04   
     PSIInit          initialize_psi              0.71   1        0.71   1.01   
     Nonlocal         getvnl                      0.55   420      0.00   0.78   
     pp_cell_vnl      getvnl                      0.55   420      0.00   0.78   
     Structure_Factor get_sk                      0.06   420      0.00   0.08   
     DiagoIterAssist  diagH_subspace              10.83  400      0.03   15.44  
     Operator         hPsi                        52.89  52341    0.00   75.40  
     Operator         EkineticPW                  0.31   52341    0.00   0.44   
     Operator         VeffPW                      37.25  52341    0.00   53.11  
     PW_Basis_K       recip2real                  21.66  83941    0.00   30.88  
     PW_Basis_K       gathers_scatterp            8.35   83941    0.00   11.90  
     PW_Basis_K       real2recip                  15.73  67941    0.00   22.43  
     PW_Basis_K       gatherp_scatters            4.28   67941    0.00   6.10   
     Operator         NonlocalPW                  15.21  52341    0.00   21.69  
     Nonlocal         add_nonlocal_pp             7.55   52341    0.00   10.76  
     DiagoIterAssist  diagH_LAPACK                0.26   400      0.00   0.38   
     ESolver_KS_PW    hamilt2density_single       68.35  20       3.42   97.45  
     HSolverPW        solve                       68.19  20       3.41   97.22  
     DiagoCG          diag_once                   52.02  400      0.13   74.16  
     DiagoCG_New      spsi_func                   0.38   103882   0.00   0.54   
     DiagoCG_New      hpsi_func                   43.33  51941    0.00   61.77  
     ElecStatePW      psiToRho                    5.10   20       0.25   7.27   
     Charge_Mixing    get_drho                    0.02   20       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   20       0.00   0.00   
     Charge           mix_rho                     0.04   19       0.00   0.05   
     Charge           Broyden_mixing              0.01   19       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.01   232      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.06   1        0.06   0.09   
     ModuleIO         write_rhog                  0.05   1        0.05   0.07   
     ModuleIO         write_istate_info           0.01   1        0.01   0.02   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:10:59 2025
     FINISH Time  : Sat Apr 26 12:12:09 2025
     TOTAL  Time  : 70
     SEE INFORMATION IN : OUT.Mn_ecutwfc_50/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:12:10 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_55/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0635428  SEC) : SETUP UNITCELL
     DONE(0.202813   SEC) : SYMMETRY
     DONE(0.491457   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.496378   SEC) : INIT PLANEWAVE
     DONE(0.502429   SEC) : LOCAL POTENTIAL
     DONE(0.548731   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 16.748
     DONE(0.548806   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.26712    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08161072e+04   0.00000000e+00   6.6264e-01  12.24
     CG2     -1.08166134e+04  -5.06257657e-01   6.9640e+00   5.74
     CG3     -1.08209709e+04  -4.35745595e+00   2.2034e+01   3.26
     CG4     -1.08305523e+04  -9.58141328e+00   1.5980e+00   3.24
     CG5     -1.08311956e+04  -6.43326566e-01   1.0637e-01   2.59
     CG6     -1.08315087e+04  -3.13087223e-01   3.6913e-02   3.07
     CG7     -1.08315647e+04  -5.59426236e-02   4.5612e-02   2.96
     CG8     -1.08315838e+04  -1.91637878e-02   3.0851e-03   2.59
     CG9     -1.08316042e+04  -2.03358392e-02   1.2267e-03   3.84
     CG10    -1.08316049e+04  -7.65966953e-04   1.1274e-03   2.78
     CG11    -1.08316055e+04  -5.22270943e-04   7.7906e-05   2.61
     CG12    -1.08316060e+04  -5.13270298e-04   3.4267e-05   4.43
     CG13    -1.08316059e+04   6.95766630e-05   1.2824e-05   2.76
     CG14    -1.08316060e+04  -5.58816086e-05   6.1781e-05   4.31
     CG15    -1.08316060e+04  -4.11336164e-05   6.8099e-06   2.90
     CG16    -1.08316060e+04   3.78595889e-06   2.5251e-06   2.65
     CG17    -1.08316060e+04  -9.08021472e-06   4.2354e-08   3.40
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       66.71  17       3.92   100.00 
     Driver           reading                     0.03   1        0.03   0.05   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 66.67  1        66.67  99.95  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.21   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.07   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    66.13  1        66.13  99.13  
     ESolver_KS_PW    runner                      66.13  1        66.13  99.13  
     ESolver_KS_PW    before_scf                  0.72   1        0.72   1.08   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.01   
     PW_Basis_Sup     recip2real                  0.07   127      0.00   0.11   
     PW_Basis_Sup     gathers_scatterp            0.03   127      0.00   0.04   
     Potential        init_pot                    0.02   1        0.02   0.03   
     Potential        update_from_charge          0.34   18       0.02   0.50   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.33   18       0.02   0.50   
     H_Hartree_pw     v_hartree                   0.03   18       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.08   161      0.00   0.12   
     PW_Basis_Sup     gatherp_scatters            0.03   161      0.00   0.04   
     PotXC            cal_v_eff                   0.30   18       0.02   0.45   
     XC_Functional    v_xc                        0.30   18       0.02   0.45   
     Potential        interpolate_vrs             0.00   18       0.00   0.00   
     Symmetry         rhog_symmetry               0.11   18       0.01   0.17   
     Symmetry         group fft grids             0.02   18       0.00   0.04   
     PSIInit          initialize_psi              0.68   1        0.68   1.02   
     Nonlocal         getvnl                      0.53   360      0.00   0.80   
     pp_cell_vnl      getvnl                      0.53   360      0.00   0.79   
     Structure_Factor get_sk                      0.06   360      0.00   0.08   
     DiagoIterAssist  diagH_subspace              9.51   340      0.03   14.25  
     Operator         hPsi                        49.69  47033    0.00   74.49  
     Operator         EkineticPW                  0.31   47033    0.00   0.46   
     Operator         VeffPW                      33.62  47033    0.00   50.40  
     PW_Basis_K       recip2real                  19.46  73893    0.00   29.17  
     PW_Basis_K       gathers_scatterp            7.68   73893    0.00   11.51  
     PW_Basis_K       real2recip                  14.25  60293    0.00   21.36  
     PW_Basis_K       gatherp_scatters            4.00   60293    0.00   6.00   
     Operator         NonlocalPW                  15.66  47033    0.00   23.48  
     Nonlocal         add_nonlocal_pp             7.87   47033    0.00   11.80  
     DiagoIterAssist  diagH_LAPACK                0.24   340      0.00   0.35   
     ESolver_KS_PW    hamilt2density_single       64.97  17       3.82   97.40  
     HSolverPW        solve                       64.83  17       3.81   97.18  
     DiagoCG          diag_once                   50.64  340      0.15   75.91  
     DiagoCG_New      spsi_func                   0.36   93386    0.00   0.54   
     DiagoCG_New      hpsi_func                   41.37  46693    0.00   62.02  
     ElecStatePW      psiToRho                    4.43   17       0.26   6.65   
     Charge_Mixing    get_drho                    0.02   17       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   17       0.00   0.00   
     Charge           mix_rho                     0.03   16       0.00   0.05   
     Charge           Broyden_mixing              0.01   16       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.01   184      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.06   1        0.06   0.09   
     ModuleIO         write_rhog                  0.05   1        0.05   0.08   
     ModuleIO         write_istate_info           0.02   1        0.02   0.03   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:12:10 2025
     FINISH Time  : Sat Apr 26 12:13:16 2025
     TOTAL  Time  : 66
     SEE INFORMATION IN : OUT.Mn_ecutwfc_55/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:13:17 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_60/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.070016   SEC) : SETUP UNITCELL
     DONE(0.214252   SEC) : SYMMETRY
     DONE(0.513415   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.518693   SEC) : INIT PLANEWAVE
     DONE(0.525178   SEC) : LOCAL POTENTIAL
     DONE(0.569255   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 18.8965
     DONE(0.569354   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.35589    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08166785e+04   0.00000000e+00   6.0630e-01  13.01
     CG2     -1.08159777e+04   7.00830325e-01   7.7004e+00   6.16
     CG3     -1.08206922e+04  -4.71449072e+00   2.3884e+01   3.95
     CG4     -1.08280928e+04  -7.40060143e+00   9.6692e+00   3.44
     CG5     -1.08313193e+04  -3.22649654e+00   2.7254e-01   2.76
     CG6     -1.08314258e+04  -1.06535961e-01   1.7255e-01   2.87
     CG7     -1.08314935e+04  -6.76405591e-02   8.2798e-02   2.77
     CG8     -1.08315608e+04  -6.72771590e-02   2.2931e-02   2.98
     CG9     -1.08316272e+04  -6.64664538e-02   5.4567e-03   3.56
     CG10    -1.08316381e+04  -1.08386549e-02   8.4480e-04   3.31
     CG11    -1.08316415e+04  -3.46294906e-03   9.9482e-04   3.92
     CG12    -1.08316420e+04  -4.77918392e-04   6.9938e-04   2.72
     CG13    -1.08316420e+04  -3.73843226e-05   5.1718e-04   2.68
     CG14    -1.08316413e+04   7.46570830e-04   7.0880e-04   2.76
     CG15    -1.08316400e+04   1.30612625e-03   1.4709e-03   2.80
     CG16    -1.08316399e+04   7.64970022e-05   8.6115e-04   2.81
     CG17    -1.08316404e+04  -4.89203449e-04   1.0133e-03   2.70
     CG18    -1.08316408e+04  -4.09953680e-04   1.3983e-06   2.79
     CG19    -1.08316410e+04  -1.37618882e-04   7.7258e-06   6.16
     CG20    -1.08316410e+04  -1.38520334e-05   6.6774e-07   4.40
     CG21    -1.08316410e+04  -4.32789954e-07   2.6954e-07   2.78
     CG22    -1.08316410e+04   4.50386223e-07   6.7292e-08   2.82
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       85.69  17       5.04   100.00 
     Driver           reading                     0.04   1        0.04   0.05   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 85.65  1        85.65  99.95  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.17   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.05   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    85.10  1        85.10  99.31  
     ESolver_KS_PW    runner                      85.10  1        85.10  99.31  
     ESolver_KS_PW    before_scf                  0.79   1        0.79   0.92   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.11   162      0.00   0.12   
     PW_Basis_Sup     gathers_scatterp            0.04   162      0.00   0.05   
     Potential        init_pot                    0.02   1        0.02   0.02   
     Potential        update_from_charge          0.45   23       0.02   0.52   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.45   23       0.02   0.52   
     H_Hartree_pw     v_hartree                   0.04   23       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.11   206      0.00   0.13   
     PW_Basis_Sup     gatherp_scatters            0.04   206      0.00   0.05   
     PotXC            cal_v_eff                   0.40   23       0.02   0.47   
     XC_Functional    v_xc                        0.40   23       0.02   0.47   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.16   23       0.01   0.18   
     Symmetry         group fft grids             0.03   23       0.00   0.04   
     PSIInit          initialize_psi              0.75   1        0.75   0.87   
     Nonlocal         getvnl                      0.75   460      0.00   0.87   
     pp_cell_vnl      getvnl                      0.75   460      0.00   0.87   
     Structure_Factor get_sk                      0.08   460      0.00   0.09   
     DiagoIterAssist  diagH_subspace              12.83  440      0.03   14.97  
     Operator         hPsi                        63.01  55811    0.00   73.53  
     Operator         EkineticPW                  0.41   55811    0.00   0.48   
     Operator         VeffPW                      40.93  55811    0.00   47.77  
     PW_Basis_K       recip2real                  24.04  90571    0.00   28.05  
     PW_Basis_K       gathers_scatterp            9.50   90571    0.00   11.09  
     PW_Basis_K       real2recip                  17.33  72971    0.00   20.23  
     PW_Basis_K       gatherp_scatters            4.91   72971    0.00   5.73   
     Operator         NonlocalPW                  21.54  55811    0.00   25.14  
     Nonlocal         add_nonlocal_pp             10.81  55811    0.00   12.62  
     DiagoIterAssist  diagH_LAPACK                0.30   440      0.00   0.35   
     ESolver_KS_PW    hamilt2density_single       83.62  22       3.80   97.58  
     HSolverPW        solve                       83.41  22       3.79   97.34  
     DiagoCG          diag_once                   64.23  440      0.15   74.96  
     DiagoCG_New      spsi_func                   0.49   110742   0.00   0.57   
     DiagoCG_New      hpsi_func                   51.90  55371    0.00   60.57  
     ElecStatePW      psiToRho                    5.79   22       0.26   6.76   
     Charge_Mixing    get_drho                    0.03   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.05   21       0.00   0.06   
     Charge           Broyden_mixing              0.02   21       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.02   264      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.18   1        0.18   0.21   
     ModuleIO         write_rhog                  0.17   1        0.17   0.20   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:13:17 2025
     FINISH Time  : Sat Apr 26 12:14:43 2025
     TOTAL  Time  : 86
     SEE INFORMATION IN : OUT.Mn_ecutwfc_60/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:14:44 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_65/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 40 * 40 * 40
     UNIFORM GRID DIM(BIG)   : 40 * 40 * 40
     DONE(0.0689529  SEC) : SETUP UNITCELL
     DONE(0.213073   SEC) : SYMMETRY
     DONE(0.503531   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.509991   SEC) : INIT PLANEWAVE
     DONE(0.517507   SEC) : LOCAL POTENTIAL
     DONE(0.565441   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 21.2402
     DONE(0.565515   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.57038    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08160038e+04   0.00000000e+00   6.4928e-01  16.11
     CG2     -1.08124401e+04   3.56367373e+00   1.0927e+01   7.68
     CG3     -1.08213069e+04  -8.86678879e+00   2.1346e+01   4.95
     CG4     -1.08276402e+04  -6.33325266e+00   1.1407e+01   4.48
     CG5     -1.08312438e+04  -3.60357025e+00   4.6341e-01   3.65
     CG6     -1.08313714e+04  -1.27645021e-01   2.1453e-01   3.44
     CG7     -1.08314534e+04  -8.19949934e-02   1.3293e-01   3.63
     CG8     -1.08315020e+04  -4.85910695e-02   8.7745e-02   3.56
     CG9     -1.08315388e+04  -3.68398666e-02   4.0377e-02   3.60
     CG10    -1.08316182e+04  -7.94089725e-02   2.4113e-03   4.03
     CG11    -1.08316442e+04  -2.59445368e-02   1.4536e-03   5.53
     CG12    -1.08316434e+04   7.65881060e-04   1.6453e-03   3.48
     CG13    -1.08316435e+04  -9.23979615e-05   1.0165e-03   3.43
     CG14    -1.08316404e+04   3.12441922e-03   7.2486e-04   3.49
     CG15    -1.08316414e+04  -1.06538296e-03   3.1125e-04   3.85
     CG16    -1.08316435e+04  -2.01232820e-03   9.0607e-04   4.99
     CG17    -1.08316430e+04   4.59786513e-04   8.2659e-04   3.45
     CG18    -1.08316427e+04   2.60693742e-04   4.4701e-04   3.44
     CG19    -1.08316421e+04   6.29854107e-04   3.4201e-04   3.44
     CG20    -1.08316411e+04   9.80067820e-04   1.8037e-04   3.56
     CG21    -1.08316419e+04  -8.08296649e-04   1.1727e-05   4.38
     CG22    -1.08316420e+04  -6.31871399e-05   8.6273e-06   4.97
     CG23    -1.08316420e+04   9.12645133e-06   9.4391e-06   3.44
     CG24    -1.08316420e+04   3.51802383e-05   5.5400e-06   3.49
     CG25    -1.08316420e+04  -2.20139103e-05   2.1369e-07   4.24
     CG26    -1.08316420e+04  -1.35783616e-06   6.9001e-07   6.27
     CG27    -1.08316420e+04   1.87277916e-06   3.2755e-07   3.55
     CG28    -1.08316420e+04  -6.84497346e-07   5.5638e-08   4.80
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       130.59 17       7.68   100.00 
     Driver           reading                     0.04   1        0.04   0.03   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 130.55 1        130.55 99.97  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.11   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    130.00 1        130.00 99.55  
     ESolver_KS_PW    runner                      130.00 1        130.00 99.55  
     ESolver_KS_PW    before_scf                  1.00   1        1.00   0.77   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.18   204      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.07   204      0.00   0.05   
     Potential        init_pot                    0.03   1        0.03   0.02   
     Potential        update_from_charge          0.76   29       0.03   0.58   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.76   29       0.03   0.58   
     H_Hartree_pw     v_hartree                   0.06   29       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.21   260      0.00   0.16   
     PW_Basis_Sup     gatherp_scatters            0.09   260      0.00   0.07   
     PotXC            cal_v_eff                   0.69   29       0.02   0.53   
     XC_Functional    v_xc                        0.69   29       0.02   0.53   
     Potential        interpolate_vrs             0.00   29       0.00   0.00   
     Symmetry         rhog_symmetry               0.24   29       0.01   0.18   
     Symmetry         group fft grids             0.05   29       0.00   0.04   
     PSIInit          initialize_psi              0.96   1        0.96   0.73   
     Nonlocal         getvnl                      1.12   580      0.00   0.86   
     pp_cell_vnl      getvnl                      1.12   580      0.00   0.86   
     Structure_Factor get_sk                      0.12   580      0.00   0.09   
     DiagoIterAssist  diagH_subspace              21.47  560      0.04   16.44  
     Operator         hPsi                        97.96  67119    0.00   75.02  
     Operator         EkineticPW                  0.58   67119    0.00   0.44   
     Operator         VeffPW                      68.23  67119    0.00   52.25  
     PW_Basis_K       recip2real                  40.66  111359   0.00   31.14  
     PW_Basis_K       gathers_scatterp            15.10  111359   0.00   11.56  
     PW_Basis_K       real2recip                  28.67  88959    0.00   21.96  
     PW_Basis_K       gatherp_scatters            7.24   88959    0.00   5.54   
     Operator         NonlocalPW                  29.00  67119    0.00   22.21  
     Nonlocal         add_nonlocal_pp             14.43  67119    0.00   11.05  
     DiagoIterAssist  diagH_LAPACK                0.39   560      0.00   0.30   
     ESolver_KS_PW    hamilt2density_single       128.05 28       4.57   98.06  
     HSolverPW        solve                       127.73 28       4.56   97.81  
     DiagoCG          diag_once                   95.12  560      0.17   72.84  
     DiagoCG_New      spsi_func                   0.66   133118   0.00   0.51   
     DiagoCG_New      hpsi_func                   78.96  66559    0.00   60.46  
     ElecStatePW      psiToRho                    10.21  28       0.36   7.82   
     Charge_Mixing    get_drho                    0.05   28       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   28       0.00   0.00   
     Charge           mix_rho                     0.08   27       0.00   0.06   
     Charge           Broyden_mixing              0.03   27       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   360      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.05   
     ModuleIO         write_rhog                  0.06   1        0.06   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:14:44 2025
     FINISH Time  : Sat Apr 26 12:16:55 2025
     TOTAL  Time  : 131
     SEE INFORMATION IN : OUT.Mn_ecutwfc_65/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:16:56 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_70/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 40 * 40 * 40
     UNIFORM GRID DIM(BIG)   : 40 * 40 * 40
     DONE(0.0634506  SEC) : SETUP UNITCELL
     DONE(0.20419    SEC) : SYMMETRY
     DONE(0.497928   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.504653   SEC) : INIT PLANEWAVE
     DONE(0.512632   SEC) : LOCAL POTENTIAL
     DONE(0.562356   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 23.6328
     DONE(0.562448   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.59193    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08159352e+04   0.00000000e+00   6.6690e-01  16.95
     CG2     -1.08192031e+04  -3.26783345e+00   4.7049e+00   8.13
     CG3     -1.08217707e+04  -2.56764059e+00   2.2248e+01   4.93
     CG4     -1.08301404e+04  -8.36970323e+00   3.4075e+00   4.51
     CG5     -1.08312665e+04  -1.12614276e+00   2.2126e-01   3.65
     CG6     -1.08314037e+04  -1.37208118e-01   1.9082e-01   3.76
     CG7     -1.08314910e+04  -8.72235688e-02   4.9773e-02   3.58
     CG8     -1.08316147e+04  -1.23720025e-01   2.3308e-02   4.41
     CG9     -1.08316298e+04  -1.51475370e-02   9.0768e-03   3.98
     CG10    -1.08316380e+04  -8.15841368e-03   5.0075e-03   4.13
     CG11    -1.08316411e+04  -3.12211165e-03   1.4325e-03   3.85
     CG12    -1.08316445e+04  -3.34596651e-03   4.6537e-04   4.52
     CG13    -1.08316457e+04  -1.22428724e-03   2.0480e-04   4.39
     CG14    -1.08316449e+04   8.18623639e-04   4.3016e-04   3.93
     CG15    -1.08316452e+04  -3.52422334e-04   6.1792e-04   4.81
     CG16    -1.08316454e+04  -2.15322174e-04   2.2036e-05   3.69
     CG17    -1.08316455e+04  -3.02298259e-05   1.7182e-05   4.64
     CG18    -1.08316455e+04   2.70273497e-06   3.5671e-06   3.63
     CG19    -1.08316455e+04  -8.54842255e-06   1.4899e-06   4.69
     CG20    -1.08316455e+04  -1.72545990e-06   1.5673e-07   4.14
     CG21    -1.08316455e+04  -2.04055160e-06   2.0140e-07   6.51
     CG22    -1.08316455e+04   4.51953120e-07   5.1397e-08   3.70
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       112.21 17       6.60   100.00 
     Driver           reading                     0.03   1        0.03   0.03   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 112.17 1        112.17 99.97  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.13   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.04   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    111.62 1        111.62 99.48  
     ESolver_KS_PW    runner                      111.62 1        111.62 99.48  
     ESolver_KS_PW    before_scf                  1.03   1        1.03   0.92   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.14   162      0.00   0.13   
     PW_Basis_Sup     gathers_scatterp            0.05   162      0.00   0.05   
     Potential        init_pot                    0.03   1        0.03   0.03   
     Potential        update_from_charge          0.62   23       0.03   0.56   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.62   23       0.03   0.55   
     H_Hartree_pw     v_hartree                   0.06   23       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.17   206      0.00   0.15   
     PW_Basis_Sup     gatherp_scatters            0.06   206      0.00   0.05   
     PotXC            cal_v_eff                   0.56   23       0.02   0.50   
     XC_Functional    v_xc                        0.56   23       0.02   0.49   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.20   23       0.01   0.18   
     Symmetry         group fft grids             0.04   23       0.00   0.04   
     PSIInit          initialize_psi              0.97   1        0.97   0.86   
     Nonlocal         getvnl                      1.02   460      0.00   0.91   
     pp_cell_vnl      getvnl                      1.02   460      0.00   0.90   
     Structure_Factor get_sk                      0.12   460      0.00   0.10   
     DiagoIterAssist  diagH_subspace              17.34  440      0.04   15.45  
     Operator         hPsi                        83.38  55588    0.00   74.30  
     Operator         EkineticPW                  0.50   55588    0.00   0.45   
     Operator         VeffPW                      56.58  55588    0.00   50.43  
     PW_Basis_K       recip2real                  33.35  90348    0.00   29.72  
     PW_Basis_K       gathers_scatterp            12.44  90348    0.00   11.09  
     PW_Basis_K       real2recip                  23.89  72748    0.00   21.29  
     PW_Basis_K       gatherp_scatters            6.15   72748    0.00   5.48   
     Operator         NonlocalPW                  26.16  55588    0.00   23.31  
     Nonlocal         add_nonlocal_pp             12.88  55588    0.00   11.48  
     DiagoIterAssist  diagH_LAPACK                0.29   440      0.00   0.26   
     ESolver_KS_PW    hamilt2density_single       109.81 22       4.99   97.86  
     HSolverPW        solve                       109.55 22       4.98   97.63  
     DiagoCG          diag_once                   83.37  440      0.19   74.30  
     DiagoCG_New      spsi_func                   0.59   110296   0.00   0.52   
     DiagoCG_New      hpsi_func                   68.16  55148    0.00   60.74  
     ElecStatePW      psiToRho                    8.09   22       0.37   7.21   
     Charge_Mixing    get_drho                    0.04   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.07   21       0.00   0.06   
     Charge           Broyden_mixing              0.03   21       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.02   264      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.06   
     ModuleIO         write_rhog                  0.06   1        0.06   0.05   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:16:56 2025
     FINISH Time  : Sat Apr 26 12:18:48 2025
     TOTAL  Time  : 112
     SEE INFORMATION IN : OUT.Mn_ecutwfc_70/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:18:49 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_75/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 45 * 45 * 45
     UNIFORM GRID DIM(BIG)   : 45 * 45 * 45
     DONE(0.0794071  SEC) : SETUP UNITCELL
     DONE(0.220888   SEC) : SYMMETRY
     DONE(0.513463   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.521252   SEC) : INIT PLANEWAVE
     DONE(0.530494   SEC) : LOCAL POTENTIAL
     DONE(0.580085   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 26.3672
     DONE(0.580178   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.06153    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08161501e+04   0.00000000e+00   6.1983e-01  23.19
     CG2     -1.08212264e+04  -5.07633326e+00   3.7564e+00  11.41
     CG3     -1.08223629e+04  -1.13645671e+00   2.1321e+01   7.06
     CG4     -1.08302617e+04  -7.89885879e+00   3.4283e+00   6.43
     CG5     -1.08313675e+04  -1.10581684e+00   1.7074e-01   5.14
     CG6     -1.08314913e+04  -1.23732973e-01   1.5459e-01   5.51
     CG7     -1.08315458e+04  -5.45538711e-02   3.4248e-02   5.15
     CG8     -1.08316201e+04  -7.42304319e-02   1.3452e-02   6.42
     CG9     -1.08316367e+04  -1.66020086e-02   3.6787e-03   5.82
     CG10    -1.08316456e+04  -8.97038185e-03   1.2436e-03   6.34
     CG11    -1.08316473e+04  -1.71174779e-03   7.3435e-04   6.00
     CG12    -1.08316477e+04  -3.93099549e-04   3.2505e-04   5.35
     CG13    -1.08316482e+04  -5.18423158e-04   1.5030e-04   5.79
     CG14    -1.08316476e+04   6.61561776e-04   3.4055e-04   5.77
     CG15    -1.08316479e+04  -2.81335733e-04   2.7688e-04   6.00
     CG16    -1.08316478e+04   2.24557706e-05   3.0887e-04   5.14
     CG17    -1.08316480e+04  -1.21180855e-04   1.7426e-05   5.12
     CG18    -1.08316480e+04  -1.72248651e-05   8.0604e-06   5.84
     CG19    -1.08316480e+04   6.02508723e-06   6.0746e-06   5.72
     CG20    -1.08316480e+04  -1.63745999e-06   1.1783e-06   5.35
     CG21    -1.08316480e+04  -6.73454202e-06   3.3433e-07   7.40
     CG22    -1.08316480e+04  -6.06390699e-07   1.6786e-07   6.35
     CG23    -1.08316480e+04   5.37060533e-08   1.8125e-08   5.38
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       159.84 17       9.40   100.00 
     Driver           reading                     0.04   1        0.04   0.03   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 159.80 1        159.80 99.97  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.09   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    159.24 1        159.24 99.62  
     ESolver_KS_PW    runner                      159.24 1        159.24 99.62  
     ESolver_KS_PW    before_scf                  1.48   1        1.48   0.93   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.23   169      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.08   169      0.00   0.05   
     Potential        init_pot                    0.04   1        0.04   0.03   
     Potential        update_from_charge          0.93   24       0.04   0.58   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.93   24       0.04   0.58   
     H_Hartree_pw     v_hartree                   0.09   24       0.00   0.06   
     PW_Basis_Sup     real2recip                  0.27   215      0.00   0.17   
     PW_Basis_Sup     gatherp_scatters            0.09   215      0.00   0.05   
     PotXC            cal_v_eff                   0.83   24       0.03   0.52   
     XC_Functional    v_xc                        0.83   24       0.03   0.52   
     Potential        interpolate_vrs             0.00   24       0.00   0.00   
     Symmetry         rhog_symmetry               0.23   24       0.01   0.15   
     Symmetry         group fft grids             0.05   24       0.00   0.03   
     PSIInit          initialize_psi              1.41   1        1.41   0.88   
     Nonlocal         getvnl                      1.17   480      0.00   0.73   
     pp_cell_vnl      getvnl                      1.17   480      0.00   0.73   
     Structure_Factor get_sk                      0.14   480      0.00   0.09   
     DiagoIterAssist  diagH_subspace              27.64  460      0.06   17.29  
     Operator         hPsi                        123.08 55604    0.00   77.00  
     Operator         EkineticPW                  0.56   55604    0.00   0.35   
     Operator         VeffPW                      93.02  55604    0.00   58.19  
     PW_Basis_K       recip2real                  56.34  91944    0.00   35.25  
     PW_Basis_K       gathers_scatterp            17.26  91944    0.00   10.80  
     PW_Basis_K       real2recip                  38.83  73544    0.00   24.29  
     PW_Basis_K       gatherp_scatters            8.14   73544    0.00   5.09   
     Operator         NonlocalPW                  29.36  55604    0.00   18.37  
     Nonlocal         add_nonlocal_pp             14.48  55604    0.00   9.06   
     DiagoIterAssist  diagH_LAPACK                0.31   460      0.00   0.20   
     ESolver_KS_PW    hamilt2density_single       156.62 23       6.81   97.98  
     HSolverPW        solve                       156.29 23       6.80   97.78  
     DiagoCG          diag_once                   114.47 460      0.25   71.61  
     DiagoCG_New      spsi_func                   0.65   110288   0.00   0.41   
     DiagoCG_New      hpsi_func                   97.80  55144    0.00   61.19  
     ElecStatePW      psiToRho                    13.60  23       0.59   8.51   
     Charge_Mixing    get_drho                    0.06   23       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   23       0.00   0.00   
     Charge           mix_rho                     0.09   22       0.00   0.06   
     Charge           Broyden_mixing              0.03   22       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   280      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.08   1        0.08   0.05   
     ModuleIO         write_rhog                  0.06   1        0.06   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:18:49 2025
     FINISH Time  : Sat Apr 26 12:21:29 2025
     TOTAL  Time  : 160
     SEE INFORMATION IN : OUT.Mn_ecutwfc_75/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:21:30 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_80/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 45 * 45 * 45
     UNIFORM GRID DIM(BIG)   : 45 * 45 * 45
     DONE(0.0678048  SEC) : SETUP UNITCELL
     DONE(0.205417   SEC) : SYMMETRY
     DONE(0.513803   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.522039   SEC) : INIT PLANEWAVE
     DONE(0.531996   SEC) : LOCAL POTENTIAL
     DONE(0.585297   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 28.5889
     DONE(0.585375   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.09613    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08156600e+04   0.00000000e+00   6.8073e-01  24.64
     CG2     -1.08244736e+04  -8.81365898e+00   2.4677e+00  11.76
     CG3     -1.08271101e+04  -2.63649362e+00   1.1131e+01   7.25
     CG4     -1.08312512e+04  -4.14108841e+00   5.0503e-01   5.80
     CG5     -1.08313730e+04  -1.21750159e-01   1.4536e-01   5.48
     CG6     -1.08315156e+04  -1.42617640e-01   2.4711e-02   6.02
     CG7     -1.08315866e+04  -7.10566819e-02   1.0649e-01   7.32
     CG8     -1.08316350e+04  -4.83585116e-02   1.0206e-02   5.35
     CG9     -1.08316364e+04  -1.36035862e-03   8.9872e-03   5.58
     CG10    -1.08316394e+04  -3.06498071e-03   9.7461e-04   5.42
     CG11    -1.08316502e+04  -1.08027165e-02   7.6838e-04   8.44
     CG12    -1.08316487e+04   1.51465712e-03   7.4487e-04   5.45
     CG13    -1.08316472e+04   1.54874684e-03   3.7030e-04   5.46
     CG14    -1.08316485e+04  -1.36633950e-03   9.2744e-05   7.61
     CG15    -1.08316483e+04   2.12127927e-04   1.0592e-04   5.94
     CG16    -1.08316484e+04  -1.35323950e-04   1.6145e-05   6.01
     CG17    -1.08316485e+04  -4.77739249e-05   2.1682e-06   6.88
     CG18    -1.08316485e+04  -4.14815959e-06   1.8831e-06   7.43
     CG19    -1.08316485e+04   3.12465719e-06   1.3091e-06   5.38
     CG20    -1.08316485e+04   4.68955577e-06   9.8771e-07   5.39
     CG21    -1.08316485e+04  -3.88540813e-06   6.7608e-07   7.85
     CG22    -1.08316485e+04  -5.68723299e-08   1.2452e-07   5.39
     CG23    -1.08316485e+04  -4.36602788e-07   7.6189e-09   7.64
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       171.70 17       10.10  100.00 
     Driver           reading                     0.04   1        0.04   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 171.66 1        171.66 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.08   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    171.09 1        171.09 99.64  
     ESolver_KS_PW    runner                      171.09 1        171.09 99.64  
     ESolver_KS_PW    before_scf                  1.51   1        1.51   0.88   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.24   169      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.08   169      0.00   0.05   
     Potential        init_pot                    0.04   1        0.04   0.03   
     Potential        update_from_charge          0.93   24       0.04   0.54   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.92   24       0.04   0.54   
     H_Hartree_pw     v_hartree                   0.08   24       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.26   215      0.00   0.15   
     PW_Basis_Sup     gatherp_scatters            0.08   215      0.00   0.04   
     PotXC            cal_v_eff                   0.83   24       0.03   0.49   
     XC_Functional    v_xc                        0.83   24       0.03   0.48   
     Potential        interpolate_vrs             0.00   24       0.00   0.00   
     Symmetry         rhog_symmetry               0.26   24       0.01   0.15   
     Symmetry         group fft grids             0.05   24       0.00   0.03   
     PSIInit          initialize_psi              1.44   1        1.44   0.84   
     Nonlocal         getvnl                      1.25   480      0.00   0.73   
     pp_cell_vnl      getvnl                      1.25   480      0.00   0.73   
     Structure_Factor get_sk                      0.15   480      0.00   0.09   
     DiagoIterAssist  diagH_subspace              28.48  460      0.06   16.59  
     Operator         hPsi                        131.57 57228    0.00   76.63  
     Operator         EkineticPW                  0.67   57228    0.00   0.39   
     Operator         VeffPW                      97.42  57228    0.00   56.74  
     PW_Basis_K       recip2real                  58.72  93568    0.00   34.20  
     PW_Basis_K       gathers_scatterp            18.12  93568    0.00   10.56  
     PW_Basis_K       real2recip                  40.74  75168    0.00   23.73  
     PW_Basis_K       gatherp_scatters            8.71   75168    0.00   5.07   
     Operator         NonlocalPW                  33.34  57228    0.00   19.42  
     Nonlocal         add_nonlocal_pp             16.55  57228    0.00   9.64   
     DiagoIterAssist  diagH_LAPACK                0.31   460      0.00   0.18   
     ESolver_KS_PW    hamilt2density_single       168.44 23       7.32   98.10  
     HSolverPW        solve                       168.08 23       7.31   97.89  
     DiagoCG          diag_once                   125.16 460      0.27   72.90  
     DiagoCG_New      spsi_func                   0.73   113536   0.00   0.42   
     DiagoCG_New      hpsi_func                   105.67 56768    0.00   61.54  
     ElecStatePW      psiToRho                    13.76  23       0.60   8.02   
     Charge_Mixing    get_drho                    0.06   23       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.01   23       0.00   0.00   
     Charge           mix_rho                     0.09   22       0.00   0.05   
     Charge           Broyden_mixing              0.03   22       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   280      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.04   
     ModuleIO         write_rhog                  0.06   1        0.06   0.03   
     ModuleIO         write_istate_info           0.02   1        0.02   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:21:30 2025
     FINISH Time  : Sat Apr 26 12:24:22 2025
     TOTAL  Time  : 172
     SEE INFORMATION IN : OUT.Mn_ecutwfc_80/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:24:23 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_85/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 45 * 45 * 45
     UNIFORM GRID DIM(BIG)   : 45 * 45 * 45
     DONE(0.0929905  SEC) : SETUP UNITCELL
     DONE(0.233617   SEC) : SYMMETRY
     DONE(0.524704   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.533586   SEC) : INIT PLANEWAVE
     DONE(0.543904   SEC) : LOCAL POTENTIAL
     DONE(0.599605   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 31.9702
     DONE(0.599676   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.21097    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08164594e+04   0.00000000e+00   6.1226e-01  25.43
     CG2     -1.08236136e+04  -7.15414391e+00   2.7477e+00  12.88
     CG3     -1.08235045e+04   1.09105222e-01   1.7825e+01   7.05
     CG4     -1.08309778e+04  -7.47334241e+00   5.2038e-01   6.95
     CG5     -1.08312329e+04  -2.55098942e-01   1.2642e-01   5.60
     CG6     -1.08315763e+04  -3.43371411e-01   1.8092e-02   6.65
     CG7     -1.08316041e+04  -2.78309754e-02   9.4908e-02   7.65
     CG8     -1.08316499e+04  -4.57540360e-02   6.0481e-03   5.58
     CG9     -1.08316527e+04  -2.78969106e-03   3.4406e-02   5.87
     CG10    -1.08316493e+04   3.35956727e-03   1.0067e-02   5.58
     CG11    -1.08316479e+04   1.41901917e-03   3.1897e-03   5.59
     CG12    -1.08316436e+04   4.24517065e-03   2.2724e-03   5.78
     CG13    -1.08316440e+04  -3.11641099e-04   6.1977e-04   5.89
     CG14    -1.08316481e+04  -4.15806280e-03   1.1477e-03   9.60
     CG15    -1.08316484e+04  -2.69014536e-04   4.6851e-04   5.89
     CG16    -1.08316484e+04  -1.01373845e-05   1.7949e-04   5.57
     CG17    -1.08316483e+04   1.36279344e-04   1.8227e-04   6.75
     CG18    -1.08316487e+04  -4.00200896e-04   3.8040e-05   6.18
     CG19    -1.08316488e+04  -1.04806652e-04   3.9439e-06   6.64
     CG20    -1.08316488e+04   8.97740122e-06   8.1091e-06   7.96
     CG21    -1.08316488e+04  -6.99184414e-06   1.4956e-06   6.16
     CG22    -1.08316488e+04  -2.47576383e-06   3.3156e-07   6.70
     CG23    -1.08316488e+04  -1.85188977e-06   1.5686e-08   8.17
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       178.45 17       10.50  100.00 
     Driver           reading                     0.06   1        0.06   0.03   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 178.39 1        178.39 99.97  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.08   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    177.83 1        177.83 99.65  
     ESolver_KS_PW    runner                      177.83 1        177.83 99.65  
     ESolver_KS_PW    before_scf                  1.61   1        1.61   0.90   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.02   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.25   169      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.08   169      0.00   0.05   
     Potential        init_pot                    0.04   1        0.04   0.02   
     Potential        update_from_charge          0.96   24       0.04   0.54   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.95   24       0.04   0.53   
     H_Hartree_pw     v_hartree                   0.10   24       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.29   215      0.00   0.16   
     PW_Basis_Sup     gatherp_scatters            0.09   215      0.00   0.05   
     PotXC            cal_v_eff                   0.85   24       0.04   0.48   
     XC_Functional    v_xc                        0.85   24       0.04   0.47   
     Potential        interpolate_vrs             0.00   24       0.00   0.00   
     Symmetry         rhog_symmetry               0.28   24       0.01   0.16   
     Symmetry         group fft grids             0.06   24       0.00   0.03   
     PSIInit          initialize_psi              1.54   1        1.54   0.86   
     Nonlocal         getvnl                      1.37   480      0.00   0.77   
     pp_cell_vnl      getvnl                      1.37   480      0.00   0.77   
     Structure_Factor get_sk                      0.15   480      0.00   0.08   
     DiagoIterAssist  diagH_subspace              29.31  460      0.06   16.42  
     Operator         hPsi                        135.80 57300    0.00   76.10  
     Operator         EkineticPW                  0.69   57300    0.00   0.39   
     Operator         VeffPW                      98.44  57300    0.00   55.17  
     PW_Basis_K       recip2real                  59.43  93640    0.00   33.30  
     PW_Basis_K       gathers_scatterp            18.26  93640    0.00   10.23  
     PW_Basis_K       real2recip                  41.31  75240    0.00   23.15  
     PW_Basis_K       gatherp_scatters            8.96   75240    0.00   5.02   
     Operator         NonlocalPW                  36.52  57300    0.00   20.47  
     Nonlocal         add_nonlocal_pp             18.21  57300    0.00   10.21  
     DiagoIterAssist  diagH_LAPACK                0.30   460      0.00   0.17   
     ESolver_KS_PW    hamilt2density_single       175.02 23       7.61   98.08  
     HSolverPW        solve                       174.63 23       7.59   97.86  
     DiagoCG          diag_once                   130.59 460      0.28   73.18  
     DiagoCG_New      spsi_func                   0.77   113680   0.00   0.43   
     DiagoCG_New      hpsi_func                   109.27 56840    0.00   61.23  
     ElecStatePW      psiToRho                    13.94  23       0.61   7.81   
     Charge_Mixing    get_drho                    0.07   23       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   23       0.00   0.00   
     Charge           mix_rho                     0.11   22       0.01   0.06   
     Charge           Broyden_mixing              0.04   22       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   280      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.08   1        0.08   0.05   
     ModuleIO         write_rhog                  0.07   1        0.07   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:24:23 2025
     FINISH Time  : Sat Apr 26 12:27:21 2025
     TOTAL  Time  : 178
     SEE INFORMATION IN : OUT.Mn_ecutwfc_85/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:27:22 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_90/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 45 * 45 * 45
     UNIFORM GRID DIM(BIG)   : 45 * 45 * 45
     DONE(0.0677761  SEC) : SETUP UNITCELL
     DONE(0.208854   SEC) : SYMMETRY
     DONE(0.487851   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.497887   SEC) : INIT PLANEWAVE
     DONE(0.508632   SEC) : LOCAL POTENTIAL
     DONE(0.56406    SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 34.5581
     DONE(0.564145   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.22324    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08162841e+04   0.00000000e+00   6.2673e-01  26.91
     CG2     -1.08239598e+04  -7.67561755e+00   2.5800e+00  13.31
     CG3     -1.08237607e+04   1.99041510e-01   1.7339e+01   7.31
     CG4     -1.08309139e+04  -7.15319068e+00   5.7505e-01   7.12
     CG5     -1.08311795e+04  -2.65600200e-01   1.1299e-01   5.84
     CG6     -1.08315875e+04  -4.08033351e-01   1.4063e-02   7.26
     CG7     -1.08316366e+04  -4.90170640e-02   6.1865e-03   8.20
     CG8     -1.08316217e+04   1.48317772e-02   1.1307e-01   6.34
     CG9     -1.08316572e+04  -3.54964184e-02   3.4478e-02   5.85
     CG10    -1.08316530e+04   4.22545819e-03   1.1391e-02   5.76
     CG11    -1.08316508e+04   2.16441001e-03   3.1589e-03   5.91
     CG12    -1.08316469e+04   3.90841079e-03   2.0268e-03   5.83
     CG13    -1.08316440e+04   2.96901349e-03   9.3873e-04   5.91
     CG14    -1.08316498e+04  -5.88466455e-03   6.6121e-04   9.06
     CG15    -1.08316500e+04  -1.28087121e-04   3.7717e-04   5.89
     CG16    -1.08316497e+04   2.70375469e-04   1.8300e-04   5.93
     CG17    -1.08316494e+04   3.06736521e-04   7.6019e-05   6.65
     CG18    -1.08316496e+04  -2.56877250e-04   5.1646e-04   8.53
     CG19    -1.08316500e+04  -3.51839544e-04   1.0089e-05   6.83
     CG20    -1.08316500e+04   2.93472235e-07   7.1867e-06   7.57
     CG21    -1.08316500e+04   5.93335183e-06   2.2992e-06   6.08
     CG22    -1.08316500e+04  -9.45494905e-06   3.8713e-08   7.93
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       178.36 17       10.49  100.00 
     Driver           reading                     0.03   1        0.03   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 178.32 1        178.32 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.08   
     PW_Basis_K       setuptransform              0.01   1        0.01   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.01   1        0.01   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    177.77 1        177.77 99.67  
     ESolver_KS_PW    runner                      177.77 1        177.77 99.67  
     ESolver_KS_PW    before_scf                  1.66   1        1.66   0.93   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.02   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.24   162      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.08   162      0.00   0.05   
     Potential        init_pot                    0.04   1        0.04   0.02   
     Potential        update_from_charge          0.94   23       0.04   0.53   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.94   23       0.04   0.53   
     H_Hartree_pw     v_hartree                   0.10   23       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.30   206      0.00   0.17   
     PW_Basis_Sup     gatherp_scatters            0.09   206      0.00   0.05   
     PotXC            cal_v_eff                   0.84   23       0.04   0.47   
     XC_Functional    v_xc                        0.83   23       0.04   0.47   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.28   23       0.01   0.16   
     Symmetry         group fft grids             0.06   23       0.00   0.03   
     PSIInit          initialize_psi              1.59   1        1.59   0.89   
     Nonlocal         getvnl                      1.42   460      0.00   0.80   
     pp_cell_vnl      getvnl                      1.42   460      0.00   0.80   
     Structure_Factor get_sk                      0.14   460      0.00   0.08   
     DiagoIterAssist  diagH_subspace              28.78  440      0.07   16.13  
     Operator         hPsi                        135.11 54805    0.00   75.75  
     Operator         EkineticPW                  0.74   54805    0.00   0.42   
     Operator         VeffPW                      95.39  54805    0.00   53.48  
     PW_Basis_K       recip2real                  57.54  89565    0.00   32.26  
     PW_Basis_K       gathers_scatterp            17.75  89565    0.00   9.95   
     PW_Basis_K       real2recip                  39.95  71965    0.00   22.40  
     PW_Basis_K       gatherp_scatters            8.71   71965    0.00   4.88   
     Operator         NonlocalPW                  38.84  54805    0.00   21.78  
     Nonlocal         add_nonlocal_pp             19.35  54805    0.00   10.85  
     DiagoIterAssist  diagH_LAPACK                0.30   440      0.00   0.17   
     ESolver_KS_PW    hamilt2density_single       174.93 22       7.95   98.08  
     HSolverPW        solve                       174.54 22       7.93   97.86  
     DiagoCG          diag_once                   131.51 440      0.30   73.73  
     DiagoCG_New      spsi_func                   0.80   108730   0.00   0.45   
     DiagoCG_New      hpsi_func                   109.26 54365    0.00   61.26  
     ElecStatePW      psiToRho                    13.35  22       0.61   7.49   
     Charge_Mixing    get_drho                    0.07   22       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.10   21       0.00   0.06   
     Charge           Broyden_mixing              0.04   21       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   264      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.08   1        0.08   0.05   
     ModuleIO         write_rhog                  0.07   1        0.07   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:27:22 2025
     FINISH Time  : Sat Apr 26 12:30:21 2025
     TOTAL  Time  : 179
     SEE INFORMATION IN : OUT.Mn_ecutwfc_90/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:30:22 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_95/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 45 * 45 * 45
     UNIFORM GRID DIM(BIG)   : 45 * 45 * 45
     DONE(0.0708866  SEC) : SETUP UNITCELL
     DONE(0.21378    SEC) : SYMMETRY
     DONE(0.516102   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.526637   SEC) : INIT PLANEWAVE
     DONE(0.538241   SEC) : LOCAL POTENTIAL
     DONE(0.594484   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 37.3413
     DONE(0.594586   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.35465    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08163655e+04   0.00000000e+00   6.1900e-01  28.75
     CG2     -1.08156998e+04   6.65774657e-01   7.5948e+00  13.46
     CG3     -1.08217481e+04  -6.04832767e+00   2.1635e+01   8.69
     CG4     -1.08281547e+04  -6.40663490e+00   1.0215e+01   7.97
     CG5     -1.08314062e+04  -3.25149606e+00   3.9249e-01   6.18
     CG6     -1.08314308e+04  -2.45820580e-02   3.0648e-01   6.09
     CG7     -1.08314895e+04  -5.87507669e-02   8.6027e-02   6.18
     CG8     -1.08315736e+04  -8.40916438e-02   3.5762e-02   6.96
     CG9     -1.08316265e+04  -5.28162400e-02   1.3571e-02   6.92
     CG10    -1.08316451e+04  -1.86014777e-02   2.2119e-03   6.98
     CG11    -1.08316508e+04  -5.77499931e-03   6.8563e-04   7.91
     CG12    -1.08316517e+04  -8.41460450e-04   1.5578e-04   7.30
     CG13    -1.08316517e+04  -4.24702229e-05   2.0160e-04   8.40
     CG14    -1.08316520e+04  -2.72652146e-04   2.6992e-05   6.65
     CG15    -1.08316520e+04  -5.72450023e-05   2.6786e-05   7.83
     CG16    -1.08316521e+04  -4.81228882e-06   1.5599e-05   6.05
     CG17    -1.08316521e+04  -1.01312855e-05   2.4944e-05   6.13
     CG18    -1.08316521e+04  -1.01566776e-05   2.4261e-07   6.10
     CG19    -1.08316521e+04  -5.16570686e-06   2.3404e-07  11.52
     CG20    -1.08316521e+04   5.99611124e-07   7.1390e-08   6.35
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       174.90 17       10.29  100.00 
     Driver           reading                     0.04   1        0.04   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 174.86 1        174.86 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.08   
     PW_Basis_K       setuptransform              0.01   1        0.01   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.01   1        0.01   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    174.28 1        174.28 99.64  
     ESolver_KS_PW    runner                      174.28 1        174.28 99.64  
     ESolver_KS_PW    before_scf                  1.76   1        1.76   1.01   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.02   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.24   148      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.09   148      0.00   0.05   
     Potential        init_pot                    0.05   1        0.05   0.03   
     Potential        update_from_charge          0.87   21       0.04   0.50   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.87   21       0.04   0.50   
     H_Hartree_pw     v_hartree                   0.08   21       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.28   188      0.00   0.16   
     PW_Basis_Sup     gatherp_scatters            0.10   188      0.00   0.06   
     PotXC            cal_v_eff                   0.78   21       0.04   0.45   
     XC_Functional    v_xc                        0.78   21       0.04   0.45   
     Potential        interpolate_vrs             0.00   21       0.00   0.00   
     Symmetry         rhog_symmetry               0.28   21       0.01   0.16   
     Symmetry         group fft grids             0.06   21       0.00   0.03   
     PSIInit          initialize_psi              1.68   1        1.68   0.96   
     Nonlocal         getvnl                      1.41   420      0.00   0.80   
     pp_cell_vnl      getvnl                      1.40   420      0.00   0.80   
     Structure_Factor get_sk                      0.14   420      0.00   0.08   
     DiagoIterAssist  diagH_subspace              26.80  400      0.07   15.32  
     Operator         hPsi                        132.14 51976    0.00   75.55  
     Operator         EkineticPW                  0.75   51976    0.00   0.43   
     Operator         VeffPW                      91.74  51976    0.00   52.46  
     PW_Basis_K       recip2real                  55.35  83576    0.00   31.64  
     PW_Basis_K       gathers_scatterp            17.16  83576    0.00   9.81   
     PW_Basis_K       real2recip                  38.42  67576    0.00   21.97  
     PW_Basis_K       gatherp_scatters            8.64   67576    0.00   4.94   
     Operator         NonlocalPW                  39.52  51976    0.00   22.60  
     Nonlocal         add_nonlocal_pp             19.69  51976    0.00   11.26  
     DiagoIterAssist  diagH_LAPACK                0.26   400      0.00   0.15   
     ESolver_KS_PW    hamilt2density_single       171.42 20       8.57   98.01  
     HSolverPW        solve                       171.04 20       8.55   97.80  
     DiagoCG          diag_once                   130.90 400      0.33   74.84  
     DiagoCG_New      spsi_func                   0.82   103152   0.00   0.47   
     DiagoCG_New      hpsi_func                   108.11 51576    0.00   61.82  
     ElecStatePW      psiToRho                    12.56  20       0.63   7.18   
     Charge_Mixing    get_drho                    0.06   20       0.00   0.04   
     Charge_Mixing    inner_product_recip_rho     0.00   20       0.00   0.00   
     Charge           mix_rho                     0.10   19       0.01   0.06   
     Charge           Broyden_mixing              0.04   19       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.03   232      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.09   1        0.09   0.05   
     ModuleIO         write_rhog                  0.07   1        0.07   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:30:22 2025
     FINISH Time  : Sat Apr 26 12:33:17 2025
     TOTAL  Time  : 175
     SEE INFORMATION IN : OUT.Mn_ecutwfc_95/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:33:18 2025
     MAKE THE DIR         : OUT.Mn_ecutwfc_100/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 48 * 48 * 48
     UNIFORM GRID DIM(BIG)   : 48 * 48 * 48
     DONE(0.0722388  SEC) : SETUP UNITCELL
     DONE(0.212517   SEC) : SYMMETRY
     DONE(0.51891    SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.529906   SEC) : INIT PLANEWAVE
     DONE(0.5457     SEC) : LOCAL POTENTIAL
     DONE(0.604199   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 40.8936
     DONE(0.604298   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.49993    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08157533e+04   0.00000000e+00   6.4699e-01  29.66
     CG2     -1.08149209e+04   8.32471797e-01   8.1580e+00  14.34
     CG3     -1.08225367e+04  -7.61582341e+00   1.9525e+01   8.95
     CG4     -1.08295175e+04  -6.98079157e+00   7.1023e+00   8.31
     CG5     -1.08314294e+04  -1.91189427e+00   3.7767e-01   6.31
     CG6     -1.08314484e+04  -1.89798834e-02   3.7326e-01   6.25
     CG7     -1.08315144e+04  -6.60798406e-02   8.9703e-02   6.31
     CG8     -1.08315686e+04  -5.41829980e-02   3.6583e-02   6.71
     CG9     -1.08316247e+04  -5.61075809e-02   5.9551e-03   7.41
     CG10    -1.08316530e+04  -2.83139504e-02   1.5831e-03   8.73
     CG11    -1.08316570e+04  -3.92393368e-03   1.2021e-03   7.30
     CG12    -1.08316543e+04   2.70491629e-03   1.4176e-03   6.40
     CG13    -1.08316542e+04   6.14301302e-05   3.9004e-04   6.34
     CG14    -1.08316541e+04   6.01142244e-05   1.7403e-04   6.89
     CG15    -1.08316543e+04  -2.12472925e-04   9.3094e-05   7.45
     CG16    -1.08316544e+04  -7.06816564e-05   3.3041e-04   6.83
     CG17    -1.08316543e+04   7.87935906e-05   6.4576e-05   6.29
     CG18    -1.08316544e+04  -3.53840880e-05   1.5282e-05   6.53
     CG19    -1.08316544e+04  -2.28525374e-05   4.2288e-06   7.86
     CG20    -1.08316544e+04  -2.67546478e-05   1.5220e-07   9.05
     CG21    -1.08316544e+04  -1.59649174e-06   3.0992e-07  11.15
     CG22    -1.08316544e+04   7.82461666e-07   8.7247e-08   6.62
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       194.32 17       11.43  100.00 
     Driver           reading                     0.04   1        0.04   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 194.29 1        194.29 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.07   
     PW_Basis_K       setuptransform              0.01   1        0.01   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.01   1        0.01   0.00   
     ppcell_vnl       init                        0.01   1        0.01   0.00   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    193.69 1        193.69 99.68  
     ESolver_KS_PW    runner                      193.69 1        193.69 99.68  
     ESolver_KS_PW    before_scf                  1.90   1        1.90   0.98   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.02   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.27   162      0.00   0.14   
     PW_Basis_Sup     gathers_scatterp            0.11   162      0.00   0.06   
     Potential        init_pot                    0.05   1        0.05   0.03   
     Potential        update_from_charge          1.04   23       0.05   0.53   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   1.03   23       0.04   0.53   
     H_Hartree_pw     v_hartree                   0.09   23       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.28   206      0.00   0.14   
     PW_Basis_Sup     gatherp_scatters            0.11   206      0.00   0.06   
     PotXC            cal_v_eff                   0.93   23       0.04   0.48   
     XC_Functional    v_xc                        0.93   23       0.04   0.48   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.34   23       0.01   0.17   
     Symmetry         group fft grids             0.07   23       0.00   0.04   
     PSIInit          initialize_psi              1.80   1        1.80   0.92   
     Nonlocal         getvnl                      1.64   460      0.00   0.84   
     pp_cell_vnl      getvnl                      1.64   460      0.00   0.84   
     Structure_Factor get_sk                      0.16   460      0.00   0.08   
     DiagoIterAssist  diagH_subspace              30.28  440      0.07   15.58  
     Operator         hPsi                        145.40 55268    0.00   74.82  
     Operator         EkineticPW                  0.86   55268    0.00   0.44   
     Operator         VeffPW                      98.22  55268    0.00   50.54  
     PW_Basis_K       recip2real                  59.59  90028    0.00   30.67  
     PW_Basis_K       gathers_scatterp            21.55  90028    0.00   11.09  
     PW_Basis_K       real2recip                  40.01  72428    0.00   20.59  
     PW_Basis_K       gatherp_scatters            10.65  72428    0.00   5.48   
     Operator         NonlocalPW                  46.18  55268    0.00   23.76  
     Nonlocal         add_nonlocal_pp             23.03  55268    0.00   11.85  
     DiagoIterAssist  diagH_LAPACK                0.29   440      0.00   0.15   
     ESolver_KS_PW    hamilt2density_single       190.51 22       8.66   98.04  
     HSolverPW        solve                       190.04 22       8.64   97.80  
     DiagoCG          diag_once                   144.46 440      0.33   74.34  
     DiagoCG_New      spsi_func                   0.92   109656   0.00   0.48   
     DiagoCG_New      hpsi_func                   118.51 54828    0.00   60.99  
     ElecStatePW      psiToRho                    14.23  22       0.65   7.32   
     Charge_Mixing    get_drho                    0.06   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.12   21       0.01   0.06   
     Charge           Broyden_mixing              0.05   21       0.00   0.03   
     Charge_Mixing    inner_product_recip_hartree 0.05   264      0.00   0.03   
     ESolver_KS_PW    after_scf                   0.09   1        0.09   0.05   
     ModuleIO         write_rhog                  0.08   1        0.08   0.04   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:33:18 2025
     FINISH Time  : Sat Apr 26 12:36:32 2025
     TOTAL  Time  : 194
     SEE INFORMATION IN : OUT.Mn_ecutwfc_100/
    
    统计结果已保存到 ./Mn/final_etot_vs_ecutwfc.csv
    Figure(800x600)
    

### 输出图片
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/475862/e010d61ade7245c499cac194228dc856/5fccda54-31b2-4bad-b5ef-a64355a51f45.png)

## K点收敛性测试
将k点从2增大到10，步长为1，ecutwfc设置为60。同样利用python脚本自动提交
```python
import os
import re
import matplotlib.pyplot as plt

# 文件路径
kpt_file_path = "./Mn/KPT"
input_file_path = "./Mn/INPUT"
abacus_command = "mpirun -np 2 abacus"  # 假设 abacus 可执行文件在 PATH 中
working_directory = "./Mn"  # 指定工作目录

# 定义网格细分范围
k_points_range = range(2, 11)  # 网格细分从 2 到 10
energy_data = []  # 用于存储网格细分和对应的总能量

# 读取原始 KPT 文件内容
with open(kpt_file_path, "r") as file:
    kpt_content = file.readlines()

# 读取原始 INPUT 文件内容
with open(input_file_path, "r") as file:
    input_content = file.readlines()

# 修改网格细分和输出后缀并运行 abacus
for k_points in k_points_range:
    # 修改 KPT 文件
    new_kpt_content = kpt_content[:]  # 复制原始内容
    new_kpt_content[3] = f"{k_points} {k_points} {k_points} 0 0 0\n"  # 修改第四行的网格细分

    # 写入修改后的内容到 KPT 文件
    with open(kpt_file_path, "w") as file:
        file.writelines(new_kpt_content)

    # 修改 INPUT 文件
    new_input_content = []
    for line in input_content:
        if line.startswith("suffix"):
            new_input_content.append(f"suffix       Mn_kpoints_{k_points}                     # 输出后缀\n")
        else:
            new_input_content.append(line)

    # 写入修改后的内容到 INPUT 文件
    with open(input_file_path, "w") as file:
        file.writelines(new_input_content)
    
    os.system(f"cd {working_directory} && {abacus_command}")

    # 提取 running_scf.log 中的 "final etot" 总能量
    log_file = os.path.join(working_directory, f"OUT.Mn_kpoints_{k_points}/running_scf.log")
    if not os.path.exists(log_file):
        print(f"警告: 找不到文件 {log_file}")
        continue

    # 读取 running_scf.log 文件
    with open(log_file, "r") as file:
        lines = file.readlines()

    # 提取 "final etot" 的总能量值
    final_etot = None
    for line in lines:
        match = re.search(r"final etot is ([\-\d\.]+) eV", line)
        if match:
            final_etot = float(match.group(1))
            break

    if final_etot is not None:
        energy_data.append((k_points, final_etot))
    else:
        print(f"警告: 在文件 {log_file} 中未找到 'final etot' 数据")

# 保存统计数据为 CSV 文件
output_csv = os.path.join(working_directory, "final_etot_vs_kpoints.csv")
with open(output_csv, "w") as file:
    file.write("k-points,final_etot (eV)\n")
    for k_points, final_etot in energy_data:
        file.write(f"{k_points},{final_etot:.6f}\n")

print(f"\n统计结果已保存到 {output_csv}")

# 绘制折线图
k_points_values, final_etot_values = zip(*energy_data)
plt.figure(figsize=(8, 6))
plt.plot(k_points_values, final_etot_values, marker='o', linestyle='-', color='b', label='Final Total Energy')
plt.xlabel('k-points')
plt.ylabel('Final Total Energy (eV)')
plt.title('Final Total Energy vs k-points')
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(working_directory, "final_etot_vs_kpoints.png"))
plt.show()
```


```
#运行脚本
! python kpoint.py
```

                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:49:35 2025
     MAKE THE DIR         : OUT.Mn_kpoints_2/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0704018  SEC) : SETUP UNITCELL
     DONE(0.210953   SEC) : SYMMETRY
     DONE(0.497485   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       4               2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.499564   SEC) : INIT PLANEWAVE
     DONE(0.50607    SEC) : LOCAL POTENTIAL
     DONE(0.555051   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 3.7793
     DONE(0.555152   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.768357   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08261816e+04   0.00000000e+00   3.5147e-01   3.57
     CG2     -1.08216966e+04   4.48500323e+00   5.0510e+00   0.94
     CG3     -1.08234983e+04  -1.80171043e+00   1.9549e+01   0.80
     CG4     -1.08295062e+04  -6.00789385e+00   4.9243e+00   0.69
     CG5     -1.08307893e+04  -1.28313472e+00   2.4853e-01   0.57
     CG6     -1.08308337e+04  -4.43575740e-02   1.5613e-01   0.57
     CG7     -1.08308930e+04  -5.92807755e-02   5.7399e-02   0.56
     CG8     -1.08309286e+04  -3.55967091e-02   3.6243e-02   0.69
     CG9     -1.08309374e+04  -8.85201382e-03   1.0592e-02   0.56
     CG10    -1.08309583e+04  -2.08416637e-02   3.1817e-03   0.67
     CG11    -1.08309638e+04  -5.50233443e-03   1.9087e-04   0.70
     CG12    -1.08309654e+04  -1.61136232e-03   1.4608e-04   0.94
     CG13    -1.08309653e+04   7.08787977e-05   1.8880e-04   0.56
     CG14    -1.08309650e+04   2.81529335e-04   1.7243e-03   0.63
     CG15    -1.08309655e+04  -5.08070631e-04   7.7122e-06   0.57
     CG16    -1.08309656e+04  -2.99436947e-05   1.2077e-05   0.79
     CG17    -1.08309656e+04  -5.64815801e-06   3.9597e-06   0.57
     CG18    -1.08309656e+04  -2.37390160e-06   3.4006e-07   0.56
     CG19    -1.08309656e+04  -2.34653272e-06   4.3023e-08   0.92
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       16.73  17       0.98   100.00 
     Driver           reading                     0.04   1        0.04   0.25   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 16.68  1        16.68  99.75  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.01   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.01   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.84   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.01   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.01   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.02   
     ppcell_vnl       init                        0.00   1        0.00   0.02   
     ppcell_vnl       init_vnl                    0.05   1        0.05   0.27   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    16.15  1        16.15  96.57  
     ESolver_KS_PW    runner                      16.15  1        16.15  96.56  
     ESolver_KS_PW    before_scf                  0.21   1        0.21   1.27   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.07   
     PW_Basis_Sup     recip2real                  0.09   142      0.00   0.52   
     PW_Basis_Sup     gathers_scatterp            0.03   142      0.00   0.19   
     Potential        init_pot                    0.02   1        0.02   0.12   
     Potential        update_from_charge          0.37   20       0.02   2.24   
     Potential        cal_fixed_v                 0.00   1        0.00   0.01   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.01   
     Potential        cal_v_eff                   0.37   20       0.02   2.22   
     H_Hartree_pw     v_hartree                   0.03   20       0.00   0.20   
     PW_Basis_Sup     real2recip                  0.10   182      0.00   0.61   
     PW_Basis_Sup     gatherp_scatters            0.04   182      0.00   0.24   
     PotXC            cal_v_eff                   0.34   20       0.02   2.01   
     XC_Functional    v_xc                        0.34   20       0.02   2.00   
     Potential        interpolate_vrs             0.00   20       0.00   0.01   
     Symmetry         rhog_symmetry               0.15   21       0.01   0.88   
     Symmetry         group fft grids             0.03   21       0.00   0.18   
     PSIInit          initialize_psi              0.18   1        0.18   1.06   
     Nonlocal         getvnl                      0.14   84       0.00   0.82   
     pp_cell_vnl      getvnl                      0.14   84       0.00   0.82   
     Structure_Factor get_sk                      0.01   84       0.00   0.09   
     DiagoIterAssist  diagH_subspace              2.24   76       0.03   13.39  
     Operator         hPsi                        11.56  10311    0.00   69.12  
     Operator         EkineticPW                  0.08   10311    0.00   0.47   
     Operator         VeffPW                      7.49   10311    0.00   44.75  
     PW_Basis_K       recip2real                  4.35   16475    0.00   26.02  
     PW_Basis_K       gathers_scatterp            1.71   16475    0.00   10.22  
     PW_Basis_K       real2recip                  3.22   13275    0.00   19.22  
     PW_Basis_K       gatherp_scatters            0.99   13275    0.00   5.89   
     Operator         NonlocalPW                  3.97   10311    0.00   23.75  
     Nonlocal         add_nonlocal_pp             1.98   10311    0.00   11.86  
     DiagoIterAssist  diagH_LAPACK                0.05   76       0.00   0.30   
     ESolver_KS_PW    hamilt2density_single       15.44  20       0.77   92.30  
     HSolverPW        solve                       15.25  20       0.76   91.16  
     DiagoCG          diag_once                   11.90  80       0.15   71.16  
     DiagoCG_New      spsi_func                   0.09   20470    0.00   0.54   
     DiagoCG_New      hpsi_func                   9.64   10235    0.00   57.65  
     ElecStatePW      psiToRho                    1.03   20       0.05   6.18   
     Charge_Mixing    get_drho                    0.03   20       0.00   0.17   
     Charge_Mixing    inner_product_recip_rho     0.00   20       0.00   0.01   
     Charge           mix_rho                     0.04   18       0.00   0.25   
     Charge           Broyden_mixing              0.02   18       0.00   0.10   
     Charge_Mixing    inner_product_recip_hartree 0.01   216      0.00   0.09   
     ESolver_KS_PW    after_scf                   0.06   1        0.06   0.38   
     ModuleIO         write_rhog                  0.06   1        0.06   0.34   
     ModuleIO         write_istate_info           0.01   1        0.01   0.06   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:49:35 2025
     FINISH Time  : Sat Apr 26 12:49:52 2025
     TOTAL  Time  : 17
     SEE INFORMATION IN : OUT.Mn_kpoints_2/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:49:53 2025
     MAKE THE DIR         : OUT.Mn_kpoints_3/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0641928  SEC) : SETUP UNITCELL
     DONE(0.21437    SEC) : SYMMETRY
     DONE(0.504046   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       4               2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.506079   SEC) : INIT PLANEWAVE
     DONE(0.512414   SEC) : LOCAL POTENTIAL
     DONE(0.557033   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 3.74512
     DONE(0.557122   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.743834   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08158356e+04   0.00000000e+00   6.8286e-01   2.67
     CG2     -1.08084646e+04   7.37094502e+00   1.8431e+01   1.21
     CG3     -1.08257074e+04  -1.72427518e+01   1.2147e+01   0.81
     CG4     -1.08305111e+04  -4.80368384e+00   2.4796e+00   0.66
     CG5     -1.08314397e+04  -9.28646405e-01   1.6002e-01   0.59
     CG6     -1.08316881e+04  -2.48365399e-01   3.4272e-02   0.63
     CG7     -1.08317482e+04  -6.00954027e-02   1.2983e-01   0.74
     CG8     -1.08317941e+04  -4.59093559e-02   9.1103e-02   0.64
     CG9     -1.08317726e+04   2.14441924e-02   5.8992e-02   0.57
     CG10    -1.08317738e+04  -1.17189476e-03   1.1138e-02   0.55
     CG11    -1.08317676e+04   6.19488253e-03   7.2863e-03   0.56
     CG12    -1.08317464e+04   2.12498475e-02   5.3329e-03   0.56
     CG13    -1.08317522e+04  -5.80910832e-03   6.7015e-04   0.62
     CG14    -1.08317616e+04  -9.45076210e-03   1.7784e-03   1.04
     CG15    -1.08317597e+04   1.90969452e-03   5.8694e-04   0.63
     CG16    -1.08317610e+04  -1.24898805e-03   1.1584e-04   0.68
     CG17    -1.08317617e+04  -7.84047784e-04   2.5023e-05   0.84
     CG18    -1.08317617e+04   4.53738269e-06   1.4337e-04   0.69
     CG19    -1.08317618e+04  -6.25427663e-05   4.4214e-06   0.56
     CG20    -1.08317618e+04  -1.17395329e-05   1.2322e-06   0.71
     CG21    -1.08317618e+04  -4.61887990e-06   3.3794e-08   0.76
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       17.58  17       1.03   100.00 
     Driver           reading                     0.04   1        0.04   0.21   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 17.54  1        17.54  99.79  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.01   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.01   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.15   1        0.15   0.85   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.01   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.01   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.02   
     ppcell_vnl       init                        0.00   1        0.00   0.02   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.23   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    17.00  1        17.00  96.73  
     ESolver_KS_PW    runner                      17.00  1        17.00  96.73  
     ESolver_KS_PW    before_scf                  0.19   1        0.19   1.06   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.06   
     PW_Basis_Sup     recip2real                  0.10   155      0.00   0.57   
     PW_Basis_Sup     gathers_scatterp            0.04   155      0.00   0.21   
     Potential        init_pot                    0.02   1        0.02   0.13   
     Potential        update_from_charge          0.43   22       0.02   2.45   
     Potential        cal_fixed_v                 0.00   1        0.00   0.01   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.43   22       0.02   2.44   
     H_Hartree_pw     v_hartree                   0.03   22       0.00   0.20   
     PW_Basis_Sup     real2recip                  0.12   197      0.00   0.65   
     PW_Basis_Sup     gatherp_scatters            0.05   197      0.00   0.27   
     PotXC            cal_v_eff                   0.39   22       0.02   2.22   
     XC_Functional    v_xc                        0.39   22       0.02   2.22   
     Potential        interpolate_vrs             0.00   22       0.00   0.01   
     Symmetry         rhog_symmetry               0.15   22       0.01   0.87   
     Symmetry         group fft grids             0.03   22       0.00   0.19   
     PSIInit          initialize_psi              0.15   1        0.15   0.84   
     Nonlocal         getvnl                      0.15   88       0.00   0.84   
     pp_cell_vnl      getvnl                      0.15   88       0.00   0.83   
     Structure_Factor get_sk                      0.02   88       0.00   0.09   
     DiagoIterAssist  diagH_subspace              2.43   84       0.03   13.80  
     Operator         hPsi                        12.12  10832    0.00   68.98  
     Operator         EkineticPW                  0.08   10832    0.00   0.46   
     Operator         VeffPW                      7.84   10832    0.00   44.58  
     PW_Basis_K       recip2real                  4.67   17468    0.00   26.55  
     PW_Basis_K       gathers_scatterp            1.87   17468    0.00   10.64  
     PW_Basis_K       real2recip                  3.28   14108    0.00   18.68  
     PW_Basis_K       gatherp_scatters            0.93   14108    0.00   5.28   
     Operator         NonlocalPW                  4.18   10832    0.00   23.80  
     Nonlocal         add_nonlocal_pp             2.11   10832    0.00   12.03  
     DiagoIterAssist  diagH_LAPACK                0.06   84       0.00   0.36   
     ESolver_KS_PW    hamilt2density_single       16.25  21       0.77   92.48  
     HSolverPW        solve                       16.06  21       0.76   91.35  
     DiagoCG          diag_once                   12.38  84       0.15   70.43  
     DiagoCG_New      spsi_func                   0.09   21496    0.00   0.53   
     DiagoCG_New      hpsi_func                   10.02  10748    0.00   57.02  
     ElecStatePW      psiToRho                    1.14   21       0.05   6.50   
     Charge_Mixing    get_drho                    0.03   21       0.00   0.15   
     Charge_Mixing    inner_product_recip_rho     0.00   21       0.00   0.01   
     Charge           mix_rho                     0.05   20       0.00   0.29   
     Charge           Broyden_mixing              0.02   20       0.00   0.12   
     Charge_Mixing    inner_product_recip_hartree 0.02   248      0.00   0.10   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.37   
     ModuleIO         write_rhog                  0.06   1        0.06   0.33   
     ModuleIO         write_istate_info           0.01   1        0.01   0.05   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:49:53 2025
     FINISH Time  : Sat Apr 26 12:50:11 2025
     TOTAL  Time  : 18
     SEE INFORMATION IN : OUT.Mn_kpoints_3/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:50:12 2025
     MAKE THE DIR         : OUT.Mn_kpoints_4/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0650619  SEC) : SETUP UNITCELL
     DONE(0.228519   SEC) : SYMMETRY
     DONE(0.522255   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       10              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.525593   SEC) : INIT PLANEWAVE
     DONE(0.531914   SEC) : LOCAL POTENTIAL
     DONE(0.576481   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 9.44824
     DONE(0.576561   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.983404   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08171065e+04   0.00000000e+00   6.1678e-01   6.49
     CG2     -1.08169196e+04   1.86900700e-01   6.6019e+00   3.13
     CG3     -1.08206942e+04  -3.77460219e+00   2.5579e+01   1.96
     CG4     -1.08299570e+04  -9.26288379e+00   4.5453e+00   1.70
     CG5     -1.08313547e+04  -1.39769073e+00   2.1870e-01   1.44
     CG6     -1.08314568e+04  -1.02099697e-01   1.0960e-01   1.44
     CG7     -1.08315233e+04  -6.65047102e-02   4.6532e-02   1.45
     CG8     -1.08316065e+04  -8.31463157e-02   5.6485e-03   1.56
     CG9     -1.08316369e+04  -3.03834011e-02   4.2370e-03   2.09
     CG10    -1.08316378e+04  -9.21537414e-04   2.3236e-03   1.41
     CG11    -1.08316387e+04  -9.50361011e-04   3.5835e-04   1.40
     CG12    -1.08316400e+04  -1.25803925e-03   2.3274e-05   1.96
     CG13    -1.08316402e+04  -2.44767192e-04   2.8222e-05   2.37
     CG14    -1.08316402e+04   3.70796206e-06   1.2139e-05   1.37
     CG15    -1.08316403e+04  -1.01453551e-05   1.4001e-07   1.40
     CG16    -1.08316403e+04  -9.60021881e-06   8.3494e-07   3.16
     CG17    -1.08316403e+04   1.79095248e-06   5.1852e-07   1.41
     CG18    -1.08316403e+04  -8.96317695e-08   7.8070e-08   1.97
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       38.78  17       2.28   100.00 
     Driver           reading                     0.04   1        0.04   0.09   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 38.75  1        38.75  99.91  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.16   1        0.16   0.42   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.01   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.01   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.10   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    38.19  1        38.19  98.46  
     ESolver_KS_PW    runner                      38.19  1        38.19  98.46  
     ESolver_KS_PW    before_scf                  0.41   1        0.41   1.05   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.02   
     PW_Basis_Sup     recip2real                  0.08   134      0.00   0.22   
     PW_Basis_Sup     gathers_scatterp            0.03   134      0.00   0.08   
     Potential        init_pot                    0.02   1        0.02   0.06   
     Potential        update_from_charge          0.36   19       0.02   0.93   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.36   19       0.02   0.92   
     H_Hartree_pw     v_hartree                   0.03   19       0.00   0.08   
     PW_Basis_Sup     real2recip                  0.09   170      0.00   0.24   
     PW_Basis_Sup     gatherp_scatters            0.03   170      0.00   0.09   
     PotXC            cal_v_eff                   0.32   19       0.02   0.83   
     XC_Functional    v_xc                        0.32   19       0.02   0.83   
     Potential        interpolate_vrs             0.00   19       0.00   0.00   
     Symmetry         rhog_symmetry               0.14   19       0.01   0.35   
     Symmetry         group fft grids             0.03   19       0.00   0.07   
     PSIInit          initialize_psi              0.37   1        0.37   0.95   
     Nonlocal         getvnl                      0.32   190      0.00   0.83   
     pp_cell_vnl      getvnl                      0.32   190      0.00   0.83   
     Structure_Factor get_sk                      0.05   190      0.00   0.12   
     DiagoIterAssist  diagH_subspace              5.26   180      0.03   13.57  
     Operator         hPsi                        28.09  25304    0.00   72.43  
     Operator         EkineticPW                  0.19   25304    0.00   0.49   
     Operator         VeffPW                      18.17  25304    0.00   46.86  
     PW_Basis_K       recip2real                  10.54  39524    0.00   27.18  
     PW_Basis_K       gathers_scatterp            4.17   39524    0.00   10.75  
     PW_Basis_K       real2recip                  7.73   32324    0.00   19.92  
     PW_Basis_K       gatherp_scatters            2.20   32324    0.00   5.66   
     Operator         NonlocalPW                  9.67   25304    0.00   24.94  
     Nonlocal         add_nonlocal_pp             4.89   25304    0.00   12.62  
     DiagoIterAssist  diagH_LAPACK                0.12   180      0.00   0.31   
     ESolver_KS_PW    hamilt2density_single       37.30  18       2.07   96.18  
     HSolverPW        solve                       37.13  18       2.06   95.73  
     DiagoCG          diag_once                   29.25  180      0.16   75.42  
     DiagoCG_New      spsi_func                   0.22   50248    0.00   0.57   
     DiagoCG_New      hpsi_func                   23.55  25124    0.00   60.72  
     ElecStatePW      psiToRho                    2.42   18       0.13   6.25   
     Charge_Mixing    get_drho                    0.02   18       0.00   0.06   
     Charge_Mixing    inner_product_recip_rho     0.00   18       0.00   0.00   
     Charge           mix_rho                     0.04   17       0.00   0.10   
     Charge           Broyden_mixing              0.02   17       0.00   0.04   
     Charge_Mixing    inner_product_recip_hartree 0.01   200      0.00   0.04   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.18   
     ModuleIO         write_rhog                  0.06   1        0.06   0.16   
     ModuleIO         write_istate_info           0.01   1        0.01   0.03   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:50:12 2025
     FINISH Time  : Sat Apr 26 12:50:50 2025
     TOTAL  Time  : 38
     SEE INFORMATION IN : OUT.Mn_kpoints_4/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:50:51 2025
     MAKE THE DIR         : OUT.Mn_kpoints_5/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0712467  SEC) : SETUP UNITCELL
     DONE(0.211703   SEC) : SYMMETRY
     DONE(0.49674    SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       10              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.500011   SEC) : INIT PLANEWAVE
     DONE(0.506366   SEC) : LOCAL POTENTIAL
     DONE(0.549932   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 9.31396
     DONE(0.550008   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.948809   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08157652e+04   0.00000000e+00   6.6584e-01   6.54
     CG2     -1.08180932e+04  -2.32800091e+00   5.4270e+00   3.09
     CG3     -1.08212683e+04  -3.17512510e+00   2.3056e+01   1.80
     CG4     -1.08303634e+04  -9.09507802e+00   2.3861e+00   1.74
     CG5     -1.08311799e+04  -8.16501043e-01   1.7944e-01   1.39
     CG6     -1.08315105e+04  -3.30638498e-01   4.6712e-02   1.56
     CG7     -1.08315786e+04  -6.80697579e-02   1.1944e-01   1.70
     CG8     -1.08316342e+04  -5.56218438e-02   1.3268e-02   1.37
     CG9     -1.08316326e+04   1.64796397e-03   5.4989e-02   1.47
     CG10    -1.08316338e+04  -1.22685445e-03   1.3614e-02   1.47
     CG11    -1.08316298e+04   4.01962765e-03   7.4397e-03   1.35
     CG12    -1.08316260e+04   3.73826583e-03   4.7919e-03   1.40
     CG13    -1.08316224e+04   3.66644969e-03   1.6262e-03   1.44
     CG14    -1.08316326e+04  -1.02600145e-02   7.6450e-03   2.37
     CG15    -1.08316344e+04  -1.71107921e-03   1.6753e-03   1.49
     CG16    -1.08316366e+04  -2.24943280e-03   9.9134e-04   1.42
     CG17    -1.08316375e+04  -8.83263124e-04   1.8001e-03   1.38
     CG18    -1.08316344e+04   3.11188305e-03   1.6694e-03   1.41
     CG19    -1.08316332e+04   1.21343336e-03   9.8671e-04   1.35
     CG20    -1.08316320e+04   1.20118012e-03   3.6470e-04   1.41
     CG21    -1.08316339e+04  -1.91807037e-03   3.0786e-05   1.84
     CG22    -1.08316340e+04  -1.17600827e-04   7.2797e-06   1.90
     CG23    -1.08316340e+04  -1.03265320e-05   3.1729e-06   1.67
     CG24    -1.08316340e+04   6.80773452e-06   2.3893e-06   1.46
     CG25    -1.08316340e+04  -8.03711450e-09   7.0556e-07   1.41
     CG26    -1.08316340e+04  -2.39817535e-06   9.3892e-08   1.88
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       48.36  17       2.84   100.00 
     Driver           reading                     0.04   1        0.04   0.09   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 48.32  1        48.32  99.91  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.29   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.01   
     ppcell_vnl       init                        0.00   1        0.00   0.01   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.08   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    47.79  1        47.79  98.82  
     ESolver_KS_PW    runner                      47.79  1        47.79  98.82  
     ESolver_KS_PW    before_scf                  0.40   1        0.40   0.82   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.03   
     PW_Basis_Sup     recip2real                  0.12   190      0.00   0.25   
     PW_Basis_Sup     gathers_scatterp            0.05   190      0.00   0.10   
     Potential        init_pot                    0.02   1        0.02   0.04   
     Potential        update_from_charge          0.52   27       0.02   1.07   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.52   27       0.02   1.07   
     H_Hartree_pw     v_hartree                   0.04   27       0.00   0.09   
     PW_Basis_Sup     real2recip                  0.14   242      0.00   0.29   
     PW_Basis_Sup     gatherp_scatters            0.05   242      0.00   0.11   
     PotXC            cal_v_eff                   0.47   27       0.02   0.97   
     XC_Functional    v_xc                        0.47   27       0.02   0.97   
     Potential        interpolate_vrs             0.00   27       0.00   0.00   
     Symmetry         rhog_symmetry               0.19   27       0.01   0.39   
     Symmetry         group fft grids             0.04   27       0.00   0.08   
     PSIInit          initialize_psi              0.36   1        0.36   0.74   
     Nonlocal         getvnl                      0.43   270      0.00   0.89   
     pp_cell_vnl      getvnl                      0.43   270      0.00   0.89   
     Structure_Factor get_sk                      0.05   270      0.00   0.10   
     DiagoIterAssist  diagH_subspace              7.54   260      0.03   15.60  
     Operator         hPsi                        35.00  30616    0.00   72.38  
     Operator         EkineticPW                  0.23   30616    0.00   0.48   
     Operator         VeffPW                      22.94  30616    0.00   47.44  
     PW_Basis_K       recip2real                  13.60  51156    0.00   28.13  
     PW_Basis_K       gathers_scatterp            5.35   51156    0.00   11.06  
     PW_Basis_K       real2recip                  9.72   40756    0.00   20.09  
     PW_Basis_K       gatherp_scatters            2.81   40756    0.00   5.81   
     Operator         NonlocalPW                  11.76  30616    0.00   24.31  
     Nonlocal         add_nonlocal_pp             5.93   30616    0.00   12.27  
     DiagoIterAssist  diagH_LAPACK                0.17   260      0.00   0.36   
     ESolver_KS_PW    hamilt2density_single       46.71  26       1.80   96.59  
     HSolverPW        solve                       46.46  26       1.79   96.07  
     DiagoCG          diag_once                   35.14  260      0.14   72.67  
     DiagoCG_New      spsi_func                   0.26   60712    0.00   0.53   
     DiagoCG_New      hpsi_func                   28.48  30356    0.00   58.89  
     ElecStatePW      psiToRho                    3.38   26       0.13   6.99   
     Charge_Mixing    get_drho                    0.04   26       0.00   0.07   
     Charge_Mixing    inner_product_recip_rho     0.00   26       0.00   0.01   
     Charge           mix_rho                     0.06   25       0.00   0.13   
     Charge           Broyden_mixing              0.03   25       0.00   0.05   
     Charge_Mixing    inner_product_recip_hartree 0.02   328      0.00   0.05   
     ESolver_KS_PW    after_scf                   0.06   1        0.06   0.13   
     ModuleIO         write_rhog                  0.06   1        0.06   0.12   
     ModuleIO         write_istate_info           0.01   1        0.01   0.02   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:50:51 2025
     FINISH Time  : Sat Apr 26 12:51:40 2025
     TOTAL  Time  : 49
     SEE INFORMATION IN : OUT.Mn_kpoints_5/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:51:41 2025
     MAKE THE DIR         : OUT.Mn_kpoints_6/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0630941  SEC) : SETUP UNITCELL
     DONE(0.203272   SEC) : SYMMETRY
     DONE(0.503455   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.508866   SEC) : INIT PLANEWAVE
     DONE(0.519537   SEC) : LOCAL POTENTIAL
     DONE(0.565509   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 18.8965
     DONE(0.565583   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.33162    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08166785e+04   0.00000000e+00   6.0630e-01  13.05
     CG2     -1.08159777e+04   7.00830325e-01   7.7004e+00   6.32
     CG3     -1.08206922e+04  -4.71449072e+00   2.3884e+01   3.98
     CG4     -1.08280928e+04  -7.40060143e+00   9.6692e+00   3.45
     CG5     -1.08313193e+04  -3.22649654e+00   2.7254e-01   2.79
     CG6     -1.08314258e+04  -1.06535961e-01   1.7255e-01   2.90
     CG7     -1.08314935e+04  -6.76405591e-02   8.2798e-02   2.87
     CG8     -1.08315608e+04  -6.72771590e-02   2.2931e-02   3.07
     CG9     -1.08316272e+04  -6.64664538e-02   5.4567e-03   3.62
     CG10    -1.08316381e+04  -1.08386549e-02   8.4480e-04   3.32
     CG11    -1.08316415e+04  -3.46294906e-03   9.9482e-04   3.88
     CG12    -1.08316420e+04  -4.77918392e-04   6.9938e-04   2.75
     CG13    -1.08316420e+04  -3.73843226e-05   5.1718e-04   2.78
     CG14    -1.08316413e+04   7.46570830e-04   7.0880e-04   2.74
     CG15    -1.08316400e+04   1.30612625e-03   1.4709e-03   2.74
     CG16    -1.08316399e+04   7.64970022e-05   8.6115e-04   2.79
     CG17    -1.08316404e+04  -4.89203449e-04   1.0133e-03   2.77
     CG18    -1.08316408e+04  -4.09953680e-04   1.3983e-06   2.74
     CG19    -1.08316410e+04  -1.37618882e-04   7.7258e-06   6.24
     CG20    -1.08316410e+04  -1.38520334e-05   6.6774e-07   4.46
     CG21    -1.08316410e+04  -4.32789954e-07   2.6954e-07   2.76
     CG22    -1.08316410e+04   4.50386223e-07   6.7292e-08   2.93
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       86.36  17       5.08   100.00 
     Driver           reading                     0.03   1        0.03   0.04   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 86.33  1        86.33  99.96  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.16   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.05   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    85.77  1        85.77  99.32  
     ESolver_KS_PW    runner                      85.77  1        85.77  99.32  
     ESolver_KS_PW    before_scf                  0.77   1        0.77   0.89   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.10   162      0.00   0.12   
     PW_Basis_Sup     gathers_scatterp            0.04   162      0.00   0.04   
     Potential        init_pot                    0.02   1        0.02   0.02   
     Potential        update_from_charge          0.43   23       0.02   0.50   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.43   23       0.02   0.50   
     H_Hartree_pw     v_hartree                   0.04   23       0.00   0.04   
     PW_Basis_Sup     real2recip                  0.11   206      0.00   0.13   
     PW_Basis_Sup     gatherp_scatters            0.04   206      0.00   0.05   
     PotXC            cal_v_eff                   0.39   23       0.02   0.46   
     XC_Functional    v_xc                        0.39   23       0.02   0.45   
     Potential        interpolate_vrs             0.00   23       0.00   0.00   
     Symmetry         rhog_symmetry               0.16   23       0.01   0.19   
     Symmetry         group fft grids             0.03   23       0.00   0.04   
     PSIInit          initialize_psi              0.73   1        0.73   0.85   
     Nonlocal         getvnl                      0.74   460      0.00   0.86   
     pp_cell_vnl      getvnl                      0.74   460      0.00   0.85   
     Structure_Factor get_sk                      0.09   460      0.00   0.10   
     DiagoIterAssist  diagH_subspace              12.95  440      0.03   14.99  
     Operator         hPsi                        63.80  55811    0.00   73.87  
     Operator         EkineticPW                  0.42   55811    0.00   0.48   
     Operator         VeffPW                      41.72  55811    0.00   48.30  
     PW_Basis_K       recip2real                  24.94  90571    0.00   28.88  
     PW_Basis_K       gathers_scatterp            9.70   90571    0.00   11.24  
     PW_Basis_K       real2recip                  17.36  72971    0.00   20.11  
     PW_Basis_K       gatherp_scatters            4.90   72971    0.00   5.67   
     Operator         NonlocalPW                  21.52  55811    0.00   24.92  
     Nonlocal         add_nonlocal_pp             10.85  55811    0.00   12.57  
     DiagoIterAssist  diagH_LAPACK                0.29   440      0.00   0.34   
     ESolver_KS_PW    hamilt2density_single       84.43  22       3.84   97.77  
     HSolverPW        solve                       84.23  22       3.83   97.53  
     DiagoCG          diag_once                   64.81  440      0.15   75.04  
     DiagoCG_New      spsi_func                   0.50   110742   0.00   0.58   
     DiagoCG_New      hpsi_func                   52.57  55371    0.00   60.87  
     ElecStatePW      psiToRho                    5.91   22       0.27   6.84   
     Charge_Mixing    get_drho                    0.03   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   22       0.00   0.00   
     Charge           mix_rho                     0.05   21       0.00   0.06   
     Charge           Broyden_mixing              0.02   21       0.00   0.02   
     Charge_Mixing    inner_product_recip_hartree 0.02   264      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.08   
     ModuleIO         write_rhog                  0.06   1        0.06   0.07   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:51:41 2025
     FINISH Time  : Sat Apr 26 12:53:07 2025
     TOTAL  Time  : 86
     SEE INFORMATION IN : OUT.Mn_kpoints_6/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:53:08 2025
     MAKE THE DIR         : OUT.Mn_kpoints_7/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0654115  SEC) : SETUP UNITCELL
     DONE(0.207508   SEC) : SYMMETRY
     DONE(0.503179   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       20              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.508696   SEC) : INIT PLANEWAVE
     DONE(0.515375   SEC) : LOCAL POTENTIAL
     DONE(0.561448   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 18.8599
     DONE(0.561517   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.33365    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08159829e+04   0.00000000e+00   6.5125e-01  13.19
     CG2     -1.08231310e+04  -7.14806780e+00   2.9097e+00   6.21
     CG3     -1.08232563e+04  -1.25279501e-01   1.8737e+01   3.51
     CG4     -1.08308272e+04  -7.57091848e+00   1.0041e+00   3.33
     CG5     -1.08312094e+04  -3.82181572e-01   1.2469e-01   2.72
     CG6     -1.08315512e+04  -3.41834126e-01   5.6608e-02   3.28
     CG7     -1.08316025e+04  -5.13303961e-02   3.2475e-02   2.90
     CG8     -1.08316188e+04  -1.62466802e-02   6.1489e-03   2.81
     CG9     -1.08316465e+04  -2.77301216e-02   1.2306e-02   3.94
     CG10    -1.08316468e+04  -3.33933051e-04   5.3277e-03   2.73
     CG11    -1.08316457e+04   1.15532517e-03   2.4950e-03   2.68
     CG12    -1.08316456e+04   7.75102602e-05   8.1246e-04   2.74
     CG13    -1.08316464e+04  -8.28967258e-04   2.5490e-04   3.20
     CG14    -1.08316474e+04  -1.02037009e-03   9.5017e-05   3.70
     CG15    -1.08316475e+04  -7.74164552e-05   1.0019e-04   3.16
     CG16    -1.08316476e+04  -2.78179939e-05   1.8348e-05   2.69
     CG17    -1.08316476e+04  -2.12260302e-05   1.2414e-05   3.26
     CG18    -1.08316476e+04   2.31046547e-05   1.0700e-05   2.81
     CG19    -1.08316475e+04   1.13039309e-05   4.0694e-06   2.86
     CG20    -1.08316476e+04  -2.16727939e-05   1.4729e-06   4.02
     CG21    -1.08316476e+04   3.44689505e-06   5.7736e-07   2.94
     CG22    -1.08316476e+04  -3.35171495e-06   4.7694e-07   4.25
     CG23    -1.08316476e+04   4.19380842e-07   7.8746e-08   2.76
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       87.15  17       5.13   100.00 
     Driver           reading                     0.04   1        0.04   0.04   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 87.11  1        87.11  99.96  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.16   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.05   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    86.56  1        86.56  99.33  
     ESolver_KS_PW    runner                      86.56  1        86.56  99.33  
     ESolver_KS_PW    before_scf                  0.77   1        0.77   0.89   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.10   169      0.00   0.12   
     PW_Basis_Sup     gathers_scatterp            0.04   169      0.00   0.04   
     Potential        init_pot                    0.02   1        0.02   0.02   
     Potential        update_from_charge          0.45   24       0.02   0.52   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.45   24       0.02   0.52   
     H_Hartree_pw     v_hartree                   0.04   24       0.00   0.05   
     PW_Basis_Sup     real2recip                  0.12   215      0.00   0.14   
     PW_Basis_Sup     gatherp_scatters            0.05   215      0.00   0.05   
     PotXC            cal_v_eff                   0.41   24       0.02   0.47   
     XC_Functional    v_xc                        0.41   24       0.02   0.47   
     Potential        interpolate_vrs             0.00   24       0.00   0.00   
     Symmetry         rhog_symmetry               0.17   24       0.01   0.20   
     Symmetry         group fft grids             0.04   24       0.00   0.04   
     PSIInit          initialize_psi              0.73   1        0.73   0.84   
     Nonlocal         getvnl                      0.78   480      0.00   0.90   
     pp_cell_vnl      getvnl                      0.78   480      0.00   0.89   
     Structure_Factor get_sk                      0.08   480      0.00   0.10   
     DiagoIterAssist  diagH_subspace              13.29  460      0.03   15.25  
     Operator         hPsi                        63.95  56508    0.00   73.38  
     Operator         EkineticPW                  0.42   56508    0.00   0.48   
     Operator         VeffPW                      41.68  56508    0.00   47.83  
     PW_Basis_K       recip2real                  24.72  92848    0.00   28.36  
     PW_Basis_K       gathers_scatterp            9.73   92848    0.00   11.17  
     PW_Basis_K       real2recip                  17.62  74448    0.00   20.22  
     PW_Basis_K       gatherp_scatters            5.03   74448    0.00   5.77   
     Operator         NonlocalPW                  21.72  56508    0.00   24.93  
     Nonlocal         add_nonlocal_pp             10.80  56508    0.00   12.40  
     DiagoIterAssist  diagH_LAPACK                0.30   460      0.00   0.34   
     ESolver_KS_PW    hamilt2density_single       85.18  23       3.70   97.75  
     HSolverPW        solve                       84.97  23       3.69   97.50  
     DiagoCG          diag_once                   64.95  460      0.14   74.53  
     DiagoCG_New      spsi_func                   0.51   112096   0.00   0.58   
     DiagoCG_New      hpsi_func                   52.46  56048    0.00   60.20  
     ElecStatePW      psiToRho                    6.09   23       0.26   6.99   
     Charge_Mixing    get_drho                    0.03   23       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   23       0.00   0.00   
     Charge           mix_rho                     0.05   22       0.00   0.06   
     Charge           Broyden_mixing              0.02   22       0.00   0.03   
     Charge_Mixing    inner_product_recip_hartree 0.02   280      0.00   0.02   
     ESolver_KS_PW    after_scf                   0.08   1        0.08   0.09   
     ModuleIO         write_rhog                  0.07   1        0.07   0.08   
     ModuleIO         write_istate_info           0.01   1        0.01   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:53:08 2025
     FINISH Time  : Sat Apr 26 12:54:35 2025
     TOTAL  Time  : 87
     SEE INFORMATION IN : OUT.Mn_kpoints_7/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:54:37 2025
     MAKE THE DIR         : OUT.Mn_kpoints_8/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0670229  SEC) : SETUP UNITCELL
     DONE(0.208006   SEC) : SYMMETRY
     DONE(0.504694   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       35              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.51358    SEC) : INIT PLANEWAVE
     DONE(0.520485   SEC) : LOCAL POTENTIAL
     DONE(0.566188   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 33.0688
     DONE(0.566265   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.87583    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08166330e+04   0.00000000e+00   6.3267e-01  22.91
     CG2     -1.08183169e+04  -1.68392841e+00   5.4734e+00  10.41
     CG3     -1.08212245e+04  -2.90751812e+00   2.4439e+01   6.63
     CG4     -1.08292460e+04  -8.02155522e+00   6.8538e+00   5.96
     CG5     -1.08313668e+04  -2.12076313e+00   2.5243e-01   4.75
     CG6     -1.08314638e+04  -9.70081264e-02   9.3038e-02   4.90
     CG7     -1.08315541e+04  -9.03229771e-02   6.4738e-02   5.56
     CG8     -1.08315838e+04  -2.96371172e-02   1.8611e-02   4.82
     CG9     -1.08316349e+04  -5.11336448e-02   4.1032e-03   6.03
     CG10    -1.08316463e+04  -1.13879483e-02   9.5538e-04   5.97
     CG11    -1.08316505e+04  -4.20200232e-03   1.2421e-03   6.08
     CG12    -1.08316466e+04   3.92103623e-03   1.6943e-03   4.74
     CG13    -1.08316472e+04  -5.99624790e-04   1.8846e-04   4.73
     CG14    -1.08316470e+04   1.45905627e-04   7.2424e-04   6.26
     CG15    -1.08316475e+04  -4.93075945e-04   2.5910e-05   5.28
     CG16    -1.08316476e+04  -1.13417283e-04   5.5517e-06   6.59
     CG17    -1.08316476e+04   1.48751228e-05   1.0694e-05   5.82
     CG18    -1.08316476e+04  -1.57296756e-05   7.4958e-07   5.09
     CG19    -1.08316476e+04  -2.72578058e-06   3.0056e-07   7.07
     CG20    -1.08316476e+04   3.15775383e-07   6.1594e-08   5.06
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       136.72 17       8.04   100.00 
     Driver           reading                     0.04   1        0.04   0.03   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 136.68 1        136.68 99.97  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.10   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    136.12 1        136.12 99.56  
     ESolver_KS_PW    runner                      136.12 1        136.12 99.56  
     ESolver_KS_PW    before_scf                  1.31   1        1.31   0.96   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.01   0.01   
     PW_Basis_Sup     recip2real                  0.10   148      0.00   0.07   
     PW_Basis_Sup     gathers_scatterp            0.04   148      0.00   0.03   
     Potential        init_pot                    0.02   1        0.02   0.01   
     Potential        update_from_charge          0.38   21       0.02   0.28   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.38   21       0.02   0.28   
     H_Hartree_pw     v_hartree                   0.03   21       0.00   0.02   
     PW_Basis_Sup     real2recip                  0.10   188      0.00   0.07   
     PW_Basis_Sup     gatherp_scatters            0.04   188      0.00   0.03   
     PotXC            cal_v_eff                   0.34   21       0.02   0.25   
     XC_Functional    v_xc                        0.34   21       0.02   0.25   
     Potential        interpolate_vrs             0.00   21       0.00   0.00   
     Symmetry         rhog_symmetry               0.14   21       0.01   0.10   
     Symmetry         group fft grids             0.03   21       0.00   0.02   
     PSIInit          initialize_psi              1.27   1        1.27   0.93   
     Nonlocal         getvnl                      1.23   735      0.00   0.90   
     pp_cell_vnl      getvnl                      1.23   735      0.00   0.90   
     Structure_Factor get_sk                      0.13   735      0.00   0.09   
     DiagoIterAssist  diagH_subspace              20.31  700      0.03   14.85  
     Operator         hPsi                        101.15 90004    0.00   73.98  
     Operator         EkineticPW                  0.66   90004    0.00   0.48   
     Operator         VeffPW                      65.44  90004    0.00   47.86  
     PW_Basis_K       recip2real                  38.84  145304   0.00   28.41  
     PW_Basis_K       gathers_scatterp            15.36  145304   0.00   11.24  
     PW_Basis_K       real2recip                  27.45  117304   0.00   20.08  
     PW_Basis_K       gatherp_scatters            8.07   117304   0.00   5.90   
     Operator         NonlocalPW                  34.85  90004    0.00   25.49  
     Nonlocal         add_nonlocal_pp             17.49  90004    0.00   12.80  
     DiagoIterAssist  diagH_LAPACK                0.46   700      0.00   0.34   
     ESolver_KS_PW    hamilt2density_single       134.24 20       6.71   98.18  
     HSolverPW        solve                       134.05 20       6.70   98.05  
     DiagoCG          diag_once                   103.56 700      0.15   75.75  
     DiagoCG_New      spsi_func                   0.81   178608   0.00   0.59   
     DiagoCG_New      hpsi_func                   83.62  89304    0.00   61.16  
     ElecStatePW      psiToRho                    9.33   20       0.47   6.82   
     Charge_Mixing    get_drho                    0.03   20       0.00   0.02   
     Charge_Mixing    inner_product_recip_rho     0.00   20       0.00   0.00   
     Charge           mix_rho                     0.05   19       0.00   0.03   
     Charge           Broyden_mixing              0.02   19       0.00   0.01   
     Charge_Mixing    inner_product_recip_hartree 0.02   232      0.00   0.01   
     ESolver_KS_PW    after_scf                   0.13   1        0.13   0.10   
     ModuleIO         write_rhog                  0.12   1        0.12   0.09   
     ModuleIO         write_istate_info           0.02   1        0.02   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:54:37 2025
     FINISH Time  : Sat Apr 26 12:56:53 2025
     TOTAL  Time  : 136
     SEE INFORMATION IN : OUT.Mn_kpoints_8/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:56:54 2025
     MAKE THE DIR         : OUT.Mn_kpoints_9/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0643955  SEC) : SETUP UNITCELL
     DONE(0.208784   SEC) : SYMMETRY
     DONE(0.508511   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       35              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.518014   SEC) : INIT PLANEWAVE
     DONE(0.524581   SEC) : LOCAL POTENTIAL
     DONE(0.567941   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 32.9834
     DONE(0.568037   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.85836    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08165651e+04   0.00000000e+00   6.1975e-01  23.24
     CG2     -1.08232747e+04  -6.70967029e+00   2.8400e+00  10.68
     CG3     -1.08247024e+04  -1.42765831e+00   1.7210e+01   6.60
     CG4     -1.08298923e+04  -5.18993054e+00   4.9494e+00   5.24
     CG5     -1.08314182e+04  -1.52589123e+00   1.6765e-01   4.74
     CG6     -1.08315117e+04  -9.34737231e-02   5.3121e-02   5.00
     CG7     -1.08315588e+04  -4.70773787e-02   1.5138e-01   5.87
     CG8     -1.08316116e+04  -5.28303222e-02   1.3207e-02   4.78
     CG9     -1.08316290e+04  -1.74223471e-02   9.2069e-03   5.54
     CG10    -1.08316359e+04  -6.85837933e-03   2.5341e-03   4.90
     CG11    -1.08316448e+04  -8.89179610e-03   1.1216e-03   5.89
     CG12    -1.08316424e+04   2.38116915e-03   1.2056e-03   5.21
     CG13    -1.08316420e+04   3.97880579e-04   2.1763e-04   4.77
     CG14    -1.08316432e+04  -1.21723853e-03   7.1908e-05   7.04
     CG15    -1.08316428e+04   4.08611127e-04   3.1392e-04   5.46
     CG16    -1.08316431e+04  -2.80151802e-04   4.3491e-05   5.49
     CG17    -1.08316431e+04  -3.95016846e-05   6.0118e-06   4.70
     CG18    -1.08316432e+04  -1.94783832e-05   5.2354e-06   6.31
     CG19    -1.08316431e+04   1.60574912e-05   5.5432e-06   4.79
     CG20    -1.08316431e+04   1.21484002e-05   2.3210e-06   4.72
     CG21    -1.08316431e+04  -1.27984490e-05   4.3235e-07   6.51
     CG22    -1.08316431e+04   1.19054786e-07   2.1792e-07   5.92
     CG23    -1.08316431e+04  -2.55764618e-07   2.0828e-08   5.60
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       150.97 17       8.88   100.00 
     Driver           reading                     0.04   1        0.04   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 150.93 1        150.93 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.14   1        0.14   0.10   
     PW_Basis_K       setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.03   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    150.37 1        150.37 99.60  
     ESolver_KS_PW    runner                      150.37 1        150.37 99.60  
     ESolver_KS_PW    before_scf                  1.29   1        1.29   0.85   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.01   
     PW_Basis_Sup     recip2real                  0.10   169      0.00   0.07   
     PW_Basis_Sup     gathers_scatterp            0.04   169      0.00   0.02   
     Potential        init_pot                    0.02   1        0.02   0.02   
     Potential        update_from_charge          0.44   24       0.02   0.29   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.44   24       0.02   0.29   
     H_Hartree_pw     v_hartree                   0.04   24       0.00   0.03   
     PW_Basis_Sup     real2recip                  0.14   215      0.00   0.09   
     PW_Basis_Sup     gatherp_scatters            0.06   215      0.00   0.04   
     PotXC            cal_v_eff                   0.40   24       0.02   0.27   
     XC_Functional    v_xc                        0.40   24       0.02   0.27   
     Potential        interpolate_vrs             0.00   24       0.00   0.00   
     Symmetry         rhog_symmetry               0.17   24       0.01   0.11   
     Symmetry         group fft grids             0.03   24       0.00   0.02   
     PSIInit          initialize_psi              1.25   1        1.25   0.83   
     Nonlocal         getvnl                      1.36   840      0.00   0.90   
     pp_cell_vnl      getvnl                      1.36   840      0.00   0.90   
     Structure_Factor get_sk                      0.15   840      0.00   0.10   
     DiagoIterAssist  diagH_subspace              23.04  805      0.03   15.26  
     Operator         hPsi                        111.78 99353    0.00   74.04  
     Operator         EkineticPW                  0.75   99353    0.00   0.50   
     Operator         VeffPW                      72.55  99353    0.00   48.06  
     PW_Basis_K       recip2real                  42.76  162948   0.00   28.33  
     PW_Basis_K       gathers_scatterp            16.78  162948   0.00   11.11  
     PW_Basis_K       real2recip                  30.76  130748   0.00   20.37  
     PW_Basis_K       gatherp_scatters            8.64   130748   0.00   5.72   
     Operator         NonlocalPW                  38.26  99353    0.00   25.34  
     Nonlocal         add_nonlocal_pp             19.16  99353    0.00   12.69  
     DiagoIterAssist  diagH_LAPACK                0.54   805      0.00   0.35   
     ESolver_KS_PW    hamilt2density_single       148.47 23       6.46   98.34  
     HSolverPW        solve                       148.25 23       6.45   98.20  
     DiagoCG          diag_once                   113.66 805      0.14   75.28  
     DiagoCG_New      spsi_func                   0.86   197096   0.00   0.57   
     DiagoCG_New      hpsi_func                   91.89  98548    0.00   60.87  
     ElecStatePW      psiToRho                    10.43  23       0.45   6.91   
     Charge_Mixing    get_drho                    0.05   23       0.00   0.03   
     Charge_Mixing    inner_product_recip_rho     0.00   23       0.00   0.00   
     Charge           mix_rho                     0.05   22       0.00   0.03   
     Charge           Broyden_mixing              0.02   22       0.00   0.01   
     Charge_Mixing    inner_product_recip_hartree 0.02   280      0.00   0.01   
     ESolver_KS_PW    after_scf                   0.08   1        0.08   0.05   
     ModuleIO         write_rhog                  0.06   1        0.06   0.04   
     ModuleIO         write_istate_info           0.02   1        0.02   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:56:54 2025
     FINISH Time  : Sat Apr 26 12:59:25 2025
     TOTAL  Time  : 151
     SEE INFORMATION IN : OUT.Mn_kpoints_9/
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 12:59:26 2025
     MAKE THE DIR         : OUT.Mn_kpoints_10/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 36 * 36 * 36
     UNIFORM GRID DIM(BIG)   : 36 * 36 * 36
     DONE(0.0670745  SEC) : SETUP UNITCELL
     DONE(0.22563    SEC) : SYMMETRY
     DONE(0.560802   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Self-consistent calculations for electrons
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       56              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      4           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.575317   SEC) : INIT PLANEWAVE
     DONE(0.581894   SEC) : LOCAL POTENTIAL
     DONE(0.628295   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 52.9102
     DONE(0.628374   SEC) : INIT BASIS
     -------------------------------------------
     SELF-CONSISTENT : 
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.6985     SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -1.08167866e+04   0.00000000e+00   6.0593e-01  36.68
     CG2     -1.08198087e+04  -3.02212052e+00   4.5312e+00  17.59
     CG3     -1.08233077e+04  -3.49894408e+00   1.9338e+01  11.06
     CG4     -1.08295987e+04  -6.29099309e+00   5.9707e+00   9.47
     CG5     -1.08313705e+04  -1.77179831e+00   3.1017e-01   7.77
     CG6     -1.08314437e+04  -7.32686708e-02   9.5129e-02   7.62
     CG7     -1.08315108e+04  -6.70225518e-02   1.7292e-01   9.03
     CG8     -1.08315813e+04  -7.05015810e-02   1.2738e-02   7.71
     CG9     -1.08316385e+04  -5.72790227e-02   4.1022e-03  10.91
     CG10    -1.08316414e+04  -2.89818285e-03   4.2335e-03   8.38
     CG11    -1.08316425e+04  -1.04426514e-03   6.2833e-04   7.62
     CG12    -1.08316438e+04  -1.36832813e-03   5.5389e-06   9.68
     CG13    -1.08316443e+04  -4.87560118e-04   3.2941e-06  18.21
     CG14    -1.08316443e+04   1.09045020e-05   4.0592e-05   7.89
     CG15    -1.08316443e+04  -1.54924075e-05   2.6106e-06   7.61
     CG16    -1.08316443e+04  -2.33198053e-06   3.7952e-07   7.68
     CG17    -1.08316443e+04  -4.72701745e-07   1.2330e-06  10.69
     CG18    -1.08316443e+04  -4.57871134e-07   1.7463e-07   7.62
     CG19    -1.08316443e+04   4.73793778e-07   5.8480e-08   7.72
    TIME STATISTICS
    ----------------------------------------------------------------------------
        CLASS_NAME               NAME             TIME/s  CALLS   AVG/s  PER/%  
    ----------------------------------------------------------------------------
                      total                       213.80 17       12.58  100.00 
     Driver           reading                     0.04   1        0.04   0.02   
     Input_Conv       Convert                     0.00   1        0.00   0.00   
     Driver           driver_line                 213.76 1        213.76 99.98  
     UnitCell         check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup     setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup     distributeg                 0.00   1        0.00   0.00   
     mymath           heapsort                    0.00   453      0.00   0.00   
     Charge_Mixing    init_mixing                 0.00   2        0.00   0.00   
     Symmetry         analy_sys                   0.16   1        0.16   0.07   
     PW_Basis_K       setuptransform              0.01   1        0.01   0.00   
     PW_Basis_K       distributeg                 0.00   1        0.00   0.00   
     PW_Basis         setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl        init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl       init                        0.00   1        0.00   0.00   
     ppcell_vnl       init_vnl                    0.04   1        0.04   0.02   
     WF_atomic        init_at_1                   0.00   1        0.00   0.00   
     wavefunc         wfcinit                     0.00   1        0.00   0.00   
     Ions             opt_ions                    213.09 1        213.09 99.67  
     ESolver_KS_PW    runner                      213.09 1        213.09 99.67  
     ESolver_KS_PW    before_scf                  2.07   1        2.07   0.97   
     H_Ewald_pw       compute_ewald               0.00   1        0.00   0.00   
     Charge           set_rho_core                0.00   1        0.00   0.00   
     Charge           atomic_rho                  0.01   2        0.00   0.00   
     PW_Basis_Sup     recip2real                  0.09   141      0.00   0.04   
     PW_Basis_Sup     gathers_scatterp            0.04   141      0.00   0.02   
     Potential        init_pot                    0.02   1        0.02   0.01   
     Potential        update_from_charge          0.38   20       0.02   0.18   
     Potential        cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal         cal_fixed_v                 0.00   1        0.00   0.00   
     Potential        cal_v_eff                   0.38   20       0.02   0.18   
     H_Hartree_pw     v_hartree                   0.03   20       0.00   0.01   
     PW_Basis_Sup     real2recip                  0.10   179      0.00   0.05   
     PW_Basis_Sup     gatherp_scatters            0.04   179      0.00   0.02   
     PotXC            cal_v_eff                   0.35   20       0.02   0.16   
     XC_Functional    v_xc                        0.35   20       0.02   0.16   
     Potential        interpolate_vrs             0.00   20       0.00   0.00   
     Symmetry         rhog_symmetry               0.14   20       0.01   0.06   
     Symmetry         group fft grids             0.03   20       0.00   0.01   
     PSIInit          initialize_psi              2.03   1        2.03   0.95   
     Nonlocal         getvnl                      1.81   1120     0.00   0.84   
     pp_cell_vnl      getvnl                      1.80   1120     0.00   0.84   
     Structure_Factor get_sk                      0.20   1120     0.00   0.09   
     DiagoIterAssist  diagH_subspace              31.08  1064     0.03   14.54  
     Operator         hPsi                        159.09 141737   0.00   74.41  
     Operator         EkineticPW                  1.05   141737   0.00   0.49   
     Operator         VeffPW                      103.08 141737   0.00   48.22  
     PW_Basis_K       recip2real                  60.29  225793   0.00   28.20  
     PW_Basis_K       gathers_scatterp            23.90  225793   0.00   11.18  
     PW_Basis_K       real2recip                  43.62  183233   0.00   20.40  
     PW_Basis_K       gatherp_scatters            12.57  183233   0.00   5.88   
     Operator         NonlocalPW                  54.64  141737   0.00   25.56  
     Nonlocal         add_nonlocal_pp             27.49  141737   0.00   12.86  
     DiagoIterAssist  diagH_LAPACK                0.73   1064     0.00   0.34   
     ESolver_KS_PW    hamilt2density_single       210.51 19       11.08  98.46  
     HSolverPW        solve                       210.33 19       11.07  98.38  
     DiagoCG          diag_once                   163.88 1064     0.15   76.65  
     DiagoCG_New      spsi_func                   1.25   281346   0.00   0.59   
     DiagoCG_New      hpsi_func                   132.24 140673   0.00   61.85  
     ElecStatePW      psiToRho                    14.20  19       0.75   6.64   
     Charge_Mixing    get_drho                    0.02   19       0.00   0.01   
     Charge_Mixing    inner_product_recip_rho     0.00   19       0.00   0.00   
     Charge           mix_rho                     0.05   18       0.00   0.02   
     Charge           Broyden_mixing              0.02   18       0.00   0.01   
     Charge_Mixing    inner_product_recip_hartree 0.01   216      0.00   0.01   
     ESolver_KS_PW    after_scf                   0.07   1        0.07   0.03   
     ModuleIO         write_rhog                  0.06   1        0.06   0.03   
     ModuleIO         write_istate_info           0.02   1        0.02   0.01   
    ----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 12:59:26 2025
     FINISH Time  : Sat Apr 26 13:03:00 2025
     TOTAL  Time  : 214
     SEE INFORMATION IN : OUT.Mn_kpoints_10/
    
    统计结果已保存到 ./Mn/final_etot_vs_kpoints.csv
    Figure(800x600)
    

### 输出图片
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/475862/e010d61ade7245c499cac194228dc856/ef69d74b-83ac-488c-88f9-c623f86e70d3.png)

# 结构弛豫
将<https://legacy.materialsproject.org/materials/mp-8634/#>中的晶格参数做一个小偏移，输入ABACUS进行计算

网站给出晶格参数为
```
_cell_length_a   2.482
_cell_length_b   2.482
_cell_length_c   2.482
_cell_angle_alpha  60.00000000
_cell_angle_beta  60.00000000
_cell_angle_gamma  60.00000000
```
设置STRU文件为
```
ATOMIC_SPECIES
Mn 54.9380 Mn_ONCV_PBE-1.2.upf        # 名称; 相对原子质量; 赝势文件名

NUMERICAL_ORBITAL
Mn_gga_10au_100Ry_2s1p1d.orb         # 轨道文件名

LATTICE_CONSTANT
1.8897259886                         # 晶格常数的单位 (Bohr), 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
1.700 1.700 0.000
1.700 0.000 1.700
0.000 1.700 1.700

ATOMIC_POSITIONS
Direct                             # 笛卡尔坐标或直接坐标
Mn                                 # 元素名称
0.0                                # 磁矩
1                                  # 原子个数
0.000 0.000 0.000 0 0 0            # 原子位置和自由度
```
根据收敛性测试结果，取kpoints=6，ecutwfc=60


```
%cd ../cell_realx/Mn/
!mpirun -np 2 abacus
```

    /opt/mamba/lib/python3.10/site-packages/IPython/core/magics/osm.py:393: UserWarning: This is now an optional IPython functionality, using bookmarks requires you to install the `pickleshare` library.
      bkms = self.shell.db.get('bookmarks', {})
    [Errno 2] No such file or directory: '../cell_realx/Mn/'
    /personal/cell_realx/Mn
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:12:52 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 24 * 24 * 24
     UNIFORM GRID DIM(BIG)   : 24 * 24 * 24
     DONE(0.0828629  SEC) : SETUP UNITCELL
     DONE(0.15694    SEC) : SYMMETRY
     DONE(0.449675   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Cell relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       16              2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.451045   SEC) : INIT PLANEWAVE
     DONE(0.452514   SEC) : LOCAL POTENTIAL
     DONE(0.524099   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 1.2041
     DONE(0.524162   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.599013   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777139e+03   0.00000000e+00   8.9373e-02   1.20
     CG2     -2.70832177e+03  -5.50384589e-01   5.2569e-02   0.23
     CG3     -2.70845356e+03  -1.31790256e-01   2.9681e-04   0.20
     CG4     -2.70846373e+03  -1.01665129e-02   9.9074e-05   0.45
     CG5     -2.70846427e+03  -5.45388212e-04   3.9052e-06   0.22
     CG6     -2.70846430e+03  -3.16680575e-05   7.6166e-08   0.35
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           307.2561366273        -0.0000000000         0.0000000000 
            -0.0000000000       307.2561366273         0.0000000000 
             0.0000000000         0.0000000000       307.2561366273 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 307.256137 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(3.347575   SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 2
     -------------------------------------------
     DONE(3.372930   SEC) : LOCAL POTENTIAL
     DONE(3.441230   SEC) : SYMMETRY
     DONE(3.441513   SEC) : INIT K-POINTS
     DONE(3.508896   SEC) : NON-LOCAL POTENTIAL
     DONE(3.519254   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70759766e+03   0.00000000e+00   4.1081e-04   0.51
     CG2     -2.70820635e+03  -6.08691380e-01   5.2023e-04   0.33
     CG3     -2.70820819e+03  -1.84063005e-03   6.2583e-04   0.25
     CG4     -2.70820969e+03  -1.49640914e-03   9.4359e-07   0.22
     CG5     -2.70820973e+03  -4.60840257e-05   4.7684e-07   0.39
     CG6     -2.70820973e+03  -7.70102437e-07   2.4864e-08   0.20
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
          -456.9016277594        -0.0000000000         0.0000000000 
            -0.0000000000      -456.9016277594         0.0000000000 
            -0.0000000000         0.0000000000      -456.9016277594 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -456.901628 KBAR
    
     ETOT DIFF (eV)       : 0.254572
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(5.516353   SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 3
     -------------------------------------------
     DONE(5.536941   SEC) : LOCAL POTENTIAL
     DONE(5.602479   SEC) : SYMMETRY
     DONE(5.602659   SEC) : INIT K-POINTS
     DONE(5.671687   SEC) : NON-LOCAL POTENTIAL
     DONE(5.683002   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70853857e+03   0.00000000e+00   1.9770e-03   0.58
     CG2     -2.70854089e+03  -2.32533285e-03   3.0151e-04   0.20
     CG3     -2.70854247e+03  -1.57699595e-03   8.4623e-05   0.25
     CG4     -2.70854272e+03  -2.56001347e-04   7.9111e-07   0.23
     CG5     -2.70854276e+03  -3.57805397e-05   5.1348e-08   0.39
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           -19.3687430452        -0.0000000000         0.0000000000 
            -0.0000000000       -19.3687430452         0.0000000000 
            -0.0000000000         0.0000000000       -19.3687430452 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -19.368743 KBAR
    
     ETOT DIFF (eV)       : -0.333028
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(7.455805   SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 4
     -------------------------------------------
     DONE(7.476328   SEC) : LOCAL POTENTIAL
     DONE(7.548210   SEC) : SYMMETRY
     DONE(7.548471   SEC) : INIT K-POINTS
     DONE(7.618281   SEC) : NON-LOCAL POTENTIAL
     DONE(7.630786   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854153e+03   0.00000000e+00   4.4278e-04   0.39
     CG2     -2.70854218e+03  -6.46320704e-04   1.4326e-04   0.20
     CG3     -2.70854273e+03  -5.49510093e-04   1.0384e-05   0.20
     CG4     -2.70854276e+03  -3.01596399e-05   2.0408e-07   0.23
     CG5     -2.70854276e+03  -4.11218631e-06   4.0511e-08   0.31
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           -20.4101212769        -0.0000000000         0.0000000000 
            -0.0000000000       -20.4101212769         0.0000000000 
            -0.0000000000        -0.0000000000       -20.4101212769 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -20.410121 KBAR
    
     ETOT DIFF (eV)       : -0.000001
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(9.082641   SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 5
     -------------------------------------------
     DONE(9.102454   SEC) : LOCAL POTENTIAL
     DONE(9.169622   SEC) : SYMMETRY
     DONE(9.169911   SEC) : INIT K-POINTS
     DONE(9.238118   SEC) : NON-LOCAL POTENTIAL
     DONE(9.247405   SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70892142e+03   0.00000000e+00   8.0934e-06   0.57
     CG2     -2.70854186e+03   3.79554565e-01   5.4434e-05   0.40
     CG3     -2.70854208e+03  -2.19198002e-04   4.2548e-05   0.33
     CG4     -2.70854219e+03  -1.02601697e-04   1.0175e-07   0.25
     CG5     -2.70854219e+03  -3.52650479e-06   2.2343e-08   0.34
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            30.8834924581         0.0000000000        -0.0000000000 
             0.0000000000        30.8834924581        -0.0000000000 
            -0.0000000000        -0.0000000000        30.8834924581 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 30.883492 KBAR
    
     ETOT DIFF (eV)       : 0.000571
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(11.236804  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 6
     -------------------------------------------
     DONE(11.257319  SEC) : LOCAL POTENTIAL
     DONE(11.326649  SEC) : SYMMETRY
     DONE(11.326942  SEC) : INIT K-POINTS
     DONE(11.395406  SEC) : NON-LOCAL POTENTIAL
     DONE(11.404358  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70847576e+03   0.00000000e+00   3.3982e-06   0.54
     CG2     -2.70854328e+03  -6.75178453e-02   6.4957e-06   0.33
     CG3     -2.70854331e+03  -2.56522758e-05   4.5704e-06   0.32
     CG4     -2.70854332e+03  -1.14301225e-05   6.1830e-09   0.20
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -2.3966684106        -0.0000000000         0.0000000000 
            -0.0000000000        -2.3966684106         0.0000000000 
             0.0000000000         0.0000000000        -2.3966684106 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -2.396668 KBAR
    
     ETOT DIFF (eV)       : -0.001127
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(12.891239  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 7
     -------------------------------------------
     DONE(12.911182  SEC) : LOCAL POTENTIAL
     DONE(12.983335  SEC) : SYMMETRY
     DONE(12.983598  SEC) : INIT K-POINTS
     DONE(13.052088  SEC) : NON-LOCAL POTENTIAL
     DONE(13.061244  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854332e+03   0.00000000e+00   6.4679e-07   0.44
     CG2     -2.70854332e+03  -9.36046578e-07   7.0637e-09   0.21
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             1.9429793075        -0.0000000000         0.0000000000 
            -0.0000000000         1.9429793075         0.0000000000 
            -0.0000000000        -0.0000000000         1.9429793075 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 1.942979 KBAR
    
     ETOT DIFF (eV)       : -0.000000
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(13.809429  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 8
     -------------------------------------------
     DONE(13.834842  SEC) : LOCAL POTENTIAL
     DONE(13.903775  SEC) : SYMMETRY
     DONE(13.904054  SEC) : INIT K-POINTS
     DONE(13.972694  SEC) : NON-LOCAL POTENTIAL
     DONE(13.981825  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854332e+03   0.00000000e+00   1.4718e-07   0.40
     CG2     -2.70854332e+03  -8.12099686e-07   1.2505e-07   0.26
     CG3     -2.70854332e+03  -2.32647860e-07   1.8716e-09   0.21
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             0.8232798932         0.0000000000        -0.0000000000 
             0.0000000000         0.8232798932        -0.0000000000 
             0.0000000000         0.0000000000         0.8232798932 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 0.823280 KBAR
    
     ETOT DIFF (eV)       : -0.000000
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(14.949644  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 9
     -------------------------------------------
     DONE(14.970365  SEC) : LOCAL POTENTIAL
     DONE(15.038850  SEC) : SYMMETRY
     DONE(15.039140  SEC) : INIT K-POINTS
     DONE(15.108133  SEC) : NON-LOCAL POTENTIAL
     DONE(15.117156  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854332e+03   0.00000000e+00   1.4388e-08   0.49
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -1.6221457371         0.0000000000        -0.0000000000 
             0.0000000000        -1.6221457371        -0.0000000000 
             0.0000000000        -0.0000000000        -1.6221457371 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -1.622146 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(15.705965  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 10
     -------------------------------------------
     DONE(15.729451  SEC) : LOCAL POTENTIAL
     DONE(15.797987  SEC) : SYMMETRY
     DONE(15.798276  SEC) : INIT K-POINTS
     DONE(15.867137  SEC) : NON-LOCAL POTENTIAL
     DONE(15.876673  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854332e+03   0.00000000e+00   2.8044e-07   0.39
     CG2     -2.70854332e+03  -1.32765869e-06   2.3172e-07   0.24
     CG3     -2.70854332e+03  -4.84296860e-07   2.3895e-09   0.20
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -0.7636457021        -0.0000000000         0.0000000000 
            -0.0000000000        -0.7636457021         0.0000000000 
             0.0000000000         0.0000000000        -0.7636457021 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -0.763646 KBAR
    
     ETOT DIFF (eV)       : -0.000000
     LARGEST GRAD (eV/A)  : 0.000000
     DONE(16.810766  SEC) : SETUP UNITCELL
     -------------------------------------------
     STEP OF RELAXATION : 11
     -------------------------------------------
     DONE(16.831827  SEC) : LOCAL POTENTIAL
     DONE(16.903282  SEC) : SYMMETRY
     DONE(16.903504  SEC) : INIT K-POINTS
     DONE(16.986392  SEC) : NON-LOCAL POTENTIAL
     DONE(16.995978  SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70854332e+03   0.00000000e+00   9.1977e-10   0.58
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -0.0465403970        -0.0000000000         0.0000000000 
            -0.0000000000        -0.0465403970         0.0000000000 
             0.0000000000         0.0000000000        -0.0465403970 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -0.046540 KBAR
    
     ETOT DIFF (eV)       : -0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       17.71  137      0.13   100.00 
     Driver            reading                     0.06   1        0.06   0.34   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 17.65  1        17.65  99.66  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   13       0.00   0.01   
     Charge_Mixing     init_mixing                 0.00   12       0.00   0.00   
     Symmetry          analy_sys                   0.77   11       0.07   4.32   
     PW_Basis_K        setuptransform              0.00   1        0.00   0.00   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   11       0.00   0.01   
     ppcell_vl         init_vloc                   0.01   11       0.00   0.06   
     ppcell_vnl        init                        0.00   1        0.00   0.01   
     ppcell_vnl        init_vnl                    0.77   11       0.07   4.34   
     WF_atomic         init_at_1                   0.00   11       0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    17.17  1        17.17  96.92  
     ESolver_KS_PW     runner                      16.34  11       1.49   92.26  
     ESolver_KS_PW     before_scf                  1.58   11       0.14   8.93   
     H_Ewald_pw        compute_ewald               0.00   11       0.00   0.01   
     Charge            set_rho_core                0.00   11       0.00   0.00   
     Charge            atomic_rho                  0.03   22       0.00   0.18   
     PW_Basis_Sup      recip2real                  0.07   416      0.00   0.41   
     PW_Basis_Sup      gathers_scatterp            0.04   416      0.00   0.20   
     Potential         init_pot                    0.06   11       0.01   0.35   
     Potential         update_from_charge          0.28   52       0.01   1.60   
     Potential         cal_fixed_v                 0.00   11       0.00   0.01   
     PotLocal          cal_fixed_v                 0.00   11       0.00   0.01   
     Potential         cal_v_eff                   0.28   52       0.01   1.58   
     H_Hartree_pw      v_hartree                   0.02   52       0.00   0.12   
     PW_Basis_Sup      real2recip                  0.08   556      0.00   0.48   
     PW_Basis_Sup      gatherp_scatters            0.04   556      0.00   0.22   
     PotXC             cal_v_eff                   0.26   52       0.00   1.45   
     XC_Functional     v_xc                        0.26   52       0.00   1.45   
     Potential         interpolate_vrs             0.00   52       0.00   0.00   
     Symmetry          rhog_symmetry               0.08   63       0.00   0.44   
     Symmetry          group fft grids             0.03   63       0.00   0.15   
     PSIInit           initialize_psi              0.07   11       0.01   0.37   
     Nonlocal          getvnl                      0.26   848      0.00   1.48   
     pp_cell_vnl       getvnl                      0.26   848      0.00   1.48   
     Structure_Factor  get_sk                      0.02   1200     0.00   0.10   
     DiagoIterAssist   diagH_subspace              2.42   784      0.00   13.67  
     Operator          hPsi                        10.46  48254    0.00   59.03  
     Operator          EkineticPW                  0.09   48254    0.00   0.51   
     Operator          VeffPW                      9.15   48254    0.00   51.68  
     PW_Basis_K        recip2real                  5.48   76558    0.00   30.93  
     PW_Basis_K        gathers_scatterp            2.80   76558    0.00   15.80  
     PW_Basis_K        real2recip                  3.65   61582    0.00   20.62  
     PW_Basis_K        gatherp_scatters            1.60   61582    0.00   9.02   
     Operator          NonlocalPW                  1.13   48254    0.00   6.38   
     Nonlocal          add_nonlocal_pp             0.44   48254    0.00   2.48   
     DiagoIterAssist   diagH_LAPACK                0.14   784      0.00   0.77   
     ESolver_KS_PW     hamilt2density_single       13.99  52       0.27   78.96  
     HSolverPW         solve                       13.88  52       0.27   78.37  
     DiagoCG           diag_once                   10.03  832      0.01   56.61  
     DiagoCG_New       spsi_func                   0.14   94940    0.00   0.77   
     DiagoCG_New       hpsi_func                   8.37   47470    0.00   47.24  
     ElecStatePW       psiToRho                    1.15   52       0.02   6.47   
     Charge_Mixing     get_drho                    0.02   52       0.00   0.09   
     Charge_Mixing     inner_product_recip_rho     0.00   52       0.00   0.01   
     Charge            mix_rho                     0.01   27       0.00   0.05   
     Charge            Broyden_mixing              0.00   27       0.00   0.01   
     Charge_Mixing     inner_product_recip_hartree 0.00   68       0.00   0.01   
     ESolver_KS_PW     after_scf                   0.51   11       0.05   2.89   
     ModuleIO          write_rhog                  0.48   11       0.04   2.70   
     Forces            cal_force                   0.11   11       0.01   0.64   
     Forces            cal_force_loc               0.01   11       0.00   0.04   
     Forces            cal_force_ew                0.00   11       0.00   0.01   
     Forces            cal_force_nl                0.09   11       0.01   0.49   
     FS_Nonlocal_tools cal_becp                    0.23   176      0.00   1.27   
     Forces            cal_force_cc                0.00   11       0.00   0.00   
     Forces            cal_force_scc               0.02   11       0.00   0.11   
     Stress_PW         cal_stress                  0.48   11       0.04   2.74   
     Stress_Func       stress_kin                  0.01   11       0.00   0.08   
     Stress_Func       stress_har                  0.00   11       0.00   0.02   
     Stress_Func       stress_ewa                  0.00   11       0.00   0.02   
     Stress_Func       stress_gga                  0.03   11       0.00   0.18   
     Stress_Func       stress_loc                  0.04   11       0.00   0.23   
     Stress_Func       stress_cc                   0.00   11       0.00   0.00   
     Stress_Func       stress_nl                   0.39   11       0.04   2.20   
     Charge_Extra      extrapolate_charge          0.02   10       0.00   0.09   
     ModuleIO          write_istate_info           0.01   1        0.01   0.07   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:12:52 2025
     FINISH Time  : Sat Apr 26 13:13:09 2025
     TOTAL  Time  : 17
     SEE INFORMATION IN : OUT.Mn/
    


```
#计算结果
!cat OUT.Mn/STRU_NOW.cif
```

    # Generated by ABACUS ModuleIO::CifParser
    data_?
    _symmetry_space_group_name_H-M   'P 1'
    _cell_length_a   2.47442649
    _cell_length_b   2.47442649
    _cell_length_c   2.47442649
    _cell_angle_alpha  60.00000000
    _cell_angle_beta  60.00000000
    _cell_angle_gamma  60.00000000
    _symmetry_Int_Tables_number   1
    _chemical_formula_structural   Mn
    _chemical_formula_sum   Mn
    _cell_volume    10.71293998
    _cell_formula_units_Z   1
    loop_
     _symmetry_equiv_pos_site_id
     _symmetry_equiv_pos_as_xyz
      1  'x, y, z'
    loop_
     _atom_site_type_symbol
     _atom_site_label
     _atom_site_symmetry_multiplicity
     _atom_site_fract_x
     _atom_site_fract_y
     _atom_site_fract_z
     _atom_site_occupancy
     Mn  Mn0   1   0.00000000   0.00000000   0.00000000 1.0
    
    

弛豫结果与网站给出结果较为接近，误差量级为0.01 angstrom

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/475862/ac0af3a696b345a2960c252c944204ca/8b02806a-99d4-44af-a5a7-c717baa3b765.png)

# 弹性模量计算
参考教程文档中[ABACUS+pymatgen 计算弹性常数](https://mcresearch.github.io/abacus-user-guide/abacus-elastic.html)，将晶格常数设置为上个部分弛豫结果。INPUT文件计算类型改为relax，并添加
```
cal_stress   1
cal_force    1
gamma_only       0
smearing_method  gaussian
smearing_sigma   0.002
mixing_type      broyden
mixing_beta      0.7
```



```
#首先进行结构弛豫
%cd ../../elasticity/Mn
!mpirun -np 2 abacus
```

                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:54:03 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 24 * 24 * 24
     UNIFORM GRID DIM(BIG)   : 24 * 24 * 24
     DONE(0.0739865  SEC) : SETUP UNITCELL
    WARNING: Symmetry cannot be kept when not all atoms are movable.
     Continue with symmetry=0 ... 
     DONE(0.147609   SEC) : SYMMETRY
     DONE(0.162015   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       112             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.167538   SEC) : INIT PLANEWAVE
     DONE(0.169136   SEC) : LOCAL POTENTIAL
     DONE(0.214089   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 9.22852
     DONE(0.214159   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(0.66708    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70769755e+03   0.00000000e+00   9.5752e-02   8.14
     CG2     -2.70841883e+03  -7.21285536e-01   3.8305e-02   1.56
     CG3     -2.70851775e+03  -9.89166180e-02   5.6837e-04   1.46
     CG4     -2.70853104e+03  -1.32946488e-02   2.4995e-04   2.84
     CG5     -2.70853200e+03  -9.58237618e-04   1.2271e-05   1.46
     CG6     -2.70853204e+03  -3.59810584e-05   7.6494e-06   2.32
     CG7     -2.70853205e+03  -9.48615010e-06   6.2281e-07   1.42
     CG8     -2.70853205e+03  -2.79322521e-06   3.0746e-08   2.27
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             5.4751821203        -0.0045050402        -0.0199059453 
            -0.0045050402         5.4161701747         0.0238729395 
            -0.0199059453         0.0238729395         5.5793122921 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 5.490222 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       22.56  17       1.33   100.00 
     Driver            reading                     0.05   1        0.05   0.22   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 22.51  1        22.51  99.78  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.01   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.01   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.33   
     PW_Basis_K        setuptransform              0.00   1        0.00   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.00   1        0.00   0.00   
     ppcell_vnl        init                        0.00   1        0.00   0.01   
     ppcell_vnl        init_vnl                    0.04   1        0.04   0.19   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    22.32  1        22.32  98.91  
     ESolver_KS_PW     runner                      21.97  1        21.97  97.36  
     ESolver_KS_PW     before_scf                  0.45   1        0.45   2.01   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.00   2        0.00   0.01   
     PW_Basis_Sup      recip2real                  0.01   58       0.00   0.04   
     PW_Basis_Sup      gathers_scatterp            0.00   58       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.05   9        0.01   0.20   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.05   9        0.01   0.20   
     H_Hartree_pw      v_hartree                   0.00   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.01   79       0.00   0.05   
     PW_Basis_Sup      gatherp_scatters            0.00   79       0.00   0.02   
     PotXC             cal_v_eff                   0.04   9        0.00   0.18   
     XC_Functional     v_xc                        0.04   9        0.00   0.18   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     PSIInit           initialize_psi              0.44   1        0.44   1.97   
     Nonlocal          getvnl                      0.37   1120     0.00   1.62   
     pp_cell_vnl       getvnl                      0.36   1120     0.00   1.62   
     Structure_Factor  get_sk                      0.02   1344     0.00   0.08   
     DiagoIterAssist   diagH_subspace              2.78   896      0.00   12.32  
     Operator          hPsi                        16.54  80853    0.00   73.30  
     Operator          EkineticPW                  0.18   80853    0.00   0.78   
     Operator          VeffPW                      14.33  80853    0.00   63.52  
     PW_Basis_K        recip2real                  8.33   114229   0.00   36.90  
     PW_Basis_K        gathers_scatterp            4.18   114229   0.00   18.53  
     PW_Basis_K        real2recip                  5.63   96085    0.00   24.95  
     PW_Basis_K        gatherp_scatters            2.42   96085    0.00   10.73  
     Operator          NonlocalPW                  1.90   80853    0.00   8.41   
     Nonlocal          add_nonlocal_pp             0.75   80853    0.00   3.33   
     DiagoIterAssist   diagH_LAPACK                0.16   896      0.00   0.70   
     ESolver_KS_PW     hamilt2density_single       21.40  9        2.38   94.86  
     HSolverPW         solve                       21.40  9        2.38   94.86  
     DiagoCG           diag_once                   17.16  1008     0.02   76.04  
     DiagoCG_New       spsi_func                   0.24   159914   0.00   1.07   
     DiagoCG_New       hpsi_func                   14.16  79957    0.00   62.75  
     ElecStatePW       psiToRho                    1.38   9        0.15   6.11   
     Charge_Mixing     get_drho                    0.00   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.00   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.28   
     ModuleIO          write_rhog                  0.05   1        0.05   0.23   
     Forces            cal_force                   0.06   1        0.06   0.27   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.06   1        0.06   0.27   
     FS_Nonlocal_tools cal_becp                    0.16   112      0.00   0.69   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.00   1        0.00   0.01   
     Stress_PW         cal_stress                  0.27   1        0.27   1.18   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.00   1        0.00   0.01   
     Stress_Func       stress_loc                  0.00   1        0.00   0.02   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.25   1        0.25   1.11   
     ModuleIO          write_istate_info           0.02   1        0.02   0.10   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:54:03 2025
     FINISH Time  : Sat Apr 26 13:54:26 2025
     TOTAL  Time  : 23
     SEE INFORMATION IN : OUT.Mn/
    


```
#产生应变结构
%cd ..
!python gene_dfm.py abacus

```

    /personal/elasticity
    /opt/mamba/lib/python3.10/site-packages/pymatgen/core/structure.py:3175: EncodingWarning: We strongly encourage explicit `encoding`, and we would use UTF-8 by default as per PEP 686
      with zopen(filename, mode="rt", errors="replace") as file:
    gen with norm [-0.01, -0.005, 0.005, 0.01]
    gen with shear [-0.01, -0.005, 0.005, 0.01]
    /opt/mamba/lib/python3.10/site-packages/pymatgen/io/vasp/inputs.py:660: EncodingWarning: We strongly encourage explicit `encoding`, and we would use UTF-8 by default as per PEP 686
      with zopen(filename, mode="wt") as file:
    


```
#计算应力
!sh run_task.sh
```

    /personal/elasticity/task.000
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:55:47 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0718281  SEC) : SETUP UNITCELL
     DONE(0.148901   SEC) : SYMMETRY
     DONE(0.460564   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.47116    SEC) : INIT PLANEWAVE
     DONE(0.478963   SEC) : LOCAL POTENTIAL
     DONE(0.535851   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.535919   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.43712    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70774764e+03   0.00000000e+00   8.8636e-02  15.87
     CG2     -2.70842871e+03  -6.81072499e-01   2.6892e-02   2.98
     CG3     -2.70850905e+03  -8.03334595e-02   2.3222e-04   3.08
     CG4     -2.70852034e+03  -1.12915203e-02   1.4838e-05   5.73
     CG5     -2.70852049e+03  -1.45880457e-04   4.5368e-06   4.15
     CG6     -2.70852050e+03  -1.00931205e-05   2.7811e-07   3.00
     CG7     -2.70852050e+03  -1.76922667e-06   1.0107e-07   4.64
     CG8     -2.70852050e+03  -2.36782040e-07   3.0313e-08   3.18
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            49.3364439604         5.0881121499        -7.2322966682 
             5.0881121499        43.5399299341        -4.2036870519 
            -7.2322966682        -4.2036870519        46.5576902560 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 46.478021 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       45.06  17       2.65   100.00 
     Driver            reading                     0.05   1        0.05   0.10   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 45.02  1        45.02  99.90  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.25  1        44.25  98.21  
     ESolver_KS_PW     runner                      43.61  1        43.61  96.77  
     ESolver_KS_PW     before_scf                  0.90   1        0.90   2.00   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.03   68       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.21   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   9        0.01   0.21   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.08   9        0.01   0.19   
     XC_Functional     v_xc                        0.08   9        0.01   0.19   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.06   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.87   1        0.87   1.93   
     Nonlocal          getvnl                      0.63   1020     0.00   1.41   
     pp_cell_vnl       getvnl                      0.63   1020     0.00   1.40   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.55   816      0.01   12.32  
     Operator          hPsi                        34.18  76488    0.00   75.86  
     Operator          EkineticPW                  0.23   76488    0.00   0.51   
     Operator          VeffPW                      30.36  76488    0.00   67.37  
     PW_Basis_K        recip2real                  17.43  106884   0.00   38.69  
     PW_Basis_K        gathers_scatterp            6.79   106884   0.00   15.06  
     PW_Basis_K        real2recip                  12.77  90360    0.00   28.34  
     PW_Basis_K        gatherp_scatters            4.06   90360    0.00   9.00   
     Operator          NonlocalPW                  3.46   76488    0.00   7.68   
     Nonlocal          add_nonlocal_pp             1.28   76488    0.00   2.85   
     DiagoIterAssist   diagH_LAPACK                0.16   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       42.54  9        4.73   94.41  
     HSolverPW         solve                       42.51  9        4.72   94.34  
     DiagoCG           diag_once                   33.93  918      0.04   75.29  
     DiagoCG_New       spsi_func                   0.31   151344   0.00   0.68   
     DiagoCG_New       hpsi_func                   29.19  75672    0.00   64.78  
     ElecStatePW       psiToRho                    3.00   9        0.33   6.65   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.14   
     ModuleIO          write_rhog                  0.05   1        0.05   0.10   
     Forces            cal_force                   0.12   1        0.12   0.27   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.25   
     FS_Nonlocal_tools cal_becp                    0.28   102      0.00   0.63   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.49   1        0.49   1.10   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.02   1        0.02   0.05   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.45   1        0.45   1.00   
     ModuleIO          write_istate_info           0.26   1        0.26   0.58   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:55:47 2025
     FINISH Time  : Sat Apr 26 13:56:32 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.001
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:56:33 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.120216   SEC) : SETUP UNITCELL
     DONE(0.1949     SEC) : SYMMETRY
     DONE(0.508353   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.518703   SEC) : INIT PLANEWAVE
     DONE(0.52672    SEC) : LOCAL POTENTIAL
     DONE(0.582877   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.6215
     DONE(0.582944   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.48453    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70771740e+03   0.00000000e+00   9.0534e-02  15.49
     CG2     -2.70842637e+03  -7.08962785e-01   2.8679e-02   3.04
     CG3     -2.70851068e+03  -8.43151229e-02   1.9928e-04   3.01
     CG4     -2.70852195e+03  -1.12675110e-02   1.9406e-05   6.03
     CG5     -2.70852208e+03  -1.27463663e-04   2.7790e-06   3.88
     CG6     -2.70852209e+03  -1.33503827e-05   9.0316e-08   3.51
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            32.2403219132         2.6585333730        -3.7692162723 
             2.6585333730        29.1910287876        -2.1834483615 
            -3.7692162723        -2.1834483615        30.7466318513 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 30.725994 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       37.19  17       2.19   100.00 
     Driver            reading                     0.09   1        0.09   0.25   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 37.09  1        37.09  99.75  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.20   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.01   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.15   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    36.57  1        36.57  98.35  
     ESolver_KS_PW     runner                      35.91  1        35.91  96.58  
     ESolver_KS_PW     before_scf                  0.90   1        0.90   2.42   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.05   
     PW_Basis_Sup      recip2real                  0.02   54       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   54       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   7        0.01   0.26   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   7        0.01   0.26   
     H_Hartree_pw      v_hartree                   0.01   7        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   71       0.00   0.08   
     PW_Basis_Sup      gatherp_scatters            0.01   71       0.00   0.03   
     PotXC             cal_v_eff                   0.09   7        0.01   0.24   
     XC_Functional     v_xc                        0.09   7        0.01   0.24   
     Potential         interpolate_vrs             0.00   7        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   8        0.00   0.05   
     Symmetry          group fft grids             0.01   8        0.00   0.02   
     PSIInit           initialize_psi              0.88   1        0.88   2.36   
     Nonlocal          getvnl                      0.53   816      0.00   1.43   
     pp_cell_vnl       getvnl                      0.53   816      0.00   1.42   
     Structure_Factor  get_sk                      0.03   1020     0.00   0.08   
     DiagoIterAssist   diagH_subspace              4.13   612      0.01   11.10  
     Operator          hPsi                        28.18  64350    0.00   75.78  
     Operator          EkineticPW                  0.20   64350    0.00   0.53   
     Operator          VeffPW                      25.00  64350    0.00   67.23  
     PW_Basis_K        recip2real                  14.25  87606    0.00   38.31  
     PW_Basis_K        gathers_scatterp            5.52   87606    0.00   14.85  
     PW_Basis_K        real2recip                  10.45  74754    0.00   28.09  
     PW_Basis_K        gatherp_scatters            3.36   74754    0.00   9.03   
     Operator          NonlocalPW                  2.86   64350    0.00   7.68   
     Nonlocal          add_nonlocal_pp             1.07   64350    0.00   2.88   
     DiagoIterAssist   diagH_LAPACK                0.12   612      0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       34.85  7        4.98   93.73  
     HSolverPW         solve                       34.83  7        4.98   93.66  
     DiagoCG           diag_once                   28.54  714      0.04   76.75  
     DiagoCG_New       spsi_func                   0.26   127476   0.00   0.71   
     DiagoCG_New       hpsi_func                   24.49  63738    0.00   65.85  
     ElecStatePW       psiToRho                    2.28   7        0.33   6.12   
     Charge_Mixing     get_drho                    0.01   7        0.00   0.02   
     Charge_Mixing     inner_product_recip_rho     0.00   7        0.00   0.00   
     Charge            mix_rho                     0.00   5        0.00   0.01   
     Charge            Broyden_mixing              0.00   5        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   20       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.16   
     ModuleIO          write_rhog                  0.04   1        0.04   0.11   
     Forces            cal_force                   0.13   1        0.13   0.34   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.32   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.79   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.51   1        0.51   1.37   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.02   
     Stress_Func       stress_loc                  0.03   1        0.03   0.07   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.46   1        0.46   1.25   
     ModuleIO          write_istate_info           0.02   1        0.02   0.05   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:56:33 2025
     FINISH Time  : Sat Apr 26 13:57:10 2025
     TOTAL  Time  : 37
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.002
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:57:11 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0731544  SEC) : SETUP UNITCELL
     DONE(0.143973   SEC) : SYMMETRY
     DONE(0.446591   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.456788   SEC) : INIT PLANEWAVE
     DONE(0.464682   SEC) : LOCAL POTENTIAL
     DONE(0.519004   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.7896
     DONE(0.519092   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.41032    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70771958e+03   0.00000000e+00   9.0595e-02  15.82
     CG2     -2.70842513e+03  -7.05550121e-01   3.0667e-02   3.09
     CG3     -2.70851044e+03  -8.53114367e-02   3.1739e-04   3.00
     CG4     -2.70852288e+03  -1.24457171e-02   3.0071e-05   5.71
     CG5     -2.70852314e+03  -2.58405780e-04   9.7824e-06   3.87
     CG6     -2.70852316e+03  -2.16609640e-05   2.9331e-07   3.12
     CG7     -2.70852317e+03  -3.96593820e-06   1.3381e-07   5.07
     CG8     -2.70852317e+03  -2.62680698e-07   5.1505e-08   3.06
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -0.3067522456        -2.1640409811         3.0528314189 
            -2.1640409811         2.2086843727         1.7567070783 
             3.0528314189         1.7567070783         0.9757480829 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 0.959227 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       44.93  17       2.64   100.00 
     Driver            reading                     0.05   1        0.05   0.10   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 44.89  1        44.89  99.90  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.01   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.16   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.39  1        44.39  98.78  
     ESolver_KS_PW     runner                      43.71  1        43.71  97.29  
     ESolver_KS_PW     before_scf                  0.89   1        0.89   1.98   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.01   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   68       0.00   0.05   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.22   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.07   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.09   9        0.01   0.20   
     XC_Functional     v_xc                        0.09   9        0.01   0.20   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.86   1        0.86   1.92   
     Nonlocal          getvnl                      0.67   1020     0.00   1.49   
     pp_cell_vnl       getvnl                      0.67   1020     0.00   1.48   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.51   816      0.01   12.27  
     Operator          hPsi                        34.13  75778    0.00   75.96  
     Operator          EkineticPW                  0.24   75778    0.00   0.53   
     Operator          VeffPW                      30.31  75778    0.00   67.46  
     PW_Basis_K        recip2real                  17.40  106174   0.00   38.73  
     PW_Basis_K        gathers_scatterp            6.85   106174   0.00   15.24  
     PW_Basis_K        real2recip                  12.79  89650    0.00   28.46  
     PW_Basis_K        gatherp_scatters            4.14   89650    0.00   9.21   
     Operator          NonlocalPW                  3.44   75778    0.00   7.66   
     Nonlocal          add_nonlocal_pp             1.32   75778    0.00   2.94   
     DiagoIterAssist   diagH_LAPACK                0.16   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       42.65  9        4.74   94.92  
     HSolverPW         solve                       42.62  9        4.74   94.85  
     DiagoCG           diag_once                   34.01  918      0.04   75.69  
     DiagoCG_New       spsi_func                   0.31   149924   0.00   0.69   
     DiagoCG_New       hpsi_func                   29.18  74962    0.00   64.94  
     ElecStatePW       psiToRho                    3.02   9        0.34   6.71   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.02   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.15   
     ModuleIO          write_rhog                  0.05   1        0.05   0.11   
     Forces            cal_force                   0.13   1        0.13   0.28   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.26   
     FS_Nonlocal_tools cal_becp                    0.30   102      0.00   0.66   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.52   1        0.52   1.16   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.47   1        0.47   1.06   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:57:11 2025
     FINISH Time  : Sat Apr 26 13:57:56 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.003
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:57:57 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0936169  SEC) : SETUP UNITCELL
     DONE(0.167333   SEC) : SYMMETRY
     DONE(0.482583   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.493319   SEC) : INIT PLANEWAVE
     DONE(0.501162   SEC) : LOCAL POTENTIAL
     DONE(0.557845   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.6776
     DONE(0.557915   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.44197    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70775109e+03   0.00000000e+00   8.9536e-02  15.83
     CG2     -2.70842588e+03  -6.74786829e-01   2.9301e-02   3.10
     CG3     -2.70851083e+03  -8.49522805e-02   2.4896e-04   3.00
     CG4     -2.70852254e+03  -1.17084084e-02   1.8658e-05   5.88
     CG5     -2.70852273e+03  -1.84366438e-04   5.6970e-06   4.03
     CG6     -2.70852273e+03  -7.95614166e-06   6.3497e-07   3.05
     CG7     -2.70852274e+03  -3.68601159e-06   6.6679e-08   4.13
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           -17.8175836810        -4.1352423822         5.8193780167 
            -4.1352423822       -12.9792663852         3.3376421959 
             5.8193780167         3.3376421959       -15.3044861544 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -15.367112 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       41.24  17       2.43   100.00 
     Driver            reading                     0.07   1        0.07   0.17   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 41.17  1        41.17  99.83  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.18   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.01   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.13   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    40.65  1        40.65  98.58  
     ESolver_KS_PW     runner                      39.96  1        39.96  96.89  
     ESolver_KS_PW     before_scf                  0.88   1        0.88   2.14   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.05   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.09   8        0.01   0.21   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.08   8        0.01   0.20   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.18   
     XC_Functional     v_xc                        0.08   8        0.01   0.18   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.05   
     Symmetry          group fft grids             0.01   9        0.00   0.02   
     PSIInit           initialize_psi              0.86   1        0.86   2.09   
     Nonlocal          getvnl                      0.58   918      0.00   1.41   
     pp_cell_vnl       getvnl                      0.58   918      0.00   1.41   
     Structure_Factor  get_sk                      0.03   1122     0.00   0.07   
     DiagoIterAssist   diagH_subspace              4.91   714      0.01   11.90  
     Operator          hPsi                        31.25  69939    0.00   75.78  
     Operator          EkineticPW                  0.22   69939    0.00   0.54   
     Operator          VeffPW                      27.74  69939    0.00   67.27  
     PW_Basis_K        recip2real                  16.05  96765    0.00   38.92  
     PW_Basis_K        gathers_scatterp            6.23   96765    0.00   15.10  
     PW_Basis_K        real2recip                  11.52  82077    0.00   27.94  
     PW_Basis_K        gatherp_scatters            3.68   82077    0.00   8.93   
     Operator          NonlocalPW                  3.15   69939    0.00   7.63   
     Nonlocal          add_nonlocal_pp             1.19   69939    0.00   2.88   
     DiagoIterAssist   diagH_LAPACK                0.14   714      0.00   0.33   
     ESolver_KS_PW     hamilt2density_single       38.92  8        4.87   94.38  
     HSolverPW         solve                       38.89  8        4.86   94.31  
     DiagoCG           diag_once                   31.32  816      0.04   75.94  
     DiagoCG_New       spsi_func                   0.30   138450   0.00   0.72   
     DiagoCG_New       hpsi_func                   26.84  69225    0.00   65.08  
     ElecStatePW       psiToRho                    2.69   8        0.34   6.52   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.02   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.15   
     ModuleIO          write_rhog                  0.05   1        0.05   0.11   
     Forces            cal_force                   0.12   1        0.12   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.28   
     FS_Nonlocal_tools cal_becp                    0.32   102      0.00   0.79   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.55   1        0.55   1.32   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.50   1        0.50   1.21   
     ModuleIO          write_istate_info           0.02   1        0.02   0.05   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:57:57 2025
     FINISH Time  : Sat Apr 26 13:58:39 2025
     TOTAL  Time  : 42
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.004
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:58:40 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0666628  SEC) : SETUP UNITCELL
     DONE(0.140573   SEC) : SYMMETRY
     DONE(0.454379   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.465309   SEC) : INIT PLANEWAVE
     DONE(0.473173   SEC) : LOCAL POTENTIAL
     DONE(0.531051   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.531119   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.40777    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70775090e+03   0.00000000e+00   9.1211e-02  15.53
     CG2     -2.70842292e+03  -6.72021855e-01   2.9010e-02   3.07
     CG3     -2.70850848e+03  -8.55662679e-02   2.1230e-04   3.06
     CG4     -2.70852033e+03  -1.18455053e-02   1.8966e-05   5.91
     CG5     -2.70852048e+03  -1.52624718e-04   4.3919e-06   3.90
     CG6     -2.70852050e+03  -1.34023904e-05   1.8731e-07   3.24
     CG7     -2.70852050e+03  -1.69347695e-06   2.0178e-08   4.72
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            46.9125401132        -5.2130823517         7.4099306408 
            -5.2130823517        40.9736566763        -4.3069347013 
             7.4099306408        -4.3069347013        44.0655368496 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 43.983911 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       41.57  17       2.45   100.00 
     Driver            reading                     0.04   1        0.04   0.10   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 41.53  1        41.53  99.90  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.01   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.18   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.13   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    41.01  1        41.01  98.65  
     ESolver_KS_PW     runner                      40.37  1        40.37  97.10  
     ESolver_KS_PW     before_scf                  0.88   1        0.88   2.11   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.09   8        0.01   0.22   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.19   
     XC_Functional     v_xc                        0.08   8        0.01   0.19   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.05   
     Symmetry          group fft grids             0.01   9        0.00   0.02   
     PSIInit           initialize_psi              0.85   1        0.85   2.05   
     Nonlocal          getvnl                      0.59   918      0.00   1.43   
     pp_cell_vnl       getvnl                      0.59   918      0.00   1.43   
     Structure_Factor  get_sk                      0.03   1122     0.00   0.06   
     DiagoIterAssist   diagH_subspace              4.87   714      0.01   11.72  
     Operator          hPsi                        31.72  71470    0.00   76.29  
     Operator          EkineticPW                  0.25   71470    0.00   0.59   
     Operator          VeffPW                      28.15  71470    0.00   67.71  
     PW_Basis_K        recip2real                  16.09  98296    0.00   38.69  
     PW_Basis_K        gathers_scatterp            6.31   98296    0.00   15.17  
     PW_Basis_K        real2recip                  11.78  83608    0.00   28.33  
     PW_Basis_K        gatherp_scatters            3.81   83608    0.00   9.16   
     Operator          NonlocalPW                  3.19   71470    0.00   7.67   
     Nonlocal          add_nonlocal_pp             1.21   71470    0.00   2.92   
     DiagoIterAssist   diagH_LAPACK                0.14   714      0.00   0.33   
     ESolver_KS_PW     hamilt2density_single       39.33  8        4.92   94.61  
     HSolverPW         solve                       39.30  8        4.91   94.54  
     DiagoCG           diag_once                   31.81  816      0.04   76.52  
     DiagoCG_New       spsi_func                   0.29   141512   0.00   0.69   
     DiagoCG_New       hpsi_func                   27.33  70756    0.00   65.73  
     ElecStatePW       psiToRho                    2.62   8        0.33   6.30   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.16   
     ModuleIO          write_rhog                  0.05   1        0.05   0.12   
     Forces            cal_force                   0.12   1        0.12   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.27   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.69   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.50   1        0.50   1.21   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.02   1        0.02   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.46   1        0.46   1.10   
     ModuleIO          write_istate_info           0.02   1        0.02   0.05   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:58:40 2025
     FINISH Time  : Sat Apr 26 13:59:21 2025
     TOTAL  Time  : 41
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.005
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 13:59:22 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0662017  SEC) : SETUP UNITCELL
     DONE(0.140037   SEC) : SYMMETRY
     DONE(0.45203    SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.462531   SEC) : INIT PLANEWAVE
     DONE(0.470493   SEC) : LOCAL POTENTIAL
     DONE(0.526772   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.526851   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.45036    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70773689e+03   0.00000000e+00   9.2332e-02  15.70
     CG2     -2.70842236e+03  -6.85471924e-01   3.0227e-02   3.10
     CG3     -2.70850989e+03  -8.75272370e-02   2.2133e-04   2.98
     CG4     -2.70852193e+03  -1.20440864e-02   1.6611e-05   5.88
     CG5     -2.70852208e+03  -1.47179671e-04   3.0895e-06   4.02
     CG6     -2.70852209e+03  -1.08451030e-05   1.1366e-07   3.39
     CG7     -2.70852209e+03  -8.70478588e-07   2.6235e-08   4.75
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            29.6369177580        -2.5604547848         3.6301623809 
            -2.5604547848        26.7001191257        -2.1028966051 
             3.6301623809        -2.1028966051        28.1983328959 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 28.178457 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       42.08  17       2.48   100.00 
     Driver            reading                     0.04   1        0.04   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 42.04  1        42.04  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.18   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.13   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    41.52  1        41.52  98.68  
     ESolver_KS_PW     runner                      40.81  1        40.81  97.00  
     ESolver_KS_PW     before_scf                  0.92   1        0.92   2.19   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.05   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.09   8        0.01   0.21   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.21   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.19   
     XC_Functional     v_xc                        0.08   8        0.01   0.19   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.05   
     Symmetry          group fft grids             0.01   9        0.00   0.02   
     PSIInit           initialize_psi              0.90   1        0.90   2.14   
     Nonlocal          getvnl                      0.60   918      0.00   1.43   
     pp_cell_vnl       getvnl                      0.60   918      0.00   1.42   
     Structure_Factor  get_sk                      0.03   1122     0.00   0.07   
     DiagoIterAssist   diagH_subspace              4.88   714      0.01   11.59  
     Operator          hPsi                        32.12  72242    0.00   76.33  
     Operator          EkineticPW                  0.22   72242    0.00   0.53   
     Operator          VeffPW                      28.47  72242    0.00   67.67  
     PW_Basis_K        recip2real                  16.19  99068    0.00   38.47  
     PW_Basis_K        gathers_scatterp            6.29   99068    0.00   14.95  
     PW_Basis_K        real2recip                  11.96  84380    0.00   28.43  
     PW_Basis_K        gatherp_scatters            3.75   84380    0.00   8.91   
     Operator          NonlocalPW                  3.29   72242    0.00   7.82   
     Nonlocal          add_nonlocal_pp             1.23   72242    0.00   2.93   
     DiagoIterAssist   diagH_LAPACK                0.14   714      0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       39.74  8        4.97   94.44  
     HSolverPW         solve                       39.71  8        4.96   94.38  
     DiagoCG           diag_once                   32.26  816      0.04   76.68  
     DiagoCG_New       spsi_func                   0.30   143056   0.00   0.71   
     DiagoCG_New       hpsi_func                   27.72  71528    0.00   65.89  
     ElecStatePW       psiToRho                    2.60   8        0.33   6.19   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.15   
     ModuleIO          write_rhog                  0.04   1        0.04   0.11   
     Forces            cal_force                   0.12   1        0.12   0.29   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.27   
     FS_Nonlocal_tools cal_becp                    0.34   102      0.00   0.82   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.56   1        0.56   1.34   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.52   1        0.52   1.23   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 13:59:22 2025
     FINISH Time  : Sat Apr 26 14:00:04 2025
     TOTAL  Time  : 42
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.006
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:00:05 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0670509  SEC) : SETUP UNITCELL
     DONE(0.145337   SEC) : SYMMETRY
     DONE(0.457887   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.468367   SEC) : INIT PLANEWAVE
     DONE(0.476574   SEC) : LOCAL POTENTIAL
     DONE(0.531757   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.531859   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.45475    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70780287e+03   0.00000000e+00   8.9182e-02  16.21
     CG2     -2.70842940e+03  -6.26531709e-01   2.8462e-02   3.09
     CG3     -2.70851153e+03  -8.21275835e-02   2.5188e-04   3.04
     CG4     -2.70852299e+03  -1.14624567e-02   1.6816e-05   5.71
     CG5     -2.70852315e+03  -1.58918055e-04   6.1721e-06   4.12
     CG6     -2.70852316e+03  -1.33927384e-05   3.1564e-07   2.98
     CG7     -2.70852317e+03  -2.36454623e-06   1.7932e-07   4.73
     CG8     -2.70852317e+03  -4.50180114e-07   3.9479e-08   3.02
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             0.3849927710         2.3381398588        -3.2984342188 
             2.3381398588         3.1027983352         1.8980356084 
            -3.2984342188         1.8980356084         1.7706713064 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 1.752821 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       45.11  17       2.65   100.00 
     Driver            reading                     0.04   1        0.04   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 45.07  1        45.07  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.55  1        44.55  98.76  
     ESolver_KS_PW     runner                      43.89  1        43.89  97.28  
     ESolver_KS_PW     before_scf                  0.92   1        0.92   2.05   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.01   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.05   
     PW_Basis_Sup      recip2real                  0.03   68       0.00   0.07   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.03   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.22   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.09   9        0.01   0.20   
     XC_Functional     v_xc                        0.09   9        0.01   0.20   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.89   1        0.89   1.98   
     Nonlocal          getvnl                      0.65   1020     0.00   1.45   
     pp_cell_vnl       getvnl                      0.65   1020     0.00   1.44   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.61   816      0.01   12.43  
     Operator          hPsi                        34.36  76489    0.00   76.16  
     Operator          EkineticPW                  0.23   76489    0.00   0.52   
     Operator          VeffPW                      30.52  76489    0.00   67.64  
     PW_Basis_K        recip2real                  17.63  106885   0.00   39.08  
     PW_Basis_K        gathers_scatterp            6.84   106885   0.00   15.17  
     PW_Basis_K        real2recip                  12.70  90361    0.00   28.16  
     PW_Basis_K        gatherp_scatters            4.13   90361    0.00   9.15   
     Operator          NonlocalPW                  3.46   76489    0.00   7.67   
     Nonlocal          add_nonlocal_pp             1.30   76489    0.00   2.89   
     DiagoIterAssist   diagH_LAPACK                0.16   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       42.79  9        4.75   94.86  
     HSolverPW         solve                       42.76  9        4.75   94.78  
     DiagoCG           diag_once                   34.15  918      0.04   75.71  
     DiagoCG_New       spsi_func                   0.31   151346   0.00   0.69   
     DiagoCG_New       hpsi_func                   29.31  75673    0.00   64.98  
     ElecStatePW       psiToRho                    2.96   9        0.33   6.56   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.15   
     ModuleIO          write_rhog                  0.05   1        0.05   0.11   
     Forces            cal_force                   0.12   1        0.12   0.27   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.25   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.64   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.52   1        0.52   1.15   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.47   1        0.47   1.05   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:00:05 2025
     FINISH Time  : Sat Apr 26 14:00:50 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.007
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:00:51 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0812697  SEC) : SETUP UNITCELL
     DONE(0.156582   SEC) : SYMMETRY
     DONE(0.467452   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.477997   SEC) : INIT PLANEWAVE
     DONE(0.485765   SEC) : LOCAL POTENTIAL
     DONE(0.542835   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.8737
     DONE(0.542911   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.44696    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70770202e+03   0.00000000e+00   9.2424e-02  15.57
     CG2     -2.70841568e+03  -7.13658852e-01   3.2386e-02   3.05
     CG3     -2.70850966e+03  -9.39783014e-02   2.4608e-04   2.94
     CG4     -2.70852254e+03  -1.28859130e-02   1.9255e-05   5.98
     CG5     -2.70852273e+03  -1.81484601e-04   5.9299e-06   3.96
     CG6     -2.70852273e+03  -8.46435497e-06   5.9391e-07   3.03
     CG7     -2.70852274e+03  -3.90280022e-06   1.9268e-07   4.25
     CG8     -2.70852274e+03  -5.52401197e-07   2.9348e-08   3.11
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           -14.6578543880         4.0873897801        -5.7520367691 
             4.0873897801        -9.8755256053         3.2990193416 
            -5.7520367691         3.2990193416       -12.1738381691 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -12.235739 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       44.08  17       2.59   100.00 
     Driver            reading                     0.04   1        0.04   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 44.04  1        44.04  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    43.50  1        43.50  98.70  
     ESolver_KS_PW     runner                      42.85  1        42.85  97.22  
     ESolver_KS_PW     before_scf                  0.90   1        0.90   2.05   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.03   68       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.02   1        0.02   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.23   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.09   9        0.01   0.20   
     XC_Functional     v_xc                        0.09   9        0.01   0.20   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.88   1        0.88   1.99   
     Nonlocal          getvnl                      0.66   1020     0.00   1.49   
     pp_cell_vnl       getvnl                      0.66   1020     0.00   1.49   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.45   816      0.01   12.37  
     Operator          hPsi                        33.53  74945    0.00   76.07  
     Operator          EkineticPW                  0.23   74945    0.00   0.53   
     Operator          VeffPW                      29.72  74945    0.00   67.43  
     PW_Basis_K        recip2real                  17.14  105341   0.00   38.90  
     PW_Basis_K        gathers_scatterp            6.67   105341   0.00   15.13  
     PW_Basis_K        real2recip                  12.48  88817    0.00   28.32  
     PW_Basis_K        gatherp_scatters            4.07   88817    0.00   9.22   
     Operator          NonlocalPW                  3.44   74945    0.00   7.80   
     Nonlocal          add_nonlocal_pp             1.31   74945    0.00   2.97   
     DiagoIterAssist   diagH_LAPACK                0.15   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       41.78  9        4.64   94.79  
     HSolverPW         solve                       41.75  9        4.64   94.71  
     DiagoCG           diag_once                   33.23  918      0.04   75.38  
     DiagoCG_New       spsi_func                   0.31   148258   0.00   0.70   
     DiagoCG_New       hpsi_func                   28.62  74129    0.00   64.92  
     ElecStatePW       psiToRho                    3.01   9        0.33   6.83   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.15   
     ModuleIO          write_rhog                  0.05   1        0.05   0.11   
     Forces            cal_force                   0.12   1        0.12   0.28   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.26   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.65   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.50   1        0.50   1.14   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.46   1        0.46   1.04   
     ModuleIO          write_istate_info           0.02   1        0.02   0.05   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:00:51 2025
     FINISH Time  : Sat Apr 26 14:01:36 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.008
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:01:37 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0674633  SEC) : SETUP UNITCELL
     DONE(0.139098   SEC) : SYMMETRY
     DONE(0.442178   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.455963   SEC) : INIT PLANEWAVE
     DONE(0.463788   SEC) : LOCAL POTENTIAL
     DONE(0.520357   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.520429   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.49332    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777838e+03   0.00000000e+00   8.8453e-02  15.84
     CG2     -2.70843058e+03  -6.52203801e-01   2.6534e-02   3.00
     CG3     -2.70850977e+03  -7.91930956e-02   2.3176e-04   3.04
     CG4     -2.70852035e+03  -1.05766893e-02   1.3566e-05   5.85
     CG5     -2.70852049e+03  -1.38111484e-04   3.0615e-06   4.19
     CG6     -2.70852050e+03  -7.15944180e-06   2.1510e-07   3.17
     CG7     -2.70852050e+03  -1.19029998e-06   5.8053e-08   4.44
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            36.6082840362         0.0000000000        -0.0000000000 
             0.0000000000        48.6792008137         8.6220834125 
            -0.0000000000         8.6220834125        42.7669150451 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 42.684800 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       41.92  17       2.47   100.00 
     Driver            reading                     0.04   1        0.04   0.10   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 41.88  1        41.88  99.90  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.13   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    41.37  1        41.37  98.69  
     ESolver_KS_PW     runner                      40.66  1        40.66  97.00  
     ESolver_KS_PW     before_scf                  0.97   1        0.97   2.32   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.09   8        0.01   0.21   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.21   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.19   
     XC_Functional     v_xc                        0.08   8        0.01   0.19   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.05   
     Symmetry          group fft grids             0.01   9        0.00   0.02   
     PSIInit           initialize_psi              0.95   1        0.95   2.26   
     Nonlocal          getvnl                      0.59   918      0.00   1.41   
     pp_cell_vnl       getvnl                      0.59   918      0.00   1.40   
     Structure_Factor  get_sk                      0.03   1122     0.00   0.07   
     DiagoIterAssist   diagH_subspace              4.89   714      0.01   11.68  
     Operator          hPsi                        31.88  72029    0.00   76.04  
     Operator          EkineticPW                  0.22   72029    0.00   0.52   
     Operator          VeffPW                      28.33  72029    0.00   67.59  
     PW_Basis_K        recip2real                  16.13  98855    0.00   38.49  
     PW_Basis_K        gathers_scatterp            6.22   98855    0.00   14.84  
     PW_Basis_K        real2recip                  11.87  84167    0.00   28.32  
     PW_Basis_K        gatherp_scatters            3.77   84167    0.00   8.99   
     Operator          NonlocalPW                  3.19   72029    0.00   7.61   
     Nonlocal          add_nonlocal_pp             1.21   72029    0.00   2.88   
     DiagoIterAssist   diagH_LAPACK                0.14   714      0.00   0.33   
     ESolver_KS_PW     hamilt2density_single       39.44  8        4.93   94.08  
     HSolverPW         solve                       39.41  8        4.93   94.01  
     DiagoCG           diag_once                   32.02  816      0.04   76.37  
     DiagoCG_New       spsi_func                   0.29   142630   0.00   0.68   
     DiagoCG_New       hpsi_func                   27.48  71315    0.00   65.55  
     ElecStatePW       psiToRho                    2.60   8        0.33   6.21   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.16   1        0.16   0.39   
     ModuleIO          write_rhog                  0.15   1        0.15   0.35   
     Forces            cal_force                   0.12   1        0.12   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.27   
     FS_Nonlocal_tools cal_becp                    0.30   102      0.00   0.71   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.56   1        0.56   1.33   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.02   1        0.02   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.51   1        0.51   1.22   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:01:37 2025
     FINISH Time  : Sat Apr 26 14:02:19 2025
     TOTAL  Time  : 42
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.009
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:02:20 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0690262  SEC) : SETUP UNITCELL
     DONE(0.149485   SEC) : SYMMETRY
     DONE(0.467355   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.47775    SEC) : INIT PLANEWAVE
     DONE(0.485499   SEC) : LOCAL POTENTIAL
     DONE(0.539381   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5935
     DONE(0.539483   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.44751    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70774389e+03   0.00000000e+00   9.0004e-02  15.78
     CG2     -2.70842771e+03  -6.83819460e-01   2.8925e-02   3.03
     CG3     -2.70851038e+03  -8.26668421e-02   2.9584e-04   2.95
     CG4     -2.70852186e+03  -1.14874188e-02   2.1972e-05   5.73
     CG5     -2.70852206e+03  -1.96529662e-04   7.3389e-06   3.99
     CG6     -2.70852209e+03  -2.79496585e-05   2.8765e-07   3.06
     CG7     -2.70852209e+03  -1.83439598e-06   1.0449e-07   4.68
     CG8     -2.70852209e+03  -3.21343049e-07   1.7010e-08   3.25
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            27.4725852353         0.0000000000        -0.0000000000 
             0.0000000000        33.4400789815         4.2409131385 
            -0.0000000000         4.2409131385        30.4864709657 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 30.466378 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       44.69  17       2.63   100.00 
     Driver            reading                     0.04   1        0.04   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 44.65  1        44.65  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.01   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.13  1        44.13  98.73  
     ESolver_KS_PW     runner                      43.45  1        43.45  97.22  
     ESolver_KS_PW     before_scf                  0.91   1        0.91   2.03   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.01   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   68       0.00   0.05   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.23   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.06   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.09   9        0.01   0.20   
     XC_Functional     v_xc                        0.09   9        0.01   0.20   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.88   1        0.88   1.97   
     Nonlocal          getvnl                      0.64   1020     0.00   1.44   
     pp_cell_vnl       getvnl                      0.64   1020     0.00   1.43   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.55   816      0.01   12.42  
     Operator          hPsi                        34.11  76351    0.00   76.32  
     Operator          EkineticPW                  0.23   76351    0.00   0.51   
     Operator          VeffPW                      30.42  76351    0.00   68.07  
     PW_Basis_K        recip2real                  17.45  106747   0.00   39.05  
     PW_Basis_K        gathers_scatterp            6.77   106747   0.00   15.15  
     PW_Basis_K        real2recip                  12.69  90223    0.00   28.39  
     PW_Basis_K        gatherp_scatters            3.96   90223    0.00   8.86   
     Operator          NonlocalPW                  3.32   76351    0.00   7.43   
     Nonlocal          add_nonlocal_pp             1.29   76351    0.00   2.89   
     DiagoIterAssist   diagH_LAPACK                0.15   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       42.38  9        4.71   94.82  
     HSolverPW         solve                       42.34  9        4.70   94.74  
     DiagoCG           diag_once                   33.85  918      0.04   75.73  
     DiagoCG_New       spsi_func                   0.31   151070   0.00   0.70   
     DiagoCG_New       hpsi_func                   29.10  75535    0.00   65.12  
     ElecStatePW       psiToRho                    2.92   9        0.32   6.52   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.14   
     ModuleIO          write_rhog                  0.04   1        0.04   0.10   
     Forces            cal_force                   0.12   1        0.12   0.27   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.11   1        0.11   0.25   
     FS_Nonlocal_tools cal_becp                    0.30   102      0.00   0.68   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.53   1        0.53   1.18   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.48   1        0.48   1.08   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:02:20 2025
     FINISH Time  : Sat Apr 26 14:03:04 2025
     TOTAL  Time  : 44
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.010
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:03:05 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0647823  SEC) : SETUP UNITCELL
     DONE(0.141607   SEC) : SYMMETRY
     DONE(0.454725   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.465724   SEC) : INIT PLANEWAVE
     DONE(0.473824   SEC) : LOCAL POTENTIAL
     DONE(0.530337   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.5375
     DONE(0.530407   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.43076    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777414e+03   0.00000000e+00   9.0801e-02  15.85
     CG2     -2.70842507e+03  -6.50927711e-01   3.0297e-02   3.10
     CG3     -2.70851061e+03  -8.55463164e-02   2.8659e-04   2.97
     CG4     -2.70852297e+03  -1.23581755e-02   2.1136e-05   5.79
     CG5     -2.70852315e+03  -1.79560604e-04   3.6175e-06   4.02
     CG6     -2.70852316e+03  -1.21469184e-05   4.5818e-07   3.48
     CG7     -2.70852317e+03  -1.84694005e-06   1.1767e-07   3.89
     CG8     -2.70852317e+03  -2.83746801e-07   4.1346e-08   3.44
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             5.1425258529         0.0000000000        -0.0000000000 
             0.0000000000        -0.3262509002        -3.8478178962 
            -0.0000000000        -3.8478178962         2.4352106286 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 2.417162 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       44.73  17       2.63   100.00 
     Driver            reading                     0.04   1        0.04   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 44.69  1        44.69  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.17   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.17  1        44.17  98.75  
     ESolver_KS_PW     runner                      43.52  1        43.52  97.29  
     ESolver_KS_PW     before_scf                  0.90   1        0.90   2.01   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.03   68       0.00   0.06   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.11   9        0.01   0.24   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.11   9        0.01   0.24   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.08   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.03   
     PotXC             cal_v_eff                   0.10   9        0.01   0.21   
     XC_Functional     v_xc                        0.10   9        0.01   0.21   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.88   1        0.88   1.96   
     Nonlocal          getvnl                      0.64   1020     0.00   1.42   
     pp_cell_vnl       getvnl                      0.63   1020     0.00   1.42   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.53   816      0.01   12.36  
     Operator          hPsi                        34.15  75758    0.00   76.36  
     Operator          EkineticPW                  0.24   75758    0.00   0.53   
     Operator          VeffPW                      30.27  75758    0.00   67.68  
     PW_Basis_K        recip2real                  17.48  106154   0.00   39.08  
     PW_Basis_K        gathers_scatterp            6.82   106154   0.00   15.25  
     PW_Basis_K        real2recip                  12.61  89630    0.00   28.19  
     PW_Basis_K        gatherp_scatters            4.06   89630    0.00   9.08   
     Operator          NonlocalPW                  3.50   75758    0.00   7.83   
     Nonlocal          add_nonlocal_pp             1.32   75758    0.00   2.96   
     DiagoIterAssist   diagH_LAPACK                0.16   816      0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       42.44  9        4.72   94.89  
     HSolverPW         solve                       42.41  9        4.71   94.82  
     DiagoCG           diag_once                   33.88  918      0.04   75.74  
     DiagoCG_New       spsi_func                   0.31   149884   0.00   0.69   
     DiagoCG_New       hpsi_func                   29.17  74942    0.00   65.21  
     ElecStatePW       psiToRho                    2.97   9        0.33   6.64   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.06   1        0.06   0.14   
     ModuleIO          write_rhog                  0.05   1        0.05   0.10   
     Forces            cal_force                   0.13   1        0.13   0.28   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.26   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.64   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.50   1        0.50   1.12   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.46   1        0.46   1.02   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:03:05 2025
     FINISH Time  : Sat Apr 26 14:03:50 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.011
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:03:51 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.119546   SEC) : SETUP UNITCELL
     DONE(0.21359    SEC) : SYMMETRY
     DONE(0.523484   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       102             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.533807   SEC) : INIT PLANEWAVE
     DONE(0.541561   SEC) : LOCAL POTENTIAL
     DONE(0.598009   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 17.7896
     DONE(0.598076   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(1.49764    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70772397e+03   0.00000000e+00   9.0218e-02  15.84
     CG2     -2.70842388e+03  -6.99910385e-01   3.0304e-02   3.08
     CG3     -2.70850974e+03  -8.58626765e-02   2.6038e-04   3.01
     CG4     -2.70852253e+03  -1.27856537e-02   2.0718e-05   6.12
     CG5     -2.70852273e+03  -1.97774104e-04   5.0100e-06   4.00
     CG6     -2.70852273e+03  -8.67843788e-06   5.9282e-07   3.23
     CG7     -2.70852274e+03  -3.65264890e-06   1.4321e-07   4.15
     CG8     -2.70852274e+03  -4.81027722e-07   1.7560e-08   3.28
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            -7.7857229591        -0.0000000000         0.0000000000 
            -0.0000000000       -17.1164087210        -6.5327867214 
             0.0000000000        -6.5327867214       -12.3595885287 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: -12.420573 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       44.98  17       2.65   100.00 
     Driver            reading                     0.05   1        0.05   0.12   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 44.93  1        44.93  99.88  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.09   1        0.09   0.21   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.01   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.12   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    44.36  1        44.36  98.61  
     ESolver_KS_PW     runner                      43.69  1        43.69  97.14  
     ESolver_KS_PW     before_scf                  0.90   1        0.90   2.00   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.02   2        0.01   0.04   
     PW_Basis_Sup      recip2real                  0.02   68       0.00   0.05   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.03   
     Potential         update_from_charge          0.10   9        0.01   0.23   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.22   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.07   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.02   
     PotXC             cal_v_eff                   0.09   9        0.01   0.20   
     XC_Functional     v_xc                        0.09   9        0.01   0.20   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   10       0.00   0.05   
     Symmetry          group fft grids             0.01   10       0.00   0.02   
     PSIInit           initialize_psi              0.88   1        0.88   1.95   
     Nonlocal          getvnl                      0.67   1020     0.00   1.48   
     pp_cell_vnl       getvnl                      0.67   1020     0.00   1.48   
     Structure_Factor  get_sk                      0.03   1224     0.00   0.07   
     DiagoIterAssist   diagH_subspace              5.60   816      0.01   12.44  
     Operator          hPsi                        34.08  75220    0.00   75.76  
     Operator          EkineticPW                  0.24   75220    0.00   0.53   
     Operator          VeffPW                      30.29  75220    0.00   67.34  
     PW_Basis_K        recip2real                  17.39  105616   0.00   38.66  
     PW_Basis_K        gathers_scatterp            6.79   105616   0.00   15.10  
     PW_Basis_K        real2recip                  12.78  89092    0.00   28.42  
     PW_Basis_K        gatherp_scatters            4.09   89092    0.00   9.10   
     Operator          NonlocalPW                  3.41   75220    0.00   7.58   
     Nonlocal          add_nonlocal_pp             1.30   75220    0.00   2.89   
     DiagoIterAssist   diagH_LAPACK                0.15   816      0.00   0.34   
     ESolver_KS_PW     hamilt2density_single       42.62  9        4.74   94.74  
     HSolverPW         solve                       42.58  9        4.73   94.67  
     DiagoCG           diag_once                   33.84  918      0.04   75.23  
     DiagoCG_New       spsi_func                   0.31   148808   0.00   0.70   
     DiagoCG_New       hpsi_func                   29.04  74404    0.00   64.56  
     ElecStatePW       psiToRho                    3.09   9        0.34   6.86   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.02   
     Charge            Broyden_mixing              0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.16   
     ModuleIO          write_rhog                  0.05   1        0.05   0.12   
     Forces            cal_force                   0.13   1        0.13   0.28   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.12   1        0.12   0.26   
     FS_Nonlocal_tools cal_becp                    0.29   102      0.00   0.64   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.01   1        0.01   0.02   
     Stress_PW         cal_stress                  0.51   1        0.51   1.13   
     Stress_Func       stress_kin                  0.01   1        0.01   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.03   1        0.03   0.06   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.46   1        0.46   1.03   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:03:51 2025
     FINISH Time  : Sat Apr 26 14:04:36 2025
     TOTAL  Time  : 45
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.012
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:04:37 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0700802  SEC) : SETUP UNITCELL
     DONE(0.14431    SEC) : SYMMETRY
     DONE(0.457222   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.47449    SEC) : INIT PLANEWAVE
     DONE(0.489305   SEC) : LOCAL POTENTIAL
     DONE(0.547166   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.495
     DONE(0.547238   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.03514    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70772751e+03   0.00000000e+00   8.9183e-02  26.46
     CG2     -2.70842391e+03  -6.96397559e-01   2.8957e-02   5.17
     CG3     -2.70850917e+03  -8.52615149e-02   2.4974e-04   5.04
     CG4     -2.70852136e+03  -1.21909023e-02   1.4245e-05   9.87
     CG5     -2.70852152e+03  -1.61513624e-04   3.4971e-06   7.02
     CG6     -2.70852153e+03  -9.57196478e-06   9.9499e-08   5.37
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            46.3576101459         0.0000000000        -0.0000000000 
             0.0000000000        10.2729060537        12.8628707138 
            -0.0000000000        12.8628707138         1.4472114286 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 19.359243 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       62.22  17       3.66   100.00 
     Driver            reading                     0.04   1        0.04   0.07   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 62.18  1        62.18  99.93  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.12   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.09   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    61.62  1        61.62  99.04  
     ESolver_KS_PW     runner                      60.50  1        60.50  97.24  
     ESolver_KS_PW     before_scf                  1.49   1        1.49   2.39   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   54       0.00   0.04   
     PW_Basis_Sup      gathers_scatterp            0.01   54       0.00   0.01   
     Potential         init_pot                    0.02   1        0.02   0.02   
     Potential         update_from_charge          0.08   7        0.01   0.13   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.08   7        0.01   0.13   
     H_Hartree_pw      v_hartree                   0.01   7        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   71       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   71       0.00   0.01   
     PotXC             cal_v_eff                   0.07   7        0.01   0.12   
     XC_Functional     v_xc                        0.07   7        0.01   0.12   
     Potential         interpolate_vrs             0.00   7        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   8        0.00   0.03   
     Symmetry          group fft grids             0.01   8        0.00   0.01   
     PSIInit           initialize_psi              1.45   1        1.45   2.33   
     Nonlocal          getvnl                      0.88   1368     0.00   1.41   
     pp_cell_vnl       getvnl                      0.87   1368     0.00   1.40   
     Structure_Factor  get_sk                      0.04   1710     0.00   0.07   
     DiagoIterAssist   diagH_subspace              6.99   1026     0.01   11.23  
     Operator          hPsi                        47.53  107432   0.00   76.38  
     Operator          EkineticPW                  0.33   107432   0.00   0.52   
     Operator          VeffPW                      42.26  107432   0.00   67.92  
     PW_Basis_K        recip2real                  24.23  146420   0.00   38.94  
     PW_Basis_K        gathers_scatterp            9.52   146420   0.00   15.30  
     PW_Basis_K        real2recip                  17.59  124874   0.00   28.27  
     PW_Basis_K        gatherp_scatters            5.63   124874   0.00   9.05   
     Operator          NonlocalPW                  4.73   107432   0.00   7.60   
     Nonlocal          add_nonlocal_pp             1.72   107432   0.00   2.77   
     DiagoIterAssist   diagH_LAPACK                0.20   1026     0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       58.86  7        8.41   94.59  
     HSolverPW         solve                       58.83  7        8.40   94.55  
     DiagoCG           diag_once                   48.08  1197     0.04   77.28  
     DiagoCG_New       spsi_func                   0.43   212812   0.00   0.70   
     DiagoCG_New       hpsi_func                   41.26  106406   0.00   66.31  
     ElecStatePW       psiToRho                    3.93   7        0.56   6.32   
     Charge_Mixing     get_drho                    0.01   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   7        0.00   0.00   
     Charge            mix_rho                     0.01   5        0.00   0.01   
     Charge            Broyden_mixing              0.00   5        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   20       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.13   
     ModuleIO          write_rhog                  0.05   1        0.05   0.08   
     Forces            cal_force                   0.22   1        0.22   0.35   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.20   1        0.20   0.32   
     FS_Nonlocal_tools cal_becp                    0.59   171      0.00   0.95   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.88   1        0.88   1.42   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.09   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.80   1        0.80   1.29   
     ModuleIO          write_istate_info           0.04   1        0.04   0.06   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:04:37 2025
     FINISH Time  : Sat Apr 26 14:05:39 2025
     TOTAL  Time  : 62
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.013
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:05:40 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0653176  SEC) : SETUP UNITCELL
     DONE(0.140493   SEC) : SYMMETRY
     DONE(0.457016   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.473932   SEC) : INIT PLANEWAVE
     DONE(0.489277   SEC) : LOCAL POTENTIAL
     DONE(0.545698   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.1193
     DONE(0.545814   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.00842    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777017e+03   0.00000000e+00   8.9086e-02  26.67
     CG2     -2.70842593e+03  -6.55755383e-01   2.8032e-02   5.07
     CG3     -2.70851149e+03  -8.55564765e-02   2.1747e-04   5.17
     CG4     -2.70852262e+03  -1.11385015e-02   1.4224e-05   9.93
     CG5     -2.70852276e+03  -1.36479666e-04   2.2336e-06   7.00
     CG6     -2.70852277e+03  -8.38202443e-06   9.5519e-08   5.89
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            24.7558938148         0.0000000000        -0.0000000000 
            -0.0000000000        12.8560048475         5.4959819862 
            -0.0000000000         5.4959819862         9.0277150104 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 15.546538 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       62.97  17       3.70   100.00 
     Driver            reading                     0.04   1        0.04   0.06   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 62.93  1        62.93  99.94  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.01   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.12   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.09   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    62.39  1        62.39  99.08  
     ESolver_KS_PW     runner                      61.28  1        61.28  97.32  
     ESolver_KS_PW     before_scf                  1.46   1        1.46   2.32   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   54       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   54       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.10   7        0.01   0.16   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   7        0.01   0.16   
     H_Hartree_pw      v_hartree                   0.01   7        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   71       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   71       0.00   0.01   
     PotXC             cal_v_eff                   0.09   7        0.01   0.14   
     XC_Functional     v_xc                        0.09   7        0.01   0.14   
     Potential         interpolate_vrs             0.00   7        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   8        0.00   0.03   
     Symmetry          group fft grids             0.01   8        0.00   0.01   
     PSIInit           initialize_psi              1.43   1        1.43   2.27   
     Nonlocal          getvnl                      0.88   1368     0.00   1.39   
     pp_cell_vnl       getvnl                      0.88   1368     0.00   1.39   
     Structure_Factor  get_sk                      0.04   1710     0.00   0.07   
     DiagoIterAssist   diagH_subspace              6.98   1026     0.01   11.09  
     Operator          hPsi                        48.22  109204   0.00   76.58  
     Operator          EkineticPW                  0.34   109204   0.00   0.53   
     Operator          VeffPW                      42.91  109204   0.00   68.14  
     PW_Basis_K        recip2real                  24.32  148192   0.00   38.63  
     PW_Basis_K        gathers_scatterp            9.47   148192   0.00   15.04  
     PW_Basis_K        real2recip                  18.03  126646   0.00   28.63  
     PW_Basis_K        gatherp_scatters            5.85   126646   0.00   9.29   
     Operator          NonlocalPW                  4.78   109204   0.00   7.59   
     Nonlocal          add_nonlocal_pp             1.74   109204   0.00   2.76   
     DiagoIterAssist   diagH_LAPACK                0.19   1026     0.00   0.31   
     ESolver_KS_PW     hamilt2density_single       59.63  7        8.52   94.70  
     HSolverPW         solve                       59.61  7        8.52   94.66  
     DiagoCG           diag_once                   48.91  1197     0.04   77.67  
     DiagoCG_New       spsi_func                   0.45   216356   0.00   0.72   
     DiagoCG_New       hpsi_func                   41.96  108178   0.00   66.63  
     ElecStatePW       psiToRho                    3.88   7        0.55   6.16   
     Charge_Mixing     get_drho                    0.01   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   7        0.00   0.00   
     Charge            mix_rho                     0.00   5        0.00   0.01   
     Charge            Broyden_mixing              0.00   5        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   20       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.13   
     ModuleIO          write_rhog                  0.05   1        0.05   0.08   
     Forces            cal_force                   0.21   1        0.21   0.34   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.20   1        0.20   0.31   
     FS_Nonlocal_tools cal_becp                    0.59   171      0.00   0.93   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.86   1        0.86   1.36   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.78   1        0.78   1.23   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:05:40 2025
     FINISH Time  : Sat Apr 26 14:06:43 2025
     TOTAL  Time  : 63
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.014
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:06:44 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0768838  SEC) : SETUP UNITCELL
     DONE(0.152881   SEC) : SYMMETRY
     DONE(0.473234   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.490347   SEC) : INIT PLANEWAVE
     DONE(0.505417   SEC) : LOCAL POTENTIAL
     DONE(0.56477    SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.5889
     DONE(0.564855   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.08836    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70774425e+03   0.00000000e+00   9.0798e-02  26.66
     CG2     -2.70842283e+03  -6.78589311e-01   3.0107e-02   5.18
     CG3     -2.70851011e+03  -8.72748559e-02   2.6155e-04   5.00
     CG4     -2.70852257e+03  -1.24636030e-02   1.7634e-05   9.85
     CG5     -2.70852275e+03  -1.81275516e-04   6.9635e-06   6.94
     CG6     -2.70852277e+03  -1.13748411e-05   4.4418e-07   4.94
     CG7     -2.70852277e+03  -3.34913143e-06   1.0241e-07   7.69
     CG8     -2.70852277e+03  -1.98906210e-07   1.0008e-08   5.28
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
             5.1226064901         0.0000000000        -0.0000000000 
             0.0000000000        18.2860785465        -3.1929177276 
            -0.0000000000        -3.1929177276        20.5778774474 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 14.662187 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       74.91  17       4.41   100.00 
     Driver            reading                     0.05   1        0.05   0.07   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 74.86  1        74.86  99.93  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.10   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    74.30  1        74.30  99.18  
     ESolver_KS_PW     runner                      73.21  1        73.21  97.73  
     ESolver_KS_PW     before_scf                  1.52   1        1.52   2.03   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.03   68       0.00   0.04   
     PW_Basis_Sup      gathers_scatterp            0.01   68       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.10   9        0.01   0.13   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   9        0.01   0.13   
     H_Hartree_pw      v_hartree                   0.01   9        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.03   89       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   89       0.00   0.01   
     PotXC             cal_v_eff                   0.09   9        0.01   0.12   
     XC_Functional     v_xc                        0.09   9        0.01   0.12   
     Potential         interpolate_vrs             0.00   9        0.00   0.00   
     Symmetry          rhog_symmetry               0.03   10       0.00   0.03   
     Symmetry          group fft grids             0.01   10       0.00   0.01   
     PSIInit           initialize_psi              1.49   1        1.49   1.99   
     Nonlocal          getvnl                      1.09   1710     0.00   1.46   
     pp_cell_vnl       getvnl                      1.09   1710     0.00   1.45   
     Structure_Factor  get_sk                      0.05   2052     0.00   0.07   
     DiagoIterAssist   diagH_subspace              9.46   1368     0.01   12.63  
     Operator          hPsi                        57.54  126891   0.00   76.81  
     Operator          EkineticPW                  0.39   126891   0.00   0.53   
     Operator          VeffPW                      51.13  126891   0.00   68.25  
     PW_Basis_K        recip2real                  29.30  177849   0.00   39.11  
     PW_Basis_K        gathers_scatterp            11.38  177849   0.00   15.19  
     PW_Basis_K        real2recip                  21.38  150147   0.00   28.54  
     PW_Basis_K        gatherp_scatters            6.77   150147   0.00   9.04   
     Operator          NonlocalPW                  5.78   126891   0.00   7.71   
     Nonlocal          add_nonlocal_pp             2.21   126891   0.00   2.95   
     DiagoIterAssist   diagH_LAPACK                0.26   1368     0.00   0.35   
     ESolver_KS_PW     hamilt2density_single       71.44  9        7.94   95.36  
     HSolverPW         solve                       71.40  9        7.93   95.31  
     DiagoCG           diag_once                   56.95  1539     0.04   76.02  
     DiagoCG_New       spsi_func                   0.52   251046   0.00   0.69   
     DiagoCG_New       hpsi_func                   49.00  125523   0.00   65.41  
     ElecStatePW       psiToRho                    4.93   9        0.55   6.58   
     Charge_Mixing     get_drho                    0.01   9        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   9        0.00   0.00   
     Charge            mix_rho                     0.01   7        0.00   0.01   
     Charge            Broyden_mixing              0.00   7        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   42       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.15   1        0.15   0.20   
     ModuleIO          write_rhog                  0.12   1        0.12   0.15   
     Forces            cal_force                   0.22   1        0.22   0.29   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.20   1        0.20   0.27   
     FS_Nonlocal_tools cal_becp                    0.58   171      0.00   0.78   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.85   1        0.85   1.13   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.07   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.77   1        0.77   1.02   
     ModuleIO          write_istate_info           0.04   1        0.04   0.05   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:06:44 2025
     FINISH Time  : Sat Apr 26 14:07:59 2025
     TOTAL  Time  : 75
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.015
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:08:00 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 32 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 32 * 30 * 30
     DONE(0.0862476  SEC) : SETUP UNITCELL
     DONE(0.163558   SEC) : SYMMETRY
     DONE(0.488504   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.50597    SEC) : INIT PLANEWAVE
     DONE(0.521151   SEC) : LOCAL POTENTIAL
     DONE(0.577479   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.3071
     DONE(0.577553   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.13391    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70776736e+03   0.00000000e+00   9.0454e-02  25.89
     CG2     -2.70842230e+03  -6.54946068e-01   2.9709e-02   4.96
     CG3     -2.70850933e+03  -8.70320130e-02   2.4144e-04   4.90
     CG4     -2.70852136e+03  -1.20272024e-02   1.6069e-05   9.62
     CG5     -2.70852153e+03  -1.63562342e-04   4.5820e-06   6.76
     CG6     -2.70852153e+03  -7.73910322e-06   4.0426e-07   5.02
     CG7     -2.70852154e+03  -2.44255698e-06   5.6110e-08   6.98
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
           -13.8211591080         0.0000000000        -0.0000000000 
            -0.0000000000        23.1995945629       -12.9123893108 
            -0.0000000000       -12.9123893108        32.6071924893 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 13.995209 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       67.47  17       3.97   100.00 
     Driver            reading                     0.06   1        0.06   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 67.41  1        67.41  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    66.86  1        66.86  99.09  
     ESolver_KS_PW     runner                      65.76  1        65.76  97.46  
     ESolver_KS_PW     before_scf                  1.56   1        1.56   2.31   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.10   8        0.01   0.14   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   8        0.01   0.14   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.05   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.12   
     XC_Functional     v_xc                        0.08   8        0.01   0.12   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.52   1        1.52   2.26   
     Nonlocal          getvnl                      0.98   1539     0.00   1.46   
     pp_cell_vnl       getvnl                      0.98   1539     0.00   1.45   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              7.99   1197     0.01   11.85  
     Operator          hPsi                        51.31  118741   0.00   76.04  
     Operator          EkineticPW                  0.37   118741   0.00   0.55   
     Operator          VeffPW                      45.29  118741   0.00   67.12  
     PW_Basis_K        recip2real                  25.87  163714   0.00   38.33  
     PW_Basis_K        gathers_scatterp            10.76  163714   0.00   15.95  
     PW_Basis_K        real2recip                  18.58  139090   0.00   27.54  
     PW_Basis_K        gatherp_scatters            6.34   139090   0.00   9.40   
     Operator          NonlocalPW                  5.42   118741   0.00   8.03   
     Nonlocal          add_nonlocal_pp             2.09   118741   0.00   3.10   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.34   
     ESolver_KS_PW     hamilt2density_single       64.02  8        8.00   94.89  
     HSolverPW         solve                       63.99  8        8.00   94.84  
     DiagoCG           diag_once                   51.74  1368     0.04   76.68  
     DiagoCG_New       spsi_func                   0.49   235088   0.00   0.73   
     DiagoCG_New       hpsi_func                   44.16  117544   0.00   65.45  
     ElecStatePW       psiToRho                    4.38   8        0.55   6.49   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.12   
     ModuleIO          write_rhog                  0.05   1        0.05   0.07   
     Forces            cal_force                   0.22   1        0.22   0.32   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.20   1        0.20   0.29   
     FS_Nonlocal_tools cal_becp                    0.59   171      0.00   0.87   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.87   1        0.87   1.28   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.78   1        0.78   1.16   
     ModuleIO          write_istate_info           0.03   1        0.03   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:08:00 2025
     FINISH Time  : Sat Apr 26 14:09:08 2025
     TOTAL  Time  : 68
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.016
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:09:09 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0857603  SEC) : SETUP UNITCELL
     DONE(0.160287   SEC) : SYMMETRY
     DONE(0.478991   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.496151   SEC) : INIT PLANEWAVE
     DONE(0.511424   SEC) : LOCAL POTENTIAL
     DONE(0.569006   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.4011
     DONE(0.569119   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.05705    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70774409e+03   0.00000000e+00   9.0208e-02  26.85
     CG2     -2.70842261e+03  -6.78527099e-01   2.9519e-02   5.08
     CG3     -2.70850950e+03  -8.68876556e-02   2.4141e-04   5.04
     CG4     -2.70852137e+03  -1.18629406e-02   1.5316e-05   9.83
     CG5     -2.70852152e+03  -1.57734088e-04   4.6974e-06   7.01
     CG6     -2.70852153e+03  -1.09237940e-05   1.9791e-07   5.16
     CG7     -2.70852154e+03  -1.57835023e-06   5.5219e-08   8.02
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            13.3600226269        15.9093943179        11.4030747388 
            15.9093943179        32.2254785007        -6.4957834319 
            11.4030747388        -6.4957834319        -4.7072920171 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 13.626070 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       70.31  17       4.14   100.00 
     Driver            reading                     0.06   1        0.06   0.08   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 70.25  1        70.25  99.92  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    69.69  1        69.69  99.11  
     ESolver_KS_PW     runner                      68.56  1        68.56  97.51  
     ESolver_KS_PW     before_scf                  1.49   1        1.49   2.12   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.01   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.09   8        0.01   0.13   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.13   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.01   
     PotXC             cal_v_eff                   0.08   8        0.01   0.12   
     XC_Functional     v_xc                        0.08   8        0.01   0.11   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.45   1        1.45   2.06   
     Nonlocal          getvnl                      0.98   1539     0.00   1.39   
     pp_cell_vnl       getvnl                      0.98   1539     0.00   1.39   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              8.15   1197     0.01   11.59  
     Operator          hPsi                        53.91  120747   0.00   76.67  
     Operator          EkineticPW                  0.38   120747   0.00   0.54   
     Operator          VeffPW                      47.94  120747   0.00   68.18  
     PW_Basis_K        recip2real                  27.53  165720   0.00   39.15  
     PW_Basis_K        gathers_scatterp            10.60  165720   0.00   15.07  
     PW_Basis_K        real2recip                  19.96  141096   0.00   28.39  
     PW_Basis_K        gatherp_scatters            6.36   141096   0.00   9.05   
     Operator          NonlocalPW                  5.36   120747   0.00   7.62   
     Nonlocal          add_nonlocal_pp             1.96   120747   0.00   2.79   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       66.91  8        8.36   95.15  
     HSolverPW         solve                       66.88  8        8.36   95.11  
     DiagoCG           diag_once                   54.30  1368     0.04   77.22  
     DiagoCG_New       spsi_func                   0.50   239100   0.00   0.70   
     DiagoCG_New       hpsi_func                   46.58  119550   0.00   66.25  
     ElecStatePW       psiToRho                    4.48   8        0.56   6.37   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.11   
     ModuleIO          write_rhog                  0.05   1        0.05   0.06   
     Forces            cal_force                   0.21   1        0.21   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.27   
     FS_Nonlocal_tools cal_becp                    0.60   171      0.00   0.86   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.88   1        0.88   1.26   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.06   1        0.06   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.80   1        0.80   1.14   
     ModuleIO          write_istate_info           0.05   1        0.05   0.07   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:09:09 2025
     FINISH Time  : Sat Apr 26 14:10:19 2025
     TOTAL  Time  : 70
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.017
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:10:20 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.106251   SEC) : SETUP UNITCELL
     DONE(0.182315   SEC) : SYMMETRY
     DONE(0.51804    SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.534464   SEC) : INIT PLANEWAVE
     DONE(0.549498   SEC) : LOCAL POTENTIAL
     DONE(0.60515    SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.448
     DONE(0.605249   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.11655    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777455e+03   0.00000000e+00   8.8109e-02  26.95
     CG2     -2.70842979e+03  -6.55239502e-01   2.7410e-02   5.10
     CG3     -2.70851175e+03  -8.19614634e-02   2.3043e-04   5.15
     CG4     -2.70852262e+03  -1.08744916e-02   1.3358e-05   9.88
     CG5     -2.70852276e+03  -1.37688932e-04   2.6583e-06   7.10
     CG6     -2.70852277e+03  -8.45961251e-06   7.9447e-08   5.70
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            17.5183508358         4.0805015017         4.3340993479 
             4.0805015017        22.2932484798        -2.4856112831 
             4.3340993479        -2.4856112831        11.6979708779 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 17.169857 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       63.29  17       3.72   100.00 
     Driver            reading                     0.08   1        0.08   0.13   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 63.20  1        63.20  99.87  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.12   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.01   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    62.64  1        62.64  98.98  
     ESolver_KS_PW     runner                      61.56  1        61.56  97.28  
     ESolver_KS_PW     before_scf                  1.51   1        1.51   2.39   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.04   2        0.02   0.06   
     PW_Basis_Sup      recip2real                  0.02   54       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   54       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.08   7        0.01   0.13   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.08   7        0.01   0.13   
     H_Hartree_pw      v_hartree                   0.01   7        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   71       0.00   0.03   
     PW_Basis_Sup      gatherp_scatters            0.01   71       0.00   0.01   
     PotXC             cal_v_eff                   0.07   7        0.01   0.12   
     XC_Functional     v_xc                        0.07   7        0.01   0.12   
     Potential         interpolate_vrs             0.00   7        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   8        0.00   0.03   
     Symmetry          group fft grids             0.01   8        0.00   0.01   
     PSIInit           initialize_psi              1.48   1        1.48   2.34   
     Nonlocal          getvnl                      0.87   1368     0.00   1.37   
     pp_cell_vnl       getvnl                      0.87   1368     0.00   1.37   
     Structure_Factor  get_sk                      0.04   1710     0.00   0.07   
     DiagoIterAssist   diagH_subspace              7.07   1026     0.01   11.17  
     Operator          hPsi                        48.31  108205   0.00   76.33  
     Operator          EkineticPW                  0.33   108205   0.00   0.52   
     Operator          VeffPW                      42.93  108205   0.00   67.83  
     PW_Basis_K        recip2real                  24.65  147193   0.00   38.95  
     PW_Basis_K        gathers_scatterp            9.54   147193   0.00   15.07  
     PW_Basis_K        real2recip                  17.91  125647   0.00   28.30  
     PW_Basis_K        gatherp_scatters            5.61   125647   0.00   8.86   
     Operator          NonlocalPW                  4.85   108205   0.00   7.66   
     Nonlocal          add_nonlocal_pp             1.84   108205   0.00   2.91   
     DiagoIterAssist   diagH_LAPACK                0.19   1026     0.00   0.31   
     ESolver_KS_PW     hamilt2density_single       59.80  7        8.54   94.49  
     HSolverPW         solve                       59.77  7        8.54   94.45  
     DiagoCG           diag_once                   48.85  1197     0.04   77.19  
     DiagoCG_New       spsi_func                   0.44   214358   0.00   0.70   
     DiagoCG_New       hpsi_func                   41.96  107179   0.00   66.30  
     ElecStatePW       psiToRho                    4.05   7        0.58   6.41   
     Charge_Mixing     get_drho                    0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   7        0.00   0.00   
     Charge            mix_rho                     0.00   5        0.00   0.01   
     Charge            Broyden_mixing              0.00   5        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   20       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.17   1        0.17   0.27   
     ModuleIO          write_rhog                  0.13   1        0.13   0.21   
     Forces            cal_force                   0.21   1        0.21   0.34   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.31   
     FS_Nonlocal_tools cal_becp                    0.58   171      0.00   0.92   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.85   1        0.85   1.34   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.76   1        0.76   1.21   
     ModuleIO          write_istate_info           0.03   1        0.03   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:10:20 2025
     FINISH Time  : Sat Apr 26 14:11:23 2025
     TOTAL  Time  : 63
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.018
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:11:24 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0698768  SEC) : SETUP UNITCELL
     DONE(0.142617   SEC) : SYMMETRY
     DONE(0.461878   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.47914    SEC) : INIT PLANEWAVE
     DONE(0.494746   SEC) : LOCAL POTENTIAL
     DONE(0.552508   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.6829
     DONE(0.55258    SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.06137    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70773898e+03   0.00000000e+00   8.9592e-02  26.45
     CG2     -2.70842633e+03  -6.87347034e-01   2.9073e-02   5.12
     CG3     -2.70851111e+03  -8.47845012e-02   2.5703e-04   5.08
     CG4     -2.70852258e+03  -1.14690959e-02   1.7250e-05   9.75
     CG5     -2.70852276e+03  -1.77369623e-04   4.5944e-06   6.88
     CG6     -2.70852277e+03  -9.21834257e-06   2.8135e-07   5.23
     CG7     -2.70852277e+03  -2.16976144e-06   5.3446e-08   7.62
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            13.1656567943        -5.9659368294        -2.8359375815 
            -5.9659368294         6.3681827222         1.6482449758 
            -2.8359375815         1.6482449758        18.9874530103 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 12.840431 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       69.43  17       4.08   100.00 
     Driver            reading                     0.04   1        0.04   0.06   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 69.38  1        69.38  99.94  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.10   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    68.84  1        68.84  99.15  
     ESolver_KS_PW     runner                      67.71  1        67.71  97.53  
     ESolver_KS_PW     before_scf                  1.51   1        1.51   2.17   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.04   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.09   8        0.01   0.13   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.13   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.08   8        0.01   0.12   
     XC_Functional     v_xc                        0.08   8        0.01   0.12   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.47   1        1.47   2.12   
     Nonlocal          getvnl                      1.00   1539     0.00   1.45   
     pp_cell_vnl       getvnl                      1.00   1539     0.00   1.44   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              8.15   1197     0.01   11.74  
     Operator          hPsi                        53.25  119938   0.00   76.69  
     Operator          EkineticPW                  0.38   119938   0.00   0.55   
     Operator          VeffPW                      47.29  119938   0.00   68.12  
     PW_Basis_K        recip2real                  27.02  164911   0.00   38.92  
     PW_Basis_K        gathers_scatterp            10.57  164911   0.00   15.23  
     PW_Basis_K        real2recip                  19.78  140287   0.00   28.49  
     PW_Basis_K        gatherp_scatters            6.37   140287   0.00   9.18   
     Operator          NonlocalPW                  5.35   119938   0.00   7.70   
     Nonlocal          add_nonlocal_pp             1.96   119938   0.00   2.82   
     DiagoIterAssist   diagH_LAPACK                0.22   1197     0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       66.03  8        8.25   95.11  
     HSolverPW         solve                       66.00  8        8.25   95.07  
     DiagoCG           diag_once                   53.39  1368     0.04   76.90  
     DiagoCG_New       spsi_func                   0.48   237482   0.00   0.69   
     DiagoCG_New       hpsi_func                   45.91  118741   0.00   66.13  
     ElecStatePW       psiToRho                    4.50   8        0.56   6.48   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.11   
     ModuleIO          write_rhog                  0.05   1        0.05   0.07   
     Forces            cal_force                   0.21   1        0.21   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.28   
     FS_Nonlocal_tools cal_becp                    0.60   171      0.00   0.86   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.89   1        0.89   1.28   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.81   1        0.81   1.16   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:11:24 2025
     FINISH Time  : Sat Apr 26 14:12:34 2025
     TOTAL  Time  : 70
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.019
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:12:35 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 32 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 32 * 30
     DONE(0.0860983  SEC) : SETUP UNITCELL
     DONE(0.161782   SEC) : SYMMETRY
     DONE(0.473677   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.491569   SEC) : INIT PLANEWAVE
     DONE(0.506932   SEC) : LOCAL POTENTIAL
     DONE(0.567981   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.495
     DONE(0.568049   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.04827    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70776554e+03   0.00000000e+00   8.9504e-02  26.37
     CG2     -2.70842371e+03  -6.58177042e-01   2.8751e-02   5.05
     CG3     -2.70851041e+03  -8.66939413e-02   2.1965e-04   5.04
     CG4     -2.70852139e+03  -1.09804343e-02   1.3874e-05   9.86
     CG5     -2.70852153e+03  -1.39631010e-04   2.2465e-06   6.82
     CG6     -2.70852153e+03  -7.54911347e-06   1.1738e-07   5.64
     CG7     -2.70852154e+03  -7.44770677e-07   3.3470e-08   7.46
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            13.8883252020       -15.9291404097       -11.0702452741 
           -15.9291404097        -4.0193906173         6.4766316396 
           -11.0702452741         6.4766316396        32.5520428753 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 14.140326 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       69.54  17       4.09   100.00 
     Driver            reading                     0.06   1        0.06   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 69.48  1        69.48  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.01   1        0.01   0.01   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    68.94  1        68.94  99.13  
     ESolver_KS_PW     runner                      67.81  1        67.81  97.51  
     ESolver_KS_PW     before_scf                  1.48   1        1.48   2.13   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.10   8        0.01   0.14   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.10   8        0.01   0.14   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.02   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.05   
     PW_Basis_Sup      gatherp_scatters            0.02   80       0.00   0.02   
     PotXC             cal_v_eff                   0.09   8        0.01   0.12   
     XC_Functional     v_xc                        0.09   8        0.01   0.12   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.45   1        1.45   2.08   
     Nonlocal          getvnl                      0.99   1539     0.00   1.43   
     pp_cell_vnl       getvnl                      0.99   1539     0.00   1.42   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              8.12   1197     0.01   11.68  
     Operator          hPsi                        53.20  121607   0.00   76.50  
     Operator          EkineticPW                  0.39   121607   0.00   0.56   
     Operator          VeffPW                      47.15  121607   0.00   67.80  
     PW_Basis_K        recip2real                  26.93  166580   0.00   38.72  
     PW_Basis_K        gathers_scatterp            10.99  166580   0.00   15.80  
     PW_Basis_K        real2recip                  19.12  141956   0.00   27.50  
     PW_Basis_K        gatherp_scatters            6.31   141956   0.00   9.07   
     Operator          NonlocalPW                  5.43   121607   0.00   7.80   
     Nonlocal          add_nonlocal_pp             2.08   121607   0.00   3.00   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.33   
     ESolver_KS_PW     hamilt2density_single       66.15  8        8.27   95.12  
     HSolverPW         solve                       66.12  8        8.26   95.07  
     DiagoCG           diag_once                   53.67  1368     0.04   77.17  
     DiagoCG_New       spsi_func                   0.50   240820   0.00   0.72   
     DiagoCG_New       hpsi_func                   45.91  120410   0.00   66.01  
     ElecStatePW       psiToRho                    4.35   8        0.54   6.25   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.12   
     ModuleIO          write_rhog                  0.05   1        0.05   0.08   
     Forces            cal_force                   0.21   1        0.21   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.27   
     FS_Nonlocal_tools cal_becp                    0.62   171      0.00   0.89   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.90   1        0.90   1.29   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.81   1        0.81   1.17   
     ModuleIO          write_istate_info           0.02   1        0.02   0.03   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:12:35 2025
     FINISH Time  : Sat Apr 26 14:13:44 2025
     TOTAL  Time  : 69
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.020
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:13:45 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0982801  SEC) : SETUP UNITCELL
     DONE(0.176136   SEC) : SYMMETRY
     DONE(0.502467   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.521621   SEC) : INIT PLANEWAVE
     DONE(0.537192   SEC) : LOCAL POTENTIAL
     DONE(0.593918   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.6829
     DONE(0.593995   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.12735    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70777112e+03   0.00000000e+00   8.8454e-02  26.96
     CG2     -2.70842567e+03  -6.54553067e-01   2.7917e-02   5.17
     CG3     -2.70850965e+03  -8.39729555e-02   2.2939e-04   5.10
     CG4     -2.70852138e+03  -1.17362276e-02   1.3406e-05   9.81
     CG5     -2.70852152e+03  -1.42675086e-04   3.1050e-06   7.01
     CG6     -2.70852153e+03  -1.00982840e-05   9.4744e-08   5.55
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            17.1857799440        15.8613770647        10.8251767588 
            15.8613770647        -0.6457555734         6.3332546447 
            10.8251767588         6.3332546447        35.6029903852 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 17.381005 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       62.95  17       3.70   100.00 
     Driver            reading                     0.07   1        0.07   0.12   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 62.88  1        62.88  99.88  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.12   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.02   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.09   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    62.32  1        62.32  99.00  
     ESolver_KS_PW     runner                      61.21  1        61.21  97.23  
     ESolver_KS_PW     before_scf                  1.53   1        1.53   2.44   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   54       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   54       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.08   7        0.01   0.12   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.08   7        0.01   0.12   
     H_Hartree_pw      v_hartree                   0.01   7        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   71       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   71       0.00   0.01   
     PotXC             cal_v_eff                   0.07   7        0.01   0.11   
     XC_Functional     v_xc                        0.07   7        0.01   0.11   
     Potential         interpolate_vrs             0.00   7        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   8        0.00   0.03   
     Symmetry          group fft grids             0.01   8        0.00   0.01   
     PSIInit           initialize_psi              1.50   1        1.50   2.39   
     Nonlocal          getvnl                      0.90   1368     0.00   1.44   
     pp_cell_vnl       getvnl                      0.90   1368     0.00   1.43   
     Structure_Factor  get_sk                      0.05   1710     0.00   0.08   
     DiagoIterAssist   diagH_subspace              7.07   1026     0.01   11.23  
     Operator          hPsi                        48.17  107994   0.00   76.52  
     Operator          EkineticPW                  0.33   107994   0.00   0.52   
     Operator          VeffPW                      42.70  107994   0.00   67.82  
     PW_Basis_K        recip2real                  24.66  146982   0.00   39.17  
     PW_Basis_K        gathers_scatterp            9.54   146982   0.00   15.15  
     PW_Basis_K        real2recip                  17.56  125436   0.00   27.90  
     PW_Basis_K        gatherp_scatters            5.54   125436   0.00   8.79   
     Operator          NonlocalPW                  4.95   107994   0.00   7.86   
     Nonlocal          add_nonlocal_pp             1.90   107994   0.00   3.02   
     DiagoIterAssist   diagH_LAPACK                0.19   1026     0.00   0.31   
     ESolver_KS_PW     hamilt2density_single       59.52  7        8.50   94.55  
     HSolverPW         solve                       59.49  7        8.50   94.51  
     DiagoCG           diag_once                   48.70  1197     0.04   77.36  
     DiagoCG_New       spsi_func                   0.45   213936   0.00   0.71   
     DiagoCG_New       hpsi_func                   41.82  106968   0.00   66.44  
     ElecStatePW       psiToRho                    3.93   7        0.56   6.25   
     Charge_Mixing     get_drho                    0.00   7        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   7        0.00   0.00   
     Charge            mix_rho                     0.00   5        0.00   0.01   
     Charge            Broyden_mixing              0.00   5        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   20       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.07   1        0.07   0.12   
     ModuleIO          write_rhog                  0.04   1        0.04   0.07   
     Forces            cal_force                   0.21   1        0.21   0.34   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.31   
     FS_Nonlocal_tools cal_becp                    0.60   171      0.00   0.95   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.88   1        0.88   1.40   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.80   1        0.80   1.27   
     ModuleIO          write_istate_info           0.03   1        0.03   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:13:45 2025
     FINISH Time  : Sat Apr 26 14:14:48 2025
     TOTAL  Time  : 63
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.021
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:14:49 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0744831  SEC) : SETUP UNITCELL
     DONE(0.151667   SEC) : SYMMETRY
     DONE(0.4734     SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.490718   SEC) : INIT PLANEWAVE
     DONE(0.506342   SEC) : LOCAL POTENTIAL
     DONE(0.563776   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.5889
     DONE(0.563843   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.03753    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70775219e+03   0.00000000e+00   8.9782e-02  26.90
     CG2     -2.70842629e+03  -6.74097391e-01   2.8766e-02   5.14
     CG3     -2.70851138e+03  -8.50839396e-02   2.2044e-04   5.16
     CG4     -2.70852262e+03  -1.12427085e-02   1.4976e-05  10.00
     CG5     -2.70852276e+03  -1.41317321e-04   2.8725e-06   6.87
     CG6     -2.70852277e+03  -7.48142058e-06   1.9367e-07   5.63
     CG7     -2.70852277e+03  -1.26583161e-06   3.1669e-08   7.34
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            13.7348820041         5.9620395520         2.8747916435 
             5.9620395520         6.9418484151         1.6708269300 
             2.8747916435         1.6708269300        19.5866698009 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 13.421133 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       70.35  17       4.14   100.00 
     Driver            reading                     0.05   1        0.05   0.07   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 70.30  1        70.30  99.93  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.08   1        0.08   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.01   1        0.01   0.02   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    69.75  1        69.75  99.15  
     ESolver_KS_PW     runner                      68.59  1        68.59  97.50  
     ESolver_KS_PW     before_scf                  1.47   1        1.47   2.09   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.09   8        0.01   0.12   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.12   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   80       0.00   0.03   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.01   
     PotXC             cal_v_eff                   0.08   8        0.01   0.11   
     XC_Functional     v_xc                        0.08   8        0.01   0.11   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.44   1        1.44   2.05   
     Nonlocal          getvnl                      0.99   1539     0.00   1.41   
     pp_cell_vnl       getvnl                      0.99   1539     0.00   1.41   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              8.08   1197     0.01   11.49  
     Operator          hPsi                        54.02  120512   0.00   76.79  
     Operator          EkineticPW                  0.38   120512   0.00   0.54   
     Operator          VeffPW                      47.86  120512   0.00   68.03  
     PW_Basis_K        recip2real                  27.38  165485   0.00   38.93  
     PW_Basis_K        gathers_scatterp            10.64  165485   0.00   15.12  
     PW_Basis_K        real2recip                  20.02  140861   0.00   28.45  
     PW_Basis_K        gatherp_scatters            6.41   140861   0.00   9.11   
     Operator          NonlocalPW                  5.56   120512   0.00   7.90   
     Nonlocal          add_nonlocal_pp             2.14   120512   0.00   3.05   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       66.95  8        8.37   95.17  
     HSolverPW         solve                       66.92  8        8.37   95.13  
     DiagoCG           diag_once                   54.38  1368     0.04   77.30  
     DiagoCG_New       spsi_func                   0.49   238630   0.00   0.69   
     DiagoCG_New       hpsi_func                   46.74  119315   0.00   66.44  
     ElecStatePW       psiToRho                    4.48   8        0.56   6.37   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.11   
     ModuleIO          write_rhog                  0.05   1        0.05   0.07   
     Forces            cal_force                   0.24   1        0.24   0.34   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.22   1        0.22   0.32   
     FS_Nonlocal_tools cal_becp                    0.63   171      0.00   0.89   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.89   1        0.89   1.27   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.81   1        0.81   1.15   
     ModuleIO          write_istate_info           0.02   1        0.02   0.03   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:14:49 2025
     FINISH Time  : Sat Apr 26 14:16:00 2025
     TOTAL  Time  : 71
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.022
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:16:01 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 30
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 30
     DONE(0.0920556  SEC) : SETUP UNITCELL
     DONE(0.166858   SEC) : SYMMETRY
     DONE(0.506581   SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.529922   SEC) : INIT PLANEWAVE
     DONE(0.551335   SEC) : LOCAL POTENTIAL
     DONE(0.610603   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.6359
     DONE(0.610698   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.12734    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70771439e+03   0.00000000e+00   9.0997e-02  26.80
     CG2     -2.70842324e+03  -7.08853044e-01   3.0301e-02   5.18
     CG3     -2.70851079e+03  -8.75465061e-02   2.6811e-04   5.03
     CG4     -2.70852258e+03  -1.17845655e-02   1.7660e-05   9.77
     CG5     -2.70852275e+03  -1.78810990e-04   4.7603e-06   6.92
     CG6     -2.70852277e+03  -1.30086344e-05   1.1760e-07   5.31
     CG7     -2.70852277e+03  -1.50627877e-06   7.5846e-08   8.62
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            12.3305814299        -4.5377716687        -4.5242776667 
            -4.5377716687        17.6405648330        -2.5946787818 
            -4.5242776667        -2.5946787818         6.0952465471 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 12.022131 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       70.97  17       4.17   100.00 
     Driver            reading                     0.07   1        0.07   0.10   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 70.91  1        70.91  99.90  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.02   1        0.02   0.03   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.06   1        0.06   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    70.33  1        70.33  99.09  
     ESolver_KS_PW     runner                      69.23  1        69.23  97.54  
     ESolver_KS_PW     before_scf                  1.52   1        1.52   2.14   
     H_Ewald_pw        compute_ewald               0.02   1        0.02   0.02   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.04   2        0.02   0.06   
     PW_Basis_Sup      recip2real                  0.03   61       0.00   0.04   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.02   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.09   8        0.01   0.12   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.12   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.02   80       0.00   0.03   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.01   
     PotXC             cal_v_eff                   0.08   8        0.01   0.11   
     XC_Functional     v_xc                        0.08   8        0.01   0.11   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.04   
     Symmetry          group fft grids             0.01   9        0.00   0.02   
     PSIInit           initialize_psi              1.46   1        1.46   2.05   
     Nonlocal          getvnl                      0.99   1539     0.00   1.39   
     pp_cell_vnl       getvnl                      0.99   1539     0.00   1.39   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              8.19   1197     0.01   11.54  
     Operator          hPsi                        54.32  122134   0.00   76.53  
     Operator          EkineticPW                  0.37   122134   0.00   0.52   
     Operator          VeffPW                      48.21  122134   0.00   67.92  
     PW_Basis_K        recip2real                  27.63  167107   0.00   38.93  
     PW_Basis_K        gathers_scatterp            10.82  167107   0.00   15.24  
     PW_Basis_K        real2recip                  20.27  142483   0.00   28.56  
     PW_Basis_K        gatherp_scatters            6.43   142483   0.00   9.07   
     Operator          NonlocalPW                  5.50   122134   0.00   7.76   
     Nonlocal          add_nonlocal_pp             2.07   122134   0.00   2.92   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.32   
     ESolver_KS_PW     hamilt2density_single       67.54  8        8.44   95.16  
     HSolverPW         solve                       67.51  8        8.44   95.12  
     DiagoCG           diag_once                   54.68  1368     0.04   77.04  
     DiagoCG_New       spsi_func                   0.50   241874   0.00   0.70   
     DiagoCG_New       hpsi_func                   46.94  120937   0.00   66.14  
     ElecStatePW       psiToRho                    4.69   8        0.59   6.60   
     Charge_Mixing     get_drho                    0.01   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.11   
     ModuleIO          write_rhog                  0.05   1        0.05   0.07   
     Forces            cal_force                   0.21   1        0.21   0.30   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.27   
     FS_Nonlocal_tools cal_becp                    0.58   171      0.00   0.82   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.02   
     Stress_PW         cal_stress                  0.87   1        0.87   1.22   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.07   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.79   1        0.79   1.11   
     ModuleIO          write_istate_info           0.02   1        0.02   0.03   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:16:01 2025
     FINISH Time  : Sat Apr 26 14:17:12 2025
     TOTAL  Time  : 71
     SEE INFORMATION IN : OUT.Mn/
    /personal/elasticity/task.023
                                                                                         
                                  ABACUS v3.9.0
    
                   Atomic-orbital Based Ab-initio Computation at UStc                    
    
                         Website: http://abacus.ustc.edu.cn/                             
                   Documentation: https://abacus.deepmodeling.com/                       
                      Repository: https://github.com/abacusmodeling/abacus-develop       
                                  https://github.com/deepmodeling/abacus-develop         
                          Commit: 68735ed (Fri Dec 27 15:05:38 2024 +0800)
    
     Sat Apr 26 14:17:13 2025
     MAKE THE DIR         : OUT.Mn/
     RUNNING WITH DEVICE  : CPU / Intel(R) Xeon(R) Platinum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Warning: the number of valence electrons in pseudopotential > 7 for Mn: [Ar] 3d5 4s2
     Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient.
     If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     UNIFORM GRID DIM        : 30 * 30 * 32
     UNIFORM GRID DIM(BIG)   : 30 * 30 * 32
     DONE(0.0881077  SEC) : SETUP UNITCELL
     DONE(0.161483   SEC) : SYMMETRY
     DONE(0.4844     SEC) : INIT K-POINTS
     ---------------------------------------------------------
     Ion relaxation calculations
     ---------------------------------------------------------
     SPIN    KPOINTS         PROCESSORS  THREADS     
     1       171             2           2           
     ---------------------------------------------------------
     Use plane wave basis
     ---------------------------------------------------------
     ELEMENT NATOM       XC          
     Mn      1           
     ---------------------------------------------------------
     Initial plane wave basis and FFT box
     ---------------------------------------------------------
     DONE(0.500976   SEC) : INIT PLANEWAVE
     DONE(0.518972   SEC) : LOCAL POTENTIAL
     DONE(0.575749   SEC) : NON-LOCAL POTENTIAL
     MEMORY FOR PSI (MB)  : 29.4011
     DONE(0.575815   SEC) : INIT BASIS
     -------------------------------------------
     STEP OF RELAXATION : 1
     -------------------------------------------
     START CHARGE      : atomic
     DONE(2.06802    SEC) : INIT SCF
     ITER       ETOT/eV          EDIFF/eV         DRHO     TIME/s
     CG1     -2.70776648e+03   0.00000000e+00   8.8033e-02  25.72
     CG2     -2.70842685e+03  -6.60371305e-01   2.7501e-02   4.89
     CG3     -2.70851041e+03  -8.35608282e-02   1.9941e-04   4.99
     CG4     -2.70852139e+03  -1.09787617e-02   1.3232e-05   9.45
     CG5     -2.70852153e+03  -1.36087234e-04   2.3432e-06   6.71
     CG6     -2.70852153e+03  -6.50756700e-06   1.5074e-07   5.40
     CG7     -2.70852154e+03  -9.60741834e-07   1.2024e-08   7.08
    ----------------------------------------------------------------
     TOTAL-STRESS (KBAR)                                            
    ----------------------------------------------------------------
            14.7787559175       -15.9030022123       -11.2440506757 
           -15.9030022123        33.6366319937        -6.4051950688 
           -11.2440506757        -6.4051950688        -3.1593432681 
    ----------------------------------------------------------------
     TOTAL-PRESSURE: 15.085348 KBAR
    
     ETOT DIFF (eV)       : 0.000000
     LARGEST GRAD (eV/A)  : 0.000000
    TIME STATISTICS
    -----------------------------------------------------------------------------
        CLASS_NAME                NAME             TIME/s  CALLS   AVG/s  PER/%  
    -----------------------------------------------------------------------------
                       total                       67.53  17       3.97   100.00 
     Driver            reading                     0.06   1        0.06   0.09   
     Input_Conv        Convert                     0.00   1        0.00   0.00   
     Driver            driver_line                 67.47  1        67.47  99.91  
     UnitCell          check_tau                   0.00   1        0.00   0.00   
     PW_Basis_Sup      setuptransform              0.00   1        0.00   0.00   
     PW_Basis_Sup      distributeg                 0.00   1        0.00   0.00   
     mymath            heapsort                    0.00   3        0.00   0.00   
     Charge_Mixing     init_mixing                 0.00   2        0.00   0.00   
     Symmetry          analy_sys                   0.07   1        0.07   0.11   
     PW_Basis_K        setuptransform              0.01   1        0.01   0.01   
     PW_Basis_K        distributeg                 0.00   1        0.00   0.00   
     PW_Basis          setup_struc_factor          0.00   1        0.00   0.00   
     ppcell_vl         init_vloc                   0.02   1        0.02   0.03   
     ppcell_vnl        init                        0.00   1        0.00   0.00   
     ppcell_vnl        init_vnl                    0.05   1        0.05   0.08   
     WF_atomic         init_at_1                   0.00   1        0.00   0.00   
     wavefunc          wfcinit                     0.00   1        0.00   0.00   
     Ions              opt_ions                    66.92  1        66.92  99.09  
     ESolver_KS_PW     runner                      65.83  1        65.83  97.47  
     ESolver_KS_PW     before_scf                  1.49   1        1.49   2.21   
     H_Ewald_pw        compute_ewald               0.00   1        0.00   0.00   
     Charge            set_rho_core                0.00   1        0.00   0.00   
     Charge            atomic_rho                  0.03   2        0.02   0.05   
     PW_Basis_Sup      recip2real                  0.02   61       0.00   0.03   
     PW_Basis_Sup      gathers_scatterp            0.01   61       0.00   0.01   
     Potential         init_pot                    0.01   1        0.01   0.02   
     Potential         update_from_charge          0.10   8        0.01   0.14   
     Potential         cal_fixed_v                 0.00   1        0.00   0.00   
     PotLocal          cal_fixed_v                 0.00   1        0.00   0.00   
     Potential         cal_v_eff                   0.09   8        0.01   0.14   
     H_Hartree_pw      v_hartree                   0.01   8        0.00   0.01   
     PW_Basis_Sup      real2recip                  0.03   80       0.00   0.04   
     PW_Basis_Sup      gatherp_scatters            0.01   80       0.00   0.02   
     PotXC             cal_v_eff                   0.09   8        0.01   0.13   
     XC_Functional     v_xc                        0.09   8        0.01   0.13   
     Potential         interpolate_vrs             0.00   8        0.00   0.00   
     Symmetry          rhog_symmetry               0.02   9        0.00   0.03   
     Symmetry          group fft grids             0.01   9        0.00   0.01   
     PSIInit           initialize_psi              1.46   1        1.46   2.16   
     Nonlocal          getvnl                      0.99   1539     0.00   1.47   
     pp_cell_vnl       getvnl                      0.99   1539     0.00   1.46   
     Structure_Factor  get_sk                      0.05   1881     0.00   0.07   
     DiagoIterAssist   diagH_subspace              7.77   1197     0.01   11.51  
     Operator          hPsi                        51.26  120969   0.00   75.90  
     Operator          EkineticPW                  0.38   120969   0.00   0.56   
     Operator          VeffPW                      45.17  120969   0.00   66.89  
     PW_Basis_K        recip2real                  25.91  165942   0.00   38.37  
     PW_Basis_K        gathers_scatterp            11.03  165942   0.00   16.33  
     PW_Basis_K        real2recip                  18.36  141318   0.00   27.19  
     PW_Basis_K        gatherp_scatters            6.38   141318   0.00   9.45   
     Operator          NonlocalPW                  5.48   120969   0.00   8.11   
     Nonlocal          add_nonlocal_pp             2.11   120969   0.00   3.13   
     DiagoIterAssist   diagH_LAPACK                0.23   1197     0.00   0.34   
     ESolver_KS_PW     hamilt2density_single       64.16  8        8.02   95.00  
     HSolverPW         solve                       64.13  8        8.02   94.95  
     DiagoCG           diag_once                   52.04  1368     0.04   77.06  
     DiagoCG_New       spsi_func                   0.50   239544   0.00   0.75   
     DiagoCG_New       hpsi_func                   44.31  119772   0.00   65.61  
     ElecStatePW       psiToRho                    4.35   8        0.54   6.44   
     Charge_Mixing     get_drho                    0.00   8        0.00   0.01   
     Charge_Mixing     inner_product_recip_rho     0.00   8        0.00   0.00   
     Charge            mix_rho                     0.01   6        0.00   0.01   
     Charge            Broyden_mixing              0.00   6        0.00   0.00   
     Charge_Mixing     inner_product_recip_hartree 0.00   30       0.00   0.00   
     ESolver_KS_PW     after_scf                   0.08   1        0.08   0.12   
     ModuleIO          write_rhog                  0.05   1        0.05   0.08   
     Forces            cal_force                   0.21   1        0.21   0.31   
     Forces            cal_force_loc               0.00   1        0.00   0.00   
     Forces            cal_force_ew                0.00   1        0.00   0.00   
     Forces            cal_force_nl                0.19   1        0.19   0.29   
     FS_Nonlocal_tools cal_becp                    0.59   171      0.00   0.87   
     Forces            cal_force_cc                0.00   1        0.00   0.00   
     Forces            cal_force_scc               0.02   1        0.02   0.03   
     Stress_PW         cal_stress                  0.86   1        0.86   1.28   
     Stress_Func       stress_kin                  0.02   1        0.02   0.03   
     Stress_Func       stress_har                  0.00   1        0.00   0.00   
     Stress_Func       stress_ewa                  0.00   1        0.00   0.00   
     Stress_Func       stress_gga                  0.01   1        0.01   0.01   
     Stress_Func       stress_loc                  0.05   1        0.05   0.08   
     Stress_Func       stress_cc                   0.00   1        0.00   0.00   
     Stress_Func       stress_nl                   0.78   1        0.78   1.16   
     ModuleIO          write_istate_info           0.02   1        0.02   0.04   
    -----------------------------------------------------------------------------
    
    
     START  Time  : Sat Apr 26 14:17:13 2025
     FINISH Time  : Sat Apr 26 14:18:20 2025
     TOTAL  Time  : 67
     SEE INFORMATION IN : OUT.Mn/
    


```
#计算弹性常数
!python compute_dfm.py abacus
```

    # Elastic Constants in GPa
     333.71  280.04  306.99  -38.05  -65.85   46.54 
     304.79  250.59  277.81  -38.43   66.51  -47.00 
     222.24  330.72  276.61   76.80    0.00    0.00 
       6.22  -79.26   91.02   29.74   51.54   74.03 
       3.30   88.41  -81.81  -30.08   52.12   73.72 
     139.99  -31.28  -73.87   60.24    0.00    0.00 
    # Bulk   Modulus BV = 287.05 GPa
    # Shear  Modulus GV = 16.11 GPa
    # Youngs Modulus EV = 47.44 GPa
    # Poission Ratio uV = 0.47 
    
