# 第一章：计算准备与理论概览

这是一篇关于 ABACUS 结合 ShengBTE 计算晶格热导率的实战教程第一章。本章将为您构建坚实的理论框架和计算环境，这是后续所有微扰计算的基石。

---

# 第一章：计算准备与理论概览

在涉足晶格热导率（Lattice Thermal Conductivity, $\kappa_L$）计算之前，我们必须达成一个共识：**热导率计算属于高阶微扰计算，对基态结构的精度极其敏感。** 任何在基态结构优化（Relax）阶段残留的微小受力，都会在后续的三阶力常数计算中被放大为巨大的数值噪音，导致出现虚频或完全错误的热导率数值。

本章将带您完成工具链的组装，并执行一个达到“科研级精度”的结构优化。

## 1.1 工具链环境配置与工作流逻辑

ABACUS 计算热导率并非单打独斗，而是一个多软件协作的接力赛。理解数据在不同软件间的流转至关重要。

### 1.1.1 核心工作流图解

我们将整个计算流程划分为三个独立的阶段，请务必在您的工作目录下建立相应的文件夹结构（如 `2nd_order`, `3rd_order`, `shengbte`）：

1.  **阶段一：二阶谐波性质 (Harmonic Properties)**
    *   **核心任务**：计算声子谱（Phonon Dispersion）和二阶力常数（2nd IFCs）。
    *   **工具链**：`ABACUS` $\leftrightarrow$ `Phonopy`
    *   **关键桥梁**：
        *   `phonopy --abacus`: 生成超胞位移。
        *   `au2si.py`: **[关键]** 将 ABACUS 的原子单位制力常数转换为 ShengBTE 需要的 eV/Å² 单位。

2.  **阶段二：三阶非谐性质 (Anharmonic Properties)**
    *   **核心任务**：计算三阶力常数（3rd IFCs），这是热阻的主要来源。
    *   **工具链**：`ABACUS` $\leftrightarrow$ `thirdorder.py` (ShengBTE套件)
    *   **关键桥梁**：
        *   `pos2stru.py`: 将 `thirdorder` 生成的 POSCAR 转换为 ABACUS 的 STRU 文件。
        *   `aba2vasp.py`: **[关键]** 将 ABACUS 的输出伪装成 `vasprun.xml`，以便 `thirdorder` 读取受力。

3.  **阶段三：玻尔兹曼输运方程求解 (BTE Solver)**
    *   **核心任务**：汇总二阶和三阶信息，计算散射率和热导率。
    *   **工具链**：`ShengBTE`
    *   **输入**：`CONTROL`, `FORCE_CONSTANTS_2ND`, `FORCE_CONSTANTS_3RD`。

### 1.1.2 必备辅助脚本 (The "Missing Links")

ABACUS 与 ShengBTE/thirdorder 的原生接口主要依赖社区提供的 Python 脚本进行格式转换。**请注意，以下脚本并非 ABACUS 二进制程序内置命令**，您通常需要在 ABACUS 官方案例库（`examples/interface_ShengBTE`）中获取它们，并放置于您的工作目录或系统 PATH 中：

*   **`pos2stru.py`**: 依赖 `ASE` 库。用于将 VASP 格式结构转换为 ABACUS STRU 格式。
*   **`aba2vasp.py`**: 用于解析 ABACUS 的输出文件，生成假的 `vasprun.xml`。
*   **`au2si.py`**: 用于二阶力常数的单位转换。

> **环境检查清单**：
> - [ ] ABACUS (v3.2.0 或更高版本，推荐 PW 基组支持)
> - [ ] Phonopy
> - [ ] ShengBTE & thirdorder.py
> - [ ] Python 环境 (需安装 `ase`, `numpy`)
> - [ ] 上述 3 个辅助脚本已就位

---

## 1.2 高精度结构优化 (Relax)

一切计算始于完美的晶体结构。对于热导率计算，普通的结构优化精度是不够的。我们需要消除所有残余应力，确保原子处于势能面的绝对极小值点。

### 1.2.1 输入文件编写 (INPUT)

以下是一个用于高精度结构优化的 `INPUT` 模板（以平面波 PW 基组为例）。

**文件名**: `INPUT`

```bash
INPUT_PARAMETERS
# 1. 基础参数
suffix          Si_relax      # 输出文件后缀
calculation     relax         # 计算类型：结构优化
basis_type      pw            # 基组类型：平面波 (推荐用于高精度声子计算)
ntype           1             # 元素种类数量

# 2. 精度控制 (至关重要)
ecutwfc         100           # [需测试] 平面波截断能量 (Ry)。建议比常规计算高 20-30%
scf_thr         1e-9          # 自洽迭代收敛阈值。Relax 阶段 1e-9 足够，但后续需更高
force_thr_ev    1e-5          # [关键] 力收敛阈值 (eV/Ang)。必须足够小！
stress_thr      1e-2          # 应力收敛阈值 (kbar)。若需优化晶胞参数，此项需严格

# 3. 迭代控制
relax_nmax      100           # 最大离子步数
cal_stress      1             # 计算应力，用于优化晶胞常数 (vc-relax)

# 4. 路径设置
pseudo_dir      ./            # 赝势目录
```

### 1.2.2 核心参数详解与避坑指南

#### 1. `scf_thr` 的精度陷阱
在结构优化阶段，`1e-9` 的电子步收敛标准通常是可以接受的。但是，**请务必在心中拉响警报**：
> **警告**：在后续进行**三阶力常数（3rd Order）**计算时，由于我们要捕捉极微小的非谐效应，`scf_thr` **必须** 设置为 **`1e-12`** (对于 PW 基组)。
>
> 许多初学者沿用 `1e-6` 或 `1e-8` 进行三阶力计算，结果得到的受力数据全是数值噪音，导致热导率计算彻底失败。

#### 2. K 点与超胞的收敛性 (The Convergence Risk)
在官方教学案例中，为了演示速度，通常使用 `2x2x2` 的 K 点网格和较小的超胞。
*   **教学演示**：Si (2x2x2 K-points) $\rightarrow$ 计算值 ~100 W/mK
*   **实验真值**：Si (300K) $\approx$ 150 W/mK
*   **科研标准**：实际研究中，您必须对 **K 点密度** (`K_POINTS`) 和 **截断能量** (`ecutwfc`) 进行严格的收敛性测试。对于热导率，K 点采样不足会直接导致声子寿命计算不准，从而低估热导率。

#### 3. 赝势选择
建议使用模守恒赝势（Norm-Conserving, ONCVPSP），它在声子计算中通常比超软赝势表现更稳定，且与 ABACUS 的 PW 基组配合极佳。

### 1.2.3 结构文件准备 (STRU)

确保您的 `STRU` 文件中晶格常数和原子位置是初始猜测值。如果需要同时优化晶胞常数（Variable Cell Relax），请在 `INPUT` 中设置 `calculation cell-relax` (或旧版本的 `vc-relax`，请查阅您版本的文档)，并确保 `cal_stress 1` 已开启。

**文件名**: `STRU` (示例)
```text
ATOMIC_SPECIES
Si 28.0855 Si_ONCV_PBE-1.0.upf

LATTICE_CONSTANT
1.88972612546  # Bohr to Angstrom conversion factor usually handled internally or set to 1.0/lat_const

LATTICE_VECTORS
0.0 2.7 2.7    # 初始猜测
2.7 0.0 2.7
2.7 2.7 0.0

ATOMIC_POSITIONS
Direct
Si
0.0
2
0.00 0.00 0.00 1 1 1
0.25 0.25 0.25 1 1 1
```

### 1.2.4 运行与检查
提交任务运行 ABACUS。计算结束后，检查输出文件（如 `OUT.Si_relax/running_relax.log`），确认：
1.  最后一步的受力是否小于 `force_thr_ev`。
2.  最后一步的应力是否小于 `stress_thr`。

获得优化后的 `STRU` 文件（通常位于 `OUT.Si_relax/STRU_ION_D` 或类似目录，取决于版本），将其重命名为 `STRU`，作为后续所有微扰计算的**唯一基准结构**。

---

**本章小结**：
我们已经搭建好了 ABACUS + ShengBTE 的工具链，并获得了一个受力极小的基态结构。下一章，我们将进入**二阶谐波性质的计算**，利用 Phonopy 结合 ABACUS 提取声子谱。