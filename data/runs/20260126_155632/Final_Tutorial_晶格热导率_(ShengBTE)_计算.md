# ABACUS 实战进阶：基于第一性原理的晶格热导率计算教程

## 前言

欢迎来到《ABACUS 实战进阶：基于第一性原理的晶格热导率计算教程》。在材料科学迈向“理性设计”的今天，理解并精确预测材料的热输运性质已成为半导体热管理、热电材料研发等领域的核心课题。

**为什么选择 ABACUS？**
ABACUS 作为国产第一性原理计算软件的佼佼者，在处理此类高阶微扰计算时展现了独特的优势：其支持的数值原子轨道（LCAO）基组在保证精度的前提下，显著降低了大规模超胞计算的内存消耗与计算时间；而其平面波（PW）基组则为追求极致精度的基态结构优化提供了坚实保障。这种灵活的基组选择，使得 ABACUS 成为研究声子输运过程的理想工具。

**学习路线图 (Roadmap)**
本教程采用循序渐进的模块化设计，共分为四个核心章节：
1. **理论与准备**：构建从力常数到玻尔兹曼输运方程（BTE）的理论框架，强调基态优化对结果敏感度的重要性。
2. **谐性性质 (Phonopy)**：掌握二阶力常数的提取，这是构建声子能带结构与群速度的地基。
3. **非谐性性质 (thirdorder)**：引入声子-声子相互作用，处理计算量最庞大且最关键的三阶力常数环节。
4. **热导率求解 (ShengBTE)**：完成数据的“总装”，利用 BTE 求解器获得最终的热导率数值并进行科学分析。

**知识图谱与定位**
本教程处于 ABACUS/DFT 知识体系的“高级物性计算”层级。它不仅要求用户掌握基础的电子结构计算，更要求用户跨越平衡态统计力学的门槛，进入激发态与非平衡态输运的研究领域。它是连接微观原子间作用力与宏观热物理性质的桥梁。

**前置知识 (Prerequisites)**
在开始本教程之前，我们建议您已具备以下基础：
- 熟练使用 Linux 命令行环境。
- 掌握第一性原理计算的基本概念（如 K 点收敛性、赝势、交换相关泛函）。
- 具备基本的 Python 脚本操作能力，用于处理 Phonopy 等后处理工具的数据。

---

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

# 第二章：谐性性质与二阶力常数 (Phonopy)

在计算晶格热导率的宏大工程中，二阶力常数（2nd Order Force Constants）是地基。它决定了声子谱（Phonon Spectrum）的形状、声子的群速度以及比热容。如果地基打不牢——例如力常数中包含数值噪音——后续的三阶散射计算将毫无意义。

本章我们将利用 **Phonopy** 结合 **ABACUS**，完成从微扰结构生成、受力计算到力常数矩阵提取的全过程。

> **⚠️ 流程全景提示**
> 计算晶格热导率涉及多个软件的接力，请务必保持清醒的头脑：
> 1.  **Phase 1 (本章)**: ABACUS + Phonopy $\rightarrow$ **二阶力常数** (需单位转换)。
> 2.  **Phase 2 (下一章)**: ABACUS + thirdorder $\rightarrow$ **三阶力常数**。
> 3.  **Phase 3**: ShengBTE $\rightarrow$ **热导率**。

---

## 2.1 超胞建立与微扰生成 (Phonopy)

为了计算力常数，我们需要采用有限位移法（Finite Displacement Method）。即在一个扩大的超胞（Supercell）中，将原子稍微移动一点点，计算它受到的回复力。

### 2.1.1 准备工作
请确保你已安装 Phonopy，并准备好优化过的晶胞结构文件 `STRU`。
在工作目录下创建一个名为 `2nd` 的文件夹，并进入该目录。

### 2.1.2 配置文件 setting.conf
编写 Phonopy 的配置文件 `setting.conf`。这里最关键的参数是 `DIM`，它定义了超胞的大小。

```ini
# setting.conf
DIM = 2 2 2
ATOM_NAME = Si
```

*   **DIM**: 指定在 $x, y, z$ 三个方向上扩胞的倍数。
    *   *教学演示*: `2 2 2` (为了快速跑通流程)。
    *   *科研实战*: **必须进行收敛性测试**。对于热导率计算，声子平均自由程较长，通常需要较大的超胞（如 4x4x4 或更大）来消除周期性边界条件引入的非物理相互作用。

### 2.1.3 生成微扰结构
使用 Phonopy 的 ABACUS 接口生成带有微扰的结构文件：

```bash
phonopy setting.conf --abacus -d
```

**执行结果**:
Phonopy 会根据晶体对称性，生成最少数量的必要微扰结构。
*   生成文件名为 `STRU-001`, `STRU-002` ...
*   生成 `phonopy_disp.yaml`：记录了具体的位移信息（后续处理必须保留此文件）。

---

## 2.2 ABACUS 力常数计算 (SCF)

现在我们需要计算这些微扰结构中原子的受力。这是一个标准的 DFT 自洽计算（SCF）过程。

### 2.2.1 输入文件 INPUT 设置
编写 `INPUT` 文件。这里有两个至关重要的细节，直接决定了你计算出的热导率是物理结果还是“随机数”。

```bash
INPUT_PARAMETERS
# ... 基础参数 (ecutwfc, basis_type 等) ...

calculation     scf         # 进行自洽计算
stru_file       STRU-001    # 指定读取的结构文件

# --- 精度控制 (生死攸关) ---
scf_thr         1e-8        # LCAO 基组推荐值
# scf_thr       1e-12       # PW (平面波) 基组推荐值
force_thr_ev    1e-7        # 力的收敛标准
cal_force       1           # 显式开启受力计算输出
```

### 🛑 专家锦囊：关于精度的血泪教训
在常规的能带计算中，`scf_thr` 设置为 `1e-6` 通常就够了。但在**声子和热导率计算**中，这远远不够！
*   **微扰极其微小**: 原子位移通常只有 0.01 Å 量级，产生的回复力非常微弱。
*   **噪音放大**: 如果 SCF 收敛精度不高，电荷密度的数值噪音会导致受力计算出现误差。这些误差在后续计算力常数矩阵时会被放大，导致声子谱出现虚频（软模），或者热导率数值完全错误。
*   **推荐标准**:
    *   **LCAO (数值原子轨道)**: 至少 `1e-8`。
    *   **PW (平面波)**: 建议 `1e-12`。

### 2.2.2 批量计算
你需要对 Phonopy 生成的每一个 `STRU-XXX` 文件运行一次 ABACUS。
*   **技巧**: 利用 `stru_file` 参数。你不需要把 `STRU-001` 重命名为 `STRU`，只需在 `INPUT` 中修改 `stru_file = STRU-001` 即可。
*   **脚本建议**: 编写一个简单的 Shell 脚本循环运行，或者使用作业调度系统提交数组作业。

---

## 2.3 力常数提取与单位转换 (核心坑点)

当所有 SCF 计算完成后，我们收集受力并计算力常数。这里存在一个**极易踩中的单位陷阱**。

### 2.3.1 收集受力 (FORCE_SETS)
假设你的计算结果保存在对应的文件夹或日志中（例如 `OUT.STRU-001/running_scf.log`），使用以下命令收集受力：

```bash
# 假设只有一个微扰结构，日志文件路径需根据实际情况调整
phonopy -f OUT.STRU-001/running_scf.log
```
这将生成 `FORCE_SETS` 文件。

### 2.3.2 计算二阶力常数矩阵
编写 `band.conf` 用于生成力常数。

```ini
# band.conf
ATOM_NAME = Si
DIM = 2 2 2
# ... 其他声子谱相关参数 ...

# --- 关键参数 ---
FORCE_CONSTANTS = WRITE        # 输出力常数文件
FULL_FORCE_CONSTANTS = .TRUE.  # 【必须开启】输出完整的力常数矩阵
```

运行 Phonopy：
```bash
phonopy -p band.conf --abacus
```
这将生成 `FORCE_CONSTANTS` 文件。

### 2.3.3 单位转换 (ShengBTE 兼容性处理)

**这是本章最大的坑点，请仔细阅读。**

*   **ABACUS/Phonopy 输出**: Phonopy 读取 ABACUS 结果时，默认输出的力常数单位通常包含原子单位制（Bohr/au），即 `eV/(Å·au)` (其中 1 au $\approx$ 0.529 Å)。
*   **ShengBTE 输入要求**: ShengBTE 严格要求二阶力常数单位为 `eV/Å²`。

如果不进行转换，直接将 `FORCE_CONSTANTS` 喂给 ShengBTE，计算出的热导率将是错误的量级。

**解决方案**:
我们需要使用一个辅助脚本 `au2si.py` 来完成这个转换。
> **注**: `au2si.py` 并非 ABACUS 内置命令，它通常随本教程的案例文件提供（位于 `examples/interface_ShengBTE/LCAO/2nd/` 目录下）。

确保 `au2si.py` 在当前目录下，运行：

```bash
python au2si.py
```

该脚本会读取当前的 `FORCE_CONSTANTS`，进行单位换算，并生成一个新的文件（通常命名为 `FORCE_CONSTANTS_2ND` 或覆盖原文件，具体视脚本实现而定）。**请务必将转换后的文件重命名为 `FORCE_CONSTANTS_2ND`**，这是 ShengBTE 识别的标准文件名。

---

## 2.4 本章总结与风险提示

至此，我们已经获得了符合 ShengBTE 要求的二阶力常数文件 `FORCE_CONSTANTS_2ND`。

**再次强调风险**:
在本教程的案例中，为了演示速度，我们使用了 `2x2x2` 的超胞和较稀疏的 K 点。
*   **案例结果**: Si 的热导率可能算出 ~100 W/mK。
*   **实验真值**: Si 的热导率约为 150 W/mK。
*   **原因**: 超胞尺寸不足以截断长程力常数，且 K 点采样未收敛。

**科研建议**: 在正式计算中，请务必测试超胞大小（如对比 3x3x3, 4x4x4）对声子谱和热导率的影响，直至结果收敛。

# 第三章：非谐性性质与三阶力常数 (thirdorder)

在上一章中，我们通过 Phonopy 获得了材料的谐性性质（声子谱）。然而，谐波近似下的声子是具有无限寿命的准粒子，这意味着材料的热导率将是无穷大。为了计算真实的晶格热导率，我们必须引入**非谐性（Anharmonicity）**，即声子之间的相互作用。

三阶力常数（3rd Order Force Constants）是描述三声子散射过程的核心物理量。本章将指导你如何结合 ABACUS 与 `thirdorder` 程序（ShengBTE 套件的一部分）完成这一计算中最昂贵且最易出错的步骤。

## 3.0 全流程概览与工具准备

计算晶格热导率是一个多软件耦合的复杂流程，为了不迷失方向，请务必牢记以下三个阶段：

1.  **2nd Order (谐性)**: ABACUS + Phonopy $\rightarrow$ `FORCE_CONSTANTS_2ND`
    *   *注意*: 需使用案例脚本 `au2si.py` 将单位转换为 ShengBTE 要求的 eV/Å²。
2.  **3rd Order (非谐性)**: ABACUS + thirdorder $\rightarrow$ `FORCE_CONSTANTS_3RD`
    *   *本章重点*: 涉及结构微扰、高精度受力计算及格式转换。
3.  **Transport (输运)**: ShengBTE
    *   *汇总*: 读取上述两个文件及 `CONTROL` 文件计算热导率 $\kappa$。

> **⚠️ 脚本依赖警告**
> 本章中提到的 `pos2stru.py` 和 `aba2vasp.py` 并非 ABACUS 内置命令，而是官方案例库（`examples/interface_ShengBTE`）中提供的 Python 辅助脚本。请务必确保这些脚本已下载并放置在你的工作目录中，且你的环境中已安装 `ASE` (Atomic Simulation Environment) 库。

---

## 3.1 三阶微扰结构生成

`thirdorder` 程序原生支持 VASP 和 QE，但暂未直接支持 ABACUS。因此，我们需要采用“曲线救国”的策略：将 ABACUS 的结构伪装成 POSCAR，生成微扰后再转回 STRU。

### 步骤 1: 准备超胞 POSCAR
首先，将优化好的晶胞（Unit Cell）扩胞为超胞（Supercell）。虽然你可以手动操作，但推荐先将优化后的 `STRU` 转换为 VASP 格式的 `POSCAR`。

### 步骤 2: 生成微扰构型 (Sow)
使用 `thirdorder` 提供的脚本生成一系列微扰结构。假设我们使用 $2\times2\times2$ 的超胞，且考虑第 nearest-neighbor 相互作用（参数 `-2`）：

```bash
# 语法: thirdorder_vasp.py sow [Nx] [Ny] [Nz] [Cutoff]
thirdorder_vasp.py sow 2 2 2 -2
```

运行后，目录下会生成大量的 `3RD.POSCAR.*` 文件（例如 `3RD.POSCAR.01`, `3RD.POSCAR.02` ...）。这些文件包含了原子微小的位移。

### 步骤 3: 格式回转 (POSCAR $\rightarrow$ STRU)
这是最关键的一步。我们需要将这些 POSCAR 转换回 ABACUS 的 STRU 格式。

> **⛔ 严禁使用 dpdata 进行转换**
> 请绝对**不要**使用 `dpdata` 将这些 POSCAR 转为 STRU。
> **原因**: `dpdata` 在处理晶胞时，会强制将晶格矢量旋转为下三角矩阵形式。然而，`thirdorder` 生成的微扰是基于原始晶格方向定义的。如果晶格发生了旋转，原子的位移方向与晶格的相对关系就会被破坏，导致后续计算出的力常数完全错误。

**正确做法**: 使用基于 ASE 的 `pos2stru.py` 脚本，它能保持晶格取向不变。

```bash
# 批量转换脚本示例
for file in 3RD.POSCAR.*; do
    # 提取编号，如 01, 02
    id=${file##*.}
    # 创建对应的计算文件夹
    mkdir SCF-$id
    # 将 POSCAR 转为 STRU 并放入文件夹
    # 假设 pos2stru.py 接受输入文件名和输出文件名作为参数
    python pos2stru.py $file > SCF-$id/STRU
done
```

---

## 3.2 极高精度受力计算 (High-Precision SCF)

现在你需要对这几十甚至上百个微扰结构进行 SCF 计算以获取受力。**这是整个热导率计算中对数值精度要求最高的环节。**

三阶力常数本质上是能量对原子位移的三阶导数（或力对位移的二阶导数）。由于微扰位移量很小，由此产生的受力变化也非常微小。如果 SCF 收敛精度不够，数值噪音将淹没真实的物理信号，导致最终的热导率结果出现巨大误差（甚至出现负值）。

### INPUT 参数设置核心原则

在 `INPUT` 文件中，除了常规的 `calculation = scf` 和 `cal_force = 1` 外，必须严格设置 `scf_thr`。

#### 1. LCAO 基组 (数值原子轨道)
对于 LCAO 基组，推荐精度如下：

```javascript
INPUT_PARAMETERS
{
    calculation     scf
    basis_type      lcao
    cal_force       1
    
    // 核心精度控制
    scf_thr         1e-8    // 必须达到 1e-8 或更小
    force_thr_ev    1e-4    // 辅助判断，但 scf_thr 是关键
}
```

#### 2. PW 基组 (平面波)
平面波基组对噪音更敏感，且通常能达到更高的收敛精度。**强烈建议**设置为 `1e-12`：

```javascript
INPUT_PARAMETERS
{
    calculation     scf
    basis_type      pw
    cal_force       1
    ecutwfc         100     // 确保能量截断足够高
    
    // 核心精度控制
    scf_thr         1e-12   // 粗体警告：1e-8 对于 PW 计算三阶力常数通常是不够的！
}
```

> **专家经验**: 
> 很多初学者发现计算出的热导率与实验值相差甚远，或者声子散射率呈现无规律的杂乱分布，90% 的原因都是 `scf_thr` 设置过大。不要为了节省一点计算时间而牺牲精度，这一步的“脏数据”会导致后续所有工作作废。

---

## 3.3 格式伪装与力常数生成 (Reap)

当所有 `SCF-*` 文件夹中的计算都正常结束后，我们需要收集受力数据。由于 `thirdorder` 只能读取 VASP 的 `vasprun.xml`，我们需要再次进行“伪装”。

### 步骤 1: 生成 vasprun.xml
使用辅助脚本 `aba2vasp.py` 进入每个文件夹，读取 ABACUS 的输出文件（如 `running_scf.log` 或二进制输出），并生成符合 VASP XML 规范的 `vasprun.xml` 文件。

```bash
# 批量处理示例
for dir in SCF-*; do
    cd $dir
    # 运行转换脚本，生成 vasprun.xml
    python ../aba2vasp.py 
    cd ..
done
```
*注：生成的 `vasprun.xml` 不需要包含完整的电子结构信息，只要包含正确的原子坐标和受力 (`<varray name="forces">`) 即可被 `thirdorder` 识别。*

### 步骤 2: 提取三阶力常数
最后，使用 `thirdorder_vasp.py` 的 `reap` 模式收集所有数据。确保命令行参数与 `sow` 阶段完全一致（超胞大小和 cutoff）。

```bash
# 语法: thirdorder_vasp.py reap [Nx] [Ny] [Nz] [Cutoff]
# 利用 find 命令将所有 vasprun.xml 排序后喂给程序
find SCF-* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -2
```

如果一切顺利，当前目录下将生成一个名为 `FORCE_CONSTANTS_3RD` 的文件。

---

## 3.4 最终汇总与风险提示

至此，你已经准备好了 ShengBTE 所需的所有核心文件：
1.  `FORCE_CONSTANTS_2ND` (来自 Phonopy + `au2si.py`)
2.  `FORCE_CONSTANTS_3RD` (来自 thirdorder)
3.  `CONTROL` (ShengBTE 的输入参数文件)

将它们放入同一个文件夹（通常命名为 `shengbte`），运行 ShengBTE 程序即可得到热导率。

### ⚠️ 风险提示：收敛性测试
在官方提供的教学案例中，为了演示速度，通常使用 $2\times2\times2$ 的 K 点网格和较小的超胞。这会导致计算出的硅（Si）热导率约为 100 W/(m·K)，而实验值约为 150 W/(m·K)。

**在实际科研中，你必须进行以下收敛性测试：**
1.  **K 点收敛**: 必须测试更密的 K 点网格（如 $4\times4\times4$ 或更高）。
2.  **超胞尺寸**: 三阶力常数的截断半径对热导率影响很大，需测试更大超胞（如 $3\times3\times3$ 或 $4\times4\times4$）以包含更长程的相互作用。

切记：**Tutorial 只是教你流程，Research 必须追求收敛。**

# 第四章：ShengBTE 热导率计算与分析

万事俱备，只欠东风。在前面的章节中，我们已经利用 ABACUS 结合 Phonopy 计算了二阶力常数，并结合 thirdorder 程序计算了三阶力常数。本章将作为“总装车间”，指导你如何配置 ShengBTE 的核心输入文件 `CONTROL`，运行玻尔兹曼输运方程（BTE）求解器，并对最终的热导率结果进行科学分析。

## 4.0 全流程回顾与数据准备

在正式运行 ShengBTE 之前，让我们通过一个逻辑图解回顾整个工作流，确保你手头的文件是正确且完整的。这是一个多软件耦合的复杂流程，清晰的阶段划分至关重要。

```mermaid
graph TD
    subgraph Stage 1: 2nd Order [ABACUS + Phonopy]
        A[结构优化 (Relax)] --> B[Phonopy 扩胞]
        B --> C[ABACUS SCF 计算受力]
        C --> D[Phonopy 生成 FORCE_CONSTANTS]
        D -->|关键: 单位转换| E(FORCE_CONSTANTS_2ND)
        style E fill:#f9f,stroke:#333,stroke-width:2px
    end

    subgraph Stage 2: 3rd Order [ABACUS + thirdorder]
        F[thirdorder 生成微扰构型] --> G[ABACUS SCF 计算受力]
        G -->|关键: 格式转换| H[aba2vasp.py -> vasprun.xml]
        H --> I[thirdorder 收集数据]
        I --> J(FORCE_CONSTANTS_3RD)
        style J fill:#f9f,stroke:#333,stroke-width:2px
    end

    subgraph Stage 3: Solver [ShengBTE]
        E --> K[ShengBTE]
        J --> K
        L[CONTROL 文件] --> K
        K --> M[BTE.kappa 热导率结果]
    end
```

**检查清单：**
进入你的 `shengbte` 工作目录，确保包含以下三个核心文件：
1.  **`FORCE_CONSTANTS_2ND`**: 由 Phonopy 产生并经过 `au2si.py` 转换单位后的二阶力常数。
2.  **`FORCE_CONSTANTS_3RD`**: 由 thirdorder 产生的原始三阶力常数。
3.  **`CONTROL`**: 我们即将编写的 ShengBTE 配置文件。

---

## 4.1 CONTROL 文件配置

`CONTROL` 是 ShengBTE 的主输入文件，采用 Fortran 的 Namelist 格式。它定义了晶体结构、超胞大小、Q 点采样网格以及求解参数。

以下是针对硅（Si）案例的标准 `CONTROL` 文件示例及详细解析：

```fortran
&allocations
    nelements = 1,          ! 元素种类的数量
    natoms = 2,             ! 单胞（Unit Cell）中的原子总数
    ngrid(:) = 10 10 10,    ! Q 点网格采样 (Q-grid)，用于积分声子玻尔兹曼方程
    norientations = 0,      ! 纳米线/薄膜方向设置，体材料设为 0
/

&crystal
    lfactor = 0.1,          ! 长度单位转换因子。ShengBTE 内部使用 nm，若输入为 Angstrom，此处设为 0.1
    lattvec(:,1) = 0.0 2.81594778072 2.81594778072,  ! 晶格矢量 a1 (Angstrom)
    lattvec(:,2) = 2.81594778072 0.0 2.81594778072,  ! 晶格矢量 a2
    lattvec(:,3) = 2.81594778072 2.81594778072 0.0,  ! 晶格矢量 a3
    elements = "Si",        ! 元素名称
    types = 1 1,            ! 每个原子对应的元素类型索引（对应 elements 列表）
    positions(:,1) = 0.00 0.00 0.00,  ! 原子 1 的分数坐标
    positions(:,2) = 0.25 0.25 0.25,  ! 原子 2 的分数坐标
    scell(:) = 2 2 2,       ! 计算力常数时使用的超胞大小 (Supercell size)
/

&parameters
    T = 300,                ! 计算热导率的温度 (Kelvin)
    scalebroad = 1.0,       ! 高斯展宽的缩放因子，通常设为 1.0
/

&flags
    nonanalytic = .FALSE.,  ! 是否包含非解析项校正（极性材料需设为 .TRUE. 并提供 BORN 文件）
    nanowires = .FALSE.,    ! 是否计算纳米线
/
```

### 关键参数详解

1.  **`scell` (Supercell Size)**:
    *   **含义**: 这里填写的必须是你**计算二阶和三阶力常数时所用的超胞大小**。
    *   **注意**: 在本教程的 Si 案例中，我们为了演示速度使用了 2x2x2 的超胞，因此这里填 `2 2 2`。如果二阶和三阶使用了不同大小的超胞，ShengBTE 会尝试进行插值，但强烈建议两者保持一致以保证精度。

2.  **`ngrid` (Q-grid)**:
    *   **含义**: 这是求解 BTE 时在布里渊区采样的网格密度。
    *   **策略**: 该参数直接决定结果的收敛性。`10 10 10` 仅为演示用。实际科研中，你需要测试 `15 15 15`, `20 20 20` 等，直到热导率数值不再剧烈变化。

3.  **`lattvec` & `positions`**:
    *   **来源**: 这些数据应直接取自你进行结构优化（Relax）后的 `STRU` 文件或转换后的 `POSCAR`。务必保证精度，不要随意截断小数位。

4.  **`lfactor`**:
    *   ABACUS 和大多数 DFT 软件使用 Å (Angstrom) 为单位，而 ShengBTE 内部使用 nm。因此，当 `lattvec` 单位为 Å 时，必须设置 `lfactor = 0.1`。

---

## 4.2 运行与结果分析

### 4.2.1 提交计算任务

ShengBTE 是一个支持 MPI 并行的程序。在准备好 `CONTROL`, `FORCE_CONSTANTS_2ND`, `FORCE_CONSTANTS_3RD` 后，使用以下命令运行：

```bash
# 假设使用 10 个核心并行计算
mpirun -n 10 ShengBTE
```

程序运行过程中会输出详细的日志。如果配置正确，计算通常在几分钟内完成（取决于 `ngrid` 的密度）。

### 4.2.2 结果分析

运行结束后，目录下会生成多个输出文件。最核心的文件是 **`BTE.kappa`**。

**`BTE.kappa` 文件结构示例：**
```text
# T          kappa_xx      kappa_yy      kappa_zz      ... (张量分量)
  300.0000   102.4512      102.4512      102.4512      ...
```

**数据解读与对比：**

*   **计算值**: 在本案例（LCAO 基组，2x2x2 超胞，2x2x2 K点）中，300K 下的计算结果大约在 **100 W/(m·K)** 左右。
*   **实验值**: 硅单晶在 300K 的实验热导率约为 **140 - 150 W/(m·K)**。

**为什么计算值偏小？（风险提示）**
你可能会发现计算结果（~100）明显低于实验值（~150）。这**不是**软件的错误，而是由于我们在教学案例中为了节省计算时间，采用了极低的参数设置：
1.  **超胞过小**: 2x2x2 的超胞截断了长程力常数，无法捕捉长波声子的贡献。
2.  **K 点稀疏**: 电子结构计算时的 2x2x2 K 点采样对于半导体 Si 来说远远不够。
3.  **Q 点稀疏**: `CONTROL` 中的 `ngrid` 仅为 10x10x10。

**科研级建议**: 在实际研究中，你必须进行严格的收敛性测试。通常对于 Si 这样的体系，需要 4x4x4 或 5x5x5 的超胞，以及更密的 K 点网格，才能得到与实验吻合的结果。

---

## 附录：常见问题与进阶建议

### 1. 常见报错与排查

*   **报错**: `Error reading FORCE_CONSTANTS`
    *   **原因 1**: Phonopy 计算二阶力常数时，未在 `band.conf` 或 `setting.conf` 中设置 `FULL_FORCE_CONSTANTS = .TRUE.`。ShengBTE 需要完整的力常数矩阵，而不仅仅是简化的。
    *   **原因 2**: 单位未转换。ABACUS+Phonopy 输出的单位是 eV/(Å·au)，ShengBTE 需要 eV/Å²。
    *   **解决**: 确保运行了 `python au2si.py` 生成 `FORCE_CONSTANTS_2ND`。

*   **现象**: 热导率结果包含大量噪音或数值异常（如负值、极大值）。
    *   **原因**: 三阶力常数计算精度不足。这是新手最容易踩的坑。
    *   **解决**: 检查计算三阶力常数时的 SCF 收敛精度。
        *   **LCAO 基组**: `scf_thr` 至少设为 **1e-8**。
        *   **PW (平面波) 基组**: `scf_thr` 必须设为 **1e-12**。
        *   **警告**: 默认的 `1e-6` 对于微扰计算完全不够，会导致数值微分结果全是数值噪音。

### 2. 辅助脚本清单

本教程中提到的脚本并非 ABACUS 内置命令，而是随案例提供的 Python 辅助工具。请确保它们在你的工作目录中：

| 脚本名 | 功能描述 | 依赖库 |
| :--- | :--- | :--- |
| **`au2si.py`** | 将 Phonopy 输出的力常数单位从原子单位制转换为 ShengBTE 所需的 eV/Å²。 | NumPy |
| **`pos2stru.py`** | 将 `thirdorder` 生成的 `POSCAR` 格式微扰结构转换为 ABACUS 的 `STRU` 格式。 | ASE |
| **`aba2vasp.py`** | 将 ABACUS 的 SCF 输出（力、能量、应力）封装成 `vasprun.xml` 格式，以便 `thirdorder` 读取。 | 无 |

### 3. 进阶：如何测试收敛性

在发表论文前，请务必完成以下测试：
1.  **K 点收敛**: 固定超胞大小，增加 DFT 计算时的 K 点密度（如 2x2x2 -> 4x4x4 -> 6x6x6），观察二阶声子谱是否稳定。
2.  **超胞收敛**: 这是最耗时的。计算 3x3x3, 4x4x4, 5x5x5 超胞的三阶力常数。通常三阶力常数的截断半径影响最大。
3.  **Q 网格收敛**: 在 `CONTROL` 文件中增加 `ngrid`（如 20x20x20 -> 30x30x30），直到 `BTE.kappa` 收敛。这一步计算成本极低，应首先保证收敛。

---

## 附录：进阶学习指南

恭喜你完成了本教程的学习！计算晶格热导率是声子物理研究的起点而非终点。为了进一步提升研究深度，我们建议在以下方向进行探索：

### 1. 进阶学习主题 (Extended Reading)
- **温度相关有效势 (TDEP)**：本教程基于冻结声子法，适用于低温或中温。对于强非谐性材料或高温工况，建议学习如何利用分子动力学轨迹结合 TDEP 方法修正力常数。
- **电子-声子耦合**：在金属或高掺杂半导体中，电子对声子的散射不可忽略。后续可探索如何利用相关工具计算电声耦合对热导率的贡献。
- **纳米尺度效应**：通过 ShengBTE 的边界散射模块，可以研究薄膜、纳米线的尺寸效应对热导率的影响，这在微纳电子器件散热研究中至关重要。

### 2. 通用调试建议 (Troubleshooting)
- **收敛性是生命线**：如果发现声子谱存在虚频（虚部频率），请优先检查基态结构优化（Relax）的受力标准（Force Convergence）是否足够严格（建议优于 1e-4 eV/Å）。
- **超胞尺寸效应**：三阶力常数的截断距离对结果影响极大。如果热导率随超胞增大而显著波动，请务必进行截断半径的收敛性测试。
- **数值噪音过滤**：在处理低对称性体系时，确保利用对称性约束来减少力常数矩阵中的数值噪音。

### 3. 官方资源与文档
- **ABACUS 官方文档**：[https://abacus.deepmodeling.com/](https://abacus.deepmodeling.com/) —— 获取最新的参数说明与功能更新。
- **开源算例库**：访问 ABACUS 在 Gitee 或 GitHub 的 `abacus-user-guide` 仓库，下载更多关于 LCAO 和 PW 基组的实战案例。
- **社区支持**：积极参与 DeepModeling 社区讨论，与全球开发者共同解决计算中的疑难杂症。

**探索科学的道路永无止境，愿本教程能成为你科研征途中的有力阶梯！**
