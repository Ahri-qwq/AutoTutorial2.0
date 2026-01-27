<!-- META_START -->
# ABACUS 实战教程：晶格热导率 (ShengBTE) 计算

## 前言
- **本书逻辑**: 本教程采用“分而治之”的教学策略，将复杂的晶格热导率计算拆解为三个核心模块：谐性声子性质（二阶力常数）、非谐性散射（三阶力常数）以及玻尔兹曼输运方程求解（ShengBTE）。教程特别强调多软件接口之间的数据流转、单位换算以及高精度计算的必要性。
- **核心知识点**: 晶格动力学、二阶/三阶力常数计算、Phonopy 接口、thirdorder 接口、ShengBTE 输入配置、单位制转换、高精度 SCF 收敛技巧。

<!-- CHAPTER_START -->
## 第一章：计算准备与理论概览
**本章逻辑**: 在开始繁杂的计算流程前，必须明确软件生态的依赖关系，并准备好高质量的基态结构。这是所有后续微扰计算的基石。

### Section 1.1: 工具链环境配置
- **内容**: 介绍 ABACUS 与外部软件（Phonopy, ShengBTE, thirdorder, ASE）的协作关系。重点说明辅助脚本（如 `pos2stru.py`, `aba2vasp.py`, `au2si.py`）的作用与获取方式。
- **关键参数**: 无（侧重环境配置）。

### Section 1.2: 高精度结构优化 (Relax)
- **内容**: 对初始晶胞进行严格的结构弛豫。强调晶格热导率对结构微扰极其敏感，因此必须保证基态结构的受力极小。
- **关键参数**: 
    - `calculation`: `relax` (结构优化)
    - `ecut`: 需收敛测试 (如 100 Ry)
    - `kspacing` 或 `k_points`: 需高密度采样 (如 2x2x2 仅为示例，实际需更高)

<!-- CHAPTER_START -->
## 第二章：谐性性质与二阶力常数 (Phonopy)
**本章逻辑**: 二阶力常数决定了声子谱和群速度。本章通过 Phonopy 接口生成微扰结构，利用 ABACUS 计算受力，最终获得符合 ShengBTE 格式要求的二阶力常数矩阵。

### Section 2.1: 超胞建立与微扰生成
- **内容**: 使用 Phonopy 生成计算所需的超胞（Supercell）及位移微扰结构。介绍如何配置 `setting.conf`。
- **关键参数**: 
    - Phonopy 参数: `DIM` (超胞尺寸)

### Section 2.2: ABACUS 力常数计算
- **内容**: 对生成的微扰结构（如 `STRU-001`）进行 SCF 自洽计算以获取原子受力。介绍如何通过 `stru_file` 参数批量读取结构。
- **关键参数**: 
    - `calculation`: `scf`
    - `stru_file`: `STRU-001` (指定结构文件)

### Section 2.3: 力常数提取与单位转换 (核心坑点)
- **内容**: 使用 Phonopy 收集受力并计算二阶力常数。**重点讲解** ABACUS/Phonopy 输出单位与 ShengBTE 输入单位的不兼容性，以及如何使用 `au2si.py` 进行单位转换。
- **关键参数**: 
    - Phonopy 参数: `FULL_FORCE_CONSTANTS = .TRUE.` (**必须开启**)
    - 脚本操作: `python au2si.py` (转换单位 eV/Å·au -> eV/Å²)

<!-- CHAPTER_START -->
## 第三章：非谐性性质与三阶力常数 (thirdorder)
**本章逻辑**: 三阶力常数描述了声子-声子散射，是热导率计算中最昂贵且最易出错的步骤。本章重点解决 ABACUS 与 `thirdorder` 程序（主要支持 VASP/QE）的接口兼容问题及精度控制。

### Section 3.1: 三阶微扰结构生成
- **内容**: 利用 `thirdorder_vasp.py` 生成大量微扰结构。**特别警告**：禁止使用 `dpdata` 进行格式转换（防止晶格旋转导致的受力错误），推荐使用基于 ASE 的 `pos2stru.py` 将 POSCAR 转为 STRU。
- **关键参数**: 无（侧重脚本操作流程）。

### Section 3.2: 极高精度受力计算 (High-Precision SCF)
- **内容**: 对数十/上百个微扰结构进行受力计算。**核心强调**：三阶力常数对数值噪音极度敏感，必须设置极严苛的收敛精度，区分 LCAO 与 PW 基组的要求。
- **关键参数**: 
    - `scf_thr`: `1e-8` (LCAO 推荐值) 或 `1e-12` (PW 推荐值)
    - `basis_type`: `lcao` 或 `pw`

### Section 3.3: 格式伪装与力常数生成
- **内容**: 将 ABACUS 的计算结果“伪装”成 VASP 格式供 `thirdorder` 读取。使用 `aba2vasp.py` 生成 `vasprun.xml`，最后运行 `thirdorder_vasp.py` 提取三阶力常数。
- **关键参数**: 无（侧重数据流处理）。

<!-- CHAPTER_START -->
## 第四章：ShengBTE 热导率计算与分析
**本章逻辑**: 万事俱备，汇总二阶和三阶力常数，配置控制文件，运行 ShengBTE 求解玻尔兹曼输运方程，得到最终热导率。

### Section 4.1: CONTROL 文件配置
- **内容**: 详解 ShengBTE 的输入文件 `CONTROL`，包括晶格参数、超胞设置、Q点网格采样等。
- **关键参数**: 
    - ShengBTE 参数: `nelements`, `natoms`, `ngrid`, `scell`, `lattvec`

### Section 4.2: 运行与结果分析
- **内容**: 提交 ShengBTE 并行计算任务。分析输出文件（如 `BTE.kappa`），对比理论值与实验值。讨论超胞大小和 K 点/Q 点网格对结果收敛性的影响（解释为何示例结果偏小）。
- **关键参数**: 无。

<!-- APPENDIX_START -->
## 附录：常见问题与进阶建议
- **常见报错**: ShengBTE 读取力常数失败（检查 `FULL_FORCE_CONSTANTS` 和单位转换）、热导率数值异常（检查 `scf_thr`）。
- **脚本清单**: `au2si.py`, `pos2stru.py`, `aba2vasp.py` 的功能速查。
- **进阶**: 如何测试超胞收敛性、磁性体系的注意事项。