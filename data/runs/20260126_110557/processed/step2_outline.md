<!-- META_START -->
# ABACUS 实战教程：结构文件 (STRU) 构建

## 前言
- **本书逻辑**: 本教程遵循“物理概念 -> 文件解剖 -> 工具实战 -> 进阶技巧”的教学逻辑。首先确立晶体结构在 ABACUS 中的物理描述方式，随后深入解析 STRU 文件的语法细节，接着介绍三种主流的构建工具（ASE, ATOMKIT, abacustest）以适应不同场景，最后处理磁性设置与结构优化约束等进阶需求。
- **核心知识点**: 
    - STRU 文件结构与核心关键词解析。
    - 玻尔 (Bohr) 与埃 (Angstrom) 的单位换算机制。
    - 赝势 (Pseudopotential) 与数值轨道 (Numerical Orbital) 的关联配置。
    - 基于 Python (ASE) 和命令行工具 (ATOMKIT/abacustest) 的格式转换工作流。
    - 初始磁矩与原子动力学约束的设置。

<!-- CHAPTER_START -->
## 第一章：STRU 文件基础与物理内涵
**本章逻辑**: 在动手写文件之前，必须先理解 ABACUS 如何描述物理世界。本章通过与 VASP (POSCAR) 类比，建立概念映射，并重点阐述 ABACUS 特有的单位制和基组依赖性，为后续操作扫清认知障碍。

### Section 1.1: STRU 文件在计算流程中的角色
- **内容**: 解释 STRU 文件作为几何初始态定义的物理意义（实空间描述）。类比 VASP 的 POSCAR，说明其包含的晶格、原子位置及元素信息。辨析 STRU 与 INPUT（控制参数）、KPT（倒空间采样）的关系。
- **关键参数**: 无

### Section 1.2: 核心语法与单位制详解
- **内容**: 逐行解析 STRU 文件的标准格式。重点讲解 `LATTICE_CONSTANT` 的作用机制，阐明为什么该值通常设为 1.8897（Bohr to Angstrom conversion factor），以及如何通过缩放因子正确定义 `LATTICE_VECTORS`。
- **关键参数**: `ATOMIC_SPECIES`, `LATTICE_CONSTANT`, `LATTICE_VECTORS`, `ATOMIC_POSITIONS`

### Section 1.3: 基组依赖性：LCAO 与 PW 的区别
- **内容**: 阐述在原子轨道 (LCAO) 和平面波 (PW) 两种基组下 STRU 文件的写法差异。重点说明 LCAO 模式下必须包含 `NUMERICAL_ORBITAL` 模块，且需与 `ATOMIC_SPECIES` 顺序严格对应。
- **关键参数**: `NUMERICAL_ORBITAL`

<!-- CHAPTER_START -->
## 第二章：基于 ASE 的 Python 脚本化构建
**本章逻辑**: ASE 是计算材料学最强大的脚本工具，也是 ABACUS 官方推荐的接口。本章从环境配置入手，教授如何利用 Python 的灵活性进行从零建模和格式转换，适合批量处理和复杂结构操作。

### Section 2.1: 环境配置与 ASE-ABACUS 接口安装
- **内容**: 强调官方 ASE (`pip install ase`) 与 ABACUS 分支版本的区别。演示通过 git clone 安装 `ase-abacus` 接口的正确步骤，以及配置 `ABACUS_PP_PATH` 和 `ABACUS_ORBITAL_PATH` 环境变量的重要性。
- **关键参数**: 无 (涉及 Shell 命令)

### Section 2.2: 从零构建：以 FCC 铝为例
- **内容**: 演示如何使用 `ase.build.bulk` 创建简单晶体结构，并利用 `write` 函数导出标准 STRU 文件。展示最基础的脚本编写流程。
- **关键参数**: 无 (涉及 Python API: `ase.io.write`, `format='abacus'`)

### Section 2.3: 格式转换：CIF/POSCAR 转 STRU
- **内容**: 这是最常见的应用场景。详细讲解如何读取 CIF 或 POSCAR 文件，并通过 Python 字典 (`dict`) 显式指定赝势 (`pp`) 和轨道 (`basis`) 文件映射，最终生成完整的 STRU 文件。
- **关键参数**: 无 (涉及 Python API: `pp={...}`, `basis={...}`)

<!-- CHAPTER_START -->
## 第三章：高效工具箱：ATOMKIT 与 abacustest
**本章逻辑**: 对于不想编写 Python 脚本的用户，或者需要高通量准备任务的场景，现成的工具更加高效。本章介绍交互式的 ATOMKIT 和工程化的 abacustest，提供“开箱即用”的解决方案。

### Section 3.1: ATOMKIT：交互式快速转换
- **内容**: 介绍 ATOMKIT 的安装与交互式菜单操作。演示如何通过简单的数字选择（如 113->101）将 CIF 文件“一键”转换为 STRU 文件，适合初学者快速上手。
- **关键参数**: 无 (涉及 CLI 交互)

### Section 3.2: abacustest：工程化输入准备
- **内容**: 介绍 abacustest 工具在准备 ABACUS 输入文件中的应用。讲解如何利用 `model inputs` 命令自动扫描目录下的赝势和轨道文件，实现批量的结构准备。
- **关键参数**: 无 (涉及 CLI 参数: `-f`, `--pp`, `--orb`)

<!-- CHAPTER_START -->
## 第四章：进阶设置：磁性与动力学约束
**本章逻辑**: 基础结构构建完成后，实际科研中往往需要处理更复杂的物理情况，如磁性体系计算和结构优化中的原子固定。本章深入 `ATOMIC_POSITIONS` 内部，讲解高级属性的设置方法。

### Section 4.1: 初始磁矩的设置 (Magnetic Moments)
- **内容**: 详细说明与 VASP (INCAR 中设置 MAGMOM) 的不同，ABACUS 的磁矩直接定义在 STRU 文件中。讲解 `mag` 关键词的使用，以及如何为反铁磁或亚铁磁体系设置不同原子的磁矩。
- **关键参数**: `mag` (位于 `ATOMIC_POSITIONS` 块内)

### Section 4.2: 结构优化中的原子约束 (Constraints)
- **内容**: 讲解如何在坐标后设置 `0` 或 `1` 来控制原子在 x, y, z 方向的移动自由度。这对于表面计算（固定底层原子）或特定方向的弛豫至关重要。
- **关键参数**: `move_x`, `move_y`, `move_z` (位于 `ATOMIC_POSITIONS` 块内)

<!-- APPENDIX_START -->
## 附录：常见问题与排错指南
- **常见报错**: 
    - 晶格常数单位混淆导致的体积错误（Bohr vs Angstrom）。
    - 缺少 `NUMERICAL_ORBITAL` 导致的 LCAO 计算崩溃。
    - 赝势/轨道文件路径未找到。
- **最佳实践**: 推荐的文件命名规范与目录结构建议。