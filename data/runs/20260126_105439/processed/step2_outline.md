<!-- META_START -->
# ABACUS 实战教程：结构文件 (STRU) 构建

## 前言
- **本书逻辑**: 本教程遵循“知其然 -> 知其所以然 -> 及其法”的教学逻辑。首先深入解析 STRU 文件的物理含义与底层格式，建立对晶体结构描述的根本认知；随后介绍基于 ASE 的编程构建方法，满足灵活定制需求；最后引入 ATOMKIT 和 abacustest 等高效工具，实现从快速转换到高通量自动化的进阶操作。
- **核心知识点**: 晶体结构描述（晶格与原子坐标）、单位制转换（Bohr 与 Angstrom）、STRU 文件核心关键词、ASE 接口配置（特定分支）、格式转换工具（ATOMKIT）、高通量输入文件准备（abacustest）。

<!-- CHAPTER_START -->
## 第一章：STRU 文件的物理图景与格式解析
**本章逻辑**: 在使用工具生成文件之前，用户必须理解文件内容的物理意义。本章将解构 STRU 文件，从模拟单元的定义到原子的微观状态，为后续排查错误打下基础。

### Section 1.1: 模拟单元的定义（晶格与基组）
- **内容**: 阐述如何定义模拟体系的边界条件和计算基组。重点解释 ABACUS 特有的单位制处理方式（晶格常数缩放因子），以及在 LCAO 基组下所需的轨道文件配置。
- **关键参数**: 
    - `ATOMIC_SPECIES`: 定义元素、质量、假势。
    - `NUMERICAL_ORBITAL`: 定义数值原子轨道文件。
    - `LATTICE_CONSTANT`: 晶格常数/缩放因子（重点解释 1.889726 的来源）。
    - `LATTICE_VECTORS`: 晶格矢量。

### Section 1.2: 原子构型的描述（位置与自由度）
- **内容**: 详解原子位置的两种表达方式（分数坐标与笛卡尔坐标），以及如何设定原子的动力学性质（固定或弛豫）和磁学性质（初始磁矩）。
- **关键参数**: 
    - `ATOMIC_POSITIONS`: 原子位置块的起始标识。
    - `Direct` / `Cartesian`: 坐标类型。
    - `move_x/y/z` (数值 0 或 1): 原子移动自由度控制。
    - `mag`: 原子初始磁矩设定。

<!-- CHAPTER_START -->
## 第二章：基于 ASE 的编程构建与定制
**本章逻辑**: Python 是计算材料学的通用语言。本章介绍如何利用 ASE（Atomic Simulation Environment）的强大生态来构建或转换结构，重点解决官方版 ASE 不支持 ABACUS 的痛点。

### Section 2.1: 环境配置与 ASE-ABACUS 接口
- **内容**: 明确指出 PyPI 官方 ASE 的局限性，详细指导用户安装带有 ABACUS 接口的特定 ASE 分支（GitLab 源），并验证安装环境。
- **关键参数**: 无（主要涉及 `git clone` 和 `pip install` 命令）。

### Section 2.2: 从零构建与格式转换
- **内容**: 演示如何使用 Python 脚本从头构建简单晶体（如 FCC Al），以及如何读取 CIF/POSCAR 等通用格式并输出为 STRU 文件。
- **关键参数**: 
    - `ase.build.bulk`: 构建体材料。
    - `ase.io.read` / `ase.io.write`: 文件读写。
    - `format='abacus'`: 指定输出格式。

### Section 2.3: 进阶：自定义赝势与轨道映射
- **内容**: 解决 ASE 默认输出中赝势和轨道文件名可能不匹配的问题，演示如何通过 Python 字典显式指定特定元素的 `upf` 和 `orb` 文件路径。
- **关键参数**: 
    - `pp` (dict): 赝势映射字典。
    - `basis` (dict): 轨道映射字典。

<!-- CHAPTER_START -->
## 第三章：高效工具链：快速转换与自动化
**本章逻辑**: 针对不同场景提供最高效的解决方案。对于简单的单次转换，介绍交互式的 ATOMKIT；对于批量任务或需要自动关联库文件的场景，介绍 abacustest。

### Section 3.1: ATOMKIT：交互式结构转换
- **内容**: 介绍 ATOMKIT 工具的获取与安装，演示如何通过命令行交互界面，快速将 CIF 或 VASP POSCAR 文件转换为标准的 ABACUS STRU 文件，适合初学者快速上手。
- **关键参数**: 
    - 功能代码 `101` (Export to ABACUS)。
    - 功能代码 `113` / `175` (Import CIF / VASP)。

### Section 3.2: abacustest：自动化输入生成
- **内容**: 介绍 abacustest 库在准备输入文件时的优势，特别是其通过 `element.json` 自动匹配赝势和轨道文件的能力，彻底解决手动修改 STRU 文件路径繁琐且易错的问题。
- **关键参数**: 
    - `abacustest model inputs`: 核心命令。
    - `-f` / `--file`: 指定结构文件。
    - `--pp`: 赝势库路径。
    - `--orb`: 轨道库路径。

<!-- APPENDIX_START -->
## 附录：常见问题与进阶建议
- **常见报错**: 
    - `LATTICE_CONSTANT` 导致的单位混淆（结构放大/缩小 ~1.89 倍）。
    - 赝势/轨道文件路径错误（File not found）。
    - ASE 安装版本错误导致 `format not supported`。
- **进阶建议**: 
    - 推荐在生产环境中使用 `abacustest` 进行标准化管理。
    - 复杂磁性结构（非共线磁性）在 STRU 中的进一步设置指引。