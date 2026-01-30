<!-- META_START -->
# ABACUS 实战教程：基于 abacustest 的自动化输入文件构建

## 前言
- **本书逻辑**: 本教程旨在解决计算材料学初学者在ABACUS输入文件准备过程中面临的“繁琐”与“易错”痛点。我们将摒弃传统的手动编辑模式，直接引入 `abacustest` 这一高效工具，通过“环境配置 -> 基础案例 -> 进阶电子结构 -> 几何优化与定制 -> 批量处理”的实战路径，让读者掌握从结构文件（CIF/POSCAR）一键生成标准化ABACUS输入文件夹（INPUT, STRU, KPT, 赝势/轨道）的核心技能。
- **核心知识点**: 
    - 自动化工作流工具 `abacustest model inputs` 的使用
    - 赝势与数值原子轨道库（APNS-pp-orb-v1）的配置
    - 基础SCF计算、磁性设置与DFT+U参数的自动化生成
    - 几何优化（Relax/Cell-Relax）与自定义参数模板
    - 高通量批量任务的输入文件准备

<!-- CHAPTER_START -->
## 第一章：环境准备与资源配置
**本章逻辑**: “工欲善其事，必先利其器”。在进行具体的计算输入文件生成前，必须先准备好核心工具 `abacustest` 以及ABACUS计算所需的赝势和轨道库。这是所有后续操作的基础。

### Section 1.1: 赝势与轨道库的获取
- **内容**: 介绍如何使用 `abacustest` 命令行工具下载官方推荐的 `APNS-pp-orb-v1` 库。解析下载后的目录结构（包含 `apns-pseudopotentials-v1` 和 `apns-orbitals-efficiency-v1`），说明其包含的元素范围（H到Bi）及轨道精度（DZP）。
- **关键参数**: 
    - 命令参数：`--download-pporb`

### Section 1.2: 环境变量配置
- **内容**: 讲解如何配置环境变量，使 `abacustest` 能够自动识别并调用下载好的赝势和轨道文件，避免每次命令都需要手动指定路径。同时说明手动指定路径的备选方案。
- **关键参数**: 
    - 环境变量：`ABACUS_PP_PATH`, `ABACUS_ORB_PATH`
    - 备选参数：`--pp`, `--orb`

<!-- CHAPTER_START -->
## 第二章：基础实战——MgO的SCF计算准备
**本章逻辑**: 以最简单的非磁性绝缘体MgO为例，演示从CIF文件到完整ABACUS输入文件夹的“Hello World”流程。本章重点在于理解标准工作流和生成文件的结构。

### Section 2.1: 一键生成输入文件夹
- **内容**: 演示使用 `MgO.cif` 文件生成SCF计算任务的完整命令。详细解析命令中各个Flag的含义，包括指定输入文件、格式、基组类型以及输出文件夹命名规则。
- **关键参数**: 
    - 命令参数：`-f`, `--ftype`, `--lcao`, `--folder-syntax`

### Section 2.2: 产物解析与文件结构
- **内容**: 深入剖析生成的 `MgO` 文件夹内部结构。解释 `INPUT` 文件的默认参数设置（如 `ecutwfc`, `mixing_type` 等），`STRU` 文件的自动生成逻辑，以及赝势和轨道文件的软链接机制。
- **关键参数**: 
    - INPUT参数：`calculation` (scf), `basis_type` (lcao), `kspacing`
    - 辅助文件：`run.sh`, `setting.json` (简要说明)

<!-- CHAPTER_START -->
## 第三章：进阶实战——Fe2O3的磁性与强关联体系
**本章逻辑**: 实际科研中常涉及磁性和强关联体系。本章以反铁磁性材料Fe2O3为例，讲解如何通过命令行参数自动化处理自旋极化（Spin）、初始磁矩设置以及DFT+U参数，避免手动修改 `STRU` 和 `INPUT` 的繁琐过程。

### Section 3.1: 自旋极化与初始磁矩设置
- **内容**: 演示如何在命令中开启共线自旋极化（nspin=2），并为Fe原子设置初始磁矩（4 μB）。分析生成的 `STRU` 文件中原子磁矩的具体格式，以及 `INPUT` 文件中混合参数（mixing_beta）的自动调整。
- **关键参数**: 
    - 命令参数：`--nspin`, `--init_mag`
    - INPUT参数：`nspin`, `mixing_beta`, `onsite_radius`, `out_mul`

### Section 3.2: 自动化DFT+U配置
- **内容**: 讲解如何通过命令行参数直接开启DFT+U功能，并为特定元素（如Fe）指定Hubbard U值（Ueff）。检查生成的 `INPUT` 文件中关于Hubbard U和uramping的设置。
- **关键参数**: 
    - 命令参数：`--dftu`, `--dftu_param`
    - INPUT参数：`dft_plus_u`, `hubbard_u`, `orbital_corr`, `uramping`

<!-- CHAPTER_START -->
## 第四章：任务类型切换与参数定制
**本章逻辑**: 除了默认的SCF计算，用户常需进行结构优化或使用特定的计算参数。本章介绍如何通过 `abacustest` 快速切换计算任务类型（如Relax），以及如何通过模板文件自定义 `INPUT` 和 `KPT` 设置。

### Section 4.1: 结构优化任务（Relax/Cell-Relax）
- **内容**: 以MgO为例，演示如何将任务类型切换为晶胞优化（cell-relax）。解析生成的 `INPUT` 文件中关于力收敛、应力收敛及优化算法的默认设置。
- **关键参数**: 
    - 命令参数：`--jtype`
    - INPUT参数：`calculation` (cell-relax), `relax_method`, `cal_force`, `cal_stress`, `force_thr_ev`, `stress_thr`

### Section 4.2: 自定义INPUT模板与K点设置
- **内容**: 以Fe2O3为例，演示当默认参数不满足需求时（例如需要更密的K点或更严格的收敛标准），如何通过 `INPUT_template` 文件和 `--kpt` 参数来覆盖默认设置。
- **关键参数**: 
    - 命令参数：`--input`, `--kpt`
    - INPUT参数：`smearing_sigma`, `kspacing` (被KPT文件替代)

<!-- CHAPTER_START -->
## 第五章：高通量场景——批量结构处理
**本章逻辑**: 在材料计算中，经常需要对一系列相似结构进行相同的计算（如表面能测试）。本章介绍如何利用通配符和Python语法，一次性为多个结构文件生成对应的ABACUS任务文件夹。

### Section 5.1: 批量生成任务文件夹
- **内容**: 以Pd(100)面不同层数的表面能计算为例（Pd100_*layer.vasp），演示如何使用通配符匹配多个文件。重点讲解 `--folder-syntax` 中使用的Python字符串切片语法（如 `"x[:-5]"`），实现文件夹的自动化规范命名。
- **关键参数**: 
    - 命令参数：`-f` (配合通配符 *), `--folder-syntax` (配合Python语法)

<!-- APPENDIX_START -->
## 附录：常见问题与进阶建议
- **内容**: 
    - 赝势/轨道文件软链接失效的排查（--copy-pp-orb 选项的使用）
    - 不同结构格式（CIF, POSCAR等）的兼容性说明
    - 如何查看 `abacustest` 的帮助文档获取更多高级参数