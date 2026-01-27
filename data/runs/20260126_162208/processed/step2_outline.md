<!-- META_START -->
# ABACUS 实战教程：AIMD 轨迹分析 (RDF/MSD)

## 前言
- **本书逻辑**: 本教程遵循“理论认知 -> 数据生产 -> 质量控制 -> 数据分析 -> 进阶处理”的科研工作流。旨在帮助读者从物理本质出发，掌握使用 ABACUS 生成高质量 AIMD 轨迹，并利用 Candela 及 Python 工具提取结构（RDF）与动力学（MSD）信息的全套技能。
- **核心知识点**: 第一性原理分子动力学 (AIMD)、径向分布函数 (RDF)、均方差位移 (MSD)、轨迹采样频率控制、热力学平衡判定、Candela 参数配置。

<!-- CHAPTER_START -->
## 第一章：物理基础与分析目标
**本章逻辑**: 在动手计算前，必须明确“算什么”以及“为什么算”。本章简述 RDF 和 MSD 的物理意义，建立结构与动力学的基本认知框架。

### Section 1.1: 结构有序性的度量——RDF
- **内容**: 解释径向分布函数 (RDF/PDF) 的定义及其在区分晶体、液体、非晶体中的作用。介绍如何通过 RDF 第一峰积分计算配位数。
- **关键参数**: 无

### Section 1.2: 扩散行为的表征——MSD
- **内容**: 解释均方差位移 (MSD) 的时间演化规律，以及如何通过爱因斯坦关系式（Einstein Relation）将 MSD 斜率关联到扩散系数。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第二章：AIMD 轨迹生产 (ABACUS)
**本章逻辑**: 高质量的分析始于正确的模拟。本章聚焦于如何在 ABACUS 中设置 MD 参数，特别是如何平衡“采样精度”与“存储开销”，为后续分析打下基础。

### Section 2.1: 基础 MD 参数设置
- **内容**: 讲解如何在 `INPUT` 文件中开启 MD 模式，选择系综（NVT/NVE），并设置合理的步长与总步数以保证统计显著性。
- **关键参数**: `calculation` (md), `md_type`, `md_nstep`, `md_dt`

### Section 2.2: 关键采样策略：输出频率控制
- **内容**: **重点章节**。详细阐述 `md_dumpfreq` 的物理含义及其对后续 RDF/MSD 分析精度的影响。讨论过稀疏采样（分析粗糙）与过密集采样（存储爆炸）的权衡。同时提及断点续算的重要性。
- **关键参数**: `md_dumpfreq`, `md_restartfreq`

<!-- CHAPTER_START -->
## 第三章：预分析与质量控制 (Quality Control)
**本章逻辑**: 拿到轨迹后的第一步绝不是直接分析，而是检查数据的可靠性。本章强调“热力学平衡”的重要性，避免对未弛豫的垃圾数据进行分析。

### Section 3.1: 热力学平衡判定
- **内容**: 讲解如何通过能量、温度随时间的波动曲线判断体系是否达到平衡。
- **关键参数**: 无 (需分析 `running_scf.log` 或 `OUT.ABACUS` 目录下的能量文件)

### Section 3.2: 轨迹截断策略
- **内容**: 确定从哪一帧开始统计数据。引入“平衡期截断”的概念，为 Candela 的 `geo_ignore` 参数设置提供依据。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第四章：结构分析实战 (Candela RDF)
**本章逻辑**: 进入具体的后处理环节。首先进行静态结构分析，讲解 Candela 工具的具体配置与运行。

### Section 4.1: Candela 输入文件配置 (RDF)
- **内容**: 详解 Candela 的 `INPUT` 文件设置。重点讲解如何指定 ABACUS 轨迹路径、截断半径的选择以及如何跳过非平衡帧。
- **关键参数**: `calculation` (pdf), `geo_in_type` (ABACUS), `geo_directory`, `geo_ignore`, `rcut`, `dr`

### Section 4.2: 多组分体系的 RDF 挑战
- **内容**: 针对多元素体系（如 SiO2），讨论如何计算特定原子对（如 Si-O）的 RDF。
- **关键参数**: `ntype`, `natom`, `geo_1`/`geo_2` (需查阅文档确认具体原子对选择方式)

<!-- CHAPTER_START -->
## 第五章：动力学分析实战 (Candela MSD)
**本章逻辑**: 进行动态性质分析。本章的核心难点在于“时间单位的换算”，这是新手最容易出错的地方。

### Section 5.1: MSD 参数配置与单位陷阱
- **内容**: 详解 Candela 计算 MSD 的设置。**核心难点**：演示如何根据 ABACUS 的 `md_dt` (fs) 和 `md_dumpfreq` 计算出 Candela 所需的 `msd_dt` (ps)。
- **关键参数**: `calculation` (msd), `msd_dt`, `msd_natom`

### Section 5.2: 结果解读与扩散系数计算
- **内容**: 分析 Candela 输出的 `MSD_total.txt`，识别线性扩散区，排除弹道输运区和统计噪声区，估算扩散系数。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第六章：进阶分析 (Python/MDTraj)
**本章逻辑**: Candela 功能固定，为了满足个性化绘图和复杂分析需求，引入 Python 生态工具作为补充。

### Section 6.1: 使用 MDTraj 处理 ABACUS 轨迹
- **内容**: 介绍如何利用 Python 的 `mdtraj` 库读取轨迹。讨论格式兼容性（直接读取目录 vs 转换格式），并展示计算 RDF 的 Python 代码示例。
- **关键参数**: 无 (Python 脚本编写)

<!-- APPENDIX_START -->
## 附录：常见问题与避坑指南
- **PBC 边界效应**: 解释为什么 `rcut` 不能超过盒长一半。
- **内存溢出处理**: 针对超大轨迹文件的处理建议。
- **统计不足**: 为什么 100 步的 MD 跑不出平滑的 RDF？
- **参数速查表**: ABACUS 与 Candela 对应参数映射表。