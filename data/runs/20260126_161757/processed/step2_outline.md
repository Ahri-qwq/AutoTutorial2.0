<!-- META_START -->
# 《ABACUS 实战教程：NEB 过渡态搜索与反应动力学模拟》

## 前言
- **本书逻辑**: 本教程遵循“理论基础 -> 前处理与插值 -> 计算策略执行 -> 结果验证与分析”的科研工作流。旨在帮助读者理解如何利用 ASE 的 NEB 算法驱动 ABACUS 进行高精度的过渡态搜索，特别强调了从传统的 NEB 到高效的 AutoNEB 的进阶，以及必不可少的振动频率验证环节。
- **核心知识点**: 过渡态理论 (TST)、最小能量路径 (MEP)、CI-NEB (爬坡图像)、AutoNEB (动态自动 NEB)、IDPP 插值、振动分析 (Vibrational Analysis)、自由能校正。

<!-- CHAPTER_START -->
## 第一章：过渡态理论与 NEB 算法原理
**本章逻辑**: 在动手计算前，必须建立正确的物理图像。本章从化学反应的势能面出发，解释过渡态的本质，并引出寻找过渡态的核心算法——NEB 及其变体，为后续参数设置提供理论依据。

### Section 1.1: 势能面与过渡态本质
- **内容**: 阐述化学反应过程中的势能面 (PES)、反应坐标、最小能量路径 (MEP) 以及过渡态 (TS) 作为一阶鞍点 (Saddle Point) 的数学定义。
- **关键参数**: 无

### Section 1.2: NEB 及其进阶算法 (CI-NEB & AutoNEB)
- **内容**: 详解 NEB (Nudged Elastic Band) 的工作原理（弹簧力与切线力的投影）。重点区分 IT-NEB (改进切线法)、CI-NEB (爬坡图像法，用于消除鞍点误差) 以及 AutoNEB (动态加点策略) 的适用场景与优势。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第二章：计算准备与初猜链生成 (Pre-processing)
**本章逻辑**: “垃圾输入导致垃圾输出”。NEB 计算的成功率高度依赖于初末态的结构精度和插值质量。本章讲解如何准备高质量的输入数据，这是防止计算发散的关键一步。

### Section 2.1: ATST-Tools 环境与工作流概览
- **内容**: 介绍 ATST-Tools 脚本集 (`neb_make.py`, `neb_run.py`, `neb_post.py`) 的功能定位，以及 ASE 如何作为驱动器调用 ABACUS 进行计算。
- **关键参数**: 无

### Section 2.2: 初末态结构优化 (IS & FS Relaxation)
- **内容**: 强调在进行 NEB 之前，必须对反应物 (Initial State) 和产物 (Final State) 进行严格的结构优化。原子映射 (Atom Mapping) 的重要性，确保初末态原子顺序一致。
- **关键参数**: 
    - `calculation`: 'relax' (初末态优化)
    - `scf_thr`: 1e-6 或 1e-7 (高精度收敛)
    - `cal_force`: 1

### Section 2.3: 初始路径插值与 IDPP 方法
- **内容**: 使用 `neb_make.py` 生成初始猜测链。对比线性插值与 IDPP (Image Dependent Pair Potential) 插值的差异，演示如何使用 IDPP 避免原子重叠和不合理的键长。
- **关键参数**: 
    - `n_images` (映像数量)
    - `--method`: IDPP (推荐插值方法)
    - `--fix`: (固定原子设置)

<!-- CHAPTER_START -->
## 第三章：标准 NEB 计算实战 (Standard NEB)
**本章逻辑**: 以简单的扩散案例（如 Li 在 Si 中扩散）入手，讲解最基础的 NEB 运行方式。涵盖串行 (DyNEB) 与并行计算的配置，重点解析 ASE-ABACUS 接口参数的设置。

### Section 3.1: 构建 ASE-ABACUS 计算接口
- **内容**: 编写或修改 `neb_run.py` 脚本。详细讲解如何在 Python 脚本中配置 ABACUS 的 Calculator 参数，使其仅作为能量和力的计算引擎。
- **关键参数**: 
    - `calculation`: 'scf' (NEB 过程中的每一步仅做 SCF)
    - `cal_force`: 1 (必须开启，ASE 需要读取力)
    - `cal_stress`: 1 (推荐开启)
    - `kpts`, `pp`, `basis` (需与初末态一致)

### Section 3.2: 运行串行与并行 NEB
- **内容**: 演示如何运行 NEB 计算。区分串行 DyNEB (适合资源受限或简单体系) 与 MPI 并行 NEB (适合服务器环境)。讲解 `climb` 和 `optimizer` 的选择。
- **关键参数**: 
    - `climb`: True (开启 CI-NEB)
    - `k`: 0.10 (弹簧常数)
    - `fmax`: 0.05 (收敛判据)
    - `algorism`: 'improvedtangent' (IT-NEB)

### Section 3.3: 过程监控与初步可视化
- **内容**: 解读 `running_neb.out` 输出文件，监控收敛过程（Energy & fmax）。使用 ASE-GUI 初步查看 `neb.traj` 轨迹，判断路径是否合理。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第四章：进阶策略：AutoNEB 自动化搜索
**本章逻辑**: 针对复杂反应（如表面催化解离），传统 NEB 容易因分辨率不足或计算资源分配不均而失败。本章引入 AutoNEB 方法，演示如何通过动态加点策略高效解决复杂体系的过渡态搜索。

### Section 4.1: AutoNEB 的工作逻辑与脚本配置
- **内容**: 解析 AutoNEB 的多阶段流程（粗糙计算 -> 寻找最大能隙/曲率 -> 插入新点 -> 精细计算）。配置 `autoneb_run.py`。
- **关键参数**: 
    - `n_simul`: (同时计算的映像数)
    - `n_images`: (最终目标映像总数)
    - `fmax`: [0.20, 0.05] (分阶段收敛判据)

### Section 4.2: 实战案例：复杂表面反应
- **内容**: 以 Cy-Pt@graphene 案例为例，展示 AutoNEB 如何自动修正路径，对比传统 NEB 可能出现的错误结果。
- **关键参数**: 无

<!-- CHAPTER_START -->
## 第五章：结果验证与热力学分析 (Post-processing)
**本章逻辑**: 找到能量最高点并不代表找到了真正的过渡态。本章强调科学严谨性，介绍如何通过振动分析确认虚频，并进一步计算反应能垒和自由能。

### Section 5.1: 数据后处理与能垒计算
- **内容**: 使用 `neb_post.py` 提取最终路径 `neb_latest.traj`，生成势能曲线图 (`nebplots.pdf`)，读取正逆反应能垒。
- **关键参数**: 无

### Section 5.2: 振动分析与虚频确认 (Vibrational Analysis)
- **内容**: 这里的核心是验证。使用 `vib_analysis.py` 对过渡态结构进行有限差分计算。如何通过结果确认**仅存在一个沿反应坐标方向的虚频**。
- **关键参数**: 
    - `vib_indices`: (指定需要振动的原子索引，通常仅包含反应中心原子)
    - `delta`: 0.01 (有限差分步长)

### Section 5.3: 热力学数据校正
- **内容**: 基于振动频率数据，计算零点能 (ZPE) 校正，并结合温度计算熵 (Entropy) 和 自由能 (Free Energy)，获得更接近实验条件的反应动力学数据。
- **关键参数**: 
    - `T`: (温度，单位 K)

<!-- APPENDIX_START -->
## 附录：常见问题与进阶建议
- **常见报错**: 原子乱序导致的能量爆炸、SCF 不收敛的解决技巧、并行核数分配错误。
- **进阶方向**: 单端搜索方法 (Dimer/Sella) 简介、从 NEB 到 Dimer 的联用策略 (D2S)。