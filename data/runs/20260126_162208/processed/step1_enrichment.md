根据提供的资料，以下是关于 **AIMD 轨迹分析 (RDF/MSD)** 的结构化元数据报告。

## 1. 物理本质 (Physics Concepts)
- **核心物理概念**:
    - **径向分布函数 (RDF / Pair Distribution Function, PDF)**: 描述系统中粒子密度的空间相关性，即给定一个参考粒子，在距离 $r$ 处找到另一个粒子的概率。它是表征液体、非晶体短程有序和长程无序结构的重要参数。
    - **均方差位移 (MSD)**: 描述粒子随时间偏离其初始位置的程度（$\langle |r(t) - r(0)|^2 \rangle$）。MSD 与时间的关系可用于计算扩散系数。
    - **第一性原理分子动力学 (AIMD)**: 基于量子力学（DFT）计算力场，模拟原子在有限温度下的运动轨迹。
- **科学问题**:
    - 区分固体（晶体）、液体和非晶态结构。
    - 计算配位数（通过对 RDF 第一峰积分）。
    - 研究物质的动力学性质（如扩散系数、离子电导率）。

## 2. 关键输入参数 (Key Parameters)

### 2.1 ABACUS 主程序 (生成轨迹)
*需在 `INPUT` 文件中设置*

| 参数名 | 推荐值 | 物理意义 |
| :--- | :--- | :--- |
| `calculation` | `md` | **必须**。指定计算类型为分子动力学。 |
| `md_type` | `nvt` / `nve` 等 | 指定系综类型（如正则系综 NVT，微正则系综 NVE）。 |
| `md_nstep` | > 1000 (视体系而定) | MD 模拟的总步数。需保证足够长以进行统计分析。 |
| `md_dt` | 1.0 ~ 2.0 (fs) | 时间步长。 |
| `md_dumpfreq` | 1 ~ 10 | **关键参数**。控制原子坐标信息的输出频率。若设置过大，轨迹帧数太少，分析结果（RDF/MSD）粗糙；若设置过小，磁盘占用过大。 |
| `md_restartfreq`| 需查阅文档确认 | 资料7提及此参数，用于控制 `STRU` 和重启文件的输出频率，防止计算中断数据丢失。 |

### 2.2 Candela 后处理工具 (分析轨迹)
*需在 Candela 的输入文件（通常命名为 `INPUT`）中设置*

| 参数名 | 推荐值 | 物理意义 |
| :--- | :--- | :--- |
| `calculation` | `pdf` 或 `msd` | `pdf`: 计算径向分布函数；`msd`: 计算均方差位移。 |
| `geo_in_type` | `ABACUS` | **必须**。指定读取的轨迹格式为 ABACUS 输出格式。 |
| `geo_directory` | `../MD_dump` | ABACUS 输出轨迹所在的目录路径（资料中显示 ABACUS 输出文件夹通常为 `MD_dump`）。 |
| `geo_ignore` | > 0 (如 50) | **平衡期截断**。跳过轨迹开头未达到热力学平衡的帧数。 |
| `rcut` | < 0.5 * Lattice | (仅 RDF) 计算截断半径，单位 Angstrom。通常取晶格常数的一半以下。 |
| `dr` | 0.01 | (仅 RDF) 径向距离的直方图间隔，单位 Angstrom。 |
| `msd_dt` | = `md_dt` * `md_dumpfreq` / 1000 | (仅 MSD) 轨迹两帧之间的时间间隔，单位 **ps**。需根据 ABACUS 的步长和输出频率换算。 |
| `msd_natom` | 体系原子数 | (仅 MSD) 需要计算 MSD 的原子数量。 |

### **【重要】知识缺口处理**
- **多组分 RDF**: 资料中 `Candela` 的输入参数只展示了 `ntype` 和 `natom`，未详细说明如何指定计算特定原子对（如 Si-O 或 O-O）的 RDF。**撰写时需提醒读者查阅 Candela 官方文档关于 `geo_1`, `geo_2` 或原子类型选择的具体定义**，不要盲目猜测是自动识别。
- **MDTraj 读取 ABACUS**: 资料1提及 `MDtraj` 支持多种格式，但未明确指出其是否**直接**支持 ABACUS 的 `MD_dump` 目录格式，还是需要先转换为 `.xyz` 或 `.dcd`。建议以 Candela 为主流程，Python/MDtraj 为高级自定义选项。

## 3. 体系与接口配置 (System & Interfaces)

- **结构 (STRU)**:
    - 标准 ABACUS `STRU` 文件。
    - 需注意 `ATOMIC_POSITIONS` 后的初始速度 `init_vel` 设置（资料7提及），这对 MD 的初始弛豫有影响。
- **外部接口**:
    - **Candela**: 专门用于处理 ABACUS、QE、VASP 等软件 MD 轨迹的后处理工具。
    - **Python (MDTraj)**: 通用 MD 分析包，灵活性更高，适合自定义绘图。
- **接口注意事项**:
    - **ABACUS -> Candela**: ABACUS 计算完成后，会生成 `MD_dump` 文件夹（或由 `geo_directory` 指定）。Candela 直接读取该文件夹下的结构信息。
    - **单位换算**: ABACUS 的 `md_dt` 通常单位为 fs，但 Candela 计算 MSD 时输出的时间单位通常为 ps（资料4提及 `msd_dt` 单位为 ps）。**必须强调这一单位转换**。

## 4. 教程编写特殊指令 (Special Instructions for Writer)

- **Critical (流程分段)**: 教程必须清晰地分为两个阶段：
    1.  **生产阶段 (Production Run)**: 运行 ABACUS 进行 AIMD，强调 `md_dumpfreq` 的设置，否则没有轨迹可分析。
    2.  **分析阶段 (Analysis)**: 使用 Candela 读取轨迹。
- **Critical (平衡性检查)**: 在计算 RDF/MSD 之前，必须强调检查系统的能量/温度曲线，确认系统已达到平衡。在 Candela 中使用 `geo_ignore` 参数去除平衡前的轨迹是标准操作，**务必在教程中高亮此步骤**。
- **锦囊妙计 (MSD 时间步长)**: 在设置 Candela 的 `msd_dt` 时，新手极易出错。请给出一个计算公式示例：
    > 如果 ABACUS 设置 `md_dt = 2.0` (fs) 且 `md_dumpfreq = 10`，则 Candela 中的 `msd_dt` 应设为 $2.0 \times 10 / 1000 = 0.02$ ps。
- **风险提示**: 资料中关于 Candela 的 `geo_1` 和 `geo_2` 参数描述为“MD 轨迹起始/结束索引”，这与 `geo_ignore` 功能似乎有重叠或互补。建议提醒读者：`geo_ignore` 是跳过帧数，而 `geo_1`/`geo_2` 可能是指定读取文件的编号范围（如果 ABACUS 输出多个分段文件）。

## 5. 常见报错与注意事项 (Pitfalls)

- **轨迹文件缺失**: 如果 ABACUS `INPUT` 中未设置 `md_dumpfreq` 或设置过大，`MD_dump` 目录可能为空或帧数不足，导致 Candela 报错或 RDF 曲线全是噪点。
- **内存溢出**: 对于超大体系或超长轨迹，一次性读取所有帧可能导致内存溢出。
- **PBC 问题**: 计算 RDF 时，截断半径 `rcut` 不能超过模拟盒子最小边长的一半，否则会引入周期性镜像的错误统计。
- **统计不足**: AIMD 步数太少（例如仅跑了 100 步）会导致 RDF 峰形不平滑，MSD 曲线无法呈现线性扩散行为。