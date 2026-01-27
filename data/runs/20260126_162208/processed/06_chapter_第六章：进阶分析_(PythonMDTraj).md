# 第六章：进阶分析 (Python/MDTraj)

虽然 Candela 提供了便捷的快速可视化和基础分析功能，但在面对复杂的结构分析（如特定原子对的径向分布函数、配位数分析、二面角统计）或需要出版级的高质量绘图时，Python 生态系统——特别是 `MDTraj`、`ASE` 和 `Matplotlib`——是不可或缺的工具。

本章将带你从 ABACUS 的原始输出出发，通过 Python 脚本实现高自由度的轨迹分析。

---

## 6.1 使用 MDTraj 处理 ABACUS 轨迹

要进行高质量的分析，首先必须确保你的 ABACUS 模拟产出了足够的数据。我们将流程分为**生产阶段**和**分析阶段**。

### 第一阶段：生产运行 (Production Run)

在运行 ABACUS 分子动力学（AIMD）之前，必须在 `INPUT` 文件中正确设置输出频率。很多初学者遇到的最大问题是：跑了一周的 MD，结果发现轨迹文件里只有寥寥几帧，无法进行统计分析。

**关键参数设置 (INPUT 文件示例)：**

```bash
INPUT_PARAMETERS
# ... (基础 DFT 参数)

# --- MD 核心参数 ---
calculation     md          # 计算类型：分子动力学
md_type         nvt         # 系综类型：nve, nvt, npt 等
md_nstep        10000       # 总步数：必须足够长以保证统计显著性
md_dt           1.0         # 时间步长：单位 fs (通常 0.5 - 2.0)

# --- 关键输出控制 ---
md_dumpfreq     10          # 轨迹输出频率：每 10 步输出一次坐标
md_out_file     1           # 显式开启输出 (默认通常开启)
```

> **专家提示**：
> *   **`md_dumpfreq`**: 这是一个权衡参数。
>     *   设为 `1`：文件巨大，I/O 瓶颈，包含过多高频噪声。
>     *   设为 `1000`：数据太稀疏，无法捕捉短时动力学行为，RDF 曲线会充满噪点。
>     *   **推荐值**：通常设为 `10` 到 `50` 之间。
> *   **轨迹文件位置**：运行结束后，轨迹通常保存在 `OUT.${suffix}/MD_dump` 文件中（格式类似 XYZ）。

### 第二阶段：平衡性检查 (Pre-Analysis Check)

**Critical (至关重要)**：在计算任何性质（RDF, MSD）之前，**必须**确认系统已经达到热力学平衡。

1.  **检查能量与温度**：提取 `OUT.${suffix}/running_md.log` 中的温度和能量数据。
2.  **判定平衡点**：观察曲线何时不再漂移并在均值附近波动。
3.  **切片处理**：假设前 2000 步是平衡期，分析时必须丢弃这些数据。

### 第三阶段：Python 分析实战 (RDF 计算)

ABACUS 的 `MD_dump` 文件通常兼容标准的 XYZ 格式（包含元素符号和坐标）。我们可以直接使用 `MDTraj` 读取。

**Python 脚本示例：计算径向分布函数 (RDF)**

```python
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

# 1. 参数设置
traj_file = 'OUT.ABACUS/MD_dump'  # ABACUS 轨迹文件路径
topology_file = None              # 对于 XYZ，MDTraj 通常能自动推断拓扑
equilibration_steps = 200         # 需要丢弃的帧数 (注意是帧数，不是 MD 步数)

# 2. 读取轨迹
print(f"正在读取轨迹: {traj_file} ...")
# ABACUS 的 MD_dump 格式通常为 XYZ，可以直接加载
traj = md.load_xyz(traj_file)

# 3. 数据预处理：去除平衡前的轨迹 (Slicing)
# 这一步对应 Candela 中的 geo_ignore
if len(traj) > equilibration_steps:
    traj_production = traj[equilibration_steps:]
    print(f"原始帧数: {len(traj)}, 分析帧数: {len(traj_production)}")
else:
    raise ValueError("轨迹过短，不足以去除平衡步！")

# 4. 拓扑修正 (如果 MDTraj 无法识别元素)
# 有时需要手动指定元素，或者加载一个 PDB 文件作为拓扑
# table, bonds = traj_production.topology.to_dataframe()
# print(table.head()) # 检查元素识别是否正确

# 5. 计算 RDF (以所有原子对为例)
# select_pairs 接受原子索引对。这里我们选取所有原子对。
pairs = traj_production.topology.select_pairs('all', 'all')
r_min, r_max = 0.0, 1.0 # 单位 nm (MDTraj 默认使用纳米)
r, g_r = md.compute_rdf(traj_production, pairs, r_range=(r_min, r_max))

# 6. 绘图
plt.figure(figsize=(8, 5))
plt.plot(r * 10, g_r, label='Total RDF', linewidth=2) # 转换为 Angstrom
plt.xlabel('r ($\AA$)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.title('Radial Distribution Function (from ABACUS MD)', fontsize=16)
plt.grid(alpha=0.3)
plt.legend()
plt.savefig('rdf_analysis.png', dpi=300)
print("RDF 图像已保存为 rdf_analysis.png")
```

---

## 附录：常见问题与避坑指南

在处理 ABACUS 轨迹时，以下是新手最容易踩的坑，请务必仔细阅读。

### 1. PBC 边界效应与 `rcut`
*   **问题**：计算 RDF 时，为什么截断半径 (`rcut` 或 `r_max`) 不能超过盒子长度的一半？
*   **原理**：在周期性边界条件 (PBC) 下，如果统计半径超过 $L/2$，粒子可能会统计到它自己的镜像（Self-interaction），导致 $g(r)$ 在远端出现错误的峰值。
*   **建议**：始终设置 `r_max` 小于最小盒长的一半。

### 2. 内存溢出 (Memory Error)
*   **场景**：`MD_dump` 文件达到数 GB，Python 脚本直接卡死。
*   **解决方案**：
    *   **MDTraj**: 使用 `md.iterload()` 分块读取。
    *   **降采样**: 在读取时使用 `stride` 参数（例如每 10 帧读 1 帧），虽然牺牲了部分统计量，但能显著降低内存占用。

### 3. 统计不足 (Noisy RDF)
*   **现象**：RDF 曲线毛刺极多，不平滑。
*   **原因**：采样点太少。100 步 MD 产生的构型不足以遍历相空间。
*   **对策**：增加 `md_nstep`。对于液体或非晶体系，通常需要数千甚至数万帧的数据才能得到平滑的 RDF。

### 4. 参数映射速查表 (ABACUS vs. Candela)

在使用 Candela 进行快速分析时，参数设置必须与 ABACUS 的 `INPUT` 保持物理一致性。

| 参数类别 | ABACUS (INPUT) | Candela (分析工具) | 物理含义/备注 |
| :--- | :--- | :--- | :--- |
| **时间步长** | `md_dt` (fs) | `msd_dt` (ps) | **极易出错！** 见下方锦囊 |
| **输出频率** | `md_dumpfreq` | (隐式影响) | 决定了 Candela 读取的帧间隔时间 |
| **平衡去除** | (手动检查日志) | `geo_ignore` | **跳过**文件开头的帧数 (Equilibration) |
| **文件分段** | (自动分段) | `geo_1`, `geo_2` | 指定读取的文件**编号**范围 |

#### 💡 锦囊妙计：MSD 时间步长计算公式

在 Candela 或 Python 中计算均方位移 (MSD) 时，必须输入正确的时间间隔 (`dt`)，否则扩散系数计算将完全错误。

**公式**：
$$ \text{Analysis\_dt (ps)} = \frac{\text{md\_dt (fs)} \times \text{md\_dumpfreq}}{1000} $$

**示例**：
如果你在 ABACUS 中设置：
*   `md_dt = 2.0` (fs)
*   `md_dumpfreq = 10`

那么每一帧之间的时间间隔是 20 fs。
在 Candela 中设置 `msd_dt` 时，应填入：
$$ 2.0 \times 10 / 1000 = \mathbf{0.02} \text{ ps} $$

#### ⚠️ 风险提示：`geo_ignore` vs `geo_1/2`
*   **`geo_ignore`**: 用于**物理上的去平衡**。例如 `geo_ignore = 100` 表示忽略读取到的前 100 帧数据，防止非平衡态数据污染结果。
*   **`geo_1` / `geo_2`**: 用于**文件读取控制**。如果 ABACUS 运行时间很长，输出了 `MD_dump_1`, `MD_dump_2` ... `MD_dump_5`。
    *   `geo_1 = 2`, `geo_2 = 4` 表示只读取 `MD_dump_2` 到 `MD_dump_4` 这些文件。
    *   **切勿混淆**：不要试图用 `geo_1` 来去除平衡，除非你的平衡过程恰好占满了一整个文件。标准做法是读取所有文件，然后用 `geo_ignore` 切掉前段。