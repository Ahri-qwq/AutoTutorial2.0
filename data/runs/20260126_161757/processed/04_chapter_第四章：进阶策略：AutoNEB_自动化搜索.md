# 第四章：进阶策略：AutoNEB 自动化搜索

在上一章中，我们了解了基础的 NEB（Nudged Elastic Band）计算流程。然而，在面对复杂的表面催化反应（如大分子解离、柔性表面重构）时，传统 NEB 往往面临两个棘手问题：一是**分辨率不足**，预设的 Image 数量可能不足以捕捉复杂的过渡态特征；二是**计算资源浪费**，大量的计算资源被消耗在远离过渡态的平坦路径上。

本章我们将引入 **AutoNEB** 方法。这是一种动态、自动化的策略，它不需要预先确定 Image 的总数，而是像“探针”一样，自动在能垒变化最剧烈的区域“加密”采样点。

> **⚠️ 核心概念预警**
>
> 请务必牢记：**ABACUS 在此模式下仅作为“Server（计算引擎）”**。
> *   **不要**在 ABACUS 的 `INPUT` 文件中寻找 `IMAGES` 或 `NEB` 相关的参数（这与 VASP 等软件不同）。
> *   所有的路径插值、弹簧力计算、收敛判断均由 **ASE** 和 **ATST-Tools** 在 Python 层面控制。
> *   ABACUS 的任务仅仅是：接收一个结构 -> 算力（Force）和能量（Energy） -> 返回给 Python 脚本。

---

## 4.1 AutoNEB 的工作逻辑与脚本配置

AutoNEB 的核心哲学是 **“先粗糙，后精细；动态加点，有的放矢”**。它不再一次性并行计算所有 Image，而是分阶段执行：

1.  **初始评估**：仅计算首末态和极少量中间点（如 1-2 个）。
2.  **动态加点**：在当前路径中，寻找**能量间隙最大**或**曲率最大**的相邻 Image 之间插入新点。
3.  **粗糙收敛**：对新插入的点进行低精度优化（较大的 `fmax`）。
4.  **循环迭代**：重复步骤 2-3，直到 Image 总数达到设定上限（`n_images`）。
5.  **最终精修**：当节点数量满足要求后，对能量最高的节点（Climbing Image）及其邻近节点进行高精度收敛计算（较小的 `fmax`）。

### 4.1.1 配置 `autoneb_run.py`

在 ATST-Tools 工作流中，`autoneb_run.py` 是控制这一过程的大脑。我们需要重点关注以下参数配置：

```python
# autoneb_run.py 核心配置片段

# --- 并行控制 ---
# n_simul: 同时计算的 Image 数量。
# 如果你的计算资源有限，不需要一次性算完所有点，AutoNEB 会自动排队。
n_simul = 4 

# --- 终止条件 ---
# n_images: 最终生成的 Image 总数（包含初末态）。
# 对于复杂反应，通常建议设为 10-16，AutoNEB 会自动加点直到达到此数值。
n_images = 12

# --- 核心收敛策略 (关键!) ---
# fmax: [粗糙收敛判据, 精细收敛判据]
# 列表的第一个值用于中间过程（节省时间），第二个值用于最终 CI-NEB 确定过渡态。
fmax = [0.20, 0.05]  # 单位: eV/Ang

# --- 算法选择 ---
algorism = "improvedtangent" # 推荐使用 IT-NEB 算法，切线估算更准
climb = True                 # 开启爬坡模式 (Climbing Image)

# --- ABACUS 计算参数 (Server 端配置) ---
# 注意：这里定义的是单个 Image 的计算精度
parameters = {
    'calculation': 'scf',
    'cal_force': 1,          # 必须开启力计算
    'cal_stress': 1,         # 推荐开启应力计算，虽不用于NEB优化，但有助于监控结构合理性
    'out_stru': 1,           # 输出结构文件
    'scf_thr': 1e-6,         # 电子步收敛精度
    'kpts': [2, 2, 1],       # K点设置
    # ... 其他 DFT 常规参数
}
```

### 4.1.2 并行资源分配策略

AutoNEB 的并行模式与传统 NEB 不同，它涉及两层并行：
1.  **Image 级并行**：由 `n_simul` 控制。
2.  **DFT 级并行**：由 `mpi` 参数（或提交脚本中的 MPI 设置）控制。

**资源计算公式**：
$$ \text{总核数} = \text{n\_simul} \times (\text{单个 ABACUS 任务所需的核数}) $$

例如，若单个单点能计算需要 16 核，你希望同时运行 4 个 Image，则提交作业时需要申请 $16 \times 4 = 64$ 个核。

---

## 4.2 实战案例：复杂表面反应 (Cy-Pt@graphene)

本案例演示环己烷（Cyclohexane）在 Pt 掺杂石墨烯表面的脱氢反应。这是一个典型的复杂体系：基底柔性大、分子构象多，传统线性插值极易导致原子重叠或路径不合理。

### Step 1: 路径插值 (Make)

在这一步，我们生成 AutoNEB 的初始猜想。

**专家建议**：
*   **插值方法**：强烈推荐使用 **IDPP** (Image Dependent Pair Potential) 方法，而非线性插值（Linear）。IDPP 能有效避免插值过程中原子距离过近导致的“核爆炸”高能结构。
*   **磁性处理**：如果体系有磁性（如本例中的 Pt-Graphene），**必须在插值阶段显式指定磁矩**。否则，中间 Image 的磁矩可能会初始化为 0，导致收敛困难或跳到非磁性态。

**操作命令**：
```bash
# -i: 指定初态和末态的输出文件 (ABACUS log 或 STRU)
# -n: 初始 Image 数量 (AutoNEB 建议从少量开始，如 4 个)
# --mag: 指定原子磁矩 (格式 元素:磁矩)
python neb_make.py -i IS/OUT.ABACUS/running_scf.log FS/OUT.ABACUS/running_scf.log \
                   -n 4 \
                   --method IDPP \
                   --mag Pt:1.5,C:0.0
```
*注：此命令会生成 `init_neb_chain.traj` 文件。*

### Step 2: 运行 AutoNEB (Run)

将配置好的 `autoneb_run.py` 提交到计算集群。

**提交脚本示例 (Slurm)**：
```bash
#!/bin/bash
#SBATCH -N 2
#SBATCH -n 64  # 假设 n_simul=4, 单个任务16核

# 必须加载 GPAW 环境，因为它是 ASE 并行计算的后端驱动
module load gpaw 
module load abacus

# 启动命令
# 注意：这里是用 mpirun 启动 python 脚本，脚本内部会再次调用 ABACUS
mpirun -np 64 gpaw python autoneb_run.py
```

**运行过程监控**：
AutoNEB 运行时会生成一系列 `run_autonebXXX.traj` 文件。
*   `run_autoneb000.traj` 到 `run_autoneb003.traj`：初始的 4 个 Image。
*   随着计算进行，你会看到 `run_autoneb004.traj` 等新文件生成，这代表 AutoNEB 正在动态插入新点。
*   查看日志 `running_autoneb.out` 可以看到类似 "Adding image between 2 and 3" 的提示。

### Step 3: 后处理与验证 (Post & Verify)

计算完成后，我们需要提取结果并验证过渡态的有效性。

**1. 生成能垒图**
```bash
python neb_post.py --autoneb run_autoneb???.traj
```
此命令会生成 `nebplots.pdf`。请重点观察：
*   **能量曲线**：是否平滑？
*   **切线力（Tangent Force）**：在最高点（TS），切线力投影应接近于 0。

**2. 核心验证：振动分析 (Vibrational Analysis)**
找到能量最高的 Image（假设是 Image 6），必须进行频率计算以确认它是真正的过渡态（鞍点）。

*   **判据**：**仅有一个虚频**（Imaginary Frequency），且该虚频对应的振动模式（Vibrational Mode）方向与反应路径方向一致。
*   **操作**：使用 ATST-Tools 的 `vib_analysis.py`。

```python
# vib_analysis.py 配置片段
# 从 NEB 结果中读取过渡态结构
neb_traj = read('neb_latest.traj', index=':')
# 自动识别过渡态索引，或者手动指定，例如 atoms = neb_traj[6]
atoms, vib_indices = neb2vib(neb_traj) 

# 只对吸附分子和表面附近的原子进行振动分析（减少计算量）
# vib_indices = [ ... ] 
```

运行振动分析后，检查输出文件：
```text
# running_vib.out 示例
---------------------
  #    meV     cm^-1
---------------------
  0   87.4i    705.0i   <-- 唯一的虚频
  1    2.4      19.3
  ...
```
如果出现多个虚频，说明该结构是高阶鞍点（Hilltop），需要微调结构后重新收敛；如果没有虚频，说明掉入了亚稳态陷阱。

---

## 本章小结

AutoNEB 策略通过动态分配计算资源，完美解决了复杂反应路径搜索中的“分辨率”与“效率”矛盾。
1.  **Input 层面**：ABACUS 仅作为力学计算器，不涉及 NEB 参数。
2.  **Workflow 层面**：`neb_make.py` (IDPP/Mag) -> `autoneb_run.py` (Dynamic) -> `vib_analysis.py` (Validation)。
3.  **策略层面**：利用 `fmax` 的分级收敛，避免在早期路径猜测较差时浪费大量 SCF 循环。

掌握了 AutoNEB，你已经具备了处理绝大多数表面催化反应过渡态搜索的能力。下一章，我们将探讨如何利用 **Dimer 方法** 进行无路径的单端过渡态搜索。