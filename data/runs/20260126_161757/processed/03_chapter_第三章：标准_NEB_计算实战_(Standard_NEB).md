# 第三章：标准 NEB 计算实战 (Standard NEB)

在掌握了 ASE-ABACUS 接口的基础之后，我们将进入计算材料学中最具挑战性也最迷人的领域之一：**过渡态搜索（Transition State Search）**。

本章将以经典的 **Li 原子在 Si 晶格中的扩散**为例，带你跑通一个完整的 NEB（Nudged Elastic Band）计算流程。

:::danger 核心观念转变：ABACUS 仅作为 Server
在开始之前，请务必摒弃使用 VASP 等传统软件时的习惯：**不要试图在 ABACUS 的 `INPUT` 文件中寻找任何与 NEB 相关的参数**（如 `IMAGES`, `SPRING_K` 等）。

在 ASE-ABACUS 工作流中：
1.  **ASE (Client)**：负责 NEB 算法的核心逻辑（弹簧力计算、切线方向投影、构型更新）。
2.  **ABACUS (Server)**：退化为一个纯粹的“力与能量计算器”。它只负责接收原子坐标，算出能量和力，然后返回给 ASE。

因此，ABACUS 的任务永远是单点的 `scf` 计算。
:::

我们将采用 **ATST-Tools** 的标准工作流：**Make (插值) -> Run (计算) -> Post (后处理)**。

---

## 3.1 Step 1: 构建初始路径 (Make)

NEB 计算的第一步是在初态（Initial State, IS）和末态（Final State, FS）之间生成一系列中间图像（Images）。

### 3.1.1 准备初末态
确保你已经完成了初态和末态的结构优化（Relaxation）。假设文件结构如下：
```text
Li-diffu-Si/
├── IS/
│   └── OUT.ABACUS/running_scf.log  (或 STRU)
├── FS/
│   └── OUT.ABACUS/running_scf.log  (或 STRU)
└── ... (赝势与轨道文件)
```

### 3.1.2 使用 IDPP 进行插值
我们将使用 `neb_make.py` 脚本生成初猜路径。强烈推荐使用 **IDPP (Image Dependent Pair Potential)** 插值方法，而非简单的线性插值。线性插值往往会导致原子距离过近（甚至重叠），引发 SCF 不收敛或巨大的非物理力。

**运行命令：**
```bash
# -n 3 表示在初末态之间插入 3 个 Image (总共 5 个 Image)
# -m IDPP 启用 IDPP 插值优化
python neb_make.py -i IS/OUT.ABACUS/running_scf.log FS/OUT.ABACUS/running_scf.log -n 3 -m IDPP
```

**输出文件：**
- `init_neb_chain.traj`: 这是一个包含 5 个构型的二进制轨迹文件，将作为下一步计算的输入。

:::warning 磁性体系特别注意
如果你的体系包含磁性（如 Fe, Co, Ni 或缺陷态），简单的几何插值会丢失磁矩信息，导致中间图像收敛到非磁性基态，计算出错误的能垒。

**必须**在 `neb_make.py` 阶段显式指定磁矩初猜：
```bash
# 示例：指定 Fe 原子初始磁矩为 3.0 μB
python neb_make.py ... --mag Fe:3.0
```
这将确保生成的 `init_neb_chain.traj` 中携带了磁矩信息，ASE 会将其传递给 ABACUS。
:::

---

## 3.2 Step 2: 配置计算脚本 (Run)

这是本章的核心。我们需要编写 `neb_run.py` 脚本。这个 Python 脚本将定义 ABACUS 如何计算能量（Calculator）以及 ASE 如何驱动 NEB（Optimizer）。

### 3.2.1 编写 `neb_run.py`

请将以下代码保存为 `neb_run.py`，并根据注释修改路径：

```python
from ase.io import read
from ase.optimize import FIRE
from abacus_neb import AbacusNEB

# ================= 配置区域 =================
# 1. 并行设置
# 串行 NEB (DyNEB) 设为 False; 并行 NEB 设为 True
parallel = True 
mpi = 4       # 每个 Image 分配的 MPI 核数
omp = 1       # OpenMP 线程数

# 2. NEB 算法参数
climb = True             # 开启 CI-NEB (爬坡图像)，这对精确定位过渡态至关重要
algorism = "improvedtangent" # 推荐使用改进切线法 (IT-NEB)，稳定性更好
k = 0.10                 # 弹簧常数 (eV/Ang^2)，通常 0.1 即可
fmax = 0.05              # 收敛判据 (eV/Ang)，最大原子受力小于此值即停止
neb_optimizer = FIRE     # 优化器，FIRE 对 NEB 通常比 BFGS 更稳健

# 3. ABACUS Calculator 参数 (物理核心)
# 注意：这里不需要设置 INPUT 文件，ASE 会自动生成
lib_dir = "/path/to/your/PP_ORB" 
pp = {'Li': 'Li_ONCV_PBE-1.2.upf', 'Si': 'Si_ONCV_PBE-1.2.upf'}
basis = {'Li': 'Li_gga_8au_100Ry_4s1p.orb', 'Si': 'Si_gga_8au_100Ry_2s2p1d.orb'}
kpts = [2, 2, 2]  # K点设置需与初末态保持一致

parameters = {
    'calculation': 'scf',     # 核心！NEB 的每一步只需做 SCF，不要设为 relax
    'nspin': 1,               # 1: 非磁, 2: 磁性
    'xc': 'pbe',
    'ecutwfc': 100,           # 平面波截断能 (Ry)
    'ks_solver': 'genelpa',
    'smearing_method': 'gaussian',
    'smearing_sigma': 0.001,
    'mixing_type': 'broyden',
    'scf_thr': 1e-6,
    
    # 关键接口参数
    'cal_force': 1,           # 必须开启！ASE 需要读取力来更新位置
    'cal_stress': 1,          # 推荐开启，有助于避免某些晶胞应力导致的计算异常
    'out_stru': 1,            # 输出结构文件
    
    # 路径设置
    'pp': pp,
    'basis': basis,
    'pseudo_dir': lib_dir,
    'basis_dir': lib_dir,
}

# ================= 运行逻辑 =================
if __name__ == "__main__":
    # 读取插值生成的轨迹
    init_chain = read("init_neb_chain.traj", index=':')
    
    # 初始化 ABACUS NEB 对象
    neb = AbacusNEB(init_chain, parameters=parameters, parallel=parallel,
                    directory="NEBrun", mpi=mpi, omp=omp, 
                    algorism=algorism, k=k)
    
    # 开始运行
    neb.run(optimizer=neb_optimizer, climb=climb, fmax=fmax)
```

### 3.2.2 关键参数解析

*   **`calculation`: 'scf'**: 再次强调，ABACUS 在这里不能进行几何优化。ASE 会根据 ABACUS 算出的力，计算出合力（势能面梯度 + 弹簧力），然后移动原子位置，再让 ABACUS 算下一步的 SCF。
*   **`cal_force`: 1**: 如果不开启，ABACUS 不会输出力，ASE 将无法获得梯度信息，计算会报错或卡死。
*   **`climb`: True**: 开启 CI-NEB。在计算初期，Image 会受弹簧力约束；当受力收敛到一定程度后，能量最高的 Image 会“爬坡”，即沿着反应路径反向受力，从而精确收敛到鞍点（Saddle Point）。

---

## 3.3 Step 3: 运行计算 (Execute)

根据你的硬件资源，选择串行或并行模式。

### 3.3.1 模式 A：串行 NEB (DyNEB)
**适用场景**：计算资源有限，或者体系非常简单（如本例的 Li 扩散）。
**配置**：在 `neb_run.py` 中设置 `parallel = False`。
**运行命令**：
```bash
python neb_run.py
```
**原理**：ASE 会依次计算 Image 1, Image 2, Image 3... 的能量和力。虽然慢，但 DyNEB (Dynamic NEB) 算法会优先优化受力最大的 Image，节省总计算量。

### 3.3.2 模式 B：并行 NEB (MPI Parallel)
**适用场景**：生产环境标准做法。所有 Image 同时计算。
**配置**：在 `neb_run.py` 中设置 `parallel = True`。
**资源分配逻辑**：
$$ \text{总核数} = \text{Image 数量} \times \text{单次 SCF 核数 (mpi)} $$
在本例中，我们有 3 个中间 Image（初末态固定不计算），如果在脚本中设置 `mpi = 4`，则总共需要 $3 \times 4 = 12$ 个核。

**运行命令**：
你需要使用支持 MPI 的 Python 解释器（通常是 `gpaw-python` 或正确配置的 `mpi4py`）：
```bash
# -np 12: 总核数
mpirun -np 12 gpaw python neb_run.py
```
此时，ASE 会自动将 MPI 通信域分裂，为每个 Image 分配 4 个核，并行执行 3 个 ABACUS 实例。

---

## 3.4 Step 4: 监控与后处理 (Post)

计算开始后，你可以实时监控进度。

### 3.4.1 监控收敛过程
查看 `running_neb.out` 文件（由 ATST-Tools 生成）：
```text
      Step     Time          Energy          fmax
FIRE:    0 18:28:24    -7051.845345         0.733359
FIRE:    1 18:34:09    -7051.879062         0.594320
...
```
- **Energy**: 当前最高能量 Image 的能量。
- **fmax**: 当前所有 Image 中最大的受力。当 `fmax` 小于你在脚本中设置的 0.05 时，计算结束。

### 3.4.2 可视化与能垒计算
计算完成后，使用 `neb_post.py` 进行一键后处理：

```bash
python neb_post.py neb.traj
```

**输出产物**：
1.  **`nebplots.pdf`**: 包含能量-反应坐标曲线图。请检查曲线是否平滑，最高点切线是否水平（CI-NEB 特征）。
2.  **`neb_latest.traj`**: 收敛后的最终路径。

**可视化检查**：
使用 ASE-GUI 查看原子运动轨迹，确保没有发生非物理的键断裂或原子碰撞：
```bash
ase gui neb_latest.traj
```
在 GUI 中，点击 `Tools` -> `NEB` 可以看到简易的能垒图。

### 3.4.3 进阶验证：振动分析 (预告)
找到过渡态结构（能量最高点）后，工作并没有结束。你必须验证它是否是**真实的一阶鞍点**。
- **标准**：对过渡态结构进行频率计算（Vibrational Analysis），必须**有且仅有一个虚频**（Imaginary Frequency），且该虚频对应的振动模式方向与反应路径方向一致。
- 这一步可以通过 ATST-Tools 的 `vib_analysis.py` 完成，我们将在后续章节详细讲解。

---

## 3.5 专家提示：AutoNEB 的优势

如果你不知道该插入多少个 Image，或者路径非常复杂，可以使用 **AutoNEB**。
它采用“先粗糙、后精细”的策略：
1.  先用少量 Image 跑粗糙的 NEB。
2.  自动在能隙大或几何变化大的地方插入新 Image。
3.  最后进行 CI-NEB 精修。

这通常比固定 Image 数量的标准 NEB 更节省计算资源，且不易陷入伪极小值。在 ATST-Tools 中，只需使用 `autoneb_run.py` 即可启用此模式。