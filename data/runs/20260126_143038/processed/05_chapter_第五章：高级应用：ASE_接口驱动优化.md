# 第五章：高级应用：ASE 接口驱动优化

在前面的章节中，我们已经掌握了如何通过直接编辑 `INPUT` 文件来运行 ABACUS 的结构优化。然而，在面对复杂的科研场景时——例如高通量计算、过渡态搜索（NEB）、或者需要严格保持晶体对称性的相变研究——仅靠 ABACUS 原生的优化器可能不够灵活。

本章将介绍如何利用 Python 的 **ASE (Atomic Simulation Environment)** 库作为“驱动器（Driver）”，调用 ABACUS 作为“计算引擎（Calculator）”。这种模式不仅能极大扩展计算的自由度，还能有效解决原生输入难以处理的复杂对称性约束问题。

---

## 5.1 ASE 与 ABACUS 的交互机制

### 5.1.1 Calculator 模式：谁在掌舵？

理解 ASE 接口的核心在于区分“谁在移动原子”。

1.  **原生模式 (Native Mode)**：
    *   **掌舵者**：ABACUS 程序本身。
    *   **设置**：在 `INPUT` 中设置 `calculation relax` 或 `cell-relax`。
    *   **流程**：ABACUS 读取结构 -> 计算力 -> ABACUS 移动原子 -> 循环直到收敛。
    *   **优点**：简单，无需编写 Python 脚本，效率略高（无 I/O 开销）。

2.  **ASE 驱动模式 (ASE Driver Mode)**：
    *   **掌舵者**：Python 脚本中的 ASE 优化器（如 `BFGS`, `FIRE`）。
    *   **引擎**：ABACUS。
    *   **设置**：在 Python 脚本中，ABACUS 的 `calculation` 参数通常设为 `scf`（单点能计算）。
    *   **流程**：
        1. ASE 给定原子位置。
        2. ASE 调用 ABACUS 计算单点能（SCF），获取能量、力和应力。
        3. ASE 根据力，利用自身的算法（如 BFGS）计算下一步原子位置。
        4. ASE 更新位置，再次调用 ABACUS。
    *   **优点**：可以使用 ASE 强大的生态（如声子谱、NEB、遗传算法），且能通过 Python 代码对每一步进行精细控制（如施加特殊的对称性约束）。

### 5.1.2 基础调用示例

以下是一个典型的 Python 脚本框架，展示了如何使用 ASE 调用 ABACUS 进行结构优化。

> **注意**：运行此脚本前，请确保环境已安装 `ase` 库，并配置好 ABACUS 的环境变量。

```python
import os
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.abacus import Abacus, AbacusProfile

# 1. 设置 ABACUS 运行环境
# 请根据实际机器修改 command，例如 mpirun -np 4 abacus
abacus_cmd = "mpirun -np 8 abacus" 
# 指向包含赝势和轨道文件的目录
pp_orb_dir = "./" 

profile = AbacusProfile(command=abacus_cmd)

# 2. 准备结构 (可以从文件读取，也可以用 ASE 建模)
atoms = read("STRU", format="abacus")

# 3. 设置计算器参数 (对应 INPUT 文件)
# 注意：这里 calculation 设为 'scf'，因为优化循环由 ASE 控制
calc = Abacus(
    profile=profile,
    directory='run_ase_opt',
    pp_dir=pp_orb_dir,
    orb_dir=pp_orb_dir,
    suffix='ase_driver',
    calculation='scf',       # 关键：只让 ABACUS 算单点能和力
    ecutwfc=100,             # Ry
    scf_thr=1e-6,
    basis_type='lcao',
    ks_solver='genelpa',
    smearing_method='gaussian',
    smearing_sigma=0.01,
    kpts=[4, 4, 4]           # 警告：实际计算请务必进行 K 点收敛测试
)

atoms.calc = calc

# 4. 设置 ASE 优化器
# logfile 记录优化过程，trajectory 保存轨迹
opt = BFGS(atoms, logfile='opt.log', trajectory='opt.traj')

# 5. 运行优化
# fmax 单位为 eV/Ang
opt.run(fmax=0.01)

# 6. 保存最终结构
atoms.write('STRU_optimized', format='abacus')
```

---

## 5.2 强对称性约束 (FixSymmetry)

在材料计算中，一个常见且棘手的问题是**对称性破缺**。

### 5.2.1 物理背景：BCC 到 FCC 的“意外”演化
假设你正在研究体心立方（BCC）铁在极端条件下的性质。如果你使用可变晶胞优化（Variable Cell Relaxation），由于数值噪音或初始结构的微小扰动，BCC 结构可能会沿着某个软模方向“滑移”，最终演化为面心立方（FCC）结构或更低对称性的结构。

虽然这在物理上可能是正确的（找到了更低能量的基态），但如果你**仅仅想研究 BCC 结构下的性质**（例如计算 BCC 下的声子谱或弹性常数），这种自动演化就是灾难性的。

### 5.2.2 解决方案：ASE FixSymmetry
ASE 提供了 `FixSymmetry` 类，它可以在优化的每一步强制清洗掉破坏对称性的力和应力分量，确保结构始终保持在指定的空间群中。

### 5.2.3 实战代码

我们需要结合 `UnitCellFilter`（用于优化晶胞）和 `FixSymmetry`（用于约束对称性）。

```python
from ase.constraints import FixSymmetry, UnitCellFilter
from ase.spacegroup import get_spacegroup

# ... (前文的 atoms 和 calc 设置代码同上) ...

# 1. 分析并施加对称性约束
# symprec 是对称性识别精度，建议设为 1e-5 或更小
atoms.set_constraint(FixSymmetry(atoms, symprec=1e-5))

# 打印当前识别出的空间群，确认是否符合预期
sg = get_spacegroup(atoms, symprec=1e-5)
print(f"Current Spacegroup: {sg.symbol} ({sg.no})")

# 2. 设置晶胞过滤器
# 如果要优化晶胞常数 (相当于 cell-relax)，必须使用 UnitCellFilter
# scalar_pressure 设置外压，单位通常为 eV/Ang^3 (ASE单位)，需注意换算
ucf = UnitCellFilter(atoms, scalar_pressure=0.0)

# 3. 将过滤器传入优化器，而不是直接传入 atoms
opt = BFGS(ucf, logfile='opt_sym.log')

# 4. 运行优化
# 此时，无论优化多少步，BCC 结构都不会变成 FCC，
# 只有保持 BCC 对称性的晶格常数和原子位置会被优化。
opt.run(fmax=0.01)
```

**专家提示**：
*   **初始结构很重要**：如果初始结构本身已经严重偏离了高对称性点，`FixSymmetry` 可能无法正确识别你想要的空间群。建议在建模时保证坐标的精确性。
*   **UnitCellFilter**：在 ASE 中，普通的 `BFGS(atoms)` 默认只优化原子坐标（固定晶胞）。如果你需要优化晶格常数（对应 ABACUS 的 `cell-relax`），必须使用 `UnitCellFilter` 或 `ExpCellFilter` 包装 atoms 对象。

---

<!-- APPENDIX_START -->
## 附录：常见报错与速查表

在进行结构优化时，参数设置失误是导致计算失败或结果不可信的主要原因。以下是针对 ABACUS 原生优化参数的速查表与避坑指南。

### 1. `relax` vs `cell-relax` 核心区分

这是初学者最容易混淆的概念，选错会导致物理意义完全不同的结果。

| 参数值 (`calculation`) | 物理含义 | 优化对象 | 典型应用场景 |
| :--- | :--- | :--- | :--- |
| **`relax`** | 几何优化 (Geometry Optimization) | **仅原子位置** (晶胞形状和体积固定) | 表面计算、缺陷体系、分子吸附、NEB 初态末态构建。 |
| **`cell-relax`** | 变胞优化 (Variable Cell Relaxation) | **原子位置 + 晶胞形状 + 体积** | 块体材料基态搜索、相变研究、寻找平衡晶格常数。 |

> **场景引导**：如果你在做一个表面吸附计算（构建了真空层），千万**不要**使用 `cell-relax`，否则真空层会因为没有原子支撑而坍缩，或者晶胞会发生非物理的畸变。

### 2. 常见“坑点”与解决方案

#### (1) `relax_nmax` 设为 1 的陷阱
*   **现象**：计算运行了一步离子步后就提示“收敛”或直接结束，但受力依然很大。
*   **原因**：在某些测试脚本或教程中，为了快速演示流程，会将 `relax_nmax`（最大离子迭代步数）设为 1。
*   **修正**：在生产计算中，**务必将 `relax_nmax` 设置得足够大**。
    *   推荐值：`relax_nmax 100` 或 `200`。
    *   如果 100 步仍未收敛，通常意味着初始结构太差或收敛标准过高，而不是步数不够。

#### (2) 应力单位 (`stress_thr`)
*   **现象**：设置了 `stress_thr 0.01`，结果计算永远不收敛，或者精度极差。
*   **原因**：ABACUS 的 `INPUT` 文件中，`stress_thr` 的默认单位是 **kBar**，而不是 eV/Ang³ 或 GPa。
*   **换算**：1 GPa = 10 kBar。
*   **推荐值**：
    *   粗略优化：`stress_thr 10` (1 GPa)
    *   精确优化：`stress_thr 1` (0.1 GPa)

#### (3) 赝势与 LCAO 轨道文件不匹配
*   **现象**：程序启动时报错，或者 SCF 第一步能量就是 `NaN`。
*   **原因**：使用 LCAO 基组时，`orbital_dir` 下的 `.orb` 文件必须与 `pseudo_dir` 下的 `.upf/.vwr` 文件是同一套生成参数制作的。
*   **修正**：始终从官方网站或可靠来源下载**成对**的赝势和轨道文件，不要混用不同版本的库。

#### (4) K 点采样的收敛性
*   **风险**：教程示例为了运行速度快，常使用 `2 2 2` 或 `3 3 3` 的 K 点。
*   **修正**：结构优化（尤其是 `cell-relax`）对 K 点密度非常敏感。K 点不足会导致晶格常数计算误差巨大。实际计算中，请务必进行 K 点收敛测试（例如对比 `4x4x4`, `6x6x6`, `8x8x8` 的结果）。