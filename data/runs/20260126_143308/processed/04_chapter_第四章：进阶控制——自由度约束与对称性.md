# 第四章：进阶控制——自由度约束与对称性

在材料模拟的实际战场中，"全自由度弛豫"（Full Relaxation）往往只是理想情况。更多时候，我们需要像外科手术一样精确控制体系的自由度：比如在计算形变势时固定晶胞体积、在模拟异质结时固定底层原子、或者在相变研究中强制保持某种对称性。

本章将深入 ABACUS 的自由度控制机制，并为习惯了 VASP `ISIF` 参数的用户提供一份详尽的迁移指南。

---

## 4.1 限定自由度优化 (Constrained Optimization)

ABACUS 的结构优化逻辑与 VASP 略有不同。VASP 通过 `ISIF` 一个参数控制所有行为，而 ABACUS 将**原子弛豫**和**晶胞弛豫**拆分为不同的计算类型，并通过额外的参数施加约束。

### 4.1.1 核心策略：分步优化 (Step-by-Step)

在开始任何晶胞变动之前，我必须强调一条**黄金法则**：

> **⚠️ 专家建议：永远不要直接运行 `cell-relax`**
>
> 即使你的最终目标是优化晶胞常数，也请务必遵循 **"先原子，后晶胞"** 的策略：
> 1.  先设置 `calculation = relax`（固定晶胞，仅优化原子位置），让原子回归到当前晶格下的能量极小值。
> 2.  读取上一步的结构（`STRU`），再开启 `calculation = cell-relax`（同时优化晶胞和原子）。
>
> **原因**：如果初始结构离平衡态较远，直接变动晶胞会导致应力（Stress）计算极度不稳定，甚至导致基组平面波数量剧烈波动（PW 模式下）或积分格点错位，引发计算发散。

### 4.1.2 晶胞约束：`fixed_axes`

当你开启 `calculation = cell-relax` 时，默认行为是三维全自由度优化（体积和形状均可变）。通过 `fixed_axes` 参数，我们可以锁定特定的自由度。

**常用场景配置 (`INPUT` 文件)：**

1.  **固定体积，优化形状 (Volume Conserved)**
    *   *应用场景*：计算剪切模量、模拟不可压缩流体极限下的晶体形变。
    *   *参数设置*：
        ```bash
        calculation  cell-relax
        cal_stress   1           # 必须开启应力计算
        fixed_axes   volume      # 锁定体积，允许形状改变
        ```

2.  **固定形状，优化体积 (Shape Conserved)**
    *   *应用场景*：状态方程 (EOS) 拟合、各向同性材料的快速预优化。
    *   *参数设置*：
        ```bash
        calculation  cell-relax
        cal_stress   1
        fixed_axes   shape       # 锁定晶胞角度和轴长比例，仅允许整体缩放
        # 注意：在旧版本或某些特定算法中可能被称为 r_shape，请以当前文档为准
        ```

3.  **固定特定方向 (自定义约束)**
    *   *注意*：ABACUS 目前原生支持的 `fixed_axes` 主要是 `volume`, `shape`, `None` (全放开)。如果需要像 VASP `ISIF=4` 那样只固定 z 轴但放开 xy 平面（例如层状材料优化），建议使用 **ASE 接口**（详见 4.2 节）或通过脚本迭代实现。

### 4.1.3 原子约束：`fixed_atoms`

原子位置的约束**不**在 `INPUT` 文件中设置，而是直接写在 `STRU` 结构文件中。这对于表面吸附、缺陷计算或异质结模拟至关重要。

**`STRU` 文件示例：**

```text
ATOMIC_POSITIONS
Direct
# 元素名  x  y  z  m_x  m_y  m_z
Si
0.000000 0.000000 0.000000 0 0 0  <-- 最后的三个 0 代表固定 x,y,z
0.250000 0.250000 0.250000 1 1 1  <-- 1 代表允许移动
```

*   **格式说明**：在原子坐标后的三个整数分别对应 x, y, z 方向的移动许可（1=Move, 0=Fix）。
*   **注意**：如果你的 `STRU` 文件中没有这三列，ABACUS 默认视为 `1 1 1`（全自由）。

---

## 4.2 对称性保持与 ASE 接口

在优化高对称性体系（如立方钙钛矿）时，数值噪声有时会导致结构“坍塌”到低对称性（如四方相）。防止这种情况有两种途径。

### 4.2.1 ABACUS 原生对称性控制

在 `INPUT` 文件中：

```bash
symmetry    1      # 开启对称性分析（默认）
```

*   **原理**：ABACUS 会在每一步 SCF 和结构优化中，强制将电荷密度和受力投影回该空间群的对称操作中。这通常足以维持晶体结构不发生对称性破缺。
*   **风险**：如果初始结构本身有微小的对称性破缺（例如手动微扰过的结构），`symmetry = 1` 可能会强制将其“修复”回高对称性，或者因为无法找到高对称性而报错。

### 4.2.2 进阶武器：ASE (Atomic Simulation Environment)

对于复杂的约束（例如：只允许单斜晶系的 $\beta$ 角变化，固定 $\alpha, \gamma$），ABACUS 原生参数可能不够灵活。此时，推荐使用 Python 的 ASE 库配合 ABACUS 计算器。

**实战代码：使用 ASE 的 `UnitCellFilter` 进行强约束优化**

```python
from ase.io import read
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS
from abacus_interface import AbacusCalculator # 假设你安装了 abacus-develop 的 python 接口

# 1. 读取结构
atoms = read("STRU", format="abacus")

# 2. 设置计算器
calc = AbacusCalculator(
    directory='.',
    inp_params={
        'calculation': 'scf',  # 注意：这里设为 scf，因为弛豫由 ASE 接管
        'ecutwfc': 100,
        'scf_thr': 1e-6,
        'cal_stress': 1,       # 必须计算应力
        'cal_force': 1
    }
)
atoms.calc = calc

# 3. 设置过滤器 (Filter)
# scalar_pressure: 设置外压
# mask: 6个分量对应 [xx, yy, zz, yz, xz, xy]
# 例如：只允许 z 轴 (zz) 弛豫，固定其他
ucf = UnitCellFilter(atoms, scalar_pressure=0.0, mask=[0, 0, 1, 0, 0, 0])

# 4. 运行优化
opt = BFGS(ucf)
opt.run(fmax=0.05) # 0.05 eV/Ang
```

*   **优势**：ASE 的 `UnitCellFilter` 提供了极高自由度的掩码（Mask）机制，可以实现任意维度的晶胞约束。

---

## 4.3 VASP 用户迁移指南 (ISIF 对照表)

这是 VASP 用户最关心的部分。VASP 的 `ISIF` 一个参数涵盖了原子位置、晶胞形状、晶胞体积三个维度。在 ABACUS 中，我们需要组合 `calculation`、`fixed_axes` 和 `STRU` 设置来实现同等效果。

**表 4-1: VASP ISIF 与 ABACUS 参数映射表**

| VASP ISIF | 物理含义 | ABACUS 对应设置 (`INPUT`) | 备注 |
| :--- | :--- | :--- | :--- |
| **2** | **仅优化原子**<br>固定晶胞 | `calculation = relax`<br>`cal_force = 1` | 最常用的设置。 |
| **3** | **全优化**<br>原子+形状+体积 | `calculation = cell-relax`<br>`fixed_axes = None`<br>`cal_stress = 1` | 对应全自由度弛豫。 |
| **4** | **固定体积**<br>优化原子+形状 | `calculation = cell-relax`<br>`fixed_axes = volume`<br>`cal_stress = 1` | 常用于相变路径计算。 |
| **5** | **仅优化形状**<br>固定原子+体积 | ❌ *不推荐直接对应* | ABACUS 优化晶胞时通常也会移动原子。如需固定原子，需在 `STRU` 中设为 `0 0 0`。 |
| **7** | **仅优化体积**<br>固定原子+形状 | `calculation = cell-relax`<br>`fixed_axes = shape`<br>+ `STRU` 中固定原子 | 类似于状态方程计算。 |

**关键区别提示：**
*   **VASP ISIF=3**：同时变动晶胞和原子。
*   **ABACUS cell-relax**：也是同时变动晶胞和原子。
*   **VASP ISIF=2**：等同于 ABACUS 的 `relax`。

---

## 附录：常见陷阱与调试清单

### 1. LCAO 基组的精度陷阱 (Egg box effect)
在使用 LCAO（数值原子轨道）基组进行结构优化时，切勿盲目追求极高的受力收敛精度。
*   **现象**：原子在某个位置附近来回震荡，无法收敛到 `1e-3 eV/Ang`。
*   **原因**：由于积分网格的存在，原子移动时会感受到非物理的“蛋盒势”（Egg box effect）。
*   **建议**：
    *   LCAO 模式下，`force_thr_ev` 设置为 **0.02 ~ 0.05 eV/Ang** 通常已足够发表级精度。
    *   如果必须追求更高精度（如声子谱计算），请切换到平面波（PW）基组或显著增加 `ecutwfc`。

### 2. 单位换算
*   `stress_thr`：ABACUS 中应力收敛阈值的单位通常是 **kBar**。
    *   1 kBar = 0.1 GPa。
    *   不要把 1 GPa 的标准误设为 1 kBar。

### 3. "SCF 不收敛，优化无意义"
这是一个老生常谈但极易被忽视的问题。如果你的 `out_stru` 显示结构跑飞了，或者能量震荡：
*   **第一步**：检查 `OUT.*/running_scf.log`。
*   **判断**：如果电子步（SCF）没有收敛（达到 `scf_nmax`），那么计算出的力（Force）和应力（Stress）就是错误的垃圾数据。
*   **对策**：先调整 `mixing_beta`、`smearing` 或 `ks_solver` 让 SCF 收敛，再谈结构优化。