# 第五章：进阶调控：占据矩阵控制 (OMC)

在强关联体系（如过渡金属氧化物、镧系/锕系元素化合物）的计算中，我们经常面临一个棘手的问题：**亚稳态困境**。由于 $d$ 或 $f$ 电子的高度局域性，体系往往存在多个能量极小值（亚稳态）。标准的 DFT+U 自洽迭代有时会陷入错误的轨道占据状态，导致磁矩、能带结构甚至晶格畸变与实验严重不符。

ABACUS 提供了强大的 **占据矩阵控制 (Occupation Matrix Control, OMC)** 功能。它允许我们显式地读取、输出并锁定轨道占据矩阵，从而“引导”或“强制”体系进入我们预期的电子基态。

> **⚠️ 核心前置条件**：
> OMC 功能目前**仅支持 LCAO 基组** (`basis_type lcao`)。如果你使用的是平面波基组 (`basis_type pw`)，请勿开启此功能。

---

## 5.1 OMC 的工作原理与模式

在 DFT+U 理论中，Hubbard $U$ 能量项依赖于局域轨道的占据矩阵 $n_{mm'}^\sigma$。OMC 功能通过控制这个矩阵在自洽迭代（SCF）中的行为来实现对电子态的调控。

### 核心参数详解

在 `INPUT` 文件中，通过以下参数控制 OMC：

1.  **`omc` (整数)**: 控制占据矩阵的读写模式。
    *   **`0` (默认)**: **标准模式**。程序在每一步 SCF 中根据当前的波函数实时计算并更新占据矩阵。这是常规 DFT+U 计算的模式。
    *   **`1` (引导模式)**: **读取并更新**。程序在 SCF 的第一步读取外部提供的 `initial_onsite.dm` 文件作为初始占据矩阵，但在随后的 SCF 步中，程序会根据迭代情况自由更新矩阵。这用于给体系一个正确的“初始推力”。
    *   **`2` (锁定模式)**: **读取并固定**。程序读取 `initial_onsite.dm`，并在整个 SCF 过程中**强制保持**该矩阵不变。这用于复现特定轨道序，或在能带计算（Non-SCF）中保持基态电子分布。

2.  **`out_chg` (整数)**:
    *   设置为 `1` 时，ABACUS 会在计算结束时将最终的占据矩阵输出到 `OUT.suffix/onsite.dm` 文件中。这是进行 OMC 续算的必要步骤。

---

## 5.2 实战：固定占据矩阵流程 (以反铁磁 NiO 为例)

本节将演示一套完整的 OMC 工作流：首先通过标准计算获得正确的反铁磁基态，然后利用 OMC 锁定该状态。

### 案例背景：反铁磁 NiO
NiO 是典型的强关联反铁磁绝缘体。为了正确描述反铁磁序（AFM），我们需要将晶胞中的两个 Ni 原子区分对待：一个自旋向上，一个自旋向下。

### 第一步：准备 STRU 文件 (关键技巧)

**初学者最易犯的错误**是在 `STRU` 中只定义一种 Ni 元素。为了实现反铁磁，必须在 `ATOMIC_SPECIES` 中将 Ni 拆分为 `Ni1` 和 `Ni2`，以便分别设置磁矩和 U 值参数。

**`STRU` 文件片段：**
```text
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf  // 第一种 Ni，用于自旋向上
Ni2 58.693 Ni_ONCV_PBE-1.0.upf  // 第二种 Ni，用于自旋向下
O   15.999 O_ONCV_PBE-1.0.upf

ATOMIC_POSITIONS
Direct
...
Ni1 0.0 0.0 0.0 1.0  // magmom = 1.0 (正磁矩)
Ni2 0.5 0.5 0.5 -1.0 // magmom = -1.0 (负磁矩)
O   0.25 0.25 0.25 0.0
...
```

### 第二步：运行标准 DFT+U 计算 (`omc=0`)

首先进行一次标准的自洽计算，生成基态的占据矩阵。

**`INPUT` 文件设置：**
```text
INPUT_PARAMETERS
# ... 基础参数 ...
basis_type      lcao
calculation     scf

# DFT+U 设置
dft_plus_u      1
# 对应 STRU 中元素的顺序：Ni1(d), Ni2(d), O(no U)
orbital_corr    2 2 -1 
hubbard_u       5.0 5.0 0.0 

# OMC 设置
omc             0       # 标准计算
out_chg         1       # 必须开启！否则不输出 onsite.dm
```

> **⚠️ 数组对应关系警告**：
> `orbital_corr` 和 `hubbard_u` 的数值个数必须严格等于 `STRU` 中 `ATOMIC_SPECIES` 的行数。
> *   本例有 3 种原子 (Ni1, Ni2, O)，所以必须写 3 个数。
> *   `2` 代表对 d 轨道加 U，`-1` 代表不加 U。
> *   顺序必须是：Ni1 -> Ni2 -> O。

**运行并验证结果：**
运行完成后，检查 `OUT.NiO/running_scf.log`：
1.  **检查 U 值读取**：搜索 `L(S)DA+U` 区块，确认程序正确识别了两个 Ni 原子的 d 轨道 (L=2)。
2.  **检查磁性**：搜索 `absolute magnetism`。
    *   `total magnetism` 应接近 0 (反铁磁抵消)。
    *   `absolute magnetism` 应显著大于 0 (例如约 3.3 $\mu_B$/cell)，证明局部磁矩存在且反向排列。

此时，输出目录中应生成了 `onsite.dm` 文件。

### 第三步：使用 OMC 锁定状态 (`omc=2`)

假设我们需要基于上一步的电子态进行后续计算（如能带结构），或者想测试“锁定”功能。

1.  **文件操作**：
    将上一步生成的 `onsite.dm` 复制到工作目录，并**重命名**为 `initial_onsite.dm`。ABACUS 仅识别此文件名的输入。
    ```bash
    cp OUT.NiO/onsite.dm ./initial_onsite.dm
    ```

2.  **修改 `INPUT`**：
    ```text
    # 修改以下参数
    omc             2       # 开启锁定模式
    # out_chg       0       # 此时可以关闭输出，除非你想再次检查
    ```

3.  **再次运行并验证**：
    运行计算。查看 `running_scf.log`，你会发现：
    *   程序提示读取了 `initial_onsite.dm`。
    *   在每一以 `L(S)DA+U` 开头的输出块中，占据矩阵的数据应与第一步计算的最终结果完全一致，且在整个 SCF 过程中保持不变。
    *   最终的总能量 (`FINAL_ETOT`) 和磁矩应与第一步完全吻合。

---

## 附录：常见问题与进阶建议

### 1. 常见报错与排查

*   **Error: `hubbard_u` size is not equal to `ntype`**
    *   **原因**: `INPUT` 中的 `hubbard_u` 或 `orbital_corr` 数组长度与 `STRU` 中的原子类型数量不一致。
    *   **解决**: 即使某些原子不加 U（如氧原子），也必须在数组对应位置填 `0.0` (U值) 和 `-1` (轨道)。

*   **Error: OMC only supports LCAO basis**
    *   **原因**: `INPUT` 中设置了 `basis_type pw` 却开启了 `omc > 0`。
    *   **解决**: 切换为 `basis_type lcao`。

### 2. 收敛困难处理技巧

在强关联体系中，DFT+U 有时比标准 DFT 更难收敛。如果开启 OMC 前无法收敛，可尝试调整混合参数：
*   **`mixing_type`**: 尝试使用 `pulay` (默认) 或 `broyden`。
*   **`mixing_beta`**: 降低该值（如从 0.4 降至 0.1），减缓电荷密度更新幅度。
*   **`mixing_restart`**: 如果电荷密度震荡，可设为 `0.0` 禁用预测。
*   **`mixing_dmr`**: 专门针对 LCAO 的密度矩阵混合参数，对于 d/f 电子体系，适当降低此值有助于稳定占据矩阵。

### 3. 拓展阅读：Yukawa 势方法

除了手动指定 U 值，ABACUS 还支持基于 Yukawa 势的屏蔽方法来估算相互作用。
*   **参数**: `yukawa_potential 1`
*   **原理**: 不直接设置 `hubbard_u`，而是通过 `yukawa_lambda` (屏蔽长度) 来定义相互作用衰减。这为确定 U 值提供了一种非经验的物理途径。详情请参阅官方文档。