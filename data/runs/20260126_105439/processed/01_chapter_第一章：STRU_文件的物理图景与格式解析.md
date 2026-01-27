# 第一章：STRU 文件的物理图景与格式解析

在第一性原理计算中，如果说 `INPUT` 文件是控制计算流程的“大脑”，那么 `STRU` 文件就是计算发生的“舞台”。它不仅定义了原子在哪里（几何结构），还定义了原子是谁（元素与赝势）、原子能做什么（基组与自由度）。

本章将带你像物理学家一样审视 `STRU` 文件，透过文本格式看到背后的物理图景。

---

## ⚠️ 核心预警：工欲善其事 (关于 ASE 的安装)

在深入文件细节之前，必须纠正一个新手最常犯的错误。许多教程会推荐使用 ASE (Atomic Simulation Environment) 来生成结构文件。

> **CRITICAL WARNING**
>
> **请勿直接使用 `pip install ase` 来获取 ABACUS 支持！**
>
> 官方的 ASE 库（PyPI 源）尚未完全包含 ABACUS 的最新接口。如果你直接 `pip install ase`，在导出 STRU 文件时可能会遇到格式错误或功能缺失。
>
> **正确做法**是安装专用的 `ase-abacus` 分支：
> ```bash
> # 1. 下载适配 ABACUS 的 ASE 接口仓库
> git clone https://gitlab.com/1041176461/ase-abacus.git
> 
> # 2. 进入目录并安装
> cd ase-abacus
> pip install .
> ```
> 只有这样，你才能确保 `write('STRU', atoms, format='abacus')` 生成的文件是完全合规的。

---

## 1.1 模拟单元的定义（晶格与基组）

一个量子力学模拟单元由两部分组成：**边界条件**（晶格）和**希尔伯特空间基底**（元素与轨道）。

### 1.1.1 元素与基组 (`ATOMIC_SPECIES` & `NUMERICAL_ORBITAL`)

这部分告诉 ABACUS 体系中包含哪些元素，以及用什么数学函数来描述电子。

```abacus
ATOMIC_SPECIES
Al  26.98154  Al_ONCV_PBE-1.0.upf
O   15.99940  O_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb
O_gga_7au_100Ry_2s2p1d.orb
```

*   **ATOMIC_SPECIES**:
    *   第一列：元素符号（Label）。
    *   第二列：原子质量（Mass），单位为 a.u. (amu)。
    *   第三列：赝势文件名。**注意**：ABACUS 会根据 `INPUT` 文件中 `pseudo_dir` 指定的路径去寻找这些文件。
*   **NUMERICAL_ORBITAL**:
    *   **仅在 LCAO 模式下需要**（即 `INPUT` 中 `basis_type lcao`）。
    *   列出数值原子轨道文件（.orb）。顺序必须与 `ATOMIC_SPECIES` 中的元素顺序严格对应。

> **💡 实战技巧：告别手动文件名匹配**
>
> 手动填写 `upf` 和 `orb` 文件名极易出错。推荐使用 **`abacustest`** 工具。
> 你只需在赝势库目录下准备一个 `element.json` 文件，定义元素到文件名的映射：
> ```json
> {
>     "Al": "Al_ONCV_PBE-1.0.upf",
>     "O": "O_ONCV_PBE-1.0.upf"
> }
> ```
> `abacustest` 可以自动读取该配置并生成正确的 STRU 文件，彻底避免文件名拼写错误。

### 1.1.2 晶格常数与单位制 (`LATTICE_CONSTANT` & `LATTICE_VECTORS`)

这是新手最容易困惑的地方。ABACUS 的底层计算基于**原子单位制 (Atomic Units)**，长度单位是 **Bohr (a.u.)**，而不是我们熟悉的 Angstrom (Å)。

```abacus
LATTICE_CONSTANT
1.8897261258

LATTICE_VECTORS
10.0  0.0  0.0
0.0  10.0  0.0
0.0   0.0 10.0
```

*   **LATTICE_CONSTANT (缩放因子)**:
    *   物理意义：这是一个**单位转换系数**。
    *   **为什么是 1.8897?**
        *   $1 \text{ Bohr} \approx 0.529177 \text{ \AA}$
        *   $1 \text{ \AA} \approx 1.889726125 \text{ Bohr}$
    *   **工作原理**：ABACUS 读取 `LATTICE_VECTORS` 中的数值后，会将它们**乘以** `LATTICE_CONSTANT` 得到最终的 Bohr 单位长度。
    *   **实战建议**：通常我们将 `LATTICE_CONSTANT` 设为 `1.8897261258`，这样 `LATTICE_VECTORS` 中的数值就可以直接填写 **Angstrom** 单位的数值。

*   **LATTICE_VECTORS**:
    *   定义晶胞的三个基矢量 $\vec{a}, \vec{b}, \vec{c}$。
    *   格式为 3x3 矩阵，每一行代表一个矢量分量 ($x, y, z$)。

---

## 1.2 原子构型的描述（位置与自由度）

定义好舞台后，我们需要安排演员（原子）的位置、动作（弛豫）和状态（磁性）。

### 1.2.1 坐标类型与数据块 (`ATOMIC_POSITIONS`)

```abacus
ATOMIC_POSITIONS
Direct  // 或 Cartesian

Al      // 元素标签，必须与 ATOMIC_SPECIES 对应
0.0     // 该元素的初始磁矩（低优先级，通常被覆盖，设为0即可）
2       // 该元素的原子数量

// 下面是核心数据行：坐标 + 移动限制 + 磁矩
// x          y          z          mx my mz  mag  moment
0.00000000 0.00000000 0.00000000 0  0  0   mag  0.0
0.50000000 0.50000000 0.50000000 1  1  1   mag  2.0
```

### 1.2.2 详解核心数据行

让我们放大看一行典型的原子定义：

`0.500000 0.500000 0.500000 1 1 1 mag 2.0`

这一行包含了三个物理维度的信息：

1.  **几何坐标 (前三列)**:
    *   如果类型是 **`Direct`** (推荐)：表示**分数坐标**。原子位置 $\vec{r} = x\vec{a} + y\vec{b} + z\vec{c}$。数值范围通常在 0 到 1 之间。
    *   如果类型是 **`Cartesian`**: 表示**笛卡尔坐标**。注意，这里的数值也会被 `LATTICE_CONSTANT` 缩放！
        *   *注：存在 `Cartesian_angstrom` 选项可跳过缩放，但兼容性较差，不建议初学者使用。*

2.  **动力学自由度 (中间三列 `mx my mz`)**:
    *   **`0`**: **固定 (Fixed)**。在结构弛豫（Relaxation）或分子动力学（MD）中，该原子在对应方向上不允许移动。
    *   **`1`**: **弛豫 (Relax)**。允许原子在对应方向上受力移动。
    *   *应用场景*：在计算表面吸附时，通常将底部的原子层设为 `0 0 0` 以模拟体相环境，而将表面和吸附原子设为 `1 1 1`。

3.  **磁学性质 (尾部 `mag moment`)**:
    *   **`mag`**: 关键词，标识后面跟随的是磁矩设置。
    *   **`2.0`**: **初始磁矩 (Initial Magnetic Moment)**，单位为 $\mu_B$。
    *   *物理意义*：这只是**初始猜测值**。在自洽计算（SCF）开始时，ABACUS 会根据这个值构建初始自旋密度。最终收敛后的磁矩由电子结构决定，可能与此值不同。
    *   *注意*：要使此设置生效，必须在 `INPUT` 文件中开启自旋极化计算（`nspin 2`）。

### 1.2.3 完整 STRU 文件示例

结合上述所有知识，以下是一个完整的 FCC 铝（Aluminum）的 STRU 文件示例，请仔细阅读注释：

```abacus
ATOMIC_SPECIES
Al 26.9815385 Al_ONCV_PBE-1.0.upf  // 元素 质量 赝势

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb        // 轨道文件(LCAO模式)

LATTICE_CONSTANT
1.8897261258369282                 // 将 Angstrom 转换为 Bohr 的系数

LATTICE_VECTORS
0.0000000000 2.0250000000 2.0250000000 // 基矢量 a (单位: Angstrom)
2.0250000000 0.0000000000 2.0250000000 // 基矢量 b
2.0250000000 2.0250000000 0.0000000000 // 基矢量 c

ATOMIC_POSITIONS
Direct                             // 使用分数坐标

Al                                 // 元素 Al
0.0                                // 默认磁矩 (占位符)
1                                  // 原子总数
// x    y    z    mx my mz  mag setting
0.0  0.0  0.0   1  1  1   mag 0.0  // 原点原子，允许移动，无初始磁矩
```

通过本章的学习，你应该已经能够“读懂” STRU 文件中的每一个数字。在下一章中，我们将结合 `INPUT` 文件，让这些静止的原子“动”起来，进行第一次电子结构计算。