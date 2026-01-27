# 第三章：输入文件配置实战

在上一章中，我们已经准备好了计算所需的“原材料”：结构文件、赝势文件和轨道文件。但在将它们交给 ABACUS 内核处理之前，我们需要通过精确的配置将这些文件“链接”起来，并告诉软件以何种精度进行计算。

本章将深入 ABACUS 的核心配置文件 `STRU` 和 `INPUT`，从开发者的视角解读如何正确映射物理参数。我们将重点解决两个新手最容易犯错的问题：如何在 `STRU` 中正确关联双重文件，以及在 LCAO 模式下如何科学设置截断能。

---

## 3.1 STRU 文件中的双重指定

在 ABACUS 的 LCAO（原子轨道线性组合）计算模式下，每个原子物种（Species）需要同时指定两类外部文件：
1.  **赝势文件 (`.upf` / `.vps` 等)**：描述原子核与内层电子对价电子的作用。
2.  **轨道文件 (`.orb`)**：描述价电子的数值原子轨道基组。

这两个文件的路径定义均位于 `STRU` 文件中，但分布在不同的模块块（Block）内。**严格的对应关系**是计算成功的关键。

### 3.1.1 ATOMIC_SPECIES：定义元素与赝势

`ATOMIC_SPECIES` 模块用于定义计算中涉及的元素种类、原子质量以及对应的赝势文件。

**基本语法：**
```text
ATOMIC_SPECIES
Element_Label   Atomic_Mass   Pseudopotential_File_Path
...
```

**实战示例 (STRU 文件片段)：**
```fortran
ATOMIC_SPECIES
Si  28.0855  Si_ONCV_PBE-1.0.upf
C   12.0107  C_ONCV_PBE-1.0.upf
```

**教授点评：**
*   **Element_Label**：这是你给元素起的“标签”。通常使用元素符号（如 `Si`），但在做反铁磁或特殊标记时，也可以写成 `Si1`、`Si2`，只要与后文 `ATOMIC_POSITIONS` 中的标签一致即可。
*   **路径问题**：示例中直接写了文件名，这意味着 ABACUS 会在当前运行目录下寻找该文件。如果文件在统一的库中，建议使用相对路径（如 `../../PP/Si.upf`）或绝对路径，避免拷贝文件带来的混乱。

### 3.1.2 NUMERICAL_ORBITAL：挂载数值轨道

这是 LCAO 模式特有的模块。在此模块中，你需要为 `ATOMIC_SPECIES` 中定义的每一种元素指定对应的轨道文件。

**关键规则：**
ABACUS 解析 `NUMERICAL_ORBITAL` 时，**不读取标签，仅依靠顺序**。这意味着：
> **`NUMERICAL_ORBITAL` 中列出文件的顺序，必须严格与 `ATOMIC_SPECIES` 中元素的出现顺序一一对应！**

**实战示例 (STRU 文件片段)：**
```fortran
NUMERICAL_ORBITAL
Si_TZVP-SR.orb
C_TZVP-SR.orb
```

**错误示范（将导致严重物理错误或程序崩溃）：**
如果在 `ATOMIC_SPECIES` 中先定义了 Si 后定义了 C，但在 `NUMERICAL_ORBITAL` 中先写了 C 的轨道文件，ABACUS 会尝试用 C 的轨道去描述 Si 原子，导致哈密顿量构建完全错误。

### 3.1.3 完整 STRU 文件结构图解

为了让你一目了然，以下是一个标准的 LCAO 模式 `STRU` 文件结构：

```fortran
ATOMIC_SPECIES
Si  28.0855  ./pp_orb/Si_ONCV_PBE-1.0.upf  // 第一种元素：硅
C   12.0107  ./pp_orb/C_ONCV_PBE-1.0.upf   // 第二种元素：碳

NUMERICAL_ORBITAL
./pp_orb/Si_TZVP-SR.orb                    // 对应第一种元素（Si）的轨道
./pp_orb/C_TZVP-SR.orb                     // 对应第二种元素（C）的轨道

LATTICE_CONSTANT
1.889726125  // 1.0 Angstrom in Bohr

LATTICE_VECTORS
5.43 0.00 0.00
0.00 5.43 0.00
0.00 0.00 5.43

ATOMIC_POSITIONS
Direct
Si  0.00 0.00 0.00 1 1 1
...
```

---

## 3.2 关键截断能参数设置

在 `INPUT` 文件中，有两个参数控制着计算的精度和格点密度：`ecutwfc` 和 `ecutrho`。

很多初学者有一个误区：“我使用的是 LCAO（原子轨道）基组，为什么还需要设置平面波截断能（Cutoff Energy）？”

**物理原理揭秘：**
虽然 LCAO 模式下波函数是用原子轨道展开的，但 ABACUS 在处理**哈密顿量矩阵元积分**（特别是涉及电子密度、哈特里势和交换关联势的部分）时，使用的是**实空间均匀格点（Real-space Uniform Grid）**。
这个格点的疏密程度，直接由“等效的平面波截断能”来控制。如果格点太稀疏，积分精度不够，会出现“蛋格效应”（Egg-box effect），导致力计算不准；如果格点太密，则白白浪费内存和计算时间。

### 3.2.1 如何设置 ecutwfc

`ecutwfc` (Energy CUToff for WaveFunCtion) 定义了波函数相关量的最大动能截断。在 LCAO 模式下，它决定了积分网格的基础密度。

**最佳实践：查阅轨道文件头**
ABACUS 官方发布的数值轨道文件（`.orb`）通常经过了严格测试，文件头中会包含生成该轨道时推荐的截断能。

**操作步骤：**
1.  在 Linux 终端使用 `head` 命令查看轨道文件前几行：
    ```bash
    head -n 20 Si_TZVP-SR.orb
    ```
2.  寻找类似 `Cutoff` 或 `Suggested cutoff` 的字段。
3.  **取最大值原则**：如果体系中有多种元素（如 Si 和 C），查看所有对应的 `.orb` 文件，取其中**最大**的推荐值作为 `INPUT` 中的 `ecutwfc`。

**经验参考值：**
如果文件中没有明确说明，对于常用的 DZP 或 TZVP 轨道，通常建议：
*   **标准精度**：50 - 60 Ry
*   **高精度**：80 - 100 Ry

### 3.2.2 如何设置 ecutrho

`ecutrho` (Energy CUToff for RHO/Charge Density) 定义了电荷密度和势场的截断能。

**物理约束：**
在密度泛函理论中，电荷密度 $\rho(r)$ 大致正比于波函数 $\psi(r)$ 的模方。在倒空间（平面波空间）中，这意味着描述密度所需的动量截止向量是波函数的 2 倍，对应的能量（与动量平方成正比）则是波函数的 4 倍。

**黄金法则：**
> **`ecutrho` 应设置为 `ecutwfc` 的 4 倍。**

虽然在某些极其特殊的软赝势情况下可以降低到 2.5-3 倍，但在 ABACUS 实战中，为了保证数值积分的稳定性，**始终推荐保持 4 倍关系**。

### 3.2.3 INPUT 文件配置实战

将上述逻辑落实到 `INPUT` 文件中，假设我们查阅轨道文件后决定使用 60 Ry 作为基准：

```bash
INPUT_PARAMETERS
# ... 其他参数 ...

# 电子结构计算精度控制
basis_type      lcao    # 明确指定使用 LCAO 基组
ecutwfc         60      # 单位：Ry，根据轨道文件推荐值设置
ecutrho         240     # 单位：Ry，严格遵守 4 * ecutwfc 规则

# ... 其他参数 ...
```

**开发者提示：**
*   **单位**：ABACUS 输入文件中的能量单位默认为 **Ry (Rydberg)**，而不是 eV。1 Ry $\approx$ 13.6057 eV。请务必注意单位换算。
*   **收敛性测试**：虽然有推荐值，但在进行全新的高精度研究（如声子谱、弹性常数计算）时，建议手动测试 `ecutwfc`（例如从 50 Ry 增加到 80 Ry），观察总能量和原子受力的收敛情况，以确定该体系的最佳参数。