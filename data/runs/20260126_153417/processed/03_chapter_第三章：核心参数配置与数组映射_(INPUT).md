# 第三章：核心参数配置与数组映射 (INPUT)

在上一章中，我们已经完成了 `STRU` 文件的构建，特别是针对反铁磁氧化镍（NiO）体系，我们将镍原子人为拆分成了 `Ni1`（自旋向上）和 `Ni2`（自旋向下）两种原子类型（Atomic Species）。

本章将进入 ABACUS 计算的核心控制室——`INPUT` 文件。对于 DFT+U 计算而言，最关键的逻辑在于**建立 INPUT 参数数组与 STRU 原子类型之间严格的索引映射关系**。如果这种映射错位，程序将无法正确施加 Hubbard U 修正，导致物理结果完全错误。

> **⚠️ 核心前置条件**：
> 目前 ABACUS 的 DFT+U 功能**仅支持 LCAO 基组**。请务必在 `INPUT` 文件中确认设置了 `basis_type lcao`。

---

## 3.1 激活 DFT+U 计算

在 `INPUT` 文件中，开启 +U 功能需要一个总开关。

### 关键参数：`dft_plus_u`
*   **类型**: Boolean / Integer
*   **默认值**: 0 (False)
*   **说明**:
    *   `0`：关闭功能。即使你设置了后面提到的 U 值参数，程序也会完全忽略它们，执行标准的 DFT (LDA/GGA) 计算。
    *   `1`：激活 DFT+U 模块。

**示例代码**：
```bash
INPUT_PARAMETERS
# ... 其他常规参数 ...
basis_type      lcao    # 必须是 lcao
dft_plus_u      1       # 开启 +U 总开关
```

---

## 3.2 轨道与 U 值的数组映射 (核心逻辑)

这是初学者最容易出错的地方。ABACUS 不通过元素符号（如 "Ni"）来识别加 U 的对象，而是通过**原子类型的索引顺序**。

### 映射规则
在 `STRU` 文件中 `ATOMIC_SPECIES` 块定义的原子类型数量记为 `ntype`。在 `INPUT` 文件中，`orbital_corr` 和 `hubbard_u` 这两个参数必须是**长度为 `ntype` 的数组**，且顺序必须与 `STRU` 中的定义严格一一对应。

### 案例演示：反铁磁 NiO
假设你的 `STRU` 文件定义如下（注意顺序）：

**STRU 文件片段**：
```bash
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf  # 索引 1 (程序内部通常从0计数，这里指第一类)
Ni2 58.693 Ni_ONCV_PBE-1.0.upf  # 索引 2
O   15.999 O_ONCV_PBE-1.0.upf   # 索引 3
```
这里 `ntype = 3`。我们需要对 `Ni1` 和 `Ni2` 的 d 轨道加 U，对 `O` 不加 U。

### 关键参数详解

#### 1. `orbital_corr` (轨道校正类型)
*   **描述**: 指定每种原子类型需要进行 +U 修正的角动量量子数 $l$。
*   **取值**:
    *   `-1`: 不进行修正（用于非关联原子，如氧、碳等）。
    *   `2`: 对 d 轨道进行修正（过渡金属常见设置）。
    *   `3`: 对 f 轨道进行修正（镧系/锕系元素）。
*   **数组设置**:
    *   对应 `Ni1`: d 轨道 -> `2`
    *   对应 `Ni2`: d 轨道 -> `2`
    *   对应 `O`: 无 -> `-1`
    *   **结果**: `2 2 -1`

#### 2. `hubbard_u` (U 值大小)
*   **描述**: 指定对应轨道的有效 U 值（Hubbard U effective），单位为 **eV**。
*   **注意**: 即使对应的 `orbital_corr` 为 -1，这里也必须填一个占位数值（通常填 0.0），以保持数组长度一致。
*   **数组设置**:
    *   对应 `Ni1`: 5.0 eV
    *   对应 `Ni2`: 5.0 eV
    *   对应 `O`: 0.0 eV
    *   **结果**: `5.0 5.0 0.0`

### 完整的 INPUT 配置块
将上述逻辑整合，你的 `INPUT` 文件中关于 DFT+U 的部分应如下所示：

```bash
#Parameter DFT+U
dft_plus_u      1             # 激活开关
orbital_corr    2 2 -1        # 对应 STRU 中的: Ni1(d), Ni2(d), O(none)
hubbard_u       5.0 5.0 0.0   # 对应 U 值: 5.0eV, 5.0eV, 0.0eV
```

> **专家提示 (Ni1 与 Ni2 的拆分意义)**：
> 你可能会问，为什么不直接写一个 Ni？
> 在反铁磁计算中，我们需要在 `STRU` 的 `ATOM_FILES` 块中为不同原子指定相反的初始磁矩（如 `mag 2.0` 和 `mag -2.0`）。虽然 ABACUS 允许对同一种类型的不同原子设置不同磁矩，但为了逻辑清晰以及后续可能针对不同自旋态微调 U 值（虽然本例中 U 值相同），将它们定义为不同的 `Species` (Ni1, Ni2) 是最稳健的做法。这也强制要求我们在 `INPUT` 中显式地为它们分别设置参数。

---

## 3.3 进阶工作流：占据矩阵控制 (OMC)

在强关联体系中，DFT+U 计算容易陷入亚稳态（Local Minima）。为了锁定正确的电子基态（如特定的轨道序或磁序），ABACUS 提供了 `omc` (Occupation Matrix Control) 参数。

### 推荐工作流

#### 第一步：标准自洽计算 (`omc 0`)
首先进行标准的 DFT+U 计算，让电子密度自由弛豫。
*   **参数**: `omc 0` (默认值)
*   **目的**: 获得初步的收敛态，检查磁矩方向是否正确。

#### 第二步：锁定占据矩阵 (`omc 2`)
如果第一步得到了正确的磁序（例如反铁磁态），但你担心在后续做结构弛豫或能带计算时电子态发生“跳变”，可以锁定当前的占据矩阵。

1.  **准备文件**: 第一步计算结束后，在输出目录（如 `OUT.NiO/`）中找到 `onsite.dm` 文件。将其复制到工作目录并重命名为 `initial_onsite.dm`。
    ```bash
    cp OUT.NiO/onsite.dm ./initial_onsite.dm
    ```
2.  **修改 INPUT**:
    ```bash
    omc 2  # 读入 initial_onsite.dm 并在整个计算过程中保持不变
    ```
3.  **重新运行**: 程序将强制使用你提供的占据矩阵，这对于复杂氧化物计算非常有效。

---

## 3.4 结果验证：你真的加上 U 了吗？

计算开始后，不要只看最后能量，必须检查日志文件 `running_scf.log` 确认参数是否被正确读入。

### 1. 检查 L(S)DA+U 区块
打开 `running_scf.log`，搜索 `L(S)DA+U`。你应该看到类似以下的输出：

```text
 ---------------------------------------------------------
 L(S)DA+U:
 atom_type=0  L=2  chi=0    U=5ev   <-- 对应 Ni1, d轨道, 5eV
 atom_type=1  L=2  chi=0    U=5ev   <-- 对应 Ni2, d轨道, 5eV
 ---------------------------------------------------------
```
*验证点*：确认 `atom_type` 的数量和 `U` 值与你的设想一致。如果没有这个区块，说明 `dft_plus_u` 可能没打开，或者 `basis_type` 不是 `lcao`。

### 2. 检查磁矩 (Magnetism)
对于反铁磁 NiO，总磁矩应该接近 0，但局部磁矩应该很大。
搜索 `absolute magnetism`：

```text
      total magnetism (Bohr mag/cell) = 0.00000000  <-- 反铁磁特征：总磁矩为0
   absolute magnetism (Bohr mag/cell) = 3.35321634  <-- 局域磁矩之和不为0
```

同时检查 Mulliken 分析（需设置 `out_mul 1`）：
```text
Total Magnetism on atom  Ni1           1.82...
Total Magnetism on atom  Ni2          -1.82...  <-- 确认符号相反
```

---

### 本章小结
1.  **开关**: `dft_plus_u 1` 且必须是 LCAO 基组。
2.  **映射**: `orbital_corr` 和 `hubbard_u` 的数组长度必须等于 `ntype`，顺序严格对应 `STRU`。
3.  **技巧**: 利用 `omc 2` 和 `initial_onsite.dm` 可以锁定复杂的电子基态。

下一章，我们将深入探讨如何通过 `running_scf.log` 和输出文件进行更详细的电子结构分析。