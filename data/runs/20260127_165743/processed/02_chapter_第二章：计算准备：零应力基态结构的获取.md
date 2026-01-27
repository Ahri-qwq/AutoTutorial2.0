# 第二章：计算准备：零应力基态结构的获取

在开始计算弹性常数之前，我们必须先进行一项至关重要的“调音”工作——获取零应力（Zero-stress）的基态结构。

弹性常数（Elastic Constants）本质上描述的是晶格在平衡位置附近对微小形变的响应。根据广义胡克定律，它定义在应力为零的平衡态。如果你的初始结构中存在残余应力（Residual Stress），哪怕只有 1-2 kbar，都会导致后续施加应变后的能量-应变曲线发生偏移，计算出的弹性模量可能出现巨大误差，甚至出现物理上不合理的负值。

本章将指导你如何使用 ABACUS 获得完美的基态结构，并理解这一步与后续步骤的关键区别。

---

## 2.1 晶胞与原子位置的同时优化 (Variable-Cell Relaxation)

要消除残余应力，我们需要允许晶胞形状（Cell）和内部原子位置（Ions）同时调整。在 ABACUS 中，这一过程通过 `cell-relax` 实现。

### 2.1.1 核心输入参数详解

在 `INPUT` 文件中，以下参数是获取零应力结构的关键：

```bash
INPUT_PARAMETERS
# 1. 任务类型
calculation     cell-relax  # 关键：同时优化晶胞和原子位置
relax_nmax      100         # 最大离子步数，建议设置足够大

# 2. 应力与力的计算
cal_stress      1           # 绝对必须：开启应力张量计算
cal_force       1           # 开启原子力计算

# 3. 收敛标准 (建议使用严格标准)
force_thr_ev    0.001       # 力收敛标准 (eV/Ang)，建议比默认值更小
stress_thr      1.0         # 应力收敛标准 (kbar)，确保残余应力极低

# 4. 结果输出
out_stru        1           # 输出优化后的 STRU 文件
```

#### 💡 教授的锦囊妙计 (Expert Tips)

1.  **`cal_stress 1` 是灵魂**：很多新手在从 SCF 计算转为优化计算时，容易忘记开启此选项。如果 ABACUS 不计算应力张量，它就不知道该如何调整晶胞，导致优化失败或仅优化了原子位置。
2.  **收敛标准要严苛**：对于常规能带计算，默认的收敛标准可能够用了。但对于弹性常数计算，我们是在衡量能量的二阶导数，对噪声非常敏感。建议将 `stress_thr` 设为 1 kbar 或更低（默认通常较宽松）。
3.  **基组选择**：虽然 LCAO 基组效率极高，但在处理晶胞形状变化时，平面波（PW）基组通常能提供更平滑的应力曲线。如果你使用 LCAO (`basis_type lcao`)，请务必确保基组精度足够（如使用 DZP 或 TZP），并进行充分的截断能测试以减小基组重叠误差（BSSE）带来的虚假应力。

### 2.1.2 流程示意

1.  **准备初猜结构**：从数据库（如 Materials Project）下载 CIF 文件，转化为 ABACUS 的 `STRU` 文件。
2.  **执行 `cell-relax`**：运行 ABACUS。
3.  **更新结构**：计算结束后，在输出目录（如 `OUT.suffix/`）中找到 `STRU_ION_D`（最终结构），将其复制并覆盖原有的 `STRU` 文件，作为后续计算的起点。

---

## 2.2 关键概念区分：Step 1 与 Step 2 的天壤之别

这是本教程中最容易出错的地方，请务必仔细阅读。

计算弹性常数分为两个截然不同的阶段，它们的物理目的完全不同，使用的计算类型也**绝对不能混淆**。

| 阶段 | 目的 | ABACUS `calculation` | 晶胞状态 | 原子状态 |
| :--- | :--- | :--- | :--- | :--- |
| **Step 1 (本章)** | **消除残余应力** | **`cell-relax`** | **可变 (Variable)** | 可变 |
| **Step 2 (下一章)** | **计算应力-应变响应** | **`relax`** (或 `scf`) | **固定 (Fixed)** | 可变 (或固定) |

> **⚠️ 红色警报 (Critical Warning)**：
>
> 在进入 **Step 2（施加应变）** 后，**切勿** 再次使用 `cell-relax`！
>
> 如果你在 Step 2 对一个已经施加了 +1% 应变的结构使用 `cell-relax`，ABACUS 会认为这个应变是“不舒服”的，并自动把晶胞优化回零应力的原始状态。这样你辛苦施加的应变就被“优化”掉了，计算出的弹性常数将完全错误。
>
> **口诀：先变胞 (`cell-relax`) 找基态，后定胞 (`relax`) 算应变。**

---

## 2.3 深入理解：原子弛豫 (Internal Relaxation)

在 Step 2（施加应变后的计算）中，我们通常使用 `calculation = relax`（固定晶胞，优化原子），这对应于物理上的**松弛离子弹性常数 (Relaxed-ion elastic constants)**。

### 为什么需要 Internal Relaxation？
当你将晶胞拉伸（施加应变 $\varepsilon$）时，晶格内的原子受力不再平衡。
- **Clamped-ion (夹持离子)**: 如果我们强行固定原子分数坐标不动（使用 `calculation = scf`），此时测得的应力对应于高频极限下的弹性响应（原子来不及调整）。
- **Relaxed-ion (松弛离子)**: 实际上，原子会发生微小的位移以释放内部的赫尔曼-费曼力（Hellmann-Feynman forces），达到新的亚稳态。这是我们在宏观实验中测量到的弹性常数。

**ABACUS 实战指南**：
- **默认推荐**：绝大多数情况下，你应该计算 Relaxed-ion 结果。因此，在后续使用 `abacustest` 生成变形任务时，默认会自动设置 `calculation = relax`。
- **特殊需求**：如果你确实需要计算 Clamped-ion 弹性常数（例如为了理论分析），需在 `abacustest` 命令中显式添加 `--norelax` 参数，这将强制任务类型为 `scf`。

---

## 2.4 晶体取向与标准化 (Standardization)

弹性张量是各向异性的，其数值高度依赖于坐标系的选取。

### 2.4.1 为什么要注意取向？
对于硅（Si）这样的立方晶系，$C_{11}$ 通常指沿着 $<100>$ 方向拉伸的模量。如果你的 `STRU` 文件中，晶格矢量 $a$ 并没有沿着笛卡尔坐标系的 $x$ 轴，而是沿着 $<110>$ 方向，那么你算出来的矩阵元数值将与标准数据库无法直接对比。

### 2.4.2 IEEE 标准与 Materials Project
Materials Project 数据库通常遵循 **IEEE 176-1987** 标准来定义晶体取向。
- **建议**：对于立方（Cubic）、四方（Tetragonal）、正交（Orthorhombic）晶系，尽量使用**惯用胞（Conventional Cell）**，并确保晶格矢量与笛卡尔坐标轴平行（$a \parallel x, b \parallel y, c \parallel z$）。
- **操作**：可以使用 ASE 或 `abacustest` 的辅助工具在准备结构阶段进行标准化旋转。

---

## 2.5 实战操作：使用 `abacustest` 准备基态优化

`abacustest` 是 ABACUS 的官方辅助工具，可以极大简化工作流。以下是针对硅（Si）进行零应力优化的标准命令。

### 1. 准备初始结构
假设你已经有了 `Si.cif` 或 `Si.stru`。

### 2. 生成优化任务
使用以下命令生成用于 `cell-relax` 的输入文件：

```bash
abacustest prepare -f Si.cif -p input_param.json -s 
```
*(注：此处假设你通过 json 或命令行参数指定了 `calculation: cell-relax`，或者直接手动修改生成的 INPUT 文件)*

更直接的命令行方式（参考核心案例）：

```bash
# 针对 Si 结构生成 cell-relax 任务
python3 -m abacustest.lib_prepare.main \
    -f Si.cif \
    --ftype cif \
    --jtype cell-relax \
    --lcao \
    --folder-syntax "Si_relax"
```

*   `--jtype cell-relax`: 自动设置 `calculation = cell-relax`。
*   `--lcao`: 自动配置 LCAO 基组相关参数（如 `basis_type`, `mixing_type` 等）。

### 3. 提交计算与更新结构
计算完成后，你需要提取优化后的结构用于下一步。

```bash
# 假设优化任务在 Si_relax/ 目录下运行
# 将优化后的结构 (STRU_ION_D) 复制出来，命名为 Si_opt.stru
cp Si_relax/OUT.*/STRU_ION_D Si_opt.stru
```

> **⚠️ 风险提示 (Risk Warning)**
>
> `abacustest` 的 `prepare` 命令非常强大，但它通常会**覆盖**目标文件夹中已有的文件。
> **切勿在计算完成后，再次在同一目录下重复执行输入文件准备命令！**
> 否则，你辛苦算出来的 `OUT` 文件夹和数据可能会被清空或覆盖。建议在不同阶段使用不同的文件夹名称（如 `step1_relax`, `step2_elastic`）。

---

**下一章预告**：
手握完美的零应力结构 `Si_opt.stru`，我们将在第三章正式进入核心环节——施加应变并提取应力张量。我们将详细讲解如何使用 `abacustest` 自动化生成那 24 个繁琐的变形结构。