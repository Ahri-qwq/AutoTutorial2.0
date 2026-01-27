# 第三章：核心实战：基于 abacustest 的自动化计算流程 (以 Si 为例)

在上一章中，我们熟悉了 ABACUS 的基本输入输出。本章我们将进入真正的“实战”环节。计算材料的弹性常数（Elastic Constants）是表征材料力学性质的基础，但手动施加应变、提交任务、提取应力张量并进行拟合的过程极其繁琐且容易出错。

我们将使用 ABACUS 官方配套的 Python 辅助工具库 `abacustest`，以最经典的立方晶系硅（Si）为例，演示如何通过一条命令完成从输入文件生成到结果分析的完整闭环。

---

## 3.1 自动化工作流环境配置与输入生成

计算弹性常数的核心思想是**应力-应变法（Stress-Strain Method）**：对晶胞施加一系列微小的应变（Strain, $\epsilon$），计算对应的应力（Stress, $\sigma$），利用广义胡克定律 $\sigma = C \epsilon$ 拟合得到弹性刚度矩阵 $C_{ij}$。

### 第一步：获取零应力基态结构 (Pre-requisite)

在施加应变之前，必须确保初始结构处于**零应力状态**。如果初始结构存在残余应力，会导致拟合结果出现巨大误差。

1.  **准备初始结构**：生成 Si 的标准惯用胞（Conventional Cell）。
2.  **执行结构优化**：
    *   **关键参数**：`calculation` 必须设为 `cell-relax`（同时优化原子位置和晶胞形状）。
    *   **目标**：确保应力张量对角元接近 0 kBar。

> **教授提示**：对于立方晶系，建议将晶格矢量沿笛卡尔坐标轴放置（即 $a$ 沿 $x$, $b$ 沿 $y$, $c$ 沿 $z$），这符合 IEEE 标准，能直接得到标准的 $C_{11}, C_{12}, C_{44}$ 而无需坐标变换。

### 第二步：生成变形结构输入文件

假设你已经完成了第一步，得到了优化后的结构文件 `STRU_relaxed`（请重命名为 `STRU`）以及对应的 `INPUT`（作为模板）和赝势/轨道文件。

现在，我们使用 `abacustest` 自动生成施加了不同应变的输入文件。在包含 `STRU`, `INPUT`, `*.upf`, `*.orb` 的目录下执行：

```bash
abacustest lib -m elastic --norm 0.01 --shear 0.01
```

**参数详解**：
*   `-m elastic`: 指定工作模式为弹性常数计算。
*   `--norm 0.01`: 设置最大**正应变**（Normal Strain）幅度为 1%。默认会生成 ±0.5% 和 ±1.0% 的变形点。
*   `--shear 0.01`: 设置最大**剪切应变**（Shear Strain）幅度为 1%。

**执行后的目录结构**：
执行成功后，当前目录下会自动生成一个名为 `elastic_result`（或类似命名，取决于版本）的文件夹，内部结构如下：

```text
.
├── run_elastic/
│   ├── deform-000/   (原始结构)
│   ├── deform-001/   (施加了第一种应变模式的结构)
│   ├── deform-002/   (施加了第二种应变模式的结构)
│   ├── ...
│   └── deform-xxx/
```

> **⚠️ 风险警示**：
> **不要在计算完成后重复执行上述输入文件准备命令！**
> `abacustest` 在生成文件时会清空目标文件夹。如果你已经跑完了计算，重复执行此命令会**直接删除所有计算结果**，导致前功尽弃。

---

## 3.2 提交计算与应力张量提取

这一步是物理意义最关键的环节。我们需要批量提交 `deform-xxx` 文件夹中的任务。

### 核心物理配置：Relaxed-ion vs Clamped-ion

在 `deform-xxx` 文件夹中，`abacustest` 自动生成的 `INPUT` 文件通常会包含以下关键设置。请务必检查确认：

```bash
# INPUT file snippet for deformation tasks

calculation     relax       # 关键点 1
cal_stress      1           # 关键点 2
cal_force       1
```

#### 关键参数解析

1.  **`calculation relax` (至关重要)**
    *   **含义**：固定晶胞形状（Fixed Cell），仅优化原子位置（Relax Atoms）。
    *   **物理意义**：这对应于计算 **"Relaxed-ion elastic constants"**（松弛离子弹性常数）。当我们施加应变时，原子内部会发生微小的相对位移以降低能量（Internal Relaxation）。
    *   **常见错误**：**切勿**在此阶段使用 `cell-relax`！如果你使用了 `cell-relax`，ABACUS 会自动优化晶胞形状，把你辛苦施加的应变“优化”回零应力状态，导致计算失败。
    *   *进阶*：如果你想计算 "Clamped-ion"（高频极限）弹性常数，需在生成输入文件时加上 `--norelax` 参数，此时计算类型为 `scf`。

2.  **`cal_stress 1`**
    *   **含义**：开启应力张量计算。
    *   **必要性**：这是绝对必须的。默认情况下，`relax` 计算可能不输出应力。如果没有这个参数，你将得到一堆优化好的结构，但没有应力数据用于拟合。

### 批量提交

你可以编写一个简单的 Shell 脚本遍历所有文件夹提交任务，或者使用 `abacustest submit` 功能（如果配置了作业调度系统）。

```bash
# 简单的串行提交示例 (仅供参考，建议使用作业调度系统)
for dir in run_elastic/deform-*; do
    cd $dir
    echo "Running in $dir"
    mpirun -np 4 abacus > output.log
    cd ../..
done
```

---

## 3.3 结果后处理与模量拟合

当所有 `deform-xxx` 文件夹中的计算都显示 "Job Finish" 后，我们就可以进行最后的数据收集与拟合。

### 执行分析命令

回到根目录，使用 `abacustest` 的分析模块：

```bash
abacustest analyze -m elastic
```

### 结果解析

程序会自动读取每个文件夹下的 `OUT.*/running_*.log` 文件，提取应力张量，并根据应变模式进行线性拟合。终端将输出类似如下的结果：

```text
-----------------------------------------
Elastic Constants (GPa):
C11 = 155.5
C12 = 58.2
C44 = 76.2
...
Bulk Modulus (Voigt) = 90.6 GPa
Shear Modulus (Voigt) = 68.1 GPa
-----------------------------------------
```

### 结果验证标准

1.  **对称性检查**：对于立方晶系 Si，理论上 $C_{11}=C_{22}=C_{33}$，$C_{12}=C_{13}=C_{23}$，$C_{44}=C_{55}=C_{66}$。如果计算结果严重偏离此规律，说明计算精度不足（如 K 点过少或截断能过低）或初始结构未充分弛豫。
2.  **对比实验值**：
    *   ABACUS 计算值 (PBE): $C_{11} \approx 155$ GPa, $C_{12} \approx 58$ GPa, $C_{44} \approx 76$ GPa。
    *   Materials Project (mp-149): $C_{11} \approx 153$ GPa, $C_{12} \approx 57$ GPa, $C_{44} \approx 74$ GPa。
    *   结果非常接近，证明计算流程正确。

### 常见问题排查 (Troubleshooting)

*   **问题**：输出的弹性常数为 0 或 NaN。
    *   **原因**：检查 `INPUT` 文件中是否遗漏了 `cal_stress 1`，或者计算任务未正常结束。
*   **问题**：$C_{11}$ 和 $C_{12}$ 差异巨大且不合理。
    *   **原因**：通常是因为第一步的“零应力基态”没有优化好。请重新做高精度的 `cell-relax`，确保残余应力小于 0.1 kBar。

---

通过本章的学习，你已经掌握了使用 ABACUS 配合 `abacustest` 进行自动化弹性常数计算的标准流程。这个流程不仅适用于 Si，也适用于其他复杂的晶体体系，只需注意晶体对称性对独立弹性常数个数的影响即可。