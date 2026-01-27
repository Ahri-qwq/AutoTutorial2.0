# 第四章：数据分析与收敛判据

在前几章中，我们已经成功配置了输入文件并运行了 ABACUS。然而，计算的结束只是科研工作的开始。面对海量的日志文件，如何提取关键物理量？如何判断计算结果是否可靠？

本章将带你掌握计算材料学中最基础也最重要的技能：**收敛性测试（Convergence Test）**。我们将通过自动化脚本批量处理数据，并制定科学的收敛判据。

---

## 4.1 自动化测试与数据提取

在进行截断能（`ecutwfc`）测试时，初学者常犯的错误是手动修改 `INPUT` 文件并反复提交任务。这不仅效率低下，而且极易出错。

作为资深开发者，我必须强调：**请永远使用脚本自动化你的工作流**。

### 4.1.1 编写自动化测试脚本

我们需要编写一个 Bash 脚本，完成以下逻辑：
1.  设定一系列 `ecutwfc` 数值（例如从 50 Ry 到 100 Ry）。
2.  为每个数值创建一个独立的工作目录，防止文件混乱。
3.  使用 `sed` 命令自动修改 `INPUT` 文件中的参数。
4.  **关键步骤**：修改 `suffix` 参数，确保输出文件不会相互覆盖。
5.  提交计算并提取最终能量。

请在案例目录下创建名为 `run_conv_test.sh` 的脚本：

```bash
#!/bin/bash

# 定义要测试的截断能列表 (单位: Ry)
# 提示：对于软赝势可从 40 起步，对于硬赝势或模守恒赝势建议从 60 或 80 起步
ecut_list="50 60 70 80 90 100"

# 原始输入文件所在目录（假设当前脚本就在该目录下）
template_dir="."

# 循环遍历每一个截断能
for ecut in $ecut_list
do
    # 1. 创建独立的工作目录，例如 run_60Ry
    work_dir="run_${ecut}Ry"
    if [ -d "$work_dir" ]; then
        rm -rf $work_dir
    fi
    mkdir $work_dir
    
    # 2. 将必要的输入文件复制到工作目录
    # 注意：根据你的实际情况，可能还需要复制 KPT 或 orbital 文件夹
    cp $template_dir/INPUT $template_dir/STRU $template_dir/KPT $work_dir/
    # 如果赝势在上一级目录，不需要复制，只需确保 INPUT 中 pseudo_dir 正确即可
    
    # 3. 进入工作目录
    cd $work_dir
    
    # 4. 使用 sed 修改 INPUT 文件 (Critical Step)
    # 修改 ecutwfc
    sed -i "s/^ecutwfc.*/ecutwfc ${ecut}/g" INPUT
    
    # 修改 suffix (防止输出文件覆盖，利于后续归档)
    sed -i "s/^suffix.*/suffix Si_ecut${ecut}/g" INPUT
    
    echo "Running calculation for ecutwfc = ${ecut} Ry..."
    
    # 5. 运行 ABACUS (根据你的机器环境修改 mpirun 参数)
    # 这里的 abacus 是可执行程序的名称
    mpirun -np 4 abacus > running_scf.log 2>&1
    
    # 6. 立即提取数据
    # 从日志中筛选 "final etot is"，并提取能量值（通常在第5列，单位 eV）
    energy=$(grep "final etot is" running_scf.log | awk '{print $5}')
    
    # 返回上一级目录
    cd ..
    
    # 7. 将结果写入汇总文件
    echo "${ecut} ${energy}" >> data_convergence.txt
done

echo "All calculations finished. Results saved in data_convergence.txt"
```

> **专家提示**：
> *   **`sed -i`**: 这是一个强大的流编辑器命令，用于直接修改文件内容。
> *   **`suffix` 的重要性**: ABACUS 的输出文件通常保存在 `OUT.${suffix}` 文件夹下。如果不修改 `suffix`，并行任务可能会尝试写入同一个文件夹，或者在后续分析时让你混淆不同参数的结果。

### 4.1.2 数据清洗

运行上述脚本后，你将得到一个名为 `data_convergence.txt` 的文件。其内容格式如下：

```text
50 -215.123456
60 -215.456789
70 -215.467890
80 -215.470123
90 -215.470567
100 -215.470650
```

这就是我们进行科学决策的基础数据：第一列是截断能（Ry），第二列是总能量（eV）。

---

## 4.2 确立收敛标准

拥有数据后，如何确定哪个 `ecutwfc` 是“足够好”的？

### 4.2.1 绘制收敛曲线

建议使用 Python (Matplotlib) 或 Gnuplot 绘制 **Total Energy vs. Ecut** 的曲线。
你会观察到：随着 `ecutwfc` 增加，总能量通常会单调下降（对于变分原理成立的计算），并逐渐趋于一个常数。

### 4.2.2 具体的收敛判据

**误区警示**：很多新手认为只要能量曲线“看起来平了”就可以。这是不严谨的。

**科学标准**：我们需要定义一个能量差的阈值 $\delta$。
$$ \Delta E = |E_{total}(E_{cut}) - E_{total}(E_{cut}^{max})| $$
或者考察相邻两个测试点的能量差。

通常的收敛标准如下：
1.  **粗略计算**：$\Delta E < 10 \text{ meV/atom}$
2.  **标准生产环境**：$\Delta E < 1 \text{ meV/atom}$ (最常用)
3.  **高精度计算（如声子谱）**：$\Delta E < 0.1 \text{ meV/atom}$

**操作步骤**：
1.  计算体系中的原子总数 $N$。
2.  检查 `data_convergence.txt` 中相邻数据点的能量差。
3.  找到那个使得能量变化小于 $1 \text{ meV} \times N$ 的最小 `ecutwfc`。
4.  为了保险起见，通常会在该值的基础上再增加 5-10 Ry 作为最终生产参数。

### 4.2.3 赝势的影响

**Critical**: 请绝对不要照搬本教程或任何文献中的 `60 Ry`。
*   **Ultrasoft (USPP) / PAW 赝势**：通常比较“软”，可能在 40-60 Ry 就能收敛。
*   **Norm-Conserving (NCPP) 赝势**：通常比较“硬”，可能需要 80-120 Ry 甚至更高才能收敛。

**必须针对你使用的每一个新赝势文件重新进行收敛性测试。**

---

## 附录：常见问题与进阶建议

### 1. 内存溢出 (OOM - Out Of Memory)
提高 `ecutwfc` 会显著增加内存消耗。平面波基组的数量与 $E_{cut}^{3/2}$ 成正比。如果你在测试 120 Ry 时节点崩溃，请检查内存使用情况，或尝试增加并行节点数以分摊内存。

### 2. 电荷密度截断能 (ecutrho)
本章主要关注波函数截断能 `ecutwfc`。在 ABACUS 的平面波 (PW) 计算中，还有一个重要参数 `ecutrho`（电荷密度截断能）。
*   通常情况下，`ecutrho` 应至少为 `ecutwfc` 的 4 倍（对于模守恒赝势）。
*   对于超软赝势 (USPP) 或 PAW，该比例可能需要更高（如 8-10 倍）。
*   **注意**：若未在 INPUT 中显式设置 `ecutrho`，ABACUS 会根据 `ecutwfc` 和赝势类型自动设定默认值。建议新手在初期使用默认设置，但在高阶优化时需查阅官方手册确认。

### 3. 力 (Force) 与应力 (Stress) 的收敛
如果你的任务是结构弛豫（`calculation cell-relax` 或 `relax`），仅看总能量收敛是不够的。
*   **力**和**应力**对截断能的敏感度通常高于总能量。
*   建议在确定 `ecutwfc` 后，额外检查最高截断能下的原子受力是否稳定。如果发现力在不同截断能下波动剧烈，需要进一步提高 `ecutwfc`。