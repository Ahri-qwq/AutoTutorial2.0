# 第三章：自动化测试工作流 (Workflow)

> **教授寄语**：
> 在计算材料学领域，手动修改 `INPUT` 文件并逐个提交任务不仅效率低下，更是错误的温床。作为未来的专家，你必须学会“让机器为你工作”。本章将教授如何利用 Shell 脚本结合 ABACUS 接口构建自动化工作流，这是通往高通量计算（High-Throughput Computing）的第一步。

---

## 3.1 构建 Shell 自动化脚本

在进行收敛性测试（如测试平面波截断能 `ecutwfc` 或 K 点密度）时，我们需要进行一系列重复的计算步骤：修改参数 -> 提交任务 -> 等待结束 -> 提取数据。我们将使用 Bash 脚本将这一过程自动化。

### 3.1.1 核心逻辑与 `sed` 命令
我们将使用 Linux 标准流编辑器 `sed` 来动态修改 `INPUT` 文件。
同时，为了防止不同参数的计算结果相互覆盖，**必须**利用 ABACUS 的 `suffix` 参数。设置 `suffix` 后，ABACUS 会将所有输出文件（日志、电荷密度、结构等）写入名为 `OUT.${suffix}` 的文件夹中。

### 3.1.2 自动化脚本模板 (`run_ecut_test.sh`)

请在包含 `INPUT`, `STRU`, `KPT`, `pseudopotentials/` 的目录下创建以下脚本：

```bash
#!/bin/bash

# ==========================================================
# ABACUS 自动化收敛性测试脚本 - Ecutwfc
# 作者: ABACUS Developer Team
# 功能: 遍历不同的 ecutwfc 值，提交计算并防止文件覆盖
# ==========================================================

# 1. 定义要扫描的参数范围 (单位: Ry)
# 注意：不同的赝势（Hard vs Soft）对 Ecut 需求差异巨大。
# 软赝势可能 40-60 Ry 即可，硬赝势可能需要 80-100 Ry 甚至更高。
# 请不要盲目照搬教程值，应从较低值开始扫描。
ECUT_LIST="40 50 60 70 80 90 100"

# 2. 定义并行核数
NP=4

# 备份原始 INPUT 文件，确保每次循环都从干净的模板开始
cp INPUT INPUT.template

echo "Start Convergence Test..."
echo "Ecut(Ry)  Total_Energy(eV)  Time(s)" > convergence_report.dat

for ecut in $ECUT_LIST
do
    # 3. 设置任务后缀 (suffix)
    # 这是一个好习惯：让输出目录名包含关键参数信息
    job_suffix="ecut_${ecut}"
    
    echo "Running calculation for ecutwfc = ${ecut} Ry ..."

    # 4. 动态修改 INPUT 文件
    # 使用 sed 替换 INPUT.template 中的参数并生成新的 INPUT
    # 假设模板中已有 ecutwfc 和 suffix 关键词
    cp INPUT.template INPUT
    
    # 修改 ecutwfc
    # 语法说明: s/原字符串/新字符串/g
    sed -i "s/^ecutwfc.*/ecutwfc ${ecut}/g" INPUT
    
    # 修改 suffix (关键步骤！防止 OUT 文件夹覆盖)
    sed -i "s/^suffix.*/suffix ${job_suffix}/g" INPUT

    # 5. 提交任务
    # 这里的 mpirun 用法取决于你的集群环境，详见 3.2 节
    mpirun -np $NP abacus > log_ecut_${ecut}.txt

    # 6. 简单的数据提取 (可选，详见后续章节)
    # 从标准输出 log 或 OUT 文件夹中的 running_scf.log 提取能量
    # 注意：ABACUS 的标准 SCF 日志位于 OUT.${suffix}/running_scf.log
    logfile="OUT.${job_suffix}/running_scf.log"
    
    if [ -f "$logfile" ]; then
        # 提取最后一步的 E_Fermi 或 Total Energy
        # 这里演示提取 Total Energy (根据实际 log 格式调整 grep 关键词)
        energy=$(grep "total energy" $logfile | tail -1 | awk '{print $5}')
        time_cost=$(grep "total time" $logfile | awk '{print $5}')
        echo "${ecut}      ${energy}      ${time_cost}" >> convergence_report.dat
    else
        echo "Error: ${logfile} not found!"
    fi

done

# 恢复原始 INPUT 文件
mv INPUT.template INPUT
echo "All done. Check convergence_report.dat"
```

> **风险提示 (Risk Warning)**：
> *   **`ecutrho` 的设置**：本脚本仅修改了波函数截断能 `ecutwfc`。在平面波（PW）计算中，电荷密度截断能 `ecutrho` 通常需要是 `ecutwfc` 的 4 倍（对于模守恒赝势）或更高（对于超软赝势）。
> *   ABACUS 默认会自动处理 `ecutrho`（通常默认为 `4 * ecutwfc`），但在高精度需求下，建议显式检查该行为。如果手动设置了 `ecutrho`，请确保在脚本中同步修改它，保持比例一致。

---

## 3.2 任务提交与并行策略

在脚本中的 `mpirun -np $NP abacus` 这一行看似简单，实则包含了并行计算的核心策略。

### 3.2.1 并行核数 (`NP`) 的选择原则
许多初学者认为“核数越多越快”，这是错误的。

1.  **原子数限制**：
    *   MPI 并行效率受限于进程间通讯（Communication Overhead）。
    *   **原则**：并行核数不应显著超过体系的原子数，也不应超过能带数（Bands）。
    *   对于小体系（例如 2 个原子的 Si 原胞），使用 4-8 核通常足够。使用 64 核反而可能因为通讯耗时导致计算变慢。

2.  **K 点并行 (K-point Parallelization)**：
    *   ABACUS 支持 K 点并行。如果你有 10 个 K 点，使用 10 的倍数个核（如 10, 20）通常效率最高，因为每个 K 点可以在一组核上独立计算。

### 3.2.2 运行命令详解
*   **`mpirun` / `mpiexec`**: 启动 MPI 并行环境。
*   **`-np $NP`**: 指定进程数。
*   **`abacus`**: 可执行程序名。如果未配置环境变量，可能需要写绝对路径（如 `/home/user/abacus/bin/abacus`）。

> **集群用户注意**：
> 如果你在超算中心（如 Slurm 系统）上运行，**严禁**在登录节点直接运行上述脚本。你需要将上述脚本封装在 `.slurm` 脚本中，使用 `sbatch` 提交，或者在脚本内部使用 `srun` 替代 `mpirun`。

---

## 3.3 结果分析与收敛判据

脚本运行结束后，你会得到一个 `convergence_report.dat` 文件。如何判断哪个 `ecutwfc` 是合适的？

### 3.3.1 判据标准
收敛不仅仅看能量曲线是否“变平”，必须量化误差。

1.  **能量差判据**：
    *   计算相邻两个参数点的能量差：$\Delta E = |E(ecut_{n}) - E(ecut_{n-1})|$。
    *   **黄金标准**：当 $\Delta E < 1 \text{ meV/atom}$ 时，通常认为该参数已收敛。
    *   *注意单位换算*：ABACUS 输出能量通常为 eV 或 Ry。$1 \text{ meV} = 0.001 \text{ eV}$。

2.  **力与应力判据**（进阶）：
    *   如果后续要进行结构弛豫（Relaxation）或晶格优化（Cell Relaxation），仅能量收敛是不够的，还需要检查原子受力和应力张量随 `ecutwfc` 的变化是否收敛。

### 3.3.2 数据来源确认
虽然 ABACUS 会生成 `istate.info` 等文件，但最权威、信息最全的日志文件是 **`running_scf.log`**（位于 `OUT.${suffix}/` 目录下）。

*   **ETOT (Total Energy)**: 系统的总能量。
*   **EDIFF**: 当前 SCF 步与上一步的能量差（用于判断 SCF 是否收敛，不要与 `ecut` 收敛混淆）。
*   **Time**: 关注计算耗时，在满足精度前提下选择耗时最少的参数。

---

### 课后练习
1.  下载一个简单的硅（Si）8原子超胞结构。
2.  使用本章提供的脚本，扫描 `ecutwfc` 从 20 Ry 到 80 Ry，步长 10 Ry。
3.  绘制 "Total Energy vs Ecut" 曲线，并找出满足 $1 \text{ meV/atom}$ 收敛标准的最小 `ecutwfc`。