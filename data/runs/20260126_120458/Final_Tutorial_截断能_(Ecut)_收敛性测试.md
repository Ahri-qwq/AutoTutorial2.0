# ABACUS 第一性原理计算实战：平面波基组与自动化收敛测试

## 前言

欢迎开启 ABACUS 的学习之旅。作为一款国产自主知识产权的开源第一性原理计算软件，ABACUS 凭借其高效的数值原子轨道（LCAO）与平面波（PW）双基组支持，在材料模拟领域展现了卓越的性能。尤其在处理大规模体系和高精度电子结构计算时，ABACUS 提供了灵活的配置方案。

**学习路线图**
本教程旨在通过四个核心阶段，带你完成从“按下运行键”到“科学分析结果”的蜕变：
1.  **物理溯源**：第一章深入剖析平面波截断能（Cutoff Energy）的物理本质，理解精度与成本的博弈。
2.  **实战配置**：第二章聚焦环境搭建与参数调优，明确 DFT 计算中“数值解”而非“绝对真理”的核心认知。
3.  **效率革命**：第三章教授如何利用 Shell 脚本构建自动化工作流，迈出高通量计算的第一步。
4.  **科学判据**：第四章回归数据本身，确立严谨的收敛性测试判据，确保科研数据的可靠性。

**知识图谱定位**
在整个 ABACUS/DFT 知识体系中，本教程处于“基础进阶”的关键位置。它不仅涵盖了 `ecutwfc` 等核心参数的设置，更侧重于计算方法论的培养。掌握了平面波收敛性测试，你将能够顺利衔接后续的 K 点密度测试、几何优化以及复杂的物性预测。

**前置要求**
在开始本教程之前，我们假设你已经：
*   完成了 ABACUS 软件的基础安装。
*   具备基本的 Linux Shell 操作常识（如 `ls`, `cd`, `grep` 等）。
*   对密度泛函理论（DFT）有初步的理论概念了解。

---

# 第一章：平面波截断的物理原理

欢迎来到 ABACUS 实战教程。作为一名计算材料学研究者，在按下“运行”按钮之前，你必须清楚计算机内部发生了什么。本章我们将探讨平面波计算中最基础、也是极其关键的一个参数——**截断能（Cutoff Energy）**。

这不仅仅是一个数字的设定，它是计算精度与计算成本之间的一场博弈。

## 1.1 平面波基组与截断近似

在 ABACUS 中，当你设置 `basis_type pw` 时，你选择了一条优雅的道路：使用平面波（Plane Waves）作为基组来展开电子波函数。

### 为什么是平面波？
根据布洛赫定理（Bloch's Theorem），周期性晶体势场中的电子波函数 $\psi_{n\mathbf{k}}(\mathbf{r})$ 可以写成一个平面波 $e^{i\mathbf{k}\cdot\mathbf{r}}$ 乘以一个具有晶格周期性的函数 $u_{n\mathbf{k}}(\mathbf{r})$。而这个周期性函数又可以展开为倒格矢 $\mathbf{G}$ 的傅里叶级数。

最终，波函数被表示为无数个平面波的叠加：
$$ \psi_{n\mathbf{k}}(\mathbf{r}) = \sum_{\mathbf{G}} c_{n,\mathbf{k}+\mathbf{G}} \cdot e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}} $$

### 截断的必要性
理论上，求和符号 $\sum_{\mathbf{G}}$ 需要包含无穷多个 $\mathbf{G}$ 矢量才能精确描述波函数（特别是原子核附近的剧烈震荡）。但在计算机中，我们无法处理无穷大的矩阵。

幸运的是，系数 $c_{n,\mathbf{k}+\mathbf{G}}$ 的大小通常随着动能 $\frac{\hbar^2}{2m}|\mathbf{k}+\mathbf{G}|^2$ 的增加而迅速衰减。因此，我们可以引入一个**截断能（Cutoff Energy, $E_{cut}$）**，只保留动能小于该值的平面波分量：

$$ \frac{\hbar^2}{2m}|\mathbf{k}+\mathbf{G}|^2 < E_{cut} $$

这个 $E_{cut}$ 就是你在 ABACUS `INPUT` 文件中看到的 `ecutwfc` 参数。

---

## 1.2 截断能 (Cutoff Energy) 的权衡艺术

选择 `ecutwfc` 是一门平衡的艺术。

*   **过低**：基组太小，无法正确描述波函数的细节（尤其是原子核附近的“硬”区域），导致总能计算错误、力（Force）计算不准，甚至优化出的晶格结构严重偏离实验值。
*   **过高**：基组数量庞大，哈密顿矩阵维度激增。计算内存需求（Memory）和 CPU 时间（Time）呈指数级或幂律级增长，造成资源浪费。

### 我们的目标
我们要寻找一个**收敛甜点（Convergence Sweet Spot）**：在这个数值之上，继续增加截断能，系统的物理量（如总能量、原子受力）变化微乎其微（例如能量差 $< 1 \text{ meV/atom}$），满足科研精度要求。

### ⚠️ 关键风险提示
1.  **赝势决定截断能**：不同的元素、甚至同一种元素的不同赝势（Hard vs. Soft），对 `ecutwfc` 的需求差异巨大。**切勿直接照搬**教程中的 60 Ry 或 100 Ry，必须针对你所使用的赝势进行测试。
2.  **电荷密度截断能**：本章仅讨论波函数截断能 `ecutwfc`。在平面波计算中，描述电荷密度通常需要更高的截断能（`ecutrho`）。ABACUS 内部通常有默认倍数关系（请查阅官方手册），若需手动设置，通常 `ecutrho` $\ge 4 \times$ `ecutwfc`。

---

## 1.3 实战：自动化收敛性测试

在实际科研中，我们**绝不**手动修改 `INPUT` 文件运行十几次。这不仅效率低下，而且容易出错（例如忘记改文件名导致结果覆盖）。

我们将编写一个 Bash 脚本，自动扫描 `ecutwfc`，并提取结果。

### 准备工作
假设你已经准备好了以下文件（以硅为例）：
*   `INPUT`: 包含基本的计算参数。
*   `STRU`: 包含晶格结构和原子位置。
*   `KPT`: K点设置。
*   `Si.upf`: 硅的赝势文件。

**初始 `INPUT` 模板 (文件名: `INPUT_template`)**:
请注意，我们将 `ecutwfc` 和 `suffix` 留作占位符或初始值。

```bash
INPUT_PARAMETERS
# ... 其他参数 ...
basis_type      pw
calculation     scf
# 待扫描参数
ecutwfc         50 
# 输出后缀，防止覆盖
suffix          Si_test
# 电子步收敛精度
scf_thr         1.0e-7
# ... 其他参数 ...
```

### 自动化脚本 (Auto-Convergence Script)

请创建一个名为 `run_ecut_test.sh` 的文件，并写入以下内容。

**脚本逻辑解析**：
1.  **循环**：遍历一系列 `ecutwfc` 值（如 40 到 100 Ry）。
2.  **修改**：使用 `sed` 命令精准修改 `INPUT` 文件中的参数。
3.  **保护**：修改 `suffix` 参数，确保每次计算的输出目录（如 `OUT.Si_ecut40`）独立，互不覆盖。
4.  **运行**：调用 `mpirun` 执行 ABACUS。
5.  **提取**：从标准的 `running_scf.log` 中提取最终能量。

```bash
#!/bin/bash

# 定义要测试的截断能列表 (单位: Ry)
# 建议从较低值开始，逐步增加
ecut_list="40 50 60 70 80 90 100"

# 准备数据记录文件
output_file="ecut_convergence.dat"
echo "# Ecut(Ry)  Total_Energy(eV)" > $output_file

# 循环测试
for ecut in $ecut_list
do
    echo "Running test for ecutwfc = $ecut Ry..."

    # 1. 修改 INPUT 文件
    # 使用 sed 替换 ecutwfc 的值
    # 假设原文件中格式为 "ecutwfc 50" 或类似
    sed -i "s/^ecutwfc.*/ecutwfc ${ecut}/g" INPUT
    
    # 2. 修改 suffix (关键步骤！)
    # 确保输出文件夹名为 OUT.Si_ecut40, OUT.Si_ecut50 等
    # 避免不同任务的日志文件相互覆盖
    suffix_name="Si_ecut${ecut}"
    sed -i "s/^suffix.*/suffix ${suffix_name}/g" INPUT

    # 3. 运行 ABACUS
    # 请根据你的机器环境调整 mpirun 参数
    mpirun -np 4 abacus > log_ecut_${ecut}.txt

    # 4. 数据提取
    # 从 running_scf.log 中提取最终能量
    # grep 查找 "final etot", awk 打印相关列
    # 注意：ABACUS 的输出目录是 OUT.${suffix}
    log_path="OUT.${suffix_name}/running_scf.log"
    
    if [ -f "$log_path" ]; then
        # 提取包含 "final etot" 的行，通常格式为: !FINAL_ETOT_IS  -xxx.xxxx eV
        energy=$(grep "final etot" $log_path | awk '{print $4}')
        
        echo "${ecut}  ${energy}" >> $output_file
        echo "  -> Done. Energy = ${energy} eV"
    else
        echo "  -> Error: $log_path not found!"
    fi
done

echo "All tests finished. Results saved in $output_file"
```

### 运行与分析

1.  赋予脚本执行权限：`chmod +x run_ecut_test.sh`
2.  运行脚本：`./run_ecut_test.sh`
3.  查看结果文件 `ecut_convergence.dat`。

**收敛判据 (Convergence Criteria)**：
不要只看能量是否“变平”。你需要计算相邻两个点的能量差 $\Delta E$。
*   **粗略标准**：$\Delta E < 10 \text{ meV/atom}$
*   **科研标准**：$\Delta E < 1 \text{ meV/atom}$ (甚至更低，取决于研究性质，如声子谱计算需要极高精度)

**示例数据分析**：
```text
# Ecut(Ry)  Energy(eV)      Delta_E (meV)
40          -200.100        -
50          -200.500        400 (未收敛)
60          -200.580        80  (未收敛)
70          -200.595        15  (接近)
80          -200.599        4   (良好)
90          -200.600        1   (收敛甜点)
100         -200.600        0   (浪费资源)
```
在这个假设的例子中，**90 Ry** 是性价比最高的选择。虽然 100 Ry 更准，但 90 Ry 的误差已在 1 meV 级别，且计算成本更低。
通过这种严谨的测试流程，你才能确信你的计算结果是建立在坚实的物理基础之上的。

# 第二章：测试环境与参数配置

作为计算材料学的从业者，我们必须建立一个核心认知：**DFT 计算的结果并非绝对真理，而是基于一系列近似参数的数值解。** 其中，平面波截断能（Energy Cutoff）是决定计算精度与成本平衡的最关键参数之一。

本章将手把手教你如何科学地配置测试环境，区分“变量”与“控制变量”，并编写自动化脚本完成测试。

---

## 2.1 关键输入参数设置

在进行截断能测试时，我们必须确保物理模型的基组类型是正确的。ABACUS 支持多种基组（如 LCAO 原子轨道基组），但在进行收敛性测试，特别是为 LCAO 生成参考基准时，我们通常使用**平面波（Plane Wave, PW）**基组。

请打开你的 `INPUT` 文件，关注以下核心参数：

### 1. 基组类型 (`basis_type`)
必须显式设置为 `pw`。这是进行截断能测试的前提。
```bash
basis_type pw
```

### 2. 波函数截断能 (`ecutwfc`)
这是本章唯一的“自变量”。
*   **物理含义**：限制平面波基组动能的最大值。动能小于该值的平面波将被纳入基组。
*   **单位**：ABACUS 默认单位为 **Ry (Rydberg)**。注意：1 Ry $\approx$ 13.6 eV。
*   **推荐扫描范围**：
    *   不要盲目照搬教程中的 60 Ry。不同的元素和赝势（Pseudopotential）对截断能的要求差异巨大。
    *   **软赝势（Soft）**：可能在 30-40 Ry 收敛。
    *   **硬赝势（Hard）**：可能需要 80-120 Ry 甚至更高。
    *   **建议策略**：从 **20 Ry** 开始，以 **10 Ry** 为步长，一直扫描到 **120 Ry**（或更高，视具体体系而定）。

**INPUT 文件示例片段**：
```bash
INPUT_PARAMETERS
# ... 其他参数 ...
basis_type pw
ecutwfc    50   # 该值将在测试脚本中被动态修改
# ... 其他参数 ...
```

---

## 2.2 关联参数与潜在风险：电荷密度截断能

在平面波计算中，除了波函数截断能 `ecutwfc`，还存在一个极易被新手忽略的参数：**电荷密度截断能 (`ecutrho`)**。

*   **物理原理**：电荷密度 $n(\mathbf{r}) = |\psi(\mathbf{r})|^2$。在倒空间中，描述电荷密度所需的截止波矢量通常是波函数的 2 倍，对应的能量截断能则是波函数的 4 倍（对于模守恒赝势 Norm-Conserving PP）。如果使用超软赝势（USPP）或 PAW 方法，这个比例通常需要更高。
*   **风险提示**：如果 `ecutrho` 设置过低，会导致积分网格过于稀疏，产生“蛋盒效应”（Egg-box effect），导致计算结果出现虚假的震荡或不收敛。
*   **设置建议**：
    *   **本教程聚焦**：本章重点关注波函数截断能 `ecutwfc` 的收敛行为。
    *   **操作指南**：请查阅 ABACUS 官方文档关于 `ecutrho` 的默认行为。通常建议显式设置 `ecutrho` 为 `ecutwfc` 的 4 倍以上（例如 `ecutwfc 60` 时，`ecutrho` 设为 240）。
    *   *注：若未显式设置，请务必检查输出日志确认程序使用的默认值是否合理。*

---

## 2.3 控制变量：结构与 K 点

科学实验的原则是“控制变量法”。当我们测试 `ecutwfc` 时，以下两个文件中的参数必须保持**绝对静止**：

1.  **`STRU` (晶胞与原子位置)**：
    *   严禁在测试能量收敛时通过 `calculation cell-relax` 或 `relax` 改变晶胞形状或原子位置。
    *   **原因**：晶胞体积的变化会直接改变倒空间格点分布，导致能量基准漂移，使得不同截断能下的能量不可比。

2.  **`KPT` (K 点网格)**：
    *   K 点密度必须固定（例如 `K_POINTS` 文件中的 `2 2 2` 或 `4 4 4`）。
    *   **原因**：K 点采样密度直接影响布里渊区积分精度。如果同时改变 K 点和截断能，你将无法判断能量变化是由哪个因素引起的。

---

## 2.4 自动化测试脚本（实战核心）

手动修改 `INPUT` 文件并运行 10 次计算既低效又容易出错。作为专业开发者，我们强烈建议使用 Shell 脚本自动化此过程。

### 脚本逻辑解析
1.  **循环**：遍历预设的能量列表（20 Ry 到 100 Ry）。
2.  **修改参数**：使用 `sed` 命令原位修改 `INPUT` 文件。
3.  **关键保护 (`suffix`)**：**这是新手最常犯的错误**。必须修改 `suffix` 参数，使每次计算的输出目录不同（如 `OUT.test_20Ry`, `OUT.test_30Ry`）。否则，新计算会覆盖旧计算的日志，导致数据丢失。
4.  **执行与提取**：运行 MPI 任务，并从标准的 SCF 输出日志 `running_scf.log` 中提取最终能量。

### 通用测试脚本模板 (`run_ecut_test.sh`)

请在包含 `INPUT`, `STRU`, `KPT` 和赝势文件的目录下创建此脚本：

```bash
#!/bin/bash

# 设定扫描范围：从 20 Ry 到 100 Ry，步长 10 Ry
ecut_list=$(seq 20 10 100)

# 清理旧数据文件
rm -f energy_convergence.dat

echo "Ecut(Ry)   Energy(eV)" > energy_convergence.dat

for ecut in ${ecut_list}; do
    echo "Running test for ecutwfc = ${ecut} Ry..."
    
    # 1. 修改 INPUT 文件中的 ecutwfc
    # 使用 sed 正则表达式匹配以 ecutwfc 开头的行并替换
    sed -i "s/^ecutwfc.*/ecutwfc ${ecut}/g" INPUT
    
    # 2. 修改 suffix 以防止输出文件覆盖
    # 这里的 suffix 设置为 test_${ecut}，输出目录将变为 OUT.test_20 等
    sed -i "s/^suffix.*/suffix test_${ecut}/g" INPUT
    
    # 3. 提交计算任务
    # 请根据你的机器实际核数调整 -np 后的数字
    mpirun -np 8 abacus > log_parallel_${ecut}.txt
    
    # 4. 数据提取
    # 从标准 SCF 日志 running_scf.log 中提取最终能量
    # 这里的 grep 关键词需严格匹配 ABACUS 输出格式
    # 假设输出目录格式为 OUT.suffix
    log_file="OUT.test_${ecut}/running_scf.log"
    
    if [ -f "$log_file" ]; then
        # 提取包含 "final etot is" 的行，并获取能量值（通常在第4列，单位eV）
        # 示例行: !FINAL_ETOT_IS -1.2345678 eV
        energy=$(grep "final etot is" $log_file | awk '{print $4}')
        echo "${ecut}      ${energy}" >> energy_convergence.dat
    else
        echo "Error: ${log_file} not found!"
    fi
done

echo "Test finished. Data saved to energy_convergence.dat"
```

### 运行脚本
```bash
chmod +x run_ecut_test.sh
./run_ecut_test.sh
```

---

## 2.5 结果分析与收敛判据

脚本运行结束后，你将得到 `energy_convergence.dat` 文件。如何判断哪个值是“收敛”的？

1.  **不要只看绝对能量**：DFT 的绝对能量没有物理意义，有意义的是能量差。
2.  **计算能量差**：计算 $E(ecut) - E(ecut_{max})$。
3.  **收敛标准**：
    *   **粗略标准**：能量随 `ecutwfc` 增加趋于平缓。
    *   **精确标准（推荐）**：当 $\Delta E < 1 \text{ meV/atom}$（即 0.001 eV/atom）时，可认为收敛。
    *   *注意*：对于声子谱或弹性常数计算，收敛标准通常需要更严格（如 $10^{-4}$ eV/atom）。


# 第三章：自动化测试工作流 (Workflow)

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

---

## 附录：进阶学习指南

完成本教程的学习仅仅是一个开始，计算材料学的世界广阔无垠，我们鼓励你继续探索以下领域：

### 1. 深度扩展与持续学习
*   **K 点收敛测试**：平面波截断能确保了基组在实空间的完备性，而 K 点密度（K-point Mesh）则决定了倒空间积分的精度。建议参考本教程的自动化思路，对 `KPT` 文件进行类似的收敛性测试。
*   **LCAO 基组进阶**：ABACUS 的一大特色是数值原子轨道。在掌握平面波后，建议学习如何利用 `basis_type lcao` 进行大规模体系计算，并了解 `orbital` 文件的选取与优化。
*   **混合基组与性能调优**：查阅官方文档中关于 `mixing_gg0` 或 `ng_beta` 等参数的说明，这些参数对于改善金属体系或复杂磁性体系的自洽迭代（SCF）收敛速度至关重要。

### 2. 通用调试建议 (Troubleshooting)
*   **检查运行日志**：当计算意外终止时，首要任务是查看 `OUT.suffix` 文件夹下的 `running_scf.log`。大多数错误（如内存不足、对称性识别失败）都会在此给出提示。
*   **并行资源分配**：根据参考资料，建议总核数取 2 的 n 次方。一个粗略的经验法则是：体系有多少个原子，使用的总核数不宜远超原子个数，以避免通信开销超过计算增益。
*   **能量判据**：通常认为单原子能量差小于 1 meV/atom 时即达到收敛。若能量曲线波动剧烈，请检查伪势文件是否匹配或初始磁矩设置是否合理。

### 3. 参考资源
*   **ABACUS 官方文档**：获取最权威的参数定义与更新动态。
*   **Bohrium 案例库**：本教程的部分脚本逻辑参考了 Bohr 平台上的氮氧化物及钡元素计算案例，建议前往该平台获取更多预装环境下的实战练习。

保持好奇，享受计算的乐趣！
