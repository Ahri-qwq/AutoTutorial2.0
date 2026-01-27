# 第二章：测试环境与参数配置

在上一章中，我们完成了 ABACUS 的安装与基本运行测试。本章将进入实战的核心环节：**收敛性测试（Convergence Test）**。

作为计算材料学的从业者，你必须建立一个核心认知：**DFT 计算的结果并非绝对真理，而是基于一系列近似参数的数值解。** 其中，平面波截断能（Energy Cutoff）是决定计算精度与成本平衡的最关键参数之一。

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

**下一步**：
在确定了合适的 `ecutwfc` 后，我们将固定该参数，在下一章中对 **K 点网格 (KPT)** 进行同样的收敛性测试。