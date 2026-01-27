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