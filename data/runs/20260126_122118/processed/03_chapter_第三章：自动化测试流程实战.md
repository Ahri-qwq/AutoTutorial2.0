# 第三章：自动化测试流程实战

在前两章中，我们已经成功编译了 ABACUS 并手动完成了一次简单的自洽计算（SCF）。但在实际科研中，手动修改参数、提交任务、记录数据的过程不仅枯燥，而且极易出错。

本章将带你进入**计算材料学的“工业化”时代**。我们将编写 Shell 脚本，构建一个“修改-运行-提取”的自动化闭环。我们将以**K 点收敛性测试**为例，演示如何通过脚本一键完成从参数扫描到数据汇总的全过程。

---

## 3.1 编写自动化批处理脚本

在进行收敛性测试时，我们需要控制变量。对于 K 点测试，原则是：**固定截断能（`ecutwfc`）和其他参数，仅改变 K 点网格密度**。

### 3.1.1 准备工作目录
首先，建立一个干净的测试目录，并准备好基础文件：
- `STRU`：结构文件（如硅或铝）。
- `INPUT`：输入参数文件。
- `KPT`：我们将用脚本动态生成此文件，初始可以不存在。
- 赝势文件（`.upf` 或 `.orb` 等）。

**关键设置**：
在 `INPUT` 文件中，我们通常设置一个 `suffix` 参数。
- **参数名**: `suffix`
- **作用**: 指定输出目录的后缀。例如设置 `suffix Si_test`，输出文件夹将命名为 `OUT.Si_test`。
- **自动化意义**: 在批量计算中，我们需要动态修改这个参数，防止后一个任务的结果覆盖前一个任务。

### 3.1.2 编写 Shell 脚本 (`auto_test_kpt.sh`)

我们将使用 Linux 系统中最常用的 Shell 脚本配合 `sed` 命令来实现自动化。请在测试目录下创建名为 `auto_test_kpt.sh` 的文件，并写入以下内容：

```bash
#!/bin/bash

# 1. 定义要测试的 K 点网格密度列表
# 这里我们测试从 2x2x2 到 8x8x8
k_values="2 3 4 5 6 7 8"

# 2. 清理旧的数据文件 (可选)
rm -f energy_convergence.dat

# 3. 开始循环
for k in $k_values
do
    echo "Testing K-point grid: $k x $k x $k"

    # --- 步骤 A: 动态生成 KPT 文件 ---
    # ABACUS 的 KPT 格式通常如下：
    # K_POINTS
    # 0
    # Gamma
    # n n n 0 0 0
    # 我们使用 cat 命令直接重写 KPT 文件
    cat > KPT <<EOF
K_POINTS
0
Gamma
$k $k $k 0 0 0
EOF

    # --- 步骤 B: 动态修改 INPUT 文件中的 suffix ---
    # 使用 sed 命令找到以 'suffix' 开头的行，并将其替换为新的后缀
    # 这样输出目录会变为 OUT.Si_k2, OUT.Si_k3 等
    sed -i "s/^suffix.*/suffix Si_k${k}/" INPUT

    # --- 步骤 C: 运行计算 ---
    # 假设使用 mpirun 进行 4 核并行计算
    # 将标准输出重定向到 log 文件中，方便后续提取
    mpirun -np 4 abacus > running_k${k}.log

    # --- 步骤 D: 提取数据 ---
    # 从日志文件中筛选包含 "final etot" 的行
    # awk '{print $4}' 表示提取该行的第 4 列数据（即能量值）
    # 注意：根据 ABACUS 版本不同，能量可能在不同列，请先手动检查一次 log 确认位置
    energy=$(grep "final etot" running_k${k}.log | awk '{print $4}')

    # 将 K 点密度和对应的能量写入汇总文件
    echo "$k $energy" >> energy_convergence.dat
done

echo "All calculations finished. Results saved in energy_convergence.dat"
```

### 3.1.3 脚本核心解析
1.  **`cat > KPT <<EOF ... EOF`**: 这是一个非常实用的“Here Document”写法，它允许我们在脚本中直接“画”出文件的内容。每次循环，$k 的值都会改变，从而生成新的 `KPT` 文件。
2.  **`sed -i ...`**: `sed` 是流编辑器。`"s/^suffix.*/suffix Si_k${k}/"` 的意思是：查找所有以 `suffix` 开头的行（`^`表示行首），将其替换为 `suffix Si_k2`（当 k=2 时）。这是防止结果覆盖的**关键操作**。

---

## 3.2 运行计算与数据分析

### 3.2.1 执行脚本
赋予脚本执行权限并运行：
```bash
chmod +x auto_test_kpt.sh
./auto_test_kpt.sh
```

此时，终端会依次打印正在计算的 K 点信息。计算完成后，目录下会出现一系列日志文件（`running_k2.log`, `running_k3.log`...）和输出文件夹（`OUT.Si_k2`, `OUT.Si_k3`...）。

### 3.2.2 结果分析与收敛判据
查看生成的汇总文件 `energy_convergence.dat`：
```bash
cat energy_convergence.dat
```
输出示例（假设值）：
```text
2 -215.12345678
3 -215.45678901
4 -215.48901234
5 -215.49012345
6 -215.49045678
7 -215.49051234
8 -215.49052345
```

**Critical: 收敛判据**
如何判断 K 点是否取够了？
1.  **能量差标准**: 计算相邻两个 K 点网格的总能量差值 $\Delta E = E_{n} - E_{n-1}$。
2.  **归一化**: 将能量差除以体系中的原子数。
3.  **阈值**: 当能量变化小于 **1 meV/atom** (即 $0.001 \text{ eV/atom}$) 时，通常认为计算已经收敛。

在上述示例中，从 5x5x5 到 6x6x6，能量变化极小（约 0.3 meV），因此对于该体系，选取 5x5x5 或 6x6x6 作为后续计算参数是合理的。

### 3.2.3 锦囊妙计：可视化建议
强烈建议不要只看数字。将 `energy_convergence.dat` 导入绘图软件（如 Gnuplot, Python Matplotlib 或 Origin）。
- **横坐标**: K-point Grid Density (2, 3, 4...)
- **纵坐标**: Total Energy (eV)
绘制 **"Total Energy vs. K-point Grid"** 曲线图。你会看到一条能量迅速下降（或上升）并最终趋于平坦的曲线。这种可视化的直觉对于判断收敛性至关重要。

---

## 3.3 进阶技巧与风险提示

### Pro Tip 1: 扩胞后的 K 点缩减法则
很多初学者在做完单胞（Unit Cell）测试后，转去做超胞（Supercell）计算时，往往不知道如何设置 K 点。
**经验法则**: K 点网格密度与晶胞尺寸成**反比**。
> 如果你的晶胞在 x 方向扩大了 $n$ 倍（例如建立了一个 $2 \times 1 \times 1$ 的超胞），那么你在 x 方向的 K 点数目可以减少为原来的 $1/n$。
> **例子**: 单胞收敛 K 点为 $8 \times 8 \times 8$。建立 $2 \times 2 \times 2$ 超胞后，K 点只需设为 $4 \times 4 \times 4$ 即可保持同等的计算精度。

### Pro Tip 2: Gamma Only 模式
对于非常大的体系（例如超过 200 个原子），或者非周期性体系（如团簇、分子），布里渊区非常小，通常只需要计算 $\Gamma$ 点（0 0 0）。
ABACUS 针对这种情况有专门的优化参数。
- 在 `INPUT` 文件中设置：
  ```text
  gamma_only 1
  ```
- 此时 `KPT` 文件应明确写为：
  ```text
  K_POINTS
  1
  Direct
  0 0 0 1
  ```
开启 `gamma_only` 可以利用实数空间的算法显著降低内存消耗并提升计算速度。

### Risk Warning: KPT 文件的 Line Mode
本章主要讨论的是自洽计算（SCF）用的均匀 K 点网格。如果你在进行能带结构（Band Structure）计算，`KPT` 文件的格式会完全不同（通常称为 Line Mode）。
- **切勿模仿 VASP**: ABACUS 的能带路径格式与 VASP 的 `KPOINTS` 不同。
- **查阅文档**: 请务必查阅 ABACUS 官方文档中关于 `KPT` 格式的说明，不要凭直觉瞎编格式，否则程序可能会报错或计算出错误的能带路径。

---

通过本章的实战，你不仅掌握了 K 点测试的方法，更重要的是掌握了**Shell 脚本自动化**这一核心技能。在后续的截断能（Ecut）测试、晶格常数优化等任务中，你只需微调上述脚本即可轻松应对。