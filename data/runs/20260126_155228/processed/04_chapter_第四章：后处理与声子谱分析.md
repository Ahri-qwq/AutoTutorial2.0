# 第四章：后处理与声子谱分析

在前一章中，我们已经成功生成了微扰结构并完成了所有超胞的自洽场（SCF）计算。此时，每个微扰文件夹（如 `disp-001`）下的日志文件中都蕴含着该微扰引起的原子受力信息。

本章的任务是将这些分散的“力”收集起来，交回给 Phonopy 计算力常数矩阵（Force Constants），进而通过动力学矩阵的对角化得到声子谱（Phonon Dispersion），并最终进行可视化与物理分析。

---

## 4.1 构建力常数矩阵 (FORCE_SET)

这是连接 ABACUS 计算结果与 Phonopy 分析的关键一步。我们需要从 ABACUS 的输出日志 `running_scf.log` 中提取原子受力。

### 4.1.1 收集受力信息

Phonopy 提供了一个便捷的命令行工具 `-f` 来解析计算软件的输出。对于 ABACUS，Phonopy 会自动识别 `running_scf.log` 中的 `FORCE` 字段。

**操作步骤：**

请在包含 `disp-xxx` 文件夹的根目录下执行以下命令：

```bash
# 语法：phonopy -f [所有微扰任务的日志文件路径]
phonopy -f ./disp-001/OUT*/running_scf.log ./disp-002/OUT*/running_scf.log
```

**命令详解：**
*   `-f`: 告诉 Phonopy 接下来的一系列文件是用来提取受力的。
*   路径通配符: `./disp-001/OUT*/running_scf.log` 是为了匹配 ABACUS 默认生成的输出目录结构（通常是 `OUT.${suffix}`）。
*   **注意**：如果你的体系对称性较低，生成了多个微扰文件夹（如 `disp-001` 到 `disp-010`），你需要将它们全部列出，或者使用 Shell 通配符（如 `./disp-*/OUT*/running_scf.log`，但需确保顺序正确，通常 Phonopy 会自动匹配，但显式列出更为稳妥）。

### 4.1.2 生成结果

命令执行成功后，当前目录下会生成一个名为 `FORCE_SET` 的文件。
*   **文件内容**：包含了原子位移与对应受力的映射关系。
*   **检查方法**：可以使用文本编辑器查看该文件，确保其中没有 `NaN` 或全零数据（除非该方向确实无受力）。

---

## 4.2 声子谱计算配置 (band.conf)

有了 `FORCE_SET`，我们就可以计算声子谱了。为了规范化计算流程，我们通常编写一个名为 `band.conf` 的配置文件，详细定义计算参数、高对称点路径等。

### 4.2.1 编写 band.conf

新建文件 `band.conf`，并写入以下内容（以 FCC Al 为例）：

```ini
# band.conf
ATOM_NAME = Al
DIM = 2 2 2
MESH = 8 8 8

# 原胞基矢矩阵 (Primitive Axes)
# 对于 FCC 结构，从常规晶胞到原胞的转换矩阵
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0

# 高对称点路径 (K-path / Q-path)
# 路径: Gamma(0 0 0) -> X(1/2 1/2 1) -> K(3/8 3/8 3/4) -> Gamma(0 0 0) -> L(1/2 1/2 1/2)
BAND = 1 1 1  1/2 1/2 1  3/8 3/8 3/4  0 0 0   1/2 1/2 1/2

# 插值密度
BAND_POINTS = 101

# 辅助能带连接（处理能带交叉）
BAND_CONNECTION = .TRUE.
```

### 4.2.2 关键参数详解

1.  **`DIM` (Dimension)**:
    *   **必须**与第二章中生成超胞时使用的参数完全一致（如 `2 2 2`）。如果这里不匹配，Phonopy 将无法正确映射力常数。

2.  **`PRIMITIVE_AXES`**:
    *   **物理意义**：我们在计算受力时使用的是扩胞后的超胞（通常基于常规晶胞），但在分析能带结构（声子谱）时，我们通常希望在**原胞（Primitive Cell）**的布里渊区中进行。
    *   该矩阵定义了如何从输入结构（常规胞）转换到原胞。对于 FCC 金属 Al，输入的是面心立方常规胞，通过上述矩阵可转换为包含 1 个原子的原胞。

3.  **`BAND`**:
    *   定义倒空间中的高对称点路径。
    *   **获取方式**：对于未知结构，强烈推荐使用 [SeeK-path](https://www.materialscloud.org/work/tools/seekpath) 工具，上传你的 `STRU` 文件，它会给出标准的 `PRIMITIVE_AXES` 和 `BAND` 路径坐标。

4.  **`ATOM_NAME`**:
    *   指定输出图片或数据中的元素标签，方便阅读。

---

## 4.3 数据可视化与结果解读

一切准备就绪，现在开始计算并绘制声子谱。

### 4.3.1 运行计算

在终端执行以下命令：

```bash
# -p 表示直接绘图(plot)，--abacus 指定接口模式
phonopy -p band.conf --abacus
```

**输出产物：**
1.  `band.yaml`: 包含声子谱所有详细数据的 YAML 文件（频率、特征向量等）。
2.  `band.pdf` (或 png): Phonopy 默认生成的声子谱图像。

### 4.3.2 高级绘图 (Gnuplot)

默认的图像可能不满足发表级质量的要求。我们可以导出数据并使用 Gnuplot 或 Python (Matplotlib) 进行精细绘制。这里演示 Gnuplot 的流程。

**第一步：导出 Gnuplot 数据格式**
```bash
phonopy-bandplot --gnuplot > pho.dat
```
此时生成的 `pho.dat` 包含了路径距离（第一列）和各支声子的频率（后续列）。

**第二步：使用 Gnuplot 绘图**
编写脚本 `plot_pho.gp`：

```gnuplot
set terminal pngcairo size 1920, 1080 font 'Arial, 36'
set output "Al-FCC_phonon.png"

set ylabel 'Frequency (THz)'
set ytics 2
unset key

# 定义高对称点的位置 (需根据 pho.dat 或 band.yaml 中的距离信息手动调整)
# 注意：以下 x 坐标值仅为示例，实际值请查看 pho.dat 的注释行或 band.yaml
x_G = 0.0000
x_X = 0.1312
x_K = 0.1775
x_G2 = 0.3166
x_L = 0.4302

set xrange [0:x_L]
set yrange [0:12]

# 设置 X 轴标签
set xtics ("{/Symbol G}" x_G, "X" x_X, "K" x_K, "{/Symbol G}" x_G2, "L" x_L)

# 绘制竖直线标记高对称点
set arrow from x_X,0 to x_X,12 nohead lt 2
set arrow from x_K,0 to x_K,12 nohead lt 2
set arrow from x_G2,0 to x_G2,12 nohead lt 2

# 绘图
plot 'pho.dat' using 1:2 w l lw 3 lc rgb "blue", \
     'pho.dat' using 1:3 w l lw 3 lc rgb "blue", \
     'pho.dat' using 1:4 w l lw 3 lc rgb "blue"
```

运行绘图：
```bash
gnuplot plot_pho.gp
```

### 4.3.3 结果解读 (以 FCC Al 为例)

当你看到生成的声子谱时，应关注以下特征：

1.  **声学支 (Acoustic Modes)**:
    *   对于包含 $N$ 个原子的原胞，总共有 $3N$ 支声子。
    *   Al 的原胞只有 1 个原子，因此只有 **3 支声学支**。
    *   **特征**：声学支必须从 $\Gamma$ 点（0,0,0）出发，且频率为 0 THz。这对应于晶格的整体平移运动。

2.  **光学支 (Optical Modes)**:
    *   如果原胞中原子数 $N > 1$，则会有 $3N - 3$ 支光学支。
    *   在 FCC Al 的例子中，**没有光学支**。如果你看到了高频的平坦能带，请检查你的 `PRIMITIVE_AXES` 设置是否正确（是否错误地在超胞布里渊区画图）。

3.  **频率范围**:
    *   Al 的最大截止频率大约在 10-12 THz 左右。如果数量级偏差巨大（如 1000 THz 或 0.01 THz），请检查单位转换或原子质量设置。

---

## 附录：常见问题与进阶建议

### 1. 虚频排查指南 (Imaginary Frequencies)
在声子谱中，虚频通常表现为**负频率**。如果在 $\Gamma$ 点附近出现微小的负频（如 -0.05 THz），通常是数值噪音，可以忽略。但如果出现明显的负频带，说明结构不稳定，可能原因如下：

*   **结构未充分弛豫 (最常见)**:
    *   **原因**: 初始结构的原子受力不为零，导致力常数计算错误。
    *   **对策**: 提高结构优化的精度。建议 `force_thr_ev` (或 `force_thr`) 至少达到 `1e-3` eV/Ang 甚至更低。
*   **超胞尺寸过小**:
    *   **原因**: 晶格振动的长程相互作用被截断。
    *   **对策**: 增大 `DIM`，例如从 `2 2 2` 增加到 `3 3 3`。
*   **收敛精度不足**:
    *   **原因**: SCF 计算时的 `scf_thr` 太大，导致受力计算不准。
    *   **对策**: 确保微扰计算时 `scf_thr` 设置为 `1e-8` 或更严。

### 2. 磁性体系注意事项
对于磁性材料（如 Fe, Ni, Co 或磁性氧化物）：
*   **磁矩一致性**: 必须确保所有微扰结构（`disp-xxx`）收敛到了与基态相同的磁性状态（铁磁、反铁磁等）。
*   **操作技巧**: 在 `INPUT` 文件中明确指定初始磁矩，并可能需要读取基态的电荷密度/波函数作为初猜（如果 ABACUS 支持该功能），或者仔细检查每个微扰任务的输出磁矩。

### 3. 文件管理技巧
在第 3 章中我们提到，为了避免频繁重命名 `STRU` 文件，推荐在 `INPUT` 文件中使用 `stru_file` 参数：

```ini
stru_file  ./STRU-001
```
这允许你直接指向 Phonopy 生成的特定微扰结构文件，而无需将其重命名为 `STRU`。这在编写批量提交脚本时非常有用，能有效降低文件操作错误的风险。