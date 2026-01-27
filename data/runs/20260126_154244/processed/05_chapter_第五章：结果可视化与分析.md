# 第五章：结果可视化与分析

在前几章中，我们已经成功完成了 ABACUS 的自洽计算（SCF），并利用 PYATB 完成了能带反折叠（Band Unfolding）的核心计算步骤。此时，你的 `pyatb` 工作目录下应该已经生成了 `spectral_weight.dat` 和 `plot_unfold.py` 等关键文件。

本章将指导你如何将这些抽象的数值数据，转化为直观的“有效能带图”（Effective Band Structure），并学会从物理角度解读这些结果。

> **⚠️ 环境与依赖警告**
> 本章涉及 Python 绘图脚本的运行。如果你使用的是容器环境，强烈推荐使用 **`abacus:3.2.3`** 镜像，该镜像已预装了 PYATB 所需的 `matplotlib`、`numpy` 等依赖库。其他环境可能需要手动安装依赖。

---

## 5.1 绘图脚本参数调整

PYATB 在运行结束后，会自动生成一个名为 `plot_unfold.py` 的 Python 脚本。这是一个基于 Matplotlib 的绘图模板。直接运行它通常能得到一张图，但为了获得最佳的视觉效果，我们需要根据物理体系的实际情况调整能量窗口。

### 5.1.1 理解 `energy_range`

打开 `plot_unfold.py` 文件，找到如下参数设置：

```python
# plot_unfold.py 片段
# ...
energy_range = [-4, 6]  # 默认值可能不适合所有体系
# ...
```

这里的 `energy_range` 定义了 Y 轴（能量轴）的显示范围，单位通常为 eV。
*   **零点基准**：这里的 $E=0$ 对应于你在 `Input` 文件中手动填入的 `fermi_energy`。
*   **调整策略**：
    *   对于**绝缘体/半导体**：通常关注价带顶（VBM）和导带底（CBM）附近的区域，例如 `[-4, 4]` 或 `[-6, 6]`。
    *   对于**金属**：关注费米面附近的电子态，范围可以更小。
    *   **案例修正**：在 GeC 案例中，为了完整展示深层能级和高能激发态，我们将范围扩大。

### 5.1.2 修改实战

你可以使用文本编辑器（如 `vim` 或 `nano`）手动修改，也可以使用 `sed` 命令快速替换。

**方式一：手动修改**
将 `energy_range` 的数值改为你关心的区间，例如 `[-14, 12]`：
```python
energy_range = [-14, 12]
```

**方式二：命令行流编辑（推荐脚本作业使用）**
```bash
# 将默认的 [-4, 6] 替换为 [-14, 12]
sed -i 's/energy_range = \[-4, 6\]/energy_range = [-14, 12]/g' plot_unfold.py
```

### 5.1.3 运行绘图

修改完成后，执行脚本生成 PDF 图像：

```bash
python plot_unfold.py
```

如果运行成功，目录下会生成 `unfold.pdf`。

---

## 5.2 结果解读

打开生成的 `unfold.pdf`，你会看到一张类似能带图的图像，但它包含了比普通能带图更丰富的信息。

### 5.2.1 谱权重（Spectral Weight）的物理意义

在反折叠能带图中，线条不再是单一的黑色实线，而是具有**颜色深浅**或**粗细**变化的“云图”。

*   **颜色/深浅**：代表 **谱权重 (Spectral Weight)** 的大小。
    *   $$P_{Km}(k_i)$$：表示超胞中的本征态 $|K, m\rangle$ 投影到原胞布洛赫态 $|k_i\rangle$ 上的概率。
*   **深色/高亮线条**：代表谱权重接近 1。这意味该能带保持了良好的原胞平移对称性，电子态是离域的（Bloch-like）。
*   **浅色/模糊线条（影子带）**：代表谱权重较小。这些通常是由于超胞引入的微扰（如缺陷、替位原子）导致对称性破缺，使得原胞的 $k$ 点与 $k+G$ 点发生了耦合，产生的“折叠”假象。
*   **能带展宽（Smearing）**：如果在某些能量区间线条变得模糊不清，说明该能量处的电子态受到了强烈的散射（例如无序合金或高浓度掺杂），动量 $k$ 不再是一个好的量子数。

### 5.2.2 典型特征识别

1.  **主能带 (Primitive Bands)**：图中清晰、连续且颜色最深的轨迹。这对应于未受扰动的理想晶体的能带结构。
2.  **能隙 (Band Gap)**：观察费米能级（$E=0$）附近。如果 $E=0$ 处没有谱权重分布，且上下有清晰的带边，则为半导体/绝缘体。
3.  **缺陷态 (Defect States)**：如果在禁带中间出现了平直的、色散很小（flat band）的线条，这往往对应于局域的缺陷能级。

---

## 附录：常见问题与进阶建议

在 ABACUS + PYATB 的实战流程中，以下环节是“深坑”高发区，请务必仔细核对。

### A.1 目录结构与文件名陷阱 (Critical)

**强烈建议**将 ABACUS 自洽计算与 PYATB 后处理分在两个不同文件夹进行！

*   **原因**：ABACUS 的输入文件名为 `INPUT` (全大写)，而 PYATB 的输入文件名为 `Input` (首字母大写)。在 Linux 系统下虽然区分大小写，但极易在 Tab 补全或脚本编写时混淆。此外，保持目录整洁有助于管理大量输出文件。
*   **推荐结构**：
    ```text
    Project_Root/
    ├── abacus_run/       # 运行 SCF
    │   ├── INPUT
    │   ├── STRU
    │   ├── KPT
    │   └── OUT.suffix/   # ABACUS 输出目录
    └── pyatb_run/        # 运行 Unfolding
        ├── Input         # PYATB 输入
        ├── STRU          # 复制自 abacus_run
        ├── data-HR-sparse_SPIN0.csr  # 复制自 OUT.suffix
        ├── data-SR-sparse_SPIN0.csr  # 复制自 OUT.suffix
        └── data-rR-sparse.csr        # 复制自 OUT.suffix
    ```

### A.2 费米能级的获取 (Manual Checkpoint)

这是全流程中唯一无法自动化的**断点**，必须手动操作。

1.  **获取数值**：在 `abacus_run` 目录下，从 SCF 日志中提取费米能级。
    ```bash
    grep "EFERMI" OUT.suffix/running_scf.log
    # 输出示例: EFERMI = 11.08291763890926 eV
    ```
2.  **填入 Input**：将上述数值**精确**填入 `pyatb_run/Input` 文件中：
    ```json
    INPUT_PARAMETERS
    {
        ...
        fermi_energy    11.08291763890926
        fermi_energy_unit eV
        ...
    }
    ```
    > **警告**：如果此数值填写错误（例如填了 0），最终画出的能带图 $E=0$ 位置将不是费米面，导致物理分析完全错误。

### A.3 扩胞矩阵 (m_matrix) 的正确设置

在 `Input` 文件中，`m_matrix` 参数定义了**超胞晶格矢量**与**原胞晶格矢量**之间的变换关系。

$$ \vec{A}_{super} = \vec{A}_{prim} \times M $$

*   **不要盲目照抄案例**：案例中的矩阵（如 `-2 2 2 ...`）是针对特定晶体结构（如 FCC 原胞转立方超胞）的。
*   **简单扩胞**：如果你只是简单地将原胞沿 $x, y, z$ 方向分别扩大 $N_x, N_y, N_z$ 倍（例如 `3x3x1` 超胞），则矩阵应为对角阵：
    ```text
    m_matrix  3 0 0  0 3 0  0 0 1
    ```
    如果是 `2x2x2` 扩胞，则是：
    ```text
    m_matrix  2 0 0  0 2 0  0 0 2
    ```

### A.4 磁性体系的处理 (Spin Polarization)

如果你的体系是磁性的（在 ABACUS `INPUT` 中设置了 `nspin 2`），ABACUS 会输出两套哈密顿量文件：
*   `data-HR-sparse_SPIN0.csr` (自旋向上 / Up)
*   `data-HR-sparse_SPIN1.csr` (自旋向下 / Down)

**处理方法**：
你需要分别运行两次 PYATB（建议分两个子文件夹），分别在 `Input` 文件中指定对应的 `HR_route` 和 `SR_route`，最后分别画出 Spin Up 和 Spin Down 的能带图进行对比。

### A.5 常见报错：MKL 链接错误

如果在运行 `pyatb` 命令时遇到类似 `undefined symbol: mkl_sparse_optimize_bsr_trsm_i8` 的错误，这是因为环境变量未正确加载 Intel MKL 库。

**解决方案**：
在运行前执行以下命令（可视情况加入 `.bashrc`）：
```bash
export LD_PRELOAD='/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so'
```
*(注：具体路径可能随镜像版本不同而微调，上述路径适用于 `abacus:3.2.3` 镜像)*