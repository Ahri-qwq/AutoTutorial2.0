# 第四章：pyatb 后处理与能带反折叠计算

在上一章中，我们已经成功完成了超胞（Supercell）的自洽计算（SCF），并生成了构建紧束缚模型所需的稀疏矩阵文件。然而，超胞的能带图由于布里渊区折叠（Band Folding），通常杂乱无章，无法直接与实验（如 ARPES）或原胞能带进行对比。

本章将是整个流程的“核心枢纽”。我们将使用 `pyatb` 软件，将 ABACUS 计算得到的电子结构信息“反折叠”（Unfold）回原胞的布里渊区。

> **⚠️ 环境预警**：
> 本章涉及 `pyatb` 的编译与运行，强烈建议使用官方推荐镜像 **`abacus:3.2.3`**。其他环境极易出现 `pybind11` 或 `MKL` 库链接错误。

---

## 4.1 工作目录规划与数据迁移

在进行后处理之前，**强烈建议**将 ABACUS 的计算目录与 `pyatb` 的运行目录物理隔离。这不仅是为了保持整洁，更是为了防止输入文件（ABACUS 的 `INPUT` 与 pyatb 的 `Input`）发生混淆。

### 4.1.1 建立独立目录
假设你的 ABACUS 自洽计算目录为 `abacus_run`，请在同级目录下新建一个名为 `pyatb_run` 的文件夹。

```bash
# 假设当前在项目根目录
mkdir pyatb_run
```

### 4.1.2 关键文件迁移
我们需要将 ABACUS 计算生成的**稀疏矩阵文件**和**结构文件**复制到新目录。

需要复制的文件清单如下：
1.  **`STRU`**: 结构文件（pyatb 需要读取晶格信息）。
2.  **`data-HR-sparse_SPIN0.csr`**: 哈密顿量矩阵（若开启自旋可能有 SPIN1）。
3.  **`data-SR-sparse_SPIN0.csr`**: 重叠矩阵。
4.  **`data-rR-sparse.csr`**: 位置算符偶极矩阵。

> **注意**：上述 `.csr` 文件通常位于 ABACUS 输出目录（如 `OUT.suffix/`）下。

**操作示例**：
```bash
# 进入 ABACUS 输出目录
cd abacus_run/OUT.GeC/

# 复制文件到 pyatb 目录
cp data-HR-sparse_SPIN0.csr ../../pyatb_run/
cp data-SR-sparse_SPIN0.csr ../../pyatb_run/
cp data-rR-sparse.csr       ../../pyatb_run/

# 别忘了复制 STRU 文件（通常在 abacus_run 根目录）
cd ../
cp STRU ../pyatb_run/
```

完成后的 `pyatb_run` 目录结构应如下所示：
```text
pyatb_run/
├── STRU
├── data-HR-sparse_SPIN0.csr
├── data-SR-sparse_SPIN0.csr
└── data-rR-sparse.csr
```

---

## 4.2 关键步骤：费米能级的手动提取

这是流程中**唯一无法自动化且至关重要**的一步。

ABACUS 计算得到的费米能级（Fermi Energy）存储在日志文件中，而 `pyatb` 无法自动读取该值。如果这一步数值填错，最终能带图的能量零点（$E-E_f$）将完全错误。

### 操作方法
请回到 ABACUS 的运行目录，从 `running_scf.log`（或屏幕输出日志）中提取 `EFERMI`。

```bash
# 在 ABACUS 输出目录下执行
grep "EFERMI" running_scf.log
```

**输出示例**：
```text
EFERMI = 11.08291763890926 eV
```

> **📝 记录下来**：请复制这个数值（例如 `11.08291763890926`），我们将在下一节的输入文件中使用它。

---

## 4.3 pyatb 输入文件 (Input) 详解

在 `pyatb_run` 目录下，我们需要创建一个名为 `Input` 的文件（注意首字母大写，区分于 ABACUS 的 `INPUT`）。

### 4.3.1 Input 文件完整示例

```json
INPUT_PARAMETERS
{
    "nspin": 1,
    "package": "ABACUS",
    "fermi_energy": 11.08291763890926,  
    "fermi_energy_unit": "eV",
    "HR_route": "data-HR-sparse_SPIN0.csr",
    "SR_route": "data-SR-sparse_SPIN0.csr",
    "rR_route": "data-rR-sparse.csr",
    "HR_unit": "Ry",
    "rR_unit": "Bohr"
}

LATTICE
{
    "lattice_constant": 1.889726,
    "lattice_constant_unit": "Bohr",
    "lattice_vector": [
        11.3149995804, 0, 0,
        0, 11.3149995804, 0,
        0, 0, 11.3149995804
    ]
}

BANDUNFOLDING
{
    "stru_file": "STRU",
    "ecut": 50,
    "band_range": [1, 500],
    "m_matrix": [
        -2, 2, 2,
        2, -2, 2,
        2, 2, -2
    ],
    "kpoint_mode": "line",
    "kpoint_num": 50,
    "high_symmetry_kpoint": [
        0.5, 0.5, 0.5, 50,  # L
        0.0, 0.0, 0.0, 50,  # Gamma
        0.5, 0.0, 0.5, 50,  # X
        0.0, 0.0, 0.0, 1    # Gamma
    ]
}
```

### 4.3.2 核心参数深度解析

#### 1. 基础参数 (`INPUT_PARAMETERS`)
*   **`package`**: 必须指定为 `ABACUS`。
*   **`fermi_energy`**: **此处填入 4.2 节中提取的数值**。
*   **`*_route`**: 指定矩阵文件的路径。由于我们已经将文件复制到了当前目录，直接写文件名即可。

#### 2. 晶格信息 (`LATTICE`)
*   建议直接从 ABACUS 的 `STRU` 文件或 `running_scf.log` 中复制晶格常数和基矢量。注意单位通常为 Bohr。

#### 3. 扩胞矩阵 (`m_matrix`) —— **易错点**
这是能带反折叠的物理核心。`m_matrix` 定义了**超胞基矢量 ($A_{super}$)** 与 **原胞基矢量 ($A_{prim}$)** 之间的变换关系。
关系式为：
$$ A_{super} = A_{prim} \times M $$

*   **不要照抄案例**：案例中的 `[-2, 2, 2, ...]` 是针对特定晶体结构的 $2\times2\times2$ 扩胞变换。
*   **如何设置**：
    *   如果你做的是简单的 $2\times2\times2$ 正交扩胞，矩阵通常是 `2 0 0 0 2 0 0 0 2`（即对角矩阵）。
    *   如果你使用了 `atomkit` 或其他工具生成了旋转扩胞，必须使用生成超胞时对应的变换矩阵。
    *   **格式**：在 Input 文件中按行优先顺序填入 9 个数字。

#### 4. 能带范围 (`band_range`)
*   指定需要计算并反折叠的能带索引范围（如 `[1, 500]`）。确保该范围覆盖了费米面附近的价带和导带。

---

## 4.4 运行 pyatb 与 MKL 环境配置

配置好 `Input` 文件后，就可以运行计算了。但在执行命令前，必须解决一个常见的 Intel MKL 库链接问题。

### 4.4.1 解决 MKL `undefined symbol` 报错
在 `abacus:3.2.3` 容器或大多数 Linux 环境中，直接运行可能会遇到如下报错：
`undefined symbol: mkl_sparse_optimize_bsr_trsm_i8`

这是由于系统预加载库顺序导致的。**解决方法**是在运行前设置 `LD_PRELOAD` 环境变量。

请在终端执行以下命令（或将其写入提交脚本）：

```bash
export LD_PRELOAD=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so
```

> **提示**：如果你使用的是不同的镜像或自定义环境，上述路径可能需要根据实际 MKL 安装位置进行调整。

### 4.4.2 执行计算
设置好环境变量后，使用 MPI 并行运行 pyatb：

```bash
# -n 指定核数，建议与物理核数一致
mpirun -n 16 pyatb
```

计算通常很快（几分钟内）。运行结束后，当前目录下会生成一个 `Out` 文件夹。

---

## 4.5 结果可视化

进入输出目录查看结果：

```bash
cd Out/Bandunfolding
ls
```

你会看到以下关键文件：
*   `spectral_weight.dat`: 包含反折叠后的能带权重数据。
*   `plot_unfold.py`: 自动生成的绘图脚本。

### 4.5.1 调整绘图范围
默认生成的 `plot_unfold.py` 脚本中，能量显示范围 (`energy_range`) 可能不符合你的需求（例如默认可能是 `[-4, 6]`）。

**必须手动修改该脚本**：
1.  打开 `plot_unfold.py`。
2.  找到 `energy_range` 变量。
3.  根据你的体系修改数值。例如，为了看清深层能级，可以改为：
    ```python
    # 修改前
    # energy_range = [-4, 6]
    
    # 修改后
    energy_range = [-14, 12]
    ```

### 4.5.2 生成能带图
运行修改后的脚本：

```bash
python plot_unfold.py
```

成功运行后，将生成 `unfold.pdf`。这张图展示了超胞能带投影回原胞布里渊区后的“有效能带结构”（Effective Band Structure），其中颜色的深浅代表了谱权重（Spectral Weight），清晰地还原了被杂质或缺陷破坏的能带色散关系。