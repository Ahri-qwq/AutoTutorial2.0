# 第三章：ABACUS 自洽计算与矩阵生成

在上一章中，我们完成了能带反折叠的理论准备。本章我们将进入实战的核心环节：使用 ABACUS 进行自洽场（SCF）计算。

在此流程中，ABACUS 的角色不仅仅是计算能量，更是一个**“哈密顿量生成器”**。我们的目标是通过 LCAO（原子轨道线性组合）基组，构建并输出体系的稀疏矩阵（Hamiltonian, Overlap, Position），这些矩阵是后续 PYATB 进行能带反折叠分析的原材料。

> **⚠️ 环境配置警告**
> 本教程涉及编译安装 PYATB，环境依赖极其敏感。强烈建议使用官方推荐的镜像 **`abacus:3.2.3`**。使用其他版本（特别是较旧版本或不兼容的 GCC 环境）极易导致后续 PYATB 安装失败或运行时出现 `undefined symbol: mkl_...` 等数学库错误。

---

## 3.0 工作目录规划（关键）

在开始计算前，必须强调一个工程习惯：**流程分离**。

由于 ABACUS 的输入文件名为 `INPUT`，而 PYATB 的输入文件通常命名为 `Input`（或 `input.json`），在同一目录下操作极易导致文件覆盖或混淆。此外，保持计算目录（生成数据）与分析目录（处理数据）的分离，有助于数据的长期管理。

**建议的文件结构：**

```text
Project_Root/
├── abacus_run/          <-- 本章工作目录
│   ├── INPUT            # ABACUS 控制参数
│   ├── STRU             # 超胞结构文件
│   ├── KPT              # K点设置
│   ├── *.upf & *.orb    # 赝势与轨道文件
│   └── OUT.suffix/      # 计算输出目录
└── pyatb_run/           <-- 下一章工作目录
    ├── Input            # PYATB 控制参数
    └── (待复制的矩阵文件)
```

请立即建立这两个文件夹，并进入 `abacus_run` 目录开始本章操作。

---

## 3.1 输入参数设置 (INPUT)

`INPUT` 文件控制着 ABACUS 的计算行为。为了适配 PYATB，我们需要开启特定的矩阵输出开关。

### 核心参数详解

1.  **`basis_type`: `lcao`**
    *   **必须设置**。能带反折叠基于紧束缚模型（Tight-Binding），必须使用原子轨道基组。平面波（`pw`）模式无法生成所需的轨道矩阵。

2.  **`calculation`: `scf`**
    *   进行自洽计算以获得收敛的电子密度和波函数。

3.  **`out_mat_hs2`: `1`**
    *   **关键开关**。设置为 `1` 时，ABACUS 会在计算结束后输出稀疏格式的哈密顿量矩阵（$H$）和重叠矩阵（$S$）。这是 PYATB 读取电子结构的核心数据源。

4.  **`out_mat_r`: `1`**
    *   **关键开关**。输出位置算符矩阵（$r$）。这对于计算偶极矩、光学性质或处理某些涉及贝里相位的计算是必要的。

5.  **`nbands`**
    *   **需特别注意**。默认情况下，ABACUS 仅计算被占据的能带（加上少量空带）。但在能带反折叠分析中，我们通常也关心导带（Conduction Bands）的结构。
    *   **建议**：显式设置一个较大的 `nbands` 值（例如总电子数的一半 + 20~50，或者更多），确保覆盖你感兴趣的高能级范围。

### INPUT 文件示例

新建 `INPUT` 文件并写入以下内容（以 GeC 体系为例，请根据实际体系调整）：

```bash
INPUT_PARAMETERS
# System variables
suffix                GeC_scf
ntype                 2
calculation           scf
symmetry              1         # 开启对称性分析可加速计算
init_chg              atomic

# Electronic structure
basis_type            lcao      # 【必须】使用 LCAO 基组
ks_solver             genelpa   # 推荐使用的对角化求解器
nspin                 1         # 1为非磁性，2为自旋极化
ecutwfc               50        # 截断能 (Ry)
scf_nmax              100       # 最大自洽步数
scf_thr               1e-7      # 自洽收敛阈值
nbands                100       # 【重要】确保包含足够的空带

# Variables related to output information for PYATB
out_mat_hs2           1         # 【关键】输出 H 和 S 稀疏矩阵
out_mat_r             1         # 【关键】输出位置矩阵
out_chg               0         # 可选：是否输出电荷密度
out_band              0         # scf计算通常不需要输出能带文件
```

---

## 3.2 结构与基组文件准备

在此步骤中，我们需要准备超胞（Supercell）的结构文件。

### 1. 结构文件 (STRU)
**注意**：这里的 `STRU` 必须是**扩胞后**的结构（Large Cell），而不是原胞（Primitive Cell）。能带反折叠的物理本质，就是将这个超胞的能带“投影”回原胞的布里渊区。

*   **逻辑检查**：如果你计划研究 3x3x1 的超胞掺杂效应，这里的 `STRU` 就必须包含 3x3x1 倍数的原子。
*   **矩阵变换预警**：请记住你从原胞扩胞到超胞使用的变换矩阵（例如 2x2x2 或 3x3x1）。在下一章配置 PYATB 的 `m_matrix` 参数时，必须与此处的结构严格对应，**切勿照抄案例中的矩阵**。

### 2. 赝势与轨道文件 (.upf & .orb)
确保 `INPUT` 文件中的 `pseudo_dir` 和 `orbital_dir` 指向正确的路径。
*   `.upf`: 赝势文件。
*   `.orb`: 数值原子轨道文件（LCAO 专用）。

### 3. K点设置 (KPT)
对于超胞计算，由于实空间晶格变大，倒空间布里渊区变小，因此 K 点网格可以相应稀疏一些（例如原胞用 12x12x12，2x2x2 超胞可能只需要 6x6x6）。

```text
K_POINTS
0
Gamma
2 2 2 0 0 0
```

---

## 3.3 运行自洽计算与结果检查

准备好所有文件后，提交任务运行 ABACUS。

```bash
# 假设使用 16 核并行
mpirun -n 16 abacus
```

### 1. 检查运行状态
计算完成后，首先检查屏幕输出或日志文件，确保 `SCF convergence reached`，即自洽收敛。

### 2. 验证矩阵输出
进入输出目录（例如 `OUT.GeC_scf/`），检查是否生成了以下关键文件：

*   `data-HR-sparse_SPIN0.csr`: 哈密顿量矩阵（稀疏格式）
*   `data-SR-sparse_SPIN0.csr`: 重叠矩阵
*   `data-rR-sparse.csr`: 位置算符矩阵

如果缺少这些文件，请回头检查 `INPUT` 中 `out_mat_hs2` 和 `out_mat_r` 是否正确设置为 `1`。

### 3. 获取费米能级（关键手动步骤）

这是一个**非自动化**的断点，至关重要。PYATB 需要知道体系的费米能级位置以对齐能量零点，但它不会自动从矩阵文件中读取。

你需要查看 ABACUS 的日志文件（通常是 `running_scf.log` 或标准输出），找到 `EFERMI` 的值。

```bash
# 在输出目录执行 grep 命令
grep "EFERMI" running_scf.log
```

**输出示例：**
```text
EFERMI = 6.54321 eV
```

> **📝 记录下来**：请务必记录这个数值（例如 `6.54321`）。在下一章编写 PYATB 的 `Input` 文件时，你需要手动将这个值填入 `fermi_energy` 字段。如果数值错误，最终能带图的费米面将发生偏移。

---

## 3.4 为 PYATB 准备数据

为了保持目录整洁，我们将生成的核心文件转移到 `pyatb_run` 文件夹中，为下一章做准备。

```bash
# 假设当前在 abacus_run/OUT.GeC_scf/ 目录下
# 且 pyatb_run 文件夹位于上两级目录

cp data-HR-sparse_SPIN0.csr ../../pyatb_run/
cp data-SR-sparse_SPIN0.csr ../../pyatb_run/
cp data-rR-sparse.csr       ../../pyatb_run/

# 此外，PYATB 也需要读取结构文件和轨道/赝势信息
cp ../STRU                  ../../pyatb_run/
# 建议将赝势和轨道文件也软链接或复制过去，或者在 pyatb Input 中指定绝对路径
```

至此，我们已经成功利用 ABACUS 生成了包含体系全部量子力学信息的稀疏矩阵。下一章，我们将切换到 PYATB，利用这些矩阵“反推”出有效能带结构。