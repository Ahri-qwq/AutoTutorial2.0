基于您提供的核心案例（ABACUS+pyatb能带反折叠计算.md）及相关知识库，以下是关于**能带反折叠流程**的结构化元数据报告。

## 1. 物理本质 (Physics Concepts)
- **核心物理概念**: **能带反折叠 (Band Unfolding)**。
- **解决的科学问题**: 
    - 当研究掺杂、缺陷或合金体系时，通常需要构建**超胞 (Supercell)**。
    - 超胞的使用会导致布里渊区 (BZ) 变小，能带发生折叠，使得能带图变得极其复杂且难以分析。
    - 实验手段（如 **ARPES**，角分辨光电子能谱）通常观测到的是基于**原胞 (Primitive Cell)** 对称性的能带特征。
    - 能带反折叠通过将超胞的波函数投影回原胞的布里渊区，恢复能带的“有效”色散关系和谱权重，从而实现理论计算与实验结果的直接对比。

## 2. 关键输入参数 (Key Parameters)

### 2.1 ABACUS 自洽计算 (INPUT 文件)
此阶段目的是生成紧束缚哈密顿量矩阵。

| 参数名 | 推荐值 | 物理意义 |
| :--- | :--- | :--- |
| `calculation` | `scf` | 进行自洽场计算。 |
| `basis_type` | `lcao` | **必须设置**。PYATB 基于数值原子轨道 (LCAO) 产生的稀疏矩阵进行后处理。 |
| `out_mat_hs2` | `1` | **核心参数**。输出稀疏格式的哈密顿量 (H) 和重叠 (S) 矩阵（生成 `data-HR-sparse_SPINx.csr` 和 `data-SR-sparse_SPINx.csr`）。 |
| `out_mat_r` | `1` | **核心参数**。输出位置算符矩阵（生成 `data-rR-sparse.csr`）。 |
| `nbands` | (根据体系设置) | 需要包含足够多的空带，以覆盖反折叠分析所需的能量范围。 |

### 2.2 PYATB 后处理 (Input 文件)
注意：PYATB 的输入文件通常命名为 `Input` (区分大小写)，以区别于 ABACUS 的 `INPUT`。

| 参数名 | 推荐值/示例 | 物理意义 |
| :--- | :--- | :--- |
| `package` | `ABACUS` | 指定接口软件为 ABACUS。 |
| `fermi_energy` | (手动填入) | **关键手动步骤**。需从 ABACUS 的 `running_scf.log` 中读取 `EFERMI` 的值并填入此处。 |
| `HR_route` | `data-HR-sparse_SPIN0.csr` | 指定 ABACUS 生成的哈密顿量矩阵文件路径。 |
| `SR_route` | `data-SR-sparse_SPIN0.csr` | 指定 ABACUS 生成的重叠矩阵文件路径。 |
| `rR_route` | `data-rR-sparse.csr` | 指定 ABACUS 生成的位置算符矩阵文件路径。 |
| `m_matrix` | (整数矩阵，如 2 0 0 0 2 0 0 0 2) | **核心参数**。描述超胞相对于原胞的扩胞矩阵 (Transformation Matrix)。例如 2x2x2 扩胞则对角线为 2。 |
| `band_range` | `1 500` | 指定计算/绘图的能带索引范围。 |
| `stru_file` | `STRU` | 指定结构文件（需与 ABACUS 计算时一致）。 |

## 3. 体系与接口配置 (System & Interfaces)

- **结构 (STRU)**: 
    - ABACUS 计算使用的是**超胞 (Supercell)** 结构。
    - PYATB 需要读取该 `STRU` 文件以及对应的轨道文件 (`.orb`) 和赝势文件 (`.upf`)。
- **外部接口 (PYATB)**:
    - **依赖库**: 需要预先安装 `pybind11`, `mpi4py`, `eigen3`。
    - **数据流**: ABACUS (生成 `.csr` 矩阵) -> 用户 (复制文件 & 填写 Fermi 能) -> PYATB (读取矩阵与结构，输出能带数据)。
- **接口注意事项**:
    - **文件命名**: ABACUS 的输入通常为 `INPUT`，PYATB 的输入通常为 `Input`，两者容易混淆，建议分文件夹管理。
    - **矩阵文件**: ABACUS 输出的 `.csr` 文件位于 `OUT.suffix/` 目录下，必须手动复制到 PYATB 的工作目录。

## 4. 教程编写特殊指令 (Special Instructions for Writer)

- **Critical (流程分离)**: 教程必须明确建议用户建立两个文件夹（例如 `abacus_run` 和 `pyatb_run`）。
    - 原因：避免输入文件 (`INPUT` vs `Input`) 混淆，保持目录整洁。
- **Critical (手动操作)**: 必须高亮强调 **“获取费米能级”** 这一步。
    - 这是一个非自动化的断点：用户必须 `grep EFERMI` 从 ABACUS 日志中获取数值，然后**手动**写入 PYATB 的 `Input` 文件。如果这一步数值错误，最终能带图的能量零点将是错误的。
- **Critical (扩胞矩阵)**: 在解释 `m_matrix` 时，需说明这是超胞晶格矢量与原胞晶格矢量的变换关系。如果用户做的是 3x3x1 的扩胞，矩阵应对应修改，不能照抄案例中的 2x2x2。
- **Plotting (绘图脚本)**: 提醒用户 PYATB 生成的 `plot_unfold.py` 脚本中，`energy_range` 参数通常需要手动修改，以匹配实际关心的能量窗口（如案例中从 `[-4, 6]` 改为 `[-14, 12]`）。
- **环境配置**: 教程开头必须警告镜像选择（推荐 `abacus:3.2.3`），因为涉及编译安装 PYATB，环境依赖极易出错。

## 5. 常见报错与注意事项 (Pitfalls)

- **MKL 库链接错误**: 
    - **现象**: 运行 PYATB 时报错 `undefined symbol: mkl_sparse_optimize_bsr_trsm_i8`。
    - **解决**: 必须在运行前设置 `LD_PRELOAD` 环境变量，加载 Intel MKL 库（参考案例中的 `os.environ['LD_PRELOAD'] = ...` 代码块）。
- **基组选择错误**:
    - 如果 ABACUS `INPUT` 中设置 `basis_type pw` (平面波)，将无法生成 PYATB 所需的稀疏矩阵，流程无法进行。必须是 `lcao`。
- **文件缺失**:
    - 忘记将 `data-HR-sparse...` 等三个 `.csr` 文件从 ABACUS 输出目录复制到 PYATB 目录，会导致 PYATB 找不到输入文件。
- **自旋设置**:
    - 案例中是 `nspin 1`。如果是自旋极化计算 (`nspin 2`)，ABACUS 会输出 `SPIN0` (up) 和 `SPIN1` (down) 的矩阵，PYATB 的 `Input` 文件中 `HR_route` 等参数可能需要对应调整或分别计算（需查阅 PYATB 文档确认多自旋处理方式，本案例仅展示了单自旋）。