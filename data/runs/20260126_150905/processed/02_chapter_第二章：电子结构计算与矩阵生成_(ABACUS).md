# 第二章：电子结构计算与矩阵生成 (ABACUS)

在第一章中，我们介绍了线性光学性质计算的基本物理图像。从本章开始，我们将进入实战环节。

计算材料的光学性质（如介电函数、吸收系数）是一个典型的**两步走**流程：
1.  **第一步（ABACUS）**：进行基态 SCF 计算。这一步不仅仅是为了得到能量，更核心的目标是**“生产数据”**——即生成描述电子体系的哈密顿量矩阵（Hamiltonian）、重叠矩阵（Overlap）和位置算符矩阵（Position Matrix）。
2.  **第二步（pyatb）**：读取上述矩阵，利用 Kubo-Greenwood 公式进行后处理，计算光学电导率和介电函数。

本章将聚焦于**第一步**。我们将以晶态二氧化硅 ($\text{SiO}_2$) 为例，详细讲解如何配置 ABACUS 的 `INPUT` 文件以正确输出稀疏矩阵，并教会你如何从输出日志中提取连接两步计算的关键“桥梁”参数。

---

## 2.1 输入文件准备与关键参数设置

要使用 pyatb 进行后处理，ABACUS 的计算必须基于 **LCAO（线性组合原子轨道）** 基组。平面波（Plane Wave）模式下的波函数无法直接被 pyatb 读取。

### 2.1.1 文件结构准备
建议为每个案例建立清晰的文件夹结构。对于 $\text{SiO}_2$，我们需要准备以下核心文件：

```text
work_dir/
├── INPUT               # 控制参数
├── STRU                # 晶体结构 (SiO2)
├── KPT                 # K点设置
├── Si_ONCV_PBE-1.0.upf # Si 赝势
├── O_ONCV_PBE-1.0.upf  # O 赝势
├── Si_gga_7au_100Ry_2s2p1d.orb # Si 轨道文件
└── O_gga_7au_100Ry_2s2p1d.orb  # O 轨道文件
```

### 2.1.2 构建 INPUT 文件 (核心)

这是本章的重中之重。为了让 ABACUS 输出 pyatb 所需的 `.csr` (Compressed Sparse Row) 格式矩阵文件，必须显式开启特定的输出开关。

以下是针对 $\text{SiO}_2$ 光学性质计算的推荐 `INPUT` 文件设置：

```bash
INPUT_PARAMETERS
# 1. 基础计算控制
suffix              silica      # 输出文件夹名称，如 OUT.silica
calculation         scf         # 自洽场计算
esolver_type        ksdft       # Kohn-Sham DFT
basis_type          lcao        # 【关键】必须使用 LCAO 基组
symmetry            0           # 关闭对称性（推荐），避免矩阵输出时的旋转问题
ks_solver           genelpa     # 求解器，genelpa 效率较高

# 2. 电子结构与收敛参数
ecutwfc             100         # 能量截断 (Ry)，需根据赝势和轨道精度测试
scf_thr             1e-8        # 自洽收敛精度 (Ry)，光学计算建议高精度
smearing_method     gaussian    # 展宽方法
smearing_sigma      0.01        # 展宽宽度 (Ry)
mixing_type         broyden     # 电荷密度混合方法
mixing_beta         0.1         # 混合参数

# 3. 矩阵输出控制 (Pyatb 接口核心)
out_chg             1           # 输出电荷密度 (建议开启)
out_mat_hs2         1           # 【关键】输出二中心哈密顿量(H)和重叠矩阵(S)
out_mat_r           1           # 【关键】输出位置算符/偶极矩阵(r)

# 4. 赝势与轨道路径
pseudo_dir          ./
orbital_dir         ./
```

#### 💡 专家解读：为什么需要这些参数？

*   **`basis_type lcao`**: pyatb 的核心算法基于紧束缚模型（Tight-Binding），依赖于原子轨道基组下的矩阵元。
*   **`out_mat_hs2 1`**: 默认情况下 ABACUS 不输出稀疏矩阵。开启此项后，软件会将实空间中的 $H_{\mu\nu}(\mathbf{R})$ 和 $S_{\mu\nu}(\mathbf{R})$ 写入磁盘。这是计算能带结构的基础。
*   **`out_mat_r 1`**: 这是**光学性质计算的灵魂**。根据 Kubo-Greenwood 公式，光电导率取决于速度矩阵元 $\langle \psi_{n\mathbf{k}} | \mathbf{v} | \psi_{m\mathbf{k}} \rangle$。在 LCAO 表象下，这需要通过位置算符矩阵 $\langle \phi_\mu | \mathbf{r} | \phi_\nu \rangle$ 转换得到。**如果忘记开启此项，后续 pyatb 将因缺少偶极矩阵而无法计算介电函数。**

### 2.1.3 STRU 与 KPT 文件简述
*   **STRU**: 包含晶格常数、晶格矢量、原子质量、赝势/轨道文件名及原子坐标。请确保坐标单位与 `LATTICE_CONSTANT` 一致。
*   **KPT**: SCF 计算的 K 点密度。例如 `6 6 6 0 0 0`。
    *   *注意*：这里 K 点主要用于收敛基态电荷密度。pyatb 后处理时会使用更密的 K 点网格（在 pyatb 的输入中定义），但基态波函数的质量取决于这一步的 K 点设置。

---

## 2.2 运行 SCF 与关键信息提取

配置好文件后，即可运行 ABACUS。

### 2.2.1 提交任务
使用 MPI 并行运行（假设使用 16 核）：
```bash
export OMP_NUM_THREADS=1
mpirun -np 16 abacus
```

### 2.2.2 验证输出矩阵
计算完成后，首先检查输出目录（例如 `OUT.silica/`）下是否生成了以下关键文件：

1.  **`data-HR-sparse_SPIN0.csr`**: 实空间哈密顿量矩阵（稀疏格式）。
2.  **`data-SR-sparse_SPIN0.csr`**: 实空间重叠矩阵。
3.  **`data-rR-sparse.csr`**: 位置算符矩阵。

> **风险提示**：如果缺少 `data-rR-sparse.csr`，请检查 `INPUT` 中是否设置了 `out_mat_r 1`。如果缺少 `HR` 或 `SR` 文件，检查 `out_mat_hs2 1`。

### 2.2.3 提取“桥梁”参数 (Critical Step)

在进入下一章使用 pyatb 之前，我们需要从 ABACUS 的日志文件（`running_scf.log` 或标准输出）中手动提取两个至关重要的物理量。**pyatb 无法自动读取这些值，必须由用户填入其输入文件。**

我们需要提取：
1.  **Fermi Energy ($E_F$)**: 费米能级。
2.  **Occupied Bands**: 占据能带的数目（即价带数）。

#### 操作指南
进入输出目录，使用 `grep` 命令查找：

```bash
# 进入输出目录
cd OUT.silica

# 1. 查找占据带数 (Occupied bands)
grep "occupied bands" running_scf.log
# 输出示例: occupied bands = 64

# 2. 查找费米能级 (Fermi Energy)
grep "E_Fermi" running_scf.log | tail -n 1
# 输出示例: E_Fermi        0.4070749075         5.5385382545
```

#### 数据解读与记录
*   **Occupied Bands**: 示例中为 `64`。请记录此数值。
*   **E_Fermi**: ABACUS 通常输出两列数值。
    *   第一列单位通常为 **Rydberg (Ry)**。
    *   第二列单位通常为 **eV**。
    *   *实战建议*：记录 **eV** 单位的数值（例如 `5.5385382545`），因为 pyatb 默认支持 eV 输入，这样更符合直觉。

### 2.2.4 准备 pyatb 所需的晶格信息

除了上述电子参数，pyatb 的 `Input` 文件还需要手动填写晶格常数和晶格矢量。**pyatb 目前不会自动从 `.csr` 文件中解析晶格信息。**

请打开你的 `STRU` 文件或查看 `running_scf.log` 开头的 `Lattice vectors` 部分，记录下：
1.  **Lattice Constant**: 晶格常数（Bohr）。
2.  **Lattice Vectors**: 三个晶格矢量（通常是归一化的或以晶格常数为单位的）。

---

## 本章小结

至此，我们已经完成了“数据生产”阶段。你现在应该手握以下“原料”：
1.  **三个稀疏矩阵文件** (`HR`, `SR`, `rR`) —— 位于 `OUT.silica/` 目录下。
2.  **两个关键参数** ($E_F \approx 5.5385 \text{ eV}$, $N_{occ} = 64$) —— 记录在案。
3.  **晶格几何信息** —— 来源于 `STRU`。

在下一章中，我们将扮演“工匠”的角色，利用 pyatb 将这些原料加工成精美的光学性质图谱。请务必保存好上述文件，不要删除 `OUT.silica` 目录。