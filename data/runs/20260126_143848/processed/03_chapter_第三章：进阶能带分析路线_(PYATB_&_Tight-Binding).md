# 第三章：进阶能带分析路线 (PYATB & Tight-Binding)

在上一章中，我们通过传统的 **NSCF（非自洽场）** 路线，利用 `calculation = nscf` 和 Line 模式的 K 点，成功绘制了基础能带图。这种方法对于完美的晶体结构非常有效，能直接给出本征值。

然而，当我们面对更复杂的物理场景时，传统 NSCF 路线会显得力不从心：
1.  **轨道成分分析**：如果你想知道某条能带是由 $d_{xy}$ 轨道贡献还是 $p_z$ 轨道贡献（即 Fat Band 投影能带）。
2.  **大体系反折叠**：当你计算掺杂、缺陷或界面体系时，必须使用超胞（Supercell）。超胞会导致布里渊区折叠，产生极其密集且难以辨认的“意大利面”能带。此时需要 **能带反折叠 (Band Unfolding)** 技术将其还原回原胞的布里渊区。

为了解决这些问题，我们需要引入 **PYATB (Python Ab-initio Tight-Binding)** 路线。

> **核心区分 (Critical Distinction)**：
> *   **传统 NSCF 路线**：`SCF` -> `NSCF` (读取电荷密度) -> 获得本征值。
> *   **PYATB 紧束缚路线**：`SCF` (开启矩阵输出) -> `PYATB` (后处理) -> 获得能带/反折叠/投影。
>
> **注意**：本章的计算流程与第二章完全独立，请勿混淆。

---

## 3.1 构建紧束缚哈密顿量 (ABACUS 端)

PYATB 的工作原理是基于 ABACUS 计算出的**局域轨道（LCAO）哈密顿量矩阵**。因此，第一步是在 ABACUS 中进行一次 SCF 计算，并命令其将这些矩阵写出到磁盘。

### 3.1.1 准备输入文件

我们需要进行一次标准的自洽计算（`calculation = scf`），但在 `INPUT` 文件中需要增加特定的参数以输出矩阵。

**关键 `INPUT` 参数设置**：

```bash
INPUT_PARAMETERS
# ... 基础参数 (ecutwfc, scf_nmax 等) ...

# 1. 必须使用 LCAO 基组
basis_type      lcao

# 2. 计算模式为自洽计算
calculation     scf

# 3. 开启紧束缚矩阵输出 (关键步骤)
# 注意：不同版本的 ABACUS 参数名可能略有不同，请以官方文档 "Tight-Binding Matrices" 章节为准
out_mat_hs2     1        # 输出稀疏格式的哈密顿量(H)和重叠矩阵(S)
out_mat_r       1        # 输出位置算符偶极矩阵(rR)，用于计算 Berry Phase 或光学性质

# 4. 推荐设置
symmetry        1        # 生成矩阵时通常保持对称性，PYATB 会处理
cal_force       0        # 能带分析通常不需要算力
cal_stress      0        # 不需要算应力
```

### 3.1.2 运行计算与输出检查

运行 ABACUS 可执行文件：
```bash
mpirun -np 16 abacus | tee out.log
```

计算结束后，检查输出目录（默认为 `OUT.ABACUS/`），你必须找到以下关键文件，它们是连接 ABACUS 与 PYATB 的桥梁：

1.  **`data-HR-sparse_SPIN0.csr`** (或类似命名): 哈密顿量矩阵 $H(R)$。
2.  **`data-SR-sparse_SPIN0.csr`**: 重叠矩阵 $S(R)$。
3.  **`data-rR-sparse_SPIN0.csr`**: 偶极矩阵 $\vec{r}(R)$。
4.  **`STRU`**: 结构文件（PYATB 需要读取晶格信息）。
5.  **`*.orb`**: 原子轨道文件（**极度重要**，PYATB 需要读取轨道信息来构建基组，请确保该文件在当前目录或路径正确）。

> **风险提示**：如果未生成 `.csr` 文件，请检查 `out_mat_hs2` 参数是否被正确识别，或查阅你所使用的 ABACUS 版本的官方文档中关于 "Tight Binding Interface" 的说明。

---

## 3.2 投影能带 (Fat Band) 计算

得到矩阵后，我们脱离 ABACUS 主程序，转而使用 Python 脚本调用 PYATB 库进行后处理。Fat Band 可以直观地展示不同原子轨道对能带的贡献权重。

### 3.2.1 准备 PYATB 脚本

创建一个 Python 脚本（例如 `run_fatband.py`）。以下是核心逻辑的伪代码示例：

```python
from pyatb import TB_model, Solver

# 1. 初始化 Tight-Binding 模型
# 需要指定 ABACUS 的输出目录和轨道文件
model = TB_model.from_abacus(
    path='./OUT.ABACUS',           # ABACUS 输出目录
    binary=False,                  # 是否为二进制格式，csr 通常为文本
    orb_files=['Si_gga_8au_60Ry_2s2p1d.orb'] # 必须与 STRU 中使用的轨道文件一致
)

# 2. 设置求解器
solver = Solver(model)

# 3. 设置 K 点路径 (高对称点)
# 坐标需参考晶体结构，例如面心立方
kpath = [
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.0, 0.5],  # X
    # ... 添加更多点
]
labels = ['G', 'X']
solver.set_kpath(kpath, npts=100, labels=labels)

# 4. 计算能带与轨道投影
# projection=True 开启 Fat Band 计算
results = solver.solve_bands(projection=True)

# 5. 绘图与保存
# PYATB 通常提供内置绘图功能，或提取 results 数据自行绘制
# results['band'] 包含本征值
# results['proj'] 包含投影权重
print("Fat band calculation done.")
```

### 3.2.2 结果分析
运行脚本后，你将得到包含轨道权重的能带数据。在绘图中，能带的粗细（或颜色）代表了特定轨道（如 Fe 的 $3d$ 轨道）的贡献。这对于分析磁性材料、强关联体系的电子结构至关重要。

---

## 3.3 能带反折叠 (Band Unfolding)

当你计算一个 $2\times2\times2$ 的超胞（例如为了引入一个空位缺陷）时，原本清晰的能带会被折叠进变小了 8 倍的布里渊区中，变得杂乱无章。Band Unfolding 技术可以将这些能带“展开”回原胞的布里渊区，使其与实验 ARPES 结果可比。

### 3.3.1 物理图像
*   **Supercell (SC)**: 实际计算的体系（大实空间，小倒空间）。
*   **Primitive Cell (PC)**: 理想的参考体系（小实空间，大倒空间）。
*   **目标**: 计算 SC 的波函数在 PC 波函数基组上的投影权重（Spectral Weight）。

### 3.3.2 操作流程

1.  **ABACUS 计算**: 对**超胞**结构进行 SCF 计算，并按照 3.1 节的方法输出 HR/SR 矩阵。
2.  **PYATB 设置**: 在 Python 脚本中，除了加载模型，还需要定义从原胞到超胞的变换矩阵。

```python
# ... (加载模型步骤同上，注意这里加载的是超胞的矩阵) ...

# 定义超胞变换矩阵 (Supercell Matrix)
# 例如 2x2x2 超胞
sc_matrix = [
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2]
]

# 设置反折叠求解
# 注意：K 点路径应该是针对“原胞”布里渊区的高对称点
solver.set_kpath(kpath_primitive, npts=100)

# 运行反折叠
# unfolding=True 激活反折叠算法
unfold_data = solver.solve_unfolding(sc_matrix=sc_matrix)

# 绘图
# 结果通常是 (E, k) 空间的谱权重图 (Spectral Function)
# 颜色深的地方代表“有效能带”，模糊的地方代表对称性破缺导致的展宽
```

---

## 附录：常见问题与避坑指南

### 1. 对称性陷阱 (Symmetry Trap)
*   **现象**: PYATB 报错提示矩阵维度不匹配或 K 点映射错误。
*   **原因**: 在进行能带反折叠时，如果 ABACUS 计算超胞时自动通过对称性简化了 K 点或调整了结构，可能导致映射失败。
*   **对策**: 虽然 SCF 计算通常开启对称性，但在处理复杂的反折叠任务时，如果遇到问题，尝试在 ABACUS `INPUT` 中设置 `symmetry 0` 以关闭对称性操作，确保输出的矩阵完全对应输入的 `STRU` 几何。

### 2. 文件依赖 (File Dependencies)
*   **现象**: PYATB 运行时报错 `FileNotFound` 或 `Orbital mismatch`。
*   **原因**: 很多人只拷贝了 `OUT.ABACUS` 文件夹，却忘记了 `.orb` 轨道文件。
*   **对策**: PYATB 需要读取 `.orb` 文件中的截断半径和轨道形状来计算重叠积分修正。**务必确保 `.orb` 文件与 Python 脚本可见，且文件名与 `STRU` 中定义的一模一样。**

### 3. 流程中断 (Workflow Break)
*   **误区**: 试图用 `nscf` 的结果去跑 PYATB。
*   **纠正**: PYATB 需要的是实空间哈密顿量矩阵（HR）。这通常是在 `scf` 计算结束时通过 `out_mat_hs2` 输出的。`nscf` 主要用于在倒空间计算本征值，通常不输出实空间稀疏矩阵。请务必分清这两个流程。

### 4. 物理图像：原胞 vs 超胞
*   **提示**: 在做反折叠时，你的 K 点路径（`kpath`）必须是写在**原胞**的倒格矢基底上的，而不是超胞的。否则你画出来的能带图横坐标会完全错误。