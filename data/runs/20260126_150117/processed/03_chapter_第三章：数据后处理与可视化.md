# 第三章：数据后处理与可视化

在完成了繁重的 SCF 迭代与几何优化后，我们手中的 `.cube` 文件就像是埋藏着物理规律的“原石”。本章的任务，就是通过数学切割（矩阵运算）和光学打磨（可视化），将这些数据转化为能够发表在顶级期刊上的精美图表。

作为开发者，我必须坦诚地告诉你：**ABACUS 内核专注于高效求解薛定谔方程，而非图形渲染。** 因此，本章将重点讲解如何构建“ABACUS + Python/VESTA”的黄金工作流。

---

## Section 3.1: 数据的矩阵运算（差分电荷密度）

差分电荷密度（Charge Density Difference）是分析电荷转移（Charge Transfer）最直观的工具。其核心定义为：
$$ \Delta \rho = \rho_{AB} - \rho_A - \rho_B $$
其中 $\rho_{AB}$ 是总体系电荷密度，$\rho_A$ 和 $\rho_B$ 是保持原子位置不变的孤立部分电荷密度。

### 3.1.1 核心流程：严谨的“三步走”

**特别注意**：ABACUS **不会**直接输出差分电荷密度文件。你必须手动设计三次独立的计算任务。

#### 第一步：计算总体系 (AB)
在 `INPUT` 文件中开启电荷密度输出：
```bash
INPUT_PARAMETERS
# ... 其他参数 ...
out_chg  1       # 开启电荷密度输出，生成 SPIN1_CHG.cube
symmetry 0       # 【关键】关闭对称性，防止格点还原时的旋转误差
```
运行后，将生成的 `SPIN1_CHG.cube` 重命名为 `rho_AB.cube`。

#### 第二步：计算分体系 (A)
修改 `STRU` 文件，**删除** B 部分的原子，但**必须保留** A 部分原子的坐标和晶胞（Lattice）参数完全不变。
*   **LCAO 基组用户的特殊提示**：
    在原子轨道（LCAO）基组下，直接删除原子意味着同时也删除了该原子携带的基函数（Basis Set）。这会导致基组重叠误差（BSSE）。
    *   **常规做法**：直接删除原子。这在定性分析成键时通常是可以接受的。
    *   **高精度做法**：引入“鬼原子”（Ghost Atom）。在 ABACUS 中，这通常需要通过修改元素类型或赝势来实现（即保留基组但不包含核电荷与电子），操作较为高阶。对于初学者，建议先采用“直接删除法”，但需知晓其物理局限性。

#### 第三步：计算分体系 (B)
同理，删除 A 部分原子，计算得到 `rho_B.cube`。

### 3.1.2 风险管控：FFT 网格陷阱

这是新手最容易“翻车”的地方。
$$ \rho_{AB}, \rho_A, \rho_B \text{ 必须定义在完全相同的 FFT 网格上！} $$

ABACUS 会根据晶胞大小和 `ecutwfc` 自动生成 FFT 网格。虽然晶胞不变通常网格也不变，但为了万无一失：
1.  检查第一步计算的日志文件（如 `OUT.ABACUS/running_scf.log`）。
2.  搜索 `FFT GRID` 关键字，找到类似 `108 108 144` 的数值。
3.  （可选但推荐）在第二、三步的 `INPUT` 中显式固定该网格（如果 ABACUS 版本支持直接指定 FFT 参数），或者务必确保 `ecutwfc` 和晶胞参数与第一步**完全一致**。

### 3.1.3 实战脚本：Python 矩阵相减

由于 `.cube` 文件本质上是带头信息的 3D 矩阵，我们可以使用 Python 脚本轻松完成运算。

**脚本示例 (`calc_diff.py`)**:
```python
import numpy as np

def read_cube(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # 解析 Header (前6行是标准Cube格式头信息)
    natoms = int(lines[2].split()[0])
    # 数据从 header + natoms 行之后开始
    start_idx = 6 + natoms
    
    # 提取网格数据
    data = []
    for line in lines[start_idx:]:
        data.extend([float(x) for x in line.split()])
    
    return np.array(data), lines[:start_idx]

def write_cube(fname, data, header):
    with open(fname, 'w') as f:
        f.writelines(header)
        # 格式化输出，每行6个数据
        for i, val in enumerate(data):
            f.write(f"{val:13.5E}")
            if (i + 1) % 6 == 0:
                f.write("\n")
            elif (i + 1) == len(data): # 最后一行
                f.write("\n")

# 主程序
print("Reading rho_AB...")
rho_AB, header = read_cube("rho_AB.cube")
print("Reading rho_A...")
rho_A, _ = read_cube("rho_A.cube")
print("Reading rho_B...")
rho_B, _ = read_cube("rho_B.cube")

# 核心运算：矩阵相减
if len(rho_AB) != len(rho_A) or len(rho_AB) != len(rho_B):
    print("Error: Grid sizes do not match! Check your FFT grids.")
    exit()

rho_diff = rho_AB - rho_A - rho_B

print("Writing rho_diff.cube...")
write_cube("rho_diff.cube", rho_diff, header)
print("Done!")
```

---

## Section 3.2: VESTA 可视化实战

拿到 `rho_diff.cube` 后，我们使用 VESTA 进行可视化。

### 3.2.1 导入与初步设置
1.  打开 VESTA，拖入 `rho_diff.cube`。
2.  此时你可能会看到一个充满整个晶胞的混乱色块，不要慌，这是因为默认的等值面（Isosurface）阈值太低。

### 3.2.2 调整等值面 (Isosurfaces)
1.  点击左侧面板的 **Properties** -> **Isosurfaces**。
2.  **设置正值（电荷积聚）**:
    *   点击列表中的第一项。
    *   **Isosurface level**: 设置为一个合理的正值，例如 `0.002`（具体数值取决于体系大小和电荷转移量，需尝试）。
    *   **Color**: 建议设置为 **黄色** (Yellow)，代表电子获得区域。
3.  **设置负值（电荷耗散）**:
    *   点击 "New" 增加一个等值面。
    *   **Isosurface level**: 设置为对应的负值，例如 `-0.002`。
    *   **Color**: 建议设置为 **青色** (Cyan) 或 **蓝色**，代表电子失去区域。
4.  点击 OK，你将看到清晰的电荷转移图像：电子从蓝色区域流向了黄色区域。

---

## Section 3.3: 定量分析与单位换算

漂亮的图片用于展示，严谨的数据用于论证。在处理 ABACUS 数据时，单位是最大的“坑”。

### 3.3.1 单位陷阱：1/Bohr³
ABACUS 的 `.cube` 文件输出单位默认为 **原子单位 (Atomic Units)**，即：
$$ \text{Unit} = \frac{e}{\text{Bohr}^3} $$

而在文献中，我们通常习惯使用 $e/\text{\AA}^3$。
*   **换算关系**: 1 Bohr $\approx$ 0.529177 $\AA$。
*   **体积因子**: $1 \text{ Bohr}^3 \approx 0.1482 \text{ \AA}^3$。
*   **密度换算公式**:
    $$ \rho (e/\text{\AA}^3) = \rho (e/\text{Bohr}^3) \times 6.748 $$

**警示**：当你进行平面平均电荷密度（Planar Average）或线性电荷分布分析时，务必乘上这个系数，否则你的电荷密度数值将比文献值小一个数量级！

### 3.3.2 Bader 电荷分析前置
虽然差分电荷密度给出了空间的分布，但有时我们需要知道“每个原子具体得失了多少电子”。这就需要 **Bader 电荷分析**。
*   **原理**: 基于电子密度的零通量面（Zero-flux surface）划分原子体积。
*   **操作**:
    1.  使用 Section 3.1 得到的高精度 `rho_AB.cube`（总电荷密度）。
    2.  使用 Henkelman 组开发的 `bader` 程序。
    3.  命令：`bader rho_AB.cube`。
    4.  输出文件 `ACF.dat` 中包含了每个原子的价电子总数。
    5.  **电荷转移量** = `Z_valence` (赝势价电子数) - `Q_bader` (Bader分析出的电子数)。

---

**本章小结**：
从 `.cube` 到差分图，再到 Bader 电荷，我们完成了从“算出数”到“看懂物理图像”的跨越。下一章，我们将进入更激动人心的领域：**能带结构与态密度（DOS）的计算与分析**。