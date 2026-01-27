# 第二章：计算前的关键预处理

> **本章核心逻辑**：“垃圾进，垃圾出（Garbage In, Garbage Out）”。
> 声子谱（Phonon Spectrum）计算本质上是求解体系势能面在平衡位置附近的二阶导数（力常数）。如果初始结构没有处于势能面的极小值点（即受力不为零），一阶导数非零将导致简谐近似失效，直接表现为声子谱中 Gamma 点附近出现巨大的虚频（Imaginary Frequency）。
>
> 本章将带你完成两项核心任务：
> 1.  **极致的几何弛豫**：将原子受力压低到声子计算所需的精度。
> 2.  **微扰构型生成**：利用 Phonopy 生成包含位移的超胞结构。

---

## Section 2.1: 初始结构的几何优化 (Relaxation)

在进行任何声子计算之前，必须对晶胞进行严格的几何优化。对于固体材料，通常建议同时优化原子位置和晶胞参数（Variable Cell Relaxation）。

### 2.1.1 准备 INPUT 文件

为了获得高质量的平衡态结构，我们需要在 `INPUT` 文件中设置比常规电子结构计算更严格的收敛标准。

以下是针对 FCC Al 优化的推荐 `INPUT` 参数配置：

```bash
INPUT_PARAMETERS
# 1. General Parameters
suffix          Al_relax      # 任务后缀名
calculation     cell-relax    # 关键参数：同时优化晶胞形状/体积和原子位置
                              # 若只想固定晶胞优化原子，可设为 'relax'

# 2. Relaxation Control
relax_nmax      100           # 最大离子步数
cal_force       1             # 显式开启力计算（虽然 relax 模式默认开启，但建议显式指定）
cal_stress      1             # 计算应力，cell-relax 必须开启
force_thr_ev    0.001         # 关键参数：力的收敛阈值 (eV/Angstrom)
                              # 声子计算建议设为 0.001 或更低，常规计算通常为 0.01-0.05
stress_thr      1             # 应力收敛阈值 (kbar)，通常设为 1-10 kbar

# 3. Electronic Iteration (SCF)
ecutwfc         60            # 平面波截断能 (Ry)，需根据赝势硬度测试收敛性
scf_thr         1.0e-8        # SCF 收敛阈值，需比 ionic 步精度高 1-2 个数量级
basis_type      lcao          # 基组类型：lcao (数值原子轨道) 或 pw (平面波)
smearing_method mp            # 金属体系推荐 Methfessel-Paxton
smearing_sigma  0.02          # 展宽宽度 (Ry)

# 4. Pseudopotential & Orbitals
pseudo_dir      ./psp         # 赝势目录
orbital_dir     ./orb         # 轨道目录 (仅 lcao 需要)
```

**专家点拨：**
*   **`calculation`**: 初始结构通常来自数据库（如 MP 或 ICSD），其实验晶格常数对应的是有限温度，而 DFT 计算对应 0K。因此，务必使用 `cell-relax` 找到 DFT 理论下的平衡晶格常数，否则会引入内应力，导致声子谱分裂或虚频。
*   **`force_thr_ev`**: 这是声子计算成败的关键。常规计算 0.02 eV/Å 足够，但声子计算对力极其敏感，建议至少收敛到 **0.001 eV/Å** 甚至更低。
*   **`scf_thr`**: 电子步收敛精度必须高于离子步精度，否则计算出的力是不准确的“噪音”，导致几何优化无法收敛。

### 2.1.2 准备 STRU 文件 (FCC Al 示例)

`STRU` 文件定义了初始晶格和原子位置。以下是 FCC Al 的标准格式：

```text
ATOMIC_SPECIES
Al 26.982 Al_ONCV_PBE-1.0.upf upf201  # 元素名 质量 赝势文件名

NUMERICAL_ORBITAL
Al_gga_7au_100Ry_4s4p1d.orb           # 轨道文件名 (仅 lcao 需要)

LATTICE_CONSTANT
1.88972612546                         # 晶格常数缩放因子 (Bohr)

LATTICE_VECTORS
2.0173 2.0173 0.0000                  # 晶格矢量 v1
2.0173 0.0000 2.0173                  # 晶格矢量 v2
0.0000 2.0173 2.0173                  # 晶格矢量 v3

ATOMIC_POSITIONS
Direct                                # 坐标类型：Direct (分数坐标) 或 Cartesian

Al                                    # 元素标签
0                                     # 磁矩 (0 表示无磁性)
1                                     # 原子数量
0.0000 0.0000 0.0000 1 1 1            # 坐标 及 移动限制(1=动, 0=定)
```

### 2.1.3 提取优化后的结构

计算完成后，ABACUS 会将优化后的结构写入输出目录。你需要执行以下操作为下一步做准备：

1.  检查 `OUT.Al_relax/running_cell-relax.log` 确认计算是否收敛（搜索 "Converge" 关键词）。
2.  找到优化后的结构文件：`OUT.Al_relax/STRU_ION_D`。
3.  将其复制并重命名为当前工作目录下的 `STRU`，作为生成超胞的基准。

---

## Section 2.2: 超胞与微扰构型的生成

在有限位移法（Finite Displacement Method）中，我们需要人为地移动原子，计算由此产生的回复力。为了避免周期性边界条件导致原子与其“镜像”发生错误的相互作用，必须建立**超胞（Supercell）**。

我们将使用 **Phonopy** 软件自动完成这一过程。

### 2.2.1 确定扩胞倍数

扩胞的大小决定了力常数计算的截断距离（Cutoff distance）。
*   **原则**：超胞的尺寸应至少大于两倍的原子间相互作用截断半径。
*   **经验值**：对于金属或共价材料，超胞边长通常需要在 **10 Å - 15 Å** 以上。
*   **FCC Al 案例**：原胞较小，我们选择 2x2x2 的扩胞。

### 2.2.2 执行 Phonopy 生成微扰

确保你已经安装了 Phonopy，并且当前目录下有优化好的 `STRU` 文件。在终端执行以下命令：

```bash
phonopy -d --dim="2 2 2" --abacus
```

**参数详解：**
*   **`-d`**: 生成位移（Displacement）。
*   **`--dim="2 2 2"`**: 指定超胞扩胞倍数，分别为 x, y, z 方向。
*   **`--abacus`**: 关键参数！指定输入输出格式为 ABACUS 接口格式。Phonopy 会自动读取当前目录下的 `STRU` 文件。

### 2.2.3 产物解析

执行上述命令后，Phonopy 会根据晶体对称性分析，生成最少数量的必要微扰结构。你会看到以下文件：

1.  **`STRU-001`, `STRU-002`, ...**
    *   这是生成的微扰结构文件。
    *   **FCC Al**：由于对称性极高，通常只会生成一个文件 `STRU-001`。
    *   **低对称性体系**：可能会生成几十个文件。每个文件代表一种独立的位移模式。
    *   **后续任务**：你需要对每一个 `STRU-XXX` 文件进行一次单点能自洽计算（SCF）以获取受力。

2.  **`phonopy_disp.yaml`**
    *   这是一个核心记录文件，记录了哪些原子被移动了、移动的方向和距离。
    *   **警告**：请妥善保存此文件，后续处理力常数（Force Sets）时必须使用它，**切勿修改或删除**。

3.  **`SPOSCAR`** (可选)
    *   这是 Phonopy 生成的完美超胞（无位移），通常用于参考，ABACUS 流程中主要使用 `STRU-XXX`。

### 2.2.4 目录结构规划建议

为了避免文件混乱，建议采用以下目录结构管理后续计算：

```text
.
├── 01_relax/             # 几何优化目录
│   ├── INPUT
│   └── STRU              # 原始结构
├── 02_phonon/            # 声子计算主目录
│   ├── STRU              # 从 01_relax 复制来的优化后结构
│   ├── phonopy_disp.yaml # phonopy -d 生成的记录
│   ├── disp-001/         # 第一个微扰结构计算目录
│   │   ├── INPUT         # 需修改 calculation=scf, cal_force=1
│   │   └── STRU          # 对应 STRU-001
│   ├── disp-002/         # 第二个微扰结构 (如有)
│   │   ├── INPUT
│   │   └── STRU          # 对应 STRU-002
│   └── ...
```

至此，我们已经完成了所有计算前的准备工作。下一章，我们将进入核心环节：批量计算微扰结构的原子受力。