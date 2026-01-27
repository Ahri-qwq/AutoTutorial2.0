# 第四章：ShengBTE 热导率计算与分析

万事俱备，只欠东风。在前面的章节中，我们已经利用 ABACUS 结合 Phonopy 计算了二阶力常数，并结合 thirdorder 程序计算了三阶力常数。本章将作为“总装车间”，指导你如何配置 ShengBTE 的核心输入文件 `CONTROL`，运行玻尔兹曼输运方程（BTE）求解器，并对最终的热导率结果进行科学分析。

## 4.0 全流程回顾与数据准备

在正式运行 ShengBTE 之前，让我们通过一个逻辑图解回顾整个工作流，确保你手头的文件是正确且完整的。这是一个多软件耦合的复杂流程，清晰的阶段划分至关重要。

```mermaid
graph TD
    subgraph Stage 1: 2nd Order [ABACUS + Phonopy]
        A[结构优化 (Relax)] --> B[Phonopy 扩胞]
        B --> C[ABACUS SCF 计算受力]
        C --> D[Phonopy 生成 FORCE_CONSTANTS]
        D -->|关键: 单位转换| E(FORCE_CONSTANTS_2ND)
        style E fill:#f9f,stroke:#333,stroke-width:2px
    end

    subgraph Stage 2: 3rd Order [ABACUS + thirdorder]
        F[thirdorder 生成微扰构型] --> G[ABACUS SCF 计算受力]
        G -->|关键: 格式转换| H[aba2vasp.py -> vasprun.xml]
        H --> I[thirdorder 收集数据]
        I --> J(FORCE_CONSTANTS_3RD)
        style J fill:#f9f,stroke:#333,stroke-width:2px
    end

    subgraph Stage 3: Solver [ShengBTE]
        E --> K[ShengBTE]
        J --> K
        L[CONTROL 文件] --> K
        K --> M[BTE.kappa 热导率结果]
    end
```

**检查清单：**
进入你的 `shengbte` 工作目录，确保包含以下三个核心文件：
1.  **`FORCE_CONSTANTS_2ND`**: 由 Phonopy 产生并经过 `au2si.py` 转换单位后的二阶力常数。
2.  **`FORCE_CONSTANTS_3RD`**: 由 thirdorder 产生的原始三阶力常数。
3.  **`CONTROL`**: 我们即将编写的 ShengBTE 配置文件。

---

## 4.1 CONTROL 文件配置

`CONTROL` 是 ShengBTE 的主输入文件，采用 Fortran 的 Namelist 格式。它定义了晶体结构、超胞大小、Q 点采样网格以及求解参数。

以下是针对硅（Si）案例的标准 `CONTROL` 文件示例及详细解析：

```fortran
&allocations
    nelements = 1,          ! 元素种类的数量
    natoms = 2,             ! 单胞（Unit Cell）中的原子总数
    ngrid(:) = 10 10 10,    ! Q 点网格采样 (Q-grid)，用于积分声子玻尔兹曼方程
    norientations = 0,      ! 纳米线/薄膜方向设置，体材料设为 0
/

&crystal
    lfactor = 0.1,          ! 长度单位转换因子。ShengBTE 内部使用 nm，若输入为 Angstrom，此处设为 0.1
    lattvec(:,1) = 0.0 2.81594778072 2.81594778072,  ! 晶格矢量 a1 (Angstrom)
    lattvec(:,2) = 2.81594778072 0.0 2.81594778072,  ! 晶格矢量 a2
    lattvec(:,3) = 2.81594778072 2.81594778072 0.0,  ! 晶格矢量 a3
    elements = "Si",        ! 元素名称
    types = 1 1,            ! 每个原子对应的元素类型索引（对应 elements 列表）
    positions(:,1) = 0.00 0.00 0.00,  ! 原子 1 的分数坐标
    positions(:,2) = 0.25 0.25 0.25,  ! 原子 2 的分数坐标
    scell(:) = 2 2 2,       ! 计算力常数时使用的超胞大小 (Supercell size)
/

&parameters
    T = 300,                ! 计算热导率的温度 (Kelvin)
    scalebroad = 1.0,       ! 高斯展宽的缩放因子，通常设为 1.0
/

&flags
    nonanalytic = .FALSE.,  ! 是否包含非解析项校正（极性材料需设为 .TRUE. 并提供 BORN 文件）
    nanowires = .FALSE.,    ! 是否计算纳米线
/
```

### 关键参数详解

1.  **`scell` (Supercell Size)**:
    *   **含义**: 这里填写的必须是你**计算二阶和三阶力常数时所用的超胞大小**。
    *   **注意**: 在本教程的 Si 案例中，我们为了演示速度使用了 2x2x2 的超胞，因此这里填 `2 2 2`。如果二阶和三阶使用了不同大小的超胞，ShengBTE 会尝试进行插值，但强烈建议两者保持一致以保证精度。

2.  **`ngrid` (Q-grid)**:
    *   **含义**: 这是求解 BTE 时在布里渊区采样的网格密度。
    *   **策略**: 该参数直接决定结果的收敛性。`10 10 10` 仅为演示用。实际科研中，你需要测试 `15 15 15`, `20 20 20` 等，直到热导率数值不再剧烈变化。

3.  **`lattvec` & `positions`**:
    *   **来源**: 这些数据应直接取自你进行结构优化（Relax）后的 `STRU` 文件或转换后的 `POSCAR`。务必保证精度，不要随意截断小数位。

4.  **`lfactor`**:
    *   ABACUS 和大多数 DFT 软件使用 Å (Angstrom) 为单位，而 ShengBTE 内部使用 nm。因此，当 `lattvec` 单位为 Å 时，必须设置 `lfactor = 0.1`。

---

## 4.2 运行与结果分析

### 4.2.1 提交计算任务

ShengBTE 是一个支持 MPI 并行的程序。在准备好 `CONTROL`, `FORCE_CONSTANTS_2ND`, `FORCE_CONSTANTS_3RD` 后，使用以下命令运行：

```bash
# 假设使用 10 个核心并行计算
mpirun -n 10 ShengBTE
```

程序运行过程中会输出详细的日志。如果配置正确，计算通常在几分钟内完成（取决于 `ngrid` 的密度）。

### 4.2.2 结果分析

运行结束后，目录下会生成多个输出文件。最核心的文件是 **`BTE.kappa`**。

**`BTE.kappa` 文件结构示例：**
```text
# T          kappa_xx      kappa_yy      kappa_zz      ... (张量分量)
  300.0000   102.4512      102.4512      102.4512      ...
```

**数据解读与对比：**

*   **计算值**: 在本案例（LCAO 基组，2x2x2 超胞，2x2x2 K点）中，300K 下的计算结果大约在 **100 W/(m·K)** 左右。
*   **实验值**: 硅单晶在 300K 的实验热导率约为 **140 - 150 W/(m·K)**。

**为什么计算值偏小？（风险提示）**
你可能会发现计算结果（~100）明显低于实验值（~150）。这**不是**软件的错误，而是由于我们在教学案例中为了节省计算时间，采用了极低的参数设置：
1.  **超胞过小**: 2x2x2 的超胞截断了长程力常数，无法捕捉长波声子的贡献。
2.  **K 点稀疏**: 电子结构计算时的 2x2x2 K 点采样对于半导体 Si 来说远远不够。
3.  **Q 点稀疏**: `CONTROL` 中的 `ngrid` 仅为 10x10x10。

**科研级建议**: 在实际研究中，你必须进行严格的收敛性测试。通常对于 Si 这样的体系，需要 4x4x4 或 5x5x5 的超胞，以及更密的 K 点网格，才能得到与实验吻合的结果。

---

## 附录：常见问题与进阶建议

### 1. 常见报错与排查

*   **报错**: `Error reading FORCE_CONSTANTS`
    *   **原因 1**: Phonopy 计算二阶力常数时，未在 `band.conf` 或 `setting.conf` 中设置 `FULL_FORCE_CONSTANTS = .TRUE.`。ShengBTE 需要完整的力常数矩阵，而不仅仅是简化的。
    *   **原因 2**: 单位未转换。ABACUS+Phonopy 输出的单位是 eV/(Å·au)，ShengBTE 需要 eV/Å²。
    *   **解决**: 确保运行了 `python au2si.py` 生成 `FORCE_CONSTANTS_2ND`。

*   **现象**: 热导率结果包含大量噪音或数值异常（如负值、极大值）。
    *   **原因**: 三阶力常数计算精度不足。这是新手最容易踩的坑。
    *   **解决**: 检查计算三阶力常数时的 SCF 收敛精度。
        *   **LCAO 基组**: `scf_thr` 至少设为 **1e-8**。
        *   **PW (平面波) 基组**: `scf_thr` 必须设为 **1e-12**。
        *   **警告**: 默认的 `1e-6` 对于微扰计算完全不够，会导致数值微分结果全是数值噪音。

### 2. 辅助脚本清单

本教程中提到的脚本并非 ABACUS 内置命令，而是随案例提供的 Python 辅助工具。请确保它们在你的工作目录中：

| 脚本名 | 功能描述 | 依赖库 |
| :--- | :--- | :--- |
| **`au2si.py`** | 将 Phonopy 输出的力常数单位从原子单位制转换为 ShengBTE 所需的 eV/Å²。 | NumPy |
| **`pos2stru.py`** | 将 `thirdorder` 生成的 `POSCAR` 格式微扰结构转换为 ABACUS 的 `STRU` 格式。 | ASE |
| **`aba2vasp.py`** | 将 ABACUS 的 SCF 输出（力、能量、应力）封装成 `vasprun.xml` 格式，以便 `thirdorder` 读取。 | 无 |

### 3. 进阶：如何测试收敛性

在发表论文前，请务必完成以下测试：
1.  **K 点收敛**: 固定超胞大小，增加 DFT 计算时的 K 点密度（如 2x2x2 -> 4x4x4 -> 6x6x6），观察二阶声子谱是否稳定。
2.  **超胞收敛**: 这是最耗时的。计算 3x3x3, 4x4x4, 5x5x5 超胞的三阶力常数。通常三阶力常数的截断半径影响最大。
3.  **Q 网格收敛**: 在 `CONTROL` 文件中增加 `ngrid`（如 20x20x20 -> 30x30x30），直到 `BTE.kappa` 收敛。这一步计算成本极低，应首先保证收敛。