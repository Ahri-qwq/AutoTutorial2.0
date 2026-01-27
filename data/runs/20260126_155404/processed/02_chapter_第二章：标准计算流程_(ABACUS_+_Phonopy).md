# 第二章：标准计算流程 (ABACUS + Phonopy)

本章将深入 ABACUS 与 Phonopy 联合使用的核心实战环节。我们将按照工业标准的“前处理 -> 自洽计算 -> 后处理”逻辑，详细拆解如何获取材料的声子谱及热力学性质。

**本章核心逻辑：**
1.  **前处理**：基于高精度弛豫结构，使用 Phonopy 生成一系列微扰后的超胞（Supercells）。
2.  **力场计算**：配置 ABACUS 计算这些超胞的原子受力（Forces）。
3.  **后处理**：收集力数据，构建力常数矩阵，计算热容、熵等热力学函数。

---

## 2.1 结构优化与前处理

在进行声子计算之前，**必须**确保晶胞处于受力极小的基态。任何残留的内应力都会导致声子谱出现虚频（软模），尤其是在 $\Gamma$ 点附近。

### 2.1.1 高精度结构弛豫
在上一章的结构优化中，我们通常使用默认收敛标准。但在声子计算前，建议使用更严格的标准：
*   **`force_thr_ev`**: 建议设为 `1.0e-4` eV/Å 或更小。
*   **`stress_thr`**: 建议设为 `1.0e-2` GPa（如果是变胞优化）。

### 2.1.2 生成超胞与微扰文件
我们将使用 Phonopy 的有限位移法（Finite Displacement Method）来生成输入结构。

**步骤 1：准备单位原胞**
准备好优化后的结构文件 `STRU`。

**步骤 2：执行 Phonopy 命令**
在终端运行以下命令生成超胞和微扰文件：

```bash
# 假设扩胞倍数为 2x2x2
phonopy --abacus -d --dim="2 2 2" -c STRU
```

**参数详解：**
*   `--abacus`: 指定输入输出格式为 ABACUS 接口。
*   `-d`: 生成位移（Displacement）。
*   `--dim="2 2 2"`: 指定超胞在 x, y, z 方向的扩胞倍数。
*   `-c STRU`: 指定输入的初始结构文件名为 `STRU`。

**输出产物：**
执行成功后，目录下会生成以下关键文件：
1.  `SPOSCAR`: 完美的超胞结构（未微扰），用于检查扩胞是否正确。
2.  `STRU-001`, `STRU-002`, ... : 包含微扰的 ABACUS 结构文件。
3.  `phonopy_disp.yaml`: 记录了原子位移与文件对应关系的映射表（**千万不要删除**）。

---

## 2.2 ABACUS 力场计算配置

这一步是计算量最大的环节。我们需要对每一个生成的 `STRU-xxx` 文件运行一次 ABACUS 的 SCF（自洽场）计算，以获取精确的原子受力。

### 2.2.1 目录管理建议
为了避免文件混乱，建议采用如下目录结构：

```text
.
├── disp-001/
│   ├── STRU      <-- (原 STRU-001 重命名而来)
│   ├── INPUT
│   └── KPT
├── disp-002/
│   ├── STRU      <-- (原 STRU-002 重命名而来)
│   ├── INPUT
│   └── KPT
...
```

### 2.2.2 INPUT 文件核心配置
在计算声子时，我们只需要“力”的信息，不需要进行结构优化或分子动力学模拟。

**关键参数清单 (INPUT):**

```bash
INPUT_PARAMETERS
# 1. 基础计算模式
calculation     scf        # 进行自洽计算
basis_type      lcao       # 或 pw，需与赝势匹配
symmetry        0          # 【重要】关闭对称性分析，防止ABACUS对微扰结构进行错误的对称性约化

# 2. 力与应力输出 (Phonopy 接口核心)
cal_force       1          # 【必须】开启力计算，输出到 running_scf.log
cal_stress      1          # 推荐开启，用于辅助判断收敛性

# 3. 精度控制 (需与结构优化时保持一致或更高)
ecutwfc         60         # 平面波截断能 (Ry)
scf_thr         1.0e-8     # 自洽收敛阈值
scf_nmax        100        # 最大迭代步数

# 4. 涂抹设置 (金属体系必填)
smearing_method gaussian
smearing_sigma  0.01
```

**风险提示：**
*   **`symmetry 0`**: 微扰后的结构对称性已被破坏，强制关闭对称性分析可以避免潜在的坐标系旋转问题，确保输出的力与 Phonopy 预期的原子索引一一对应。
*   **`cal_force 1`**: 如果遗漏此参数，日志中将不会包含 `FORCE` 字段，导致后续步骤失败。

### 2.2.3 批量提交任务
对于包含数十个微扰结构的任务，建议编写 Shell 脚本批量提交：

```bash
#!/bin/bash
# 假设生成的文件夹名为 disp-001, disp-002 ...
for dir in disp-*; do
    cd $dir
    # 提交任务命令，例如使用 slurm:
    # sbatch submit.sh
    # 或者直接运行 (仅供测试):
    # OMP_NUM_THREADS=1 mpirun -np 8 abacus > output.log &
    cd ..
done
```

---

## 2.3 后处理与热力学性质计算

当所有 SCF 任务完成后，我们需要从 ABACUS 的日志文件 `running_scf.log` 中提取力数据，并计算热力学性质。

### 2.3.1 提取力常数 (FORCE_SETS)
Phonopy 需要读取每个微扰结构的受力情况。

**命令：**
```bash
phonopy --abacus -f disp-001/running_scf.log disp-002/running_scf.log ...
```
*注意：可以使用通配符 `disp-*/running_scf.log`，但必须确保文件顺序与 `phonopy_disp.yaml` 中的顺序一致（通常 shell 的默认排序即可，但需小心 1 和 10 的排序问题）。*

**检查点：**
*   命令执行后，当前目录下应生成 `FORCE_SETS` 文件。
*   **常见报错**：如果报错提示找不到力数据，请检查 `running_scf.log` 结尾是否显示 `JOB DONE`，以及文件中是否包含 `FORCE (eV/Angstrom)` 数据块。

### 2.3.2 计算热力学性质
我们需要编写一个配置文件 `mesh.conf` 来指定计算参数。

**配置文件示例 (`mesh.conf`):**
```ini
DIM = 2 2 2           # 必须与前处理时的扩胞倍数一致
MP = 31 31 31         # Q点网格密度 (Mesh Sampling)，越大越精确
TMIN = 0              # 起始温度 (K)
TMAX = 1000           # 结束温度 (K)
TSTEP = 10            # 温度步长 (K)
```

**执行计算：**
```bash
phonopy --abacus -p mesh.conf
```
*   `-p`: Plot，直接画图并输出数据。

### 2.3.3 结果文件深度解读
运行结束后，重点关注 `thermal_properties.yaml` 文件。

**关键物理量与单位：**

| 字段 (Field) | 物理意义 | 单位 (Unit) | 备注 |
| :--- | :--- | :--- | :--- |
| **temperature** | 温度 | **K** | |
| **free_energy** | 亥姆霍兹自由能 ($F$) | **kJ/mol** | 注意：包含零点能 (ZPE) 和振动熵贡献 ($F = E_{ZPE} - TS_{vib}$)。<br>**不包含** DFT 基态能量 $E_{el}$。 |
| **entropy** | 振动熵 ($S$) | **J/K/mol** | |
| **heat_capacity** | 等容热容 ($C_v$) | **J/K/mol** | 在固体物理近似下，$C_v \approx C_p$ (低温时)。 |

**专家提示：**
这里的 Free Energy **不包含** $PV$ 项。在固相计算中，通常假设 $P \approx 0$，因此吉布斯自由能 $G \approx F$。如果你需要计算相图，记得手动加上 DFT 计算得到的基态能量 $E_{el}$。

---

## 特别专题：处理 Selective Dynamics (手动方案)

**警告：高风险操作**
目前的 ABACUS+Phonopy 接口**不支持**像 VASP 那样通过 `selective dynamics` 标签自动固定特定原子。如果你只需要计算表面吸附分子或特定局域原子的振动模式，必须采用以下**手动替代方案**。此过程极易出错，请务必核对原子索引。

### 操作流程

1.  **手动提取子结构**：
    从你的超胞结构中，手动提取出你**关心**的原子（例如吸附分子 A），记为集合 $S_A$。其余固定不动的原子记为 $S_B$。

2.  **生成微扰 (针对 $S_A$)**：
    创建一个只包含 $S_A$ 原子的临时结构文件 `STRU_A`。使用 Phonopy 对 `STRU_A` 生成微扰文件。
    *注意：此时 `DIM` 通常设为 `1 1 1`，因为你是在操作子集。*

3.  **重组结构 (Re-embedding)**：
    这是最繁琐的一步。对于生成的每一个 `STRU_A-xxx`：
    *   读取其坐标。
    *   将 $S_B$ 原子的坐标（保持基态位置）**手动追加**回文件中。
    *   **关键**：必须保证原子在文件中的顺序与原始超胞完全一致！如果 $S_A$ 原本在 $S_B$ 之后，重组时也要保持这个顺序。

4.  **运行 ABACUS**：
    对重组后的所有结构运行 SCF 计算。

5.  **过滤力数据 (Filtering)**：
    ABACUS 会输出所有原子的力。但 Phonopy 基于 `STRU_A` 生成的 `phonopy_disp.yaml` 只期待看到 $S_A$ 原子的力。
    *   你需要编写脚本，读取 `running_scf.log`。
    *   **只保留** $S_A$ 对应原子的力数据。
    *   构造一个新的符合 Phonopy 格式的力文件（或手动修改 `FORCE_SETS`）。

6.  **生成 FORCE_SETS**：
    使用修改后的数据让 Phonopy 计算力常数。

**建议**：除非计算资源极其受限，否则建议对整个体系进行完整的声子计算，然后在后处理阶段通过分析本征矢（Eigenvectors）来提取局部振动模式，这比上述手动方案更安全。