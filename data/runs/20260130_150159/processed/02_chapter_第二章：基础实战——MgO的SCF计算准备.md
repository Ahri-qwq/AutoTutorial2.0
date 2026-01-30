# 第二章：基础实战——MgO的SCF计算准备

在计算材料学领域，MgO（氧化镁）常被视作“Hello World”级别的标准测试体系。它结构简单（岩盐矿结构）、非磁性且为绝缘体，非常适合初学者用来熟悉从结构文件到第一性原理计算的完整工作流。

本章将演示如何利用 `abacustest` 工具，从一个通用的 CIF 晶体结构文件出发，一键生成符合 ABACUS 标准的输入文件夹。这一过程将自动处理最令人头疼的赝势（Pseudopotential）与数值原子轨道（Numerical Atomic Orbitals）匹配问题。

## 2.1 一键生成输入文件夹

在传统的 DFT 计算准备中，用户通常需要手动转化结构文件格式、去赝势库寻找对应的元素文件、编写输入参数，并确保轨道文件与赝势匹配。而在 ABACUS 的现代工作流中，我们可以通过 `abacustest` 极大地简化这一过程。

假设你当前目录下已经准备好了一个名为 `MgO.cif` 的结构文件，且已经配置好了推荐的 `APNS-pp-orb-v1` 赝势轨道库（即环境变量 `ABACUS_PP_PATH` 和 `ABACUS_ORB_PATH` 已设置），只需执行以下一条命令即可生成完整的 SCF 计算任务：

```bash
abacustest model inputs -f MgO.cif --ftype cif --lcao --folder-syntax MgO
```

### 命令参数详解

这条命令虽然简短，但包含了几个定义计算核心的关键参数：

*   **`-f MgO.cif`**：指定输入的原始结构文件。
*   **`--ftype cif`**：明确告诉程序输入文件的格式是 CIF。ABACUS 的原生结构文件格式是 STRU，但通过此工具我们可以直接使用 CIF 或 VASP POSCAR 等通用格式。
*   **`--lcao`**：**这是最关键的参数之一**。它指定计算将使用 **LCAO（线性组合原子轨道）** 基组。如果不加此选项，工具默认会准备平面波（PW）基组的输入文件。由于 LCAO 计算需要额外的 `.orb` 轨道文件，加上此标记后，程序会自动去库中检索并链接对应的轨道文件。
*   **`--folder-syntax MgO`**：指定生成的 ABACUS 输入文件夹的名称。这里我们将文件夹命名为 `MgO`。如果不设置此参数，程序将默认使用数字编号（如 `000000`）作为文件夹名。

> **专家提示**：如果你的环境中没有设置赝势库的环境变量，你也可以在上述命令中通过 `--pp /path/to/pp` 和 `--orb /path/to/orb` 显式指定路径。

---

## 2.2 产物解析与文件结构

执行上述命令后，当前目录下会生成一个名为 `MgO` 的文件夹。理解这个文件夹内的结构对于掌握 ABACUS 的运行机制至关重要。

让我们进入该文件夹查看其结构：

```bash
MgO
├── INPUT
├── STRU
├── Mg_gga_10au_100Ry_2s1p.orb -> /path/to/apns-orbitals-efficiency-v1/Mg_gga_10au_100Ry_2s1p.orb
├── Mg.PD04.PBE.UPF -> /path/to/apns-pseudopotentials-v1/Mg.PD04.PBE.UPF
├── O_gga_6au_100Ry_2s2p1d.orb -> /path/to/apns-orbitals-efficiency-v1/O_gga_6au_100Ry_2s2p1d.orb
├── O.upf -> /path/to/apns-pseudopotentials-v1/O.upf
├── struinfo.txt
├── run.sh
└── setting.json
```

### 1. 核心输入文件：INPUT

`INPUT` 文件是 ABACUS 的总控中心。`abacustest` 自动为我们生成了一套适合绝大多数 LCAO-SCF 计算的默认参数：

```python
INPUT_PARAMETERS
calculation     scf           # 计算类型：自洽场计算
basis_type      lcao          # 基组类型：原子轨道
ecutwfc         100           # 波函数截断能 (Ry)
scf_thr         1e-07         # SCF收敛阈值
scf_nmax        100           # 最大SCF迭代步数
mixing_type     broyden       # 电荷密度混合方法
mixing_beta     0.8           # 混合参数，0.8为经验推荐值
kspacing        0.14          # K点间距 (1/bohr)，自动生成K点网格
smearing_method gauss         # 展宽方法
smearing_sigma  0.015         # 展宽宽度 (Ry)
ks_solver       genelpa       # 对角化求解器
precision       double        # 计算精度
```

*   **`calculation scf`**: 明确本次任务是进行电子密度的自洽迭代，不涉及离子步（即不移动原子位置）。
*   **`kspacing 0.14`**: 在这个自动生成的配置中，你可能没有看到独立的 `KPT` 文件。这是因为使用了 `kspacing` 参数，ABACUS 会根据倒空间距离自动生成 K 点网格。对于绝缘体 MgO，0.14 (1/bohr) 是一个相对密集的采样，足以保证精度。
*   **`mixing_type` & `mixing_beta`**: 控制电荷密度混合的策略。对于宽带隙绝缘体 MgO，默认的 `broyden` 和 `0.8` 的混合系数通常能实现快速收敛。

### 2. 结构文件：STRU

`STRU` 文件是 ABACUS 能够识别的结构格式，由 CIF 文件转化而来。它包含了晶格常数、晶格矢量、原子种类（及其对应的赝势文件）、数值轨道文件以及原子的坐标信息。

### 3. 赝势与轨道文件的软链接

你会注意到目录下有 `.upf`（赝势）和 `.orb`（轨道）文件，它们以**软链接（Symbolic Link）**的形式存在，指向了你系统中的 `APNS` 库。
*   **机制优势**：这种机制避免了在每个计算文件夹中重复拷贝庞大的数据文件，既节省磁盘空间，又确保了所用参数的一致性。
*   **文件对应**：
    *   `Mg.PD04.PBE.UPF`：镁的 PBE 泛函赝势。
    *   `Mg_gga_10au_100Ry_2s1p.orb`：镁的配套数值轨道，文件名中的 `10au` 代表截断半径，`2s1p` 代表轨道构型（DZP）。

### 4. 辅助文件

*   **`run.sh` & `setting.json`**：这是为使用 `dpdispatcher` 提交任务准备的脚本。如果你是在本地直接运行或使用简单的 Slurm 脚本，可以忽略或删除它们。
*   **`struinfo.txt`**：记录了生成该 `STRU` 文件的原始结构路径等元数据，用于追溯。

至此，一个标准的 MgO SCF 计算任务文件夹已经准备就绪，无需手动编辑任何文件即可直接提交计算。