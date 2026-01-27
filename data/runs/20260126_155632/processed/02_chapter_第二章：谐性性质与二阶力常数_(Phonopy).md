# 第二章：谐性性质与二阶力常数 (Phonopy)

在计算晶格热导率的宏大工程中，二阶力常数（2nd Order Force Constants）是地基。它决定了声子谱（Phonon Spectrum）的形状、声子的群速度以及比热容。如果地基打不牢——例如力常数中包含数值噪音——后续的三阶散射计算将毫无意义。

本章我们将利用 **Phonopy** 结合 **ABACUS**，完成从微扰结构生成、受力计算到力常数矩阵提取的全过程。

> **⚠️ 流程全景提示**
> 计算晶格热导率涉及多个软件的接力，请务必保持清醒的头脑：
> 1.  **Phase 1 (本章)**: ABACUS + Phonopy $\rightarrow$ **二阶力常数** (需单位转换)。
> 2.  **Phase 2 (下一章)**: ABACUS + thirdorder $\rightarrow$ **三阶力常数**。
> 3.  **Phase 3**: ShengBTE $\rightarrow$ **热导率**。

---

## 2.1 超胞建立与微扰生成 (Phonopy)

为了计算力常数，我们需要采用有限位移法（Finite Displacement Method）。即在一个扩大的超胞（Supercell）中，将原子稍微移动一点点，计算它受到的回复力。

### 2.1.1 准备工作
请确保你已安装 Phonopy，并准备好优化过的晶胞结构文件 `STRU`。
在工作目录下创建一个名为 `2nd` 的文件夹，并进入该目录。

### 2.1.2 配置文件 setting.conf
编写 Phonopy 的配置文件 `setting.conf`。这里最关键的参数是 `DIM`，它定义了超胞的大小。

```ini
# setting.conf
DIM = 2 2 2
ATOM_NAME = Si
```

*   **DIM**: 指定在 $x, y, z$ 三个方向上扩胞的倍数。
    *   *教学演示*: `2 2 2` (为了快速跑通流程)。
    *   *科研实战*: **必须进行收敛性测试**。对于热导率计算，声子平均自由程较长，通常需要较大的超胞（如 4x4x4 或更大）来消除周期性边界条件引入的非物理相互作用。

### 2.1.3 生成微扰结构
使用 Phonopy 的 ABACUS 接口生成带有微扰的结构文件：

```bash
phonopy setting.conf --abacus -d
```

**执行结果**:
Phonopy 会根据晶体对称性，生成最少数量的必要微扰结构。
*   生成文件名为 `STRU-001`, `STRU-002` ...
*   生成 `phonopy_disp.yaml`：记录了具体的位移信息（后续处理必须保留此文件）。

---

## 2.2 ABACUS 力常数计算 (SCF)

现在我们需要计算这些微扰结构中原子的受力。这是一个标准的 DFT 自洽计算（SCF）过程。

### 2.2.1 输入文件 INPUT 设置
编写 `INPUT` 文件。这里有两个至关重要的细节，直接决定了你计算出的热导率是物理结果还是“随机数”。

```bash
INPUT_PARAMETERS
# ... 基础参数 (ecutwfc, basis_type 等) ...

calculation     scf         # 进行自洽计算
stru_file       STRU-001    # 指定读取的结构文件

# --- 精度控制 (生死攸关) ---
scf_thr         1e-8        # LCAO 基组推荐值
# scf_thr       1e-12       # PW (平面波) 基组推荐值
force_thr_ev    1e-7        # 力的收敛标准
cal_force       1           # 显式开启受力计算输出
```

### 🛑 专家锦囊：关于精度的血泪教训
在常规的能带计算中，`scf_thr` 设置为 `1e-6` 通常就够了。但在**声子和热导率计算**中，这远远不够！
*   **微扰极其微小**: 原子位移通常只有 0.01 Å 量级，产生的回复力非常微弱。
*   **噪音放大**: 如果 SCF 收敛精度不高，电荷密度的数值噪音会导致受力计算出现误差。这些误差在后续计算力常数矩阵时会被放大，导致声子谱出现虚频（软模），或者热导率数值完全错误。
*   **推荐标准**:
    *   **LCAO (数值原子轨道)**: 至少 `1e-8`。
    *   **PW (平面波)**: 建议 `1e-12`。

### 2.2.2 批量计算
你需要对 Phonopy 生成的每一个 `STRU-XXX` 文件运行一次 ABACUS。
*   **技巧**: 利用 `stru_file` 参数。你不需要把 `STRU-001` 重命名为 `STRU`，只需在 `INPUT` 中修改 `stru_file = STRU-001` 即可。
*   **脚本建议**: 编写一个简单的 Shell 脚本循环运行，或者使用作业调度系统提交数组作业。

---

## 2.3 力常数提取与单位转换 (核心坑点)

当所有 SCF 计算完成后，我们收集受力并计算力常数。这里存在一个**极易踩中的单位陷阱**。

### 2.3.1 收集受力 (FORCE_SETS)
假设你的计算结果保存在对应的文件夹或日志中（例如 `OUT.STRU-001/running_scf.log`），使用以下命令收集受力：

```bash
# 假设只有一个微扰结构，日志文件路径需根据实际情况调整
phonopy -f OUT.STRU-001/running_scf.log
```
这将生成 `FORCE_SETS` 文件。

### 2.3.2 计算二阶力常数矩阵
编写 `band.conf` 用于生成力常数。

```ini
# band.conf
ATOM_NAME = Si
DIM = 2 2 2
# ... 其他声子谱相关参数 ...

# --- 关键参数 ---
FORCE_CONSTANTS = WRITE        # 输出力常数文件
FULL_FORCE_CONSTANTS = .TRUE.  # 【必须开启】输出完整的力常数矩阵
```

运行 Phonopy：
```bash
phonopy -p band.conf --abacus
```
这将生成 `FORCE_CONSTANTS` 文件。

### 2.3.3 单位转换 (ShengBTE 兼容性处理)

**这是本章最大的坑点，请仔细阅读。**

*   **ABACUS/Phonopy 输出**: Phonopy 读取 ABACUS 结果时，默认输出的力常数单位通常包含原子单位制（Bohr/au），即 `eV/(Å·au)` (其中 1 au $\approx$ 0.529 Å)。
*   **ShengBTE 输入要求**: ShengBTE 严格要求二阶力常数单位为 `eV/Å²`。

如果不进行转换，直接将 `FORCE_CONSTANTS` 喂给 ShengBTE，计算出的热导率将是错误的量级。

**解决方案**:
我们需要使用一个辅助脚本 `au2si.py` 来完成这个转换。
> **注**: `au2si.py` 并非 ABACUS 内置命令，它通常随本教程的案例文件提供（位于 `examples/interface_ShengBTE/LCAO/2nd/` 目录下）。

确保 `au2si.py` 在当前目录下，运行：

```bash
python au2si.py
```

该脚本会读取当前的 `FORCE_CONSTANTS`，进行单位换算，并生成一个新的文件（通常命名为 `FORCE_CONSTANTS_2ND` 或覆盖原文件，具体视脚本实现而定）。**请务必将转换后的文件重命名为 `FORCE_CONSTANTS_2ND`**，这是 ShengBTE 识别的标准文件名。

---

## 2.4 本章总结与风险提示

至此，我们已经获得了符合 ShengBTE 要求的二阶力常数文件 `FORCE_CONSTANTS_2ND`。

**再次强调风险**:
在本教程的案例中，为了演示速度，我们使用了 `2x2x2` 的超胞和较稀疏的 K 点。
*   **案例结果**: Si 的热导率可能算出 ~100 W/mK。
*   **实验真值**: Si 的热导率约为 150 W/mK。
*   **原因**: 超胞尺寸不足以截断长程力常数，且 K 点采样未收敛。

**科研建议**: 在正式计算中，请务必测试超胞大小（如对比 3x3x3, 4x4x4）对声子谱和热导率的影响，直至结果收敛。