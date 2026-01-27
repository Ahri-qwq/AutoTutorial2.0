# 第二章：计算准备与初猜链生成 (Pre-processing)

> **教授寄语**：
> “垃圾输入导致垃圾输出 (Garbage In, Garbage Out)”。
> 在我多年的计算材料学经验中，超过 80% 的 NEB 计算失败并非源于算法本身，而是源于糟糕的初末态结构或不合理的初始路径插值。NEB 就像是在迷雾中寻找山脊线的最低点，如果你连起点和终点都定错了，或者初始路线直接穿过了一座高山（原子重叠），那么再强大的算法也无能为力。
> 本章我们将学习如何为 ABACUS 准备“完美”的输入数据。

---

## 2.1 核心概念：ASE-ABACUS 的“客户-服务器”模式

在开始操作之前，必须纠正一个许多从 VASP 转过来的用户最常见的**认知误区**：

<div style="background-color: #ffe6e6; padding: 15px; border-left: 5px solid #ff0000; margin-bottom: 20px;">
<strong>⚠️ 核心警告 (Critical Warning)</strong><br>
请不要在 ABACUS 的 <code>INPUT</code> 文件中寻找任何与 NEB 相关的参数（如 <code>IMAGES</code>, <code>SPRING</code> 等）！
</div>

**ABACUS 在此模式下的角色**：
*   **ABACUS (Server)**：它只负责做一件事——给定一个原子结构，计算出能量 (Energy) 和力 (Forces)。它并不知道自己在做 NEB，它以为自己只是在做普通的 SCF 计算。
*   **ASE/ATST-Tools (Client/Driver)**：这是“大脑”。它负责生成一系列结构（Images），调用 ABACUS 算力，收集数据，计算弹簧力，并更新结构位置。

因此，我们的工作流是围绕 **ATST-Tools** 脚本集展开的：
1.  **Make (`neb_make.py`)**: 准备数据，生成初猜链。
2.  **Run (`neb_run.py`)**: 驱动计算，ASE 调用 ABACUS。
3.  **Post (`neb_post.py`)**: 数据分析，可视化。

---

## 2.2 初末态结构优化 (IS & FS Relaxation)

NEB 计算的前提是：反应物 (Initial State, IS) 和产物 (Final State, FS) 必须是势能面上的局部极小点。

### 2.2.1 ABACUS 优化参数设置
我们需要对 IS 和 FS 进行高精度的结构优化。精度不够会导致 NEB 在端点附近产生虚假的“受力”，导致收敛困难。

**推荐的 `INPUT` 参数 (IS/FS 共用)：**

```bash
INPUT_PARAMETERS
# 核心计算控制
calculation     relax       # 必须进行结构优化
relax_nmax      100         # 最大优化步数
cal_force       1           # 必须计算力
cal_stress      1           # 建议计算应力（如果是固体体系）

# 精度控制 (至关重要)
scf_thr         1.0e-7      # 自洽场收敛精度，建议比普通计算高一个量级
force_thr_ev    0.01        # 力收敛判据 (eV/Ang)，越小越好，建议 0.01 或 0.02

# 基组与泛函 (根据体系调整)
basis_type      lcao        # 或 pw
ecutwfc         100         # 截断能
ks_solver       genelpa     # 求解器
...
```

### 2.2.2 原子映射 (Atom Mapping) —— 失败的头号杀手
**原子映射**是指初态的第 $i$ 个原子，必须对应末态的第 $i$ 个原子。
*   **错误案例**：初态中第 1 号原子是左边的氢，末态中第 1 号原子变成了右边的氢。
*   **后果**：NEB 会试图让这两个原子“互换位置”，导致路径能量极高，计算发散。

**操作建议**：
1.  先搭建并优化好 **IS** 结构。
2.  复制 IS 的结构文件 (`STRU`) 作为 FS 的起点。
3.  在可视化软件（如 VESTA, MS）中，**手动移动**关键原子到产物位置，**不要删除或重新添加原子**，以保持原子 ID 顺序不变。
4.  对修改后的 FS 进行结构优化。

---

## 2.3 初始路径插值与 IDPP 方法

当 IS 和 FS 都优化收敛后（假设分别位于文件夹 `IS` 和 `FS` 中），我们需要生成中间的过渡图像（Images）。

### 2.3.1 为什么选择 IDPP？
*   **线性插值 (Linear Interpolation)**：直接连接初末态坐标。
    *   *缺点*：对于旋转或复杂位移，原子可能会直线“穿过”彼此，导致中间图像原子重叠，能量爆炸，SCF 无法收敛。
*   **IDPP (Image Dependent Pair Potential)**：
    *   *原理*：通过构建成对势函数，让插值路径尽量保持键长合理，自动绕开原子重叠区域。
    *   *结论*：**始终推荐使用 IDPP**。

### 2.3.2 使用 `neb_make.py` 生成初猜
我们将使用 ATST-Tools 提供的 `neb_make.py` 脚本。

**基本命令格式**：
```bash
python3 neb_make.py -i [IS_LOG] [FS_LOG] -n [N_IMAGES] --method IDPP
```
*   `[IS_LOG]`: 初态 ABACUS 输出文件路径 (如 `IS/OUT.ABACUS/running_scf.log` 或 `running_relax.log`)。
*   `[FS_LOG]`: 末态 ABACUS 输出文件路径。
*   `[N_IMAGES]`: 中间插入的图像数量（不含首尾）。例如插 5 个点，总共就是 7 个 Image。

**实战案例**：
假设目录结构如下：
```text
.
├── IS/
│   └── OUT.ABACUS/running_relax.log
├── FS/
│   └── OUT.ABACUS/running_relax.log
```

运行命令：
```bash
python3 /opt/ATST-Tools/neb/neb_make.py \
    -i IS/OUT.ABACUS/running_relax.log FS/OUT.ABACUS/running_relax.log \
    -n 5 \
    --method IDPP
```
**输出**：生成文件 `init_neb_chain.traj`。这是 ASE 的标准轨迹文件，包含了插值后的所有结构。

### 2.3.3 处理特殊情况（磁性与固定原子）

**场景 A：磁性体系 (Magnetism)**
ABACUS 的 `INPUT` 文件虽然可以设置磁性，但在 NEB 插值过程中，中间图像的磁矩初始化往往需要显式指定，否则可能导致中间态磁矩丢失（变为非磁），导致能量曲线不连续。

*   **解决方案**：使用 `--mag` 参数。
*   **格式**：`元素:磁矩`。
*   **示例**：铁原子初始磁矩设为 3.0，氧原子设为 0。
    ```bash
    python3 neb_make.py ... --mag Fe:3.0,O:0.0
    ```
    *注意：这会为轨迹文件中的原子添加初始磁矩信息，后续 ASE 调用 ABACUS 时会尝试将此信息传递给计算器（具体依赖于接口实现，建议在后续 `neb_run.py` 中也再次确认磁性设置）。*

**场景 B：表面催化中的原子固定 (Constraints)**
在表面反应中，我们通常固定底部的几层原子。
*   **解决方案**：使用 `--fix` 参数。
*   **格式**：`高度阈值:方向`。
*   **示例**：固定 Z 轴方向 (2) 上，相对坐标小于 0.5 的所有原子。
    ```bash
    python3 neb_make.py ... --fix 0.5:2
    ```

---

## 2.4 进阶策略：AutoNEB 的准备

如果你不知道该插多少个点，或者反应路径非常复杂，**AutoNEB** 是最佳选择。
AutoNEB 的策略是：
1.  先用少量的点（如 3-4 个）和较低的精度（粗糙的 K 点或截断能）跑通一条大概的路径。
2.  程序自动检测哪里能量变化剧烈，就在哪里**动态加点**。
3.  最后对关键区域进行高精度 CI-NEB 计算。

**准备工作的区别**：
对于 AutoNEB，`neb_make.py` 的步骤是一样的，但你通常只需要生成较少的初始点（例如 `-n 3`），后续的加点工作交给 `autoneb_run.py` 脚本在运行时自动处理。

---

## 2.5 检查与验证 (Checklist)

在进入下一章“运行计算”之前，请务必完成以下检查：

1.  [ ] **可视化检查**：使用 `ase gui init_neb_chain.traj` 或将 traj 转换为 cif/xyz 文件查看。
    *   *检查点*：播放动画，确认原子从 IS 运动到 FS 的过程中没有发生“瞬移”或严重的原子重叠。
2.  [ ] **原子映射**：确认 IS 和 FS 的原子数量、种类顺序完全一致。
3.  [ ] **收敛性**：确认 IS 和 FS 的结构优化已经达到 `force_thr_ev` 要求。
4.  [ ] **文件就位**：确认目录下已有 `init_neb_chain.traj`。

做好这些准备，我们就为成功的过渡态搜索打下了坚实的基础。下一章，我们将编写 Python 驱动脚本，正式启动 ABACUS 引擎。