# 原子轨道基组下的能带反折叠：ABACUS 与 PYATB 实战指南

## 前言

欢迎来到《原子轨道基组下的能带反折叠：ABACUS 与 PYATB 实战指南》。

在计算材料学领域，处理掺杂、缺陷、合金以及表面态等复杂体系时，超胞（Supercell）方法是不可或缺的工具。然而，超胞计算带来的布里渊区折叠效应往往使能带结构变得杂乱无章，难以与实验中的角分辨光电子能谱（ARPES）直接对比。**ABACUS** 作为一款优秀的国产开源密度泛函理论（DFT）软件，其基于数值原子轨道（LCAO）基组的特性，为构建紧束缚哈密顿量提供了天然的便利。结合 **PYATB**（Python Ab-initio Tight-Binding）后处理工具，我们可以高效地实现能带反折叠，还原出清晰的物理图像。

### 学习路线图
本教程采用循序渐进的模块化设计，旨在带你从零完成一次完整的反折叠分析：
1.  **物理基础（第一章）**：理解能带折叠的本质及反折叠的数学图像，明确标准工作流。
2.  **环境准备（第二章）**：攻克复杂的依赖库安装难题，建立基于容器化的稳定计算环境。
3.  **矩阵生成（第三章）**：掌握 ABACUS 的自洽计算，导出关键的哈密顿量与重叠矩阵。
4.  **后处理核心（第四章）**：深入 PYATB 的操作逻辑，完成从超胞到原胞投影的计算。
5.  **可视化分析（第五章）**：将计算数据转化为具备物理深度的有效能带图。

### 知识体系定位
本教程处于 ABACUS/DFT 知识体系的“高级应用”层级。它向上承接了基础的自洽场（SCF）计算与能带结构分析，向下延伸至拓扑材料表面态研究、磁性无序态模拟等前沿领域。相较于传统的平面波反折叠方法，本教程介绍的基于 LCAO 基组的方法在处理大规模体系时具有显著的效率优势。

### 前置要求
在开始学习前，建议读者具备以下基础知识：
-   **固体物理基础**：熟悉倒空间、布里渊区及能带理论。
-   **ABACUS 基础操作**：了解 INPUT、STRU、KPT 等核心输入文件的编写。
-   **Linux 与 Python 基础**：能够使用命令行执行程序，并具备基本的脚本运行能力。

---

# 第一章：能带反折叠的物理基础与应用场景

欢迎来到《ABACUS 实战教程》。在开始敲击键盘运行代码之前，作为开发者，我必须先带你理解“能带反折叠”（Band Unfolding）背后的物理图像。很多初学者能够跑通流程，但往往因为对物理概念理解偏差，导致 `m_matrix` 设置错误或对结果图产生误读。

此外，本章还将明确 ABACUS 配合 PYATB 进行反折叠计算的**标准工作流规范**。遵守这些规范将为你节省大量的排错时间。

---

## 1.1 超胞计算与能带折叠困境

### 为什么我们需要超胞（Supercell）？
在计算材料学中，我们经常需要处理破坏了完美晶体平移对称性的体系，例如：
- **掺杂（Doping）**：在硅中掺杂磷原子。
- **缺陷（Defects）**：空位或间隙原子。
- **合金（Alloys）**：如 SiGe 随机合金。
- **表面与界面**：为了模拟表面重构。

为了在周期性边界条件（Bloch 定理）下处理这些问题，我们通常构建**超胞（Supercell）**。即把 $N \times N \times N$ 个原胞（Primitive Cell, PC）拼成一个大盒子，然后在其中引入杂质或缺陷。

### 什么是能带折叠（Band Folding）？
根据固体物理的基本原理，实空间晶格矢量变大，倒空间布里渊区（Brillouin Zone, BZ）就会相应变小。
- **实空间**：超胞体积是原胞的 $M$ 倍（$V_{SC} = M \times V_{PC}$）。
- **倒空间**：超胞布里渊区体积是原胞的 $1/M$。

这就导致了一个严重的后果：**能带折叠**。
原胞中原本舒展清晰的能带，被迫“折叠”进了超胞那个极小的布里渊区内。这就好比把一张原本画在 A4 纸上的图，折叠了无数次塞进了一个火柴盒里。结果就是能带图变得极度密集、杂乱无章（Spaghetti bands），原本清晰的物理特征（如带隙性质、有效质量、能谷位置）完全无法辨认。

更重要的是，实验手段（如角分辨光电子能谱 **ARPES**）测量的是基于材料**有效原胞**的电子结构。如果我们直接拿超胞的一团乱麻的能带去和 ARPES 对比，是完全对不上号的。

---

## 1.2 反折叠（Unfolding）原理简介

为了解决上述问题，我们需要一种数学手段，将折叠进超胞布里渊区的波函数，“投影”回原胞的布里渊区。这就是**能带反折叠**。

### 核心思想：谱权重（Spectral Weight）
在超胞计算中，我们得到的本征态 $|\Psi_{KN}\rangle$（$K$ 为超胞 k 点，$N$ 为能带索引）虽然是超胞的本征态，但它依然包含了原胞的物理信息。

反折叠的核心思想是：将超胞的波函数 $|\Psi_{KN}\rangle$，用原胞的一组完备的布洛赫波函数 $|\psi_{kn}\rangle$（$k$ 为原胞 k 点，$n$ 为能带索引）进行展开。

我们需要计算一个物理量：**谱权重 $W_{Kn}(k)$**。
$$ P_{KN}(k) = \sum_{n} | \langle \psi_{kn} | \Psi_{KN} \rangle |^2 $$

这个权重的物理意义非常直观：**它代表了超胞的第 $N$ 条能带在 $K$ 点的量子态中，有多少成分是来自于原胞的 $k$ 点的。**

### 结果呈现：有效能带结构（EBS）
经过反折叠处理后，我们不再画单一的线条，而是画带有“深浅”或“粗细”的能带图：
- **权重高**：表示该能带保持了较好的原胞对称性（如未受杂质干扰的深层能级），在图中显示为深色/粗线。
- **权重低**：表示该能带被杂质散射严重，平移对称性破坏较大，在图中显示为浅色/模糊。
- **消失**：如果某处权重为零，说明该能带被完全破坏（如带隙中的杂质态）。

这与 ARPES 实验得到的谱函数强度图是直接对应的。

---

## 1.3 关键工作流规范（开发者建议）

在进入下一章的实操前，作为 ABACUS 和 PYATB 的使用者，你必须遵守以下**工程规范**。这是无数用户踩坑后总结出的血泪经验。

### 1. 环境配置警告
> **Warning**: PYATB 依赖特定的数学库（Eigen, MKL）和 Python 环境（mpi4py, pybind11）。环境搭建极易出错。

*   **推荐方案**：请务必使用官方验证过的镜像 `abacus:3.2.3`。
*   **原因**：不同版本的 GCC 和 Python 库极易导致 PYATB 编译失败或运行时出现 `undefined symbol` 错误。

### 2. 目录结构分离（Critical）
**绝对不要**在同一个文件夹下同时运行 ABACUS 的自洽计算（SCF）和 PYATB 的反折叠计算。

*   **ABACUS 输入文件**：名为 `INPUT` (全大写)。
*   **PYATB 输入文件**：名为 `Input` (首字母大写)。
*   **风险**：在不区分大小写的文件系统（如某些配置下的 macOS 或 Windows 挂载卷）中，或者当你手误使用 `rm` 命令时，这两个文件极易混淆或覆盖。

**推荐目录结构**：
```text
Work_Dir/
├── abacus_scf/       <-- 运行 ABACUS 自洽计算
│   ├── INPUT
│   ├── STRU
│   ├── KPT
│   └── ...
└── pyatb_run/        <-- 运行 PYATB 后处理
    ├── Input         <-- PYATB 的输入文件
    ├── STRU          <-- 必须与 abacus_scf 中完全一致
    ├── *.csr         <-- 从 abacus_scf 拷贝过来的稀疏矩阵文件
    └── ...
```

### 3. 必须手动获取费米能级（Manual Breakpoint）
这是全流程中唯一无法自动化的“断点”。
PYATB 需要知道体系的费米能级（Fermi Energy）来确定能量零点，但它目前**不会**自动去读取 ABACUS 的日志文件。

*   **操作**：自洽计算结束后，你必须执行：
    ```bash
    grep "EFERMI" OUT.suffix/running_scf.log
    ```
*   **动作**：将获取的数值（例如 `11.0829 eV`）**手动复制粘贴**到 PYATB 的 `Input` 文件中的 `fermi_energy` 字段。
*   **后果**：如果这一步填错或没填，画出来的能带图能量零点会偏移，导致费米面位置错误，所有物理分析将毫无意义。

### 4. 理解扩胞矩阵（m_matrix）
在 PYATB 的 `Input` 文件中，有一个参数叫 `m_matrix`。这**不是**随便抄一个案例就能用的。
它描述了**超胞晶格矢量**与**原胞晶格矢量**的变换关系：
$$ \vec{A}_{super} = M \times \vec{A}_{prim} $$

*   **案例**：如果你做的是 $2 \times 2 \times 2$ 的扩胞：
    ```text
    m_matrix  2 0 0  0 2 0  0 0 2
    ```
*   **注意**：如果你做的是 $3 \times 3 \times 1$ 的表面板层计算，矩阵必须相应修改为：
    ```text
    m_matrix  3 0 0  0 3 0  0 0 1
    ```
    *请务必根据你的实际 `STRU` 结构来设置此参数！*

### 5. 绘图脚本的微调
PYATB 运行结束后会生成一个 `plot_unfold.py` 脚本。
*   **默认行为**：它通常预设了一个较小的能量窗口（如 `[-4, 6]` eV）。
*   **操作**：你通常需要打开该 Python 脚本，找到 `energy_range` 变量，将其修改为你关心的范围（例如包含深层价带的 `[-14, 12]`），否则生成的 PDF 图中能带可能会被截断。

---

准备好了吗？理解了物理图像并记住了这些“避坑指南”后，请进入第二章，我们将开始实战操作。

# 第二章：计算环境搭建与依赖库安装

在进行能带反折叠（Band Unfolding）计算之前，构建一个稳定、兼容的计算环境是至关重要的。由于 **PYATB** (Python Ab-initio Tight-Binding) 是一个独立的 Python 后处理包，它并不包含在 ABACUS 的核心二进制文件中，且对底层的数学库（如 Eigen3、MKL）有严格的依赖。

本章将手把手带你完成环境搭建。请务必严格遵循以下步骤，特别是关于**镜像选择**和**文件目录管理**的建议，这能帮你规避 90% 以上的常见报错。

---

## 2.1 基础依赖库安装

### 2.1.1 镜像环境选择（Critical）

> **⚠️ 警告**：本教程强烈推荐使用 **`abacus:3.2.3`** 镜像环境。
>
> PYATB 需要编译 C++ 扩展模块，对编译器版本和系统库路径非常敏感。使用其他版本的镜像（特别是精简版或旧版）极易导致 `pybind11` 编译失败或运行时出现 `MKL link error`。

如果你使用的是 Bohrium 平台或其他容器化环境，请首先确认镜像版本。

### 2.1.2 Python 依赖与数学库安装

PYATB 的核心算法依赖于 `Eigen3` 线性代数库，并通过 `pybind11` 实现 Python 与 C++ 的接口。我们需要依次安装这些组件。

#### 1. 安装 Python 基础包
在终端中执行以下命令，安装 `pybind11`（用于 C++/Python 绑定）和 `mpi4py`（用于并行计算）：

```bash
pip install pybind11
pip install mpi4py
```

#### 2. 安装 Eigen3 库
`Eigen3` 是一个纯头文件的 C++ 库，不需要编译，但必须放置在编译器能找到的路径下。

```bash
# 进入临时目录
cd ~

# 下载 Eigen 3.4.0 源码包（推荐版本）
# 如果网络不佳，请使用国内镜像或手动上传
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz

# 解压
tar -xzf eigen-3.4.0.tar.gz

# 创建标准包含目录（如果不存在）
mkdir -p /usr/local/include

# 将 Eigen 库移动到系统标准路径
# 注意：PYATB 默认会在 /usr/local/include/eigen3 寻找头文件
mv eigen-3.4.0 /usr/local/include/eigen3
```

#### 3. 配置 MKL 环境变量（解决链接错误）
在 `abacus:3.2.3` 等环境中，直接运行 PYATB 可能会遇到 `undefined symbol: mkl_sparse_optimize_bsr_trsm_i8` 错误。这是由于 MKL 动态库加载顺序导致的。

请在终端执行（或写入 `~/.bashrc`）：

```bash
export LD_PRELOAD=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so
```

---

## 2.2 pyatb 编译与安装

环境准备就绪后，我们从源码编译安装 PYATB。

### 2.2.1 获取源码
```bash
cd ~
git clone https://github.com/pyatb/pyatb.git
cd pyatb
```

### 2.2.2 编译配置（siteconfig.py）
PYATB 使用 `siteconfig.py` 来定位依赖库。通常情况下，如果你按照 2.1 节将 Eigen3 放在了 `/usr/local/include/eigen3`，则**无需修改**该文件。

*如果你的 Eigen3 安装在自定义路径，请编辑 `siteconfig.py` 中的 `include_dirs` 参数。*

### 2.2.3 执行安装
在 `pyatb` 目录下执行：

```bash
pip install ./
```

安装完成后，可以通过以下命令验证：
```bash
python -c "import pyatb; print('PYATB installed successfully')"
```

---

## 2.3 关键工作流架构：流程分离（Critical）

在正式开始计算前，必须强调**文件目录管理的规范性**。

> **❌ 常见错误**：在同一个文件夹下既跑 ABACUS 自洽计算，又跑 PYATB 反折叠计算。
> **💥 后果**：ABACUS 的输入文件名为 `INPUT`，而 PYATB 的输入文件名为 `Input`。在 Linux 系统中虽然区分大小写，但极易造成用户混淆，且输出文件会相互覆盖，导致目录一片狼藉。

**✅ 推荐的最佳实践**：
请为每个反折叠任务建立如下目录结构：

```text
Project_Folder/
├── abacus_scf/          <-- 步骤1：运行 ABACUS 自洽计算
│   ├── INPUT
│   ├── STRU (超胞结构)
│   ├── KPT
│   └── ...
└── pyatb_run/           <-- 步骤2：运行 PYATB 后处理
    ├── Input            <-- 注意大小写！这是 PYATB 的输入文件
    ├── STRU (软链接或复制自 abacus_scf)
    ├── *.orb / *.upf    (需要复制过来)
    └── ...
```

---

## 2.4 核心参数与手动操作详解

在从 `abacus_scf` 转向 `pyatb_run` 的过程中，存在无法自动化的**手动断点**。如果忽略这些步骤，计算结果将完全错误。

### 1. 必须手动获取费米能级 (The "EFERMI" Check)
PYATB 需要知道体系的费米能级来对齐能带图的零点。这个数值必须从 ABACUS 的自洽计算日志中提取。

**操作步骤**：
1.  在 `abacus_scf` 目录下，自洽计算完成后，执行：
    ```bash
    grep "EFERMI" OUT.*/running_scf.log
    ```
    *输出示例：`EFERMI = 11.08291763890926 eV`*

2.  **手动复制**该数值，填入 `pyatb_run/Input` 文件的 `fermi_energy` 参数中：
    ```javascript
    // pyatb_run/Input
    INPUT_PARAMETERS
    {
        ...
        fermi_energy    11.08291763890926  // <--- 必须精确匹配
        fermi_energy_unit eV
        ...
    }
    ```

### 2. 必须传递稀疏矩阵文件
ABACUS 计算结束后，会生成包含哈密顿量信息的 `.csr` 文件。你需要将以下三个文件从 `abacus_scf/OUT.suffix/` 目录复制到 `pyatb_run/` 目录：
*   `data-HR-sparse_SPIN0.csr` (哈密顿量矩阵)
*   `data-SR-sparse_SPIN0.csr` (重叠矩阵)
*   `data-rR-sparse.csr` (位置算符矩阵)

### 3. 理解并修改扩胞矩阵 (m_matrix)
在 `pyatb_run/Input` 中，`m_matrix` 参数定义了**超胞基矢量与原胞基矢量之间的变换关系**。

> **⚠️ 警告**：不要盲目照抄案例中的矩阵！

*   **物理含义**：如果超胞基矢 $A$ 与原胞基矢 $a$ 的关系为 $A = M \cdot a$，则 `m_matrix` 对应矩阵 $M$ 的转置或其展开形式（具体取决于定义，通常为整数矩阵）。
*   **案例分析**：
    如果你的超胞是简单的 $2 \times 2 \times 2$ 扩胞（即 $A_1=2a_1, A_2=2a_2, A_3=2a_3$），矩阵通常是对角线上为 2。
    **但是**，如果你做的是 $3 \times 3 \times 1$ 的表面板层扩胞，你需要相应修改矩阵，例如：
    ```javascript
    // 针对 3x3x1 扩胞的示例（仅作示意，需根据实际变换填写）
    m_matrix  3  0  0
              0  3  0
              0  0  1
    ```
    *注：PYATB 输入通常将其写为一行，如 `3 0 0 0 3 0 0 0 1`。请务必根据你实际构建超胞的方式来设置此参数。*

### 4. 绘图脚本的能量范围修正
PYATB 运行结束后会生成 `plot_unfold.py`。这是一个自动生成的 Python 脚本。
默认生成的 `energy_range` 往往不符合你的研究需求（例如默认范围太窄或太宽）。

**操作建议**：
打开 `plot_unfold.py`，找到以下行并手动修改：
```python
# 修改前
# energy_range = [-4, 6] 

# 修改后 (根据你的体系，例如半导体通常关心带边附近)
energy_range = [-14, 12] 
```
修改后再次运行 `python plot_unfold.py` 即可得到清晰的能带反折叠图。

# 第三章：ABACUS 自洽计算与矩阵生成

在上一章中，我们完成了能带反折叠的理论准备。本章我们将进入实战的核心环节：使用 ABACUS 进行自洽场（SCF）计算。

在此流程中，ABACUS 的角色不仅仅是计算能量，更是一个**“哈密顿量生成器”**。我们的目标是通过 LCAO（原子轨道线性组合）基组，构建并输出体系的稀疏矩阵（Hamiltonian, Overlap, Position），这些矩阵是后续 PYATB 进行能带反折叠分析的原材料。

> **⚠️ 环境配置警告**
> 本教程涉及编译安装 PYATB，环境依赖极其敏感。强烈建议使用官方推荐的镜像 **`abacus:3.2.3`**。使用其他版本（特别是较旧版本或不兼容的 GCC 环境）极易导致后续 PYATB 安装失败或运行时出现 `undefined symbol: mkl_...` 等数学库错误。

---

## 3.0 工作目录规划（关键）

在开始计算前，必须强调一个工程习惯：**流程分离**。

由于 ABACUS 的输入文件名为 `INPUT`，而 PYATB 的输入文件通常命名为 `Input`（或 `input.json`），在同一目录下操作极易导致文件覆盖或混淆。此外，保持计算目录（生成数据）与分析目录（处理数据）的分离，有助于数据的长期管理。

**建议的文件结构：**

```text
Project_Root/
├── abacus_run/          <-- 本章工作目录
│   ├── INPUT            # ABACUS 控制参数
│   ├── STRU             # 超胞结构文件
│   ├── KPT              # K点设置
│   ├── *.upf & *.orb    # 赝势与轨道文件
│   └── OUT.suffix/      # 计算输出目录
└── pyatb_run/           <-- 下一章工作目录
    ├── Input            # PYATB 控制参数
    └── (待复制的矩阵文件)
```

请立即建立这两个文件夹，并进入 `abacus_run` 目录开始本章操作。

---

## 3.1 输入参数设置 (INPUT)

`INPUT` 文件控制着 ABACUS 的计算行为。为了适配 PYATB，我们需要开启特定的矩阵输出开关。

### 核心参数详解

1.  **`basis_type`: `lcao`**
    *   **必须设置**。能带反折叠基于紧束缚模型（Tight-Binding），必须使用原子轨道基组。平面波（`pw`）模式无法生成所需的轨道矩阵。

2.  **`calculation`: `scf`**
    *   进行自洽计算以获得收敛的电子密度和波函数。

3.  **`out_mat_hs2`: `1`**
    *   **关键开关**。设置为 `1` 时，ABACUS 会在计算结束后输出稀疏格式的哈密顿量矩阵（$H$）和重叠矩阵（$S$）。这是 PYATB 读取电子结构的核心数据源。

4.  **`out_mat_r`: `1`**
    *   **关键开关**。输出位置算符矩阵（$r$）。这对于计算偶极矩、光学性质或处理某些涉及贝里相位的计算是必要的。

5.  **`nbands`**
    *   **需特别注意**。默认情况下，ABACUS 仅计算被占据的能带（加上少量空带）。但在能带反折叠分析中，我们通常也关心导带（Conduction Bands）的结构。
    *   **建议**：显式设置一个较大的 `nbands` 值（例如总电子数的一半 + 20~50，或者更多），确保覆盖你感兴趣的高能级范围。

### INPUT 文件示例

新建 `INPUT` 文件并写入以下内容（以 GeC 体系为例，请根据实际体系调整）：

```bash
INPUT_PARAMETERS
# System variables
suffix                GeC_scf
ntype                 2
calculation           scf
symmetry              1         # 开启对称性分析可加速计算
init_chg              atomic

# Electronic structure
basis_type            lcao      # 【必须】使用 LCAO 基组
ks_solver             genelpa   # 推荐使用的对角化求解器
nspin                 1         # 1为非磁性，2为自旋极化
ecutwfc               50        # 截断能 (Ry)
scf_nmax              100       # 最大自洽步数
scf_thr               1e-7      # 自洽收敛阈值
nbands                100       # 【重要】确保包含足够的空带

# Variables related to output information for PYATB
out_mat_hs2           1         # 【关键】输出 H 和 S 稀疏矩阵
out_mat_r             1         # 【关键】输出位置矩阵
out_chg               0         # 可选：是否输出电荷密度
out_band              0         # scf计算通常不需要输出能带文件
```

---

## 3.2 结构与基组文件准备

在此步骤中，我们需要准备超胞（Supercell）的结构文件。

### 1. 结构文件 (STRU)
**注意**：这里的 `STRU` 必须是**扩胞后**的结构（Large Cell），而不是原胞（Primitive Cell）。能带反折叠的物理本质，就是将这个超胞的能带“投影”回原胞的布里渊区。

*   **逻辑检查**：如果你计划研究 3x3x1 的超胞掺杂效应，这里的 `STRU` 就必须包含 3x3x1 倍数的原子。
*   **矩阵变换预警**：请记住你从原胞扩胞到超胞使用的变换矩阵（例如 2x2x2 或 3x3x1）。在下一章配置 PYATB 的 `m_matrix` 参数时，必须与此处的结构严格对应，**切勿照抄案例中的矩阵**。

### 2. 赝势与轨道文件 (.upf & .orb)
确保 `INPUT` 文件中的 `pseudo_dir` 和 `orbital_dir` 指向正确的路径。
*   `.upf`: 赝势文件。
*   `.orb`: 数值原子轨道文件（LCAO 专用）。

### 3. K点设置 (KPT)
对于超胞计算，由于实空间晶格变大，倒空间布里渊区变小，因此 K 点网格可以相应稀疏一些（例如原胞用 12x12x12，2x2x2 超胞可能只需要 6x6x6）。

```text
K_POINTS
0
Gamma
2 2 2 0 0 0
```

---

## 3.3 运行自洽计算与结果检查

准备好所有文件后，提交任务运行 ABACUS。

```bash
# 假设使用 16 核并行
mpirun -n 16 abacus
```

### 1. 检查运行状态
计算完成后，首先检查屏幕输出或日志文件，确保 `SCF convergence reached`，即自洽收敛。

### 2. 验证矩阵输出
进入输出目录（例如 `OUT.GeC_scf/`），检查是否生成了以下关键文件：

*   `data-HR-sparse_SPIN0.csr`: 哈密顿量矩阵（稀疏格式）
*   `data-SR-sparse_SPIN0.csr`: 重叠矩阵
*   `data-rR-sparse.csr`: 位置算符矩阵

如果缺少这些文件，请回头检查 `INPUT` 中 `out_mat_hs2` 和 `out_mat_r` 是否正确设置为 `1`。

### 3. 获取费米能级（关键手动步骤）

这是一个**非自动化**的断点，至关重要。PYATB 需要知道体系的费米能级位置以对齐能量零点，但它不会自动从矩阵文件中读取。

你需要查看 ABACUS 的日志文件（通常是 `running_scf.log` 或标准输出），找到 `EFERMI` 的值。

```bash
# 在输出目录执行 grep 命令
grep "EFERMI" running_scf.log
```

**输出示例：**
```text
EFERMI = 6.54321 eV
```

> **📝 记录下来**：请务必记录这个数值（例如 `6.54321`）。在下一章编写 PYATB 的 `Input` 文件时，你需要手动将这个值填入 `fermi_energy` 字段。如果数值错误，最终能带图的费米面将发生偏移。

---

## 3.4 为 PYATB 准备数据

为了保持目录整洁，我们将生成的核心文件转移到 `pyatb_run` 文件夹中，为下一章做准备。

```bash
# 假设当前在 abacus_run/OUT.GeC_scf/ 目录下
# 且 pyatb_run 文件夹位于上两级目录

cp data-HR-sparse_SPIN0.csr ../../pyatb_run/
cp data-SR-sparse_SPIN0.csr ../../pyatb_run/
cp data-rR-sparse.csr       ../../pyatb_run/

# 此外，PYATB 也需要读取结构文件和轨道/赝势信息
cp ../STRU                  ../../pyatb_run/
# 建议将赝势和轨道文件也软链接或复制过去，或者在 pyatb Input 中指定绝对路径
```

至此，我们已经成功利用 ABACUS 生成了包含体系全部量子力学信息的稀疏矩阵。下一章，我们将切换到 PYATB，利用这些矩阵“反推”出有效能带结构。

# 第四章：pyatb 后处理与能带反折叠计算

在上一章中，我们已经成功完成了超胞（Supercell）的自洽计算（SCF），并生成了构建紧束缚模型所需的稀疏矩阵文件。然而，超胞的能带图由于布里渊区折叠（Band Folding），通常杂乱无章，无法直接与实验（如 ARPES）或原胞能带进行对比。

本章将是整个流程的“核心枢纽”。我们将使用 `pyatb` 软件，将 ABACUS 计算得到的电子结构信息“反折叠”（Unfold）回原胞的布里渊区。

> **⚠️ 环境预警**：
> 本章涉及 `pyatb` 的编译与运行，强烈建议使用官方推荐镜像 **`abacus:3.2.3`**。其他环境极易出现 `pybind11` 或 `MKL` 库链接错误。

---

## 4.1 工作目录规划与数据迁移

在进行后处理之前，**强烈建议**将 ABACUS 的计算目录与 `pyatb` 的运行目录物理隔离。这不仅是为了保持整洁，更是为了防止输入文件（ABACUS 的 `INPUT` 与 pyatb 的 `Input`）发生混淆。

### 4.1.1 建立独立目录
假设你的 ABACUS 自洽计算目录为 `abacus_run`，请在同级目录下新建一个名为 `pyatb_run` 的文件夹。

```bash
# 假设当前在项目根目录
mkdir pyatb_run
```

### 4.1.2 关键文件迁移
我们需要将 ABACUS 计算生成的**稀疏矩阵文件**和**结构文件**复制到新目录。

需要复制的文件清单如下：
1.  **`STRU`**: 结构文件（pyatb 需要读取晶格信息）。
2.  **`data-HR-sparse_SPIN0.csr`**: 哈密顿量矩阵（若开启自旋可能有 SPIN1）。
3.  **`data-SR-sparse_SPIN0.csr`**: 重叠矩阵。
4.  **`data-rR-sparse.csr`**: 位置算符偶极矩阵。

> **注意**：上述 `.csr` 文件通常位于 ABACUS 输出目录（如 `OUT.suffix/`）下。

**操作示例**：
```bash
# 进入 ABACUS 输出目录
cd abacus_run/OUT.GeC/

# 复制文件到 pyatb 目录
cp data-HR-sparse_SPIN0.csr ../../pyatb_run/
cp data-SR-sparse_SPIN0.csr ../../pyatb_run/
cp data-rR-sparse.csr       ../../pyatb_run/

# 别忘了复制 STRU 文件（通常在 abacus_run 根目录）
cd ../
cp STRU ../pyatb_run/
```

完成后的 `pyatb_run` 目录结构应如下所示：
```text
pyatb_run/
├── STRU
├── data-HR-sparse_SPIN0.csr
├── data-SR-sparse_SPIN0.csr
└── data-rR-sparse.csr
```

---

## 4.2 关键步骤：费米能级的手动提取

这是流程中**唯一无法自动化且至关重要**的一步。

ABACUS 计算得到的费米能级（Fermi Energy）存储在日志文件中，而 `pyatb` 无法自动读取该值。如果这一步数值填错，最终能带图的能量零点（$E-E_f$）将完全错误。

### 操作方法
请回到 ABACUS 的运行目录，从 `running_scf.log`（或屏幕输出日志）中提取 `EFERMI`。

```bash
# 在 ABACUS 输出目录下执行
grep "EFERMI" running_scf.log
```

**输出示例**：
```text
EFERMI = 11.08291763890926 eV
```

> **📝 记录下来**：请复制这个数值（例如 `11.08291763890926`），我们将在下一节的输入文件中使用它。

---

## 4.3 pyatb 输入文件 (Input) 详解

在 `pyatb_run` 目录下，我们需要创建一个名为 `Input` 的文件（注意首字母大写，区分于 ABACUS 的 `INPUT`）。

### 4.3.1 Input 文件完整示例

```json
INPUT_PARAMETERS
{
    "nspin": 1,
    "package": "ABACUS",
    "fermi_energy": 11.08291763890926,  
    "fermi_energy_unit": "eV",
    "HR_route": "data-HR-sparse_SPIN0.csr",
    "SR_route": "data-SR-sparse_SPIN0.csr",
    "rR_route": "data-rR-sparse.csr",
    "HR_unit": "Ry",
    "rR_unit": "Bohr"
}

LATTICE
{
    "lattice_constant": 1.889726,
    "lattice_constant_unit": "Bohr",
    "lattice_vector": [
        11.3149995804, 0, 0,
        0, 11.3149995804, 0,
        0, 0, 11.3149995804
    ]
}

BANDUNFOLDING
{
    "stru_file": "STRU",
    "ecut": 50,
    "band_range": [1, 500],
    "m_matrix": [
        -2, 2, 2,
        2, -2, 2,
        2, 2, -2
    ],
    "kpoint_mode": "line",
    "kpoint_num": 50,
    "high_symmetry_kpoint": [
        0.5, 0.5, 0.5, 50,  # L
        0.0, 0.0, 0.0, 50,  # Gamma
        0.5, 0.0, 0.5, 50,  # X
        0.0, 0.0, 0.0, 1    # Gamma
    ]
}
```

### 4.3.2 核心参数深度解析

#### 1. 基础参数 (`INPUT_PARAMETERS`)
*   **`package`**: 必须指定为 `ABACUS`。
*   **`fermi_energy`**: **此处填入 4.2 节中提取的数值**。
*   **`*_route`**: 指定矩阵文件的路径。由于我们已经将文件复制到了当前目录，直接写文件名即可。

#### 2. 晶格信息 (`LATTICE`)
*   建议直接从 ABACUS 的 `STRU` 文件或 `running_scf.log` 中复制晶格常数和基矢量。注意单位通常为 Bohr。

#### 3. 扩胞矩阵 (`m_matrix`) —— **易错点**
这是能带反折叠的物理核心。`m_matrix` 定义了**超胞基矢量 ($A_{super}$)** 与 **原胞基矢量 ($A_{prim}$)** 之间的变换关系。
关系式为：
$$ A_{super} = A_{prim} \times M $$

*   **不要照抄案例**：案例中的 `[-2, 2, 2, ...]` 是针对特定晶体结构的 $2\times2\times2$ 扩胞变换。
*   **如何设置**：
    *   如果你做的是简单的 $2\times2\times2$ 正交扩胞，矩阵通常是 `2 0 0 0 2 0 0 0 2`（即对角矩阵）。
    *   如果你使用了 `atomkit` 或其他工具生成了旋转扩胞，必须使用生成超胞时对应的变换矩阵。
    *   **格式**：在 Input 文件中按行优先顺序填入 9 个数字。

#### 4. 能带范围 (`band_range`)
*   指定需要计算并反折叠的能带索引范围（如 `[1, 500]`）。确保该范围覆盖了费米面附近的价带和导带。

---

## 4.4 运行 pyatb 与 MKL 环境配置

配置好 `Input` 文件后，就可以运行计算了。但在执行命令前，必须解决一个常见的 Intel MKL 库链接问题。

### 4.4.1 解决 MKL `undefined symbol` 报错
在 `abacus:3.2.3` 容器或大多数 Linux 环境中，直接运行可能会遇到如下报错：
`undefined symbol: mkl_sparse_optimize_bsr_trsm_i8`

这是由于系统预加载库顺序导致的。**解决方法**是在运行前设置 `LD_PRELOAD` 环境变量。

请在终端执行以下命令（或将其写入提交脚本）：

```bash
export LD_PRELOAD=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so
```

> **提示**：如果你使用的是不同的镜像或自定义环境，上述路径可能需要根据实际 MKL 安装位置进行调整。

### 4.4.2 执行计算
设置好环境变量后，使用 MPI 并行运行 pyatb：

```bash
# -n 指定核数，建议与物理核数一致
mpirun -n 16 pyatb
```

计算通常很快（几分钟内）。运行结束后，当前目录下会生成一个 `Out` 文件夹。

---

## 4.5 结果可视化

进入输出目录查看结果：

```bash
cd Out/Bandunfolding
ls
```

你会看到以下关键文件：
*   `spectral_weight.dat`: 包含反折叠后的能带权重数据。
*   `plot_unfold.py`: 自动生成的绘图脚本。

### 4.5.1 调整绘图范围
默认生成的 `plot_unfold.py` 脚本中，能量显示范围 (`energy_range`) 可能不符合你的需求（例如默认可能是 `[-4, 6]`）。

**必须手动修改该脚本**：
1.  打开 `plot_unfold.py`。
2.  找到 `energy_range` 变量。
3.  根据你的体系修改数值。例如，为了看清深层能级，可以改为：
    ```python
    # 修改前
    # energy_range = [-4, 6]
    
    # 修改后
    energy_range = [-14, 12]
    ```

### 4.5.2 生成能带图
运行修改后的脚本：

```bash
python plot_unfold.py
```

成功运行后，将生成 `unfold.pdf`。这张图展示了超胞能带投影回原胞布里渊区后的“有效能带结构”（Effective Band Structure），其中颜色的深浅代表了谱权重（Spectral Weight），清晰地还原了被杂质或缺陷破坏的能带色散关系。

# 第五章：结果可视化与分析

在前几章中，我们已经成功完成了 ABACUS 的自洽计算（SCF），并利用 PYATB 完成了能带反折叠（Band Unfolding）的核心计算步骤。此时，你的 `pyatb` 工作目录下应该已经生成了 `spectral_weight.dat` 和 `plot_unfold.py` 等关键文件。

本章将指导你如何将这些抽象的数值数据，转化为直观的“有效能带图”（Effective Band Structure），并学会从物理角度解读这些结果。

> **⚠️ 环境与依赖警告**
> 本章涉及 Python 绘图脚本的运行。如果你使用的是容器环境，强烈推荐使用 **`abacus:3.2.3`** 镜像，该镜像已预装了 PYATB 所需的 `matplotlib`、`numpy` 等依赖库。其他环境可能需要手动安装依赖。

---

## 5.1 绘图脚本参数调整

PYATB 在运行结束后，会自动生成一个名为 `plot_unfold.py` 的 Python 脚本。这是一个基于 Matplotlib 的绘图模板。直接运行它通常能得到一张图，但为了获得最佳的视觉效果，我们需要根据物理体系的实际情况调整能量窗口。

### 5.1.1 理解 `energy_range`

打开 `plot_unfold.py` 文件，找到如下参数设置：

```python
# plot_unfold.py 片段
# ...
energy_range = [-4, 6]  # 默认值可能不适合所有体系
# ...
```

这里的 `energy_range` 定义了 Y 轴（能量轴）的显示范围，单位通常为 eV。
*   **零点基准**：这里的 $E=0$ 对应于你在 `Input` 文件中手动填入的 `fermi_energy`。
*   **调整策略**：
    *   对于**绝缘体/半导体**：通常关注价带顶（VBM）和导带底（CBM）附近的区域，例如 `[-4, 4]` 或 `[-6, 6]`。
    *   对于**金属**：关注费米面附近的电子态，范围可以更小。
    *   **案例修正**：在 GeC 案例中，为了完整展示深层能级和高能激发态，我们将范围扩大。

### 5.1.2 修改实战

你可以使用文本编辑器（如 `vim` 或 `nano`）手动修改，也可以使用 `sed` 命令快速替换。

**方式一：手动修改**
将 `energy_range` 的数值改为你关心的区间，例如 `[-14, 12]`：
```python
energy_range = [-14, 12]
```

**方式二：命令行流编辑（推荐脚本作业使用）**
```bash
# 将默认的 [-4, 6] 替换为 [-14, 12]
sed -i 's/energy_range = \[-4, 6\]/energy_range = [-14, 12]/g' plot_unfold.py
```

### 5.1.3 运行绘图

修改完成后，执行脚本生成 PDF 图像：

```bash
python plot_unfold.py
```

如果运行成功，目录下会生成 `unfold.pdf`。

---

## 5.2 结果解读

打开生成的 `unfold.pdf`，你会看到一张类似能带图的图像，但它包含了比普通能带图更丰富的信息。

### 5.2.1 谱权重（Spectral Weight）的物理意义

在反折叠能带图中，线条不再是单一的黑色实线，而是具有**颜色深浅**或**粗细**变化的“云图”。

*   **颜色/深浅**：代表 **谱权重 (Spectral Weight)** 的大小。
    *   $$P_{Km}(k_i)$$：表示超胞中的本征态 $|K, m\rangle$ 投影到原胞布洛赫态 $|k_i\rangle$ 上的概率。
*   **深色/高亮线条**：代表谱权重接近 1。这意味该能带保持了良好的原胞平移对称性，电子态是离域的（Bloch-like）。
*   **浅色/模糊线条（影子带）**：代表谱权重较小。这些通常是由于超胞引入的微扰（如缺陷、替位原子）导致对称性破缺，使得原胞的 $k$ 点与 $k+G$ 点发生了耦合，产生的“折叠”假象。
*   **能带展宽（Smearing）**：如果在某些能量区间线条变得模糊不清，说明该能量处的电子态受到了强烈的散射（例如无序合金或高浓度掺杂），动量 $k$ 不再是一个好的量子数。

### 5.2.2 典型特征识别

1.  **主能带 (Primitive Bands)**：图中清晰、连续且颜色最深的轨迹。这对应于未受扰动的理想晶体的能带结构。
2.  **能隙 (Band Gap)**：观察费米能级（$E=0$）附近。如果 $E=0$ 处没有谱权重分布，且上下有清晰的带边，则为半导体/绝缘体。
3.  **缺陷态 (Defect States)**：如果在禁带中间出现了平直的、色散很小（flat band）的线条，这往往对应于局域的缺陷能级。

---

## 附录：常见问题与进阶建议

在 ABACUS + PYATB 的实战流程中，以下环节是“深坑”高发区，请务必仔细核对。

### A.1 目录结构与文件名陷阱 (Critical)

**强烈建议**将 ABACUS 自洽计算与 PYATB 后处理分在两个不同文件夹进行！

*   **原因**：ABACUS 的输入文件名为 `INPUT` (全大写)，而 PYATB 的输入文件名为 `Input` (首字母大写)。在 Linux 系统下虽然区分大小写，但极易在 Tab 补全或脚本编写时混淆。此外，保持目录整洁有助于管理大量输出文件。
*   **推荐结构**：
    ```text
    Project_Root/
    ├── abacus_run/       # 运行 SCF
    │   ├── INPUT
    │   ├── STRU
    │   ├── KPT
    │   └── OUT.suffix/   # ABACUS 输出目录
    └── pyatb_run/        # 运行 Unfolding
        ├── Input         # PYATB 输入
        ├── STRU          # 复制自 abacus_run
        ├── data-HR-sparse_SPIN0.csr  # 复制自 OUT.suffix
        ├── data-SR-sparse_SPIN0.csr  # 复制自 OUT.suffix
        └── data-rR-sparse.csr        # 复制自 OUT.suffix
    ```

### A.2 费米能级的获取 (Manual Checkpoint)

这是全流程中唯一无法自动化的**断点**，必须手动操作。

1.  **获取数值**：在 `abacus_run` 目录下，从 SCF 日志中提取费米能级。
    ```bash
    grep "EFERMI" OUT.suffix/running_scf.log
    # 输出示例: EFERMI = 11.08291763890926 eV
    ```
2.  **填入 Input**：将上述数值**精确**填入 `pyatb_run/Input` 文件中：
    ```json
    INPUT_PARAMETERS
    {
        ...
        fermi_energy    11.08291763890926
        fermi_energy_unit eV
        ...
    }
    ```
    > **警告**：如果此数值填写错误（例如填了 0），最终画出的能带图 $E=0$ 位置将不是费米面，导致物理分析完全错误。

### A.3 扩胞矩阵 (m_matrix) 的正确设置

在 `Input` 文件中，`m_matrix` 参数定义了**超胞晶格矢量**与**原胞晶格矢量**之间的变换关系。

$$ \vec{A}_{super} = \vec{A}_{prim} \times M $$

*   **不要盲目照抄案例**：案例中的矩阵（如 `-2 2 2 ...`）是针对特定晶体结构（如 FCC 原胞转立方超胞）的。
*   **简单扩胞**：如果你只是简单地将原胞沿 $x, y, z$ 方向分别扩大 $N_x, N_y, N_z$ 倍（例如 `3x3x1` 超胞），则矩阵应为对角阵：
    ```text
    m_matrix  3 0 0  0 3 0  0 0 1
    ```
    如果是 `2x2x2` 扩胞，则是：
    ```text
    m_matrix  2 0 0  0 2 0  0 0 2
    ```

### A.4 磁性体系的处理 (Spin Polarization)

如果你的体系是磁性的（在 ABACUS `INPUT` 中设置了 `nspin 2`），ABACUS 会输出两套哈密顿量文件：
*   `data-HR-sparse_SPIN0.csr` (自旋向上 / Up)
*   `data-HR-sparse_SPIN1.csr` (自旋向下 / Down)

**处理方法**：
你需要分别运行两次 PYATB（建议分两个子文件夹），分别在 `Input` 文件中指定对应的 `HR_route` 和 `SR_route`，最后分别画出 Spin Up 和 Spin Down 的能带图进行对比。

### A.5 常见报错：MKL 链接错误

如果在运行 `pyatb` 命令时遇到类似 `undefined symbol: mkl_sparse_optimize_bsr_trsm_i8` 的错误，这是因为环境变量未正确加载 Intel MKL 库。

**解决方案**：
在运行前执行以下命令（可视情况加入 `.bashrc`）：
```bash
export LD_PRELOAD='/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_def.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_avx512.so.2:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so'
```
*(注：具体路径可能随镜像版本不同而微调，上述路径适用于 `abacus:3.2.3` 镜像)*

---

## 附录：进阶学习与调试指南

### 进阶探索建议
完成本教程的实战练习只是第一步，能带反折叠技术在实际科研中有着广阔的应用空间，鼓励读者在以下方向进一步探索：
-   **复杂体系研究**：尝试将本流程应用于掺杂半导体、高熵合金或拓扑绝缘体的表面态研究。参考资料中提到的 MnBi2Te4 顺磁态模拟是一个极佳的进阶案例。
-   **非共线磁性**：在涉及自旋-轨道耦合（SOC）的体系中，反折叠过程将涉及非共线磁矩的处理，这需要更深入地理解哈密顿矩阵的结构。
-   **算法深入**：阅读戴祖建、金敢等人的原始文献（*Computational Materials Science 213 (2022) 111656*），理解基于数值原子轨道投影到平面波子空间的底层算法逻辑。

### 通用调试建议 (Troubleshooting)
在实践过程中，你可能会遇到各种报错，以下是通用的排查思路：
1.  **环境一致性**：90% 的 PYATB 报错源于底层数学库（MKL/Eigen3）版本不匹配。强烈建议始终使用官方推荐的 Docker 镜像 `abacus:3.2.3`。
2.  **矩阵文件完整性**：确保 ABACUS 输出的 `SR`、`HR`、`rR` 文件完整且未在传输过程中损坏。LCAO 计算时必须开启对应的矩阵输出参数。
3.  **K 点路径匹配**：反折叠计算中的 K 点路径必须与你想要投影的原胞高对称路径严格对应，否则会导致谱权重分布异常。
4.  **内存管理**：对于超大规模超胞，PYATB 的矩阵处理可能消耗大量内存，建议在计算节点上监控内存使用情况。

### 官方资源与支持
-   **ABACUS 官方文档**：[https://abacus.deepmodeling.com/](https://abacus.deepmodeling.com/) 获取最新的参数说明与更新动态。
-   **PYATB 开源仓库**：关注 GitHub 上的相关项目，获取最新的后处理脚本示例。
-   **社区支持**：若在 Bohrium 平台上运行遇到问题，可联系 [bohrium@dp.tech](mailto:bohrium@dp.tech) 或在 DeepModeling 社区发帖讨论。

祝你在材料模拟的微观世界中探索愉快！
