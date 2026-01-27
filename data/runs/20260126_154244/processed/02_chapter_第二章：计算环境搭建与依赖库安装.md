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