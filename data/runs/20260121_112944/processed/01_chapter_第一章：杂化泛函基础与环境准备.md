# 第一章：杂化泛函基础与环境准备

欢迎来到《ABACUS 实战教程》。作为第一章，我们不急于马上运行庞大的计算任务，而是要先夯实物理基础并确保你的“武器”——ABACUS 编译版本——已经为杂化泛函做好了准备。

在我的开发和教学经验中，超过 50% 的用户在初次尝试杂化泛函计算时，遇到的不是收敛问题，而是直接的“未定义符号”或“缺少库”的报错。本章将帮助你理解我们为什么要使用昂贵的杂化泛函，并手把手教你排查和解决环境依赖问题。

---

## 1.1 物理背景与适用场景

### 1.1.1 为什么我们需要杂化泛函？

在标准的密度泛函理论（DFT）中，无论是局域密度近似（LDA）还是广义梯度近似（GGA，如 PBE），都存在一个著名的病态问题：**自相互作用误差（Self-Interaction Error, SIE）**。

简单来说，在经典的 Hartree 势中，电子会错误地与自己的电荷密度发生相互作用。虽然交换关联泛函（XC）试图抵消这一项，但 LDA/GGA 无法完全抵消它。这种残留的排斥力导致电子过度离域化（delocalization）。

**后果是严重的：**
1.  **带隙低估（Bandgap Underestimation）**：这是半导体物理中最痛的点。例如，PBE 计算出的硅（Si）带隙可能只有实验值的一半，甚至将某些窄带隙半导体错误地描述为金属。
2.  **强关联体系失效**：对于过渡金属氧化物（如 NiO, FeO）或 f 电子体系，电子高度局域化，LDA/GGA 往往无法正确描述其磁性基态或绝缘性质。

### 1.1.2 杂化泛函的解决方案

杂化泛函（Hybrid Functionals）的思路非常直接：既然 Hartree-Fock (HF) 理论包含精确的交换能（Exact Exchange, EXX）且完全没有自相互作用误差（但缺少关联能），而 DFT 包含关联能但有 SIE，为什么不把它们混合起来呢？

通用的杂化泛函能量公式可以写为：
$$ E_{xc}^{hybrid} = \alpha E_{x}^{HF} + (1-\alpha)E_{x}^{DFT} + E_{c}^{DFT} $$

其中 $\alpha$ 是混合系数（Mixing parameter）。
-   **PBE0**：$\alpha = 0.25$。即混合了 25% 的 HF 交换能。
-   **HSE06**：这是固体物理中最常用的泛函。为了加速计算，它使用屏蔽库伦势（Screened Coulomb Potential），只在短程引入 HF 交换项，长程使用 DFT。这极大地降低了计算量，同时保留了修正带隙的能力。

### 1.1.3 ABACUS 的独特优势

ABACUS 在处理杂化泛函时具有独特的优势，特别是结合 **数值原子轨道（LCAO）** 基组时。相比于平面波基组，LCAO 基组天然具有局域性，这使得精确交换项（EXX）的计算可以通过 **LibRI**（Resolution of Identity）库进行极其高效的处理。

**适用场景总结**：
*   **必须使用**：计算半导体/绝缘体的准确带隙、能带边位置（Band edge alignment）。
*   **推荐使用**：缺陷形成能、极化子（Polaron）研究、强关联过渡金属氧化物。
*   **慎重使用**：大体系（>500 原子）的动力学模拟（除非你有极其充裕的算力，因为计算量通常是 PBE 的 10-100 倍）。

---

## 1.2 编译依赖与安装检查

这是本章最关键的部分。ABACUS 的杂化泛函功能并非“开箱即用”的默认配置，它依赖于几个特定的外部数学库。

### 1.2.1 核心依赖库：“三驾马车”

要在 ABACUS (LCAO模式) 中高效运行杂化泛函，编译时必须链接以下三个库：

1.  **Libxc**：
    *   **作用**：提供各种交换关联泛函的定义（如 HSE06, PBE0 的具体数学形式）。
    *   **必要性**：如果没有它，ABACUS 根本不知道 `HSE` 是什么。
2.  **LibRI** (Resolution of Identity Library)：
    *   **作用**：这是 ABACUS 杂化泛函计算速度的核心。它专门用于处理 LCAO 基组下的四中心双电子积分（ERI）。
    *   **必要性**：在 LCAO 模式下计算 EXX，必须依赖此库。
3.  **LibComm**：
    *   **作用**：辅助的基础通信库，支持 LibRI 的运行。

### 1.2.2 快速自检：你的 ABACUS 能跑杂化泛函吗？

在开始复杂的计算前，请执行以下简单的“冒烟测试”。

**步骤 1：准备一个最小化的 INPUT 文件**
创建一个名为 `INPUT` 的文件，写入以下内容：
```fortran
INPUT_PARAMETERS
# 基础参数
calculation     scf
basis_type      lcao
ks_solver       genelpa  # 或者 scalapack_gvx

# 关键测试参数
dft_functional  hse      # 指定使用 HSE 杂化泛函
```
*(注：你还需要配套的 `STRU` 和 `KPT` 以及赝势/轨道文件才能运行，这里假设你已有基本的 Si 或 H2O 测试算例)*

**步骤 2：尝试运行**
执行 ABACUS 可执行文件。

**步骤 3：观察报错**

*   **情况 A：运行成功或开始迭代**
    *   🎉 恭喜，你的版本已经支持杂化泛函。

*   **情况 B：报错 "Unrecognized exchange-correlation functional"**
    *   ❌ **原因**：编译时未链接 **Libxc**。
    *   *错误示例*：
        ```text
        Unrecognized exchange-correlation functional 'HSE'.
        ```

*   **情况 C：报错 "compile with libri..."**
    *   ❌ **原因**：编译时未链接 **LibRI**。这是 LCAO 模式下最常见的错误。
    *   *错误示例*：
        ```text
        compile with libri to use hybrid functional in lcao basis
        ```

### 1.2.3 编译修复指南

如果你遇到了上述错误，需要重新编译 ABACUS。请严格按照以下步骤操作：

#### 1. 获取完整的源码（解决 LibRI 缺失的根源）
很多用户直接使用 `git clone` 下载 ABACUS 源码，却忘记了下载子模块。`LibRI` 和 `LibComm` 位于源码的 `deps` 目录下。

**检查方法**：
查看 ABACUS 源码目录下的 `deps/LibRI` 和 `deps/LibComm` 文件夹。如果它们是空的，说明你缺少了源码。

**修复命令**（在 ABACUS 源码根目录下执行）：
```bash
git submodule init
git submodule update --remote --recursive
```
*教授提示：这一步至关重要。如果 `deps` 文件夹为空，后续编译无论怎么设置参数都会失败。*

#### 2. 编译时的 CMake 设置
在编译 ABACUS 时，必须显式开启对这些库的支持。请参考 ABACUS 官方文档的 [Advanced Installation Options](https://abacus.deepmodeling.com/en/latest/advanced/install.html) 章节。

通常，你需要确保 CMake 能够找到 Libxc，并启用 LibRI。一个典型的 CMake 配置逻辑如下（仅供参考，具体路径需根据你的集群环境修改）：

```bash
cmake -B build \
      -DENABLE_LIBXC=ON \
      -DENABLE_LIBRI=ON \
      ... (其他参数)
```

**关键点**：
*   确保 `Libxc` 已经安装在系统中，或者通过 `module load` 加载。
*   `ENABLE_LIBRI=ON` 会告诉编译器去编译 `deps` 目录下的 LibRI 源码并链接。

### 1.2.4 运行时的环境准备
如果你使用的是集群上的模块（Module），运行时不仅需要 ABACUS 的可执行文件，有时还需要加载对应的库环境。

```bash
# 示例：运行时可能需要加载 Libxc
module load libxc/5.1.0  # 版本号视具体情况而定
mpirun -np 4 abacus
```

---

## 本章小结

1.  **物理意义**：杂化泛函通过引入 Hartree-Fock 交换项（EXX），有效修正了 DFT 的自相互作用误差，是计算**带隙**和**局域化电子态**的首选。
2.  **核心依赖**：ABACUS (LCAO) 运行杂化泛函必须依赖 **Libxc**（泛函定义）、**LibRI**（积分加速）和 **LibComm**。
3.  **常见坑**：源码 `deps` 目录为空。务必使用 `git submodule update --remote --recursive` 补全子模块。
4.  **启动开关**：在 INPUT 文件中设置 `dft_functional hse` (或 `pbe0`) 即可开启杂化泛函模式，但这需要编译环境的支持。

做好了这些准备，在下一章中，我们将深入具体的参数设置，教你如何优雅地驾驭 HSE06 计算。