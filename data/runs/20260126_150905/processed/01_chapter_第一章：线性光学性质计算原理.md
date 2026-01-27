# 第一章：线性光学性质计算原理

在开始具体的代码操作之前，我们需要建立清晰的物理图像。计算材料的光学性质，本质上是求解光（电磁波）与材料中电子的相互作用。

本章将带你理解这一过程的两个核心环节：
1.  **宏观层面**：如何通过介电函数（Dielectric Function）推导出折射率、吸收系数等实验可测量的光学常数。
2.  **微观层面**：ABACUS 如何利用原子轨道（LCAO）基组构建紧束缚模型，并将电子结构信息传递给后处理工具 `pyatb` 进行光导率计算。

---

## 1.1 介电函数与宏观光学常数

在线性响应理论中，材料对光的响应完全由复介电函数 $\epsilon(\omega)$ 描述：

$$
\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)
$$

### 1.1.1 物理意义
*   **虚部 $\epsilon_2(\omega)$（吸收项）**：直接对应电子在能级间的跃迁。当光子的能量 $\hbar\omega$ 等于导带与价带的能量差时，电子吸收光子能量发生跃迁，产生吸收峰。这是我们计算的**起点**。
*   **实部 $\epsilon_1(\omega)$（色散项）**：描述材料的极化程度和折射行为。它并不直接通过波函数计算，而是通过 **Kramers-Kronig (K-K) 关系** 由虚部 $\epsilon_2(\omega)$ 推导而来。

### 1.1.2 导出其他光学性质
一旦我们计算出了 $\epsilon_1(\omega)$ 和 $\epsilon_2(\omega)$，所有其他的线性光学性质都可以通过简单的代数关系得到。这也是为什么我们在计算流程中只关注介电函数的原因。

以下是常用光学常数的推导公式：

1.  **折射率 (Refractive Index, $n$)**：
    $$ n(\omega) = \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} + \epsilon_1}{2}} $$
2.  **消光系数 (Extinction Coefficient, $\kappa$)**：
    $$ \kappa(\omega) = \sqrt{\frac{\sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1}{2}} $$
3.  **吸收系数 (Absorption Coefficient, $\alpha$)**：
    $$ \alpha(\omega) = \frac{2\omega}{c}\kappa(\omega) = \sqrt{\frac{2\omega^2}{c^2} \left( \sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1 \right)} $$
    > **注意**：这里的单位换算至关重要，我们将在 1.3 节详细讨论。
4.  **反射率 (Reflectivity, $R$)**（正入射情况）：
    $$ R(\omega) = \frac{(n-1)^2 + \kappa^2}{(n+1)^2 + \kappa^2} $$
5.  **能量损失函数 (Energy Loss Function, $L$)**：
    $$ L(\omega) = \text{Im}\left(-\frac{1}{\epsilon(\omega)}\right) = \frac{\epsilon_2}{\epsilon_1^2 + \epsilon_2^2} $$

---

## 1.2 基于 LCAO 的紧束缚模型与工作流

ABACUS 的一大特色是支持数值原子轨道（LCAO）基组。在计算光学性质时，我们采用**“两步走”**策略。

### 1.2.1 核心流程：从 DFT 到 Tight-Binding
计算光学性质需要求和布里渊区中大量的 $k$ 点以保证积分精度。直接使用 DFT 进行全布里渊区计算极其昂贵。
ABACUS 的解决方案是：
1.  **Step 1 (ABACUS)**：在稀疏 $k$ 点网格上进行自洽计算（SCF），得到实空间的哈密顿量矩阵 ($H$)、重叠矩阵 ($S$) 和位置算符矩阵 ($r$)。
2.  **Step 2 (pyatb)**：利用这些矩阵构建紧束缚（Tight-Binding）模型，利用 Wannier 插值技术快速得到任意密集的 $k$ 点上的本征值和波函数，进而计算 Kubo-Greenwood 公式。

### 1.2.2 关键矩阵与参数设置
为了让 ABACUS 输出 Step 2 所需的矩阵，必须在 `INPUT` 文件中显式开启以下参数：

| 参数名 | 推荐值 | 物理含义 |
| :--- | :--- | :--- |
| `out_mat_hs2` | `1` | 输出二阶哈密顿矩阵 ($H$) 和重叠矩阵 ($S$)。生成文件通常为 `data-HR-sparse_SPIN0.csr` 和 `data-SR-sparse_SPIN0.csr`。 |
| `out_mat_r` | `1` | 输出位置算符（偶极）矩阵 ($r$)。这是计算光学跃迁矩阵元 $\langle n\mathbf{k}|\mathbf{r}|m\mathbf{k}\rangle$ 的基础。生成文件为 `data-rR-sparse.csr`。 |
| `basis_type` | `lcao` | 必须使用原子轨道基组。 |

> **风险提示**：默认情况下 `out_mat_*` 参数均为 0。如果你忘记设置它们为 1，ABACUS 计算结束后将不会生成任何矩阵文件，导致后续 pyatb 无法运行。

### 1.2.3 数据的“手动”传递 (Critical)
这是新手最容易出错的环节。`pyatb` 作为一个独立的后处理工具，它**不能**自动读取 ABACUS 的二进制输出或日志文件来获取费米能级和占据态信息。

**你必须手动完成以下数据的传递：**

1.  **查找数据**：
    打开 ABACUS 的标准输出文件（通常是屏幕输出重定向的文件，或 `OUT.*/running_scf.log`），在文件末尾寻找以下关键词：
    *   `occupied bands`: 体系的占据能带数。
    *   `E_Fermi`: 费米能级（通常单位为 eV，但也可能是 Ry，需核对 `running_scf.log` 中的单位说明，通常 ABACUS 输出 eV）。

2.  **填入 pyatb**：
    将上述数值准确填入 `pyatb` 的 `Input` 文件中：
    ```text
    OPTICAL_CONDUCTIVITY
    {
     occ_band  [填入 occupied bands 的数值]
     ...
    }
    
    INPUT_PARAMETERS
    {
     fermi_energy [填入 E_Fermi 的数值]
     fermi_energy_unit eV
     ...
    }
    ```

> **警告**：如果费米能级填写错误，会导致带间跃迁的占据数判断错误（$f_{nk}-f_{mk}$），从而计算出完全错误的介电函数谱。

### 1.2.4 晶格信息的补充
`pyatb` 目前版本不会自动从 `.csr` 稀疏矩阵文件中读取晶格常数和晶格矢量。你必须在 `pyatb` 的 `Input` 文件中显式添加 `LATTICE` 块，内容应与 ABACUS 的 `STRU` 文件保持一致（注意单位通常是 Bohr）。

---

## 1.3 后处理中的单位换算逻辑

在 Step 2 完成后，`pyatb` 会输出介电函数的实部和虚部。但在绘制**吸收系数 (Absorption Coefficient)** 图谱时，直接套用公式会遇到单位问题。

通常计算得到的能量 $E$ 单位是 eV，而实验常用的吸收系数单位是 $\text{cm}^{-1}$。

**标准换算流程**：

公式：$\alpha(E) = \frac{E}{\hbar c} \sqrt{2 \left( \sqrt{\epsilon_1^2 + \epsilon_2^2} - \epsilon_1 \right)}$

1.  **常数准备**：
    *   普朗克常数 $\hbar \approx 4.1357 \times 10^{-15} \text{ eV}\cdot\text{s}$
    *   光速 $c \approx 3.0 \times 10^{10} \text{ cm/s}$ (注意这里用 cm/s 而非 m/s)
2.  **预因子计算**：
    $$ \frac{1}{\hbar c} \approx \frac{1}{4.1357 \times 10^{-15} \times 3.0 \times 10^{10}} \approx \frac{1}{1.24 \times 10^{-4}} \approx 8065.5 \text{ cm}^{-1}/\text{eV} $$
    
    或者使用更精确的推导逻辑（如本教程案例脚本所示）：
    $$ \text{系数} = \frac{1}{c_{\text{SI}}} \times \frac{1}{\hbar_{\text{SI}}} \dots $$
    
    在 Python 脚本中，为了保持物理清晰，建议保留完整的常数项：
    ```python
    # 示例逻辑
    omega_in_inverse_seconds = energy_in_eV / (4.1357e-15) # E = h_bar * omega
    c_in_cm_per_s = 3.0e10
    prefactor = omega_in_inverse_seconds / c_in_cm_per_s
    # 最终 alpha 单位为 cm^-1
    ```

理解了这些原理后，下一章我们将进入实战环节，对 $\text{SiO}_2$ 体系进行完整的计算操作。