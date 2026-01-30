# 第一章：环境准备与资源配置

“工欲善其事，必先利其器”。在进行具体的计算输入文件生成前，必须先准备好核心工具 `abacustest` 以及 ABACUS 计算所需的赝势和轨道库。这是所有后续操作的基础，也是确保计算流程标准化与自动化的关键一步。

## 1.1 赝势与轨道库的获取

ABACUS 的 LCAO（线性组合原子轨道）基组计算高度依赖于高质量的数值原子轨道和赝势。为了简化这一准备过程，我们推荐使用 `abacustest` 工具来获取官方推荐的 `APNS-pp-orb-v1` 库。

请在终端中执行以下命令，将推荐的资源库下载到当前文件夹：

```bash
abacustest model inputs --download-pporb
```

命令执行完成后，当前目录下会自动生成两个文件夹，分别对应赝势库和轨道库：

1.  **`apns-pseudopotentials-v1`**：包含推荐使用的赝势文件。
2.  **`apns-orbitals-efficiency-v1`**：包含配套的数值原子轨道文件。

**资源库特性说明**：
这两个文件夹是对原始 `APNS-pp-orb-v1` 版本的重新打包，具有以下关键特征：
*   **元素覆盖范围**：涵盖了从氢（H）到铋（Bi）的83种元素。
*   **轨道精度**：所有轨道均为 **DZP**（Double Zeta plus Polarization，双ζ加极化）水平，在计算精度与效率之间取得了良好的平衡。
*   **特殊元素处理**：对于镧系元素，4f 电子被视为核电子（frozen core），且选用了截断半径为 8 au 的轨道，以适应常规计算需求。

## 1.2 环境变量配置

为了使 `abacustest` 在后续生成输入文件（如 `STRU` 和 `INPUT`）时能够自动识别并调用上述下载的资源，我们需要配置环境变量。这样可以避免在每次执行命令时都手动指定冗长的文件路径。

请根据你实际下载的路径，在 shell 配置文件（如 `.bashrc` 或 `.zshrc`）或当前终端中设置以下两个环境变量：

*   **`ABACUS_PP_PATH`**：指向赝势文件夹的路径。
*   **`ABACUS_ORB_PATH`**：指向轨道文件夹的路径。

**配置示例**：

假设你将文件下载到了 `/your/path/to/` 目录下，请执行：

```bash
export ABACUS_PP_PATH=/your/path/to/apns-pseudopotentials-v1
export ABACUS_ORB_PATH=/your/path/to/apns-orbitals-efficiency-v1
```

**配置生效后的行为**：
当环境变量配置完成后，使用 `abacustest model inputs` 命令处理结构文件时，程序将自动从这两个路径中检索所需的元素数据，无需额外干预。

**备选方案**：
如果你不希望设置全局环境变量，或者需要临时使用不同版本的库，也可以在生成输入文件时通过命令行参数显式指定路径：
*   `--pp`：指定赝势路径。
*   `--orb`：指定轨道路径。

例如：
```bash
abacustest model inputs ... --pp /path/to/pp --orb /path/to/orb
```

但在本教程的后续章节中，我们将默认你已完成了环境变量的配置，以保持命令行的简洁性。