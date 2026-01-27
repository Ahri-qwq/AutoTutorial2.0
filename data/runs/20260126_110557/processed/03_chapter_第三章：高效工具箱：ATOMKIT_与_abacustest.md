# 第三章：高效工具箱：ATOMKIT 与 abacustest

在前两章中，我们深入了解了 ABACUS 的输入参数和物理原理。然而，在实际科研工作中，手动编写 `STRU` 文件不仅效率低下，而且极易出错——尤其是在处理复杂超胞或需要批量测试数十种结构时。

本章将为你打开 ABACUS 的“高效工具箱”。对于不想编写 Python 脚本的初学者，我们有交互式的 **ATOMKIT**；对于需要高通量生成任务的进阶用户，我们有工程化的 **abacustest**。

但在介绍这些“黑盒”工具之前，作为开发者，我必须先澄清几个关于 `STRU` 文件的核心概念。无论你使用什么工具生成结构，理解这些底层逻辑都是排查报错（Debug）的关键。

---

## 3.0 教授的课前预警：STRU 文件的四大核心坑点

在开始使用自动化工具前，请务必理解以下四个手动修改或检查 `STRU` 时最容易混淆的概念。

### 1. ASE 安装避坑指南 (Critical)
很多用户习惯直接使用 `pip install ase` 安装 ASE，**这是错误的**。官方 ASE 尚未完全合并 ABACUS 的所有接口功能。
要使用 Python 脚本处理 ABACUS 结构，必须安装带有 ABACUS 接口的分支版本：

```bash
# ❌ 不要这样做
# pip install ase

# ✅ 正确的安装方式 (推荐在虚拟环境中进行)
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```

### 2. 神秘的 `1.8897...` 与 LATTICE_CONSTANT
打开一个生成的 `STRU` 文件，你常会看到这样的设置：
```text
LATTICE_CONSTANT
1.8897261258369282
```
**这是什么？**
ABACUS 内部计算使用原子单位制（Bohr），而晶体结构通常使用埃（Angstrom）。
*   `1 Bohr ≈ 0.529177 Angstrom`
*   `1 Angstrom ≈ 1.889726 Bohr`

**逻辑如下：**
*   **场景 A (推荐)**：你的 `LATTICE_VECTORS` 是以**埃 (Angstrom)** 为单位写的。此时 `LATTICE_CONSTANT` 必须设为 `1.8897...`，ABACUS 会自动将矢量乘以这个常数，转换为内部的 Bohr。
*   **场景 B**：你的 `LATTICE_VECTORS` 是以**Bohr** 为单位写的。此时 `LATTICE_CONSTANT` 设为 `1.0`。

> **风险提示**：虽然 `ATOMIC_POSITIONS` 支持 `Cartesian_angstrom` 选项（强制使用埃），但为了保证与后处理工具的最大兼容性，建议始终使用标准的 `LATTICE_CONSTANT` (1.8897...) 配合 `Direct`（分数坐标）或 `Cartesian`（以晶格常数为缩放因子的笛卡尔坐标）。

### 3. LCAO 与 PW 的区别
*   **平面波 (PW) 计算**：`STRU` 文件中**不需要** `NUMERICAL_ORBITAL` 这一栏。
*   **原子轨道 (LCAO) 计算**：`STRU` 文件中**必须包含** `NUMERICAL_ORBITAL`，且轨道文件的顺序必须与 `ATOMIC_SPECIES` 中的元素顺序严格对应。

### 4. 磁矩在哪里设置？
VASP 用户习惯在 `INCAR` 中通过 `MAGMOM` 设置磁矩。但在 ABACUS 中，**初始磁矩是结构的一部分**，必须在 `STRU` 文件中逐个原子设置。

```text
ATOMIC_POSITIONS
Direct
Fe
1.0  # 初始磁矩 (mag) - 低优先级，通常被忽略或作为默认值
2    # 原子数量
0.0 0.0 0.0 1 1 1 mag 2.5  # <--- 高优先级：针对该特定原子的磁矩
0.5 0.5 0.5 1 1 1 mag -2.5 # <--- 反铁磁设置示例
```

---

## 3.1 ATOMKIT：交互式快速转换

如果你手头有一个 `CIF` 文件，不想写 Python 代码，只想“一键”转成 `STRU`，**ATOMKIT** 是目前最便捷的交互式工具。它由 VASPKIT 团队开发，操作逻辑非常符合直觉。

### 3.1.1 安装与启动
下载并解压 ATOMKIT 后，直接运行：
```bash
~/software/atomkit/atomkit
```

### 3.1.2 实战演练：CIF 转 STRU
假设你当前目录下有一个 `Si.cif` 文件。启动 ATOMKIT 后，你会看到一个数字菜单。请按照以下“密码”顺序操作：

1.  **输入 `1`**：选择 `1) Export Structures`（导出结构）。
2.  **输入 `101`**：选择目标格式 `101) ABACUS (*.stru)`。
3.  **输入 `113`**：选择源文件格式 `113) CIF Format (*.cif)`。
    *   *注：你也可以直接输入 `113 Si.cif` 跳过文件名询问。*
4.  **输入文件名**：`Si.cif`。

**屏幕输出示例**：
```text
-->> (01) Reading Structure from Si.cif File...
...
-->> (02) Written Si.STRU File!
```

此时目录下会生成 `Si.STRU`。打开检查，你会发现 ATOMKIT 已经自动处理了：
*   晶格常数转换为 `1.889726`。
*   坐标转换为 `Direct` 分数坐标。
*   生成了 `ATOMIC_SPECIES` 和 `NUMERICAL_ORBITAL` 的模板。

> **注意**：ATOMKIT 无法预知你具体使用哪个版本的赝势和轨道文件。它生成的 `STRU` 文件中，赝势文件名通常是占位符（如 `Si.upf`），你需要手动修改为实际的文件名（如 `Si_ONCV_PBE-1.0.upf`）。

---

## 3.2 abacustest：工程化输入准备

当你需要对一个文件夹里的 50 个结构进行批量计算，或者希望自动匹配赝势和轨道文件时，交互式的 ATOMKIT 就显得太慢了。这时我们需要使用 **abacustest**。

它是 ABACUS 官方提供的测试与前处理框架，具备强大的 `model inputs` 功能。

### 3.2.1 核心命令：model inputs
该命令可以扫描指定目录，自动匹配结构、赝势和轨道，生成完整的输入文件集。

**基本语法**：
```bash
abacustest model inputs -f [结构文件列表] --pp [赝势目录] --orb [轨道目录]
```

### 3.2.2 实战演练：批量准备任务
假设你的工作目录结构如下：
```text
/work
├── structures/
│   ├── Si.cif
│   └── GaAs.cif
├── pp_dir/          # 存放 .upf 赝势文件
└── orb_dir/         # 存放 .orb 轨道文件
```

**步骤 1：准备环境**
确保 `pp_dir` 和 `orb_dir` 中的文件名包含元素符号（如 `Si_ONCV.upf`），或者在目录下提供 `element.json` 映射文件，以便程序能自动识别。

**步骤 2：执行命令**
```bash
abacustest model inputs \
    -f structures/*.cif \
    --pp ./pp_dir \
    --orb ./orb_dir \
    --d ./jobs  # 指定输出目录
```

**步骤 3：检查结果**
程序会自动在 `jobs` 目录下生成 `Si` 和 `GaAs` 两个子文件夹，每个文件夹内包含：
*   `STRU`：已自动填好对应的赝势和轨道文件名。
*   `INPUT`：基于默认模板生成的输入参数。
*   `KPT`：自动生成的 K 点文件。

### 3.2.3 进阶技巧：自定义参数模板
如果你不希望使用默认的 `INPUT` 参数，可以先准备一个 `INPUT_template` 文件，然后通过 `--input` 参数传递给 abacustest：

```bash
abacustest model inputs -f *.cif --pp ./pp --orb ./orb --input INPUT_template
```
这样，所有生成的任务都会继承你模板中的参数（如 `ecutwfc`, `mixing_type` 等），同时针对每个结构自动调整 `STRU`。

---

## 3.3 附录：ASE 高级脚本模板 (Power User)

对于习惯使用 Python 的用户，这里提供一份**能够同时处理赝势和轨道映射**的 ASE 转换脚本。这是实现“完美转换”的关键。

```python
from ase.io import read, write
import os

# 1. 读取结构 (支持 cif, poscar, xsd 等)
atoms = read("Si.cif")

# 2. 定义赝势映射 (必须与实际文件名一致)
# 键是元素符号，值是文件名
pp_map = {
    "Si": "Si_ONCV_PBE-1.0.upf",
    "O":  "O_ONCV_PBE-1.0.upf"
}

# 3. 定义轨道映射 (仅 LCAO 需要)
orb_map = {
    "Si": "Si_gga_7au_100Ry_4s4p1d.orb",
    "O":  "O_gga_6au_100Ry_2s2p1d.orb"
}

# 4. 设置磁矩 (可选)
# 注意：ABACUS 接口会读取 atoms 对象的初始磁矩并写入 STRU
# atoms.set_initial_magnetic_moments([1.0, -1.0]) 

# 5. 输出 STRU
# format='abacus' 是关键
write("STRU", atoms, format='abacus', pp=pp_map, basis=orb_map)

print("STRU file generated successfully with PP and ORB mappings.")
```

**总结**：
*   **新手/偶尔转换**：使用 **ATOMKIT** (113 -> 101)。
*   **批量/工程化**：使用 **abacustest** (`model inputs`)。
*   **定制化/脚本控**：使用 **ASE-ABACUS** (注意安装正确的分支)。

掌握了这些工具，你已经跨过了 ABACUS 使用门槛中最高的一级台阶。下一章，我们将正式运行计算，并解读输出文件。