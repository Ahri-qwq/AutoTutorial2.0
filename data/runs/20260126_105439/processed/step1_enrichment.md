根据您提供的核心案例（`stru.md`）及检索到的知识库，以下是关于 **结构文件 (STRU) 构建** 的结构化元数据报告。

## 1. 物理本质 (Physics Concepts)
- **核心物理概念**: 晶体结构描述 (Crystal Structure Description)。
- **科学问题**: 定义模拟体系的几何构型，包括晶格常数、晶格矢量（周期性边界条件）、原子种类、原子位置坐标、原子移动自由度以及初始磁矩。这是进行任何第一性原理计算（如 SCF、结构弛豫、MD）的物理基础。

## 2. 关键输入参数 (Key Parameters)
**注意**: 本主题主要涉及 `STRU` 文件的编写，而非 `INPUT` 文件。以下列出 `STRU` 文件中的关键关键词（Keywords）。

- **ATOMIC_SPECIES**
    - **物理意义**: 定义计算涉及的元素种类。
    - **格式**: `元素符号 摩尔质量 赝势文件名`。
    - **示例**: `Al 26.9815385 Al_ONCV_PBE-1.0.upf`。

- **NUMERICAL_ORBITAL**
    - **物理意义**: 指定数值原子轨道文件（仅在使用 LCAO 基组时需要）。
    - **格式**: `轨道文件名`。
    - **示例**: `Al_gga_7au_100Ry_4s4p1d.orb`。

- **LATTICE_CONSTANT**
    - **物理意义**: 晶格缩放系数/长度单位转换。ABACUS 内部默认单位为 Bohr。
    - **推荐值**: 通常设为 `1.8897261258...` (即 1 Angstrom 对应的 Bohr 数值)，以便后续使用 Angstrom 单位。
    - **注意**: 如果设为 1.0，则后续晶格矢量单位需直接为 Bohr。

- **LATTICE_VECTORS**
    - **物理意义**: 定义晶胞的基矢量（3x3 矩阵）。
    - **格式**: 九分量形式，每行一个矢量。

- **ATOMIC_POSITIONS**
    - **物理意义**: 定义原子的坐标、磁矩和移动性。
    - **坐标类型**: 
        - `Direct`: 分数坐标（相对于晶格矢量）。
        - `Cartesian`: 笛卡尔坐标（单位由 LATTICE_CONSTANT 决定）。
        - `Cartesian_angstrom`: 直接使用埃为单位的笛卡尔坐标（**注意**：资料提到此选项在接口软件中支持不广，需谨慎使用）。
    - **原子块格式**:
        - `元素符号`
        - `初始磁矩` (低优先级，通常设为 0.0)
        - `原子数量`
        - `x y z move_x move_y move_z mag moment`
            - `move_x/y/z`: `1` 代表可动，`0` 代表固定（用于结构弛豫）。
            - `mag`: 关键词，后接该原子的初始磁矩数值（高优先级）。

- **【知识缺口处理】**:
    - 资料中未详细说明非共线磁性（Non-collinear）在 `STRU` 中的具体写法（如是否涉及 3 分量磁矩），仅展示了标量磁矩 `mag 0.0`。如涉及复杂磁性设置，需提醒查阅官方文档。

## 3. 体系与接口配置 (System & Interfaces)
- **结构 (STRU) 特殊要求**:
    - 必须包含赝势文件（及轨道文件）的准确路径或文件名，且文件需存在。
    - 原子必须按元素分块（Block）书写，不能像 VASP 那样混排。

- **外部接口 (External Tools)**:
    1.  **ASE (Atomic Simulation Environment)**
        - **特定分支**: 必须使用 `ase-abacus` 分支 (GitLab: 1041176461/ase-abacus)，官方 PyPI 的 ASE 暂未完全支持 ABACUS。
        - **操作**: Python 脚本，使用 `write('STRU', structure, format='abacus')`。
    2.  **ATOMKIT**
        - **类型**: 命令行交互式工具。
        - **操作**: 选择功能 `101` (Export to ABACUS) 或 `113` (Import CIF) 等进行转换。
    3.  **abacustest**
        - **类型**: Python 库/命令行工具，支持高通量。
        - **操作**: `abacustest model inputs -f <file> ...`。支持通过 `element.json` 自动匹配赝势和轨道。

- **接口注意事项**:
    - **ASE**: 需要编写 Python 脚本，适合批量处理和复杂建模。需注意安装来源。
    - **ATOMKIT**: 适合快速将 CIF/POSCAR 转为 STRU，交互式操作简单直观。
    - **abacustest**: 适合自动化流程，特别是自动关联 PP/Orb 路径的功能非常实用。

## 4. 教程编写特殊指令 (Special Instructions for Writer)
- **Critical (关于 ASE 安装)**: 必须在教程中显眼位置强调，**不能**直接用 `pip install ase` 获取 ABACUS 接口支持，必须 `git clone` 指定的 `ase-abacus` 仓库并安装。这是新手最容易踩的坑。
- **Critical (关于单位)**: 务必解释 `LATTICE_CONSTANT` 的作用。新手常困惑为何它是 ~1.8897。解释清楚这是将 Angstrom 转换为 Bohr 的系数，因为 ABACUS 底层计算使用原子单位制。
- **Tip (关于赝势/轨道匹配)**: 在介绍 `abacustest` 时，强调其可以通过目录下的 `element.json` 自动查找赝势和轨道，这比手动修改 STRU 文件中的文件名要高效得多。
- **Format**: 在展示 STRU 文件示例时，请务必注释清楚每一部分的含义，特别是 `ATOMIC_POSITIONS` 下的那一行数字（坐标、移动限制、磁矩）分别代表什么。

## 5. 常见报错与注意事项 (Pitfalls)
- **单位混淆**: 用户在 `LATTICE_CONSTANT` 设为 1.8897 的情况下，在 `Cartesian` 坐标中又直接填入了以埃为单位的数值，导致结构被放大了约 1.89 倍（或反之）。
- **文件缺失**: `STRU` 中指定的 `.upf` 或 `.orb` 文件名与实际目录下的文件名不匹配，导致计算无法启动。
- **格式错误**: `ATOMIC_POSITIONS` 中原子未按元素归类，或者原子数量与列表行数不符。
- **ASE 版本**: 使用了标准版 ASE 导致 `format='abacus'` 报错不支持。