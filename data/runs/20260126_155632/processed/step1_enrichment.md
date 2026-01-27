根据提供的参考资料，以下是关于 **晶格热导率 (ShengBTE) 计算** 的结构化元数据报告。

---

# Metadata: ABACUS + ShengBTE 晶格热导率计算

## 1. 物理本质 (Physics Concepts)
- **核心物理概念**: 
    - 晶格动力学 (Lattice Dynamics)
    - 声子-声子相互作用 (Phonon-Phonon Interaction)
    - 非谐效应 (Anharmonicity)
    - 玻尔兹曼输运方程 (Boltzmann Transport Equation, BTE)
- **解决的科学问题**: 
    - 计算材料的晶格热导率 ($\kappa_L$)。
    - 获得二阶力常数 (声子谱) 和三阶力常数 (散射率)。
    - 理解材料热输运机制（如声子寿命、平均自由程）。

## 2. 关键输入参数 (Key Parameters)

### ABACUS INPUT 参数
| 参数名 | 推荐值 / 设定要求 | 物理意义 | 备注/来源 |
| :--- | :--- | :--- | :--- |
| `scf_thr` | **LCAO**: `1e-8`<br>**PW**: `1e-12` | 自洽场收敛精度 (电荷密度/能量) | **关键参数**。计算三阶力常数对微扰极其敏感，必须使用极高精度。资料明确指出 PW 需更严格。 |
| `stru_file` | `STRU-001` (示例) | 指定读取的结构文件名 | 方便批量计算 Phonopy 产生的微扰结构，无需重命名为 STRU。 |
| `basis_type` | `lcao` 或 `pw` | 基组类型 | 资料中对比了两种基组的结果，均可用于此计算。 |
| `calculation` | `scf` | 计算类型 | 用于计算微扰后的原子受力。 |
| `ecut` | 需收敛测试 (例: 100 Ry) | 平面波截断能量 | 资料提及 LCAO 中也需设置此参数。 |

### 外部软件关键配置 (非 ABACUS INPUT，但必须提及)
- **Phonopy (`band.conf`)**:
    - `FULL_FORCE_CONSTANTS = .TRUE.`: **必须开启**。ShengBTE 读取二阶力常数时需要完整的力常数矩阵，否则报错。
- **ShengBTE (`CONTROL`)**:
    - `nelements`, `natoms`, `ngrid`, `lattvec`, `scell`: 描述晶胞和超胞信息。

### 【重要】知识缺口处理
- **应力计算**: 资料未提及是否需要开启应力计算 (`cal_stress`)，通常计算力常数主要依赖 Force，但需确认是否涉及晶格松弛步骤的应力设置。
- **磁性设置**: 资料中的 Si 案例为非磁性 (`magnetism 0`)，若涉及磁性材料，需查阅文档确认磁矩初始化的影响。

## 3. 体系与接口配置 (System & Interfaces)

### 结构 (STRU)
- **预处理**: 必须先进行高精度的结构优化 (Relax)。
- **超胞**: 计算力常数需要建立超胞 (资料示例为 2x2x2)，实际科研需测试超胞大小收敛性。

### 外部接口工具链
该计算流程极其依赖外部工具的格式转换，流程如下：
1.  **Phonopy**: 用于生成二阶微扰结构和计算二阶力常数。
2.  **thirdorder (thirdorder_vasp.py)**: 用于生成三阶微扰结构和计算三阶力常数。
3.  **ASE**: 用于结构文件格式转换 (`pos2stru.py` 中调用)。
4.  **ShengBTE**: 最终求解热导率。

### 接口注意事项 (Critical)
1.  **单位转换 (Unit Conversion)**:
    - ABACUS+Phonopy 输出的力常数单位为 `eV/(Ang*au)` (1 au = 0.52918 Ang)。
    - ShengBTE 要求的二阶力常数单位为 `eV/Ang^2`。
    - **操作**: 必须运行脚本 `python au2si.py` 进行转换，否则结果错误。
2.  **格式伪装 (Format Masquerading)**:
    - `thirdorder` 程序仅支持 VASP/QE 格式。
    - **操作**: 需使用 `aba2vasp.py` 将 ABACUS 的输出 (Force) 包装成 `vasprun.xml` 格式供 `thirdorder` 读取。
3.  **结构转换避坑**:
    - **禁止使用 `dpdata`**: 资料明确警告，`dpdata` 转换时会强制将晶格转为下三角矩阵（旋转晶格），导致原子受力方向与微扰方向不匹配，计算结果错误。
    - **推荐**: 使用基于 ASE 的脚本 (`pos2stru.py`) 进行 POSCAR 到 STRU 的转换。

## 4. 教程编写特殊指令 (Special Instructions for Writer)

- **锦囊 1 (流程图解)**: 这是一个多步骤、多软件耦合的流程。建议撰稿人清晰地划分三个阶段：
    1.  **2nd Order**: ABACUS + Phonopy (注意 `FULL_FORCE_CONSTANTS` 和 `au2si.py`)。
    2.  **3rd Order**: ABACUS + thirdorder (注意 `aba2vasp.py` 和 `scf_thr`)。
    3.  **ShengBTE**: 汇总计算。
- **锦囊 2 (精度强调)**: 必须在教程中用**粗体**强调 `scf_thr` 的设置。对于 PW 基组，明确指出 `1e-8` 是不够的，必须 `1e-12`。这是新手最容易导致计算结果全是噪音的地方。
- **锦囊 3 (脚本依赖)**: 教程中提到的 `au2si.py`, `pos2stru.py`, `aba2vasp.py` 并非 ABACUS 内置命令，而是案例文件夹中提供的辅助脚本。撰写时需提醒用户确保这些脚本在工作目录中，或从案例库下载。
- **风险提示**: 提醒用户资料中的 2x2x2 K点和超胞仅为教学演示，实际科研必须进行 K 点和超胞尺寸的收敛性测试，否则热导率数值会偏低 (如案例中 Si 计算值 100 vs 实验值 150)。

## 5. 常见报错与注意事项 (Pitfalls)

- **ShengBTE 报错**: 如果 ShengBTE 运行时崩溃或报错，首先检查 `FORCE_CONSTANTS_2ND` 是否设置了 `FULL_FORCE_CONSTANTS = .TRUE.` 生成，以及是否执行了 `au2si.py` 单位转换。
- **三阶力常数全为零或噪音**: 检查 ABACUS `scf_thr` 是否足够小。如果使用了 `dpdata` 进行结构转换，也会导致力常数计算完全错误。
- **计算耗时**: 三阶力常数计算需要对数十甚至上百个微扰结构进行 SCF 计算，建议使用脚本批量提交作业。
- **文件覆盖**: 在处理大量 SCF 文件夹时，注意不要混淆不同微扰构型的输出文件。