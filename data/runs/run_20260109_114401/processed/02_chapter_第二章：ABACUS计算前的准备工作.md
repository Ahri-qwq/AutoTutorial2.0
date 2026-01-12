# 第二章：ABACUS计算前的准备工作

在开始实际计算之前，了解如何正确设置初始条件是非常重要的一步。本章将详细介绍晶胞优化、INPUT文件详解以及初始结构文件准备等关键步骤。

## Section 2.1: 晶胞优化

### 内容
在进行任何类型的模拟之前，必须先完成晶胞优化过程。这一步骤可以帮助消除内部应力，确保计算结果的准确性和可靠性。晶胞优化通过调整原子位置和晶格参数，使得体系达到能量最低状态。

### 关键参数
- `cal_stress`: 计算应力张量。开启此选项可以得到晶胞应力信息，有助于判断是否需要进一步优化。
- `cal_force`: 计算原子受力。开启此选项可以得到每个原子所受的力，用于优化过程中调整原子位置。

```markdown
示例 INPUT 文件中的相关设置：
```
```bash
cal_stress .true.
cal_force .true.
```

### 物理意义
- `cal_stress`：计算应力张量，用于评估晶胞内部应力。如果应力较大，说明晶胞结构可能不稳定，需要进一步优化。
- `cal_force`：计算原子受力，用于调整原子位置以达到能量最低状态。

## Section 2.2: INPUT文件详解

### 内容
INPUT文件是ABACUS计算的核心配置文件，包含了所有必要的参数设置。本节将逐条解析构成INPUT文件的关键项，特别是那些直接影响到弹性性质计算结果的部分。

### 关键参数
- `cal_stress`: 计算应力张量。
- `cal_force`: 计算原子受力。
- `gamma_only`: 是否只计算Γ点。
- `smearing_method`: 能带展宽方法。
- `smearing_sigma`: 能带展宽参数。
- `mixing_type`: 密度混合类型。

```markdown
示例 INPUT 文件：
```
```bash
INPUT_PARAMETERS
suffix       Mn           # 输出后缀
ntype        1            # 元素种类
ecutwfc      20           # 展开截止能量 (Ry)
scf_thr      1e-7         # 电荷密度收敛阈值
basis_type   pw           # 基函数类型 (pw/lcao)
calculation  scf          # 计算类型
cal_stress   .true.       # 计算应力张量
cal_force    .true.       # 计算原子受力
gamma_only   .false.      # 是否只计算Γ点
smearing_method fermi_dirac # 能带展宽方法
smearing_sigma 0.05       # 能带展宽参数 (Ry)
mixing_type  pulay        # 密度混合类型
```

### 物理意义
- `cal_stress` 和 `cal_force`：如前所述，用于计算应力和受力，确保晶胞优化的准确性。
- `gamma_only`：如果设置为 `.true.`，则只计算Γ点，适用于对称性较高的体系，可以减少计算时间。
- `smearing_method` 和 `smearing_sigma`：能带展宽方法和参数，用于处理部分占据态，提高自洽迭代的稳定性。
- `mixing_type`：密度混合类型，选择合适的混合方法可以加速自洽迭代的收敛。

## Section 2.3: 初始结构文件准备

### 内容
初始结构文件（STRU）描述了目标体系的原子排列方式。本节将教授如何创建或修改STRU格式的文件。

### 关键参数
- `ATOMIC_SPECIES`：定义元素名称、相对原子质量和赝势文件名。
- `LATTICE_CONSTANT`：晶格常数。
- `LATTICE_VECTORS`：晶格向量。
- `ATOMIC_POSITIONS`：原子位置。

```markdown
示例 STRU 文件：
```
```bash
ATOMIC_SPECIES
Mn 54.938045 Mn.upf

LATTICE_CONSTANT
1.8897259886  # 1.8897259886 Bohr = 1.0 Angstrom

LATTICE_VECTORS
3.52000000000000000 0.0000000000000000 0.0000000000000000
-1.7600000000000000 3.03108891324553771 0.0000000000000000
0.0000000000000000 0.0000000000000000 5.2400000000000000

ATOMIC_POSITIONS
Mn 0.0000000000000000 0.0000000000000000 0.0000000000000000
Mn 0.5000000000000000 0.5000000000000000 0.5000000000000000
```

### 物理意义
- `ATOMIC_SPECIES`：定义元素名称、相对原子质量和赝势文件名，确保计算中使用的赝势与元素匹配。
- `LATTICE_CONSTANT`：晶格常数，单位默认为Bohr。
- `LATTICE_VECTORS`：晶格向量，定义晶胞的形状和大小。
- `ATOMIC_POSITIONS`：原子位置，描述原子在晶胞中的具体坐标。

## 配置文件 `config.json`

在进行高通量计算时，通常需要准备一个 `config.json` 文件来管理多个计算任务。该文件包含了一系列参数，用于控制计算流程和输出结果。

```markdown
示例 `config.json` 文件：
```
```json
{
  "system": {
    "ntype": 1,
    "nbands": 10,
    "pseudo_dir": "./",
    "stru_file": "STRU"
  },
  "model": {
    "basis_type": "pw",
    "ecutwfc": 20
  },
  "calculations": [
    {
      "type": "scf",
      "scf_thr": 1e-7
    },
    {
      "type": "relax",
      "relax_nmax": 50
    }
  ]
}
```

### 物理意义
- `system`：系统参数，包括元素种类、能带数、赝势文件路径和结构文件名。
- `model`：模型参数，包括基组类型和平面波截断能量。
- `calculations`：计算任务列表，定义了自洽迭代和结构优化的具体参数。

## 使用 `gene_dfm.py` 和 `compute_dfm.py` 脚本

为了生成和处理形变结构，可以使用 `gene_dfm.py` 和 `compute_dfm.py` 脚本。这些脚本可以帮助自动化生成不同形变下的结构文件，并计算相应的物理性质。

```markdown
示例命令：
```
```bash
python gene_dfm.py -i INPUT -s STRU -o deformed
python compute_dfm.py -i INPUT -d deformed
```

### 物理意义
- `gene_dfm.py`：生成形变结构文件，用于计算弹性常数。
- `compute_dfm.py`：计算形变结构的物理性质，如应力和能量。

### 风险提示
- 关于 `gene_dfm.py` 和 `compute_dfm.py` 的具体使用细节，请查阅 GitHub 仓库中的文档。
- 关于 `INPUT` 文件中某些参数的具体设置，请查阅 ABACUS 官方文档。

## 应力-应变法的优势

在高通量计算中，应力-应变法相对于能量法具有明显优势。应力-应变法可以直接从应力张量中提取弹性常数，避免了能量法中需要大量计算不同形变下能量的复杂性。此外，应力-应变法在处理复杂体系时更为高效和准确。

通过以上步骤，你可以准备好进行ABACUS计算所需的所有初始条件，确保计算结果的准确性和可靠性。