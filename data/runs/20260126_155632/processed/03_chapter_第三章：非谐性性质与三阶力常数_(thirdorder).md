# 第三章：非谐性性质与三阶力常数 (thirdorder)

在上一章中，我们通过 Phonopy 获得了材料的谐性性质（声子谱）。然而，谐波近似下的声子是具有无限寿命的准粒子，这意味着材料的热导率将是无穷大。为了计算真实的晶格热导率，我们必须引入**非谐性（Anharmonicity）**，即声子之间的相互作用。

三阶力常数（3rd Order Force Constants）是描述三声子散射过程的核心物理量。本章将指导你如何结合 ABACUS 与 `thirdorder` 程序（ShengBTE 套件的一部分）完成这一计算中最昂贵且最易出错的步骤。

## 3.0 全流程概览与工具准备

计算晶格热导率是一个多软件耦合的复杂流程，为了不迷失方向，请务必牢记以下三个阶段：

1.  **2nd Order (谐性)**: ABACUS + Phonopy $\rightarrow$ `FORCE_CONSTANTS_2ND`
    *   *注意*: 需使用案例脚本 `au2si.py` 将单位转换为 ShengBTE 要求的 eV/Å²。
2.  **3rd Order (非谐性)**: ABACUS + thirdorder $\rightarrow$ `FORCE_CONSTANTS_3RD`
    *   *本章重点*: 涉及结构微扰、高精度受力计算及格式转换。
3.  **Transport (输运)**: ShengBTE
    *   *汇总*: 读取上述两个文件及 `CONTROL` 文件计算热导率 $\kappa$。

> **⚠️ 脚本依赖警告**
> 本章中提到的 `pos2stru.py` 和 `aba2vasp.py` 并非 ABACUS 内置命令，而是官方案例库（`examples/interface_ShengBTE`）中提供的 Python 辅助脚本。请务必确保这些脚本已下载并放置在你的工作目录中，且你的环境中已安装 `ASE` (Atomic Simulation Environment) 库。

---

## 3.1 三阶微扰结构生成

`thirdorder` 程序原生支持 VASP 和 QE，但暂未直接支持 ABACUS。因此，我们需要采用“曲线救国”的策略：将 ABACUS 的结构伪装成 POSCAR，生成微扰后再转回 STRU。

### 步骤 1: 准备超胞 POSCAR
首先，将优化好的晶胞（Unit Cell）扩胞为超胞（Supercell）。虽然你可以手动操作，但推荐先将优化后的 `STRU` 转换为 VASP 格式的 `POSCAR`。

### 步骤 2: 生成微扰构型 (Sow)
使用 `thirdorder` 提供的脚本生成一系列微扰结构。假设我们使用 $2\times2\times2$ 的超胞，且考虑第 nearest-neighbor 相互作用（参数 `-2`）：

```bash
# 语法: thirdorder_vasp.py sow [Nx] [Ny] [Nz] [Cutoff]
thirdorder_vasp.py sow 2 2 2 -2
```

运行后，目录下会生成大量的 `3RD.POSCAR.*` 文件（例如 `3RD.POSCAR.01`, `3RD.POSCAR.02` ...）。这些文件包含了原子微小的位移。

### 步骤 3: 格式回转 (POSCAR $\rightarrow$ STRU)
这是最关键的一步。我们需要将这些 POSCAR 转换回 ABACUS 的 STRU 格式。

> **⛔ 严禁使用 dpdata 进行转换**
> 请绝对**不要**使用 `dpdata` 将这些 POSCAR 转为 STRU。
> **原因**: `dpdata` 在处理晶胞时，会强制将晶格矢量旋转为下三角矩阵形式。然而，`thirdorder` 生成的微扰是基于原始晶格方向定义的。如果晶格发生了旋转，原子的位移方向与晶格的相对关系就会被破坏，导致后续计算出的力常数完全错误。

**正确做法**: 使用基于 ASE 的 `pos2stru.py` 脚本，它能保持晶格取向不变。

```bash
# 批量转换脚本示例
for file in 3RD.POSCAR.*; do
    # 提取编号，如 01, 02
    id=${file##*.}
    # 创建对应的计算文件夹
    mkdir SCF-$id
    # 将 POSCAR 转为 STRU 并放入文件夹
    # 假设 pos2stru.py 接受输入文件名和输出文件名作为参数
    python pos2stru.py $file > SCF-$id/STRU
done
```

---

## 3.2 极高精度受力计算 (High-Precision SCF)

现在你需要对这几十甚至上百个微扰结构进行 SCF 计算以获取受力。**这是整个热导率计算中对数值精度要求最高的环节。**

三阶力常数本质上是能量对原子位移的三阶导数（或力对位移的二阶导数）。由于微扰位移量很小，由此产生的受力变化也非常微小。如果 SCF 收敛精度不够，数值噪音将淹没真实的物理信号，导致最终的热导率结果出现巨大误差（甚至出现负值）。

### INPUT 参数设置核心原则

在 `INPUT` 文件中，除了常规的 `calculation = scf` 和 `cal_force = 1` 外，必须严格设置 `scf_thr`。

#### 1. LCAO 基组 (数值原子轨道)
对于 LCAO 基组，推荐精度如下：

```javascript
INPUT_PARAMETERS
{
    calculation     scf
    basis_type      lcao
    cal_force       1
    
    // 核心精度控制
    scf_thr         1e-8    // 必须达到 1e-8 或更小
    force_thr_ev    1e-4    // 辅助判断，但 scf_thr 是关键
}
```

#### 2. PW 基组 (平面波)
平面波基组对噪音更敏感，且通常能达到更高的收敛精度。**强烈建议**设置为 `1e-12`：

```javascript
INPUT_PARAMETERS
{
    calculation     scf
    basis_type      pw
    cal_force       1
    ecutwfc         100     // 确保能量截断足够高
    
    // 核心精度控制
    scf_thr         1e-12   // 粗体警告：1e-8 对于 PW 计算三阶力常数通常是不够的！
}
```

> **专家经验**: 
> 很多初学者发现计算出的热导率与实验值相差甚远，或者声子散射率呈现无规律的杂乱分布，90% 的原因都是 `scf_thr` 设置过大。不要为了节省一点计算时间而牺牲精度，这一步的“脏数据”会导致后续所有工作作废。

---

## 3.3 格式伪装与力常数生成 (Reap)

当所有 `SCF-*` 文件夹中的计算都正常结束后，我们需要收集受力数据。由于 `thirdorder` 只能读取 VASP 的 `vasprun.xml`，我们需要再次进行“伪装”。

### 步骤 1: 生成 vasprun.xml
使用辅助脚本 `aba2vasp.py` 进入每个文件夹，读取 ABACUS 的输出文件（如 `running_scf.log` 或二进制输出），并生成符合 VASP XML 规范的 `vasprun.xml` 文件。

```bash
# 批量处理示例
for dir in SCF-*; do
    cd $dir
    # 运行转换脚本，生成 vasprun.xml
    python ../aba2vasp.py 
    cd ..
done
```
*注：生成的 `vasprun.xml` 不需要包含完整的电子结构信息，只要包含正确的原子坐标和受力 (`<varray name="forces">`) 即可被 `thirdorder` 识别。*

### 步骤 2: 提取三阶力常数
最后，使用 `thirdorder_vasp.py` 的 `reap` 模式收集所有数据。确保命令行参数与 `sow` 阶段完全一致（超胞大小和 cutoff）。

```bash
# 语法: thirdorder_vasp.py reap [Nx] [Ny] [Nz] [Cutoff]
# 利用 find 命令将所有 vasprun.xml 排序后喂给程序
find SCF-* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -2
```

如果一切顺利，当前目录下将生成一个名为 `FORCE_CONSTANTS_3RD` 的文件。

---

## 3.4 最终汇总与风险提示

至此，你已经准备好了 ShengBTE 所需的所有核心文件：
1.  `FORCE_CONSTANTS_2ND` (来自 Phonopy + `au2si.py`)
2.  `FORCE_CONSTANTS_3RD` (来自 thirdorder)
3.  `CONTROL` (ShengBTE 的输入参数文件)

将它们放入同一个文件夹（通常命名为 `shengbte`），运行 ShengBTE 程序即可得到热导率。

### ⚠️ 风险提示：收敛性测试
在官方提供的教学案例中，为了演示速度，通常使用 $2\times2\times2$ 的 K 点网格和较小的超胞。这会导致计算出的硅（Si）热导率约为 100 W/(m·K)，而实验值约为 150 W/(m·K)。

**在实际科研中，你必须进行以下收敛性测试：**
1.  **K 点收敛**: 必须测试更密的 K 点网格（如 $4\times4\times4$ 或更高）。
2.  **超胞尺寸**: 三阶力常数的截断半径对热导率影响很大，需测试更大超胞（如 $3\times3\times3$ 或 $4\times4\times4$）以包含更长程的相互作用。

切记：**Tutorial 只是教你流程，Research 必须追求收敛。**