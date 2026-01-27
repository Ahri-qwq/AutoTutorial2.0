根据提供的资料，以下是关于“静电势与功函数”计算的结构化元数据报告：

# ABACUS Metadata: 静电势与功函数 (Electrostatic Potential & Work Function)

## 1. 物理本质 (Physics Concepts)
- **核心概念**: 
    - **静电势 (Electrostatic Potential, $V_{static}$)**: 将单位电荷从参考点移动到场内某一点所做的功。在 DFT 框架下，定义为 $V_{static} = V_H + V_{ext} + V_{efield} + V_{dipole}$（包含哈特里势、离子局域势、外场及偶极修正，**不包含交换关联势**）。
    - **功函数 (Work Function, $\Phi$)**: 把电子从固体内部移到真空中所需的最小能量。
- **科学问题**: 
    - 计算表面体系的电子逸出功。
    - 分析催化活性、化学反应机理及电子输运性质。
    - 验证偶极修正（Dipole Correction）是否正确消除表面伪电场。

## 2. 关键输入参数 (Key Parameters)

### INPUT 文件参数
| 参数名 | 推荐值 | 物理意义与说明 |
| :--- | :--- | :--- |
| `out_pot` | `2` | **核心参数**。<br>`2`: 输出静电势（Electrostatic Potential），不含交换关联势，用于计算功函数。<br>`1`: 输出总局域势（Total Local Potential），包含交换关联势。<br>**注意**：计算功函数必须设为 `2`。 |
| `efield_pos_max` | 视体系而定 (0.0~1.0) | **偶极修正相关**。<br>锯齿状电势（Sawtooth potential）的最大值位置（分数坐标）。通常置于真空层中心。 |
| `efield_pos_dec` | 视体系而定 (如 0.1) | **偶极修正相关**。<br>电势从最大值衰减到最小值的长度（分数坐标范围）。锯齿位于 `efield_pos_max` 到 `efield_pos_max + efield_pos_dec` 区域。 |

### 【重要】知识缺口处理
- **偶极修正开关**: 资料中提到了偶极修正的参数设置（`efield_pos_max` 等）和公式项 ($V_{dipole}$)，但未明确给出在 `INPUT` 中**开启**偶极修正功能的具体布尔值参数名（通常可能是 `dipole_correction`，但资料未显示）。**请撰稿人查阅官方文档确认开启该功能的具体参数名**。
- **单位转换**: 资料指出 `ElecStaticPot.cube` 单位为 a.u. (Hartree)，而费米能 $E_F$ 通常以 eV 为单位。计算功函数时的单位换算细节（是否由后处理脚本自动完成）资料未详述，需确认。

## 3. 体系与接口配置 (System & Interfaces)

- **结构 (STRU)**:
    - 必须包含**真空层**（Vacuum layer）的表面模型（Slab model）。
    - 建议原子位于超胞中部，真空层留在边缘，或者根据 `efield_pos_max` 调整原子位置，确保原子不位于电势锯齿突变区域。
- **输出文件**:
    - `OUT.${suffix}/ElecStaticPot.cube`: 3D 实空间静电势文件（单位 a.u.）。
    - `OUT.${suffix}/running_scf.log`: 提取费米能级 ($E_F$)。
    - `ElecStaticPot_AVE` (资料提及): 可能由后处理生成的平面平均静电势文件。
- **外部接口**:
    - **可视化**: [VESTA](https://jp-minerals.org/vesta/en/) (直接读取 `.cube` 文件)。
    - **数据处理**: Python 脚本（ABACUS 提供）用于将 3D cube 文件转换为沿 Z 轴的平面平均电势（Planar Average Potential），以便读取真空能级。

## 4. 教程编写特殊指令 (Special Instructions for Writer)

- **Critical (核心区分)**: 
    - 必须着重强调 `out_pot 1` 和 `out_pot 2` 的区别。计算功函数**错误地使用 `out_pot 1`（包含 XC 势）会导致结果错误**。
- **计算公式明确化**: 
    - 教程中需明确写出计算步骤：$\Phi = V_{vac} - E_F$。
    - $V_{vac}$ (真空静电势): 从平面平均静电势图的真空层“平台”处读取。
    - $E_F$ (费米能): 从 scf log 中读取。
- **锦囊妙计 (Tips)**:
    - 在绘制静电势沿 Z 轴分布图时，指导用户如何判断偶极修正是否合理：真空层区域应出现平坦的电势平台（Plateau）。如果真空层电势倾斜，说明偶极修正未正确设置或真空层太薄。
    - 提醒用户注意 `ElecStaticPot.cube` 的单位是 Hartree (a.u.)，而通常发表文章使用 eV，作图或计算时可能需要乘以 27.2114。

## 5. 常见报错与注意事项 (Pitfalls)

- **锯齿位置错误**: 如果 `efield_pos_max` 设置的位置穿过原子层，会导致严重的非物理结果。必须确保锯齿（电势突变区）位于纯真空区域。
- **真空层厚度不足**: 如果真空层太薄，静电势可能无法在真空中达到常数平台，导致 $V_{vac}$ 无法准确读取。
- **单位混淆**: 混淆 Hartree 和 eV 会导致功函数数值数量级错误。