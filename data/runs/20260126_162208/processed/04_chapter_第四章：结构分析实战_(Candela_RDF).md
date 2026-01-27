# 第四章：结构分析实战 (Candela RDF)

在完成了分子动力学（MD）模拟的“跑分”环节后，我们面对的是海量的原子轨迹数据。如何从这些看似杂乱无章的原子运动中提取出有序的物理规律？这就是后处理（Post-processing）的核心任务。

本章我们将使用 **Candela** 工具进行结构分析。Candela 是专为大规模原子模拟设计的分析工具集，能够高效处理 ABACUS 的输出轨迹。我们将以最基础也最重要的**径向分布函数（RDF）**为例，打通从“生产”到“分析”的全流程。

---

## 4.1 生产与准备：从 ABACUS 到 Candela

很多初学者容易犯的一个错误是：直接拿跑完的 MD 轨迹去计算性质，而忽略了数据的有效性和采样频率。为了确保分析结果的物理意义，我们必须严格遵守 **“生产-检查-分析”** 的三步走流程。

### 第一阶段：生产运行 (Production Run)

在 ABACUS 的 `INPUT` 文件中，有两个参数直接决定了你是否有“粮”可分析。如果设置不当，后续 Candela 将无数据可读。

**ABACUS `INPUT` 关键设置：**

```bash
# ABACUS INPUT file snippet
calculation     md
md_nstep        10000    # 总步数
md_dt           1.0      # 时间步长 (fs)

# --- 关键参数 ---
md_dumpfreq     10       # 轨迹输出频率 (每10步输出一帧)
# ----------------
```

*   **`md_dumpfreq`**: 这是重中之重。如果设为 0 或过大，你将无法获得足够的轨迹帧用于统计。
    *   **推荐值**: 对于 RDF 计算，通常每 10~100 步输出一帧即可。过密会浪费磁盘空间，过疏会导致采样不足。

### 第二阶段：平衡性检查 (Equilibrium Check)

**在运行 Candela 之前，必须进行平衡性检查！**
MD 模拟的前半段通常是系统从初始结构向平衡态弛豫的过程（Equilibration），这段数据**必须被剔除**，否则会污染统计结果。

1.  **检查能量/温度曲线**: 使用 Python 脚本或 Gnuplot 绘制 `OUT.suffix/running_md.log` 中的温度和能量随时间的变化。
2.  **确定截断点**: 找到系统能量和温度开始在均值附近波动的时刻。假设系统在第 2000 步（即第 200 帧，如果 `md_dumpfreq=10`）达到平衡，那么前 200 帧数据必须丢弃。

### 第三阶段：Candela 输入文件配置

进入分析阶段，我们需要配置 Candela 的 `INPUT` 文件。假设你的 ABACUS 输出目录为 `OUT.ABACUS/MD_dump`。

**Candela `INPUT` 示例 (计算 RDF):**

```bash
calculation     pdf            # 计算径向分布函数 (RDF)
geo_in_type     ABACUS         # 指定读取 ABACUS 格式轨迹
geo_directory   ../OUT.ABACUS/MD_dump  # 轨迹文件所在路径

# --- 轨迹筛选与截断 (核心) ---
geo_1           0              # 起始文件索引
geo_2           1000           # 结束文件索引
geo_interval    1              # 读取间隔 (每隔多少个文件读一次)
geo_ignore      200            # 【高亮】忽略前200帧 (平衡期)

# --- RDF 物理参数 ---
rcut            8.0            # 截断半径 (Angstrom)
dr              0.05           # 直方图的分辨率 (Angstrom)

# --- 体系信息 ---
ntype           1              # 元素种类数量
natom           64             # 总原子数
```

#### 关键参数详解

1.  **`calculation pdf`**:
    *   虽然我们要算的是 RDF (Radial Distribution Function)，但在 Candela 中对应的指令通常是 `pdf` (Pair Distribution Function)。
2.  **`geo_ignore` (风险控制)**:
    *   **作用**: 忽略轨迹开头的不平衡帧数。
    *   **设置原则**: 对应上文“平衡性检查”中确定的步数。如果你的系统需要很长时间平衡，务必增大此值。
3.  **`geo_1` / `geo_2` vs `geo_ignore` (易混淆警告)**:
    *   ABACUS 的 `MD_dump` 目录下通常包含一系列文件（如 `STRU_TRAJ_0`, `STRU_TRAJ_1`...）。
    *   **`geo_1` / `geo_2`**: 指定要读取的**文件编号范围**。例如 `geo_1 0` 和 `geo_2 1000` 表示读取从 0 号到 1000 号文件。
    *   **`geo_ignore`**: 指定在读取到的数据中，**跳过前多少帧**。
    *   **逻辑关系**: Candela 会先根据 `geo_1/2` 加载文件，然后从加载的数据流中切掉前 `geo_ignore` 帧。
4.  **`rcut`**:
    *   RDF 的计算范围。通常取晶胞最短边长的一半，或者足以覆盖前几个配位壳层即可（如 8-10 Å）。

---

> **💡 锦囊妙计：关于时间步长的换算 (MSD 预警)**
>
> 虽然本章主讲 RDF，但在配置 Candela 时，新手极易在后续计算 MSD（均方位移）时搞错时间单位。请牢记以下换算逻辑：
>
> Candela 中的 `msd_dt` (ps) 并不等于 ABACUS 的 `md_dt` (fs)。
>
> **计算公式**：
> $$ \text{msd\_dt (ps)} = \frac{\text{ABACUS md\_dt (fs)} \times \text{ABACUS md\_dumpfreq}}{1000} $$
>
> **示例**：
> *   ABACUS 设置: `md_dt = 2.0` (fs), `md_dumpfreq = 10`
> *   Candela 设置: `msd_dt = 0.02` (ps)
>
> *如果不进行此换算，计算出的扩散系数将会有数量级的错误！*

---

## 4.2 多组分体系的 RDF 挑战

对于单质（如纯铝），RDF 只有一条曲线。但对于多组分体系（如 $SiO_2$ 或 $Li_{10}GeP_{2}S_{12}$），我们需要区分不同元素对之间的关联（如 Si-O 键长，或 Li-Li 扩散通道）。

### 输入文件调整

在多组分体系中，`INPUT` 文件需要准确描述原子类型：

```bash
# 针对 SiO2 体系 (假设总原子数 96, 32个Si, 64个O)
calculation     pdf
geo_in_type     ABACUS
# ... (路径设置同上) ...

ntype           2              # 两种元素
natom           96             # 总原子数

# 注意：Candela 通常会自动识别 STRU 中的原子顺序
# 确保 ABACUS 的 STRU 文件中原子按类型分组排列
```

### 结果分析

运行 Candela 后，输出文件通常会包含多列数据或生成多个文件。对于双组分体系（Type 1 和 Type 2），你通常会得到以下几类 RDF：

1.  **1-1 (e.g., Si-Si)**: 描述同类原子的长程有序性。
2.  **1-2 (e.g., Si-O)**: **最重要**。第一峰的位置对应键长，第一峰的积分（Coordination Number, CN）对应配位数。
3.  **2-2 (e.g., O-O)**: 描述阴离子骨架。

**实战技巧**：
*   **配位数计算**: Candela 输出文件的第三列通常是积分值 $N(r)$。
*   要获得 Si 原子的氧配位数，需查看 **Si-O (1-2)** 的 RDF 曲线。
*   读取 RDF 第一极小值处对应的积分值，即为第一配位层的平均配位数。

### 常见报错排除

*   **Error: `box size is too small`**: 你的 `rcut` 设置超过了模拟盒子的一半。请减小 `rcut`。
*   **Error: `segmentation fault`**: 检查 `natom` 是否与轨迹文件中的实际原子数一致，或者 `geo_directory` 路径是否为空。