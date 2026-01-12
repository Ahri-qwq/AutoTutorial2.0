## 第一章：几何结构与晶格建模

**本章逻辑**: 了解材料的几何结构是进行任何计算的基础。本章将介绍如何使用ABACUS定义晶体结构，并为后续的电子结构计算打下坚实的基础。

### Section 1.1: 晶体结构基础

在开始具体的计算之前，我们需要先理解一些基本的晶体学概念。晶体是由原子、分子或离子按照一定的周期性排列而成的固体。这种周期性的排列可以通过布拉维格子（Bravais lattice）来描述。布拉维格子是一个无限重复的点阵，每个点代表一个晶胞中的等效位置。常见的布拉维格子包括简单立方（SC）、体心立方（BCC）、面心立方（FCC）、六方密堆积（HCP）等。

- **布拉维格子**：布拉维格子是描述晶体中原子排列的基本单位。它由三个基矢量 \(\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3\) 定义，这些基矢量决定了晶胞的形状和大小。
- **晶胞参数**：晶胞参数是指晶胞的边长 \(a, b, c\) 和夹角 \(\alpha, \beta, \gamma\)。这些参数完全确定了晶胞的几何形状。
- **对称性**：晶体具有高度的对称性，这使得我们可以利用对称操作（如旋转、镜像和平移）来简化计算。

### Section 1.2: 实战演练 - 金刚石结构建模

**Mapped Case ID**: `problem_0`

**Selection Reason**: 金刚石结构是半导体材料中最常见的结构之一，其简单的对称性和明确的能带特性使其成为初学者理解晶体结构的理想案例。

#### 理论铺垫

金刚石结构是一种典型的面心立方（FCC）结构，其中每个原子都被四个最近邻原子包围，形成四面体配位。金刚石结构的一个重要特点是它的间接带隙，这意味着电子从价带到导带的跃迁需要同时吸收光子并发生动量变化。

#### 关键操作

1. **定义晶胞参数**：对于硅（Si），晶胞参数 \(a = 5.43 \, \text{Å}\)。
2. **设置原子位置**：在FCC结构中，原子位于晶胞的八个顶点和六个面的中心。
3. **高对称点路径**：选择高对称点路径 \(\Gamma-X-W-K-\Gamma\) 来计算能带结构。

```markdown
STRU:
element: Si
lattice_constant: 5.43
crystal_type: diamond

INPUT:
mode: pyatb
high_symm_points:
  G: [0, 0, 0]
  X: [0.5, 0, 0.5]
  W: [0.5, 0.25, 0.75]
  K: [0.375, 0.375, 0.75]
energy_min: -10
energy_max: 10
insert_point_nums: 30
calculation: cell-relax
ecutwfc: 100
scf_thr: 1e-07
scf_nmax: 100
smearing_sigma: 0.015
mixing_beta: 0.8
ks_solver: genelpa
precision: double
cal_force: 1
cal_stress: 1
relax_method: cg
relax_nmax: 60
force_thr_ev: 0.01
stress_thr: 0.5
fixed_axes: None
relax_cell: 1
stress_thr_kbar: 0.5
max_steps: 60

KPT:
kpath: [G, X, W, K, G]
kspacing: 0.14
```

> **Note**: 在此案例中，`dft_functional` 未明确指定，因此 ABACUS 默认使用 PBE 泛函。PBE 泛函在带隙预测上通常会低估实际值，特别是在半导体材料中。例如，实验测得的硅带隙约为 1.1 eV，而 DFT-PBE 计算结果通常在 0.6-0.7 eV 之间。

> **Tip**: 选择 \(\Gamma-X-W-K-\Gamma\) 路径的原因在于，这条路径涵盖了布里渊区的主要极值点，能够全面展示硅的能带结构特征。特别是 \(\Gamma\) 点处的价带最大值和 \(X\) 点附近的导带最小值，体现了硅的间接带隙特性。

#### 计算步骤

1. **能带结构计算 (Non-SCF)**：首先进行非自洽场（Non-SCF）的能带结构计算。
2. **建模**：构建初始的金刚石结构。
3. **几何优化 (弛豫)**：通过几何优化使结构达到能量最低状态。

### Section 1.3: 拓展案例 - 二维材料石墨烯

**Mapped Case ID**: `problem_1`

**简要提及**: 石墨烯作为一种典型的二维材料，其独特的六角蜂窝结构和Dirac锥特性值得简要介绍。

#### 理论铺垫

石墨烯是一种由碳原子以 sp² 杂化轨道形成的二维蜂窝状晶格。石墨烯的独特性质来源于其特殊的电子结构：在费米能级附近，石墨烯的能带结构呈现出线性的 Dirac 锥，导致其具有非常高的电子迁移率和半金属特性。

#### 关键操作

1. **定义晶胞参数**：石墨烯的晶胞参数 \(a = 2.46 \, \text{Å}\)。
2. **设置原子位置**：石墨烯的晶胞包含两个碳原子，分别位于晶胞的两个子晶格上。
3. **高对称点路径**：选择高对称点路径 \(\Gamma-M-K-\Gamma\) 来计算能带结构。

```markdown
STRU:
element: C
lattice_constant: 2.46
crystal_type: hcp

INPUT:
relax_cell: 1
force_thr_ev: 0.01
stress_thr_kbar: 0.5
max_steps: 60
mode: pyatb
high_symm_points:
  G: [0, 0, 0]
  M: [0.5, 0, 0]
  K: [0.3333, 0.3333, 0]
energy_min: -5
energy_max: 5
insert_point_nums: 30
calculation: cell-relax
ecutwfc: 100
scf_thr: 1e-07
scf_nmax: 100
smearing_sigma: 0.015
mixing_beta: 0.8
ks_solver: genelpa
precision: double
cal_force: 1
cal_stress: 1
relax_method: cg
relax_nmax: 60
stress_thr: 0.5
fixed_axes: None

KPT:
kpath: [G, M, K, G]
kspacing: 0.14
```

> **Note**: 在此案例中，`dft_functional` 未明确指定，因此 ABACUS 默认使用 PBE 泛函。PBE 泛函在带隙预测上通常会低估实际值，特别是在二维材料中。石墨烯实际上是一个零带隙半金属，但在数值计算中可能会出现微小的带隙（如 0.001 eV）。

> **Tip**: 选择 \(\Gamma-M-K-\Gamma\) 路径的原因在于，这条路径涵盖了布里渊区的主要极值点，能够全面展示石墨烯的能带结构特征。特别是 \(K\) 点处的 Dirac 锥，体现了石墨烯的线性色散关系。

#### 计算步骤

1. **建模**：构建初始的石墨烯结构。
2. **几何优化 (弛豫)**：通过几何优化使结构达到能量最低状态。
3. **能带结构计算 (Non-SCF)**：最后进行非自洽场（Non-SCF）的能带结构计算。

通过上述步骤，读者可以逐步掌握如何使用 ABACUS 进行晶体结构建模和能带结构计算。下一章我们将进一步探讨电子基态计算的相关内容。