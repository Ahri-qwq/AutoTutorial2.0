根据您提供的核心案例（ABACUS + pyatb 计算介电函数）及相关知识库，以下是为该计算主题准备的结构化元数据。

---

# Metadata: 光学性质 (介电函数/吸收谱) 计算

## 1. 物理本质 (Physics Concepts)
- **核心概念**:
    - **线性响应理论 (Linear Response Theory)**: 基于微扰论描述体系对外场的响应。
    - **Kubo-Greenwood 公式**: 用于计算频率依赖的光导率 ($\sigma(\omega)$)，这是连接电子结构与光学性质的核心桥梁。
    - **介电函数 (Dielectric Function, $\epsilon(\omega)$)**: 复数函数 $\epsilon = \epsilon_1 + i\epsilon_2$。虚部 $\epsilon_2$ 直接关联电子跃迁（吸收），实部 $\epsilon_1$ 通过 Kramers-Kronig 关系由虚部推导得出。
    - **紧束缚模型 (Tight-Binding Model)**: ABACUS 利用 LCAO 基组生成哈密顿量 ($H$)、重叠矩阵 ($S$) 和偶极矩阵 ($r$)，pyatb 基于这些矩阵构建紧束缚模型进行后续计算。
- **解决的科学问题**:
    - 预测材料对光的响应，包括**折射率 ($n$)**、**消光系数 ($\kappa$)**、**吸收系数 ($\alpha$)**、**反射率 ($R$)** 和 **能量损失函数 ($L$)**。
    - 分析能带结构（特别是带隙）与光学吸收峰之间的对应关系。

## 2. 关键输入参数 (Key Parameters)

### A. ABACUS `INPUT` 文件 (第一步：SCF 计算)
| 参数名 | 推荐值/示例值 | 物理意义与说明 |
| :--- | :--- | :--- |
| `calculation` | `scf` | 进行自洽场计算以获得基态电子波函数和电荷密度。 |
| `basis_type` | `lcao` | **必须设置**。pyatb 依赖于原子轨道基组生成的稀疏矩阵。 |
| `out_mat_hs2` | `1` | **核心参数**。输出二中心哈密顿量 (H) 和重叠矩阵 (S) 的稀疏矩阵文件 (`data-HR-sparse_SPINx.csr`, `data-SR-sparse_SPINx.csr`)。 |
| `out_mat_r` | `1` | **核心参数**。输出位置/偶极矩阵文件 (`data-rR-sparse.csr`)，用于计算跃迁矩阵元。 |
| `out_chg` | `1` | 输出电荷密度文件。 |
| `nbands` | *(需根据能带需求设置)* | 虽然案例未显式设置，但为了包含足够的高能级空带以计算高频跃迁，通常需要设置比默认占据态更多的能带数。 |

### B. pyatb `Input` 文件 (第二步：光学性质计算)
| 参数名 | 来源/示例值 | 物理意义与说明 |
| :--- | :--- | :--- |
| `fermi_energy` | *从 ABACUS scf日志读取* | **必须精确匹配**。ABACUS SCF 计算得到的费米能级 (eV)，用于确定电子占据分布。 |
| `occ_band` | *从 ABACUS scf日志读取* | **必须精确匹配**。体系的占据能带数。 |
| `HR_route` | `data-HR-sparse_SPIN0.csr` | 指定 ABACUS 生成的哈密顿矩阵文件路径。 |
| `SR_route` | `data-SR-sparse_SPIN0.csr` | 指定 ABACUS 生成的重叠矩阵文件路径。 |
| `rR_route` | `data-rR-sparse.csr` | 指定 ABACUS 生成的偶极矩阵文件路径。 |
| `omega` | `0 30` | 计算的能量/频率范围 (eV)。 |
| `domega` | `0.01` | 能量/频率的步长 (eV)。 |
| `eta` | `0.1` | 洛伦兹展宽因子 (Smearing)，用于平滑光谱，避免狄拉克 $\delta$ 函数导致的数值尖峰。 |
| `grid` | `20 20 20` | 用于积分光导率的 K 点网格密度，通常需要比 SCF 计算的 K 点更密以获得平滑光谱。 |

## 3. 体系与接口配置 (System & Interfaces)
- **结构 (STRU)**:
    - 标准 ABACUS STRU 格式。
    - **注意**: pyatb 的输入文件中需要手动填入 `LATTICE` 信息（晶格常数和矢量），需确保与 ABACUS `STRU` 文件中的定义完全一致。
- **外部接口 (pyatb)**:
    - 本流程强依赖 **pyatb** (Python Atomic Tight-Binding) 扩展包。
    - **数据流**: ABACUS (产生 `.csr` 矩阵) $\rightarrow$ 文件复制 $\rightarrow$ pyatb (读取矩阵与参数，计算 $\epsilon_2$) $\rightarrow$ Python 脚本 (后处理计算 $\epsilon_1, n, \alpha$ 等)。
- **接口注意事项**:
    - ABACUS 计算完成后，必须将生成的 `data-HR-sparse_*.csr`, `data-SR-sparse_*.csr`, `data-rR-sparse.csr` 移动或复制到 pyatb 的工作目录下。

## 4. 教程编写特殊指令 (Special Instructions for Writer)
- **Critical (流程强调)**: 务必强调这是一个**两步走**的流程。第一步 ABACUS 计算仅仅是为了“生产数据”（矩阵），真正的光学性质计算是在第二步通过 pyatb 完成的。
- **Critical (参数一致性)**: 在编写教程步骤时，必须明确指示读者如何从 ABACUS 的输出文件（如 `running_scf.log` 或 `OUT.*/` 目录）中查找 `E_Fermi` 和 `occupied bands`，并将这两个数值准确填入 pyatb 的 `Input` 文件中。这是最容易出错的地方。
- **Critical (单位换算)**: 在后处理部分（Python 脚本计算吸收系数），必须解释单位换算逻辑。知识库中提供了详细的推导（从 eV/s 到 cm$^{-1}$），教程中应保留这一逻辑，因为直接输出的数据通常是原子单位或 eV 单位，而非实验常用的 cm$^{-1}$。
- **风险提示**: 
    - 提醒读者 pyatb 的 `LATTICE` 部分需要手动填写，不要假设它会自动从 `.csr` 文件读取晶格信息。
    - 提醒读者 `out_mat_*` 参数默认为 0，必须显式设置为 1，否则 pyatb 无数据可用。

## 5. 常见报错与注意事项 (Pitfalls)
- **K 点收敛性**: 光学性质计算（尤其是吸收谱）对 K 点密度非常敏感。pyatb 中的 `grid` 参数通常需要比 ABACUS SCF 中的 `K_POINTS` 大得多，否则光谱会出现非物理的震荡或锯齿。
- **空带数量**: 如果关注高能区的吸收谱，ABACUS SCF 计算时可能需要增加 `nbands` 以包含足够的高能级空轨道。
- **文件缺失**: 常见报错是 pyatb 找不到 `.csr` 文件。需检查 ABACUS 是否成功运行结束，以及是否正确开启了 `out_mat_hs2` 和 `out_mat_r`。
- **费米能级错误**: 如果 pyatb 输入的费米能级与 SCF 不一致，会导致占据态判断错误，进而导致错误的跃迁几率（例如允许了占据态到占据态的跃迁），光谱结果将完全错误。