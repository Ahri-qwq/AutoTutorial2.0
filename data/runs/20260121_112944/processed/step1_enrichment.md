基于提供的核心案例和知识库，以下是关于 **ABACUS 杂化泛函使用教程** 的结构化元数据：

# 1. 物理本质 (Physics Concepts)
- **核心物理概念**: 杂化泛函 (Hybrid Functional) 通过在 DFT 的交换关联泛函中引入一定比例的 Hartree-Fock (HF) 精确交换能 (Exact Exchange, EXX)，以修正局域密度近似 (LDA) 或广义梯度近似 (GGA) 中的自相互作用误差 (Self-Interaction Error)。
- **解决的科学问题**: 
    - 改善半导体和绝缘体的带隙预测（通常 LDA/GGA 会低估带隙）。
    - 更准确地描述强关联体系或局域化电子态（如 d 电子和 f 电子）。
    - 提高反应势垒、电荷转移激发等性质的计算精度。

# 2. 关键输入参数 (Key Parameters)

## 必须设置 (Mandatory)
- **`dft_functional`**
    - **推荐值**: `hse` (即 HSE06), `pbe0`, 或 `hf`。
    - **物理意义**: 指定使用的杂化泛函类型。设置后会自动配置对应的 EXX 系数。
- **`basis_type`**
    - **推荐值**: `lcao`
    - **物理意义**: 杂化泛函计算目前在数值原子轨道 (LCAO) 基组下效率最高且功能最完善。平面波 (PW) 基组下仅为实验性功能且效率较低。

## 性能优化与算法控制 (Performance & Algorithm - Core Case)
- **`exx_separate_loop`**
    - **推荐值**: `1` (双层循环，省内存) 或 `0` (单层循环，收敛性好)。
    - **物理意义**: 控制自洽迭代的循环方式。
        - `0`: 单层循环，同时更新 GGA 和 EXX。收敛性好，但内存消耗高（需存 `mixing_ndim` 步密度矩阵）。
        - `1`: 双层循环，外层更新 EXX，内层只更新 GGA。内存消耗低，但 GGA 计算次数多。
- **`exx_real_number`**
    - **推荐值**: `1` (如果确认体系具有中心反演对称性或无自旋轨道耦合等导致复数的情况)。
    - **物理意义**: 强制将 EXX 计算限制在实数域。若体系密度矩阵实部为主，开启此项可显著提速（案例中提速约 2.5 倍）。
- **`exx_ccp_rmesh_times`**
    - **推荐值**: `1.5` (针对 HSE), `5` (针对 PBE0/HF)。
    - **物理意义**: 决定计算库仑势所需的径向格点相对于原子轨道格点的倍数。减小此值可显著加速，但需测试精度。

## 泛函系数微调 (Functional Tuning - Optional)
*注：通常由 `dft_functional` 自动设置，仅在需自定义泛函时修改。*
- **`exx_fock_alpha` / `exx_hybrid_alpha`**
    - **物理意义**: Hartree-Fock 精确交换项的混合比例（如 PBE0 默认为 0.25）。
- **`exx_erfc_omega` / `exx_hse_omega`**
    - **物理意义**: 范围分离参数（Range-separation parameter），用于 HSE 等泛函中屏蔽长程交换作用（HSE06 默认为 0.11）。

## 知识缺口处理 (Knowledge Gaps)
- **`ks_solver`**: 资料中提到推荐使用 `genelpa` 或 `scalapack_gvx`，但未明确不同体系大小下的严格界限，建议撰写时提示用户根据编译情况选择。
- **平面波杂化泛函参数**: 资料明确指出平面波杂化泛函为实验性功能且不支持 K 点并行，具体高级参数文档较少，建议教程重点聚焦于 LCAO 基组。

# 3. 体系与接口配置 (System & Interfaces)

- **编译依赖 (Compilation)**:
    - **必须包含**: `Libxc` (泛函库), `LibRI` (电子积分库), `LibComm`。
    - **注意事项**: 若未正确链接 LibRI，运行时会报错 "compile with libri to use hybrid functional in lcao basis" 或 "Unrecognized exchange-correlation functional"。
- **结构文件 (STRU)**:
    - **原子轨道基组**: 建议使用 ABACUS 官网提供的默认截断半径 (Rcut) 和基组大小 (如 dzp)。
    - **警告**: EXX 计算量与 $R_{cut}^9$ 和基组数目四次方成正比，切勿盲目使用大截断半径或高角动量基组 (如 tzdp)，否则计算时间会爆炸式增长。
- **K点设置 (KPT)**:
    - 由于库伦势的长程相互作用，杂化泛函自洽计算所需的 K 点密度可能需要高于普通 GGA 计算。

# 4. 教程编写特殊指令 (Special Instructions for Writer)

- **Critical (并行策略)**: 
    - 必须花费大量篇幅讲解**混合并行策略**。
    - **核心原则**: 杂化泛函极度消耗内存。
    - **最佳实践**: 推荐“少 MPI 进程，多 OpenMP 线程”。
    - **具体示例**: 必须引用案例中的配置，例如单节点 56 核，建议 `mpirun -np 1 -env OMP_NUM_THREADS=56` (最省内存) 或 `mpirun -np 4 -env OMP_NUM_THREADS=14` (速度可能更优，视内存而定)。
    - **对比**: 严禁使用纯 MPI 并行（如单核单进程），这会导致内存溢出。

- **Critical (性能权衡)**:
    - 务必解释 `exx_separate_loop` 的两种模式。教程应提供决策树：
        - 内存不足？ -> 用 `exx_separate_loop 1`。
        - 体系难以收敛？ -> 用 `exx_separate_loop 0`。
    - 强调 `exx_real_number` 的提速效果，但要提示适用条件（实空间密度矩阵虚部为零）。

- **Risk Warning (基组陷阱)**:
    - 明确警告用户不要随意增加基组精度（如从 dzp 换到 tzdp）或截断半径（Rcut），除非他们清楚自己在做什么。案例数据显示 tzdp 比 dzp 慢 6 倍。

- **Scope Limit**:
    - 明确教程主要针对 LCAO 基组。对于 PW 基组，仅作为“实验性功能”简要提及，指出其效率低且不支持 K 点并行。

# 5. 常见报错与注意事项 (Pitfalls)

- **报错**: "Unrecognized exchange-correlation functional 'HSE'"
    - **原因**: 编译时未链接 Libxc 或 LibRI。
- **问题**: 内存溢出 (OOM)
    - **原因**: MPI 进程数开太多，或者使用了单层循环 (`exx_separate_loop 0`) 处理大体系。
    - **对策**: 减少 MPI 进程增加 OpenMP 线程；切换到双层循环 (`exx_separate_loop 1`)。
- **性能**: 计算极慢
    - **原因**: 使用了过大的原子轨道截断半径，或者使用了平面波基组。
- **收敛**: 震荡不收敛
    - **原因**: 双层循环 (`exx_separate_loop 1`) 虽然省内存，但有时收敛性不如单层循环。尝试切换回 `0` 或调整混合参数。