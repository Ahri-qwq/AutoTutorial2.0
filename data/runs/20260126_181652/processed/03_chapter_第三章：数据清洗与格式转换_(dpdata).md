# 第三章：数据清洗与格式转换 (dpdata)

在上一章中，我们成功运行了 ABACUS 的 AIMD（从头算分子动力学）任务，生成了包含原子轨迹、受力和维里（Virial）信息的原始数据。然而，DeepMD-kit 无法直接读取 ABACUS 的原始输出文件。

我们需要一位“翻译官”——**dpdata**。

本章将带你通过 Python 脚本，将 ABACUS 的 DFT 计算结果清洗、打乱并转化为 DeepMD 训练所需的压缩 NumPy 格式（`.npy`）。

---

## 3.1 数据的“翻译”：使用 dpdata 读取 ABACUS 轨迹

`dpdata` 是 DeepModeling 社区开发的一个强大的数据处理库，它能够处理几十种第一性原理软件的格式。对于 ABACUS 用户来说，最常用的类是 `LabeledSystem`，因为它不仅包含结构信息（System），还包含标签信息（Label：即能量、受力和维里）。

### 3.1.1 核心前置检查（Pre-flight Check）

在开始转换之前，请务必确认你的 ABACUS `INPUT` 文件中包含了以下关键设置。如果缺失这些参数，生成的轨迹将不包含训练所需的标签，数据将无法使用：

```bash
# INPUT 文件检查清单
calculation     md      # 必须是分子动力学或弛豫
cal_force       1       # 关键：必须输出受力 (Force)
cal_stress      1       # 关键：必须输出应力/维里 (Virial)
```

> **专家提示**：
> 这里必须严格区分 **DFT 数据生成** 与 **DPMD 模型推理**。
> *   **生成训练数据时**（当前阶段）：`esolver_type` 必须是默认的 `ks` (Kohn-Sham DFT) 或显式指定为 `ksdft`。我们是在用量子力学“教会”神经网络。
> *   **使用模型推理时**（未来阶段）：`esolver_type` 才会设置为 `dp`。
> **切勿在生成训练数据时将 `esolver_type` 设为 `dp`，那是逻辑闭环错误。**

### 3.1.2 读取数据脚本

假设你的 ABACUS MD 任务运行在文件夹 `00.data/abacus_md` 下，且包含 `INPUT`、`STRU` 以及输出目录（通常是 `OUT.ABACUS` 或 `OUT.${suffix}`）。

创建一个 Python 脚本 `convert_data.py`：

```python
import dpdata
import numpy as np

# 1. 指定 ABACUS MD 的路径
# fmt='abacus/md' 告诉 dpdata 这是一个 ABACUS 的 MD 轨迹
# 路径应指向包含 INPUT 和 OUT.xxx 文件夹的根目录
data_path = './00.data/abacus_md'

try:
    # 加载数据
    # type_map 指定原子类型的名称列表，需与 STRU 文件中的顺序一致
    ls = dpdata.LabeledSystem(data_path, fmt='abacus/md', type_map=['Al'])
    
    print(f"成功加载数据！")
    print(f"原子数: {ls.get_natoms()}")
    print(f"帧数: {ls.get_nframes()}")
    
    # 简单检查数据完整性
    if ls.get_nframes() > 0:
        print(f"第一帧能量: {ls['energies'][0]}")
        # 检查是否包含维里 (Virial)
        if 'virials' in ls.data:
            print("维里 (Virial) 数据已检测到。")
        else:
            print("警告：未检测到维里数据，请检查 INPUT 中是否设置了 cal_stress 1")
            
except Exception as e:
    print(f"读取错误: {e}")
    print("请检查路径是否正确，以及 ABACUS 任务是否已生成输出文件。")
```

---

## 3.2 数据集划分策略 (Train/Validation)

机器学习的一条铁律是：**永远不要用考试题来训练学生**。我们需要将数据划分为“训练集”和“验证集”。

### 3.2.1 为什么要随机打乱？

MD 轨迹在时间上是高度相关的（第 $t$ 帧和第 $t+1$ 帧的结构非常相似）。如果我们直接切取前 80% 做训练，后 20% 做验证，模型可能无法在验证集中表现良好，因为验证集的状态可能在训练集中从未出现过（例如温度漂移或发生了相变）。

**最佳实践**：使用随机索引（Random Indexing）抽取验证集。

### 3.2.2 划分脚本实现

在上述脚本的基础上，添加以下代码：

```python
# ... (接上文代码)

# 2. 数据集划分策略
n_frames = ls.get_nframes()
validation_ratio = 0.2  # 20% 的数据用于验证

# 生成随机索引
indexes = np.arange(n_frames)
np.random.shuffle(indexes)

# 计算切分点
n_val = int(n_frames * validation_ratio)
idx_val = indexes[:n_val]      # 验证集索引
idx_train = indexes[n_val:]    # 训练集索引

print(f"训练集帧数: {len(idx_train)}")
print(f"验证集帧数: {len(idx_val)}")

# 使用 sub_system 方法创建子数据集
ts_train = ls.sub_system(idx_train)
ts_val = ls.sub_system(idx_val)
```

---

## 3.3 导出 DeepMD 格式

DeepMD-kit 要求的输入格式是一系列 `.npy` 文件，通常组织在 `set.000`, `set.001` 等文件夹中。`dpdata` 提供了极其便捷的一键导出功能。

### 3.3.1 导出脚本

```python
# ... (接上文代码)

# 3. 导出为 DeepMD 的 npy 格式
# 训练集输出到 training_data 目录
ts_train.to_deepmd_npy('deepmd_data/training_set')

# 验证集输出到 validation_data 目录
ts_val.to_deepmd_npy('deepmd_data/validation_set')

print("数据转换完成！输出目录: ./deepmd_data")
```

### 3.3.2 产物结构解析

运行脚本后，你会看到如下目录结构。理解这些文件对于排查错误至关重要：

```text
deepmd_data/
├── training_set
│   ├── set.000/             # 压缩的数据块
│   │   ├── box.npy          # 晶胞参数 (N_frames, 9)
│   │   ├── coord.npy        # 原子坐标 (N_frames, N_atoms * 3)
│   │   ├── energy.npy       # 体系总能量 (N_frames, 1)
│   │   ├── force.npy        # 原子受力 (N_frames, N_atoms * 3)
│   │   └── virial.npy       # 维里张量 (N_frames, 9)
│   ├── type.raw             # 原子类型索引 (如 0 0 0 1 1 ...)
│   └── type_map.raw         # 原子类型名称 (如 Al Cu)
└── validation_set
    ├── set.000/
    ├── type.raw
    └── type_map.raw
```

> **注意**：
> *   `set.000` 是 DeepMD 的默认分包方式。如果数据量极大，`dpdata` 可能会自动生成 `set.001`, `set.002` 等。
> *   `type.raw` 文件将原子映射为整数（0, 1, 2...），这与 `INPUT` 文件中的 `ntype` 顺序对应。

---

## 常见问题与解决方案 (Troubleshooting)

**Q1: 报错 `KeyError: 'virial'`**
*   **原因**: ABACUS 计算时未开启应力计算。
*   **解决**: 检查 DFT 的 `INPUT` 文件，确保设置了 `cal_stress 1`。如果是弛豫任务（`calculation cell-relax`），通常默认会计算，但 MD 任务必须显式指定。

**Q2: `dpdata` 读取不到数据，显示帧数为 0**
*   **原因**: 路径错误或 ABACUS 尚未输出。
*   **解决**:
    1. 确保 `data_path` 指向的是包含 `INPUT` 文件的**目录**，而不是某个具体文件。
    2. 检查该目录下是否存在 `OUT.${suffix}` 文件夹（例如 `OUT.ABACUS`）。
    3. 确保 `OUT.${suffix}` 内部包含 `Running_MD` (MD任务) 或 `running_cell-relax.log` (弛豫任务) 等日志文件。

**Q3: 混合了多种不同原子数的体系怎么办？**
*   **解决**: DeepMD 要求同一个 system 内原子数量和类型必须一致。如果你有不同原子数的体系（例如一个 64 原子，一个 128 原子），需要将它们分别保存为不同的 `deepmd_data` 子目录（例如 `system_64/` 和 `system_128/`），并在训练时的 json 配置文件中同时指向这两个目录。

完成本章后，你已经拥有了标准的“燃料”。下一章，我们将正式启动引擎，配置并训练你的第一个 Deep Potential 模型。