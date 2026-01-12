# 软件案例｜ABACUS

本文介绍如何使用 LBG 在 Bohrium上 运行 ABACUS 任务。

## 简介

ABACUS 密度泛函理论软件是由何力新教授、任新国研究员和陈默涵研究员主导开发的，拥有完全自主知识产权的一款国产开源的密度泛函理论软件。ABACUS 可采用平面波基矢量和数值原子轨道基矢量来进行模拟计算，通过求解 Kohn-Sham 方程得到材料基态电荷密度分布，并由此计算目标材料的各项物理性质。目前，ABACUS 正在发展基于机器学习辅助的泛函模型 DeePKS，为实现跨尺度的分子动力学模拟提供了强有力的基石。此外，软件研发团队还在发展适用于多场景的密度泛函理论（如适用于大尺度计算的无轨道密度泛函理论和适用于高温高压条件的随机波函数密度泛函理论）。

[代码仓库](https://github.com/deepmodeling/abacus-develop)

[软件文档](https://abacus.deepmodeling.com/)

## 如何在 Bohrium 上运行 ABACUS 任务

> 本案例任务运行约需 40s

### 步骤一，准备输入数据

ABACUS 的输入文件（如输入参数、赝势文件、配置文件等）均已存储至 `Bohrium_ABACUS_example` 文件夹内，（你可以在左侧点击数据集查看相应文件）：


```python
# 出于安全考虑，我们没有数据集所在文件夹的写入权限，因此我们将其复制到 `/data/` 目录下:
! cp -nr /bohr/ /data/

# 我们在这里定义一些路径，并切换到工作路径，方便后续调用：
import os

bohr_dataset_url = "/bohr/bohrium-abacus-u26n/v1/"  # url 可从左侧数据集复制
work_path = os.path.join("/data", bohr_dataset_url[1:])
os.chdir(work_path)
print(f"当前路径为：{os.getcwd()}")
```

    当前路径为：/data/bohr/bohrium-abacus-u26n/v1
    

### 步骤二，准备配置文件


```python
# 进入 Bohrium_ABACUS_example 文件夹
try:
    os.chdir("./Bohrium_ABACUS_example")
    print(f"当前路径为：{os.getcwd()}")
except:
    print(f"当前路径{os.getcwd()}，请检查工作路径是否正确或数据集是否已复制。")
```

    当前路径为：/data/bohr/bohrium-abacus-u26n/v1/Bohrium_ABACUS_example
    

**文件夹内已经包含配置文件 `job.json`，您也可以通过以下 Python 代码创建或编辑 `job.json`：**


```python
import json

# 设置 job.json 的参数
job_params = {
    "job_name": "ABACUS test",                                  # Job 名称
    "command": "OMP_NUM_THREADS=1 mpirun -np 8 abacus > log",   # 计算命令
    "log_file": "log_file",                                     # 日志文件
    "backward_files": [],                                       # 后向文件
    "project_id": 12345,                                        # 项目 ID
    "platform": "ali",                                          # 计算平台（ali/paratera）
    "job_type": "container",                                    # Job 类型
    "machine_type": "c16_m64_cpu",                              # 机型
    "image_address": "registry.dp.tech/dptech/abacus:3.1.0",    # 镜像地址
}

with open("job.json", "w") as f:
    json.dump(job_params, f, indent=4, ensure_ascii=False)                     

print(f"已生成 {os.path.join(os.getcwd(), 'job.json')}, \n参数为:\n{json.dumps(job_params, indent=4, ensure_ascii=False)}")
```

    已生成 /data/bohr/bohrium-abacus-u26n/v1/Bohrium_ABACUS_example/job.json, 
    参数为:
    {
        "job_name": "ABACUS test",
        "command": "OMP_NUM_THREADS=1 mpirun -np 8 abacus > log",
        "log_file": "log_file",
        "backward_files": [],
        "project_id": 12345,
        "platform": "ali",
        "job_type": "container",
        "machine_type": "c16_m64_cpu",
        "image_address": "registry.dp.tech/dptech/abacus:3.1.0"
    }
    

**注意**： `"project_id"`需要替换为您自己的项目ID，可在“[项目管理](https://bohrium.dp.tech/projects)”页查看。建议将 MPI 进程数设置为 CPU 核心数 / 2，如本例选用的机型为 16 核 32G 内存的 CPU 机器，则使用 `mpirun -np 8` 来提交作业。

### 步骤三，提交任务

**第一次使用 Lebesgue Utility 需要配置，若您还未配置，请参阅 [Bohrium 帮助文档｜Lebesgue Utility](https://nb.bohrium.dp.tech/detail/6643676953)）**

使用 Lebesgue Utility 提交任务（Lebesgue Utility 是 Bohrium 平台的作业管理系统）：


```python
os.system("lbg job submit -i job.json -p ./")
```

其中：

- -i 指定任务的配置文件，本案例中是 job.json
- -p 指定输入文件所在的目录，Bohrium 会将指定的目录打包上传，在计算节点上解压后，将工作目录切换为该目录。本案例中是 ./

在命令行看到如下输出即表示提交成功。同时可以看到任务的 JOB ID，后续可用此 ID 追踪任务进度。

```Bash
Submit job succeed. JOB GROUP ID: <JOB GROUP ID>, JOB ID: <JOB ID>
```

## 查看任务

您可以在[监控任务文档](https://bohrium-doc.dp.tech/docs/quickstart/Status)中了解如何在 Bohrium 平台查看任务状态。

## 下载结果

您可以在[结果下载文档](https://bohrium-doc.dp.tech/docs/quickstart/Result)中了解如何在 Bohrium 平台下载任务结果。
