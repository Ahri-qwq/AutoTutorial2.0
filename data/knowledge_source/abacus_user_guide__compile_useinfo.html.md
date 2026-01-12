# https://mcresearch.github.io/abacus-user-guide/compile-useinfo.html

> Source: https://mcresearch.github.io/abacus-user-guide/compile-useinfo.html

ABACUS 的-i/-I/--info 功能使用说明

**作者：周徐源，邮箱：xy_z@pku.edu.cn**

**审核：陈默涵，邮箱：mohanchen@pku.edu.cn**

**最后更新时间：2025/12/15**

功能背景

先前 ABACUS 仅支持 `-v/-V/--version`

选项，仅输出版本号，无 `-h/--help`

选项，输出信息过少。接收完#6734 PR 后，`3.9.0.20`

版本([https://github.com/deepmodeling/abacus-develop/releases/tag/v3.9.0.20](https://github.com/deepmodeling/abacus-develop/releases/tag/v3.9.0.20))的 abacus 支持输入 `-i/-I/--info`

选项输出**编译平台**、**功能支持**及**依赖库版本**等详细信息。考虑到需要首先通过 `configure`

命令获取版本号再基于版本号编译，采用引入通过 CMake 自动填充模板头文件的技术路线。

编译注意事项

由于在 CMake 的 configure 过程或 Makefile 通过 `generate_build_info.sh`

生成 `build_info.h`

，为避免污染源代码，需要采用“源外构建”，即在单独的 `BUILD_DIR=build`

中构建

```
mkdir build && cd build
cmake ../ -D...
or
make -f ../source/Makefile ...
```


使用说明

编译完成后假设可执行文件名为 `abacus`

，则通过运行 `abacus -i`

或 `abacus -I`

或 `abacus --info`

，输出包括如下部分的信息：

核心版本信息和平台

- ABACUS 版本：v3.x.x.x
- Git commit 的哈希值
- 目标平台：CPU、GPU、DCU、DSP 或 SW
- 构建类型：Custom 或 Debug

构建基本环境

- 构建用户
- 构建主机
- 构建时间

编译器和编译标志

- C++ 编译器类型与路径
- C++ 编译器版本号
- C++ 编译标志
- 链接标志
- CUDA 编译标志

清查和调试

- 地址清查支持情况
- 调试符号支持情况

CMake 构建概况

- 构建选项
- 依赖包查找结果

并行和通信

- MPI 应用：IntelMPI、OpenMPI
- MPI 版本
- CUDA 感知 MPI：若支持则显示版本号
- OpenMP：若支持则显示版本号

核心数学库

加速器和硬件

- NVIDIA CUDA：若支持则显示版本号（或路径）
- AMD POCm：若支持则显示版本号（或路径）
- CUSolverMP：若支持则显示版本号（或路径）

杂化泛函相关

[Cereal 序列化](https://github.com/USCiLab/cereal)：若支持则显示版本号（或路径）[LibRI](https://github.com/abacusmodeling/LibRI)：若支持则显示版本号（或路径）[LibComm](https://github.com/abacusmodeling/libComm)：若支持则显示版本号（或路径）

人工智能和机器学习相关

[LibTorch](https://github.com/pytorch/pytorch)：若支持则显示版本号（或路径）[Libnpy](https://github.com/llohse/libnpy)：若支持则显示版本号（或路径）[DeePMD-kit 力场](https://github.com/deepmodeling/deepmd-kit)：若支持则显示版本号（或路径）[NEP 力场](https://github.com/brucefan1983/NEP_CPU)：若支持则显示版本号（或路径）[TensorFlow](https://github.com/tensorflow/tensorflow)：若支持则显示版本号（或路径）

测试和其它支持

[Google 测试](https://github.com/google/googletest)：若支持则显示版本号（或路径）[Google 基准测试](https://github.com/google/benchmark)：若支持则显示版本号（或路径）[RapidJSON](https://github.com/Tencent/rapidjson)：若支持则显示版本号（或路径）[PEXSI](https://github.com/HPC-AI-Team/pexsi)：若支持则显示版本号（或路径）[cnpy](https://github.com/rogersce/cnpy)：若支持则显示版本号（或路径）
