# https://mcresearch.github.io/abacus-user-guide/

> Source: https://mcresearch.github.io/abacus-user-guide/

ABACUS 中文文档主页

一、介绍

ABACUS（Atomic-orbtial Based Ab-initio Computation at UStc，中文名原子算筹）是国产开源密度泛函理论软件，相关介绍 ABACUS 的新闻可在 [ABACUS 新闻稿整理](news.html)查看，以下是一些常用地址：

ABACUS 在 DeepModeling 社区中的 GitHub 仓库地址为：

[https://github.com/deepmodeling/abacus-develop](https://github.com/deepmodeling/abacus-develop)

ABACUS 的 Gitee 镜像仓库地址为：

[https://gitee.com/deepmodeling/abacus-develop](https://gitee.com/deepmodeling/abacus-develop)

ABACUS 网站访问：

文档（包括安装方法、输入输出参数介绍、功能介绍、算例介绍、开发者须知等）：

[https://abacus.deepmodeling.com/en/latest/](https://abacus.deepmodeling.com/en/latest/)

本教程系列旨在帮助新手用户入门了解 ABACUS 的使用。秉着开源软件的理念，本文档是由开源社区的老师同学们贡献所成。如果你也想贡献一份文档，我们十分欢迎，请参考[如何贡献 ABACUS 使用教程](contribute.html)。

本教程中标有 [ Logo的部分可以直接在 ][Bohrium Notebook](https://nb.bohrium.dp.tech/) 上打开。

在 Bohrium Notebook 上快速学习，见[快速开始 ABACUS｜ 自洽 能带 态密度 结构优化](https://nb.bohrium.dp.tech/detail/4641406377)；在 Bohrium 平台上运行大任务，见[教程](https://bohrium-doc.dp.tech/docs/software/ABACUS/)。

二、**用户文档**

2.1 ABACUS 编译教程

[官方编译教程](https://abacus.deepmodeling.com/en/latest/quick_start/easy_install.md)（英文官网）[GCC 编译 ABACUS 教程](abacus-gcc.html)[Intel oneAPI 2024.x 编译 ABACUS 教程](abacus-oneapi.html)[Intel oneAPI 编译 ABACUS 教程](abacus-intel.html)[编译 Nvidia GPU 版本的 ABACUS](abacus-gpu.html)[ABACUS LCAO 基组 GPU 版本使用说明](abacus-gpu-lcao.html)[在超算环境编译 ABACUS 的建议](abacus-hpc.html)[ABACUS 在曙光 DCU 集群上的编译与使用](abacus-dcu.html)- ABACUS toolchain 脚本集
[ABACUS 编译教程系列之一：基于 Intel 编译器](https://www.bilibili.com/video/BV1ZN411L75Z/)（B 站视频）[ABACUS 编译教程系列之二：基于 CUDA](https://www.bilibili.com/video/BV1Jb4y1L7KB/)（B 站视频）[ABACUS 编译教程系列之三：docker 的使用](https://www.bilibili.com/video/BV13C4y1R7DL/)（B 站视频）

2.2 建模

- 准备晶胞和原子位置等信息的文件 STRU：如何转换 STRU 的格式
[ABACUS 如何选择晶胞朝向获得最佳并行效率？以碳纳米管为例](abacus-eff1.html)[ABACUS 如何选择晶胞朝向获得最佳并行效率？以二维氮化硼为例](abacus-eff2.html)[ABACUS 如何选择晶胞朝向获得最佳并行效率？以铜表面一氧化碳吸附为例](abacus-eff3.html)- 准备赝势：
[模守恒赝势生成方法简介](abacus-upf.html) - 数值原子轨道基组生成教程：

2.3 Kohn-Sham 密度泛函理论

[ABACUS 的平面波计算与收敛性测试](abacus-pw.html)[ABACUS 平面波基组下的杂化泛函](abacus-exx.html)- 电子自洽迭代
- ABACUS 使用教程 ｜ 结构优化
- ABACUS 磁性材料计算使用教程
- ABACUS 使用 DFT+U 计算教程 | 基础版
[ABACUS+LibRI 做杂化泛函计算教程](abacus-libri.html)[ABACUS 收敛性问题解决手册](abacus-conv.html)[ABACUS 答疑手册 v0.2 版本](abacus-question.html)[ABACUS 对比 CP2K 精度和效率测试 | Si 的状态方程（EOS）](https://www.bohrium.com/notebooks/77351186918)- 有 VASP 使用背景的用户上手 ABACUS 教程：
[ABACUS新人使用的一些注意事项](https://xmywuqhxb0.feishu.cn/docx/KN3KdqbX6o9S6xxtbtCcD5YPnue)

2.4 分子动力学

2.5 AI 辅助功能

- DeePKS 方法
[ABACUS+DPGEN 使用教程](abacus-dpgen.html)- ABACUS+DeepH 建立碳材料的哈密顿量模型

2.6 特色功能

2.7 后处理

[ABACUS+Atomkit 计算态密度和能带](abacus-dos.html)[ABACUS 计算 PDOS](abacus-pdos.html)[ABACUS 输出部分的电荷密度和波函数及可视化教程](abacus-chg.html)[ABACUS 计算电子局域函数 ELF 使用教程](abacus-elf.html)[ABACUS+Bader charge 分析教程](abacus-bader.html)[ABACUS+pymatgen 计算弹性常数](abacus-elastic.html)[ABACUS+Phonopy 计算声子谱](abacus-phonopy.html)- ABACUS+pyatb 能带反折叠计算
[ABACUS+ShengBTE 计算晶格热导率](abacus-shengbte.html)- ABACUS+Phono3py 计算晶格热导率
[ABACUS+Wannier90 使用教程](abacus-wannier.html)[ABACUS+Candela 使用教程](abacus-candela.html)[ABACUS+USPEX 接口教程](abacus-uspex.html)[ABACUS+Hefei NAMD 使用教程](abacus-namd.html)- ABACUS+ASE 做过渡态计算
[ATST-Tools: ASE-ABACUS 过渡态计算工作流套件与算例](https://github.com/QuantumMisaka/ATST-Tools)[(支持 NEB，Dimer，AutoNEB 等过渡态方法)](https://nb.bohrium.dp.tech/detail/39369325971)[ABACUS-ASE做NEB计算](https://dptechnology.feishu.cn/wiki/wikcnzar41sN8ZtGLtm3PLnarSc)（简单算例）


三、教程

3.1 基于 ABACUS 的表面计算教程

3.2 《计算材料学》采用 ABACUS 的计算模拟实例

- ABACUS 计算模拟实例 | 概述
- ABACUS 计算模拟实例 | I. 原子及小分子气体能量计算
- ABACUS 计算模拟实例 | II. C2H5OH 的振动模式与频率计算
- ABACUS 计算模拟实例 | III. 材料平衡晶格常数计算
- ABACUS 计算模拟实例 | IV. 堆垛层错能的计算
- ABACUS 计算模拟实例 | V. Al 的弹性性能指标计算
- ABACUS 计算模拟实例 | VI. 空位形成能与间隙能计算
- 2024 秋计算材料学-上机练习：ABACUS 能带和态密度计算
- ABACUS 计算模拟实例 | VIII. 基于 HSE06 的态密度与能带计算
- ABACUS 计算模拟实例 | IX. 表面能的计算
- ABACUS 计算模拟实例 | XI. Pt 表面简单物种的吸附能计算
- ABACUS 计算模拟实例 | XII. Pt(111)表面羟基解离的过渡态搜索
- ABACUS 计算模拟实例 | XIII. Pt 表面的 ORR 催化路径

3.3 公众号文章推荐

[新服务器安装ABACUS](https://mp.weixin.qq.com/s/Cc1TWrGaiYPeZJwbo28DeA)[ABACUS中坐标变换——调整真空层方向](https://mp.weixin.qq.com/s/HaI17wtxg--AmZNWrxd-UA)[ABACUS+Wannier90+WannierTools计算Bi2Se3的能带和拓扑性质](https://mp.weixin.qq.com/s/HAByRaMFqScnTsE_6kw8tg)[更新：ABACUS+Wannier90+WannierTools计算Bi2Se3的能带和拓扑性质](https://mp.weixin.qq.com/s/nfpUpt9sClRX0_lHYwkIHw)

3.4 视频推荐

四、**开发者文档**

4.1 基础规范

[ABACUS 的 Github 仓库 Issues 处理流程](develop-issue.html)[ABACUS 开源项目 C++ 代码规范](develop-C++.html)[ABACUS 注释规范：Doxygen 入门 (c++)](develop-dox.html)[ABACUS 线上文档输入参数撰写规范](develop-input.html)[ABACUS 代码存放规范](develop-rule.html)[如何在 ABACUS 中新增一个输入参数（v3.7.0 后）](develop-addinp2.html)[如何在 ABACUS 中进行异构计算](develop-hetero.html)[ABACUS formatter-2.0 版本使用说明书](develop-formatter2.html)[ABACUS 中使用格式化工具 clang-format](develop-format.html)[如何在 ABACUS 中新增一个输入参数（截至 v3.5.3）](develop-addinp.html)

4.2 性能工具

4.3 编程进阶

[ABACUS 中的测试（一）：测试的重要性](develop-test1.html)[ABACUS 中的测试（二）：测试工具 gtest](develop-test2.html)[C++ 程序设计的一些想法](develop-design.html)[文件输出功能的实现代码结构设计建议：以 ABCUS CifParser 为例](develop-cifparser.html)[以格点积分程序为例：一些代码开发习惯小贴士](develop-grid.html)[在 ABACUS 中进行差分测试](algorithm-delta.html)[ABACUS 如何处理内存 bug？](develop-memory.html)[Tensor 类文档 1：构造和使用说明](develop-tensor1.html)[Tensor 类文档 2：使用和拓展](develop-tensor2.html)

4.4 模块介绍

4.5 平面波代码介绍

[Introduction to ABACUS: Path to PW calculation - Part 1](develop-path1.html)[Introduction to ABACUS: Path to PW calculation - Part 2](develop-path2.html)[Introduction to ABACUS: Path to PW calculation - Part 3](develop-path3.html)[Introduction to ABACUS: Path to PW calculation - Part 4](develop-path4.html)[Introduction to ABACUS: Path to PW calculation - Part 5](develop-path5.html)[Introduction to ABACUS: Path to PW calculation - Summary 1](develop-sm1.html)[Introduction to ABACUS: Path to PW calculation - Part 6](develop-path6.html)[Introduction to ABACUS: Path to PW calculation - Part 7](develop-path7.html)[Introduction to ABACUS: Path to PW calculation - Part 8](develop-path8.html)[Introduction to ABACUS: Path to PW calculation - Part 9](develop-path9.html)[Introduction to ABACUS: Path to PW calculation - Part 10](develop-path10.html)[Introduction to ABACUS: Path to PW calculation - Part 11](develop-path11.html)[Introduction to ABACUS: Path to PW calculation - Summary Final](develop-sm2.html)
