# ABACUS 使用教程｜DFT+U 计算

作者：李文菲

### 1. **介绍**

本教程旨在介绍在ABACUS软件中，LCAO基组下的DFT+U功能。在计算含过渡金属体系时，DFT+U方法通常会被采用，因为它可以避免在使用LDA或GGA泛函时出现的d轨道过度分散问题。

关于DFT+U有很多参考材料，其中ABACUS中DFT+U算法细节和简单应用见 [1]。


本教程中将会展示如何在SCF计算中打开+U功能，以及如何使用 occupation matrix control 的功能。

### 2. **准备**

在示例数据集中，我们为你准备好了 DFT+U 计算的示例文件。我们可以直接访问：（你可以在左侧点击数据集查看相应文件）：


```python
! tree /bohr/
```


```python
# 出于安全考虑，我们没有数据集所在文件夹的写入权限，因此我们将其复制到 `/data/` 目录下:
! cp -nr /bohr/ /data/

# 我们在这里定义一些路径，并切换到工作路径，方便后续调用：
import os

bohr_dataset_url = "/bohr/abacus-magnetic-eu2y/v4/"  # url 可从左侧数据集复制
work_path = os.path.join("/data", bohr_dataset_url[1:])
os.chdir(work_path)
print(f"当前路径为：{os.getcwd()}")
```

这是一个反铁磁的NiO算例，其中两个Ni原子上的磁矩彼此相反。让我们来看一下 STRU 文件：


```python
! cat ./ABACUS_DFT+U/STRU
```

可以看到在STRU文件当中，共有三类“元素”：不同磁矩的Ni被定义为不同的atomic species，虽然使用的赝势和轨道文件相同，但原子磁矩不同。

```
ATOMIC_SPECIES
Ni1 58.693 Ni_ONCV_PBE-1.0.upf
Ni2 58.693 Ni_ONCV_PBE-1.0.upf
O   15.999 O_ONCV_PBE-1.0.upf

...
Ni1
2.0 //磁矩

...
Ni2
-2.0 //磁矩
```

现在查看 INPUT 文件：


```python
! cat ./ABACUS_DFT+U/INPUT
```

与DFT+U有关的控制参数主要有三个，如下：

```
#Parameter DFT+U
dft_plus_u    1
orbital_corr    2 2 -1
hubbard_u    5.0 5.0 0.0
```

这三个参数在ABACUS的[线上文档](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#dft-u-correction)中均有说明，在这里再进行简单概述：

- dft_plus_u是总控制，设为1就会在计算中打开+U功能；如果设为0，那么即使设定了+U有关的其他参数，在实际计算中也会忽略。
- orbital_corr是一个数组，长度为ntype（即atomic species的个数），那么在这个例子中就是3。每一个数字控制的是对应的原子类型会在哪一个l量子数上进行+U操作。图中的设置意味着我们会在两个Ni原子的d轨道上进行+U，而最后一位的-1则意味着不对O原子进行+U操作。
- hubbard_u的长度与orbital_corr相同，是U的数值。

此外，在线上文档输入参数介绍的DFT+U一节中还有三个参数，对应的是+U计算的两个额外功能。

- yukawa_potential和yukawa_lambda：在 [1] 文中描述了一种不预先设置U值，而是在进行SCF计算时通过将电子相互作用近似为Yukawa potential，在程序内自行计算U值的方法。如果想要尝试这种方法，可将yukawa_potential设为1。而使用Yukawa potential时需要确定screening length，既可以通过电子密度进行on-the-fly的计算，也可以通过yukawa_lambda手动设置。
- omc：occupation matrix control，可以在+U计算的过程中固定occupation matrix，即[1]文中式10表示的物理量，共三种模式：

$$n_{I, m m^{\prime}}^\sigma=\frac{1}{N_{\mathbf{k}}} \sum_{n \mathbf{k}} f_{n \mathbf{k}}^\sigma\left\langle\psi_{n \mathbf{k}}^\sigma \sigma\left|\hat{P}_{I, m m^{\prime}}^\sigma\right| \psi_{n \mathbf{k}}^\sigma \sigma\right\rangle$$

- omc = 0：无occupation matrix control，进行标准的+U计算。
  - omc = 1：在SCF第一步读入提供的initial_onsite.dm，后续计算中照常更新occupation matrix。
  - omc = 2：读入initial_onsite.dm，并在后续的计算中始终使用这个occupation matrix，不再更新。
  - initial_onsite.dm文件的格式将在下一节计算流程中进行说明。

### 3. **流程**

1. 运行标准的 DFT+U SCF 计算


```bash
%%bash
# 进入工作文件夹
cd ./ABACUS_DFT+U
# OMP_NUM_THREADS=1 表示使用单线程，如果你的机器配置比较高，可以使用多线程，比如 4 线程，就可以写成 OMP_NUM_THREADS=4
# mpirun -n 后面的数字表示计算所使用的 CPU 核心数，这里使用 2 个核心，你可以根据你的机器配置进行修改。
OMP_NUM_THREADS=1 mpirun -n 8 abacus
```

2. 计算结果与分析

在DFT+U算例的文件夹中，运行abacus程序，进入输出文件夹OUT.NiO。

- 首先查看running_scf.log，
  - 搜索FINAL_ETOT，得到总能量：
  
  ```
    --------------------------------------------
    !FINAL_ETOT_IS -9255.7279034240546025 eV
    --------------------------------------------
  ```

- 由于INPUT文件中设置了out_bandgap = 1，我们可以搜索E_bandgap，找到最后一个occurence，得到能隙：

    ```
    E_bandgap               +0.205369322748               +2.794192983776
    ```

- 搜索'absolute magnetism'，找到最后一个occurence，得到总磁矩：

    ```
          total magnetism (Bohr mag/cell) = 0.00000000
       absolute magnetism (Bohr mag/cell) = 3.35321634
    ```

- 由于INPUT文件设置了out_mul = 1，会进行Mulliken charge analysis，生成Mulliken.txt文件。在其中搜索Magnetism，得到原子磁矩：

    ```
    Total Magnetism on atom  Ni1           1.8268646
    Total Magnetism on atom  Ni2          -1.8268646
    Total Magnetism on atom  O      -3.6718263e-13
    Total Magnetism on atom  O       1.7330755e-13
    ```

可以看出，我们确实得到了一个反铁磁的体系。

- occupation matrix：进行DFT+U计算时，会在running log当中会在每一步输出occupation matrix，可以在文件中搜索以'L(S)DA+U'开头的block。这个block的格式如下：
  - 首先，会对每一类+U的“元素”进行说明，给出+U的l量子数和U值，可以看到这里是与输入文件中的设定一致的：

    ```
    atom_type=0  L=2  chi=0    U=5ev
    atom_type=1  L=2  chi=0    U=5ev
    ```

- 接下来，就是输出每个有+U的原子上的occupation matrix；由于d轨道共有五个，因此就是一系列5*5的矩阵。在本算例中，这一段的结构如下：

```
atoms  0 //原子编号，这里指的是第一个Ni原子
L  2 //+U的l channel
zeta  0 //对应l=2的基组，由于这个基组中，Ni原子上只有一个d轨道基组，因此只有zeta = 0
spin  0 //spin up
//first 5*5 matrix
spin  1 //spin down
//second 5*5 matrix
atoms  1 //第二个Ni原子，以下与第一个原子类似
L  2
zeta  0
spin  0
//3rd 5*5 matrix
spin  1
//4th 5*5 matrix
```

- 由于INPUT文件中设置了out_chg 1，在OUT.NiO中也输出了一个onsite.dm，写入了最后一步的occupation matrix。文件的格式与上一条中occupation matrix的输出一样。

3. 尝试occupation matrix control

接着上一步，将onsite.dm粘贴到work directory当中，改名initial_onsite.dm，并在INPUT文件中设置omc = 2，就可以从头开始固定使用这个收敛后的occupation matrix进行+U计算。（当然，在实际使用omc的时候，用户可以自行设置initial_onsite.dm，但这个文件的格式是与onsite.dm完全相同的。）

这样计算得到的最终结果与上一步的SCF是完全相同的。也可以打开OUT.NiO/running_scf.log确认整个过程中的occupation matrix是不变的。

### 4. **结语**

在过渡金属元素体系的计算中，往往会需要用到+U的功能。ABACUS在LCAO基组下实现了DFT+U功能，其输入参数的设置比较直观，实现上的一些细节也在[1]中进行了说明，可以进行参考。

---

[1] Qu X, Xu P, Jiang H, He L, Ren X. DFT+ U within the framework of linear combination of numerical atomic orbitals. The Journal of Chemical Physics. 2022;156(23):234104. [doi:10.1063/5.0090122](https://doi.org/10.1063/5.0090122)
