# 用ABACUS-ASE自动产生能带路径

<a href="" target="_blank"><img src="https://cdn.dp.tech/bohrium/web/static/images/open-in-bohrium.svg" alt="Open In Bohrium"/></a>

<font color="grey">


推荐镜像：`abacus-user-guide:3.3.2` \
推荐计算资源：`CPU` \
内容：本教程主要介绍如何用ABACUS-ASE自动产生能带路径。 \
使用方式：您可在 Bohrium Notebook上直接运行。您可以点击界面上方蓝色按钮 `开始连接`，选择 `abacus-user-guide:3.3.2` 镜像及`c4_m8_cpu`款节点配置，稍等片刻即可运行。如您遇到任何问题，请联系 [bohrium@dp.tech](mailto:bohrium@dp.tech) 。 \
共享协议：本作品采用[知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议](https://creativecommons.org/licenses/by-nc-sa/4.0/)进行许可。

</font>

本notebook改编自abacus使用指南，更多信息详见[这里](https://github.com/MCresearch/abacus-user-guide)。

**请注意：运行本notebook需要选择** `Kernel`**为Bash**

## 1. 安装ASE-ABACUS接口：

**本镜像中已经提前安装好ASE-ABACUS接口，可以直接使用**
自行安装可参考指令
* `git clone https://gitlab.com/1041176461/ase-abacus.git`
* `cd ase-abacus`
* `python3 setup.py install`

## 2. python脚本


```bash
cd ~
cp -r /bohr/abacus-band-path-pzlw/v2 .  
cd v2
ls
```

    [0m[01;32mBN_mp-1639_primitive.cif[0m   STRU                     [01;32mget_kpath.py[0m
    [01;32mGaN_mp-804_primitive.cif[0m   [01;32mSi_mp-149_primitive.cif[0m
    [01;32mMgO_mp-1265_primitive.cif[0m  cif2STRU.py
    

产生K点路径的code


```bash
cat get_kpath.py
```

    import os,glob
    from ase.io import read, write
    from pathlib import Path
    from ase.calculators.abacus import Abacus,AbacusTemplate
    from ase.geometry.dimensionality import (analyze_dimensionality,isolate_components)
    
    ####################################################################################################
    # This script is used to generate KPATH for 2D and 3D materials automatically.
    # Users should provide ABACUS STRU files in names of 'STRU_*'
    # the directory of STRU files is specified by `stru_dir`
    # The output KPT files are in name of KPT_STRU_* for each STRU_*, respectively.
    # insert_num determine how many kpoints to be inserted between 2 high-symmetry points
    # Author: Tianqi Zhao @2023/01/12
    ####################################################################################################
    
    def write_kpt(f1=None,kpath=None,insert_num=10):
        num_highsym_kpt = 0
        kpt_new = []
        insert_list = []
        k_path = list(kpath.path)
        for i in range(len(k_path)):
            if (k_path[i] != ','):
                kpt_new.append(kpath.special_points[k_path[i]])
                if i == len(k_path) -1:
                    insert_list.append(1)
                elif k_path[i+1] == ',':
                     insert_list.append(1)
                else:
                    insert_list.append(insert_num)
        lines = []
        lines.append('K_POINTS')
        lines.append(f'{len(kpt_new)}')
        lines.append(f'Line')
        for i in range(len(kpt_new)):
            lines.append(f'{kpt_new[i][0]:0<12f} {kpt_new[i][1]:0<12f} {kpt_new[i][2]:0<12f} {insert_list[i]}')
        lines.append('')
        f1 = open(kpt_dir+'KPT_'+stru_name,'w')
        f1.write('\n'.join(lines))
        f1.close()
    
    insert_num = 20
    stru_dir = './'
    kpt_dir = stru_dir
    stru_files = glob.glob(stru_dir+'STRU_*')
    stru_list = []
    stru_dim = {}
    atoms_dict = {}
    for i in range(len(stru_files)):
        stru_list.append(stru_files[i].split('/')[-1])
        atoms_from_stru = read(stru_files[i],format='abacus')
        atoms_dict[stru_list[i]] = atoms_from_stru
        intervals = analyze_dimensionality(atoms_from_stru,method='RDA')
        #intervals = analyze_dimensionality(atoms_from_stru,method='TSA')
        m = intervals[0]
        m.dimtype
        stru_dim[stru_list[i]] = m.dimtype
    
    for stru_name,dim_type in stru_dim.items():
        if dim_type == '3D':
            kpath = atoms_dict[stru_name].cell.bandpath(npoints=100)
            print(stru_name,dim_type,kpath)
            write_kpt(kpt_dir+'KPT_'+stru_name,kpath=kpath,insert_num=insert_num)
        elif dim_type == '2D':
            lat_length0 = atoms_dict[stru_name].cell.lengths()
            lat_length1 = []
            result = isolate_components(atoms_dict[stru_name],kcutoff=1.5)
            for dim,components in result.items():
                for atoms in components:
                    lat_length1 = atoms.cell.lengths()
            for i in range(len(lat_length0)):
                lat_length0[i] = round(lat_length0[i],5)
                lat_length1[i] = round(lat_length1[i],5)
            pbc = []
            for i in range(len(lat_length0)):
                if lat_length0[i] in lat_length1:
                    pbc.append(1)
                else:
                    pbc.append(0)
            lat2d_pbc = atoms_dict[stru_name].cell.get_bravais_lattice(pbc=pbc)
            kpath2d = atoms_dict[stru_name].cell.bandpath(npoints=100,pbc=pbc)
            print(stru_name,dim_type,kpath2d)
            write_kpt('KPT_'+stru_name,kpath=kpath2d,insert_num=insert_num)
    

## 3. 流程

针对计算的体系，可以从数据库[materials project](https://legacy.materialsproject.org/)下载对应的cif文件, 这里以Si、BN、MgO和GaN为例，cif文件提前在数据集中已下载好。


```bash
cat cif2STRU.py
```

    from ase.io import read, write
    from pathlib import Path
    import sys
    
    def cif2stru(cs_dir, cif):
        #cs_dir = '/root/v1'
        #cs_vasp = Path(cs_dir, 'Si_mp-149_primitive.cif')
        cs_vasp = Path(cs_dir, cif)
        cs_atoms = read(cs_vasp, format='cif')
        cs_stru = Path(cs_dir, 'STRU')
        pp = {}
        basis = {}
        write(cs_stru, cs_atoms, format='abacus', pp=pp, basis=basis)
    
    
    if __name__=='__main__':
        os_dir = sys.argv[1]
        cif  = sys.argv[2]
        cif2stru(os_dir, cif)
    

将cif文件转为STRU，传入两个参数，第一个是路径名，第二个是cif文件名


```bash
python cif2STRU.py . Si_mp-149_primitive.cif
```


```bash
ls
```

    [0m[01;32mBN_mp-1639_primitive.cif[0m   STRU                     [01;32mget_kpath.py[0m
    [01;32mGaN_mp-804_primitive.cif[0m   [01;32mSi_mp-149_primitive.cif[0m
    [01;32mMgO_mp-1265_primitive.cif[0m  cif2STRU.py
    


```bash
cp STRU STRU_Si
```

注意：计算能带路径往往仅限于晶体结构（结构高度对称且有序），**对于非晶、液体、孤立分子往往对称性差，结构有序度很低，一般不计算能带**，可以通过计算态密度来分析电子结构。**本code只适用于二维/三维具有对称性的晶体结构能带计算**，对于1维分子k点路径只需要0 0 0到0.5, 0 0，孤立分子用1 1 1。


```bash
python get_kpath.py && ls
```

    STRU_Si 3D BandPath(path='GXWKGLUWLK,UX', cell=[3x3], special_points={GKLUWX}, kpts=[100x3])
    [0m[01;32mBN_mp-1639_primitive.cif[0m  [01;32mMgO_mp-1265_primitive.cif[0m  [01;32mSi_mp-149_primitive.cif[0m
    [01;32mGaN_mp-804_primitive.cif[0m  STRU                       cif2STRU.py
    KPT_STRU_Si               STRU_Si                    [01;32mget_kpath.py[0m
    


```bash
cat KPT_STRU_Si
```

    K_POINTS
    12
    Line
    0.0000000000 0.0000000000 0.0000000000 20
    0.5000000000 0.0000000000 0.5000000000 20
    0.5000000000 0.2500000000 0.7500000000 20
    0.3750000000 0.3750000000 0.7500000000 20
    0.0000000000 0.0000000000 0.0000000000 20
    0.5000000000 0.5000000000 0.5000000000 20
    0.6250000000 0.2500000000 0.6250000000 20
    0.5000000000 0.2500000000 0.7500000000 20
    0.5000000000 0.5000000000 0.5000000000 20
    0.3750000000 0.3750000000 0.7500000000 1
    0.6250000000 0.2500000000 0.6250000000 20
    0.5000000000 0.0000000000 0.5000000000 1
    

类似地，我们可以产生其他体系的K点路径


```bash
python cif2STRU.py . BN_mp-1639_primitive.cif && cp STRU STRU_BN
python cif2STRU.py . GaN_mp-804_primitive.cif && cp STRU STRU_GaN
python cif2STRU.py . MgO_mp-1265_primitive.cif && cp STRU STRU_MgO
```


```bash
python get_kpath.py && ls
```

    STRU_GaN 3D BandPath(path='GMKGALHA,LM,KH', cell=[3x3], special_points={AGHKLM}, kpts=[100x3])
    STRU_BN 3D BandPath(path='GXWKGLUWLK,UX', cell=[3x3], special_points={GKLUWX}, kpts=[100x3])
    STRU_Si 3D BandPath(path='GXWKGLUWLK,UX', cell=[3x3], special_points={GKLUWX}, kpts=[100x3])
    STRU_MgO 3D BandPath(path='GXWKGLUWLK,UX', cell=[3x3], special_points={GKLUWX}, kpts=[100x3])
    [0m[01;32mBN_mp-1639_primitive.cif[0m  KPT_STRU_Si                STRU_MgO
    [01;32mGaN_mp-804_primitive.cif[0m  [01;32mMgO_mp-1265_primitive.cif[0m  STRU_Si
    KPT_STRU_BN               STRU                       [01;32mSi_mp-149_primitive.cif[0m
    KPT_STRU_GaN              STRU_BN                    cif2STRU.py
    KPT_STRU_MgO              STRU_GaN                   [01;32mget_kpath.py[0m
    

最后得到`KPT_STRU_BN`, `KPT_STRU_MgO`, `KPT_STRU_GaN`, 可以直接用于能带计算.
