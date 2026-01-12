# https://abacus.deepmodeling.com/en/latest/advanced/interface/Hefei-NAMD.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/interface/Hefei-NAMD.html

# Hefei-NAMD[#](#hefei-namd)

[Hefei-NAMD](https://github.com/QijingZheng/Hefei-NAMD) Non-adiabatic molecular dynamics applies surface hopping to incorporate quantum mechanical effects into molecular dynamics simulations. Surface hopping partially incorporates the non-adiabatic effects by including excited adiabatic surfaces in the calculations, and allowing for ‘hops’ between these surfaces.

Detailed instructions on installing and running Hefei-NAMD can be found on its official [website](http://staff.ustc.edu.cn/~zqj/posts/Hefei-NAMD-Training/).

The steps are as follows :

Add output parameters in INPUT when running MD using ABACUS .


```
out_wfc_lcao 1
out_mat_hs 1
```

Clone Hefei-NAMD codes optimized for ABACUS from

[website](https://github.com/vtzf/abacus-namd).`Args.py`

including directory of ABACUS output files and NAMD parameters. We can see detailed explanation for all parameters in`Args.py`

.Run

`NAC.py`

to prepare related files for NAMD simulations.

```
sbatch sub_nac
```

Run

`SurfHop.py`

to perform NAMD simulations.

```
sbatch sub_sh
```

And results are under directory namddir in `Args.py`

.
