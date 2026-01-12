# https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/dos.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/dos.html

# Calculating DOS and PDOS[#](#calculating-dos-and-pdos)

## DOS[#](#dos)

[examples/dos](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dos). We first, do a ground-state energy calculation * with one additional keyword “*:

[out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg)” in the INPUT file

```
out_chg 1
```

`STRU`

file, pseudopotential file and atomic orbital file (and the local density matrix file [onsite.dm](http://onsite.dm) if DFT+U is used) to do a non-self-consistent calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:

```
INPUT_PARAMETERS
#Parameters (General)
suffix Si2_diamond
ntype 1
nbands 8
calculation nscf
basis_type lcao
read_file_dir ./
#Parameters (Accuracy)
ecutwfc 60
symmetry 1
scf_nmax 50
scf_thr 1.0e-9
pw_diag_thr 1.0e-7
#Parameters (File)
init_chg file
out_dos 1
dos_sigma 0.07
```

Some parameters in the INPUT file are explained:

calculation

[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#calculation).pw_diag_thr

[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#pw_diag_thr).For LCAO calculations, this parameter will be neglected !

init_chg

[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#init_chg).out_dos

`(number of states)/(eV * unitcell)`

. For its more information please see the[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#out_dos).dos_sigma

the gaussian smearing parameter(DOS), in unit of eV. For its more information please see the

[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#dos_sigma).read_file_dir

the location of electron density file. For its more information please see the

[here](https://abacus-rtd--1282.org.readthedocs.build/en/1282/advanced/input_files/input-main.html#read_file_dir).

```
K_POINTS
0
Gamma
8 8 8 0 0 0
```

```
-5.49311 0.0518133 0.0518133
-5.48311 0.0641955 0.116009
-5.47311 0.0779299 0.193939
-5.46311 0.0926918 0.28663
-5.45311 0.108023 0.394653
-5.44311 0.123346 0.517999
...
```

## PDOS[#](#pdos)

```
<pdos>
<nspin>1</nspin>
<norbitals>26</norbitals>
<energy_values units="eV">
-5.50311
-5.49311
-5.48311
-5.47311
...
```

The rest of the fileis arranged in sections, each section with a header such as below:

```
<orbital
index=" 1"
atom_index=" 1"
species="Si"
l=" 0"
m=" 0"
z=" 1"
>
<data>
...
</data>
```

`(number of states)/(eV * unitcell)`

.
