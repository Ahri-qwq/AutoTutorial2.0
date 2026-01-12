# https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/band.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/band.html

# Extracting Band Structure[#](#extracting-band-structure)

[examples/band](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/band). Similar to the [DOS case](https://abacus-rtd.readthedocs.io/en/latest/advanced/elec_properties/dos.html), one first needs to perform a ground-state energy calculation * with one additional keyword “*:

[out_chg](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-chg)” in the INPUT file

```
out_chg 1
```

`chgs1.cube`

, `chgs2.cube`

, etc. Then, one can use the same `STRU`

file, pseudopotential files and atomic orbital files (and the local density matrix file [onsite.dm](http://onsite.dm) if DFT+U is used) to do a non-self-consistent (NSCF) calculation. In this example, the potential is constructed from the ground-state charge density from the proceeding calculation. Now the INPUT file is like:

```
INPUT_PARAMETERS
#Parameters (General)
nbands 8
calculation nscf
basis_type lcao
read_file_dir ./
#Parameters (Accuracy)
ecutwfc 60
scf_nmax 50
scf_thr 1.0e-9
pw_diag_thr 1.0e-7
#Parameters (File)
init_chg file
out_band 1
out_proj_band 1
#Parameters (Smearing)
smearing_method gaussian
smearing_sigma 0.02
```

Here is a relevant k-point file KPT (in LINE mode):

```
K_POINTS # keyword for start
6 # number of high symmetry lines
Line # line-mode
0.5 0.0 0.5 20 # X
0.0 0.0 0.0 20 # G
0.5 0.5 0.5 20 # L
0.5 0.25 0.75 20 # W
0.375 0.375 0.75 20 # K
0.0 0.0 0.0 1 # G
```

This means we are using the following k-points:

20/1 number of k points along the segment line, which is constructed by two adjacent k points.


`eigs1.txt`

in the output directory. Plot it and you will obtain the energy band structure!

`<band structure>`

is the same as data in the BANDS_1.dat file, such as:

```
<pband>
<nspin>1</nspin>
<norbitals>153</norbitals>
<band_structure nkpoints="96" nbands="50" units="eV">
...
```

The rest of the files arranged in sections, each section with a header such as below:

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

The shape of text contents of element `<data>`

is (Number of k-points, Number of bands)
