# https://abacus.deepmodeling.com/en/latest/advanced/scf/advanced.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/scf/advanced.html

# SCF in Complex Environments[#](#scf-in-complex-environments)

## Implicit Solvation Model[#](#implicit-solvation-model)

[methodology](https://aip.scitation.org/doi/10.1063/1.4865107) developed by Mathew, Sundararaman, Letchworth-Weaver, Arias, and Hennig in 2014.

[webpage](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#implicit-solvation-model):

```
INPUT_PARAMETERS
imp_sol 1
eb_k 80
tau 0.000010798
sigma_k 0.6
nc_k 0.00037
```

Example of running DFT calculation with the implicit solvation model is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/implicit_solvation_model/Pt-slab).

## External Electric Field[#](#external-electric-field)

`efield_flag`

in `INPUT`

(setting to 1 to turn on the field). Related keywords that control the external field are listed as follows with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#electric-field-and-dipole-correction):

```
INPUT_PARAMETERS
efield_flag 1
efield_dir 2
efield_pos_max 0.5
efield_pos_dec 0.1
efield_amp 0.001
```

[directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/electric_field/Pt-slab).

## Dipole Correction[#](#dipole-correction)

[methodology](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.59.12301) proposed by Bengtsson in 1999. This correction must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE. Note that the common input parameters shared between the external electric field and dipole correction, with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#electric-field-and-dipole-correction). The following keywords settings add dipole correction only without applying any external electric field:

```
INPUT_PARAMETERS
efield_flag 1
dip_cor_flag 1
efield_dir 2
efield_pos_max 0.5
efield_pos_dec 0.1
efield_amp 0
```

```
INPUT_PARAMETERS
efield_flag 1
dip_cor_flag 1
efield_dir 2
efield_pos_max 0.5
efield_pos_dec 0.1
efield_amp 0.001
```

[directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dipole_correction/Pt-slab). There are two input files, where `INPUT1`

considers only the dipole correction without no applied external field, while `INPUT2`

considers the dipole correction under an applied external field.

`INPUT`

, and run ABACUS.

## Compensating Charge[#](#compensating-charge)

[methodology](http://dx.doi.org/10.1103/PhysRevB.89.245406) developed by Brumme, Calandra, and Mauri in 2014. Input parameters that control the compensating charge are listed as follows with detailed explaination provided [here](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#gate-field-compensating-charge):

```
INPUT_PARAMETERS
gate_field 1
efield_dir 2
zgate 0.5
block 1
block_down 0.45
block_up 0.55
block_height 0.1
```

Example of running DFT calculation with the compensating charge is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/compensating_charge/Pt-slab).

## Van-der-Waals Correction[#](#van-der-waals-correction)

[D2](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20495), [D3(0)](https://aip.scitation.org/doi/10.1063/1.3382344) and [D3(BJ)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21759), to describe Van der Waals interactions. Among them, the D3 method has been implemented in ABACUS based on the dftd3 [program](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3) written by Stefan Grimme, Stephan Ehrlich and Helge Krieg.

To use VdW-correction, users need to supply value to the `vdw_method`

keyword in the `INPUT`

file:

(Default) none: no VdW correction

d2: DFT-D2 method

d3_0: DFT-D3(0) method

d3_bj: DFT-D3(BJ) method


[list of keywords](../input_files/input-main.html#vdw-correction) to control relevant parmeters used in calculating the VdW correction, such as the scale factor (s6) term. Recommended values of such parameters can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of the parameters in ABACUS are set to be the recommended values for PBE.

[directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/vdw/si2). There are two input files, where `INPUT1`

shows how to apply D2 correction with user-specified \(C_6\) parameter, and `INPUT2`

shows how to apply D3(BJ) correction with default VdW parameters.

`INPUT`

, and run ABACUS.
