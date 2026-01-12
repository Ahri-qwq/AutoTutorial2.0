# https://abacus.deepmodeling.com/en/latest/advanced/scf/construct_H.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/scf/construct_H.html

# Constructing the Hamiltonian[#](#constructing-the-hamiltonian)

## Exchange-Correlation Functionals[#](#exchange-correlation-functionals)

`dft_functional`

keyword in `INPUT`

file. If `dft_functional`

is not specified, ABACUS will use the xc functional indicated in the pseudopotential file.

[file](../../_downloads/06b3cf6c86225a8f6a3bf4b28215f1f1/xc_funcs.h) for a complete list of functionals implemented in ABACUS. Furthermore, if ABACUS is compiled with LIBXC, we also support all the LDA, GGA and meta-GGA functionals provided therein.

Here, we use a simple [example calculation](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/scf/lcao_Si2) for illustration.

**Default setting:**`INPUT`

file, there is no specification of the`dft_functional`

keyword. As a result, we use the default option, that is to use the xc functional in the pseudopotential file,`Si.pz-vbc.UPF`

. We can take a look at the first few lines of the`<PP_HEADER>`

section from the pseudopotential file:`dft_functional`

is specified, users should make sure that all pseudopotentials are using the same functional. Otherwise, the type of xc functional should be specified explicitly.**Using PBE**`dft_functional`

parameter. For example, to use PBE functional, add the following line to`INPUT`

file and rerun the calculation:dft_functional PBE

**More functionals from LIBXC**For this part, users should compile the ABACUS code with LIBXC linked (version 5.1.7 or higher).

To use SCAN functional, make the following modification to the

`INPUT`

file:dft_functional SCAN

[source code](../../_downloads/d54055d49c4276ef7a86279feb37a0f6/xc_functional.cpp).dft_functional LDA_X_YUKAWA+LDA_C_1D_CSC

The list of LIBXC keywords can be found on its

[website](https://www.tddft.org/programs/libxc/functionals/).**Temperature-dependent functional**`xc_temperature`

(unit is Rydberg) is used to specify the temperature, such as the following:dft_functional LDA_XC_CORRKSDT xc_temperature 10

**Hybrid functional**`dft_functional`

. Options are`hf`

(pure Hartree-Fock),`pbe0`

(PBE0),`hse`

, and`scan0`

(SCAN0) (Note: in order to use HSE or SCAN0 functional, LIBXC is required). Note also that only HSE has been tested while other hybrid functionals have NOT been fully tested yet, and the maximum parallel cpus for running exx is N^4, with N being the number of atoms.[Exact Exchange](../input_files/input-main.html#exact-exchange)in the list of input variables for more information.[directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/hse/lcao_Si2). Apart from the input files (`INPUT`

,`STRU`

,`KPT`

), we further provide two files: running_scf.log_ref and log_ref, which contains reference for running_scf.log and standard output from the program, respectively.

## DFT+*U*[#](#dft-u)

*d*/*f* shells. These include transition metals ™ and their oxides, rare-earth compounds, and actinides, to name a few, where L(S)DA/GGAs typically yield quantitatively or even qualitatively wrong results. To address this failure, an efficient and successful method named DFT+*U*, which inherits the efficiency of L(S)DA/GGA but gains the strength of the Hubbard model in describing the physics of strongly correlatedsystems, has been developed.

*U* method is accessible in ABACUS. The details of the DFT+*U* method could be found in this [paper](https://doi.org/10.1063/5.0090122). It should be noted that the DFT+*U* works only within the NAO scheme, which means that the value of the keyword `basis_type`

must be lcao when DFT+*U* is called. To turn on DFT+*U*, users need to set the value of the `dft_plus_u`

keyword in the `INPUT`

file to be 1. All relevant parmeters used in DFT+*U* calculations are listed in the [DFT+ U correction](../input_files/input-main.html#dft-u-correction) part of the

[list of keywords](../input_files/input-main.html).

Examples of DFT+*U* calculations are provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dft_plus_u).
