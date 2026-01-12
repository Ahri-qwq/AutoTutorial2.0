# https://abacus.deepmodeling.com/en/latest/advanced/scf/initialization.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/scf/initialization.html

# Initializing SCF[#](#initializing-scf)

## Charge Density[#](#charge-density)

`init_chg`

is used for choosing the method of charge density initialization.

`atomic`

: initial charge density by atomic charge density from pseudopotential file under keyword`PP_RHOATOM`

`file`

: initial charge density from files produced by previous calculations with.`out_chg 1`

`auto`

: Abacus first attempts to read the density from a file; if not found, it defaults to using atomic density.

## Wave function[#](#wave-function)

`init_wfc`

is used for choosing the method of wavefunction coefficient initialization.

`basis_type=pw`

, setting of `random`

and `atomic`

are supported. Atomic wave function is read from pseudopotential file under keyword `PP_PSWFC`

, if setting is `atomic`

and number of band of atomic wavefunction less than `nbands`

in INPUT file, the extra bands will be initialed by random.

`basis_type=lcao`

, we further support reading of initial wavefunction by setting `init_wfc`

to `file`

. In LCAO code, wave function is used to initialize density matrix and real-space charge density. For such purpose, a file containing wavefunction must be prepared. Such files can be generated from previous calculations with [ out_wfc_lcao 1](../elec_properties/wfc.html).
