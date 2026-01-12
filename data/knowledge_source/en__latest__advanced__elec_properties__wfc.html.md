# https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/wfc.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/wfc.html

# Extracting Wave Functions[#](#extracting-wave-functions)

[examples/11_wfc](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc).

## Wave Function in G-Space[#](#wave-function-in-g-space)

`INPUT`

file while performing SCF calculation:

**PW basis**: Setto`out_wfc_pw`

`1`

. Output file format:`wfs[spin]k[kpoint]_pw.txt`

, where`[spin]`

is the spin channel index, and`[kpoint]`

the k-point index.**LCAO basis**: Setto`out_wfc_lcao`

`1`

.**Multi-k calculations**: Generates multiple files`wfs[spin]k[kpoint]_nao.txt`

.**Gamma-only calculations**:`wfs[spin]_nao.txt`

instead.


## Wave Function in Real Space[#](#wave-function-in-real-space)

[ out_wfc_norm](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-norm) or

[.](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-re-im)

`out_wfc_re_im`

[ basis_type](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#basis-type) is

`lcao`

, only `get_wf`

[is effective. An example is](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#calculation)

`calculation`

[examples/11_wfc/lcao_ienvelope_Si2](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc/lcao_ienvelope_Si2).
