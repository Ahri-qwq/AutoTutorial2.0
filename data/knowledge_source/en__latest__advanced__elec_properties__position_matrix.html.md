# http://abacus.deepmodeling.com/en/latest/advanced/elec_properties/position_matrix.html

> Source: http://abacus.deepmodeling.com/en/latest/advanced/elec_properties/position_matrix.html

# Extracting Position Matrices[#](#extracting-position-matrices)

[out_mat_r](../input_files/input-main.html#out-mat-r) to write the position matrices into a file named `data-rR-tr`

in the directory `OUT.${suffix}`

. The position matrices is defined as:

[gamma_only](../input_files/input-main.html#gamma-only), and explicitly specifies that gamma point is the only k point in the KPT file.

```
-5 -5 -5 //R (lattice vector)
...
-5 -5 -4 //R (lattice vector)
...
-5 -5 -3 //R (lattice vector)
```

[out_app_flag](../input_files/input-main.html#out-app-flag) is set to true, then `data-rR-tr`

is written in an append manner. Otherwise, output files will be put in a separate directory, `matrix`

, and named as `$x`

_data-rR-tr, where `$x`

is the number of MD step. In addition, the output frequency is controlled by [out_freq_ion](../input_files/input-main.html#out-freq-ion). For example, if we are running a 10-step MD with out_freq_ion = 3, then `$x`

will be 0, 3, 6, and 9.

## get_s[#](#get-s)

`INPUT`

file we need to set the keyword [calculation](../input_files/input-main.html#calculation) to `get_s`

, and [out_mat_r](../input_files/input-main.html#out-mat-r) to `true`

.
