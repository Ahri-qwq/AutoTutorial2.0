# https://abacus.deepmodeling.com/en/latest/advanced/interface/deeph.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/interface/deeph.html

# DeepH[#](#deeph)

[DeepH](https://doi.org/10.1038/s43588-022-00265-6) applies meaching learning to predict the Hamiltonian in atomic basis representation. For such purpose, DeepH uses the Hamiltonian and overlap matrices from DFT calculations. Here we introduce how to extract relevant information from ABACUS for the purpose of DeepH training and prediction.

[website](https://deeph-pack.deepmodeling.com/en/latest/#deeph). An [example](https://deeph-pack.deepmodeling.com/en/latest/demo/demo3.html) for using DeepH with ABACUS is also provided.

`INPUT`

files.

Note: Use the LCAO basis for DeepH-related calculations


[README.md](http://README.md) file in the above-mentioned example, there are two stages where users need to run ABACUS calculations.

`INPUT`

file:

```
out_mat_hs2 1
```

`${x}`

.csr and data-SR-sparse_SPIN`${x}`

.csr will be generated, which contain the Hamiltonian and overlap matrices respectively in csr format. `${x}`

takes value of 0 or 1, based on the spin component. More details on this keyword can be found in the [list of input keywords](../input_files/input-main.html#out-mat-hs2).

`INPUT`

file we need to make the following specification of the keyword `calculation`

:

```
calculation get_S
```

A file named `SR.csr`

will be generated in the working directory, which contains the overlap matrix.
