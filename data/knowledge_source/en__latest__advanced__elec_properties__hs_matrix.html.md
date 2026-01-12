# https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/hs_matrix.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/hs_matrix.html

# Extracting Hamiltonian and Overlap Matrices[#](#extracting-hamiltonian-and-overlap-matrices)

\(H(k)=\sum_R H(R)e^{-ikR}\)


and

\(S(k)=\sum_R S(R)e^{-ikR}\)


## out_mat_hs[#](#out-mat-hs)

[out_mat_hs](../input_files/input-main.html#out-mat-hs) to true to print the upper triangular part of the Hamiltonian matrices and overlap matrices for each k point into files in the directory `OUT.${suffix}`

. It is available for both gamma_only and multi-k calculations.

[Basis Set](../pp_orb.html#basis-set).

As for information on the k points, one may look for the `SETUP K-POINTS`

section in the running log.

## out_mat_hs2[#](#out-mat-hs2)

[out_mat_hs2](../input_files/input-main.html#out-mat-hs2). This functionality is not available for gamma_only calculations. To generate such matrices for gamma only calculations, users should turn off [gamma_only](../input_files/input-main.html#gamma-only), and explicitly specify that gamma point is the only k point in the KPT file.

`hrs1_nao.csr`

and `sr_nao.csr`

are generated, which contain the Hamiltonian matrix \(H(R)\) and overlap matrix \(S(R)\) respectively. For nspin = 2, three files `hrs1_nao.csr`

and `hrs2_nao.csr`

and `sr_nao.csr`

are created, where the first two files correspodn to \(H(R)\) for spin up and spin down, respectively.

`R`

are in the file.

`R`

and the number of nonzero matrix elements, such as:

```
-3 1 1 1020
```

which means there are 1020 nonzero elements in the (-3,1,1) cell.

[out_freq_ion](../input_files/input-main.html#out-freq-ion) and [out_app_flag](../input_files/input-main.html#out-app-flag).

## get_s[#](#get-s)

`INPUT`

file we need to set the value keyword [calculation](../input_files/input-main.html#calculation) to be `get_s`

.

`sr_nao.csr`

will be generated in the working directory, which contains the overlap matrix.

`nspin`

is set to 1 or 2, the dimension of the overlap matrix is nlocal \(\times\) nlocal, where nlocal is the total number of numerical atomic orbitals. These numerical atomic orbitals are ordered from outer to inner loop as atom, angular quantum number \(l\), zeta (multiple radial orbitals corresponding to each \(l\)), and magnetic quantum number \(m\). When`nspin`

is set to 4, the dimension of the overlap matrix is (2 \(\times\) nlocal) \(\times\) (2 \(\times\) nlocal). In this case, the numerical atomic orbitals are ordered from outer to inner loop as atom, angular quantum number \(l\), zeta (multiple radial orbitals corresponding to each \(l\)), magnetic quantum number \(m\), and npol (index of spin, ranges from 0 to 1).

## examples[#](#examples)

We provide [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/matrix_hs) of outputting the matrices. There are four examples:

out_hs_gammaonly: writing H(k) and S(k) for gamma-only calculation

out_hs_multik: writing H(k) and S(k) for multi-k calculation

out_hs2_multik: writing H® and S® for multi-k calculation

out_s_multik: running calculation=get_s to obtain overlap matrix for multi-k calculation


Reference output files are provided in each directory.
