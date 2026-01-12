# https://abacus.deepmodeling.com/en/latest/advanced/scf/hsolver.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/scf/hsolver.html

# Solving the Hamiltonian[#](#solving-the-hamiltonian)

## Explicit Diagonalization[#](#explicit-diagonalization)

Method of explicit solving KS-equation can be chosen by variable “ks_solver” in INPUT file.

`ks_solver`

can be `cg`

, `bpcg`

or `dav`

. The default setting `cg`

is recommended, which is band-by-band conjugate gradient diagonalization method. There is a large probability that the use of setting of `dav`

, which is block Davidson diagonalization method, can be tried to improve performance.

`ks_solver`

can be `genelpa`

or `scalapack_gvx`

. The default setting `genelpa`

is recommended, which is based on ELPA (EIGENVALUE SOLVERS FOR PETAFLOP APPLICATIONS) ([https://elpa.mpcdf.mpg.de/](https://elpa.mpcdf.mpg.de/)) and the kernel is auto choosed by GENELPA([pplab/GenELPA](https://github.com/pplab/GenELPA)), usually faster than the setting of “scalapack_gvx”, which is based on ScaLAPACK(Scalable Linear Algebra PACKage)

## Stochasic DFT[#](#stochasic-dft)

[Phys. Rev. B 106, 125132 (2022)](https://doi.org/10.1103/PhysRevB.106.125132)]. Different from traditional KSDFT with the explicit diagonalization method, SDFT and MDFT calculate physical quantities with trace of the corresponding operators. The advantages of SDFT and MDFT compared to the traditional KSDFT are the ability to simulate larger sizes and higher temperatures. In our package, SDFT and MDFT can be used by setting the `esolver_type`

parameter to `sdft`

for SCF calculations or MD calculations. To start with, you can refer to two [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/stochastic) and an explanation of the [input variables](../input_files/input-main.html#electronic-structure-sdft).

When we have a hamiltonian, the electronic density can be calculated with:

\(\rho(\mathbf{r})={\rm Tr}[f(\hat{H})\ket{\mathbf{r}}\bra{\mathbf{r}}]\),

`smearing_method`

, the parameter `smearing_sigma`

is equal the temperature \(T\) (in Ry) and `nche_sto`

represents the order of the expansion.

For physical quantities represented by operator \(\hat{O}\), SDFT calculates its trace with:

\({\rm Tr}[\hat{O}]=\sum_{i=1}^{N_\chi}{\bra{\chi_i}\hat{O}\ket{\chi_i}}\),

while MDFT calculates the trace as:

\(\ket{\tilde\chi_i}=\ket{\chi_i}-\sum_{n=1}^{N_\phi}\braket{\phi_n|\chi_i}\ket{\phi_n}\).

`nbands`

while the number of stochastic orbitals \(N_\chi\) is controlled by `nbands_sto`

.

`out_dos`

and `cal_cond`

separately.
