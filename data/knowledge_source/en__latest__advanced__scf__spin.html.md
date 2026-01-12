# https://abacus.deepmodeling.com/en/latest/advanced/scf/spin.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/scf/spin.html

# Spin-polarization and SOC[#](#spin-polarization-and-soc)

## Non-spin-polarized Calculations[#](#non-spin-polarized-calculations)

**“nspin 1”** in INPUT file means calculation with non-polarized spin. In this case, electrons with spin up and spin down have same occupations at every energy states, weights of bands per k point would be double.

## Collinear Spin Polarized Calculations[#](#collinear-spin-polarized-calculations)

**“nspin 2”** in INPUT file means calculation with polarized spin along z-axis. In this case, electrons with spin up and spin down will be calculated respectively, number of k points would be doubled. Potential of electron and charge density will separate to spin-up case and spin-down case.

Magnetic moment Settings in [STRU files](../input_files/stru.html) are not ignored until **“nspin 2”** is set in INPUT file

When **“nspin 2”** is set, the screen output file will contain magnetic moment information. e.g.

```
ITER TMAG AMAG ETOT(eV) EDIFF(eV) DRHO TIME(s)
GE1 4.16e+00 4.36e+00 -6.440173e+03 0.000000e+00 6.516e-02 1.973e+01
```

[Mulliken charge analysis](../elec_properties/Mulliken.html).

### Constraint DFT for collinear spin polarized calculations[#](#constraint-dft-for-collinear-spin-polarized-calculations)

For some special need, there are two method to constrain electronic spin.

**“ocp”**and**“ocp_set”**If**“ocp=1”**and**“ocp_set”**is set in INPUT file, the occupations of states would be fixed by**“ocp_set”**, this method is often used for excited states calculation. Be careful that: when**“nspin=1”**, spin-up and spin-down electrons will both be set, and when**“nspin=2”**, you should set spin-up and spin-down respectively.**“nupdown”**If**“nupdown”**is set to non-zero, number of spin-up and spin-down electrons will be fixed, and Fermi energy level will split to E_Fermi_up and E_Fermi_down. By the way, total magnetization will also be fixed, and will be the value of**“nupdown”**.

## Noncollinear Spin Polarized Calculations[#](#noncollinear-spin-polarized-calculations)

**“noncolin 1”**, in which case the coupling between spin up and spin down will be taken into account. In this case, nspin is automatically set to 4, which is usually not required to be specified manually. The weight of each band will not change, but the number of occupied states will be double. If the nbands parameter is set manually, it is generally set to twice what it would be when nspin<4.

[SOC effects](#soc-effects). When **“lspinorb 1”** in INPUT file, “nspin” is also automatically set to 4. Note: different settings for “noncolin” and “lspinorb” correspond to different calculations:

noncolin=0 lspinorb=0 nspin<4 : Non-collinear magnetic moments and SOC effects are not considered.

noncolin=1 lspinorb=1 : The SOC effect and non-collinear magnetic moment are both calculated.


## For the continuation job[#](#for-the-continuation-job)

# SOC Effects[#](#soc-effects)

## SOC[#](#soc)

`lspinorb`

is used for control whether or not SOC(spin-orbit coupling) effects should be considered.

Both `basis_type=pw`

and `basis_type=lcao`

support `scf`

and `nscf`

calculation with SOC effects.

Atomic forces and cell stresses can not be calculated with SOC effects yet.

## Pseudopotentials and Numerical Atomic Orbitals[#](#pseudopotentials-and-numerical-atomic-orbitals)

For Norm-Conserving pseudopotentials, there are differences between SOC version and non-SOC version.

`PP_HEADER`

part, keyword `has_so=1`

and `relativistic="full"`

refer to SOC effects have been considered, which would lead to different `PP_NONLOCAL`

and `PP_PSWFC`

parts. Please be careful that `relativistic="full"`

version can be used for SOC or non-SOC calculation, but `relativistic="scalar"`

version only can be used for non-SOC calculation. When full-relativistic pseudopotential is used for non-SOC calculation, ABACUS will automatically transform it to scalar-relativistic version.

## Partial-relativistic SOC Effect[#](#partial-relativistic-soc-effect)

`soc_lambda`

, which has value range [0.0, 1.0] , is used for modulate SOC effect. In particular, `soc_lambda 0.0`

refers to scalar-relativistic case and `soc_lambda 1.0`

refers to full-relativistic case.
