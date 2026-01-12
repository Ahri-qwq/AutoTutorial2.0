# https://abacus.deepmodeling.com/en/latest/advanced/pp_orb.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/pp_orb.html

# Basis Set and Pseudopotentials[#](#basis-set-and-pseudopotentials)

## Basis Set[#](#basis-set)

ABACUS supports both PW and LCAO basis set, controlled by keyword [basis_type](input_files/input-main.html#basis-type) in INPUT file.

[kinetic energy cutoff](input_files/input-main.html#ecutwfc) of the plane wave.

[official website](http://abacus.ustc.edu.cn/pseudo/list.htm). For more information, also check the `NUMERICAL_ORBITAL`

section in the specification of the [STRU file](input_files/stru.html).

[Table of spherical harmonics - Wikipedia](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics), by a factor of \((-1)^m\).

**except for the intra-m ordering**. Specifically, orbitals are first ordered by their atomic species in accordance with the `ATOMIC_SPECIES`

section of the STRU file. For orbitals of the same species, orbitals belonging to each atom are put together, with their overall order following the `ATOMIC_POSITIONS`

section of the STRU file. Orbitals on each atom are further ascendingly ordered by their angular momentum (s,p,d,f,…), followed by an order based on their their zeta number. Finally, m is ordered as 0, 1, -1, 2, 2, \(\ldots\), l, -l, which is the only exception of the lexicographic order.

## Generating atomic orbital bases[#](#generating-atomic-orbital-bases)

Guidelines for generating atomic orbital bases are as follows:

[Numerical Atomic Orbitals 1: the nomenclature and usage of numerical atomic orbitals in ABACUS](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html)(Chinese)[Numerical Atomic Orbitals 3: generate high-precision numerical atomic orbitals](https://mcresearch.github.io/abacus-user-guide/abacus-nac1.html)(Chinese)

[Github repository of ABACUS ORBGEN project](https://github.com/kirk0830/ABACUS-ORBGEN), the usage of which can be found in README (in English) file.

*NOTE*: users are encouraged to cite the above works when numerical atomic orbitals and its generation codes are used in their research.

## BSSE Correction[#](#bsse-correction)

`STRU`

file when an element name contains the “empty” suffix, such as “H_empty”, “O_empty” and so on. Here we provide an [example](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/bsse/water) of calculating the molecular formation energy of \(H_2O\) with BSSE correction.

In the example, we provide four STRU files:

STRU_0 : used along with ntype = 2;normal calculation of water molecule (\(E(\text{H}_2\text{O})\))

obtained total energy of -466.4838149140513 eV

STRU_1 : used along with ntype = 2;calculation of single O atom (\(E_O\))

obtained total energy of -427.9084406198214 eV

STRU_2 : used along with ntype = 3;calculation of 1st H atom (\(E_{H1}\))

obtained total energy of -12.59853381731160 eV

STRU_3 : used along with ntype = 3;calculation of 2nd H atom (\(E_{H2}\))

obtained total energy of -12.59853378720844 eV


Note : Remember to adjust the parameter

`ntype`

in INPUT file

Thus, the formation energy is given by:

## Pseudopotentials[#](#pseudopotentials)

### Supported formats[#](#supported-formats)

### Usage[#](#usage)

`ATOMIC_SPECIES`

section in the specification of the [STRU file](input_files/stru.html).

### Download[#](#download)

Users can find pseudopotentials in the following links:

**Website**

[Quantum ESPRESSO](https://www.quantum-espresso.org/pseudopotentials): the official website of Quantum ESPRESSO, where you can find a large number of pseudopotential files.[Stantard Solid State Pseudopotential library](https://www.materialscloud.org/sssp): a library of**high-quality**pseudopotentials for solid-state calculations, with**a large number of tests on efficiency and precison**.[PWmat](http://www.pwmat.com/potential-download): a website that provides a large number of pseudopotential files, various kinds of semi-core constructed pseudopotentials are included.[THEOS](http://theossrv1.epfl.ch/Main/Pseudopotentials): PSlibrary 0.3.1, a library of pseudopotentials for DFT calculations, including ultrasoft, paw, norm-conserving both full-relativistic and scalar-relativistic pseudopotentials.[ABACUS@USTC](https://abacus.ustc.edu.cn/pseudo/list.htm):**ABACUS official website**where you can find a large number of pseudopotential files and numerical atomic orbital files.[BLPS](https://github.com/PrincetonUniversity/BLPSLibrary): BLPS format pseudopotential library

**Norm-conserving pseudopotentials**

[SG15](http://www.quantum-simulation.org/potentials/sg15_oncv/):**vastly used in ABACUS**DFT calculation and numerical atomic orbital generation.[PseudoDOJO](http://www.pseudo-dojo.org/): another widely used pseudopotential database, developed by Abinit group,**including Lanthanide pseudopotentials (f-electrons frozen)**.[The Rappe group](https://www.sas.upenn.edu/rappegroup/research/pseudo-potential-gga.html): a collection of GGA pseudopotentials which are generated with Opium code, several tests proves that are out-performing in alloy systems.[Matteo Giantomassi’s Github repo](https://github.com/gmatteo/pseudos_ac_she): a Github repository that contains norm-conserving pseudopotentials for**Actinides and superheavy elements to 120-th element**.

**Ultrasoft pseudopotentials**

[Vanderbilt](http://www.physics.rutgers.edu/~dhv/uspp/): a collection of ultrasoft pseudopotentials generated by Vanderbilt group.[GBRV](https://www.physics.rutgers.edu/gbrv/)by Kevin F. Garrity, Joseph W. Bennett, Karin M. Rabe, and David Vanderbilt: presently the most popular ultrasoft pseudpotentials in Quantum ESPRESSO user community.

### Pseudopotential Generation[#](#pseudopotential-generation)

For pseudopotential generation, please refer to the following links for more information:

[A brief introduction of norm-conserving pseudopotential generation](https://mcresearch.github.io/abacus-user-guide/abacus-upf.html)
