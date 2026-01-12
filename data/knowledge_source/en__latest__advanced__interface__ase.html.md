# https://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/interface/ase.html

# ASE[#](#ase)

## Introduction[#](#introduction)

[ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment) provides a set of Python tools for setting, running, and analysing atomic simulations. We have developed an ABACUS calculator ([ase-abacus](https://gitlab.com/1041176461/ase-abacus)) to be used together with the ASE tools, which exists as an external project with respect to ASE and is maintained by ABACUS developers.

## Installation[#](#installation)

```
git clone https://gitlab.com/1041176461/ase-abacus.git
cd ase-abacus
pip install .
```

Another direct way:

```
pip install git+https://gitlab.com/1041176461/ase-abacus.git
```

## Environment variables[#](#environment-variables)

[ABACUS](http://abacus.ustc.edu.cn) supports two types of basis sets: PW, LCAO. The path of pseudopotential and numerical orbital files can be set throught the environment variables `ABACUS_PP_PATH`

and `ABACUS_ORBITAL_PATH`

, respectively, e.g.:

```
PP=${HOME}/pseudopotentials
ORB=${HOME}/orbitals
export ABACUS_PP_PATH=${PP}
export ABACUS_ORBITAL_PATH=${ORB}
```

`ABACUS_PP_PATH`

is needed. For LCAO calculations, both `ABACUS_PP_PATH`

and `ABACUS_ORBITAL_PATH`

should be set.

Also, one can manally set the paths of PP and ORB when using ABACUS calculator in ASE.

## ABACUS Calculator[#](#abacus-calculator)

The default initialization command for the ABACUS calculator is

```
from ase.calculators.abacus import Abacus
```

keyword | description |
|---|---|
| dict of pseudopotentials for involved elememts, |
|
|
|
|
|
|
|
|
|
|

[here](../input_files/input-main.html).

The input parameters can be set like::

```
# for ABACUS calculator
calc = Abacus(profile=profile,
ecutwfc=100,
scf_nmax=100,
smearing_method='gaussian',
smearing_sigma=0.01,
basis_type='pw',
ks_solver='dav',
calculation='scf',
pp=pp,
basis=basis,
kpts=kpts)
```

The command to run jobs can be set by specifying `AbacusProfile`

::

```
from ase.calculators.abacus import AbacusProfile
# for OpenMP setting inside python env
import os
os.environ("OMP_NUM_THREADS") = 1
# for MPI setting used in abacus
mpi_num = 4
# for ABACUS Profile
abacus = '/usr/local/bin/abacus' # specify abacus exec
profile = AbacusProfile(command=f'mpirun -n {mpi_num} {abacus}') # directly the command for running ABACUS
```

in which `abacus`

sets the absolute path of the `abacus`

executable.

## MD Analysis[#](#md-analysis)

`running_md.log`

can be read. If the ‘STRU_MD_*’ files are not continuous (e.g. ‘STRU_MD_0’, ‘STRU_MD_5’, ‘STRU_MD_10’…), the index parameter of read should be as a slice object. For example, when using the command `read('running_md.log', index=slice(0, 15, 5), format='abacus-out')`

to parse ‘running_md.log’, ‘STRU_MD_0’, ‘STRU_MD_5’ and ‘STRU_MD_10’ will be read.

The `MD_dump`

file is also supported to be read-in by `read('MD_dump', format='abacus-md')`


## SPAP Analysis[#](#spap-analysis)

[SPAP](https://github.com/chuanxun/StructurePrototypeAnalysisPackage) (Structure Prototype Analysis Package) is written by Dr. Chuanxun Su to analyze symmetry and compare similarity of large amount of atomic structures. The coordination characterization function (CCF) is used to measure structural similarity. An unique and advanced clustering method is developed to automatically classify structures into groups.

If you use this program and method in your research, please read and cite the publication:

and you should install it first with command `pip install spap`

.
