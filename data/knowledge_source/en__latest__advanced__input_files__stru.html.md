# http://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html

> Source: http://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html

# The STRU file[#](#the-stru-file)

## Examples[#](#examples)

`STRU`

file contains the information about the lattice geometry, the name(s) and/or location(s) of the pseudopotential and numerical orbital files, as well as the structural information about the system. We supply two ways of specifying the lattice geometry. Below are two examples of the `STRU`

file for the same system:

### No latname[#](#no-latname)

`latname`

in the INPUT file. (See [input parameters](input-main.html#latname).)

```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf upf201 // label; mass; pseudo_file; pseudo_type
NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
LATTICE_CONSTANT
10.2 // lattice scaling factor (Bohr)
LATTICE_VECTORS
0.5 0.5 0.0 // latvec1
0.5 0.0 0.5 // latvec2
0.0 0.5 0.5 // latvec3
ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 // magnetism(Be careful: value 1.0 refers to 1.0 bohr mag, but not fully spin up !!!)
2 // number of atoms
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```

### latname fcc[#](#latname-fcc)

`latname="fcc"`

in the INPUT file. (See [input parameters](input-main.html#latname).) And the `STRU`

file becomes:

```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file
NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
LATTICE_CONSTANT
10.2 // lattice scaling factor (Bohr)
ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 // magnetism
2 // number of atoms
0.00 0.00 0.00 0 0 0//the position of atoms and other parameter specify by key word
0.25 0.25 0.25 1 1 1
```

The LATTICE_VECTORS section is removed.

## Structure of the file[#](#structure-of-the-file)

`STRU`

file contains several sections, and each section must start with a keyword like `ATOMIC_SPECIES`

, `NUMERICAL_ORBITAL`

, or `LATTICE_CONSTANT`

, etc. to signify what type of information that comes below.

### ATOMIC_SPECIES[#](#atomic-species)

`Si`

:

```
Si 28.00 Si_ONCV_PBE-1.0.upf upf201 // label; mass; pseudo_file; pseudo_type
```

`Si_ONCV_PBE-1.0.upf`

is the pseudopotential file. When the path is not specified, the file is assumed to be located in work directory. Otherwise, please explicitly specify the location of the pseudopotential files.

`upf201`

is the type of pseudopotential. There are five options: `upf`

(.UPF format), `upf201`

(the new .UPF format), `vwr`

(.vwr format), `blps`

(bulk-derived local pseudopotential), and `auto`

(automatically identified). If no pseudopotential type is assigned, the default value is `auto`

, and the pseudopotential type will be automatically identified.

When [esolver_type](input-main.html#esolver-type) is set to `lj`

or `dp`

, the keyword `pseudo_file`

and `pseudo_type`

is needless.

[dft_functional](input-main.html#dft-functional) keyword.

Common sources of the pseudopotential files include:

### NUMERICAL_ORBITAL[#](#numerical-orbital)

`LCAO`

calculations. Thus this section will be neglected in calcultions with plane wave basis. In the above example, numerical atomic orbitals is specified for the element `Si`

:

```
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
```

Numerical atomic orbitals may be downloaded from the [official website](http://abacus.ustc.edu.cn/pseudo/list.htm).

### LATTICE_CONSTANT[#](#lattice-constant)

The lattice constant of the system in unit of Bohr.

### LATTICE_VECTORS[#](#lattice-vectors)

*the lattice vectors given here are scaled by the lattice constant*. This section must be removed if the type Bravais lattice is specified using the input parameter `latname`

. (See [input parameters](input-main.html#latname).)

### LATTICE_PARAMETERS[#](#lattice-parameters)

`latname`

(see [input parameters](input-main.html#latname)) is used to specify the Bravais lattice type. The example above is a fcc lattice, where no additional information except the lattice constant is required to determine the geometry of the lattice.

`LATTICE_PARAMETERS`

must be present. It contains **one single line** with some parameters (separated by blank space if multiple parameters are needed), where the number of parameters required depends on specific type of lattice.

latname = “sc”: the LATTICE_PARAMETERS section is not required:

v1 = (1, 0, 0) v2 = (0, 1, 0) v3 = (0, 0, 1)

latname = “fcc”: the LATTICE_PARAMETERS section is not required:

v1 = (-0.5, 0, 0.5) v2 = (0, 0.5, 0.5) v3 = (-0.5, 0.5, 0)

latname = “bcc” : the LATTICE_PARAMETERS section is not required:

v1 = (0.5, 0.5, 0.5) v2 = (-0.5, 0.5, 0.5) v3 = (-0.5, -0.5, 0.5)

v1 = (1.0, 0, 0) v2 = (-0.5, sqrt(3)/2, 0) v3 = (0, 0, x)

v1 = (tx, -ty, tz) v2 = (0, 2ty, tz) v3 = (-tx, -ty, tz)

where tx=sqrt((1-x)/2), ty=sqrt((1-x)/6), and tz=sqrt((1+2x)/3).

v1 = (1, 0, 0) v2 = (0, 1, 0) v3 = (0, 0, x)

v1 = (0.5, -0.5, x) v2 = (0.5, 0.5, x) v3 = (-0.5, -0.5, x)

v1 = (1, 0, 0) v2 = (0, x, 0) v3 = (0, 0, y)

v1 = (0.5, x/2, 0) v2 = (-0.5, x/2, 0) v3 = (0, 0, y)

v1 = (0.5, 0, y/2) v2 = (0.5, x/2, 0) v3 = (0, x/2, y/2)

v1 = (0.5, x/2, y/2) v2 = (-0.5, x/2, y/2) v3 = (-0.5, -x/2, y/2)

v1 = (1, 0, 0) v2 = (x*z, x*sqrt(1-z^2, 0) v3 = (0, 0, y)

v1 = (0.5, 0, -y/2) v2 = (x*z, x*sqrt(1-z^2), 0) v3 = (0.5, 0, y/2)

v1 = (1, 0, 0) v2 = (x*m, x*sqrt(1-m^2), 0) v3 = (y*n, y*(l-n*m/sqrt(1-m^2)), y*fac)

where \(fac=\frac{\sqrt{1+2*m*n*l-m^2 -n^2 -l^2 }}{\sqrt{1-m^2}}\)


### ATOMIC_POSITIONS[#](#atomic-positions)

This section specifies the positions and other information of individual atoms.

The first line signifies method that atom positions are given, the following options are supported:

`Direct`

: coordinates of atom positions below would in fraction coordinates.`Cartesian`

: Cartesian coordinates in unit of ‘LATTICE_CONSTANT’.`Cartesian_au`

: Cartesian coordinates in unit of Bohr, same as setting of`Cartesian`

with`LATTICE_CONSTANT`

= 1.0 .`Cartesian_angstrom`

: Cartesian coordinates in unit of Angstrom, same as setting of`Cartesian`

with`LATTICE_CONSTANT`

= 1.889726125457828 .`Cartesian_angstrom_center_xy`

: Cartesian coordinates in unit of Angstrom, with Direct coordinate (0.5, 0.5, 0.0) as reference.`Cartesian_angstrom_center_xz`

: Cartesian coordinates in unit of Angstrom, with Direct coordinate (0.5, 0.0, 0.5) as reference…`Cartesian_angstrom_center_yz`

: Cartesian coordinates in unit of Angstrom, with Direct coordinate (0.0, 0.5, 0.5) as reference…`Cartesian_angstrom_center_xyz`

: Cartesian coordinates in unit of Angstrom, with Direct coordinate (0.5, 0.5, 0.5) as reference…

`Fe`

), the initial magnetic moment (`1.0`

), and the number of atoms for this particular element (`2`

) repsectively. Notice this magnetic moment will be a default value for every atom of this type but will be overrided if one define it for each atom by keyword(see below).

### More Key Words[#](#more-key-words)

Several other parameters could be defined after the atom position using key words :

`m`

or NO key word: three numbers, which take value in 0 or 1, control how the atom move in geometry relaxation calculations. In example below, the numbers`0 0 0`

following the coordinates of the first atom means this atom are*not allowed*to move in all three directions, and the numbers`1 1 1`

following the coordinates of the second atom means this atom*can*move in all three directions.`v`

or`vel`

or`velocity`

: set the three components of initial velocity of atoms in geometry relaxation calculations(e. g.`v 1.0 1.0 1.0`

).`mag`

or`magmom`

: set the start magnetization for each atom. In colinear case only one number should be given. In non-colinear case one have two choice:either set one number for the norm of magnetization here and specify two polar angle later(e. g. see below), or set three number for the xyz commponent of magnetization here (e. g.`mag 0.0 0.0 1.0`

). Note that if this parameter is set, the initial magnetic moment setting in the second line will be overrided.`angle1`

: in non-colinear case, specify the angle between z-axis and real spin, in angle measure instead of radian measure`angle2`

: in non-colinear case, specify angle between x-axis and real spin in projection in xy-plane , in angle measure instead of radian measuree.g.:

`nspin==2 || nspin == 4`

), e.g.:Fe 0.0 2 0.0 0.0 0.0 m 0 0 0 0.5 0.5 0.5 m 1 1 1 O 0.0 2 0.0 0.0 0.0 m 0 0 0 0.5 0.5 0.5 m 1 1 1

For

`nspin==2`

, we will autoset atomic magmon is`1.0`

:Fe 1.0 2 0.0 0.0 0.0 m 0 0 0 0.5 0.5 0.5 m 1 1 1 Fe 1.0 2 0.0 0.0 0.0 m 0 0 0 0.5 0.5 0.5 m 1 1 1

For

`nspin==4`

, we will autoset atomic magmon as follow:However, this autoset will not be vaild once

`STRU`

specalize a finite magnetic for any single atom.
