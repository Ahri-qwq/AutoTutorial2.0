# http://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html

> Source: http://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html

# The KPT file[#](#the-kpt-file)

`STRU`

file. For the input k-point (`KPT`

) file, the file should either contain the k-point coordinates and weights or the mesh size for creating the k-point gird. Both options are allowed in `ABACUS`

.

## Gamma-only Calculations[#](#gamma-only-calculations)

[gamma_only](input-main.html#gamma-only) to be 1. Due to details of implementation, gamma-only calculation will be slightly faster than running a non gamma-only calculation and explicitly setting gamma point to be the only the k-point, but the results should be consistent.

## Generate k-mesh automatically[#](#generate-k-mesh-automatically)

`KPT`

) file used in `ABACUS`

.

```
K_POINTS //keyword for start
0 //total number of k-point, `0' means generate automatically
Gamma //which kind of Monkhorst-Pack method, `Gamma' or `MP'
2 2 2 0 0 0 //first three number: subdivisions along reciprocal vectors
//last three number: shift of the mesh
```

`K_POINTS`

, or `KPOINTS`

or just `K`

. The second line is an integer, and its value determines how to get k-points. In this example, `0`

means using Monkhorst-Pack (MP) method to generate k-points automatically.

`Gamma`

or `MP`

, different Monkhorst Pack (MP) method. Monkhorst-Pack (MP) is a method which uses the uniform k-points sampling in Brillouin-zone, while `Gamma`

means the Γ-centered Monkhorst-Pack method. The first three numbers of the last line are integers, which give the MP k grid dimensions, and the rest three are real numbers, which give the offset of the k grid. In this example, the numbers `0 0 0`

means that there is no offset, and this is the a standard 2by2by2 k grid.

## Set k-points explicitly[#](#set-k-points-explicitly)

```
K_POINTS //keyword for start
8 //total number of k-point
Direct //`Direct' or `Cartesian' coordinate
0.0 0.0 0.0 0.125 //coordinates and weights
0.5 0.0 0.0 0.125
0.0 0.5 0.0 0.125
0.5 0.5 0.0 0.125
0.0 0.0 0.5 0.125
0.5 0.0 0.5 0.125
0.0 0.5 0.5 0.125
0.5 0.5 0.5 0.125
```

## Band structure calculations[#](#band-structure-calculations)

```
K_POINTS // keyword for start
6 // number of high symmetry lines
Line // line-mode
0.5 0.0 0.5 20 // X
0.0 0.0 0.0 20 // G
0.5 0.5 0.5 20 // L
0.5 0.25 0.75 20 // W
0.375 0.375 0.75 20 // K
0.0 0.0 0.0 1 // G
```
