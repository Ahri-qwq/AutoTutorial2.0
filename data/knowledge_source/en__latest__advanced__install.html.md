# https://abacus.deepmodeling.com/en/latest/advanced/install.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/install.html

# Advanced Installation Options[#](#advanced-installation-options)

[easy-installation guide](../quick_start/easy_install.html) before.

## Build with Libxc[#](#build-with-libxc)

Dependency: [Libxc](https://tddft.org/programs/libxc/) >= 5.1.7 .

`Libxc_DIR`

to the corresponding directory.

```
cmake -B build -DLibxc_DIR=~/libxc
```

## Build with ML-ALGO[#](#build-with-ml-algo)

C++ compiler, supporting

**C++14**(GCC >= 5 is sufficient)CMake >= 3.18

[LibTorch](https://pytorch.org/)with cxx11 ABI supporting CPU or GPU

```
cmake -B build -DENABLE_MLALGO=1 -DTorch_DIR=~/libtorch/share/cmake/Torch/ -Dlibnpy_INCLUDE_DIR=~/libnpy/include
```

CMake will try to download Libnpy if it cannot be found locally.


## Build with DeePMD-kit[#](#build-with-deepmd-kit)

[this page]for ABACUS interface with DP-GEN.

[TensorFlow](https://www.tensorflow.org/)(optional)[LibTorch](https://pytorch.org/)(optional)

`tensorflow_cc`

and `torch`

libraries are in the same directory as the `deepmd_c`

/`deepmd_cc`

libraries, then

```
cmake -B build -DDeePMD_DIR=/dir_to_deepmd-kit
```

DeePMD-kit supports TensorFlow backend but its libraries are placed at another directory, then

```
cmake -B build -DDeePMD_DIR=/dir_to_deepmd-kit -DTensorFlow_DIR=/dir_to_tensorflow
```

```
cmake -B build -DDeePMD_DIR=/dir_to_deepmd-kit -DTorch_DIR=/dir_to_pytorch
```

## Build with NEP[#](#build-with-nep)

[NEP_CPU](https://github.com/brucefan1983/NEP_CPU) library, which can be easily installed using toolchain as shown below:

```
./install_abacus_toolchain.sh --with-nep=install
```

To build ABACUS:

```
cmake -B build -DNEP_DIR=/path/to/nep_cpu
```

## Build with LibRI and LibComm[#](#build-with-libri-and-libcomm)

The new EXX implementation depends on two external libraries:

[deps](https://github.com/deepmodeling/abacus-develop/tree/develop/deps) folder. Set `-DENABLE_LIBRI=ON`

to build with these two libraries.

`-DLIBRI_DIR=${path to your LibRI folder} -DLIBCOMM_DIR=${path to your LibComm folder}`

.

## Build Unit Tests[#](#build-unit-tests)

`BUILD_TESTING`

flag. You can also specify path to local installation of [Googletest](https://github.com/google/googletest) by setting `GTEST_DIR`

flags. If not found in local, the configuration process will try to download it automatically.

```
cmake -B build -DBUILD_TESTING=1
```

After building and installing, unit tests can be performed with `ctest`

.

`ctest -R <test-match-pattern>`

to perform tests with name matched by given pattern.

## Build Performance Tests[#](#build-performance-tests)

`ENABLE_GOOGLEBENCH`

flag. You can also specify the path to a local installation of [Google Benchmark](https://github.com/google/benchmark.git) by setting `BENCHMARK_DIR`

flags. If not found locally, the configuration process will try to download it automatically.

```
cmake -B build -DENABLE_GOOGLEBENCH=1
```

`ENABLE_GOOGLEBENCH`

to ON, `BUILD_TESTING`

is automatically enabled. After building and installing, performance tests can be executed with `ctest`

.

## Build with CUDA support[#](#build-with-cuda-support)

### Extra prerequisites[#](#extra-prerequisites)

`USE_CUDA`

flag. You can also specify path to local installation of CUDA Toolkit by setting `CMAKE_CUDA_COMPILER`

flags.

```
cmake -B build -DUSE_CUDA=1 -DCMAKE_CUDA_COMPILER=${path to cuda toolkit}/bin/nvcc
```

`-DUSE_CUDA_MPI=ON`

. In this case, the program will directly communicate data with the CUDA hardware, rather than transferring it to the CPU first before communication. But note that if CUDA Aware is not supported, adding `-DUSE_CUDA_MPI=ON`

will cause the program to throw an error.

## Build math library from source[#](#build-math-library-from-source)

`USE_ABACUS_LIBM`

flag. It is expected to get a better performance on legacy versions of `gcc`

and `clang`

.

Currently supported math functions: `sin`

, `cos`

, `sincos`

, `exp`

, `cexp`


```
cmake -B build -DUSE_ABACUS_LIBM=1
```

## Build with PEXSI support[#](#build-with-pexsi-support)

[PEXSI Installation Guide](https://pexsi.readthedocs.io/en/latest/install.html) for more details. Note that PEXSI requires ParMETIS and SuperLU_DIST.

`ENABLE_PEXSI`

to `ON`

. If the libraries are not installed in standard paths, you can set `PEXSI_DIR`

, `ParMETIS_DIR`

and `SuperLU_DIST_DIR`

to the corresponding directories.

```
cmake -B build -DENABLE_PEXSI=ON -DPEXSI_DIR=${path to PEXSI installation directory} -DParMETIS_DIR=${path to ParMETIS installation directory} -DSuperLU_DIST_DIR=${path to SuperLU_DIST installation directory}
```

## Build ABACUS with make[#](#build-abacus-with-make)

Note: We suggest using CMake to configure and compile.


`make`

, users need to specify the location of the compilers, headers and libraries in `source/Makefile.vars`

:

```
# This is the Makefile of ABACUS API
#======================================================================
# Users set
#======================================================================
CXX = mpiicpc
# mpiicpc: compile intel parallel version
# icpc: compile intel sequential version
# make: ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: nothing need to be set except LIBXC_DIR
#
# mpicxx: compile gnu parallel version
# g++: compile gnu sequential version
# make: FFTW_DIR, OPENBLAS_LIB_DIR, SCALAPACK_LIB_DIR, ELPA_DIR, ELPA_INCLUDE_DIR, CEREAL_DIR must also be set.
# make pw: FFTW_DIR, OPENBLAS_LIB_DIR must be set.
# GPU = OFF #We do not support GPU yet
# OFF: do not use GPU
# CUDA: use CUDA
OPENMP = OFF
# the default is not to use OPENMP to accelerate.
# Change OPENMP to ON to use OPENMP.
#======================================================================
#-------------------- FOR INTEL COMPILER ----------------------------
## ELPA_DIR should contain an include folder and lib/libelpa.a
## CEREAL_DIR should contain an include folder.
#----------------------------------------------------------------------
ELPA_DIR = /usr/local/include/elpa-2021.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/elpa
CEREAL_DIR = /usr/local/include/cereal
##------------------- FOR GNU COMPILER ------------------------------
## FFTW_DIR should contain lib/libfftw3.a.
## OPENBLAS_LIB_DIR should contain libopenblas.a.
## SCALAPACK_LIB_DIR should contain libscalapack.a
## All three above will only be used when CXX=mpicxx or g++
## ELPA_DIR should contain an include folder and lib/libelpa.a
## CEREAL_DIR should contain an include folder.
##---------------------------------------------------------------------
# FFTW_DIR = /public/soft/fftw_3.3.8
# OPENBLAS_LIB_DIR = /public/soft/openblas/lib
# SCALAPACK_LIB_DIR = /public/soft/openblas/lib
# ELPA_DIR = /public/soft/elpa_21.05.002
# ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
# CEREAL_DIR = /public/soft/cereal
##------------------- OPTIONAL LIBS ---------------------------------
## To use MLALGO: set ENABLE_MLALGO = ON, and set LIBTORCH_DIR and LIBNPY_DIR
## To use LIBXC: set LIBXC_DIR which contains include and lib/libxc.a (>5.1.7)
## To use DeePMD: set DeePMD_DIR LIBTORCH_DIR and TensorFlow_DIR
## To use LibRI: set LIBRI_DIR and LIBCOMM_DIR
## To use PEXSI: set PEXSI_DIR DSUPERLU_DIR and PARMETIS_DIR
##---------------------------------------------------------------------
# LIBTORCH_DIR = /usr/local
# LIBNPY_DIR = /usr/local
ENABLE_MLALGO = OFF
# LIBXC_DIR = /public/soft/libxc
# DeePMD_DIR = ${deepmd_root}
# TensorFlow_DIR = ${tensorflow_root}
# LIBRI_DIR = /public/software/LibRI
# LIBCOMM_DIR = /public/software/LibComm
# PEXSI_DIR = /public/software/pexsi
# DSUPERLU_DIR = /public/software/superlu_dist
# PARMETIS_DIR = /public/software/parmetis
##---------------------------------------------------------------------
# NP = 14 # It is not supported. use make -j14 or make -j to parallelly compile
# DEBUG = OFF
# Only for developers
# ON: use gnu compiler and check segmental defaults
# OFF: nothing
#======================================================================
```

```
CXX = mpiicpc #(or CXX = icpc)
ELPA_DIR = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR = /public/soft/cereal
```

`CXX=mpiicpc`

, a parallel version will be compiled. When `CXX=icpc`

, a sequential version will be compiled.

Another example is where the Gnu C++ compiler, MPICH, OPENBLAS, ScaLAPACK, ELPA and CEREAL are used:

```
CXX = mpicxx/g++
FFTW_DIR = /public/soft/fftw_3.3.8
OPENBLAS_LIB_DIR = /public/soft/openblas/lib
SCALAPACK_LIB_DIR = /public/soft/openblas/lib
ELPA_DIR = /public/soft/elpa_21.05.002
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2021.05.002
CEREAL_DIR = /public/soft/cereal
```

`CXX=mpicxx`

, a parallel version will be compiled. When `CXX=g++`

, a sequential version will be compiled.

Except modifying `Makefile.vars`

, you can also directly use

```
make CXX=mpiicpc ELPA_DIR=/public/soft/elpa_21.05.002 \
ELPA_INCLUDE_DIR=${ELPA_DIR}/include/elpa-2021.05.002 \
CEREAL_DIR=/public/soft/cereal
```

`make`

or `make abacus`

to compile full version which supports LCAO calculations. Use `make pw`

to compile pw version which only supports pw calculations. For pw version, `make pw CXX=mpiicpc`

, you do not need to provide any libs. For `make pw CXX=mpicxx`

, you need provide `FFTW_DIR`

and `OPENBLAS_LIB_DIR`

.

After modifying the `Makefile.vars`

file, execute `make`

or `make -j12`

to build the program.

`ABACUS.mpi`

will be created in directory `bin/`

.

### Add Libxc Support[#](#add-libxc-support)

To compile ABACUS with LIBXC, you need to define `LIBXC_DIR`

in the file `Makefile.vars`

or use

```
make LIBXC_DIR=/pulic/soft/libxc
```

directly.

### Add ML-ALGO Support[#](#add-ml-algo-support)

`ENABLE_MLALGO = ON`

, and define `LIBTORCH_DIR`

and `LIBNPY_DIR`

in the file `Makefile.vars`

or use

```
make ENABLE_MLALGO=ON LIBTORCH_DIR=/opt/libtorch/ LIBNPY_DIR=/opt/libnpy/
```

directly.

### Add DeePMD-kit Support[#](#add-deepmd-kit-support)

[this page]for ABACUS interface with DP-GEN.

`DeePMD_DIR`

and `TensorFlow_DIR`

(TensorFlow Backend, optional) and/or `LIBTORCH_DIR`

(PyTorch Backend, optional) in the file `Makefile.vars`

.

`tensorflow_cc`

and `torch`

libraries are in the same directory as the `deepmd_c`

/`deepmd_cc`

libraries, then

```
make DeePMD_DIR=/dir_to_deepmd-kit
```

DeePMD-kit supports TensorFlow backend but its libraries are placed at another directory, then

```
make DeePMD_DIR=/dir_to_deepmd-kit TensorFlow_DIR=/dir_to_tensorflow
```

```
make DeePMD_DIR=/dir_to_deepmd-kit Torch_DIR=/dir_to_pytorch
```

### Add LibRI Support[#](#add-libri-support)

[LibRI](https://github.com/abacusmodeling/LibRI) and [LibComm](https://github.com/abacusmodeling/LibComm) and need to define `LIBRI_DIR`

and `LIBCOMM_DIR`

in the file `Makefile.vars`

or use

```
make LIBRI_DIR=/public/software/LibRI LIBCOMM_DIR=/public/software/LibComm
```

directly.
