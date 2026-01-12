# https://abacus.deepmodeling.com/en/latest/quick_start/easy_install.html

> Source: https://abacus.deepmodeling.com/en/latest/quick_start/easy_install.html

# Easy Installation[#](#easy-installation)

`make`

, please refer to [the advanced installation guide](../advanced/install.html)`cmake`

to avoid dependency issues. We recommend compiling ABACUS (and possibly its requirements) from the source code using the latest compiler for the best performace. You can use [toolchain](#install-by-toolchain) to install ABACUS and dependencies in a source-code compilation way with convience. You can also deploy ABACUS **without building** by [Docker](#container-deployment) or [conda](#install-by-conda). Please note that ABACUS only supports Linux; for Windows users, please consider using [WSL](https://learn.microsoft.com/en-us/windows/wsl/) or docker.

## Get ABACUS source code[#](#get-abacus-source-code)

ABACUS source code can be obtained via one of the following choices:

Clone the whole repo with git:

`git clone https://github.com/deepmodeling/abacus-develop.git`

`git clone https://github.com/deepmodeling/abacus-develop.git --depth=1`

`wget https://github.com/deepmodeling/abacus-develop/archive/refs/heads/develop.zip`

Get the source code of a stable version

[here](https://github.com/deepmodeling/abacus-develop/releases)[Gitee repo](https://gitee.com/deepmodeling/abacus-develop/): e.g.`git clone https://gitee.com/deepmodeling/abacus-develop.git`

. This Gitee repo is updated synchronously with GitHub.

## Update to latest release by git[#](#update-to-latest-release-by-git)

Please check the [release page](https://github.com/deepmodeling/abacus-develop/releases) for the release note of a new version.

It is OK to download the new source code from beginning following the previous step.

You can update your cloned git repo (from Github or Gitee) in-place with the following commands:

```
git remote -v
# Check if the output contains the line below
# origin https://github.com/deepmodeling/abacus-develop.git (fetch)
# The remote name is marked as "upstream" if you clone the repo from your own fork.
# Replace "origin" with "upstream" or the remote name corresponding to deepmodeling/abacus-develop if necessary
git fetch origin
git checkout v3.x.x # Replace the tag with the latest version, like v3.10.0
git describe --tags # Verify if the tag has been successfully checked out
```

[Build and Install](#build-and-install) part. If you encountered errors, try remove the `build`

directory first and reconfigure.

To use the codes under active development:

```
git checkout develop
git pull
```

## Prerequisites[#](#prerequisites)

To compile ABACUS, please make sure that the following prerequisites are present:

[CMake](https://cmake.org/)>= 3.16 .C++ compiler, supporting C++11. You can use

[Intel® C++ compiler](https://software.intel.com/enus/c-compilers)or[GCC](https://gcc.gnu.org/).

[(ref)].

MPI library. The recommended versions are

[Intel MPI](https://software.intel.com/enus/mpi-library),[MPICH](https://www.mpich.org/)or[Open MPI](https://www.open-mpi.org/).`BLAS`

,`LAPACK`

,`ScaLAPACK`

, and`ELPA`

from source file. You can use[Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)or[GFortran](https://gcc.gnu.org/fortran/).

## Install by toolchain[#](#install-by-toolchain)

[toolchain](https://github.com/deepmodeling/abacus-develop/tree/develop/toolchain) scripts to compile and install all the requirements and ABACUS itself automatically and suitable for machine characteristic in an online or offline way. The toolchain can be downloaded with ABACUS repo, and users can easily compile the requirements by running *toolchain_[gnu,intel,gcc-aocl,aocc-aocl].sh* and ABACUS itself by running *build_abacus_[gnu,intel,gcc-aocl,aocc-aocl].sh* script in the toolchain directory in `GNU`

, `Intel-oneAPI`

, `GCC-AMD AOCL`

and `AMD AOCC-AOCL`

toolchain. Sometimes, ABACUS by toolchain installation may have better efficient performance due to the suitable compiled dependencies. One should read the [README in toolchain](https://github.com/deepmodeling/abacus-develop/tree/develop/toolchain/README.md) for most of the information before use, and related tutorials can be accessed via ABACUS WeChat platform.

## Install by conda[#](#install-by-conda)

[DeepModeling conda FAQ](https://docs.deepmodeling.com/faq/conda.html) for how to setup a conda environment. A pre-built ABACUS binary with all requirements is available at [conda-forge](https://anaconda.org/conda-forge/abacus). It supports advanced features including Libxc, LibRI, and DeePKS. Conda will install the GPU-supported version of ABACUS if a valid GPU driver is present. Please refer to [the advanced installation guide](../advanced/install.html) for more details.

```
# Install
# We recommend installing ABACUS in a new environment to avoid potential conflicts:
conda create -n abacus_env abacus "libblas=*=*mkl" mpich -c conda-forge
# Run
conda activate abacus_env
OMP_NUM_THREADS=1 mpirun -n 4 abacus
# Update
conda update -n abacus_env abacus -c conda-forge
```

`"openblas=*=openmp*"`

or`"libblas=*=*mkl"`

. See[switching BLAS implementation in conda].

`OpenMPI`

and`MPICH`

variant. Install`mpich`

or`openmpi`

package to switch MPI library if required.

[conda recipe file](https://github.com/deepmodeling/abacus-develop/blob/develop/conda/meta.yaml).

Note: The

[deepmodeling conda channel]offers historical versions of ABACUS.

### Developing with conda[#](#developing-with-conda)

It is possible to build ABACUS from source based on the conda environment.

```
conda create -n abacus_env abacus -c conda-forge
conda activate abacus_env
export CMAKE_PREFIX_PATH=$CONDA_PREFIX:$CMAKE_PREFIX_PATH
# By default OpenBLAS is used; run `conda install "blas=*=mkl" mkl_fft mkl-devel -c conda-forge` to switch implementation.
export MKLROOT=$CONDA_PREFIX # If Intel MKL is required.
export CMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'`:$CMAKE_PREFIX_PATH # If DEEPKS support is required;
# usually expands to `$CONDA_PREFIX/lib/python3.1/site-packages/torch/share/cmake`
```

[Build and Install](#build-and-install) part above withou manually setting paths to dependencies. See [the advanced installation guide](../advanced/install.html) for more features. Make sure the environment variables are set before running `cmake`

. Possible command: `cmake -B build -DENABLE_MLALGO=ON -DENABLE_LIBXC=ON -DENABLE_LIBRI=ON`

.

## Install ABACUS manually[#](#install-abacus-manually)

### Install requirements[#](#install-requirements)

`apt`

and `yum`

:

```
sudo apt update && sudo apt install -y libopenblas-openmp-dev liblapack-dev libscalapack-mpi-dev libelpa-dev libfftw3-dev libcereal-dev libxc-dev g++ make cmake bc git pkgconf
```

[Intel® oneAPI toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/commercial-base-hpc.html) (former Intel® Parallel Studio) as toolchain. The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#base-kit) contains Intel® oneAPI Math Kernel Library (aka `MKL`

), including `BLAS`

, `LAPACK`

, `ScaLAPACK`

and `FFTW3`

. The [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html#hpc-kit) contains Intel® MPI Library, and C++ compiler(including MPI compiler).

`elpa`

with a different MPI library may cause conflict. Don’t forget to[set environment variables]before you start!`cmake`

will use Intel MKL if the environment variable`MKLROOT`

is set.

Please refer to our [guide](https://github.com/deepmodeling/abacus-develop/wiki/Building-and-Running-ABACUS) on installing requirements.

### Configure[#](#configure)

The basic command synopsis is:

```
cd abacus-develop
cmake -B build [-D <var>=<value>] ...
```

`CMAKE_INSTALL_PREFIX`

: the path of ABACUS binary to install;`/usr/local/bin/abacus`

by defaultCompilers

`CMAKE_CXX_COMPILER`

: C++ compiler; usually`g++`

(GNU C++ compiler) or`icpx`

(Intel C++ compiler). Can also set from environment variable`CXX`

. It is OK to use MPI compiler here.`MPI_CXX_COMPILER`

: MPI wrapper for C++ compiler; usually`mpicxx`

or`mpiicpx`

(for Intel toolkits) or`mpiicpc`

(for classic Intel Compiler Classic MPI before 2024.0).

Requirements: Unless indicated, CMake will try to find under default paths.

`MKLROOT`

: If environment variable`MKLROOT`

exists,`cmake`

will take MKL as a preference, i.e. not using`LAPACK`

,`ScaLAPACK`

and`FFTW`

. To disable MKL, unset environment variable`MKLROOT`

, or pass`-DMKLROOT=OFF`

to`cmake`

.`LAPACK_DIR`

: Path to OpenBLAS library`libopenblas.so`

(including BLAS and LAPACK)`SCALAPACK_DIR`

: Path to ScaLAPACK library`libscalapack.so`

`ELPA_DIR`

: Path to ELPA install directory; should be the folder containing ‘include’ and ‘lib’.

`ln -s elpa/include/elpa-2021.05.002/elpa elpa/include/elpa`

to help the build system find ELPA headers.`FFTW3_DIR`

: Path to FFTW3.`CEREAL_INCLUDE_DIR`

: Path to the parent folder of`cereal/cereal.hpp`

. Will download from GitHub if absent.`Libxc_DIR`

: (Optional) Path to Libxc.

`LIBRI_DIR`

: (Optional) Path to LibRI.`LIBCOMM_DIR`

: (Optional) Path to LibComm.

`ENABLE_LCAO=ON`

: Enable LCAO calculation. If SCALAPACK, ELPA or CEREAL is absent and only require plane-wave calculations, the feature of calculating LCAO basis can be turned off.`ENABLE_LIBXC=OFF`

:[Enable Libxc](../advanced/install.html#add-libxc-support)to suppport variety of functionals. If`Libxc_DIR`

is defined,`ENABLE_LIBXC`

will set to ‘ON’.`ENABLE_LIBRI=OFF`

:[Enable LibRI](../advanced/install.html#add-libri-support)to suppport variety of functionals. If`LIBRI_DIR`

and`LIBCOMM_DIR`

is defined,`ENABLE_LIBRI`

will set to ‘ON’.`USE_OPENMP=ON`

: Enable OpenMP support. Building ABACUS without OpenMP is not fully tested yet.`BUILD_TESTING=OFF`

:[Build unit tests](../advanced/install.html#build-unit-tests).`ENABLE_GOOGLEBENCH=OFF`

:[Build performance tests](../advanced/install.html#build-performance-tests)`ENABLE_MPI=ON`

: Enable MPI parallel compilation. If set to`OFF`

, a serial version of ABACUS will be compiled. It now supports both PW and LCAO.`ENABLE_COVERAGE=OFF`

: Build ABACUS executable supporting[coverage analysis](../CONTRIBUTING.html#generating-code-coverage-report). This feature has a drastic impact on performance.`ENABLE_ASAN=OFF`

: Build with Address Sanitizer. This feature would help detecting memory problems.`USE_ELPA=ON`

: Use ELPA library in LCAO calculations. If this value is set to OFF, ABACUS can be compiled without ELPA library.


Here is an example:

```
CXX=mpiicpx cmake -B build -DCMAKE_INSTALL_PREFIX=~/abacus -DELPA_DIR=~/elpa-2025.01.001/build -DCEREAL_INCLUDE_DIR=~/cereal/include
```

### Build and Install[#](#build-and-install)

After configuring, build and install by:

```
cmake --build build -j`nproc`
cmake --install build
```

`-j`

on your need: set to the number of CPU cores(`nproc`

) to reduce compilation time.

## Run ABACUS[#](#run-abacus)

### Load ABACUS[#](#load-abacus)

`CMAKE_INSTALL_PREFIX`

, please add it to your environment variable `PATH`

to locate the correct executable. `my-install-dir`

should be changed to the location of your installed abacus:`/home/your-path/abacus/bin/`

.)

```
export PATH=/my-install-dir/:$PATH
```

*abacus_env.sh*. You can source it to set the environment variables.

```
source /path/to/abacus/toolchain/abacus_env.sh
```

```
conda activate abacus_env
```

### Run with Parallelism Setting[#](#run-with-parallelism-setting)

Please set OpenMP threads by setting environment variable:

```
export OMP_NUM_THREADS=1
```

`INPUT`

file. Please make sure structure, pseudo potential, or orbital files indicated by `INPUT`

is at the correct location.

```
cd abacus-develop/examples/force/pw_Si2
```

Use 4 MPI processes to run, for example:

```
mpirun -n 4 abacus
```

`OMP_NUM_THREADS`

before running `mpirun`

:

```
OMP_NUM_THREADS=4 mpirun -n 4 abacus
```

In this case, the total thread count is 16.

`--bind-to core`

. This means that no matter how many threads you set, these threads will be restricted to run on 1 or 2 CPU cores. Therefore, setting a higher number of OpenMP threads might result in slower program execution. Hence, when using`mpirun -n`

set to 1 or 2, it is recommended to set`--bind-to none`

to avoid performance degradation. For example:`OMP_NUM_THREADS=6 mpirun --bind-to none -n 1 abacus`

. The detailed binding strategy of OpenMPI can be referred to at[https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/mpirun.1.html#quick-summary].

`OMP_NUM_THREADS`

is not set. However, it is **required** to set `OMP_NUM_THREADS`

before running `mpirun`

to avoid potential performance issues.

Please refer to [hands-on guide](hands_on.html) for more instructions.

`lscpu | grep 'per core'`

and see if ‘Thread(s) per core’ is 2.

## Container Deployment[#](#container-deployment)

[here](https://github.com/deepmodeling/abacus-develop/pkgs/container/abacus). For a quick start: pull the image, prepare the data, run container. Instructions on using the image can be accessed in [Dockerfile](#../../Dockerfile). A mirror is available by `docker pull registry.dp.tech/deepmodeling/abacus`

.

[Package Page](https://github.com/orgs/deepmodeling/packages?repo_name=abacus-develop).

[Developing inside a Container](https://code.visualstudio.com/docs/remote/containers#_quick-start-try-a-development-container). Choose `Open a Remote Window -> Clone a Repository in Container Volume`

in VS Code command palette, and put the [git address](https://github.com/deepmodeling/abacus-develop.git) of `ABACUS`

when prompted.

For online development environment, we support [GitHub Codespaces](https://github.com/codespaces): [Create a new Codespace](https://github.com/codespaces/new?machine=basicLinux32gb&repo=334825694&ref=develop&devcontainer_path=.devcontainer%2Fdevcontainer.json&location=SouthEastAsia)

We also support [Gitpod](https://www.gitpod.io/): [Open in Gitpod](https://gitpod.io/#https://github.com/deepmodeling/abacus-develop)

## Command line options[#](#command-line-options)

`abacus --version`

, the result will be like:

```
ABACUS version v3.9.0.2
```

`INPUT`

file by running the command `abacus --check-input`

, the result will be like:

```
ABACUS v3.9.0.2
Atomic-orbital Based Ab-initio Computation at UStc
Website: http://abacus.ustc.edu.cn/
Documentation: https://abacus.deepmodeling.com/
Repository: https://github.com/abacusmodeling/abacus-develop
https://github.com/deepmodeling/abacus-develop
Commit: unknown
Tue Jun 18 14:20:31 2024
MAKE THE DIR : OUT.ABACUS/
----------------------------------------------------------
INPUT parameters have been successfully checked!
----------------------------------------------------------
```

Warnings will be given if there are any errors in the `INPUT`

file.
