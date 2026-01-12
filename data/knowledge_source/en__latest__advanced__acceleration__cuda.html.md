# https://abacus.deepmodeling.com/en/latest/advanced/acceleration/cuda.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/acceleration/cuda.html

# CUDA GPU Implementations[#](#cuda-gpu-implementations)

**Full gpu implementations**: During the SCF progress,`Psi`

,`Hamilt`

,`Hsolver`

,`DiagCG`

, and`DiagoDavid`

classes are stored or calculated by the GPU devices.**Electronic state data**: (e.g. electronic density) are moved from the GPU to the CPU(s) every scf step.**Accelerated by the NVIDIA libraries**:`cuBLAS`

for common linear algebra calculations,`cuSolver`

for eigen values/vectors, and`cuFFT`

for the conversions between the real and recip spaces.**Multi GPU supprted**: Using multiple MPI tasks will often give the best performance. Note each MPI task will be bind to a GPU device with automatically computing load balancing.**Parallel strategy**: K point parallel.

## Required hardware/software[#](#required-hardware-software)

Check if you have an NVIDIA GPU: cat /proc/driver/nvidia/gpus/*/information

Install a driver and toolkit appropriate for your system (SDK is not necessary)


## Building ABACUS with the GPU support:[#](#building-abacus-with-the-gpu-support)

Check the [Advanced Installation Options](https://abacus-rtd.readthedocs.io/en/latest/advanced/install.html#build-with-cuda-support) for the installation of CUDA version support.

[enable GPU support](https://github.com/marekandreas/elpa/blob/master/documentation/INSTALL.md).

## Run with the GPU support by editing the INPUT script:[#](#run-with-the-gpu-support-by-editing-the-input-script)

`INPUT`

file we need to set the input parameter [device](../input_files/input-main.html#device) to `gpu`

. If this parameter is not set, ABACUS will try to determine if there are available GPUs.

`ks_solver`

: For the PW basis, CG, BPCG and Davidson methods are supported on GPU; set the input parameter[ks_solver](../input_files/input-main.html#ks-solver)to`cg`

,`bpcg`

or`dav`

. For the LCAO basis,`cusolver`

,`cusolvermp`

and`elpa`

is supported on GPU.**multi-card**: ABACUS allows for multi-GPU acceleration. If you have multiple GPU cards, you can run ABACUS with several MPI processes, and each process will utilize one GPU card. For example, the command`mpirun -n 2 abacus`

will by default launch two GPUs for computation. If you only have one card, this command will only start one GPU.

## Examples[#](#examples)

We provides [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/gpu) of gpu calculations.

## Known limitations[#](#known-limitations)

PW basis:

`kpar`

will be set to match the number of MPI tasks automatically.or the environmental variable`CMAKE_CUDA_ARCHITECTURES`

.`CUDAARCHS`


LCAO basis:

When using elpa on GPUs, some ELPA internal logs will be output.
