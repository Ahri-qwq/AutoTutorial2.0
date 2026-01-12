# https://abacus.deepmodeling.com/en/latest/advanced/md.html

> Source: https://abacus.deepmodeling.com/en/latest/advanced/md.html

# Molecular Dynamics[#](#molecular-dynamics)

[calculation](input_files/input-main.html#calculation) to be `md`

, ABACUS currently provides several different MD evolution methods, which is specified by keyword [md_type](input_files/input-main.html#md-type) in the `INPUT`

file:

fire: a MD-based relaxation algorithm, see

[details](#fire)herenve: NVE ensemble with velocity Verlet algorithm

nvt: NVT ensemble

npt: Nose-Hoover style NPT ensemble

langevin: NVT ensemble with Langevin thermostat

msst: MSST method


[md_type](input_files/input-main.html#md-type) is set to nvt, [md_thermostat](input_files/input-main.html#md-thermostat) is used to specify the temperature control method used in NVT ensemble.

nhc: Nose-Hoover chain

anderson: Anderson thermostat

berendsen: Berendsen thermostat

rescaling: velocity Rescaling method 1

rescale_v: velocity Rescaling method 2


[md_type](input_files/input-main.html#md-type) is set to npt, [md_pmode](input_files/input-main.html#md-pmode) is used to specify the cell fluctuation mode in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.

iso: isotropic cell fluctuations

aniso: anisotropic cell fluctuations

tri: non-orthogonal (triclinic) simulation box


[list of keywords](input_files/input-main.html#molecular-dynamics) to control relevant parmeters used in MD simulations.

`MD_dump`

， in which the atomic forces, atomic velocities, and lattice virial are controlled by keyword [dump_force](input_files/input-main.html#dump-force), [dump_vel](input_files/input-main.html#dump-vel), and [dump_virial](input_files/input-main.html#dump-virial), respectively.

[Examples](#../../examples/md/lcao_gammaonly_Si8/) of MD simulations are also provided. There are eight INPUT files corresponding to eight different MD evolution methods in the directory. For examlpe, `INPUT_0`

shows how to employ the NVE simulation.

`INPUT`

, and run ABACUS.

## FIRE[#](#fire)

[FIRE](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201) (fast inertial relaxation engine) is a MD-based minimization algorithm. It is based on conventional molecular dynamics with additional velocity modifications and adaptive time steps. The MD trajectory will descend to an energy-minimum.

## NVE[#](#nve)

Currently NVE ensemble in ABACUS is implemented based on the [velocity verlet algorithm](https://aip.scitation.org/doi/abs/10.1063/1.442716).

## Nose Hoover Chain[#](#nose-hoover-chain)

[Nose-Hoover style non-Hamiltonian equations of motion](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.31.1695) which are designed to generate positions and velocities sampled from NVT and NPT ensemble.

## Langevin[#](#langevin)

[Langevin thermostat](https://en.wikipedia.org/wiki/Langevin_dynamics) can be used for molecular dynamics equations by assuming that the atoms being simulated are embedded in a sea of much smaller fictional particles. In many instances of solute-solvent systems, the behavior of the solute is desired, and the behavior of the solvent is non-interesting(e.g. proteins, DNA, nanoparticles in solution). In these cases, the solvent influences the dynamics of the solute(typically nanoparticles) via random collisions, and by imposing a frictional drag force on the motion of the nanoparticle in the solvent. The damping factor and the random force combine to give the correct NVT ensemble.

## Anderson[#](#anderson)

[Anderson thermostat](https://aip.scitation.org/doi/abs/10.1063/1.439486) couples the system to a heat bath that imposes the desired temperature to simulate the NVT ensemble. The coupling to a heat bath is represented by stochastic collision that act occasionally on randomly selected particles.

## Berendsen[#](#berendsen)

[Berendsen thermostat](https://aip.scitation.org/doi/10.1063/1.448118), which rescales their velocities every timestep. In this scheme, the system is weakly coupled to a heat bath with some temperature. Though the thermostat does not generate a correct canonical ensemble (especially for small systems), for large systems on the order of hundreds or thousands of atoms/molecules, the approximation yields roughly correct results for most calculated properties.

## Rescaling[#](#rescaling)

[md_tolerance](input_files/input-main.html#md-tolerance) (Kelvin).

## Rescale_v[#](#rescale-v)

[md_nraise](input_files/input-main.html#md-nraise) steps the current temperature is rescaled to target temperature.

## MSST[#](#msst)

[Multi-Scale Shock Technique (MSST) integration](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.235503) to update positions and velocities each timestep to mimic a compressive shock wave passing over the system. The MSST varies the cell volume and temperature in such a way as to restrain the system to the shock Hugoniot and the Rayleigh line. These restraints correspond to the macroscopic conservation laws dictated by a shock front.

## DPMD[#](#dpmd)

Compiling ABACUS with [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit), MD calculations based on machine learning DP model is enabled.

[esolver_type](input_files/input-main.html#esolver-type) should be set to `dp`

. And the filename of DP model is specified by keyword [pot_file](input_files/input-main.html#pot-file).

## NEP[#](#nep)

[NEP](https://gpumd.org/potentials/nep.html)), MD simulations using NEP models are enabled. To use this feature, set [esolver_type](input_files/input-main.html#esolver-type) to `nep`

and specify the potential file path with the [pot_file](input_files/input-main.html#pot-file) keyword in your INPUT file.
