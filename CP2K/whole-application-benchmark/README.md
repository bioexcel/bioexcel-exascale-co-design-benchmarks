# CP2K whole-application biomolecular QM/MM benchmark

## Description

This directory provides the CP2K whole-application biomolecular QM/MM benchmark that formed the basis for the construction of the co-design kernel benchmark also provided in this repository (under `../kernel-benchmark`).  The benchmark was contributed by BioExcel researchers at the University of Jyvaskula and relates to the project's fluorescent protein use case [1]. It calculates the QM/MM forces and energies for a phytochrome dimer (CBD-PHY) in a 12x12x12 nm water box to which a Biliverdin chromophore is bound.

The chromophore, consisting of 68 atoms, is located within a designated QM box of dimensions 2.5x2.5x2.5 nm. Properties of the QM region are computed using the Gaussian and Plane Waves (GPW) method, the DZVP basis set (655 Gaussians), and the PBE correlation-exchange functional. The QM/MM coupling is described using the Gaussian expansion of the Electrostatic Potential (GEEP) method described in [2][3]. Bonds that cross the QM/MM boundary are treated using the Generalized Hybrid Orbital link method.

The phytochrome and waters (167,922 atoms all together) are treated using MM. The Amber03 forcefield is specified, however classical MM-MM interactions (both bonded and non-bonded) are disabled for this calculation since it is planned to call CP2K from GROMACS via an interface, hence we are primarily concerned with CP2K's performance computing QM and QM/MM electrostatic interactions (forces and energies). These will ultimately be used by GROMACS in combination with its internal MM interaction computation to perform time integration of atomic motion. 

A nearly-SCF-converged wave function solution for the QM region is provided as a "RESTART" input, thereby allowing the benchmark to proceed almost directly to a one-time computation of the QM/MM forces and energies, which is done with periodic boundary conditions using the method described in [2][3]. After this, the run finishes without any time integration or atomic motion.

Compared to currently common simulations of photoactive systems this benchmark is fairly large - three to four times the size of usual fluorescent proteins.


### 
[1] https://bioexcel.eu/research/projects/electronic-interaction-phenomena-proton-dynamics-and-fluorescent-proteins/
[2] An Efficient Real Space Multigrid QM/MM Electrostatic Coupling. Laino, T.; Mohamed, F.; Laio, A.; Parrinello, M.  J. Chem. Theory Comput. 2005, 1, 1176-1184. https://doi.org/10.1021/ct050123f  
[3] An Efficient Linear-Scaling Electrostatic Coupling for Treating Periodic Boundary Conditions in QM/MM Simulations. Laino, T.; Mohamed, F. Laio, A. and Parrinello, M. J. Chem. Theory Comput. 2006 2 (5), 1370-1378 https://doi.org/10.1021/ct6001169  


## Requirements

The whole-application benchmark can be run directly with a regular production build of CP2K, in serial or in parallel (see sections below on how to run and expected runtimes). It has been tested with [CP2K version 6.1](https://github.com/cp2k/cp2k/releases/tag/v6.1.0), see details on <https://www.cp2k.org/download> and <https://www.cp2k.org/howto:compile>. Memory requirements are modest, corresponding to on average around 0.6GB per MPI rank and less than 1GB for any rank. 


## How to run

The benchmark can be run directly with a regular production build of CP2K as follows:

`mpirun -n 24 cp2k.popt -in force-opt-qmmm.in`

The number of MPI ranks (24 in this example) should be adjusted to equal the number of available cores.

Or, to redirect output from stdout to file for future reference:

`mpirun -n 24 cp2k.popt -in force-opt-qmmm.in > force-opt-qmmm.log`

Alternatively, if using the hybrid MPI+OpenMP-threaded version of CP2K (`cp2k.psmp` executable), then equivalently to the above do:

```
export OMP_NUM_THREADS=1
mpirun -n 24 cp2k.psmp -in force-opt-qmmm.in > force-opt-qmmm.log
```

Note: OpenMP-threading was not found to yield any performance advantage for this benchmark over running the pure MPI executable (`cp2k.popt`).

## Expected runtimes and reference outputs

Running on 24 cores (two 12-core Intel E5-2697v2@2.7GHz processors) on a single compute node of ARCHER (a CRAY XC30 machine) using a production build of CP2K, the benchmark takes around 200 seconds to complete and uses less than 1GB of memory per rank (memory high watermark) with a total requirement of at most 16GB at any given time during the simulation. Runtimes and maximum total memory requirements on a range of numbers of cores on the same platform, including across two compute nodes for 48 cores, are given in the table below:

| cores | runtime (s) | speedup     | parallel efficiency (%) | total memory required (GB) |
| ----- | ----------- | ----------- | ----------------------- | -------------------------- |
| 1     | 2657        | 1	    | 100		      | 0.6	    		   |
| 2     | 1392        | 1.9	    | 95		      | 1.3	    		   |
| 6     | 539         | 4.9	    | 82		      | 2.7	    		   |
| 12    | 322         | 8.3	    | 69		      | 7	    		   |
| 24    | 203         | 13.0	    | 54	    	      | 16			   |
| 48    | 145         | 18.3	    | 38		      | 26			   |







Reference output log files are included in the `outputs` directory. 

