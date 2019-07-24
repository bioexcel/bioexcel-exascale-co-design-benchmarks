# CP2K whole-application biomolecular QM/MM benchmark

## Description

The benchmark system consists of a phytochrome dimer (CBD-PHY) in a 12x12x12 nm water box, to which a Biliverdin chromophore is bound. The phytochrome and waters (167922 atoms all together) are treated using MM. The Amber03 forcefield is specified, however classical MM-MM interactions (both bonded and non-bonded) are disabled for this calculation since given the planned usage of CP2K for QM/MM in BioExcel through an interface driven by GROMACS (which computes MM interactions and performs motion) we are only concerned with CP2K’s performance computing QM and QM/MM electrostatic interactions. The chromophore (68 atoms) is located within a designated QM box of dimensions 2.5x2.5x2.5 nm. QM energies and forces are computed using the DZVP basis set (655 Gaussians) and PBE DFT functional. 

A nearly-SCF-converged wave function solution for the QM region, force-opt-qmmm-RESTART.wfn, is provided as a “RESTART” input, thereby allowing the benchmark to proceed almost directly to a one-time computation of the QM/MM forces and energies, which is done with periodic boundary conditions. After this, the run finishes without any time integration or atomic motion. 

## Requirements

The whole-application benchmark can be run directly with a regular production build of CP2K, in serial or in parallel (see sections below on how to run and expected runtimes). It has been tested with [CP2K version 6.1](https://github.com/cp2k/cp2k/releases/tag/v6.1.0), see details on <https://www.cp2k.org/download> and <https://www.cp2k.org/howto:compile>. Memory requirements are modest, corresponding to on average around 0.6GB per MPI rank and less than 1GB for any rank. 


## How to run

The benchmark can be run directly with a regular production build of CP2K as follows:

`mpirun -n 24 cp2k.popt -in force-opt-qmmm.in`

Or, to redirect output from stdout to file for future reference:

`mpirun -n 24 cp2k.popt -in force-opt-qmmm.in > force-opt-qmmm.log`

The number of MPI ranks (24 in the above example) should be adjusted to be equal to the number of available cores.

## Expected runtimes and reference outputs

Running on 24 cores on a CRAY XC30 compute node (two 12-core Intel E5-2697v2@2.7GHz processors) using a production build of CP2K, the benchmark takes around 200 seconds to complete and uses less than 1GB of memory per rank (memory high watermark) with a total requirement of at most 16GB at any given time during the simulation. Runtimes and maximum total memory requirements on a range of numbers of cores on the same platform, including across two compute nodes for 48 cores, are given in the table below:

| cores | runtime (s) | memory (GB) |
| ----- | ----------- | ----------- |
| 48    | 145         | 26	    |
| 24    | 203         | 16	    |
| 12    | 322         | 7	    |
| 6     | 539         | 2.7	    |
| 2     | 1392        | 1.3	    |
| 1     | 2657        | 0.6	    |

Running the hybrid MPI+OpenMP-threaded version of CP2K (cp2k.psmp executable) was not found to yield any performance advantage for this benchmark over running pure MPI. Reference output log files are included in the `outputs` directory. 

