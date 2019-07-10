# CP2K whole-application biomolecular QM/MM benchmark

## Description

The benchmark system consists of a phytochrome dimer (CBD-PHY) in a 12x12x12 nm water box, to which a Biliverdin chromophore is bound. The phytochrome and waters (167922 atoms all together) are treated using MM. The Amber03 forcefield is specified, however classical MM-MM interactions (both bonded and non-bonded) are disabled for this calculation since given the planned usage of CP2K for QM/MM in BioExcel through an interface driven by GROMACS (which computes MM interactions and performs motion) we are only concerned with CP2K’s performance computing QM and QM/MM electrostatic interactions. The chromophore (68 atoms) is located within a designated QM box of dimensions 2.5x2.5x2.5 nm. QM energies and forces are computed using the DZVP basis set (655 Gaussians) and PBE DFT functional. 

A nearly-SCF-converged wave function solution for the QM region, force-opt-qmmm-RESTART.wfn, is provided as a “RESTART” input, thereby allowing the benchmark to proceed almost directly to a one-time computation of the QM/MM forces and energies, which is done with periodic boundary conditions. After this, the run finishes without any time integration or atomic motion. 


## How to run

The whole-application benchmark described in the parent directory's README can be run directly with a regular production build of CP2K, e.g. as follows:

mpirun -n 48 cp2k.popt -in qmmm-force-opt.in > qmmm-force-opt.out

