# CP2K QM/MM kernel benchmark

## Description

The code chosen for the kernel consists of the CP2K subroutine qmmm_forces_with_gaussian_LG, which was found to be by far the most costly single subroutine in the execution of the CP2K whole-application biomolecular QM/MM benchmark, taking almost 25% of overall runtime. The subroutine and hence the kernel computes the force contributions due to the long-range electrostatic interaction between QM and MM regions. It does this by treating the QM-MM coupling using the method of multigrids and Gaussian Expansion of the QM/MM Electrostatic Potential (GEEP) described in [2] and taking into account periodic boundary conditions as described in [3].  Computation of the long-range electrostatic QM/MM energies (subroutine qmmm_elec_with_gaussian_LG) were found to contribute another 10% of runtime running the whole-application benchmark, and since the relevant subroutine is similar in form to the first part of that computing the forces, any optimisations found to speed up the kernel benchmark are likely to also yield improvement in the QM/MM energies subroutine.
 
As well as being motivated by performance profiling results for the whole-application benchmark, the choice of benchmark kernel is also justified with reference to an established recognition in the CP2K QM/MM literature [2] that generally speaking evaluating the electrostatic interaction between QM and MM parts can be as time consuming as solving the QM part itself using DFT methods. As the latter is already the target of wider attention outside of BioExcel, this provides an additional reason for focusing here on computation of the QM/MM potential as a target for performance optimisation and co-design. 

### References
[1] https://www.cp2k.org/docs  
[2] An Efficient Real Space Multigrid QM/MM Electrostatic Coupling. Laino, T.; Mohamed, F.; Laio, A.; Parrinello, M.  J. Chem. Theory Comput. 2005, 1, 1176-1184. https://doi.org/10.1021/ct050123f  
[3] An Efficient Linear-Scaling Electrostatic Coupling for Treating Periodic Boundary Conditions in QM/MM Simulations. Laino, T.; Mohamed, F. Laio, A. and Parrinello, M. J. Chem. Theory Comput. 2006 2 (5), 1370-1378 https://doi.org/10.1021/ct6001169  
[4] https://www.cp2k.org/performance  
[5] https://repository.prace-ri.eu/git/UEABS/ueabs  
[6] https://www.cp2k.org/dev:profiling

## Requirements

The benchmark kernel is entirely self contained - the kernel code itself (`qmmm_gpw_forces.f90`), a  wrapper (`kernel_benchmark.f90`), some dependencies adapted from the CP2K codebase, and input files corresponding to the whole-application benchmark also contained in this repository are all included (more information in `./src/README.md` and `./data/README.md`). All that is required to build the kernel benchmark is a Fortran 2008-compatible compiler (gfortran is recommended, as it is for CP2K as a whole). If a Fortran 2008-compatible compiler is not available, you will have to modify a few file IO statements in `./src/kernel_benchmark.f90` to not make use of the 'newunit' feature. The code runs serially so no MPI library is required. The kernel benchmark has negligible memory requirements (~20MB memory high watermark). 

## How to build

- In the Makefile in this directory, specify a Fortran compiler that supports Fortran 2008 (in particular the `newunit` IO feature), and specify an optimisation level by setting FCFLAGS (e.g. to `-O3`)

- Type `make`

- The resulting executable kernel_benchmark will now appear in this directory

- `make clean` can be executed in this directory to clear all object and module files as well as the `kernel_benchmark` executable

## How to run

The executable `kernel_benchmark` should be run without any input parameters. It expects the `data` subdirectory included in this repository to be located in the working directory where it executes. The executable prints out times obtained using calls to `system_clock()` to time the initialisation and the time spent running the kernel. 

## Checking results and expected runtime

A reference output file `./outputs/Forces.out` is provided against which any `Forces.out` produced by running the benchmark can be compared using the `diff` command. The reference output was produced by compiling with gfortran version 8.2.0 with `-g` and default optimisation level (`-O0`) and running the benchmark on a single core of a quad-core Intel i7-3820QM@2.7GHz on a macOS laptop, for which the kernel call takes just over 3 minutes. No difference in `Forces.out` was encountered when running on the same platform with `-O3`, for which the kernel call takes around 60 seconds. Running on a single core of an Intel E5-2697v2@2.7GHz processor on an HPC compute node obtains similar timings, taking 2 minutes 34 seconds with `-O0` and 53 seconds with `-O3`.

## Adapting the kernel benchmark

One might like to run the kernel benchmark code with different inputs. This is facilitated by inclusion of the file `./tools/qmmm_gpw_forces.F`, which was used to write to file the data structures needed to run the kernel benchmark and which are included here under `./data`. For more information see `./tools/README.md`