# CP2K QM/MM kernel benchmark

## Description

The code chosen for the kernel consists of the CP2K subroutine qmmm_forces_with_gaussian_LG, which was found to be by far the most costly single subroutine in the execution of the CP2K whole-application biomolecular QM/MM benchmark, taking almost 25% of overall runtime. The subroutine and hence the kernel computes the force contributions due to the long-range electrostatic interaction between QM and MM regions. It does this by treating the QM-MM coupling using the method of multigrids and Gaussian Expansion of the QM/MM Electrostatic Potential (GEEP) described in [2] and taking into account periodic boundary conditions as described in [3].  Computation of the long-range electrostatic QM/MM energies (subroutine qmmm_elec_with_gaussian_LG) were found to contribute another 10% of runtime running the whole-application benchmark, and since the relevant subroutine is similar in form to the first part of that computing the forces, any optimisations found to speed up the kernel benchmark are likely to also yield improvement in the QM/MM energies subroutine.
 
As well as being motivated by performance profiling results for the whole-application benchmark, the choice of benchmark kernel is also justified with reference to an established recognition in the CP2K QM/MM literature [2] that generally speaking evaluating the electrostatic interaction between QM and MM parts can be as time consuming as solving the QM part itself using DFT methods. As the latter is already the target of wider attention outside of BioExcel, this provides an additional reason for focusing here on computation of the QM/MM potential as a target for performance optimisation and co-design. 

### References
[1] https://www.cp2k.org/docs
[2] An Efficient Real Space Multigrid QM/MM Electrostatic Coupling. Laino, T.; Mohamed, F.; Laio, A.; Parrinello, M.  J. Chem. Theory Comput. 2005, 1, 1176-1184. DOI: 10.1021/ct050123f
[3] An Efficient Linear-Scaling Electrostatic Coupling for Treating Periodic Boundary Conditions in QM/MM Simulations. Laino, T.; Mohamed, F. Laio, A. and Parrinello, M. J. Chem. Theory Comput. 2006 2 (5), 1370-1378 DOI: 10.1021/ct6001169
[4] https://www.cp2k.org/performance
[5] https://repository.prace-ri.eu/git/UEABS/ueabs
[6] https://www.cp2k.org/dev:profiling

## How to build

## How to run

