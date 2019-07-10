# bioexcel-exascale-co-design-benchmarks

Exascale co-design benchmarks and kernels produced by BioExcel, an EC
HPC Centre of Excellence.

The contents relates to the BioExcel core applications: GROMACS,
HADDOCK, and CP2K. It is intended to provide a simple way for people
not expert in those software packages to obtain relevant feedback.
Such feedback could be relevant for considering changes to hardware,
drivers, libraries, compilers, or system configurations would be for
users of these applications.

These are

* a version of GROMACS containing a benchmark tool of production
  versions of key CPU pair-interaction kernels
* one such GROMACS pair-interaction kernel in microbenchmark form
  suitable for running in a CPU simulator
* a script to illustrate HADDOCK I/O usage ... TODO UU
* a CP2K QM/MM kernel benchmark (code + input data)
* a set of CP2K whole-application inputs on which the kernel benchmark was based

Further details are found in the respective subfolders.

Feedback and suggestions for further collaboration are warmly invited
and are best directed to the maintainers named in the respective
subdirectories.

The content is made available under the MIT license, except where the
bundled source code declares another free-and-open-source license. It is
open and may be used freely and for any purpose.
