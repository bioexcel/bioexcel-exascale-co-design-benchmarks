# Background

Between 50% and 80% of the computational time in molecular dynamics
simulations with GROMACS is spent in the non-bonded pair interaction
kernels. Thus these kernels are the prime target for code optimization
and co-design. The two co-design benchmarks in this suite both provide
functionality to benchmark the non-bonded kernels only. They both use
the same kernel code that is used in the full GROMACS package. These
kernels rely on explicit SIMD parallelization to maxmimize performance.
The nonbonded_benchmark tool also supports OpenMP parallelization,
whereas the minimal kernel_microbenchmark tool does not.

The quantity of interest in the number of iterations per second as a
function of the number of atoms per core or thread. This can be
normalized by the number of atoms to get the efficiency of the kernel
as a function of the number of atoms. OpenMP parallization can often
improve the performance as the kernels are usually instruction latency
limited.

The relevant performance regime is that which occurs in most
biomolecular MD simulations. This is from around 400 atoms per core
or thread up to a few thousand. At lower atom per core counts
the performance of GROMACS will be limited by communication overhead
(not present in these tools). High atom counts are less relevant
because one would usally use more hardware for a larger problem to
reduce the time to solution. Of particular interest is the performance
behavior at cache size boundaries. In the low atom count per core
regime everything will usually fit in level-2 cache. For larger counts
level-3 cache can come into play which can lower the throughput
of the kernels.

# Contents

Two different benchmark cases are contained in respective independent
subdirectories. Please see the respective subdirectories for more
details

## nonbonded_benchmark

nonbonded_benchmark provides a tool for observing the performance of
production versions of key GROMACS kernels on real CPU hardware.

## kernel_microbenchmark

kernel_microbenchmark provides a much simpler set of code that runs a
single version of one such kernel. It does so in a form suitable for
running a single iteration in a CPU processor simulator without the
prohibitive cost of running the full application in the simulator.
