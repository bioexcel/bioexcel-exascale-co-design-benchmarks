GROMACS kernel_microbenchmark tool
==================================

When developing new architectures, CPU designs are first tested in
simulators. It might be too much effort to build or run the complete,
or a significant part of, the GROMACS code base in an
emulator. Therefore, we have extracted the kernels themselves as well
as their inputs into this stand-alone app. The app is built from the
minimal GROMACS code needed to run the benchmarks. This enables
quicker and simpler co-design of GROMACS for future architectures.

Only a single system size is supported and only the most important
kind of kernel is currently supported. The fixed, same system size is
not a limitation when running on an emulator, because on an emulator
one does not need to run with large numbers of threads as would be
needed on a large multicore processor to assess thermally limited
performance.

The current version of this kernel computes pair interaction distances
and looks up the pairwise N-body interaction strength from
pre-computed tables. In future, a version that supports computing the
interaction directly will be added, so that CPU designers can see the
impact of design choices on more or less CPU-bound forms of the same
workload.

In future we will consider adding SIMD-accelerated and/or
accelerator-based version of the kernel benchmark app. Although the
same limitations for accelerator benchmarking discussed for the
nonbonded_benchmark tool hold here, the situation for emulation is
slightly different. For a completely new accelerator architecture it
might be needed to write the actual kernel from scratch. What is
needed then is a framework which hands prepared input data to the
kernel. This is exactly what the benchmark app provides. So, this can
already be used when the data and pairlist is laid out for GROMACS CPU
SIMD use. With little effort the app can be extended to provide the
data and pairlist in GROMACS GPU layout. But since several details of
the optimal layout will depend on the accelerator architecture and
capabilities, it is not useful and in many cases impossible to
pre-generate all data. The setups of interest can be generated with
little effort by the GROMACS developers in a co-design effort.

Installation
^^^^^^^^^^^^

The subdirectory contains a simple standalone CMake project that
requires a standard-compliant C++14 compiler. To build and run, use

    cmake . -DCMAKE_CXX_COMPILER=/path/to/your/cc
    make
    ./gromacs_kernel_microbenchmark

The tool checks for correct results after the kernel completes,
reporting success or failure to the terminal and a normal UNIX error
code (0 for success).

Maintainers
^^^^^^^^^^^

Mark Abraham <mjab@kth.se>

