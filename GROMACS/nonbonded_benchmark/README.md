# GROMACS nonbonded_benchmark tool

This benchmark tool is capably of running a reproducible CPU-only
N-body pair-interaction kernels. These kernels are those used in the
production version of GROMACS. They support multi-threading with
OpenMP and any SIMD architecture supported in the GROMACS SIMD
module. A range of representative input sizes are available. This tool
is part of the GROMACS distribution and will be released with the 2020
and future releases. It is also available in the code distribution and
repository associated with this deliverable.

The description of the function, capabilities and the options of the
tool are provided through the help function on the tool (e.g. gmx
nonbonded_bench -h). These are copied below for convenience.

## Installation

Build exactly as you would build a production version of GROMACS.  See
the INSTALL file in the gromacs_source_code subdirectory for details.

This code is based on a development version of GROMACS between 2019
and 2020 (commit 219d24220e5e552 from gerrit.gromacs.org), so there is
no official regressiontest suite for it. If you want to run the tests
with "make check", please direct CMake to use the version of the
regressiontest suite contained in this repository, e.g.  with

```bash
cd gromacs_source_code
mkdir build
cd build
cmake .. -DREGRESSIONTEST_PATH=../../gromacs_regressiontests
make check
bin/gmx nonbonded_bench
```

## Maintainers

Berk Hess <berk at kth.se>
Mark Abraham <mjab at kth.se>

## Detailed description

```
SYNOPSIS

gmx nonbonded_bench [-size <int>] [-nt <int>] [-simd <enum>]
             [-coulomb <enum>] [-[no]table] [-combrule <enum>] [-[no]halflj]
             [-[no]energy] [-[no]all] [-cutoff <real>] [-iter <int>]
             [-warmup <int>] [-[no]cycles]

DESCRIPTION

gmx nonbonded_bench runs benchmarks for one or more so-called Nbnxm non-bonded
pair kernels. The non-bonded pair kernels are the most compute intensive part
of MD simulations and usually comprise 60 to 90 percent of the runtime. For
this reason they are highly optimized and several different setups are
available to compute the same physical interactions. In addition, there are
different physical treatments of Coulomb interactions and optimizations for
atoms without Lennard-Jones interactions. There are also different physical
treatments of Lennard-Jones interactions, but only a plain cut-off is
supported in this tool, as that is by far the most common treatment. And
finally, while force output is always necessary, energy output is only
required at certain steps. In total there are 12 relevant combinations of
options. The combinations double to 24 when two different SIMD setups are
supported. These combinations can be run with a single invocation using the
-all option. The behavior of each kernel is affected by caching behavior,
which is determined by the hardware used together with the system size and the
cut-off radius. The larger the number of atoms per thread, the more L1 cache
is needed to avoid L1 cache misses. The cut-off radius mainly affects the data
reuse: a larger cut-off results in more data reuse and makes the kernel less
sensitive to cache misses.

OpenMP parallelization is used to utilize multiple hardware threads within a
compute node. In this benchmarks there is no interaction between threads,
apart from starting and closing a single OpenMP parallel region per iteration.
Additionally thread interact through sharing and evicting data from shared
cashes. The number of threads to use is set with the -nt option. Thread
affinity is important, especially with SMT and shared cashes. Affinities can
be set through the OpenMP library using the GOMP_CPU_AFFINITY environment
variable.

The benchmark tool times one or more kernels by running them repeatedly for a
number of iterations set by the -iter option. An initial kernel call is done
to avoid additional initial cache misses. Times are recording in cycles read
from efficient, high accuracy counters in the CPU. Note that these often do
not correspond to actual clock cycles. For each kernel the total number of
cycles is reported, the cycles per iteration, as well as the total and useful
number of pair interactions per cycle. Because a cluster pair list is ued
instead atom pair list, interactions are also computed for some atom pairs
that are beyond the cut-off distance. These pairs are not useful (except for
additional buffering, but that is not of interest here), only a side effect of
the cluster-pair setup. The SIMD 2xMM kernel has a higher useful pair ratio
then the 4xM kernel due to a smaller cluster size, but a lower total pair
throughput. It is best to run this, or for that matter any, benchmark with
locked CPU clocks, as thermal throttling can significantly affect performance.
If that is not an option, the -warmup option can be used to run initial,
untimed iterations to warm up the processor.

The most relevant regime is between 0.1 to 1 millisecond per iteration. Thus
it is useful to run with system sizes that cover both ends of this regime.

The -simd and -table options select different implementations to compute the
same physics. The choice of these options should ideally be optimized for the
target hardware. Up till now, we only found tabulated Ewald correction to be
useful on 2-wide SIMD or 4-wide SIMD without FMA support. As all modern
architectures are wider and support FMA, we do not use tables by default. The
only exception are kernels without SIMD, which only support tables.Options
-coulomb, -combrule and -halflj depend on the force field and composition of
the simulated system. The optimization of computing Lennard-Jones interactions
for only half of the atoms in a cluster is useful for water, which does not
use Lennard-Jones on hydrogen atoms in most water models. In the MD engine,
any clusters where at most half of the atoms have LJ interactions will
automatically use this kernel. And finally, the -energy option selects the
computation of energies, which are usually only needed infrequently.

OPTIONS

Other options:

 -size   <int             (1)
           The system size is 3000 atoms times this value
 -nt     <int             (1)
           The number of OpenMP threads to use
 -simd   <enum            (auto)
           SIMD type, auto runs all supported SIMD setups or no SIMD when SIMD
           is not supported: auto, no, 4xm, 2xmm
 -coulomb <enum           (ewald)
           The functional form for the Coulomb interactions: ewald,
           reaction-field
 -[no]table                 (no)
           Use lookup table for Ewald correction instead of analytical
 -combrule <enum          (geometric)
           The LJ combination rule: geometric, lb, none
 -[no]halflj                (no)
           Use optimization for LJ on half of the atoms
 -[no]energy                (no)
           Compute energies in addition to forces
 -[no]all                   (no)
           Run all 12 combinations of options for coulomb, halflj, combrule
 -cutoff <real            (1)
           Pair-list and interaction cut-off distance
 -iter   <int             (100)
           The number of iterations for each kernel
 -warmup <int             (0)
           The number of iterations for initial warmup
 -[no]cycles                (no)
           Report cycles/pair instead of pairs/cycle
```
