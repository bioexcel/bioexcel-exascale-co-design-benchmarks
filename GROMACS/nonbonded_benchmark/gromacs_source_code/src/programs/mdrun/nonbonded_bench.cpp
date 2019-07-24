/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief This file contains the main function for the non-bonded kernel benchmark
 *
 * \author Berk Hess <hess@kth.se>
 */

#include "gmxpre.h"

#include "nonbonded_bench.h"

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/nbnxm/benchmark/bench_setup.h"
#include "gromacs/utility/arraysize.h"

namespace gmx
{

//! Implements C-style main function for mdrun
int gmx_nonbonded_bench(int argc, char *argv[])
{
    std::vector<const char *>        desc = {
        "[THISMODULE] runs benchmarks for one or more so-called Nbnxm",
        "non-bonded pair kernels. The non-bonded pair kernels are",
        "the most compute intensive part of MD simulations",
        "and usually comprise 60 to 90 percent of the runtime.",
        "For this reason they are highly optimized and several different",
        "setups are available to compute the same physical interactions.",
        "In addition, there are different physical treatments of Coulomb",
        "interactions and optimizations for atoms without Lennard-Jones",
        "interactions. There are also different physical treatments of",
        "Lennard-Jones interactions, but only a plain cut-off is supported",
        "in this tool, as that is by far the most common treatment.",
        "And finally, while force output is always necessary, energy output",
        "is only required at certain steps. In total there are",
        "12 relevant combinations of options. The combinations double to 24",
        "when two different SIMD setups are supported. These combinations",
        "can be run with a single invocation using the [TT]-all[tt] option.",
        "The behavior of each kernel is affected by caching behavior,",
        "which is determined by the hardware used together with the system size",
        "and the cut-off radius. The larger the number of atoms per thread,",
        "the more L1 cache is needed to avoid L1 cache misses.",
        "The cut-off radius mainly affects the data reuse: a larger cut-off",
        "results in more data reuse and makes the kernel less sensitive to cache",
        "misses.[PAR]",
        "OpenMP parallelization is used to utilize multiple hardware threads",
        "within a compute node. In these benchmarks there is no interaction",
        "between threads, apart from starting and closing a single OpenMP",
        "parallel region per iteration. Additionally, threads interact",
        "through sharing and evicting data from shared caches.",
        "The number of threads to use is set with the [TT]-nt[tt] option.",
        "Thread affinity is important, especially with SMT and shared",
        "caches. Affinities can be set through the OpenMP library using",
        "the GOMP_CPU_AFFINITY environment variable.[PAR]",
        "The benchmark tool times one or more kernels by running them",
        "repeatedly for a number of iterations set by the [TT]-iter[tt]",
        "option. An initial kernel call is done to avoid additional initial",
        "cache misses. Times are recording in cycles read from efficient,",
        "high accuracy counters in the CPU. Note that these often do not",
        "correspond to actual clock cycles. For each kernel, the tool",
        "reports the total number of cycles, cycles per iteration,",
        "and (total and useful) pair interactions per cycle.",
        "Because a cluster pair list is used instead of an atom pair list,",
        "interactions are also computed for some atom pairs that are beyond",
        "the cut-off distance. These pairs are not useful (except for",
        "additional buffering, but that is not of interest here),",
        "only a side effect of the cluster-pair setup. The SIMD 2xMM kernel",
        "has a higher useful pair ratio then the 4xM kernel due to a smaller",
        "cluster size, but a lower total pair throughput.",
        "It is best to run this, or for that matter any, benchmark",
        "with locked CPU clocks, as thermal throttling can significantly",
        "affect performance. If that is not an option, the [TT]-warmup[TT]",
        "option can be used to run initial, untimed iterations to warm up",
        "the processor.[PAR]",
        "The most relevant regime is between 0.1 to 1 millisecond per",
        "iteration. Thus it is useful to run with system sizes that cover",
        "both ends of this regime.[PAR]",
        "The [TT]-simd[tt] and [TT]-table[tt] options select different",
        "implementations to compute the same physics. The choice of these",
        "options should ideally be optimized for the target hardware.",
        "Historically, we only found tabulated Ewald correction to be useful",
        "on 2-wide SIMD or 4-wide SIMD without FMA support. As all modern",
        "architectures are wider and support FMA, we do not use tables by",
        "default. The only exceptions are kernels without SIMD, which only",
        "support tables.",
        "Options [TT]-coulomb[tt], [TT]-combrule[tt] and [TT]-halflj[tt]",
        "depend on the force field and composition of the simulated system.",
        "The optimization of computing Lennard-Jones interactions for only",
        "half of the atoms in a cluster is useful for water, which does not",
        "use Lennard-Jones on hydrogen atoms in most water models.",
        "In the MD engine, any clusters where at most half of the atoms",
        "have LJ interactions will automatically use this kernel.",
        "And finally, the [TT]-energy[tt] option selects the computation",
        "of energies, which are usually only needed infrequently."
    };

    static int                       sizeFactor     = 1;
    static Nbnxm::KernelBenchOptions options;

    const char *nbnxmsimdStrings[Nbnxm::nbnxmsimdNR + 1] =
    { nullptr, "auto", "no", "4xm", "2xmm", nullptr };

    const char *combruleStrings[Nbnxm::combruleNR + 1] =
    { nullptr, "geometric", "lb", "none", nullptr };

    const char *coulombTypeStrings[4] =
    { nullptr, "ewald", "reaction-field", nullptr };

    t_pargs           pa[] = {
        { "-size", FALSE, etINT,
          { &sizeFactor }, "The system size is 3000 atoms times this value" },
        { "-nt", FALSE, etINT,
          { &options.numThreads}, "The number of OpenMP threads to use" },
        { "-simd", FALSE, etENUM,
          { nbnxmsimdStrings }, "SIMD type, auto runs all supported SIMD setups or no SIMD when SIMD is not supported" },
        { "-coulomb", FALSE, etENUM,
          { coulombTypeStrings }, "The functional form for the Coulomb interactions" },
        { "-table", FALSE, etBOOL,
          { &options.useTabulatedEwaldCorr }, "Use lookup table for Ewald correction instead of analytical" },
        { "-combrule", FALSE, etENUM,
          { combruleStrings }, "The LJ combination rule" },
        { "-halflj", FALSE, etBOOL,
          { &options.useHalfLJOptimization }, "Use optimization for LJ on half of the atoms" },
        { "-energy", FALSE, etBOOL,
          { &options.computeVirialAndEnergy }, "Compute energies in addition to forces" },
        { "-all", FALSE, etBOOL,
          { &options.doAll }, "Run all 12 combinations of options for coulomb, halflj, combrule" },
        { "-cutoff", FALSE, etREAL,
          { &options.pairlistCutoff }, "Pair-list and interaction cut-off distance" },
        { "-iter", FALSE, etINT,
          { &options.numIterations}, "The number of iterations for each kernel" },
        { "-warmup", FALSE, etINT,
          { &options.numWarmupIterations}, "The number of iterations for initial warmup" },
        { "-cycles", FALSE, etBOOL,
          { &options.cyclesPerPair }, "Report cycles/pair instead of pairs/cycle" }
    };

    gmx_output_env_t *oenv;
    if (!parse_common_args(&argc, argv, 0,
                           0, nullptr, asize(pa), pa, desc.size(), desc.data(), 0, nullptr,
                           &oenv))
    {
        return 0;
    }
    options.nbnxmsimd         = nenum(nbnxmsimdStrings);

    options.ljCombinationRule = nenum(combruleStrings);

    options.coulombType       = nenum(coulombTypeStrings);

    options.nbnxmsimd         = nenum(nbnxmsimdStrings);

    // We compute the Ewald coefficient here to avoid a dependency of the Nbnxm on the Ewald module
    const real ewald_rtol     = 1e-5;
    options.ewaldcoeff_q      = calc_ewaldcoeff_q(options.pairlistCutoff,
                                                  ewald_rtol);

    Nbnxm::bench(sizeFactor, options);

    return 0;
}

}
