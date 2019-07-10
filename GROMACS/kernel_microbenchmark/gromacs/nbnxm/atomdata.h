/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_NBNXN_ATOMDATA_H
#define GMX_NBNXN_ATOMDATA_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace gmx
{
class MDLogger;
}

struct nbnxn_atomdata_t;
struct nonbonded_verlet_t;
struct t_mdatoms;
struct tMPI_Atomic;

namespace Nbnxm
{
class GridSet;
enum class KernelType;
}

enum {
    nbatXYZ, nbatXYZQ, nbatX4, nbatX8
};

//! Stride for coordinate/force arrays with xyz coordinate storage
static constexpr int STRIDE_XYZ  = 3;
//! Stride for coordinate/force arrays with xyzq coordinate storage
static constexpr int STRIDE_XYZQ = 4;
//! Size of packs of x, y or z with SIMD 4-grouped packed coordinates/forces
static constexpr int c_packX4    = 4;
//! Size of packs of x, y or z with SIMD 8-grouped packed coordinates/forces
static constexpr int c_packX8    = 8;
//! Stridefor a pack of 4 coordinates/forces
static constexpr int STRIDE_P4   = DIM*c_packX4;
//! Stridefor a pack of 8 coordinates/forces
static constexpr int STRIDE_P8   = DIM*c_packX8;

// Struct that holds force and energy output buffers
struct nbnxn_atomdata_output_t
{
    std::vector<real> f;      // f, size natoms*fstride
    std::vector<real> fshift; // Shift force array, size SHIFTS*DIM
    std::vector<real> Vvdw;   // Temporary Van der Waals group energy storage
    std::vector<real> Vc;     // Temporary Coulomb group energy storage
    std::vector<real>   VSvdw;  // Temporary SIMD Van der Waals group energy storage
    std::vector<real>   VSc;    // Temporary SIMD Coulomb group energy storage
};

using Outputs = std::vector<nbnxn_atomdata_output_t>;

/* Flags for telling if threads write to force output buffers */
#if false
typedef struct {
    int               nflag;       /* The number of flag blocks                         */
    gmx_bitmask_t    *flag;        /* Bit i is set when thread i writes to a cell-block */
    int               flag_nalloc; /* Allocation size of cxy_flag                       */
} nbnxn_buffer_flags_t;
#endif

/* LJ combination rules: geometric, Lorentz-Berthelot, none */
enum {
    ljcrGEOM, ljcrLB, ljcrNONE, ljcrNR
};

/* Struct that stores atom related data for the nbnxn module
 *
 * Note: performance would improve slightly when all std::vector containers
 *       in this struct would not initialize during resize().
 */
struct nbnxn_atomdata_t
{   //NOLINT(clang-analyzer-optin.performance.Padding)
    struct Params
    {
        // The number of different atom types
        int                   numTypes;
        // Lennard-Jone 6*C6 and 12*C12 parameters, size numTypes*2*2
        std::vector<real> nbfp;
        // Combination rule, see enum defined above
        int                   comb_rule;
        // LJ parameters per atom type, size numTypes*2
        std::vector<real> nbfp_comb;
        // As nbfp, but with a stride for the present SIMD architecture
        std::vector<real>   nbfp_aligned;
        // Atom types per atom
        std::vector<int>  type;
        // LJ parameters per atom for fast SIMD loading
        std::vector<real> lj_comb;
        // Charges per atom, not set with format nbatXYZQ
        std::vector<real> q;
        // The number of energy groups
        int                   nenergrp;
        // 2log(nenergrp)
        int                   neg_2log;
        // The energy groups, one int entry per cluster, only set when needed
        std::vector<int>  energrp;
    };

    // Diagonal and topology exclusion helper data for all SIMD kernels
    struct SimdMasks
    {
        // Helper data for setting up diagonal exclusion masks in the SIMD 4xN kernels
        std::vector<real>     diagonal_4xn_j_minus_i;
        // Helper data for setting up diaginal exclusion masks in the SIMD 2xNN kernels
        std::vector<real>     diagonal_2xnn_j_minus_i;
        // Filters for topology exclusion masks for the SIMD kernels
        std::vector<uint32_t> exclusion_filter;
        // Filters for topology exclusion masks for double SIMD kernels without SIMD int32 logical support
        std::vector<uint64_t> exclusion_filter64;
        // Array of masks needed for exclusions
        std::vector<real>     interaction_array;
    };

        // The LJ and charge parameters
        Params                     params_;
        // The total number of atoms currently stored
        int                        numAtoms_;
        int                        natoms_local; /* Number of local atoms                           */
        int                        XFormat;      /* The format of x (and q), enum                      */
        int                        FFormat;      /* The format of f, enum                              */
        gmx_bool                   bDynamicBox;  /* Do we need to update shift_vec every step?    */
        std::vector<gmx::RVec> shift_vec;    /* Shift vectors, copied from t_forcerec              */
        int                        xstride;      /* stride for a coordinate in x (usually 3 or 4)      */
        int                        fstride;      /* stride for a coordinate in f (usually 3 or 4)      */
        std::vector<real>      x_;           /* x and possibly q, size natoms*xstride              */

        // Masks for handling exclusions in the SIMD kernels
        const SimdMasks          simdMasks;

        /* Output data */
        std::vector<nbnxn_atomdata_output_t> out; /* Output data structures, 1 per thread */

        /* Reduction related data */
        gmx_bool                 bUseBufferFlags;     /* Use the flags or operate on all atoms     */
        //nbnxn_buffer_flags_t     buffer_flags;        /* Flags for buffer zeroing+reduc.  */
        gmx_bool                 bUseTreeReduce;      /* Use tree for force reduction */
};

enum {
    enbnxninitcombruleDETECT, enbnxninitcombruleGEOM, enbnxninitcombruleLB, enbnxninitcombruleNONE
};

#endif
