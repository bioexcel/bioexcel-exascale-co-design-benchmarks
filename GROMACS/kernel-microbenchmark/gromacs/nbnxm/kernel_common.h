/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief
 * Declares the nbnxm pair interaction kernel function types and kind counts, also declares utility functions used in nbnxm_kernel.cpp.
 *
 * \author Berk Hess <hess@kth.se>
 */

#ifndef GMX_NBXNM_KERNEL_COMMON_H
#define GMX_NBXNM_KERNEL_COMMON_H

#include "gromacs/math/vectypes.h"

struct NbnxnPairlistCpu;
struct interaction_const_t;
struct nbnxn_atomdata_t;
struct nbnxn_atomdata_output_t;

// TODO: Consider using one nbk_func type now ener and noener are identical

/*! \brief Pair-interaction kernel type that also calculates energies.
 */
typedef void (nbk_func_ener)(const NbnxnPairlistCpu     *nbl,
                             const nbnxn_atomdata_t     *nbat,
                             const interaction_const_t  *ic,
                             const rvec                 *shift_vec,
                             nbnxn_atomdata_output_t    *out);

/*! \brief Pointer to \p nbk_func_ener.
 */
typedef nbk_func_ener *p_nbk_func_ener;

/*! \brief Pair-interaction kernel type that does not calculates energies.
 */
typedef void (nbk_func_noener)(const NbnxnPairlistCpu     *nbl,
                               const nbnxn_atomdata_t     *nbat,
                               const interaction_const_t  *ic,
                               const rvec                 *shift_vec,
                               nbnxn_atomdata_output_t    *out);

/*! \brief Pointer to \p nbk_func_noener.
 */
typedef nbk_func_noener *p_nbk_func_noener;

#endif
