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
 * \brief
 * Declares the GridSet class.
 *
 * This class holds the grids for the local and non-local domain decomposition
 * zones, as well as the cell and atom data that covers all grids.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRIDSET_H
#define GMX_NBNXM_GRIDSET_H

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"

#include "grid.h"

namespace Nbnxm
{

/*! \internal
 * \brief Holds a set of search grids for the local + non-local DD zones
 */
class GridSet
{
    public:
        //! The search grids
        std::vector<Grid>     grids_;
        //! The actual cell indices for all atoms, covering all grids
        std::vector<int>      cells_;
        //! The actual array of atom indices, covering all grids
        std::vector<int>      atomIndices_;
        //! Tells whether we have perturbed non-bonded interactions
        bool                  haveFep_;
        //! The periodic unit-cell
        matrix                box_;
        //! The number of local real atoms, i.e. without padded atoms, local atoms: 0 to numAtomsLocal_
        int                   numRealAtomsLocal_;
        //! The total number of real atoms, i.e. without padded atoms
        int                   numRealAtomsTotal_;
        //! Working data for constructing a single grid, one entry per thread
        std::vector<GridWork> gridWork_;

};

} // namespace Nbnxm

#endif
