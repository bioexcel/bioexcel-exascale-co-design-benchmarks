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
 * Declares the Grid class.
 *
 * This class provides functionality for setting up and accessing atoms
 * on a grid for one domain decomposition zone. This grid is used for
 * generating cluster pair lists for computing non-bonded pair interactions.
 * The grid consists of a regular array of columns along dimensions x and y.
 * Along z the number of cells and their boundaries vary between the columns.
 * Each cell can hold one or more clusters of atoms, depending on the grid
 * geometry, which is set by the pair-list type.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GRID_H
#define GMX_NBNXM_GRID_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"


struct gmx_domdec_zones_t;
struct nbnxn_atomdata_t;
struct nbnxn_search;

namespace Nbnxm
{

/*! \internal
 * \brief Bounding box for a nbnxm atom cluster
 *
 * \note Should be aligned in memory to enable 4-wide SIMD operations.
 */
struct BoundingBox
{
    /*! \internal
     * \brief Corner for the bounding box, padded with one element to enable 4-wide SIMD operations
     */
    struct Corner
    {
        float x;       //!< x coordinate
        float y;       //!< y coordinate
        float z;       //!< z coordinate
        float padding; //!< padding, unused, but should be set to avoid operations on unitialized data
    };

    Corner lower; //!< lower, along x and y and z, corner
    Corner upper; //!< upper, along x and y and z, corner
};

/*! \internal
 * \brief Bounding box for one dimension of a grid cell
 */
struct BoundingBox1D
{
    float lower; //!< lower bound
    float upper; //!< upper bound
};

/*! \brief The number of bounds along one dimension of a bounding box */
static constexpr int c_numBoundingBoxBounds1D = 2;

} // namespace Nbnxm

namespace Nbnxm
{

/*! \internal
 * \brief Helper struct to pass data that is shared over all grids
 *
 * To enable a single coordinate and force array, a single cell range
 * is needed which covers all grids. This helper struct contains
 * references to the index lists mapping both ways, as well as
 * the free-energy boolean, which is the same for all grids.
 */
struct GridSetData
{
    //! The cell indices for all atoms
    std::vector<int> &cells;
    //! The atom indices for all atoms stored in cell order
    std::vector<int> &atomIndices;
    //! Tells whether we are have perturbed non-bonded interations
    const bool        haveFep;
};

/*! \internal
 * \brief Working arrays for constructing a grid
 */
struct GridWork
{
    //! Number of atoms for each grid column
    std::vector<int> numAtomsPerColumn;
    //! Buffer for sorting integers
    std::vector<int> sortBuffer;
};

/*! \internal
 * \brief A pair-search grid object for one domain decomposition zone
 *
 * This is a rectangular 3D grid covering a potentially non-rectangular
 * volume which is either the whole unit cell or the local zone or part
 * of a non-local zone when using domain decomposition. The grid cells
 * are even spaced along x/y and irregular along z. Each cell is sub-divided
 * into atom clusters. With a CPU geometry, each cell contains 1 or 2 clusters.
 * With a GPU geometry, each cell contains up to 8 clusters. The geometry is
 * set by the pairlist type which is the only argument of the constructor.
 *
 * When multiple grids are used, i.e. with domain decomposition, we want
 * to avoid the overhead of multiple coordinate arrays or extra indexing.
 * Therefore each grid stores a cell offset, so a contiguous cell index
 * can be used to index atom arrays. All methods returning atom indices
 * return indices which index into a full atom array.
 *
 * Note that when atom groups, instead of individual atoms, are assigned
 * to grid cells, individual atoms can be geometrically outside the cell
 * and grid that they have been assigned to (as determined by the center
 * or geometry of the atom group they belong to).
 */
class Grid
{
    public:
        /*! \internal
         * \brief The cluster and cell geometry of a grid
         */
        struct Geometry
        {
            bool isSimple;             //!< Is this grid simple (CPU) or hierarchical (GPU)
            int  numAtomsICluster;     //!< Number of atoms per cluster
            int  numAtomsJCluster;     //!< Number of atoms for list j-clusters
            int  numAtomsPerCell;      //!< Number of atoms per cell
            int  numAtomsICluster2Log; //!< 2log of na_c
        };

        // The physical dimensions of a grid
        struct Dimensions
        {
            //! The lower corner of the (local) grid
            rvec lowerCorner;
            //! The upper corner of the (local) grid
            rvec upperCorner;
            //! The physical grid size: upperCorner - lowerCorner
            rvec gridSize;
            //! An estimate for the atom number density of the region targeted by the grid
            real atomDensity;
            //! The maximum distance an atom can be outside of a cell and outside of the grid
            real maxAtomGroupRadius;
            //! Size of cell along dimension x and y
            real cellSize[DIM - 1];
            //! 1/size of a cell along dimensions x and y
            real invCellSize[DIM - 1];
            //! The number of grid cells along dimensions x and y
            int  numCells[DIM - 1];
        };

        /* Data members */
        //! The geometry of the grid clusters and cells
        Geometry   geometry_;
        //! The physical dimensions of the grid
        Dimensions dimensions_;

        //! The total number of cells in this grid
        int        numCellsTotal_;
        //! Index in nbs->cell corresponding to cell 0
        int        cellOffset_;
        //! The maximum number of cells in a column
        int        numCellsColumnMax_;

        //! The start of the source atom range mapped to this grid
        int        srcAtomBegin_;
        //! The end of the source atom range mapped to this grid
        int        srcAtomEnd_;

        /* Grid data */
        /*! \brief The number of, non-filler, atoms for each grid column.
         *
         * \todo Needs a useful name. */
        std::vector<int> cxy_na_;
        /*! \brief The grid-local cell index for each grid column
         *
         * \todo Needs a useful name. */
        std::vector<int> cxy_ind_;

        //! The number of cluster for each cell
        std::vector<int> numClusters_;

        /* Bounding boxes */
        //! Bounding boxes in z for the cells
        std::vector<BoundingBox1D>                          bbcz_;
        //! 3D bounding boxes for the sub cells
        std::vector < BoundingBox > bb_;
        //! 3D j-bounding boxes
        std::vector < BoundingBox > bbj_;

        /* Bit-flag information */
        //! Flags for properties of clusters in each cell
        std::vector<int>          flags_;
        //! Signal bits for atoms in each cell that tell whether an atom is perturbed
        std::vector<unsigned int> fep_;

        /* Statistics */
        //! Total number of clusters, used for printing
        int numClustersTotal_;
};

} // namespace Nbnxm

#endif
