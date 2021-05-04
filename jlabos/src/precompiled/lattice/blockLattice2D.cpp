/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * The dynamics of a 2D block lattice -- template instantiation.
 */
#ifdef COMPILE_2D

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/blockLattice2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

    template class BlockLattice2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
    template class BlockLatticeDataTransfer2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;

    template double getStoredAverageDensity<FLOAT_T, descriptors::DESCRIPTOR_2D>
        (BlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> const& blockLattice);
    template double getStoredAverageEnergy<FLOAT_T, descriptors::DESCRIPTOR_2D>
        (BlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> const& blockLattice);
    template double getStoredMaxVelocity<FLOAT_T, descriptors::DESCRIPTOR_2D>
        (BlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> const& blockLattice);

}

#endif  // COMPILE_2D