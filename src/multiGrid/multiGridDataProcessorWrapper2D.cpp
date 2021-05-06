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

#include "multiGrid/multiGridDataProcessorWrapper2D.h"

#include "atomicBlock/dataProcessor2D.h"
#include "core/plbDebug.h"
#include "multiGrid/multiGridOperations2D.h"

namespace plb {

/* *************** BoxProcessing2D, general case *************************** */

void applyProcessingFunctional(
    BoxProcessingFunctional2D *functional, Box2D domain, std::vector<MultiGrid2D *> multiBlocks,
    plint referenceLevel)
{
    executeDataProcessor(BoxProcessorGenerator2D(functional, domain), multiBlocks, referenceLevel);
}

void integrateProcessingFunctional(
    BoxProcessingFunctional2D *functional, Box2D domain, std::vector<MultiGrid2D *> multiBlocks,
    plint referenceLevel, plint level)
{
    addInternalProcessor(
        BoxProcessorGenerator2D(functional, domain), multiBlocks, level, referenceLevel);
}

/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(
    DotProcessingFunctional2D *functional, DotList2D const &dotList,
    std::vector<MultiGrid2D *> multiBlocks, plint referenceLevel)
{
    executeDataProcessor(DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel);
}

void integrateProcessingFunctional(
    DotProcessingFunctional2D *functional, DotList2D const &dotList,
    std::vector<MultiGrid2D *> multiBlocks, plint level, plint referenceLevel)
{
    addInternalProcessor(
        DotProcessorGenerator2D(functional, dotList), multiBlocks, referenceLevel, level);
}

/* *************** BoundedBoxProcessing2D, general case *************************** */

void applyProcessingFunctional(
    BoundedBoxProcessingFunctional2D *functional, Box2D domain,
    std::vector<MultiGrid2D *> multiBlocks, plint boundaryWidth, plint referenceLevel)
{
    std::vector<BoxProcessorGenerator2D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        executeDataProcessor(*generators[iGen], multiBlocks, referenceLevel);
        delete generators[iGen];
    }
}

void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional2D *functional, Box2D domain,
    std::vector<MultiGrid2D *> multiBlocks, plint boundaryWidth, plint referenceLevel, plint level)
{
    std::vector<BoxProcessorGenerator2D *> generators;
    functional->getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen = 0; iGen < generators.size(); ++iGen) {
        addInternalProcessor(*generators[iGen], multiBlocks, level, referenceLevel);
        delete generators[iGen];
    }
}

}  // namespace plb
