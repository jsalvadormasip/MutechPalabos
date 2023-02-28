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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H

#include "complexDynamics/advectionDiffusionBoundaryCondition2D.h"
#include "core/cell.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataInitializerFunctional2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator2D :
    public OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor> {
public:
    AdvectionDiffusionBoundaryConditionInstantiator2D();

    void addTemperatureBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    void addTemperatureCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    void addTemperatureCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

private:
    template <int direction, int orientation>
    void addTemperatureBoundary(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addTemperatureCorner(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);

    template <int direction, int orientation>
    void addTemperatureBoundary(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
    template <int normalX, int normalY>
    void addTemperatureCorner(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType);
};

///////// class AdvectionDiffusionBoundaryConditionInstantiator2D ////////////////////////

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator2D<
    T, Descriptor, BoundaryManager>::AdvectionDiffusionBoundaryConditionInstantiator2D()
{ }

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getAdvectionDiffusionBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>, domain,
            lattice, bcType);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getAdvectionDiffusionBoundaryProcessor<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCorner(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getAdvectionDiffusionCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyDensityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getAdvectionDiffusionCornerProcessor<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0 == domain.x1 || domain.y0 == domain.y1);

    setCompositeDynamics(
        lattice, domain,
        BoundaryManager::template getAdvectionDiffusionBoundaryDynamics<direction, orientation>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation>, domain,
            lattice, bcType);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getAdvectionDiffusionBoundaryProcessor<direction, orientation>();
    if (functional) {
        integrateProcessingFunctional(functional, domain, lattice, bcType);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
template <int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCorner(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    setCompositeDynamics(
        lattice, Box2D(x, x, y, y),
        BoundaryManager::template getAdvectionDiffusionCornerDynamics<xNormal, yNormal>(
            new NoDynamics<T, Descriptor>));

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType == boundary::neumann) {
        integrateProcessingFunctional(
            new CopyDensityFunctional2D<T, Descriptor, xNormal, yNormal>, Box2D(x, x, y, y),
            lattice);
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data
    // processor.
    BoxProcessingFunctional2D_L<T, Descriptor> *functional =
        BoundaryManager::template getAdvectionDiffusionCornerProcessor<xNormal, yNormal>();
    if (functional) {
        integrateProcessingFunctional(functional, Box2D(x, x, y, y), lattice);
    }
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<0, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, -1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureBoundary<1, 1>(domain, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<-1, 1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, -1>(x, y, lattice, bcType);
}

template <typename T, template <typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T, Descriptor, BoundaryManager>::
    addTemperatureCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType)
{
    addTemperatureCorner<1, 1>(x, y, lattice, bcType);
}

}  // namespace plb

#endif
