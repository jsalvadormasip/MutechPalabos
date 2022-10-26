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

#ifndef ELIULT_OFF_LATTICE_MODEL_3D_HH
#define ELIULT_OFF_LATTICE_MODEL_3D_HH

#include "offLattice/ELIULToffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "core/dynamics.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace plb {

/**
 * See (Marson et al., 2021)
 *
 * (Marson et al., 2021) Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., Latt, J., 2021.
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows. Phys. Rev. E 103, 053308.
 * https://doi.org/10.1103/PhysRevE.103.053308
 *
 * @tparam T
 * @tparam Descriptor
 */
template<typename T, template<typename U> class Descriptor>
ELIULTLocalModel3D<T,Descriptor>::ELIULTLocalModel3D (
        BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_)
    : OffLatticeModel3D<T,Array<T,3> >(shape_, flowType_)
{ }

template<typename T, template<typename U> class Descriptor>
ELIULTLocalModel3D<T,Descriptor>* ELIULTLocalModel3D<T,Descriptor>::clone() const {
    return new ELIULTLocalModel3D(*this);
}

template<typename T, template<typename U> class Descriptor>
plint ELIULTLocalModel3D<T,Descriptor>::getNumNeighbors() const {
    return 1;
}

/**
 *  Filippova-Haenel is a completion scheme for a layer of cells on the
 *  "solid" side of the boundary, like Guo.
 * @tparam T
 * @tparam Descriptor
 * @return true
 */
template<typename T, template<typename U> class Descriptor>
bool ELIULTLocalModel3D<T,Descriptor>::isExtrapolated() const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
void ELIULTLocalModel3D<T,Descriptor>::prepareCell (
        Dot3D const& cellLocation,
        AtomicContainerBlock3D& container )
{
    typedef Descriptor<T> D;
    Dot3D offset = container.getLocation();
    auto* info = dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<int> liquidNeighbors;
    std::vector<plint> ids;
    Dot3D absLoc = cellLocation+offset;
    if (this->isSolid(absLoc)) {
        for (int iPop=0; iPop<D::q; ++iPop) {
            Dot3D neighbor(cellLocation.x+D::c[iPop][0], cellLocation.y+D::c[iPop][1], cellLocation.z+D::c[iPop][2]);
            Dot3D neighborLoc = neighbor + offset;
            // If the non-fluid node has a fluid neighbor ...
            if (this->isFluid(neighborLoc)) {
                // ... check how many fluid nodes it has ahead of it ...
                plint iTriangle=-1;
                global::timer("intersect").start();
                Array<T,3> locatedPoint;
                T distance;
                Array<T,3> wallNormal;
                Array<T,3> surfaceData;
                OffBoundary::Type bdType;
#ifdef PLB_DEBUG
                bool ok =
#endif
                    this->pointOnSurface (
                            cellLocation+offset, Dot3D(D::c[iPop][0],D::c[iPop][1],D::c[iPop][2]), locatedPoint, distance,
                            wallNormal, surfaceData, bdType, iTriangle );
                    
                // In the following, the importance of directions is sorted wrt. how well they
                //   are aligned with the wall normal. It is better to take the continuous normal,
                //   because it is not sensitive to the choice of the triangle when we shoot at
                //   an edge.
                //wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
                global::timer("intersect").stop();
                PLB_ASSERT( ok );
                // ... then add this node to the list.
                liquidNeighbors.push_back(iPop);
                ids.push_back(iTriangle);
            }
        }
        if (!liquidNeighbors.empty()) {
            info->getDryNodes().push_back(cellLocation);
            info->getDryNodeFluidDirections().push_back(liquidNeighbors);
            info->getDryNodeIds().push_back(ids);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ContainerBlockData*
    ELIULTLocalModel3D<T,Descriptor>::generateOffLatticeInfo() const
{
    return new OffLatticeInfo3D;
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> ELIULTLocalModel3D<T,Descriptor>::getLocalForce (
                AtomicContainerBlock3D& container ) const
{
    auto* info =
        dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    return info->getLocalForce();
}

template<typename T, template<typename U> class Descriptor>
void ELIULTLocalModel3D<T,Descriptor>::boundaryCompletion (
        AtomicBlock3D& nonTypeLattice,
        AtomicContainerBlock3D& container,
        std::vector<AtomicBlock3D *> const& args )
{
    auto& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&> (nonTypeLattice);
    auto* info =
        dynamic_cast<OffLatticeInfo3D*>(container.getData());
    PLB_ASSERT( info );
    std::vector<Dot3D> const&
        dryNodes = info->getDryNodes();
    std::vector<std::vector<int > > const&
        dryNodeFluidDirections = info->getDryNodeFluidDirections();
    std::vector<std::vector<plint> > const&
        dryNodeIds = info->getDryNodeIds();
    PLB_ASSERT( dryNodes.size() == dryNodeFluidDirections.size() );

    Dot3D absoluteOffset = container.getLocation();

    Array<T,3>& localForce = info->getLocalForce();
    localForce.resetToZero();
    for (pluint iDry=0; iDry<dryNodes.size(); ++iDry) {
        cellCompletion (
            lattice, dryNodes[iDry], dryNodeFluidDirections[iDry],
            dryNodeIds[iDry], absoluteOffset, localForce, args );
    }
}

template<typename T, template<typename U> class Descriptor>
void ELIULTLocalModel3D<T,Descriptor>::cellCompletion (
        BlockLattice3D<T,Descriptor>& lattice, Dot3D const& guoNode,
        std::vector<int> const& dryNodeFluidDirections,
        std::vector<plint> const& dryNodeIds, Dot3D const& absoluteOffset,
        Array<T,3>& localForce, std::vector<AtomicBlock3D *> const& args )
{
    typedef Descriptor<T> D;
    Cell<T,Descriptor>& cellS =
        lattice.get( guoNode.x, guoNode.y, guoNode.z );
#ifdef PLB_DEBUG
    int noDynId =
        NoDynamics<T,Descriptor>().getId();
#endif
    PLB_ASSERT(cellS.getDynamics().getId() == noDynId
               && "Filippova-Haenel BC needs the dynamics to be set to NoDynamics.");
    for (plint iDirection=0; iDirection<(plint)dryNodeFluidDirections.size(); ++iDirection)
    {
        int iOpp = dryNodeFluidDirections[iDirection];
        int iPop = indexTemplates::opposite<Descriptor<T> >(iOpp);
        Dot3D fluidDirection(D::c[iOpp][0],D::c[iOpp][1],D::c[iOpp][2]);
        plint dryNodeId = dryNodeIds[iDirection];

        Array<T,3> wallNode, wall_vel;
        T wallDistance;
        OffBoundary::Type bdType;
        Cell<T,Descriptor> const& cellF =
            lattice.get( guoNode.x+fluidDirection.x,
                         guoNode.y+fluidDirection.y,
                         guoNode.z+fluidDirection.z );

        Cell<T,Descriptor> cellFstar(cellF);
        BlockStatistics statsCopy(lattice.getInternalStatistics());
        cellFstar.collide(statsCopy);

        T f_rhoBar, ff_rhoBar;
        Array<T,3> f_j, ff_j;
        Array<T,3> wallNormal;

        if (args.empty()) {
            Cell<T,Descriptor> const& cellFF =
            lattice.get( guoNode.x+2*fluidDirection.x,
                         guoNode.y+2*fluidDirection.y,
                         guoNode.z+2*fluidDirection.z );

            cellF.getDynamics().computeRhoBarJ(cellF, f_rhoBar, f_j);
            cellFF.getDynamics().computeRhoBarJ(cellFF, ff_rhoBar, ff_j);
        } else {
            if ((plint) args.size() == 2) {
                // 1 field for rhoBar, 1 field for j.
                // #ifdef PLB_DEBUG
                auto const* rhoBarField =
                    dynamic_cast<ScalarField3D<T> const*>( args[0] );
                // #endif
                auto const* jField =
                    dynamic_cast<TensorField3D<T,3> const*>( args[1] );
                PLB_ASSERT( rhoBarField );
                PLB_ASSERT( jField );

                Dot3D posJ = guoNode + fluidDirection+computeRelativeDisplacement(lattice, *jField);
                Dot3D posRho = guoNode + fluidDirection+computeRelativeDisplacement(lattice, *rhoBarField);

                f_rhoBar = rhoBarField->get(posRho.x, posRho.y, posRho.z);
                f_j = jField->get(posJ.x, posJ.y, posJ.z);
                Dot3D posJ2 = posJ + fluidDirection;

                ff_j = jField->get(posJ2.x, posJ2.y, posJ2.z);
            } else {
                PLB_ASSERT(false); // Not implemented for 1 arg.
            }
        }
        
        T f_rho = D::fullRho(f_rhoBar);
        T f_jSqr = normSqr(f_j);

#ifdef PLB_DEBUG
        bool ok =
#endif
        this->pointOnSurface( guoNode+absoluteOffset, fluidDirection,
                              wallNode, wallDistance, wallNormal,
                              wall_vel, bdType, dryNodeId );
        PLB_ASSERT( ok );

        Array<T,3> w_j = wall_vel*f_rho;
        T d = std::sqrt(D::cNormSqr[iOpp]);
        PLB_ASSERT( wallDistance <= d );
        T q = 1.0 - wallDistance / d;

        Array<T,3> wf_j; wf_j.resetToZero();
        T omega = cellF.getDynamics().getOmega(); // TODO: add TRT case with omega- and omega+
//
        T a3 = q / (1. + q);
        T a4 = 1. / (1. + q);
        T f_1plusq = cellFstar[iOpp];
        T f_0 = cellF.getDynamics().computeEquilibrium(iOpp,f_rhoBar,w_j,normSqr(w_j))+
                (cellF[iPop]-cellF.getDynamics().computeEquilibrium(iPop,f_rhoBar,f_j,f_jSqr));
        cellS[iOpp] = a3 * f_1plusq + a4 * f_0;

        localForce[0] += D::c[iPop][0]*(cellS[iPop] - cellS[iOpp]);
        localForce[1] += D::c[iPop][1]*(cellS[iPop] - cellS[iOpp]);
        localForce[2] += D::c[iPop][2]*(cellS[iPop] - cellS[iOpp]);
    }
}

}  // namespace plb

#endif  // ELIULT_OFF_LATTICE_MODEL_3D_HH
