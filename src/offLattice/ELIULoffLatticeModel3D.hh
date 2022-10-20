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

#ifndef ELIUL_OFF_LATTICE_MODEL_3D_HH
#define ELIUL_OFF_LATTICE_MODEL_3D_HH

#include "offLattice/ELIULoffLatticeModel3D.h"
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
ELIULlocalModel3D<T,Descriptor>::ELIULlocalModel3D (
        BoundaryShape3D<T,Array<T,3> >* shape_, int flowType_,ParametrizationType param_,NonEquilibriumType neq_t_)
    : OffLatticeModel3D<T,Array<T,3> >(shape_, flowType_), param(param_), neq_t(neq_t_)
{
   typedef Descriptor<T> D;
   invAB.resize(D::q);
   invAB[0] = T();
   for (plint iPop = 1; iPop < D::q; ++iPop) {
       invAB[iPop] =
           (T)1
           / std::sqrt(
               util::sqr(D::c[iPop][0]) + util::sqr(D::c[iPop][1]) + util::sqr(D::c[iPop][2]));
   }
}

template <typename T, template <typename U> class Descriptor>
ELIULlocalModel3D<T, Descriptor> *ELIULlocalModel3D<T, Descriptor>::clone() const
{
   return new ELIULlocalModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint ELIULlocalModel3D<T, Descriptor>::getNumNeighbors() const
{
   return 1;
}

template <typename T, template <typename U> class Descriptor>
bool ELIULlocalModel3D<T, Descriptor>::isExtrapolated() const
{
   // Bouzidi is a completion scheme for a layer of cells on the
   // "fluid" side of the boundary, unlike Guo.
   return false;
}

template <typename T, template <typename U> class Descriptor>
void ELIULlocalModel3D<T, Descriptor>::prepareCell(
   Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
   typedef Descriptor<T> D;
   Dot3D offset = container.getLocation();
   ELIULOffLatticeInfo3D *info = dynamic_cast<ELIULOffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   std::vector<int> solidDirections;
   std::vector<plint> boundaryIds;
   std::vector<bool> hasFluidNeighbor;
   if (this->isFluid(cellLocation + offset)) {
       for (plint iPop = 1; iPop < D::q; ++iPop) {
           Dot3D neighbor(
               cellLocation.x + D::c[iPop][0], cellLocation.y + D::c[iPop][1],
               cellLocation.z + D::c[iPop][2]);
           Dot3D prevNode(
               cellLocation.x - D::c[iPop][0], cellLocation.y - D::c[iPop][1],
               cellLocation.z - D::c[iPop][2]);
           // If the fluid node has a non-fluid neighbor ...
           if (this->isSolid(neighbor + offset)) {
               plint iTriangle = -1;
               global::timer("intersect").start();
               Array<T, 3> locatedPoint;
               T distance;
               Array<T, 3> wallNormal;
               Array<T, 3> surfaceData;
               OffBoundary::Type bdType;
#ifdef PLB_DEBUG
               bool ok =
#endif
                   this->pointOnSurface(
                       cellLocation + offset, Dot3D(D::c[iPop][0], D::c[iPop][1], D::c[iPop][2]),
                       locatedPoint, distance, wallNormal, surfaceData, bdType, iTriangle);
               // In the following, the importance of directions is sorted wrt. how well they
               //   are aligned with the wall normal. It is better to take the continuous normal,
               //   because it is not sensitive to the choice of the triangle when we shoot at
               //   an edge.
               global::timer("intersect").stop();
               PLB_ASSERT(ok);
               // ... then add this node to the list.
               solidDirections.push_back(iPop);
               boundaryIds.push_back(iTriangle);
               bool prevNodeIsPureFluid = this->isFluid(prevNode + offset);
               if (prevNodeIsPureFluid) {
                   hasFluidNeighbor.push_back(true);
               } else {
                   hasFluidNeighbor.push_back(false);
               }
           }
       }
       if (!solidDirections.empty()) {
           info->getBoundaryNodes().push_back(cellLocation);
           info->getSolidDirections().push_back(solidDirections);
           info->getBoundaryIds().push_back(boundaryIds);
           info->getHasFluidNeighbor().push_back(hasFluidNeighbor);
       }
   }
}

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *ELIULlocalModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
   return new ELIULOffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> ELIULlocalModel3D<T, Descriptor>::getLocalForce(
   AtomicContainerBlock3D &container) const
{
   ELIULOffLatticeInfo3D *info = dynamic_cast<ELIULOffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void ELIULlocalModel3D<T, Descriptor>::boundaryCompletion(
   AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
   std::vector<AtomicBlock3D *> const &args)
{
   BlockLattice3D<T, Descriptor> &lattice =
       dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
   ELIULOffLatticeInfo3D *info = dynamic_cast<ELIULOffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   std::vector<Dot3D> const &boundaryNodes = info->getBoundaryNodes();
   std::vector<std::vector<int> > const &solidDirections = info->getSolidDirections();
   std::vector<std::vector<plint> > const &boundaryIds = info->getBoundaryIds();
   std::vector<std::vector<bool> > const &hasFluidNeighbor = info->getHasFluidNeighbor();
   PLB_ASSERT(boundaryNodes.size() == solidDirections.size());
   PLB_ASSERT(boundaryNodes.size() == boundaryIds.size());
   PLB_ASSERT(boundaryNodes.size() == hasFluidNeighbor.size());

   Dot3D absoluteOffset = container.getLocation();

   Array<T, 3> &localForce = info->getLocalForce();
   localForce.resetToZero();
   for (pluint i = 0; i < boundaryNodes.size(); ++i) {
       cellCompletion(
           lattice, boundaryNodes[i], solidDirections[i], boundaryIds[i], hasFluidNeighbor[i],
           absoluteOffset, localForce, args);
   }
}

template <typename T, template <typename U> class Descriptor>
void ELIULlocalModel3D<T, Descriptor>::cellCompletion(
   BlockLattice3D<T, Descriptor> &lattice, Dot3D const &boundaryNode,
   std::vector<int> const &solidDirections, std::vector<plint> const &boundaryIds,
   std::vector<bool> const &hasFluidNeighbor, Dot3D const &absoluteOffset, Array<T, 3> &localForce,
   std::vector<AtomicBlock3D *> const &args)
{
   typedef Descriptor<T> D;
   Array<T, D::d> deltaJ;
   deltaJ.resetToZero();

   plint numNeumannNodes = 0;
   T neumannDensity = T();
   Cell<T, Descriptor> &cellF = lattice.get(boundaryNode.x, boundaryNode.y, boundaryNode.z);
   if (this->computesStat()) {
       for (int iPop : solidDirections) {
           deltaJ[0] += D::c[iPop][0] * cellF[iPop];
           deltaJ[1] += D::c[iPop][1] * cellF[iPop];
           deltaJ[2] += D::c[iPop][2] * cellF[iPop];
       }
   }
   for (pluint i = 0; i < solidDirections.size(); ++i) {
       int i_solid = solidDirections[i];
       int i_fluid = indexTemplates::opposite<D>(i_solid);
       Array<T, 3> wallNode, wall_vel;
       T AC;
       OffBoundary::Type bdType;
       Array<T, 3> wallNormal;
       plint id = boundaryIds[i];
#ifdef PLB_DEBUG
       bool ok =
#endif
           this->pointOnSurface(
               boundaryNode + absoluteOffset, Dot3D(D::c[i_solid][0], D::c[i_solid][1], D::c[i_solid][2]),
               wallNode, AC, wallNormal, wall_vel, bdType, id);
       PLB_ASSERT(ok);
       T q = AC * invAB[i_solid];
       Cell<T, Descriptor> &cellS = lattice.get(
           boundaryNode.x + D::c[i_solid][0], boundaryNode.y + D::c[i_solid][1],
           boundaryNode.z + D::c[i_solid][2]);
       Cell<T, Descriptor> &cellFF = lattice.get(
           boundaryNode.x - D::c[i_solid][0], boundaryNode.y - D::c[i_solid][1],
           boundaryNode.z - D::c[i_solid][2]);
       if (bdType == OffBoundary::dirichlet) {
           T u_ci = D::c[i_solid][0] * wall_vel[0] + D::c[i_solid][1] * wall_vel[1]
                    + D::c[i_solid][2] * wall_vel[2];
           plint numUnknown = 0;
           if (hasFluidNeighbor[i]) {
               const auto trt = dynamic_cast<const BaseTRTdynamics<T, Descriptor> *>(&cellF.getDynamics());
               T omega_plus  = cellF.getDynamics().getOmega(); // TODO: add TRT case with omega- and omega+
               T omega_minus = 0.;
               if(trt)
                   omega_minus = cellF.getDynamics().getParameter(dynamicParams::omega_minus);
               else
                   omega_minus = omega_plus;
               T tauPlus = 1./omega_plus;
               T tauMinus = 1./omega_minus;
               T fPlus  = 0.5*(cellS[i_solid] + cellFF[i_fluid]);
               T fMinus = 0.5*(cellS[i_solid] - cellFF[i_fluid]);
               T Kelip = q-tauPlus;
               T Kelim = q-tauMinus;
               T rhoBar;
               Array<T, 3> j;
               momentTemplates<T, Descriptor>::get_rhoBar_j(cellF, rhoBar, j);
               T jSqr = dot(j,j);
               T feq_fluid = cellF.getDynamics().computeEquilibrium(i_fluid, rhoBar, j, jSqr);
               T feq_solid = cellF.getDynamics().computeEquilibrium(i_solid, rhoBar, j, jSqr);
               T eqPlus = 0.5 * (feq_solid-feq_fluid);
               T eqMinus = 0.5 * (feq_solid-feq_fluid);
               cellF[i_fluid] = Kelip * (fPlus-eqPlus) / (tauPlus-1.)
                                + Kelim * (fMinus-eqMinus) / (tauMinus-1.)
                                + eqPlus
//                                - eqMinusWall
                   ;
//               cellS[i_fluid] = a4 * f_0 + a5 * f_1 - (this->param == PT::k2) ? k*f0neq_minus*omega_minus : 0.0;
           } else {
               ++numUnknown;
               cellF[i_fluid] = cellS[i_solid];
           }
       } else {
           // Not implemented yet.
           PLB_ASSERT(false);
       }
   }

   if (this->computesStat()) {
       Cell<T, Descriptor> collidedCell(cellF);
       BlockStatistics statsCopy(lattice.getInternalStatistics());
       collidedCell.collide(statsCopy);

       for (int iPop : solidDirections) {
           int oppPop = indexTemplates::opposite<D>(iPop);
           deltaJ[0] -= D::c[oppPop][0] * collidedCell[oppPop];
           deltaJ[1] -= D::c[oppPop][1] * collidedCell[oppPop];
           deltaJ[2] -= D::c[oppPop][2] * collidedCell[oppPop];
       }
   }

   localForce += deltaJ;
   if (numNeumannNodes > 0) {
       neumannDensity /= numNeumannNodes;
       T oldRhoBar;
       Array<T, 3> j;
       momentTemplates<T, Descriptor>::get_rhoBar_j(cellF, oldRhoBar, j);
       T newRhoBar = D::rhoBar(neumannDensity);
       T jSqr = normSqr(j);
       for (plint iPop = 0; iPop < D::q; ++iPop) {
           T oldEq = cellF.getDynamics().computeEquilibrium(iPop, oldRhoBar, j, jSqr);
           T newEq = cellF.getDynamics().computeEquilibrium(iPop, newRhoBar, j, jSqr);
           cellF[iPop] += newEq - oldEq;
       }
   }
}

}  // namespace plb

#endif  // ELIUL_OFF_LATTICE_MODEL_3D_HH

