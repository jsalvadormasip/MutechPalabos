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

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/dynamics.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/meiLuoShyyOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"

namespace plb {

/**
* This class implements the Mei-Luo-Shyy (MLS,1999) boundary condition on a BoundaryShape.
* The BoundaryShape determines whether the points of the discrete lattice are voxelFlag::inside
* or voxelFlag::outside some geometry.
*
* It can handle moving boundaries using the momentum correction of ladd (LADD, 1994).
* The wall velocity is recovered from SurfaceData stored in BoundaryShape3D<T,SurfaceData>*
*
* NOTE: this class was previously called FilippovaHaenelModel3D in Palabos before June 2020
*
* (MLS,1999) R. Mei, L.-S. Luo, and W. Shyy, “An Accurate Curved Boundary Treatment in the Lattice
* Boltzmann Method,” Journal of Computational Physics, vol. 155, no. 2, pp. 307–330, Nov. 1999,
* doi: 10.1006/jcph.1999.6334.
*
* (LADD, 1994) A. J. C. Ladd, “Numerical simulations of particulate suspensions via a discretized
* Boltzmann equation. Part 1. Theoretical foundation,” Journal of Fluid Mechanics, vol. 271, pp.
* 285–309, Jul. 1994, doi: 10.1017/S0022112094001771.
*
* @tparam T
* @tparam Descriptor
*/
template <typename T, template <typename U> class Descriptor>
ELIULModel3D<T, Descriptor>::ELIULModel3D(
   BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
   OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

template <typename T, template <typename U> class Descriptor>
ELIULModel3D<T, Descriptor> *ELIULModel3D<T, Descriptor>::clone() const
{
   return new ELIULModel3D(*this);
}

template <typename T, template <typename U> class Descriptor>
plint ELIULModel3D<T, Descriptor>::getNumNeighbors() const
{
   return 1;
}

template <typename T, template <typename U> class Descriptor>
bool ELIULModel3D<T, Descriptor>::isExtrapolated() const
{
   return true;
}

template <typename T, template <typename U> class Descriptor>
void ELIULModel3D<T, Descriptor>::prepareCell(
   Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
   typedef Descriptor<T> D;
   Dot3D offset = container.getLocation();
   OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   std::vector<int> liquidNeighbors;
   std::vector<plint> ids;
   Dot3D absLoc = cellLocation + offset;
   if (this->isSolid(absLoc)) {
       for (int iPop = 0; iPop < D::q; ++iPop) {
           Dot3D neighbor(
               cellLocation.x + D::c[iPop][0], cellLocation.y + D::c[iPop][1],
               cellLocation.z + D::c[iPop][2]);
           Dot3D neighborLoc = neighbor + offset;
           // If the non-fluid node has a fluid neighbor ...
           if (this->isFluid(neighborLoc)) {
               // ... check how many fluid nodes it has ahead of it ...
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
               // wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
               global::timer("intersect").stop();
               PLB_ASSERT(ok);
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

template <typename T, template <typename U> class Descriptor>
ContainerBlockData *ELIULModel3D<T, Descriptor>::generateOffLatticeInfo() const
{
   return new OffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> ELIULModel3D<T, Descriptor>::getLocalForce(AtomicContainerBlock3D &container) const
{
   OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void ELIULModel3D<T, Descriptor>::boundaryCompletion(
   AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
   std::vector<AtomicBlock3D *> const &args)
{
   BlockLattice3D<T, Descriptor> &lattice =
       dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
   OffLatticeInfo3D *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   std::vector<Dot3D> const &dryNodes = info->getDryNodes();
   std::vector<std::vector<int> > const &dryNodeFluidDirections =
       info->getDryNodeFluidDirections();
   std::vector<std::vector<plint> > const &dryNodeIds = info->getDryNodeIds();
   PLB_ASSERT(dryNodes.size() == dryNodeFluidDirections.size());

   Dot3D absoluteOffset = container.getLocation();

   Array<T, 3> &localForce = info->getLocalForce();
   localForce.resetToZero();
   for (pluint iDry = 0; iDry < dryNodes.size(); ++iDry) {
       cellCompletion(
           lattice, dryNodes[iDry], dryNodeFluidDirections[iDry], dryNodeIds[iDry], absoluteOffset,
           localForce, args);
   }
}

template <typename T, template <typename U> class Descriptor>
void ELIULModel3D<T, Descriptor>::cellCompletion(
   BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
   std::vector<int> const &dryNodeFluidDirections, std::vector<plint> const &dryNodeIds,
   Dot3D const &absoluteOffset, Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args)
{
   typedef Descriptor<T> D;
   Array<T, D::d> deltaJ;
   deltaJ.resetToZero();
   Cell<T, Descriptor> &cellS = lattice.get(guoNode.x, guoNode.y, guoNode.z);
//#ifdef PLB_DEBUG
//   int noDynId = NoDynamics<T, Descriptor>().getId();
//#endif
//   PLB_ASSERT(
//       cellS.getDynamics().getId() == noDynId
//       && "ELIUL BC needs the dynamics to be set to NoDynamics.");
   for (plint iDirection = 0; iDirection < (plint)dryNodeFluidDirections.size(); ++iDirection) {
       int i_fluid = dryNodeFluidDirections[iDirection];
       int i_solid = indexTemplates::opposite<Descriptor<T> >(i_fluid);
       Dot3D fluidDirection(D::c[i_fluid][0], D::c[i_fluid][1], D::c[i_fluid][2]);
       plint dryNodeId = dryNodeIds[iDirection];

       Array<T, 3> wallNode, wall_vel;
       T wallDistance;
       OffBoundary::Type bdType;
       Cell<T, Descriptor> &cellF = lattice.get(
           guoNode.x + fluidDirection.x, guoNode.y + fluidDirection.y,
           guoNode.z + fluidDirection.z);
       Cell<T, Descriptor> const &cellFF = lattice.get(
           guoNode.x + 2*fluidDirection.x, guoNode.y + 2*fluidDirection.y,
           guoNode.z + 2*fluidDirection.z);

       Cell<T, Descriptor> collidedCell(cellF);
       BlockStatistics statsCopy(lattice.getInternalStatistics());
       collidedCell.collide(statsCopy);

       T f_rhoBar, ff_rhoBar;
       Array<T, 3> f_j, ff_j;
       Array<T, 3> wallNormal;

       if (args.empty()) {
           cellF.getDynamics().computeRhoBarJ(cellF, f_rhoBar, f_j);
           cellFF.getDynamics().computeRhoBarJ(cellFF, ff_rhoBar, ff_j);
       } else {
           if ((plint)args.size() == 1) {
               auto const *macroField =
                   dynamic_cast<NTensorField3D<T> const *>(args[0]);
               PLB_ASSERT(macroField);
               // 1 Variable for rhoBar, 3 variables for j.
               PLB_ASSERT(macroField->getNdim() == 4);
               Dot3D off = computeRelativeDisplacement(lattice, *macroField);

               Dot3D macroPos = guoNode + fluidDirection + off;
               Dot3D macroPos2 = macroPos + fluidDirection;
               T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
               f_rhoBar = macroscopic[0];
               f_j.from_cArray(macroscopic + 1);
               T const *macroscopic2 = macroField->get(macroPos2.x, macroPos2.y, macroPos2.z);
               ff_j.from_cArray(macroscopic2 + 1);
           } else if ((plint)args.size() == 2) {
               // 1 field for rhoBar, 1 field for j.
               // #ifdef PLB_DEBUG
               auto const *rhoBarField =
                   dynamic_cast<ScalarField3D<T> const *>(args[0]);
               // #endif
               auto const *jField =
                   dynamic_cast<TensorField3D<T, 3> const *>(args[1]);
               PLB_ASSERT(rhoBarField);
               PLB_ASSERT(jField);

               Dot3D posJ =
                   guoNode + fluidDirection + computeRelativeDisplacement(lattice, *jField);
               Dot3D posRho =
                   guoNode + fluidDirection + computeRelativeDisplacement(lattice, *rhoBarField);

               f_rhoBar = rhoBarField->get(posRho.x, posRho.y, posRho.z);
               f_j = jField->get(posJ.x, posJ.y, posJ.z);
               Dot3D posJ2 = posJ + fluidDirection;

               ff_j = jField->get(posJ2.x, posJ2.y, posJ2.z);
           } else {
               PLB_ASSERT(false);  // Not implemented for 3 args.
           }
       }

       T f_rho = D::fullRho(f_rhoBar);
       T f_jSqr = normSqr(f_j);

#ifdef PLB_DEBUG
       bool ok =
#endif
           this->pointOnSurface(
               guoNode + absoluteOffset, fluidDirection, wallNode, wallDistance, wallNormal,
               wall_vel, bdType, dryNodeId);
       PLB_ASSERT(ok);

       Array<T, 3> w_j = wall_vel * f_rho;
       T d = std::sqrt(D::cNormSqr[i_fluid]);
       PLB_ASSERT(wallDistance <= d);
       T delta = 1.0 - wallDistance / d;
       T& q = delta;
       PLB_ASSERT(q <= 1);
       PLB_ASSERT(q >= 0);


       const auto trt = dynamic_cast<const BaseTRTdynamics<T, Descriptor> *>(&cellF.getDynamics());
       T omega_plus  = cellF.getDynamics().getOmega();
       T omega_minus = 0.;
       if(trt)
           omega_minus = cellF.getDynamics().getParameter(dynamicParams::omega_minus);
       else
           omega_minus = omega_plus;
       T tauPlus = 1./omega_plus;
       T LambdaPlus = tauPlus-0.5;
       T tauMinus = 1./omega_minus;
       T LambdaMinus = tauMinus-0.5;
       T fPlus  = 0.5*(cellS[i_solid] + cellFF[i_fluid]);
       T fMinus = 0.5*(cellS[i_solid] - cellFF[i_fluid]);
       T Kelip = q-tauPlus;
       T Kelim = q-tauMinus;
       T feq_fluid = cellF.getDynamics().computeEquilibrium(i_fluid, f_rhoBar, f_j, f_jSqr);
       T feq_solid = cellF.getDynamics().computeEquilibrium(i_solid, f_rhoBar, f_j, f_jSqr);
       T eqPlus = 0.5 * (feq_solid+feq_fluid);
       T eqMinus = 0.5 * (feq_solid-feq_fluid);


       // YLI
//       T alphaPlus = 0.;
//       T alphaMinus = 2./(1.+q);
//       T beta = (1.-q)/(1.+q);
//       Kelip = 0.;
//       Kelim = 0.;
//       eqPlus = 0.;


       // ELIUL
       T alphaPlus = -1.;
       T alphaMinus = 1.;
       T beta = 0.;
       T Kmin = Kelim;

       // ELIUL-C
       Kmin = 0.0;
//       // k1 eli
//       Kmin = 1.-alphaMinus/2.0;
//       // k4 eli
//       Kmin = 1.+alphaMinus*(LambdaMinus-0.5);

       cellS[i_fluid] =   0.5 * (alphaPlus + alphaMinus) * cellS[i_solid]
                        + (1. + 0.5 * (alphaPlus - alphaMinus)) * cellF[i_fluid]
                        + beta * cellFF[i_fluid]
                        + Kelip * (fPlus-eqPlus) / (-1.+tauPlus)
                        + Kmin * (fMinus-eqMinus) / (-1.+tauMinus)
                        - alphaPlus * eqPlus
           //                                - eqMinusWall
           ;

       localForce[0] += D::c[i_solid][0] * (cellS[i_fluid] + cellS[i_solid]);
       localForce[1] += D::c[i_solid][1] * (cellS[i_fluid] + cellS[i_solid]);
       localForce[2] += D::c[i_solid][2] * (cellS[i_fluid] + cellS[i_solid]);

   }
   for (plint iDirection = 0; iDirection < (plint)dryNodeFluidDirections.size(); ++iDirection){
       int i_fluid = dryNodeFluidDirections[iDirection];
       int i_solid = indexTemplates::opposite<Descriptor<T> >(i_fluid);
       Dot3D fluidDirection(D::c[i_fluid][0], D::c[i_fluid][1], D::c[i_fluid][2]);
       plint dryNodeId = dryNodeIds[iDirection];
       Cell<T, Descriptor> &cellF = lattice.get(
           guoNode.x + fluidDirection.x, guoNode.y + fluidDirection.y,
           guoNode.z + fluidDirection.z);
       Cell<T, Descriptor> const &cellFF = lattice.get(
           guoNode.x + 2*fluidDirection.x, guoNode.y + 2*fluidDirection.y,
           guoNode.z + 2*fluidDirection.z);
       cellF[i_fluid] = cellS[i_fluid];
   }
}

}  // namespace plb

#endif  // ELIUL_OFF_LATTICE_MODEL_3D_HH

