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

#ifndef ELI_OFF_LATTICE_MODEL_3D_HH
#define ELI_OFF_LATTICE_MODEL_3D_HH

#ifdef _DEBUG
constexpr bool debug_mode = true;
#else
constexpr bool debug_mode = false;
#endif

#include <algorithm>
#include <cmath>
#include <vector>
#include <array>

#include "core/dynamics.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/meiLuoShyyOffLatticeModel3D.h"
#include "offLattice/nextNeighbors3D.h"
#include "ELIoffLatticeModel3D.h"

namespace plb {

/**
*
*
* @tparam T
* @tparam Descriptor
*/
template <typename T, template <typename U> class Descriptor>
ELIModels3D<T, Descriptor>::ELIModels3D(
   BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_) :
   OffLatticeModel3D<T, Array<T, 3> >(shape_, flowType_)
{ }

//template <typename T, template <typename U> class Descriptor>
//ELIModels3D<T, Descriptor> *ELIModels3D<T, Descriptor>::clone() const
//{
//   return new ELIModels3D(*this);
//}

template <typename T, template <typename U> class Descriptor>
plint ELIModels3D<T, Descriptor>::getNumNeighbors() const
{
   return 1;
}

template <typename T, template <typename U> class Descriptor>
bool ELIModels3D<T, Descriptor>::isExtrapolated() const
{
   return true;
}

template <typename T, template <typename U> class Descriptor>
void ELIModels3D<T, Descriptor>::prepareCell(
   Dot3D const &cellLocation, AtomicContainerBlock3D &container)
{
   typedef Descriptor<T> D;
   Dot3D offset = container.getLocation();
   auto *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
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

               if constexpr (debug_mode) {
                   [[maybe_unused]] bool ok = this->pointOnSurface(
                       cellLocation + offset, Dot3D(D::c[iPop][0], D::c[iPop][1], D::c[iPop][2]),
                       locatedPoint, distance, wallNormal, surfaceData, bdType, iTriangle);
                   PLB_ASSERT(ok);
               } else{
                   this->pointOnSurface(
                       cellLocation + offset, Dot3D(D::c[iPop][0], D::c[iPop][1], D::c[iPop][2]),
                       locatedPoint, distance, wallNormal, surfaceData, bdType, iTriangle);
               }
               // In the following, the importance of directions is sorted wrt. how well they
               //   are aligned with the wall normal. It is better to take the continuous normal,
               //   because it is not sensitive to the choice of the triangle when we shoot at
               //   an edge.
               // wallNormal = this->computeContinuousNormal(locatedPoint, iTriangle);
               global::timer("intersect").stop();
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
ContainerBlockData *ELIModels3D<T, Descriptor>::generateOffLatticeInfo() const
{
   return new OffLatticeInfo3D;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 3> ELIModels3D<T, Descriptor>::getLocalForce(AtomicContainerBlock3D &container) const
{
   auto *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
   PLB_ASSERT(info);
   return info->getLocalForce();
}

template <typename T, template <typename U> class Descriptor>
void ELIModels3D<T, Descriptor>::boundaryCompletion(
   AtomicBlock3D &nonTypeLattice, AtomicContainerBlock3D &container,
   std::vector<AtomicBlock3D *> const &args)
{
   auto &lattice =
       dynamic_cast<BlockLattice3D<T, Descriptor> &>(nonTypeLattice);
   auto *info = dynamic_cast<OffLatticeInfo3D *>(container.getData());
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
void ELIModels3D<T, Descriptor>::cellCompletion(
   BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
   std::vector<int> const &dryNodeFluidDirections, std::vector<plint> const &dryNodeIds,
   Dot3D const &absoluteOffset, Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args)
{
   typedef Descriptor<T> D;
   Array<T, D::d> deltaJ;
   deltaJ.resetToZero();
   Cell<T, Descriptor> &cellS = lattice.get(guoNode.x, guoNode.y, guoNode.z);
#ifdef PLB_DEBUG
   int noDynId = NoDynamics<T, Descriptor>().getId();
#endif
   PLB_ASSERT(
       cellS.getDynamics().getId() == noDynId
       && "ELIUL BC needs the dynamics to be set to NoDynamics.");
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

       Array<T, 3> wallNormal;

       auto [f_rhoBar, f_j] =
           getRhoBarJ(lattice, guoNode, args, fluidDirection, cellF);

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
//       T LambdaPlus = tauPlus-0.5;
       T tauMinus = 1./omega_minus;
//       T LambdaMinus = tauMinus-0.5;
       T fPlus  = 0.5*(cellS[i_solid] + cellFF[i_fluid]);
       T fMinus = 0.5*(cellS[i_solid] - cellFF[i_fluid]);
       T feq_fluid = cellF.getDynamics().computeEquilibrium(i_fluid, f_rhoBar, f_j, f_jSqr);
       T feq_solid = cellF.getDynamics().computeEquilibrium(i_solid, f_rhoBar, f_j, f_jSqr);
       T eqPlus = 0.5 * (feq_solid+feq_fluid);
       T eqMinus = 0.5 * (feq_solid-feq_fluid);
       T c_i_w_j = D::c[i_solid][0] * w_j[0] + D::c[i_solid][1] * w_j[1] + D::c[i_solid][2] * w_j[2];
       T eqMinusWall = -2.0 * D::t[i_solid] * D::invCs2 * c_i_w_j;

       auto [alphaPlus,alphaMinus,beta,Kplus,Kmin] = eliCoefficients(q,tauPlus,tauMinus);

       cellF[i_fluid] =   0.5 * (alphaPlus + alphaMinus) * cellS[i_solid]
                        + (1. + 0.5 * (alphaPlus - alphaMinus)) * cellF[i_fluid] /* = cellFF[i_fluid]*/
                        + beta * cellF[i_solid]
                        + Kplus * (fPlus-eqPlus) / (-1.+tauPlus)
                        + Kmin * (fMinus-eqMinus) / (-1.+tauMinus)
                        - alphaPlus * eqPlus
                        - eqMinusWall
           ;

       // the following line is needed to make the algorithm more compact and avoid fetching the
       // population in cellFF[i_fluid], but simple in cellF[i_fluid]
       cellS[i_fluid] = cellF[i_fluid];

       localForce[0] += D::c[i_solid][0] * (cellF[i_fluid] + cellS[i_solid]);
       localForce[1] += D::c[i_solid][1] * (cellF[i_fluid] + cellS[i_solid]);
       localForce[2] += D::c[i_solid][2] * (cellF[i_fluid] + cellS[i_solid]);

   }
}
template <typename T, template<typename U> class Descriptor>
std::tuple<T,Array<T,3>> ELIModels3D<T, Descriptor>::getRhoBarJ(
    const BlockLattice3D<T, Descriptor> &lattice, const Dot3D &guoNode,
    const std::vector<AtomicBlock3D *> &args, const Dot3D &fluidDirection,
    const Cell<T, Descriptor> &cellF) const
{
    T f_rhoBar = 0;
    Array<T, 3> f_j;
    if (args.empty()) {
        plbLogicError(
            "ELI requires externally "
            "provided rhoBar and j fields!");
    } else {
        if ((plint)args.size() == 1) {
            auto const *macroField =
                dynamic_cast<NTensorField3D<T> const *>(args[0]);
            PLB_ASSERT(macroField);
            // 1 Variable for rhoBar, 3 variables for j.
            PLB_ASSERT(macroField->getNdim() == 4);
            Dot3D off = computeRelativeDisplacement(lattice, *macroField);

            Dot3D macroPos = guoNode + fluidDirection + off;
            T const *macroscopic = macroField->get(macroPos.x, macroPos.y, macroPos.z);
            f_rhoBar = macroscopic[0];
            f_j.from_cArray(macroscopic + 1);
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
        } else {
            plbLogicError(
                "ELI requires externally "
                "provided rhoBar and j fields!");
        }
    }
    return std::tuple<T,Array<T,3>>{f_rhoBar, f_j};
}

template <typename T, template<typename U> class Descriptor>
inline std::array<T,5> ELIUL<T, Descriptor>::eliCoefficients(T q,T tauPlus, T tauMinus) const
{
    T alphaPlus = -1.;
    T alphaMinus = 1.;
    T Kplus = q-tauPlus;
    T Kelim = q-tauMinus;
    T Kmin = Kelim;
    T beta = 0;
    return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
}

template <typename T, template <typename U> class D>
ELIUL<T, D> *ELIUL<T, D>::clone() const
{
    return new ELIUL<T,D>(*this);
}

template <typename T, template <typename U> class D>
ELIULC<T, D> *ELIULC<T, D>::clone() const
{
    return new ELIULC<T,D>(*this);
}
template <typename T, template <typename U> class D>
std::array<T, 5> ELIULC<T, D>::eliCoefficients(T q, T tauPlus, T tauMinus) const
{
    T alphaPlus = -1.;
    T alphaMinus = 1.;
    T Kplus = q - tauPlus;
//    T LambdaMinus = tauMinus-0.5;
    T Kmin = 0.0;
    T beta = 0.0;
    return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
}

template <typename T, template <typename U> class D>
ELIULK1<T, D> *ELIULK1<T, D>::clone() const
{
    return new ELIULK1<T,D>(*this);
}
template <typename T, template <typename U> class D>
std::array<T, 5> ELIULK1<T, D>::eliCoefficients(T q, T tauPlus, T tauMinus) const
{
    T alphaPlus = -1.;
    T alphaMinus = 1.;
    T Kplus = q - tauPlus;
    //    T Kelim = q - tauMinus;
    T LambdaMinus = tauMinus-0.5;
    T Kmin = 1. + alphaMinus * (LambdaMinus - 0.5);
    T beta = 0.0;
    return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
}

template <typename T, template <typename U> class D>
ELIULK4<T, D> *ELIULK4<T, D>::clone() const
{
    return new ELIULK4<T, D>(*this);
}

template <typename T, template <typename U> class D>
std::array<T, 5> ELIULK4<T, D>::eliCoefficients(T q, T tauPlus, T tauMinus) const
{
    T alphaPlus = -1.;
    T alphaMinus = 1.;
    T Kplus = q - tauPlus;
//    T Kelim = q - tauMinus;
    T LambdaMinus = tauMinus-0.5;
    T Kmin = 1. + alphaMinus * (LambdaMinus - 0.5);
    T beta = 0.0;
    return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
}


template <typename T, template <typename U> class D>
ELIULK3<T, D> *ELIULK3<T, D>::clone() const
{
    return new ELIULK3<T,D>(*this);
}


template <typename T, template <typename U> class D>
std::array<T, 5> ELIULK3<T, D>::eliCoefficients(T q, T tauPlus, T tauMinus) const
{
    T alphaPlus = -1.;
    T alphaMinus = 1.;
    T Kplus = q - tauPlus;
    //    T Kelim = q - tauMinus;
    T LambdaMinus = tauMinus-0.5;
    T LambdaPlus = tauPlus-0.5;
    T Kmin = 1. + alphaMinus * LambdaMinus - (alphaMinus*(q*q+LambdaPlus))/(2.*LambdaPlus);
    T beta = 0.0;
    return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
}

template <typename T, template <typename U> class D, typename Function>
#if ((defined(_MSVC_LANG) && _MSVC_LANG > 201703L) || __cplusplus > 201703L)
    requires std::invocable<Function, T, T, T> /*&& object<Function>*/
#endif
ELIgeneric<T, D, Function>::ELIgeneric(
    BoundaryShape3D<T, Array<T, 3>> *shape_, int flowType_, Function coefficients_):
    ELIModels3D<T, D>(shape_,flowType_),
    compute_coefficients(coefficients_)
{
//    compute_coefficients = coefficients_;
}


template <typename T, template <typename U> class D, typename Function>
#if ((defined(_MSVC_LANG) && _MSVC_LANG > 201703L) || __cplusplus > 201703L)
    requires std::invocable<Function, T, T, T> /*&& object<Function>*/
#endif
ELIgeneric<T, D, Function> *ELIgeneric<T, D, Function>::clone() const
{
    return new ELIgeneric<T, D, Function>(*this);
}

template <typename T, template <typename U> class D, typename Function>
#if ((defined(_MSVC_LANG) && _MSVC_LANG > 201703L) || __cplusplus > 201703L)
    requires std::invocable<Function, T, T, T> /*&& object<Function>*/
#endif
std::array<T, 5> ELIgeneric<T, D, Function>::eliCoefficients(T q, T tauPlus, T tauMinus) const
{
    return compute_coefficients(q,tauPlus,tauMinus);
}

}  // namespace plb

#endif  // ELI_OFF_LATTICE_MODEL_3D_HH

