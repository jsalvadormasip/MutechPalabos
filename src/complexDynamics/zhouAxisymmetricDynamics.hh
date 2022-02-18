/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
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
 * Axisymmetric LBM dynamics, refering to 
 * Zhou, J. G., 2011, “Axisymmetric Lattice Boltzmann Method Revised,” Physical review E, 84(3), p. 036704.
 * can be instantiated -- generic implementation.
 */
#ifndef ZHOU_AXISYMMETRIC_DYNAMICS_HH
#define ZHOU_AXISYMMETRIC_DYNAMICS_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "core/hierarchicSerializer.h"

namespace plb {

// ============== Zhou (2011) Axisymmetric dynamics ====================== //

template<typename T, template<typename U> class Descriptor>
int ZhouAxisymmetricDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,ZhouAxisymmetricDynamics<T,Descriptor> >("Zhou_Axisymmetric_Dynamics");

template<typename T, template<typename U> class Descriptor>
ZhouAxisymmetricDynamics<T,Descriptor>::ZhouAxisymmetricDynamics(Dynamics<T,Descriptor>* baseDynamics_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_, false)  // false is for automaticPrepareCollision.
{
    PLB_ASSERT(Descriptor<T>::d == 2);
}

template<typename T, template<typename U> class Descriptor>
ZhouAxisymmetricDynamics<T,Descriptor>::ZhouAxisymmetricDynamics(HierarchicUnserializer& unserializer)
    : CompositeDynamics<T,Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell, BlockStatistics& statistics)
{

    Array<T,Descriptor<T>::numPop> h1;         /// h1
    Array<T,Descriptor<T>::numPop> h2;         /// h2
    h1.resetToZero();
    h2.resetToZero();

    T rhoBar;
    Array<T,Descriptor<T>::d> j;

    if(absoluteR != 0) { // h1 & h2 have a non-zero value only if r != 0
        Array<T,Descriptor<T>::q>& f = cell.getRawPopulations();
        momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
        T rho = Descriptor<T>::fullRho(rhoBar);
        T invRho = Descriptor<T>::invRho(rhoBar);
        const T jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);
        Array<T,Descriptor<T>::d> u;
        u = j*invRho; // u(0) = u_x and u(1) = u_r

        for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) {
            h1[iPop] = -Descriptor<T>::t[iPop] * j[1] / (T) absoluteR; // Zhou 2011
        }

        Array<T,Descriptor<T>::q> fNeq;
        T tau = 1.0/this->getOmega();
        T nu = (2.0 * tau - 1.0) / 6.0;
        for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) {
            fNeq[iPop] = f[iPop] - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (
                    iPop, rhoBar, invRho, j, jSqr );
            h2[iPop] = - (2*tau-1) / (2*tau*(T)absoluteR) * Descriptor<T>::c[iPop][1] * fNeq[iPop]
                    + (T) 1.0 / (T) 6.0 * rho * (
                            Descriptor<T>::c[iPop][0] * (-u[0]*u[1]/(T) absoluteR) +
                            Descriptor<T>::c[iPop][1] * (-u[1]*u[1]/(T) absoluteR - 2*nu*u[1]/util::sqr((T)absoluteR))
                            );
        }
    }
    else { // Specular boundary condition for r = 0;
        Array<int,Descriptor<T>::q> reflection; reflection[0] = -1; reflection[1] = 1;
        for (plint iPop = 3; iPop <= 5; iPop++) { // populations at southwest, south, southeast
            Array<int,Descriptor<T>::d> v; 
            v[0] = -Descriptor<T>::c[indexTemplates::opposite<Descriptor<T> >(iPop)][0]; 
            v[1] = Descriptor<T>::c[indexTemplates::opposite<Descriptor<T> >(iPop)][1];
            plint jPop = indexTemplates::findVelocity<Descriptor<T> >(v);
            PLB_ASSERT(jPop != Descriptor<T>::q);
            cell[jPop] = cell[iPop];
        }
    }

    this->getBaseDynamics().collide(cell, statistics);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] += h1[iPop] + h2[iPop];
    }

    momentTemplates<T,Descriptor >::get_rhoBar_j(cell, rhoBar, j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    T uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j)*invRho*invRho;
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& statistics )
{

    Array<T,Descriptor<T>::numPop> h1;         /// h1
    Array<T,Descriptor<T>::numPop> h2;         /// h2
    h1.resetToZero();
    h2.resetToZero();

    if(absoluteR != 0) { // h1 & h2 have a non-zero value only if r != 0
        Array<T,Descriptor<T>::q>& f = cell.getRawPopulations();
        T rho = Descriptor<T>::fullRho(rhoBar);
        T invRho = Descriptor<T>::invRho(rhoBar);
        const T jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);
        Array<T,Descriptor<T>::d> u;
        u = j*invRho; // u(0) = u_x and u(1) = u_r

        for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) {
            h1[iPop] = -Descriptor<T>::t[iPop] * j[1] / (T) absoluteR; // Zhou 2011
        }

        Array<T,Descriptor<T>::q> fNeq;
        T tau = 1.0/this->getOmega();
        T nu = (2.0 * tau - 1.0) / 6.0;
        for (plint iPop=0; iPop < Descriptor<T>::q; ++iPop) {
            fNeq[iPop] = f[iPop] - dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium (
                    iPop, rhoBar, invRho, j, jSqr );
            h2[iPop] = - (2*tau-1) / (2*tau*(T)absoluteR) * Descriptor<T>::c[iPop][1] * fNeq[iPop]
                    + (T) 1.0 / (T) 6.0 * rho * (
                            Descriptor<T>::c[iPop][0] * (-u[0]*u[1]/(T) absoluteR) +
                            Descriptor<T>::c[iPop][1] * (-u[1]*u[1]/(T) absoluteR - 2*nu*u[1]/util::sqr((T)absoluteR))
                            );
        }
    }
    else { // Specular boundary condition for r = 0;
        Array<int,Descriptor<T>::q> reflection; reflection[0] = -1; reflection[1] = 1;
        for (plint iPop = 3; iPop <= 5; iPop++) {
            Array<int,Descriptor<T>::d> v; 
            v[0] = -Descriptor<T>::c[indexTemplates::opposite<Descriptor<T> >(iPop)][0]; 
            v[1] = Descriptor<T>::c[indexTemplates::opposite<Descriptor<T> >(iPop)][1];
            plint jPop = indexTemplates::findVelocity<Descriptor<T> >(v);
            PLB_ASSERT(jPop != Descriptor<T>::q);
            cell[jPop] = cell[iPop];
        }
    }

    this->getBaseDynamics().collideExternal(cell, rhoBar, j, thetaBar, statistics);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] += h1[iPop] + h2[iPop];
    }
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
{ }

template<typename T, template<typename U> class Descriptor>
int ZhouAxisymmetricDynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    CompositeDynamics<T,Descriptor>::serialize(serializer);
    serializer.addValue<plint>(absoluteR);
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue<plint>(absoluteR);
}

template<typename T, template<typename U> class Descriptor>
void ZhouAxisymmetricDynamics<T,Descriptor>::setAbsoluteR(plint absoluteR_) {
    absoluteR = absoluteR_;
}

/********************************* GetAbsoluteRFunctional ******************/
template<typename T, template<typename U> class Descriptor>
const int GetAbsoluteRFunctional<T,Descriptor>::staticId =
        meta::registerProcessor2D < GetAbsoluteRFunctional<T, Descriptor>,
                T, Descriptor> (std::string("GetAbsoluteRFunctional"));

template<typename T, template<typename U> class Descriptor>
void GetAbsoluteRFunctional<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    Dot2D absoluteOffset = lattice.getLocation();
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            ZhouAxisymmetricDynamics<T,Descriptor>& dynamics = dynamic_cast<ZhouAxisymmetricDynamics<T,Descriptor>&>(lattice.get(iX,iY).getDynamics());
            dynamics.setAbsoluteR(iY + absoluteOffset.y);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
GetAbsoluteRFunctional<T,Descriptor>*
GetAbsoluteRFunctional<T,Descriptor>::clone() const
{
    return new GetAbsoluteRFunctional<T,Descriptor>(*this);
}


}  // namespace plb

#endif  // ZHOU_AXISYMMETRIC_DYNAMICS_HH

