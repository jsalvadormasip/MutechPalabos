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
 * can be instantiated -- header file.
 */
#ifndef ZHOU_AXISYMMETRIC_DYNAMICS_H
#define ZHOU_AXISYMMETRIC_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {


/// This class implements the Axisymmetric boundary condition by Zhou (2011)
template<typename T, template<typename U> class Descriptor>
class zhouAxisymmetricDynamics : public CompositeDynamics<T,Descriptor> {
public:
    zhouAxisymmetricDynamics(Dynamics<T,Descriptor>* baseDynamics_);
    zhouAxisymmetricDynamics(HierarchicUnserializer& unserializer);
    virtual void collide(Cell<T,Descriptor>& cell, BlockStatistics& statistics_);
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);
    virtual zhouAxisymmetricDynamics<T,Descriptor>* clone() const {
        return new zhouAxisymmetricDynamics<T,Descriptor>(*this);
    }
    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);
    void setAbsoluteR(plint absoluteR_);

private:
    static int id;
    plint absoluteR;    
};

template<typename T, template<typename U> class Descriptor>
class GetAbsoluteRFunctional : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual GetAbsoluteRFunctional<T,Descriptor>* clone() const;
    virtual int getStaticId() const { return staticId; }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
private:
    static const int staticId;
};

}  // namespace plb

#endif  // ZHOU_AXISYMMETRIC_DYNAMICS_H
