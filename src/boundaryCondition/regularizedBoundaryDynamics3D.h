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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_3D_H
#define REGULARIZED_BOUNDARY_DYNAMICS_3D_H

#include "boundaryCondition/regularizedBoundaryDynamics.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
class RegularizedVelocityInnerEdgeDynamics3D : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    enum { direction1 = (plane + 1) % 3, direction2 = (plane + 2) % 3 };

public:
    /* *************** Construction and Destruction ********************** */

    RegularizedVelocityInnerEdgeDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);

    /// Clone the object, based on its dynamic type
    virtual RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2> *clone()
        const;

    virtual void replaceBaseDynamics(Dynamics<T, Descriptor> *newBaseDynamics);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u_) const;

    /// Define velocity. Stores value inside Dynamics object.
    virtual void defineVelocity(Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u_);

    /// Compute density from incoming particle populations
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /* *************** Other virtual methods ***************************** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Default completion scheme, does nothing
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    RegularizedVelocityBoundaryDynamics<T, Descriptor, direction1, normal1> dynamics1;
    RegularizedVelocityBoundaryDynamics<T, Descriptor, direction2, normal2> dynamics2;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class RegularizedVelocityInnerCornerDynamics3D : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /* *************** Construction and Destruction ********************** */

    RegularizedVelocityInnerCornerDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics_, bool automaticPrepareCollision_ = true);

    /// Clone the object, based on its dynamic type
    virtual RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, normalX, normalY, normalZ>
        *clone() const;

    virtual void replaceBaseDynamics(Dynamics<T, Descriptor> *newBaseDynamics);

    /* *************** Computation of macroscopic variables ************** */

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity(
        Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u_) const;

    /// Define velocity. Stores value inside Dynamics object.
    virtual void defineVelocity(Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u_);

    /// Compute density from incoming particle populations
    virtual T computeDensity(Cell<T, Descriptor> const &cell) const;

    /* *************** Other virtual methods ***************************** */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T, Descriptor> const &cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar_, Array<T, Descriptor<T>::d> &j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
        Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const;

    /// Default completion scheme, does nothing
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    RegularizedVelocityBoundaryDynamics<T, Descriptor, 0, normalX> xDynamics;
    RegularizedVelocityBoundaryDynamics<T, Descriptor, 1, normalY> yDynamics;
    RegularizedVelocityBoundaryDynamics<T, Descriptor, 2, normalZ> zDynamics;
};

}  // namespace plb

#endif  // REGULARIZED_BOUNDARY_DYNAMICS_3D_H
