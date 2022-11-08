/* This code is a showcase for the Palabos library.
 *
 * The Palabos software is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * You can redistribute it and/or modify this code
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#ifndef CURVED_BOUNDARIES_SETUP_POISEUILLE_H
#define CURVED_BOUNDARIES_SETUP_POISEUILLE_H
using namespace plb;
using namespace plb::descriptors;
/// Velocity on the parabolic Poiseuille profile
template <typename Real>
Real poiseuilleVelocity(plint iY, IncomprFlowParam<Real> const& parameters) {
    Real y = (Real)iY / (Real)parameters.getNy();
    return 4. * parameters.getLatticeU() * (y - y * y);
}

/// Linearly decreasing pressure profile
template <typename Real>
Real poiseuillePressure(plint iX, IncomprFlowParam<Real> const& parameters) {
    Real Lx = parameters.getNx() - 1;
    Real Ly = parameters.getNy() - 1;
    return 8. * parameters.getLatticeNu() * parameters.getLatticeU() /
           (Ly * Ly) * (Lx / (Real)2 - (Real)iX);
}

/// Convert pressure to density according to ideal gas law
template <typename Real, template <typename U> class Descriptor>
Real poiseuilleDensity(plint iX, IncomprFlowParam<Real> const& parameters) {
    return poiseuillePressure(iX, parameters) * Descriptor<Real>::invCs2 +
           (Real)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template <typename Real>
class PoiseuilleVelocity {
public:
    explicit PoiseuilleVelocity(IncomprFlowParam<Real> parameters_)
        : parameters(parameters_) {}
    void operator()(plint iX, plint iY, Array<Real, 2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = Real();
    }
    void operator()(plint iX, plint iY, plint iZ, Array<Real, 3>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = Real();
        u[2] = Real();
    }

private:
    IncomprFlowParam<Real> parameters;
};

/// A functional, used to initialize the velocity for the boundary conditions
template <typename Real>
class ConstantVelocity {
public:
    explicit ConstantVelocity(IncomprFlowParam<Real> parameters_)
        : parameters(parameters_) {}
    void operator()(plint iX, plint iY, Array<Real, 2>& u) const {
        u[0] = parameters.getLatticeU();
        u[1] = Real();
    }
    void operator()(plint iX, plint iY, plint iZ, Array<Real, 3>& u) const {
        u[0] = parameters.getLatticeU();
        u[1] = Real();
        u[2] = Real();
    }

private:
    IncomprFlowParam<Real> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template <typename T>
class ConstantDensity {
public:
    explicit ConstantDensity(T density_) : density(density_) {}
    T operator()(plint iX, plint iY) const { return density; }
    T operator()(plint iX, plint iY, plint iZ) const { return density; }

private:
    T density;
};

/// A functional, used to create an initial condition for the density and
/// velocity
template <typename Real, template <typename U> class Descriptor>
class PoiseuilleVelocityAndDensity {
public:
    explicit PoiseuilleVelocityAndDensity(IncomprFlowParam<Real> parameters_)
        : parameters(parameters_) {}
    void operator()(plint iX, plint iY, Real& rho, Array<Real, 2>& u) const {
        rho = poiseuilleDensity<Real, Descriptor>(iX, parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = Real();
    }
    void operator()(plint iX, plint iY, plint iZ, Real& rho,
                    Array<Real, 3>& u) const {
        rho = poiseuilleDensity<Real, Descriptor>(iX, parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = Real();
        u[2] = Real();
    }

private:
    IncomprFlowParam<Real> parameters;
};

/// A functional, used to create an initial condition for the density and
/// velocity
template <typename Real, template <typename U> class Descriptor>
class ConstantVelocityAndDensity {
public:
    explicit ConstantVelocityAndDensity(IncomprFlowParam<Real> parameters_)
        : parameters(parameters_) {}
    void operator()(plint iX, plint iY, Real& rho, Array<Real, 2>& u) const {
        rho = 1.0;
        u[0] = parameters.getLatticeU();
        u[1] = Real();
    }
    void operator()(plint iX, plint iY, plint iZ, Real& rho,
                    Array<Real, 3>& u) const {
        rho = 1.0;
        u[0] = parameters.getLatticeU();
        u[1] = Real();
        u[2] = Real();
    }

private:
    IncomprFlowParam<Real> parameters;
};

/// A functional, used to create an initial condition for the density and
/// velocity
template <typename Real, template <typename U> class Descriptor>
class PoiseuilleDensity {
public:
    explicit PoiseuilleDensity(IncomprFlowParam<Real> parameters_)
        : parameters(parameters_) {}
    Real operator()(plint iX, plint iY) const {
        return poiseuilleDensity<Real, Descriptor>(iX, parameters);
    }
    Real operator()(plint iX, plint iY, plint iZ) const {
        return poiseuilleDensity<Real, Descriptor>(iX, parameters);
    }

private:
    IncomprFlowParam<Real> parameters;
};

#endif  // CURVED_BOUNDARIES_SETUP_POISEUILLE_H
