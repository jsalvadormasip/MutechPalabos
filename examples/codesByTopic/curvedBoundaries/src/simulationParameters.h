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

#ifndef CURVED_BOUNDARIES_SIMULATIONPARAMETERS_H
#define CURVED_BOUNDARIES_SIMULATIONPARAMETERS_H

#undef NDEBUG
#include <cassert>  // reinclude the header to update the definition of assert()

namespace incompressible_simulation_parameters {
template <typename Real>
class NonDimensional {
public:
    NonDimensional(Real l_ref, Real rho_f);

    NonDimensional();

    [[maybe_unused]] NonDimensional<Real> initReLxLyLz(Real re_, Real lx_,
                                                         Real ly_, Real lz_);

    [[maybe_unused]] NonDimensional<Real> initReMaLxLyLz(Real re_, Real ma_,
                                                            Real lx_, Real ly_,
                                                            Real lz_);

    [[maybe_unused]] NonDimensional<Real> initReLxLyLzRhos(
        Real re_, Real lx_, Real ly_, Real lz_, Real rho_solid_);

public:
    Real getRe() const {
        assert(re > 0);
        return re;
    }

    Real getLx() const {
        assert(re > 0);
        return lx;
    }

    Real getLy() const {
        assert(re > 0);
        return ly;
    }

    Real getLz() const {
        assert(re > 0);
        return lz;
    }

    Real getMa() const {
        if (ma_is_set)
            return ma;
        else {
            pcout << "ERROR! Ma was not set!!" << std::endl;
            abort();
        }
    }

    Real getRhof() const { return rho_f; }

    Real getRhoSolid() const {
        if (rho_is_set)
            return rho_solid;
        else
            abort();
    }

    NonDimensional<Real> printParameters() const;

    [[nodiscard]] bool initialized() const { return is_initialized; }

private:
    Real re;
    Real ma;
    bool ma_is_set;
    Real rho_solid;  // nb: in rho_fluid units i.e. rho_solid/rho_fluid
    bool rho_is_set;
    Real ly;  // nb: in reference length units
    Real lx;
    Real lz;
    const Real l_ref;
    const Real rho_f;
    bool is_initialized;
};

/**
 * This class handles all the simulation parameters in lattice units. It allows
 * to change the parameters only by defining them in block, by means of init_
 * functions. It does not store dimensionless number and parameters, because
 * everything is consistently defined in lattice units. After initializing, you
 * can get the values of the attributes using the getters, that implement simple
 * assert conditions to avoid units errors. After initialization
 * @tparam Real Template type for real numbers
 * @tparam Int Template type for integer numbers
 */
template <typename Real, typename Int>
class Numerics {
public:
    Numerics();
    ~Numerics() = default;

    [[maybe_unused]] Numerics<Real, Int>& initLU(Int l_ref_lu_, Real tau_lu_,
                                                 Real u_lb_, Int lx_lu,
                                                 Int ly_lu, Int lz_lu = 0);

    [[maybe_unused]] Numerics<Real, Int>& initReTauUlb(Real re, Real tau_lu_,
                                                       Real u_lb_, Int lx_lu,
                                                       Int ly_lu,
                                                       Int lz_lu = 0);

    [[maybe_unused]] Numerics<Real, Int>& initLrefReUlb(Int l_ref_lu_, Real re,
                                                        Real u_lb_, Int lx_lu,
                                                        Int ly_lu,
                                                        Int lz_lu = 0);

    [[maybe_unused]] Numerics<Real, Int>& initLrefReTau(Int l_ref_lu_, Real re,
                                                        Real tau_, Int lx_lu,
                                                        Int ly_lu,
                                                        Int lz_lu = 0);

    [[maybe_unused]] Numerics<Real, Int>& initLrefluNodim(
        Int l_ref_lu_, NonDimensional<Real>* dimless_, Real u_lb_ = -1,
        Real tau_ = -1);

public:
    Real getUlb() const {
        assert(u_lb >= 0);
        return u_lb;
    }

    Real getLref() const {
        assert(l_ref_lu >= 0);
        return l_ref_lu;
    }

    Real getTau() const {
        assert(tau >= 0);
        return tau;
    }

    Real getOmega() const {
        assert(tau >= 0);
        return 1. / tau;
    }

    Real getNulu(Real tau_, Real cs2) const {
        assert(tau >= 0);
        return cs2 * (tau - 1. / 2.);
    }

    Real getNulu(Real tau_) const {
        assert(tau_ >= 0);
        return (tau_ - 1. / 2.) / 3.0;
    }

    Real getLx() const {
        assert(lx_domain >= 0);
        return lx_domain;
    }

    Real getLy() const {
        assert(ly_domain >= 0);
        return ly_domain;
    }

    Real getLz() const {
        assert(lz_domain >= 0);
        return lz_domain;
    }

    Real getCs2() const { return cs2; }

    Real getCs2Inv() const { return 1. / cs2; }

    Real getCs() const { return sqrt(cs2); }

    [[nodiscard]] bool initialized() const { return is_initialized; }

    Numerics<Real, Int>& printParameters();

    Numerics<Real, Int>& writeParametersLog() const;

    plb::IncomprFlowParam<Real> getIncomprFlowParam() const;

private:
    void setTau(Real tau_) {
        assert(tau_ < 3 && tau_ > 0.5);
        tau = tau_;
    }

    void setUlb(Real u_lb_) {
        assert(u_lb_ < 0.4);
        u_lb = u_lb_;
    }

private:
    Real tau;
    Real u_lb;
    Int l_ref_lu;
    Int lx_domain;
    Int ly_domain;
    Int lz_domain;
    Int max_iter;
    Real cs2;
    Real rhof0;
    bool is_initialized;
    NonDimensional<Real>* dimless;
};

}  // namespace incompressible_simulation_parameters
#endif  // CURVED_BOUNDARIES_SIMULATIONPARAMETERS_H
