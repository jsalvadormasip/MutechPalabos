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

#include "simulationParameters.h"
#undef NDEBUG
#include <cassert>  // reinclude the header to update the definition of assert()
namespace incompressible_simulation_parameters {
template <typename Real>
NonDimensional<Real>::NonDimensional(Real l_ref, Real rho_f) : l_ref(l_ref), rho_f(rho_f)
{
    re = 0;
    ma = -1;
    ma_is_set = false;
    rho_solid = -1;  // nb: in rho_fluid units i.e. rho_solid/rho_fluid
    rho_is_set = false;
    ly = 0;  // nb: in reference length units
    lx = 0;
    lz = 0;
    is_initialized = false;
}

template <typename Real>
NonDimensional<Real>::NonDimensional() : l_ref(1), rho_f(1)
{
    re = 0;
    ma = -1;
    ma_is_set = false;
    rho_solid = -1;  // nb: in rho_fluid units i.e. rho_solid/rho_fluid
    rho_is_set = false;
    ly = 0;  // nb: in reference length units
    lx = 0;
    lz = 0;
    is_initialized = false;
}

/**
 * Initializes the dimensionless parameters from Re number and dimensionless
 * domain dimension.
 * @param re_ Reynolds' number
 * @param lx_
 * @param ly_
 * @param lz_
 */
template <typename Real>
NonDimensional<Real> NonDimensional<Real>::initReLxLyLz(Real re_, Real lx_, Real ly_, Real lz_)
{
    re = re_;
    lx = lx_;
    ly = ly_;
    lz = lz_;
    is_initialized = true;
    return *this;
}

template <typename Real>
NonDimensional<Real> NonDimensional<Real>::initReMaLxLyLz(
    Real re_, Real ma_, Real lx_, Real ly_, Real lz_)
{
    re = re_;
    lx = lx_;
    ly = ly_;
    lz = lz_;
    ma = ma_;
    ma_is_set = true;
    is_initialized = true;
    return *this;
}

template <typename Real>
NonDimensional<Real> NonDimensional<Real>::initReLxLyLzRhos(
    Real re_, Real lx_, Real ly_, Real lz_, Real rho_solid_)
{
    re = re_;
    lx = lx_;
    ly = ly_;
    lz = lz_;
    rho_solid = rho_solid_;
    rho_is_set = true;
    is_initialized = true;
    return *this;
}

template <typename Real>
NonDimensional<Real> NonDimensional<Real>::printParameters() const
{
    plb::pcout << "- DIMENSIONLESS PARAMETERS ------" << std::endl;
    plb::pcout << std::setw(12) << "re:" << std::setw(12) << re << std::endl;
    plb::pcout << std::setw(12) << "ma:" << std::setw(12) << ma << std::endl;
    plb::pcout << std::setw(12) << "ma_is_set:" << std::setw(12) << ma_is_set << std::endl;
    plb::pcout << std::setw(12) << "rho_solid:" << std::setw(12) << rho_solid << std::endl;
    plb::pcout << std::setw(12) << "rho_is_set:" << std::setw(12) << rho_is_set << std::endl;
    plb::pcout << std::setw(12) << "lx:" << std::setw(12) << lx << std::endl;
    plb::pcout << std::setw(12) << "ly:" << std::setw(12) << ly << std::endl;
    plb::pcout << std::setw(12) << "lz:" << std::setw(12) << lz << std::endl;
    plb::pcout << std::setw(12) << "l_ref:" << std::setw(12) << l_ref << std::endl;
    plb::pcout << std::setw(12) << "rho_f:" << std::setw(12) << rho_f << std::endl;
    plb::pcout << std::endl;
    return *this;
}

template <typename Real, typename Int>
Numerics<Real, Int>::Numerics()
{
    tau = -1;
    u_lb = -1;
    l_ref_lu = -1;
    lx_domain = -1;
    ly_domain = -1;
    lz_domain = -1;
    max_iter = -1;
    cs2 = 1. / 3.;
    is_initialized = false;
    rhof0 = (Real)1.0;
    dimless = new NonDimensional<Real>();
}

template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::initLU(
    Int l_ref_lu_, Real tau_lu_, Real u_lb_, Int lx_lu, Int ly_lu, Int lz_lu)
{
    l_ref_lu = l_ref_lu_;
    tau = tau_lu_;
    u_lb = u_lb_;
    lx_domain = lx_lu;
    ly_domain = ly_lu;
    lz_domain = lz_lu;
    is_initialized = true;
    dimless->initReLxLyLz(
        l_ref_lu_ * u_lb_ / getNulu(tau_lu_), lx_lu / l_ref_lu_, ly_lu / l_ref_lu_,
        lz_lu / l_ref_lu_);
    return *this;
}

template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::initReTauUlb(
    Real re, Real tau_lu_, Real u_lb_, Int lx_lu, Int ly_lu, Int lz_lu)
{
    l_ref_lu = re * getNulu(tau_lu_) / u_lb_;
    setTau(tau_lu_);
    setUlb(u_lb_);
    lx_domain = lx_lu;
    ly_domain = ly_lu;
    lz_domain = lz_lu;
    is_initialized = true;
    dimless->initReLxLyLz(re, lx_lu / l_ref_lu, ly_lu / l_ref_lu, lz_lu / l_ref_lu);
    return *this;
}

template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::initLrefReUlb(
    Int l_ref_lu_, Real re, Real u_lb_, Int lx_lu, Int ly_lu, Int lz_lu)
{
    l_ref_lu = l_ref_lu_;
    setTau(3. * (u_lb_ * l_ref_lu_ / re) + 0.5);
    setUlb(u_lb_);
    lx_domain = lx_lu;
    ly_domain = ly_lu;
    lz_domain = lz_lu;
    is_initialized = true;
    dimless->initReLxLyLz(
        l_ref_lu_ * u_lb_ / getNulu(tau), lx_lu / l_ref_lu_, ly_lu / l_ref_lu_, lz_lu / l_ref_lu_);
    return *this;
}

template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::initLrefReTau(
    Int l_ref_lu_, Real re, Real tau_, Int lx_lu, Int ly_lu, Int lz_lu)
{
    l_ref_lu = l_ref_lu_;
    setTau(tau_);
    setUlb(re * getNulu(tau_) / l_ref_lu_);
    lx_domain = lx_lu;
    ly_domain = ly_lu;
    lz_domain = lz_lu;
    is_initialized = true;
    dimless->initReLxLyLz(
        l_ref_lu_ * u_lb / getNulu(tau), lx_lu / l_ref_lu_, ly_lu / l_ref_lu_, lz_lu / l_ref_lu_);
    return *this;
}

/// Generates the parameters in lattice units from an object with the non
/// dimensional parameters and the lattice resolution, i.e. l_ref_lu_. l_ref_lu_
/// is Real, but the generated box dimension are rounded to int! \param
/// l_ref_lu_ resolution (reference length in lattice units) \param dimless_
/// NonDimensional<Real> object that has already been initialized
/// @param u_lb_ if given > 0, it overrides the Ma number in the dimless_
/// @param tau_ if given > 0, it overrides the Ma number in the dimless_
/// @return *this
template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::initLrefluNodim(
    Int l_ref_lu_, NonDimensional<Real> *dimless_, Real u_lb_, Real tau_)
{
    assert(dimless_->initialized());
    l_ref_lu = l_ref_lu_;
    if (u_lb_ > getEpsilon(1)) {
        setUlb(u_lb_);
        setTau(3. * (u_lb_ * l_ref_lu_ / dimless_->getRe()) + 0.5);
        if (tau_ > 0.499)
            pcout << "initLrefluNodim() Waring: you cannot set both u_lb_ "
                     "and tau_. Giving priority to u_lb_...."
                  << std::endl;
    } else if (tau_ > 0.499) {
        setUlb(dimless_->getRe() * getNulu(tau_) / l_ref_lu_);
        setTau(tau_);
    } else {
        pcout << "Detected invalid tau value, trying to recompute it from Ma number..."
              << std::endl;
        setUlb(dimless_->getMa() * getCs());
        setTau(getCs2Inv() * (u_lb * l_ref_lu_ / dimless_->getRe()) + 0.5);
    }
    lx_domain = util::roundToInt(dimless_->getLx() * l_ref_lu_);
    ly_domain = util::roundToInt(dimless_->getLy() * l_ref_lu_);
    lz_domain = util::roundToInt(dimless_->getLz() * l_ref_lu_);
    is_initialized = true;
    delete dimless;
    dimless = dimless_;
    return *this;
}

/**
 * Writes numerical parameters in lattice units to screen
 * @tparam Real
 * @tparam Int
 * @return *this
 */
template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::printParameters()
{
    plb::pcout << "- NUMERICAL PARAMETERS ------" << std::endl;
    plb::pcout << std::setw(12) << "tau:" << std::setw(12) << tau << std::endl;
    plb::pcout << std::setw(12) << "u_lb:" << std::setw(12) << u_lb << std::endl;
    plb::pcout << std::setw(12) << "l_ref_lu:" << std::setw(12) << l_ref_lu << std::endl;
    plb::pcout << std::setw(12) << "lx_domain:" << std::setw(12) << lx_domain << std::endl;
    plb::pcout << std::setw(12) << "ly_domain:" << std::setw(12) << ly_domain << std::endl;
    plb::pcout << std::setw(12) << "lz_domain:" << std::setw(12) << lz_domain << std::endl;
    plb::pcout << std::setw(12) << "lz_domain:" << std::setw(12) << lz_domain << std::endl;
    plb::pcout << std::setw(12) << "max_iter:" << std::setw(12) << max_iter << std::endl;
    plb::pcout << std::setw(12) << "cs2:" << std::setw(12) << cs2 << std::endl;
    plb::pcout << std::endl;
    return *this;
}

/**
 * Writes numerical parameters to a log file in the output folder
 * @tparam Real
 * @tparam Int
 * @return *this
 */
template <typename Real, typename Int>
Numerics<Real, Int> &Numerics<Real, Int>::writeParametersLog() const
{
    std::string fullName = global::directories().getLogOutDir() + "lattice_units_Log.dat";
    plb_ofstream ofile(fullName.c_str());
    ofile << "- NUMERICAL PARAMETERS ------" << std::endl;
    ofile << std::setw(12) << "tau:" << std::setw(12) << tau << std::endl;
    ofile << std::setw(12) << "u_lb:" << std::setw(12) << u_lb << std::endl;
    ofile << std::setw(12) << "l_ref_lu:" << std::setw(12) << l_ref_lu << std::endl;
    ofile << std::setw(12) << "lx_domain:" << std::setw(12) << lx_domain << std::endl;
    ofile << std::setw(12) << "ly_domain:" << std::setw(12) << ly_domain << std::endl;
    ofile << std::setw(12) << "lz_domain:" << std::setw(12) << lz_domain << std::endl;
    ofile << std::setw(12) << "lz_domain:" << std::setw(12) << lz_domain << std::endl;
    ofile << std::setw(12) << "max_iter:" << std::setw(12) << max_iter << std::endl;
    ofile << std::setw(12) << "cs2:" << std::setw(12) << cs2 << std::endl;
    ofile << std::endl;
    return *this;
}

/**
 * Returns Palabos IncomprFlowParam<Real> object. Important Note:
 * IncomprFlowParam<Real> only handles integer values for the "resolution" that
 * it is actually the reference length in lattice units. Thus, the conversion is
 * not perfect because the getLref of Numerics will be rounded to int before
 * being passed to IncomprFlowParam.
 * @tparam Real
 * @tparam Int
 * @return *this
 */
template <typename Real, typename Int>
IncomprFlowParam<Real> Numerics<Real, Int>::getIncomprFlowParam() const
{
    assert(dimless->initialized() && initialized());
    return IncomprFlowParam<Real>(
        getUlb(),                     // uMax
        dimless->getRe(),             // Re
        util::roundToInt(getLref()),  // N
        dimless->getLx(),             // lx
        dimless->getLy(),             // ly
        dimless->getLz());
}
}  // namespace incompressible_simulation_parameters