///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifndef HYPERELASTICITY_H
#define HYPERELASTICITY_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Valanis - Landel Hyperelasticity
// Attention: Hill's stability criterion (aka Drucker's condition)
// See Xu 2015 for more details
///////////////////////////////////////////////////////////////////////////////
// Trianglular Elements
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar f_tr(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa){
    // Custom Energy
    //return miu * (1. / 4.) * std::pow(x - 1., 4.);

    // Skalak 1973
    //return (miu / 8.) * std::pow(x, 4.) - (miu / 4.) * std::pow(x, 2.);

    // Custom Energy Combo
    return miu*(1./4.)*std::pow(x - 1., 4.) + (lambda/8.)*std::pow(x, 4.) - (lambda/4.)*std::pow(x, 2.);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar f_prime_tr(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa){
    // Custom Energy
    //return miu * std::pow(x - 1., 3.);

    // Skalak 1973
    //return (miu / 2.) * std::pow(x, 3.) - (miu / 2.) * x;

    // Custom Energy Combo
    return miu*std::pow(x - 1., 3.) + (lambda/2.)*std::pow(x, 3.) - (lambda/2.)*x;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar g_tr(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa)
{
    // Custom Energy
    //return lambda * (1. / 4.) * std::pow(x - 1., 4.);

    // Skalak 1973
    //return (lambda / 8.) * std::pow(x, 4.) - (lambda / 4.) * std::pow(x, 2.);

    // Other Energy (PD) that conserves Area
    return 0.;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar g_prime_tr(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa)
{
    // Custom Energy
    //return lambda * std::pow(x - 1., 3.);

    // Skalak 1973
    //return (lambda / 2.) * std::pow(x, 3.) - (lambda / 2.) * x;

    // Other Energy (PD) that conserves Area
    return 0.;
}
///////////////////////////////////////////////////////////////////////////////
// Tetrahedra
///////////////////////////////////////////////////////////////////////////////
inline Scalar f_tet(const Scalar& x, const Scalar& miu, const Scalar& lambda,
    const Scalar& kappa)
{
    return miu * (1. / 4.) * std::pow(x - 1., 4.);
}
///////////////////////////////////////////////////////////////////////////////
inline Scalar f_prime_tet(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa)
{
    return miu * std::pow(x - 1., 3.);
}
///////////////////////////////////////////////////////////////////////////////
inline Scalar g_tet(const Scalar& x, const Scalar& miu, const Scalar& lambda,
    const Scalar& kappa)
{
    return 0.;
}
///////////////////////////////////////////////////////////////////////////////
inline Scalar g_prime_tet(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa)
{
    return 0.;
}
///////////////////////////////////////////////////////////////////////////////
inline Scalar h_tet(const Scalar& x, const Scalar& miu, const Scalar& lambda,
    const Scalar& kappa)
{
    return lambda * (1. / 4.) * std::pow(x - 1., 4.);
}
///////////////////////////////////////////////////////////////////////////////
inline Scalar h_prime_tet(const Scalar& x, const Scalar& miu,
    const Scalar& lambda, const Scalar& kappa)
{
    return lambda * std::pow(x - 1., 3.);
}
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // HYPERELASTICITY_H
///////////////////////////////////////////////////////////////////////////////
