#ifndef ADVECTION_DIFFUSION_FD_H
#define ADVECTION_DIFFUSION_FD_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

#include <memory>

namespace plb {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3rd order Weighted Essentially Non-Oscillatory
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*********************WENO procedure for the convective term (includes diffusion term and settlinf velocity)******************************************/

    template<typename T>
    class WENO3 : public BoxProcessingFunctional3D
    {
    public:
        WENO3(T d_, T eps_, bool neumann_, plint nx_, plint ny_, plint nz_);
        virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields );
        virtual WENO3<T>* clone() const;
        virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    private:
        T d, eps;
        bool neumann;
        plint nx, ny, nz;
    };

/******************************3rd order Runge-Kutta for the time discretization***************************************************/
/* ******** RK3_Step1_functional3D ****************************************** */

template<typename T>
class RK3_Step1_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual RK3_Step1_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

  /* ******** RK3_Step2_functional3D ****************************************** */

template<typename T>
class RK3_Step2_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual RK3_Step2_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* ******** RK3_Step3_functional3D ****************************************** */

template<typename T>
class RK3_Step3_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
  virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
  virtual RK3_Step3_functional3D<T>* clone() const;
  virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  virtual BlockDomain::DomainT appliesTo() const;
};

}

#endif  // ADVECTION_DIFFUSION_FD_H
