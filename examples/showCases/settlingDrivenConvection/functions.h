#include <cmath>
#include <vector>

#ifndef FUNCTIONS_SDC_H
#define FUNCTIONS_SDC_H
#define M_PI        3.14159265358979323846

namespace plb {


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Buoyant force term
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename T, template<typename U> class FluidDescriptor>
class ScalarBuoyanTermProcessor3D :  public BoxProcessingFunctional3D
{
public:

    ScalarBuoyanTermProcessor3D(T gravity_, T rho0_, T rhoP_, T TotalVolFrac_, T dt_,
                                 Array<T,FluidDescriptor<T>::d> dir_);

    virtual void processGenericBlocks( Box3D domain, std::vector<AtomicBlock3D*> fields );
    virtual ScalarBuoyanTermProcessor3D<T,FluidDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:
      T gravity, rho0, rhoP, TotalVolFrac, dt;
      Array<T,FluidDescriptor<T>::d> dir;
};



template< typename T, template<typename U> class FluidDescriptor>
ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::
        ScalarBuoyanTermProcessor3D(T gravity_, T rho0_, T rhoP_, T TotalVolFrac_, T dt_, Array<T,FluidDescriptor<T>::d> dir_)
    :  gravity(gravity_), rho0(rho0_), rhoP(rhoP_), TotalVolFrac(TotalVolFrac_), dt(dt_),
       dir(dir_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T,FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}


template< typename T, template<typename U> class FluidDescriptor >
void ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    typedef FluidDescriptor<T> D;
    enum {
        forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt
    };

    PLB_PRECONDITION(fields.size()==3);
    BlockLattice3D<T, FluidDescriptor>* fluid = dynamic_cast<BlockLattice3D<T, FluidDescriptor>*>(fields[0]);
    ScalarField3D<T>* volfracfield = dynamic_cast<ScalarField3D<T>*>(fields[1]);
    ScalarField3D<T>* densityfield = dynamic_cast<ScalarField3D<T>*>(fields[2]);


    Dot3D offset1 = computeRelativeDisplacement(*fluid, *volfracfield);
    Dot3D offset2 = computeRelativeDisplacement(*fluid, *densityfield);


    Array<T,D::d> gravOverrho0 (
            gravity*dir[0]/rho0,
            gravity*dir[1]/rho0,
            gravity*dir[2]/rho0 );


    T maxiT = 1.0/dt;
    T iT = fluid->getTimeCounter().getTime();
    T gain = util::sinIncreasingFunction(iT, maxiT);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ)
            {

                T localVolfrac = volfracfield->get(iX+offset1.x,iY+offset1.y,iZ+offset1.z);
                // Computation of the Boussinesq force
                T *force = fluid->get(iX,iY,iZ).getExternal(forceOffset);
                T dens = densityfield->get(iX+offset2.x,iY+offset2.y,iZ+offset2.z);
                // volfracfield is the order-0 moment of the advection-diffusion lattice.
                const T diffT = rhoP-rho0;
                for (pluint iD = 0; iD < D::d; ++iD)
                {
			                 force[iD] = - gain * gravOverrho0[iD] * (diffT * localVolfrac * TotalVolFrac + (dens - rho0)*(1 - localVolfrac*TotalVolFrac));
                }
            }
        }
    }
}

template< typename T, template<typename U> class FluidDescriptor>
ScalarBuoyanTermProcessor3D<T,FluidDescriptor>*
    ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::clone() const
{
    return new ScalarBuoyanTermProcessor3D<T,FluidDescriptor>(*this);
}

template< typename T, template<typename U> class FluidDescriptor >
void ScalarBuoyanTermProcessor3D<T,FluidDescriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
        modified[2] = modif::nothing;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Settling velocity field
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename T>
class Get_v_sedimentation :
    public BoxProcessingFunctional3D
{
public:

    Get_v_sedimentation(T rhoP_, T Dp_, T convers_,  T mu_, T g_);

    virtual void processGenericBlocks( Box3D domain,
                          std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual Get_v_sedimentation<T>* clone() const;

  private:
    T rhoP; T Dp; T convers; T mu; T g;

};

template< typename T >
Get_v_sedimentation<T>::Get_v_sedimentation(T rhoP_, T Dp_, T convers_, T mu_, T g_)
            : rhoP(rhoP_), Dp(Dp_), convers(convers_), mu(mu_), g(g_)
{
}


template< typename T >
void Get_v_sedimentation<T>::processGenericBlocks (Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    //typedef DensityDescriptor<T> D;
    PLB_PRECONDITION(atomicBlocks.size()==3);
    ScalarField3D<T>* densityField = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[0]);
    ScalarField3D<T>* volfracField = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[1]);
    ScalarField3D<T>* v_sed = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[2]);


    Dot3D offset1 = computeRelativeDisplacement(*densityField, *volfracField);
    Dot3D offset2 = computeRelativeDisplacement(*densityField, *v_sed);


    for (plint iX=domain.x0; iX<=domain.x1; ++iX){
        for (plint iY=domain.y0; iY<=domain.y1; ++iY){
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ){


                T dens = densityField->get(iX, iY, iZ);
                T localvolfrac = volfracField->get(iX+offset1.x,iY+offset1.y,iZ+offset1.z);
                T vel_sed;

                if (localvolfrac>0) vel_sed=-convers*(0.5*Dp*Dp*g*(rhoP-dens))/(9*mu);
                else vel_sed = 0;

                v_sed->get(iX+offset2.x,iY+offset2.y,iZ+offset2.z)=vel_sed;



            }
        }
    }
}

template< typename T >
Get_v_sedimentation<T>*
    Get_v_sedimentation<T>::clone() const
{
    return new Get_v_sedimentation<T>(*this);
}

template< typename T >
void Get_v_sedimentation<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
        modified[2] = modif::staticVariables;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Units conversion class
///////////////////////////////////////////////////////////////////////////////////////////////////////////


    /// A useful class for the conversion between dimensionless and lattice units.
    template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
    class RayleighTaylorFlowParam {
    public:
        /// Constructor
        /** \param Re_  Reynolds number
         * \param Ra_  Raylegh number
         *  \param Pr_  Prandtl number
         *  \param coldTemperature_  minimum temperature
         *  \param hotTemperature_  maximum temperature
         *  \param deltaT_ time discretization number
         *  \param N_  resolution (a lattice of size 1 has N_+1 cells)
         *  \param lx_ x-length in dimensionless units (e.g. 1)
         *  \param ly_ y-length in dimensionless units (e.g. 1)
         *  \param lz_ z-length in dimensionless units (e.g. 1)
         */
        RayleighTaylorFlowParam(T Ri_, T Gr_, T uMax_, T uCar_,  T resolution_,
                                T lx_, T ly_, T lz_, T Di_=T() )
            : Ri(Ri_), Gr(Gr_), uMax(uMax_), uCar(uCar_),
              resolution(resolution_), lx(lx_), ly(ly_), lz(lz_), Di(Di_)
        { }
        /// Peclet number
        //T getPe() const      { return getResolution() * getLatticeU() / getLatticeKappa(); }
        /// Rayleigh number
        T getRi() const      { return Ri; }
        /// Grasshoff number
        T getGr() const      { return Gr; }

        T getUcar() const       { return uCar;}

        T getResolution() const { return resolution; }
        /// x-length in dimensionless units
        T getLx() const      { return lx; }
        /// y-length in dimensionless units
        T getLy() const      { return ly; }
        /// z-length in dimensionless units
        T getLz() const      { return lz; }
        /// lattice spacing in dimensionless units
        T getDeltaX() const  { return 1 /(T)resolution; }
        /// caracteristic velocity
        //T getCarvel() const  { return uCar; }
        /// time step in dimensionless units
        //T getDeltaT() const  { return getDeltaX()*getDeltaX(); }
        T getDeltaT() const  { return getLatticeU()*getDeltaX()/getUcar(); }
        /// conversion from dimensionless to lattice units for space coordinate
        plint nCell(T l) const { return (plint)(l/getDeltaX()+(T)0.5); }
        /// conversion from dimensionless to lattice units for time coordinuLbate
        plint nStep(T t) const { return (plint)(t/getDeltaT()+(T)0.5); }
        /// number of lattice cells in x-direction
        plint getNx() const    { return nCell(lx)+1; }
        /// number of lattice cells in y-direction
        plint getNy() const    { return nCell(ly)+1; }
        /// number of lattice cells in z-direction
        plint getNz() const    { return nCell(lz)+1; }
        /// velocity in lattice units (proportional to Mach number)
        //T getLatticeU() const       { return (getRe()*getLatticeNu())/(getNz()-1); }
        T getLatticeU() const       { return uMax;}
        /// Reynolds number
        T getRe() const      { return std::sqrt(getGr()/getRi()); }
        /// viscosity in lattice units
        T getLatticeNu() const      { return getLatticeU()*(getNz()-1)/getRe(); }

        /// thermal conductivity in lattice units
        T getLatticeKappa() const   { return (getDeltaT()/(getDeltaX()*getDeltaX())); }
        /// viscosity in lattice units
        T getLatticeGravity() const { return getDeltaT() * getDeltaT() / getDeltaX(); }
        /// relaxation time
        T getSolventTau() const   { return nsDescriptor<T>::invCs2*getLatticeNu()+(T)0.5; }
        /// relaxation frequency
        T getSolventOmega() const { return (T)1 / getSolventTau(); }
        /// relaxation time
        T getTemperatureTau() const    { return adDescriptor<T>::invCs2*getLatticeKappa()+(T)0.5; }
        /// relaxation frequency
        T getTemperatureOmega() const  { return (T)1 / getTemperatureTau(); }
    private:
        T Ri, Gr, uMax, uCar, resolution, lx, ly, lz, Di;
    };

    template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
    void writeLogFile(RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> const& parameters,
                      std::string const& title)
    {
        std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
        std::ofstream ofile(fullName.c_str());
        ofile << title << "\n\n";
        ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
        //ofile << "Peclet number:             Pe=" << parameters.getPe() << "\n";
        ofile << "Richardson number:          Ri=" << parameters.getRi() << "\n";
        ofile << "Grasshoff number:            Gr=" << parameters.getGr() << "\n";
        ofile << "Kinematic viscosity:       Nu=" << parameters.getLatticeNu() << "\n";
        ofile << "Thermal conductivity:   Kappa=" << parameters.getLatticeKappa() << "\n";
        ofile << "Lattice resolution:         N=" << parameters.getResolution() << "\n";
        ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
        ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
        ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
        ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
        ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
        ofile << "Solvent omega:        omega_S=" << parameters.getSolventOmega() << "\n";
        ofile << "Temperature omega:    omega_T=" << parameters.getTemperatureOmega() << "\n";
        ofile << "Caracteristic vel:       uLb=" << parameters.getLatticeU() << "\n";
        ofile << "Number of cells x:       Nx=" << parameters.getNx() << "\n";
        ofile << "Number of cells y:       Ny=" << parameters.getNy() << "\n";
        ofile << "Number of cells z:       Nz=" << parameters.getNz() << "\n";
    }


}

#endif
