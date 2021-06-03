/* This file is part of the Palabos library.
* Copyright (C) 2009 Jonas Latt
* E-mail contact: jonas@lbmethod.org
* The most recent release of Palabos can be downloaded at 
* <http://www.lbmethod.org/palabos/>
*
* The library Palabos is free software: you can redistribute it and/or
* modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
*
* The library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
*
**/

#include "palabos3D.h"
#include "palabos3D.hh"   // include full template code
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
// #define DESCRIPTOR D3Q27Descriptor
#define DESCRIPTOR D3Q19Descriptor

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
template<typename T>
class Param {
public:
    plint N;                    // Resolution
    plint nx, ny, nz;           // Domain size (LB units)
    T lx, ly, lz;               // Domain size (Physical units)
    T dx, dt;                   // Space and time steps
    T Re, Ma, tc;               // Dimensionless parameters
    T u0;                       // Velocity (LB units)
    T soundSpeed;               // Speed of sound (Physical units)
    T cSmago;                   // Smagorinsky constant (turbulent flows)
    T omega, tau, nu;           // Collision parameters (LB units)
    T omega3, omega4;           // Relaxation frequencies of 3rd and 4th order moments
    T omega5, omega6;           // Relaxation frequencies of 5th and 6th order moments
    bool iniPressure, iniFneq;  // Init parameters (init density through Poisson equation,
                                //                  include nonequilibrium part)
    plint intIniFneq, intIniPressure;

    T tAdim, vtsT;              // Simulation time and frequency of .vti outputs
    std::string hoOmega;

    std::string lbm;
    std::string dynName, outDirName, logfileName;
    std::string fnameBase;

    Param() { }
    Param(std::string xmlFname) {
        //                                             1          2  3 4 
        //pcout << "Example: " << global::argv(0) << " config.xml Re N Ma " << std::endl;


        //////// Parameters from the xml file
        XMLreader document(xmlFname);
        document["lattice"]["lbm"].read(lbm);
        document["lattice"]["dynName"].read(dynName);
        document["lattice"]["hoOmega"].read(hoOmega);
        document["simuParam"]["soundSpeed"].read(soundSpeed);
        document["simuParam"]["dx"].read(dx);
        document["simuParam"]["nz"].read(nz);
        document["simuParam"]["cSmago"].read(cSmago);
        document["initialization"]["intIniFneq"].read(intIniFneq);
        document["initialization"]["intIniPressure"].read(intIniPressure);
        document["io"]["output"].read(outDirName);
        document["io"]["tAdim"].read(tAdim);
        document["io"]["vtsT"].read(vtsT);

        //////// Numerical discretization 
        // dx = 0.01;                                // Space step in physical units [m] (now read from .xml file  to avoid har coded value)
        T cs = ::sqrt(DESCRIPTOR<T>::cs2);        // Sound speed in LB units
        // T gamma = 1.4;                            // Specific heat ratio of air (diatomic gas)
        // T rGas = 287;                             // Gas constant                       
        // T Tref = 273.15 + 20.;                    // Ambient temperature (20°C)
        // soundSpeed = ::sqrt(gamma*rGas*Tref;      // Sound speed in physical units (m/s)
        //            = 343.20208332701009;
        // soundSpeed = 340.;                        // Simplified value (now read from .xml file to avoid har coded value)       
        dt = (cs/soundSpeed)*dx;                  // Time step in physical units [s]
        
        //////// Simulation domain parameters
        global::argv(3).read(N);
        nx = N; ny = N; //nz = 3; nz is given in the xml. file to avoid any hard coded value 
        lx = nx * dx; ly = ny * dx; lz = nz * dx;

        //////// Dimensionless parameters
        global::argv(4).read(Ma); 
        u0 = (T)(Ma * cs);                             // Ma = u0/cs
        tc = (T)(N/u0);
        global::argv(2).read(Re);
        nu = (u0*N)/Re;                           // Re = (u0*N)/nu
        tau = nu/DESCRIPTOR<T>::cs2;
        omega = 1./(tau + 0.5);
        if (hoOmega == "SRT") { 
            omega3 = omega;
            omega4 = omega;
            omega5 = omega;
            omega6 = omega;
        } else if (hoOmega == "REG") { 
            omega3 = 1.;
            omega4 = 1.;
            omega5 = 1.;
            omega6 = 1.;
        } else {
            pcout << "Error: Relaxation of high-order moments not correct." << std::endl;
            exit(-1);
        }

        //////// Improved initialization step
        iniPressure = intIniPressure == 0 ? false : true;
        if (iniPressure){
            pcout << "Error: Poisson solver does not work for 3D models." << std::endl;
            exit(-1);
        }
        iniFneq = intIniFneq == 0 ? false : true;

        //////// File name for log and stats
        std::stringstream fnameBaseStr;
        //pcout << "Syntax: " << global::argv(0) << " config.xml Re N Ma iniPressure iniFneq" << std::endl;
        fnameBaseStr << std::setprecision(7) << Re;
        fnameBaseStr << "_";
        fnameBaseStr << N;
        fnameBaseStr << "_0_";
        fnameBaseStr << std::setprecision(7) << (int)(Ma*100.);
        fnameBaseStr << "_";
        fnameBaseStr << iniFneq ? 1 : 0;
        fnameBaseStr << "_";
        fnameBaseStr << iniPressure ? 1 : 0;
        fnameBaseStr >> fnameBase;
    }
    
    void writeLogFile() {
        plb_ofstream fout((outDirName+"/log_"+fnameBase+".dat").c_str());//========== Numerical Parameters ===========//

        fout << " //======== LBM Parameters ===============// " << std::endl;
        fout << "Lattice  -->           "<< lbm << std::endl;
        fout << "Dynamics -->           "<< dynName << std::endl;
        fout << "HO Relaxation Rype --> "<< hoOmega << std::endl;
        fout << std::endl;

        fout << " //======== Physical Parameters ==========// " << std::endl;
        fout << "Flow properties (dimensionless):    " << std::endl;
        fout << "Re = " << Re << std::endl;
        fout << "Ma = " << Ma << std::endl;
        fout << "Flow properties (physical units):    " << std::endl;
        fout << "nu = " << nu*dx*dx/dt << " [m2/s]" << std::endl;
        fout << "c  = " << soundSpeed << " [m/s]" << std::endl;
        fout << "u0 = " << Ma * ::sqrt(DESCRIPTOR<T>::cs2) * dx/dt << " [m/s]" << std::endl;
        fout << "tc = " << tc * dt << " [s]" << std::endl;
        fout << "Geometry (physical units):    " << std::endl;
        fout << "lx = " << lx << " [m]" << std::endl;
        fout << "ly = " << ly << " [m]" << std::endl;
        fout << "lz = " << lz << " [m]" << std::endl;
        fout << std::endl;

        fout << " //======== Numerical Parameters =========// " << std::endl;
        fout << "Numerical discretization (physical units):    " << std::endl;        
        fout << "dx = " << dx << " [m]" << std::endl;
        fout << "dt = " << dt << " [s]" << std::endl;
        fout << "Geometry (LB units):    " << std::endl;
        fout << "N  = " << N << " (resolution)" << std::endl;
        fout << "nx = " << nx << std::endl;
        fout << "ny = " << ny << std::endl;
        fout << "nz = " << nz << std::endl;
        fout << "Flow properties (LB units):    " << std::endl;
        fout << "nuLB = " << nu << std::endl;
        fout << "u0LB = " << Ma * ::sqrt(DESCRIPTOR<T>::cs2) << std::endl;
        fout << "tcLB = " << round(tc) << " (" << tc << ")" << std::endl;
        fout << "Collision parameters (LB units):    " << std::endl;
        fout << "tau = " << tau << std::endl;
        fout << "omega = " << omega << std::endl;
        if (lbm == "D3Q19"){
            fout << "omega3 = " << omega3 << std::endl;
            fout << "omega4 = " << omega4 << std::endl;
        } else if (lbm == "D3Q27"){
            fout << "omega3 = " << omega3 << std::endl;
            fout << "omega4 = " << omega4 << std::endl;
            fout << "omega5 = " << omega5 << std::endl;
            fout << "omega6 = " << omega6 << std::endl;
        }
        fout << "Large Eddy Simulation parameters:    " << std::endl;
        fout << "cSmago = " << cSmago << std::endl;
        fout << std::endl;

        fout << " //======== Improved Initialization ======// " << std::endl;
        fout << "Initializes pressure (Poisson solver): ";
        if (iniPressure) fout << "True" << std::endl;
        else fout << "False" << std::endl;
        fout << "Initializes fNeq (using FD gradients): ";
        if (iniFneq) fout << "True" << std::endl;
        else fout << "False" << std::endl;
        fout << std::endl;

        fout << " //======== Simulation parameters ======// " << std::endl;
        fout << "output= " << outDirName << std::endl;
        fout << "tAdim = " << tAdim << " * tc" << std::endl;
        fout << "      = " << (int)(tAdim * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(tAdim * tc) << " [iterations]" << std::endl;
        fout << "vtsT  = " << vtsT << " * tc" << std::endl;
        fout << "      = " << (int)(vtsT * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(vtsT * tc) << " [iterations]" << std::endl;
    }
};

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> *getDynamics(Param<T> &param)
{
    Dynamics<T,Descriptor> *dyn;
    if (param.dynName == "BGK_Ma2") { 
        dyn = new BGKdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "RM") { 
        dyn = new RMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "HM") { 
        dyn = new HMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CM") { 
        dyn = new CMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CHM") { 
        dyn = new CHMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "K") { 
        dyn = new Kdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "GH") { 
        dyn = new GHdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "RR") { 
        dyn = new RRdynamics<T,Descriptor>(param.omega);
    } else {
        pcout << "Error: dynamics name does not exist." << std::endl;
        exit(-1);
    }

    return dyn;
}


//////// A functional, used to initialize the simulation
template<typename T>
class DoubleShearLayerInitialVelocityField {
public:
    DoubleShearLayerInitialVelocityField(Param<T> const& param_)
        : param(param_)
    { }
    
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;

        const plint nx = param.nx;
        const plint ny = param.ny;
        
        T ux = 0.;
        T uy = 0.;
        T uz = 0.;

        T kappa =80.;
        T delta = 0.05;
        T u0 = param.u0;

        T x = (T)iX/T(nx);
        T y = (T)iY/T(ny);
        if (y <= 0.5){
            // pcout << "top (iX,iY,iZ) = (" << iX << "," << iY << "," << iZ << ")" << std::endl;
            // pcout << "top (x,y) = (" << x << "," << y << "," << z << ")" << std::endl;
            ux   = u0*tanh(kappa*(y-0.25));
            uy   = u0*delta*sin(2.*M_PI*(x+0.25));
        }
        else{
            // pcout << "bottom (iX,iY,iZ) = (" << iX << "," << iY << "," << iZ << ")" << std::endl;
            // pcout << "bottom (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
            ux   = u0*tanh(kappa*(0.75-y));
            uy   = u0*delta*sin(2.*M_PI*(x+0.25));
        }

        u[0] = ux;
        u[1] = uy;
        u[2] = uz;
    }
    
private:
    Param<T> param;
};


//////// Initialize the simulation
void simulationSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                     const Param<T> &param)
{
    // Set periodic boundaries.
    lattice.periodicity().toggleAll(true); 

    // Initialize the simulation domain.
    initializeAtEquilibrium (
        lattice, lattice.getBoundingBox(),
        DoubleShearLayerInitialVelocityField<T>(param) );

    // Call initialize to get the lattice ready for the simulation.
    lattice.initialize();
}


//////// Post processing
/// Produce a GIF snapshot of the velocity-norm.
void writeGif(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                            *computeVelocityNorm(lattice), imSize, imSize);
}
/// Write the full velocity and the velocity-norm into a VTK file.
void writeVTS(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              Param<T> &param, plint iter)
{
    T dx = param.dx;
    T dt = param.dt;
    VtkStructuredImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "machNorm", 1.0/std::sqrt(DESCRIPTOR<T>::cs2));
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.0);
    vtkOut.writeData<3,float>(*computeGradient(*computeDensity(lattice)), "gradRho", 1.0/dx);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", (T)1/dt);
    vtkOut.writeData<6,float>(*computeStrainRate(*computeVelocity(lattice)), "strain", (T)1/dt);
}

T computeAveragedEnstrophy(MultiBlockLattice3D<T,DESCRIPTOR>& lattice)
{
    unique_ptr<MultiTensorField3D<T,3> > velocity = computeVelocity(lattice);
    unique_ptr<MultiTensorField3D<T,3> > vorticity = computeVorticity(*velocity);
    unique_ptr<MultiScalarField3D<T> > vorticityNorm = computeNorm(*vorticity);
    unique_ptr<MultiScalarField3D<T> > enstrophy = 
        multiply(0.5, *multiply(*vorticityNorm,*vorticityNorm) );
    return computeAverage(*enstrophy);
}

T computeRMSkinEnergy(MultiBlockLattice3D<T,DESCRIPTOR>& lattice)
{
    unique_ptr<MultiScalarField3D<T> > velocityNorm = computeVelocityNorm(lattice);
    unique_ptr<MultiScalarField3D<T> > kinEnergy = 
        multiply(0.5, *multiply(*velocityNorm,*velocityNorm) );
    T averagedKinEnergy = computeAverage(*kinEnergy);
    subtractInPlace(*kinEnergy, averagedKinEnergy);   
    multiplyInPlace(*kinEnergy, *kinEnergy);   
    return std::sqrt( computeAverage(*kinEnergy));
}

T computeRMSenstrophy(MultiBlockLattice3D<T,DESCRIPTOR>& lattice)
{
    unique_ptr<MultiTensorField3D<T,3> > velocity = computeVelocity(lattice);
    unique_ptr<MultiTensorField3D<T,3> > vorticity = computeVorticity(*velocity);
    unique_ptr<MultiScalarField3D<T> > vorticityNorm = computeNorm(*vorticity);
    unique_ptr<MultiScalarField3D<T> > enstrophy = 
        multiply(0.5, *multiply(*vorticityNorm,*vorticityNorm) );
    T averagedEns = computeAverage(*enstrophy);
    subtractInPlace(*enstrophy, averagedEns);
    multiplyInPlace(*enstrophy, *enstrophy);      
    return std::sqrt( computeAverage(*enstrophy));
}


std::unique_ptr<MultiScalarField3D<T> > solvePoissonEquation(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T epsilon)
{
    T beta = (T)1.0;
    Box3D domain = lattice.getBoundingBox();
    Box3D reducedDomain = domain.enlarge(-1);

    std::unique_ptr<MultiScalarField3D<T> > density = computeDensity(lattice);
    std::unique_ptr<MultiTensorField3D<T,3> > velocity = computeVelocity(lattice);
    std::unique_ptr<MultiScalarField3D<T> > pressure1 = multiply( DESCRIPTOR<T>::cs2, *add(-1., *density) );
    std::unique_ptr<MultiScalarField3D<T> > pressure2 = multiply( DESCRIPTOR<T>::cs2, *add(-1., *density) );
    std::unique_ptr<MultiScalarField3D<T> > rhs = computePoissonRHS(*velocity);  ;
    rhs = multiply( -(T)1, *rhs );
    T average = sqrt(computeAverage( *multiply(*rhs,*rhs),domain));
    pcout << "Poisson average = " << average << endl;
    
//     setToConstant(*rhs,domain,T());
    bool poissonConverged = false;
    util::ValueTracer<T> converge(1.0,(T)1000,1.0e-5);
    for (plint iPop = 0; iPop < 10000000; ++iPop)
    {
        
        poissonIterate(*pressure1,  *pressure2, *rhs, beta, domain);
        std::swap(pressure1,pressure2);
        
        T residual = computePoissonResidue(*pressure1, *rhs, reducedDomain);
//         T average  = computeAverage(*pressure1, reducedDomain);
        if (iPop % 5000 == 0) {
            pcout << "Iteration number = " << iPop << ", residual = " << residual/average << ", " << residual << endl;
        }
        
        converge.takeValue(residual/average);
        if (converge.hasConverged() || residual/average < epsilon) {
            pcout << "Poisson solver converged." << endl;
            poissonConverged = true;
            break;
        }
    }
    if (!poissonConverged) {
        pcout << "Poisson solver did not converge, continuing anyway." << endl;
    }

    return add(-computeAverage(*pressure1), *pressure1);
}

void accurateInitialCondition(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, const Param<T> &param, T epsilon)
{
    std::unique_ptr<MultiScalarField3D<T> > density;
    if (param.iniPressure) {
        density = solvePoissonEquation(lattice, epsilon);
        // Renormalize density to 1.
        *density = *multiply(*density, DESCRIPTOR<T>::invCs2);
    }
    density = computeDensity(lattice);

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("press", 0, 6), *density);

    if (param.iniFneq) {
        std::unique_ptr<MultiTensorField3D<T,3> > velocity = computeVelocity(lattice);
        std::unique_ptr<MultiTensorField3D<T,6> > S = computeStrainRate(*velocity);

        recomposeFromFlowVariables(lattice, *density, *velocity, *S);
    }
}

/// Solve the pressure from the velocity derivatives with Gauss-Seidel
template<typename T, template<typename U> class Descriptor>
void solvePressure(MultiBlockLattice3D<T,Descriptor>& lattice, T epsilon){
    std::unique_ptr<MultiScalarField3D<T> > density = computeDensity(lattice);
    std::unique_ptr<MultiTensorField3D<T,3> > velocity = computeVelocity(lattice);
    // initial pressure
    std::unique_ptr<MultiScalarField3D<T> > initialValue = multiply( DESCRIPTOR<T>::cs2, *add(-1., *density) );
    // the result holder
    std::unique_ptr<MultiScalarField3D<T> > result = multiply( DESCRIPTOR<T>::cs2, *add(-1., *density) );
    // right hand side of the equation
    std::unique_ptr<MultiScalarField3D<T> > rhs = computePoissonRHS(*velocity);  ;
    rhs = multiply( -(T)1, *rhs );
    
    GaussSeidelSolver( *initialValue, *result, *rhs, lattice.getBoundingBox(),epsilon, 10000000 );
    //     multiGridVCycle<T>(*initialValue, *result, *rhs, lattice.getBoundingBox());
    
    *result = *multiply(*result, DESCRIPTOR<T>::invCs2);
    *result = *add((T)1, *result);
    plb_ofstream poisson(createFileName("poisson",lattice.getNx(),6).c_str());
    poisson << *result;
    
    // writing the image
    ImageWriter<T> imageWriter("leeloo");
    plint imSizeX = 300;
    plint imSizeY = 300;
    imageWriter.writeScaledGif(createFileName("pressure",lattice.getNx(),6),*density,imSizeX,imSizeY);
    
    std::unique_ptr<MultiTensorField3D<T,6> > S = computeStrainRate(*velocity);
    recomposeFromFlowVariables(lattice, *density, *velocity, *S);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    if (argc != 5)
    {
        pcout << argc << std::endl;
        pcout << "Error! Wrong number of parameters." << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0) << " config.xml Re N Ma" << std::endl;
        pcout << "Example: " << (std::string)global::argv(0) << " config.xml 10000 256 0.1" << std::endl;
        exit(1);
    }

    ///// Read and print the parameters of the simulation
    Param<T> param(argv[1]);
    param.writeLogFile();
    ///// Output directory
    global::directories().setOutputDir(param.outDirName);

    ///// Initialize the diagonal relaxation matrix from the .xml file.
    // Q19 and Q27 do not have the same number of relaxation parameters!
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    if (param.lbm == "D3Q19"){
        allOmega[0] = param.omega;  // relaxation of M200 and cyclic permutations
        allOmega[1] = param.omega;  // relaxation of M110 and cyclic permutations
        allOmega[2] = param.omega3; // relaxation of M210 and cyclic permutations
        allOmega[3] = param.omega4; // relaxation of M220 and cyclic permutations
    } else if (param.lbm == "D3Q27"){
        allOmega[0] = param.omega;  // relaxation of M200 and cyclic permutations
        allOmega[1] = param.omega;  // relaxation of M110 and cyclic permutations
        allOmega[2] = param.omega3; // relaxation of M210 and cyclic permutations
        allOmega[3] = param.omega3; // relaxation of M111 and cyclic permutations
        allOmega[4] = param.omega4; // relaxation of M220 and cyclic permutations
        allOmega[5] = param.omega4; // relaxation of M211 and cyclic permutations
        allOmega[6] = param.omega5; // relaxation of M221 and cyclic permutations
        allOmega[7] = param.omega6; // relaxation of M222 and cyclic permutations
    } else {
        pcout << "Error: lbm name does not exist." << std::endl;
        exit(-1);
    }
    
    ///// Generate the dynamics and the corresponding lattice from the .xml file.
    Dynamics<T,DESCRIPTOR> *dyn = getDynamics<T,DESCRIPTOR>(param);
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
              param.nx, param.ny, param.nz, dyn);

    if (param.dynName == "RM") { 
        RMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "HM") { 
        HMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CM") { 
        CMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CHM") { 
        CHMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "K") { 
        Kdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "GH") { 
        GHdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "RR") { 
        RRdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    }

    lattice.toggleInternalStatistics(false);


    ///// Initialization from analytical profiles
    simulationSetup(lattice, param);
    ////////////////////////////////
    // Poisson not yet working... //
    ////////////////////////////////
    //// Either add fNeq, or properly initialize density (incompressible flows) or do nothing
    accurateInitialCondition(lattice, param,1.e-4);

    ///// Initial state is output
    writeVTS(lattice, param, 0);    

    ///// Simulation maximal time, and output frequency (in terms of iterations).
    plint vtsTout = param.vtsT * param.tc;
    plint tmax = param.tAdim * param.tc;
          
    plb_ofstream statsOut((param.outDirName+"/stats_" + param.fnameBase + ".dat").c_str());
    T previous_energy = std::numeric_limits<T>::max();

    ///// Main loop over time iterations.
    for (plint iT=0; iT<=tmax; ++iT) {
        // pcout << "iteration " << iT << std::endl;
        // T energy = 4*computeAverageEnergy(lattice);
        T kinEnergy = computeAverageEnergy(lattice);
        T kinEnergyAdim = kinEnergy/(0.5*param.u0*param.u0);
        T RMSkinEnergy = computeRMSkinEnergy(lattice);
        T RMSkinEnergyAdim = RMSkinEnergy/kinEnergy;

        T enstrophy = computeAveragedEnstrophy(lattice);
        T enstrophyAdim = enstrophy*(param.nx*param.nx)/(0.5*param.u0*param.u0);
        T RMSenstrophy = computeRMSenstrophy(lattice);
        T RMSenstrophyAdim = RMSenstrophy/enstrophy;

        // statsOut << iT/param.tc << " " 
        //          << kinEnergy << " " << kinEnergyAdim << " " << RMSkinEnergy << " " << RMSkinEnergyAdim << " "
        //          << enstrophy << " " << enstrophyAdim << " " << RMSenstrophy << " " << RMSenstrophyAdim << std::endl;    
        statsOut << iT/param.tc << " " 
                 << kinEnergyAdim << " " << " " << RMSkinEnergyAdim << " "
                 << enstrophyAdim << " " << " " << RMSenstrophyAdim << std::endl;    


        ///// Stability test based on the decrease of the kinetic energy.
        // The kinetic energy should be smaller than its initial value (weak condition)
        if (iT == 0) previous_energy = kinEnergy;
        if ((iT > (int)(0.01*param.tc)) and !(kinEnergy < previous_energy)) {
            pcout << "Catastrophic error: energy has increased or is NaN!" << std::endl;
            return 1;
        }
        // // The kinetic energy should decrease over time (strong condition)
        // if ((iT > (int)(0.1*param.tc)) and !(energy < previous_energy)) {
        //     pcout << "Catastrophic error: energy has increased or is NaN!" << std::endl;
        //     return 1;
        // }
        // if (iT % (int)(0.01*param.tc) == 0) previous_energy = energy;
        // pcout << "previous_energy = " << previous_energy << std::endl;
        // pcout << "energy = " << energy << std::endl;
        
        
        ///// Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        ///// Output.
        if (iT % vtsTout == 0) {
            pcout << "Writing VTS file at iteration = " << iT << std::endl;
            // pcout << "maximal velocity = " << computeMax(*computeVelocityNorm(lattice)) << std::endl;
            writeVTS(lattice, param, iT);
        }

    }
    statsOut.close();


    return 0;
}
