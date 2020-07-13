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

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor
#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

template<typename T>
class CylinderShapeDomain3D : public plb::DomainFunctional3D {
public:
    CylinderShapeDomain3D(plb::plint cx_, plb::plint cz_, plb::plint radius)
        : cx(cx_),
        cz(cz_),
        radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY, plb::plint iZ) const {
        return plb::util::sqr(iX - cx) + plb::util::sqr(iZ - cz) <= radiusSqr;
    }
    virtual CylinderShapeDomain3D<T>* clone() const {
        return new CylinderShapeDomain3D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cz;
    plb::plint radiusSqr;
};

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu

    T deltaP = - (alpha * nu);

    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }

    sum *= ((T)4 * alpha * a *a /std::pow(pi,(T)3));
    sum += (alpha / (T)2 * (x * x - a*a / (T)4));
    
    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

void simulationSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                      IncomprFlowParam<T> const& parameters,
                      OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                      Array<plint,3> &forceIds )
{
    // No periodic boundaries
    lattice.periodicity().toggleAll(false);

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    // Flow direction: z
    // Top/Bottom direction: y
    Box3D top    = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1); // Full Area
    Box3D bottom = Box3D(0, nx-1, 0, 0, 0, nz-1); // Full Area
    
    Box3D inlet  = Box3D(0, nx-1, 0, ny-1, 0, 0); // Full Area
    Box3D outlet = Box3D(1, nx-2, 1, ny-2, nz-1, nz-1); // Offset from wall boundaries by 1 lattice unit
    
    Box3D left   = Box3D(0, 0, 0, ny-1, 0, nz-1); // Full Area
    Box3D right  = Box3D(nx-1, nx-1, 0, ny-1, 0, nz-1); // Full Area
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
    boundaryCondition.addVelocityBoundary2P(outlet, lattice, boundary::neumann);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, top);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, bottom);
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, left);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, right);
    
    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));

    setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));

    // Add the obstacle: cylinder 3d
    plint cx     = nx/2 + 2; // cx is slightly offset to avoid full symmetry,
                             // and to get a Von Karman Vortex street.
    plint cz     = nz/4;
    plint radius = parameters.getResolution() / 2; // the diameter is the reference length

    lattice.toggleInternalStatistics(true);
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();
    forceIds[2] = lattice.internalStatSubscription().subscribeSum();

    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain3D<T>(cx, cz, radius),
                   new plb::MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
    initializeMomentumExchange(lattice, lattice.getBoundingBox());

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
    IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const plint xComponent = 2;

    Box3D slice(0, nx - 1, ny / 2, ny / 2, 0, nz - 1);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif(createFileName("uz", iter, 6),
        *computeVelocityComponent(lattice, slice, zComponent),
        imSize, imSize);
    imageWriter.writeScaledGif(createFileName("uNorm", iter, 6),
        *computeVelocityNorm(lattice, slice),
        imSize, imSize);
    imageWriter.writeScaledGif(createFileName("omega", iter, 6),
        *computeNorm(*computeVorticity(
            *computeVelocity(lattice)), slice),
        imSize, imSize);
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

#ifdef MSVC
	VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
#else
	ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 3, dx);
#endif

    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    if (argc != 2)
    {
        pcout << "Error: the parameters are wrong. Give resolution.";
        exit(1);
    }

    T W_ = 21.0;
    T L_ = 75.0;
    T h_ = 6.0;
    T D_ = 6.0;

    // Use the class IncomprFlowParam to convert from
    // dimensionless variables to lattice units, in the
    // context of incompressible flows.
    IncomprFlowParam<T> parameters(
        0.01, // Reference velocity (the maximum velocity in the Poiseuille profile) in lattice units.
        100., // Reynolds number
        atoi(argv[1]), // Resolution of the reference length (cylinder diameter)
        W_/D_, // dimensionless: channel lateral length
        h_/D_, // dimensionless: channel height
        L_/D_ // dimensionless: channel length
    );
    
    const T vtkSave = (T) 0.1; // Time intervals at which to save GIF VTKs, in dimensionless time units
    const T maxT    = (T)10.0; // Total simulation time, in dimensionless time units

    pcout << "omega= " << parameters.getOmega() << std::endl;
    writeLogFile(parameters, "3D square Poiseuille with Cylinder as an obstacle");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new DYNAMICS );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    Array<plint,3> forceIds;
    simulationSetup(lattice, parameters, *boundaryCondition, forceIds);

    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT)
    {
        if (iT%parameters.nStep(vtkSave)==0 && iT>0)
        {
            pcout << "step " << iT << "; t=" << iT*parameters.getDeltaT() << std::endl;
            writeVTK(lattice, parameters, iT);
        }

        // Execute a time iteration.
        lattice.collideAndStream();
    }
    
    delete boundaryCondition;
}
