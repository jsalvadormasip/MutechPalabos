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
 * A fluid constrained between a hot bottom wall (no-slip for the velocity) and a cold
 * top wall (no-slip for the velocity). The lateral walls are periodic. Under the
 * influence of gravity, convection rolls are formed. Thermal effects are modelled
 * by means of a Boussinesq approximation: the fluid is incompressible, and the influence
 * of the temperature is visible only through a body-force term, representing buoyancy
 * effects. The temperature field obeys an advection-diffusion equation.
 *
 * The simulation is first created in a fully symmetric manner. The symmetry is therefore
 * not spontaneously broken; while the temperature drops linearly between the hot and
 * and cold wall, the convection rolls fail to appear at this point. In a second stage, a
 * random noise is added to trigger the instability.
 *
 * This application is technically a bit more advanced than the other ones, because it
 * illustrates the concept of data processors. In the present case, they are used to
 * create the initial condition, and to trigger the instability.
 **/

#include <cstdlib>
#include <iostream>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define ADESCRIPTOR  descriptors::AdvectionDiffusionD3Q7Descriptor

#define ADYNAMICS  AdvectionDiffusionBGKdynamics
#define NSDYNAMICS GuoExternalForceBGKdynamics

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(
    plint iY, RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> const &parameters)
{
    T y = (T)iY / parameters.getResolution();
    return 4. * parameters.getLatticeU() * (y - y * y);
}

/// A functional, used to initialize the velocity for the boundary conditions
template <typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> parameters_) :
        parameters(parameters_)
    { }
    void operator()(
        [[maybe_unused]] plint iX, [[maybe_unused]] plint iY, plint iZ, Array<T, 3> &u) const
    {
        u[0] = poiseuilleVelocity(iZ, parameters);
        u[1] = T();
        u[2] = T();
    }

private:
    RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> parameters;
};

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template <typename T>
class OneDensityAndPoiseuilleVelocity {
public:
    OneDensityAndPoiseuilleVelocity(
        RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> parameters_) :
        parameters(parameters_)
    { }
    void operator()(
        [[maybe_unused]] plint iX, [[maybe_unused]] plint iY, plint iZ, T &rho,
        Array<T, 3> &u) const
    {
        rho = (T)1;
        u[0] = poiseuilleVelocity(iZ, parameters);
        u[1] = T();
        u[2] = T();
    }

private:
    RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> parameters;
};

/// Initialization of the temperature field.
template <
    typename T, template <typename NSU> class nsDescriptor,
    template <typename ADU> class adDescriptor>
struct IniTemperatureRayleighBenardProcessor3D :
    public BoxProcessingFunctional3D_L<T, adDescriptor> {
    IniTemperatureRayleighBenardProcessor3D(
        RayleighBenardFlowParam<T, nsDescriptor, adDescriptor> parameters_) :
        parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T, adDescriptor> &adLattice)
    {
        Dot3D absoluteOffset = adLattice.getLocation();

        plint meanRegion = parameters.getNx() / 8;
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absoluteX = absoluteOffset.x + iX;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteZ = absoluteOffset.z + iZ;

                    T temperature = T();
                    if (absoluteX < meanRegion || absoluteX > 7 * meanRegion) {
                        temperature =
                            (parameters.getHotTemperature() + parameters.getColdTemperature())
                            / (T)2;
                    } else {
                        temperature = parameters.getHotTemperature()
                                      - parameters.getDeltaTemperature()
                                            / (T)(parameters.getNz() - 1) * (T)absoluteZ;
                    }

                    Array<T, adDescriptor<T>::d> jEq(0., 0., 0.);
                    adLattice.get(iX, iY, iZ).defineDensity(temperature);
                    iniCellAtEquilibrium(adLattice.get(iX, iY, iZ), temperature, jEq);
                }
            }
        }
    }
    virtual IniTemperatureRayleighBenardProcessor3D<T, nsDescriptor, adDescriptor> *clone() const
    {
        return new IniTemperatureRayleighBenardProcessor3D<T, nsDescriptor, adDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }

private:
    RayleighBenardFlowParam<T, nsDescriptor, adDescriptor> parameters;
};

/// Perturbation of the temperature field to instantiate the instability.
template <
    typename T, template <typename NSU> class nsDescriptor,
    template <typename ADU> class adDescriptor>
struct PerturbTemperatureRayleighBenardProcessor3D :
    public BoxProcessingFunctional3D_L<T, adDescriptor> {
    PerturbTemperatureRayleighBenardProcessor3D(
        RayleighBenardFlowParam<T, nsDescriptor, adDescriptor> parameters_) :
        parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T, adDescriptor> &lattice)
    {
        Dot3D absoluteOffset = lattice.getLocation();

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteX = absoluteOffset.x + iX;
                    plint absoluteY = absoluteOffset.y + iY;
                    plint absoluteZ = absoluteOffset.z + iZ;

                    if ((absoluteX == (parameters.getNx() - 1) / 8 + 2)
                        && (absoluteY == (parameters.getNy() - 1) / 2) && (absoluteZ == 0))
                    {
                        T temperature = T();
                        temperature = parameters.getHotTemperature() * 1.1;

                        Array<T, adDescriptor<T>::d> jEq(0., 0., 0.);
                        lattice.get(iX, iY, iZ).defineDensity(temperature);
                        iniCellAtEquilibrium(lattice.get(iX, iY, iZ), temperature, jEq);
                    }
                }
            }
        }
    }
    virtual PerturbTemperatureRayleighBenardProcessor3D<T, nsDescriptor, adDescriptor> *clone()
        const
    {
        return new PerturbTemperatureRayleighBenardProcessor3D<T, nsDescriptor, adDescriptor>(
            *this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    RayleighBenardFlowParam<T, nsDescriptor, adDescriptor> parameters;
};

void rayleighBenardSetup(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiBlockLattice3D<T, ADESCRIPTOR> &adLattice,
    OnLatticeBoundaryCondition3D<T, NSDESCRIPTOR> &nsBoundaryCondition,
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, ADESCRIPTOR> &adBoundaryCondition,
    RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> &parameters)
{
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    Box3D bottom(1, nx - 2, 0, ny - 1, 0, 0);
    Box3D top(1, nx - 2, 0, ny - 1, nz - 1, nz - 1);
    Box3D inlet(0, 0, 0, ny - 1, 1, nz - 2);
    Box3D inletNN(0, 0, 0, ny - 1, 0, 0);
    Box3D inletPN(0, 0, 0, ny - 1, nz - 1, nz - 1);
    Box3D outlet(nx - 1, nx - 1, 0, ny - 1, 1, nz - 2);
    Box3D outletPP(nx - 1, nx - 1, 0, ny - 1, nz - 1, nz - 1);
    Box3D outletNP(nx - 1, nx - 1, 0, ny - 1, 0, 0);

    nsBoundaryCondition.addVelocityBoundary2N(bottom, nsLattice);
    nsBoundaryCondition.addVelocityBoundary2P(top, nsLattice);
    nsBoundaryCondition.addVelocityBoundary0N(inlet, nsLattice);
    nsBoundaryCondition.addVelocityBoundary0P(outlet, nsLattice, boundary::neumann);
    nsBoundaryCondition.addExternalVelocityEdge1NN(inletNN, nsLattice);
    nsBoundaryCondition.addExternalVelocityEdge1PN(inletPN, nsLattice);
    nsBoundaryCondition.addExternalVelocityEdge1PP(outletPP, nsLattice);
    nsBoundaryCondition.addExternalVelocityEdge1NP(outletNP, nsLattice);

    setBoundaryVelocity(nsLattice, nsLattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters));
    initializeAtEquilibrium(
        nsLattice, nsLattice.getBoundingBox(), OneDensityAndPoiseuilleVelocity<T>(parameters));

    adBoundaryCondition.addTemperatureBoundary2N(bottom, adLattice);
    adBoundaryCondition.addTemperatureBoundary2P(top, adLattice);
    adBoundaryCondition.addTemperatureBoundary0N(inlet, adLattice);
    adBoundaryCondition.addTemperatureBoundary0P(outlet, adLattice, boundary::neumann);
    adBoundaryCondition.addTemperatureEdge1NN(inletNN, adLattice);
    adBoundaryCondition.addTemperatureEdge1PN(inletPN, adLattice);
    adBoundaryCondition.addTemperatureEdge1PP(outletPP, adLattice);
    adBoundaryCondition.addTemperatureEdge1NP(outletNP, adLattice);

    // adBoundaryCondition.addTemperatureBoundary2N(bottom, adLattice);
    // adBoundaryCondition.addTemperatureBoundary2P(top, adLattice);

    applyProcessingFunctional(
        new IniTemperatureRayleighBenardProcessor3D<T, NSDESCRIPTOR, ADESCRIPTOR>(parameters),
        adLattice.getBoundingBox(), adLattice);

    nsLattice.initialize();
    adLattice.initialize();
}

void writeVTK(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiBlockLattice3D<T, ADESCRIPTOR> &adLattice,
    RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> const &parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(nsLattice), "velocityNorm", dx / dt);
    vtkOut.writeData<3, float>(*computeVelocity(nsLattice), "velocity", dx / dt);
    // Temperature is the order-0 moment of the advection-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    vtkOut.writeData<float>(*computeDensity(adLattice), "temperature", (T)1);
}

void writeGif(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiBlockLattice3D<T, ADESCRIPTOR> &adLattice,
    int iT)
{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    const plint nz = nsLattice.getNz();
    Box3D slice(0, nx - 1, (ny - 1) / 2, (ny - 1) / 2, 0, nz - 1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(
        createFileName("u", iT, 6), *computeVelocityNorm(nsLattice, slice), imSize, imSize);
    // Temperature is the order-0 moment of the advection-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(
        createFileName("temperature", iT, 6), *computeDensity(adLattice, slice), imSize, imSize);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    global::timer("simTime").start();

    T Ra = 0.;
    try {
        global::argv(1).read(Ra);
    } catch (PlbIOException &exception) {
        pcout << exception.what() << endl;
        pcout << "The structure of the input parameters should be : " << (string)global::argv(0)
              << " Ra" << endl;
        exit(1);
    }

    const T lx = 5.0;
    const T ly = 4.0;
    const T lz = 1.0;
    const T uMax = 0.1;
    const T Pr = 1.0;

    const T hotTemperature = 1.0;
    const T coldTemperature = 0.0;
    const plint resolution = 20;

    global::directories().setOutputDir("./tmp/");

    RayleighBenardFlowParam<T, NSDESCRIPTOR, ADESCRIPTOR> parameters(
        Ra, Pr, uMax, coldTemperature, hotTemperature, resolution, lx, ly, lz);

    writeLogFile(parameters, "palabos.log");

    const double rayleigh = parameters.getResolution() * parameters.getResolution()
                            * parameters.getResolution() * parameters.getDeltaTemperature()
                            * parameters.getLatticeGravity()
                            / (parameters.getLatticeNu() * parameters.getLatticeKappa());

    const double prandtl = parameters.getLatticeNu() / parameters.getLatticeKappa();

    pcout << rayleigh << " " << prandtl << endl;

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    T nsOmega = parameters.getSolventOmega();
    T adOmega = parameters.getTemperatureOmega();

    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice(
        nx, ny, nz, new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega));
    // Use periodic boundary conditions.
    nsLattice.periodicity().toggleAll(true);

    MultiBlockLattice3D<T, ADESCRIPTOR> adLattice(
        nx, ny, nz, new ADYNAMICS<T, ADESCRIPTOR>(adOmega));
    // Use periodic boundary conditions.
    adLattice.periodicity().toggleAll(true);

    OnLatticeBoundaryCondition3D<T, NSDESCRIPTOR> *nsBoundaryCondition =
        createLocalBoundaryCondition3D<T, NSDESCRIPTOR>();

    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, ADESCRIPTOR> *adBoundaryCondition =
        createLocalAdvectionDiffusionBoundaryCondition3D<T, ADESCRIPTOR>();

    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);

    rayleighBenardSetup(
        nsLattice, adLattice, *nsBoundaryCondition, *adBoundaryCondition, parameters);

    Array<T, NSDESCRIPTOR<T>::d> forceOrientation(T(), T(), (T)1);
    plint processorLevel = 1;
    integrateProcessingFunctional(
        new BoussinesqThermalProcessor3D<T, NSDESCRIPTOR, ADESCRIPTOR>(
            parameters.getLatticeGravity(), parameters.getAverageTemperature(),
            parameters.getDeltaTemperature(), forceOrientation),
        nsLattice.getBoundingBox(), nsLattice, adLattice, processorLevel);

#ifndef PLB_REGRESSION
    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for rayleighBenardSetup:" << tIni << endl;
#endif
    global::timer("simTime").start();

#ifndef PLB_REGRESSION
    plint evalTime = 10000;
#endif
    plint iT = 0;
    plint maxT = 1000000000;
    plint statIter = 10;
#ifndef PLB_REGRESSION
    plint saveIter = 1000;
#endif
    util::ValueTracer<T> converge((T)1, (T)100, 1.0e-3);
    bool convergedOnce = false;

    // Main loop over time iterations.
    for (iT = 0; iT <= maxT; ++iT) {
#ifndef PLB_REGRESSION
        if (iT == (evalTime)) {
            T tEval = global::timer("simTime").stop();
            T remainTime = (tEval - tIni) / (T)evalTime * (T)maxT / (T)3600;
            global::timer("simTime").start();
            pcout << "Remaining " << (plint)remainTime << " hours, and ";
            pcout << (plint)((T)60 * (remainTime - (T)((plint)remainTime)) + 0.5) << " minutes."
                  << endl;
        }
#endif
        if (iT % statIter == 0) {
            int zDirection = 2;
            T nusselt = computeNusseltNumber(
                nsLattice, adLattice, nsLattice.getBoundingBox(), zDirection,
                parameters.getDeltaX(), parameters.getLatticeKappa(),
                parameters.getDeltaTemperature());
            converge.takeValue(nusselt, true);
        }
        if (converge.hasConverged()) {
            if (!convergedOnce) {
                convergedOnce = true;
                converge.resetValues();
                converge.setEpsilon(1.0e-11);
                applyProcessingFunctional(
                    new PerturbTemperatureRayleighBenardProcessor3D<T, NSDESCRIPTOR, ADESCRIPTOR>(
                        parameters),
                    adLattice.getBoundingBox(), adLattice);
                pcout << "Intermetiate convergence.\n";
#ifndef PLB_REGRESSION
                pcout << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
                writeVTK(nsLattice, adLattice, parameters, iT);
#endif
            } else {
                pcout << "Simulation is over.\n";
                break;
            }
        }
#ifndef PLB_REGRESSION
        if (iT % saveIter == 0) {
            pcout << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
            writeVTK(nsLattice, adLattice, parameters, iT);

            pcout << iT << " : Writing gif." << endl;
            writeGif(nsLattice, adLattice, iT);
        }
#endif

        // Lattice Boltzmann iteration step.
        adLattice.collideAndStream();
        nsLattice.collideAndStream();
    }

#ifndef PLB_REGRESSION
    writeGif(nsLattice, adLattice, iT);

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd - tIni;
    T nx100 = nsLattice.getNx() / (T)100;
    T ny100 = nsLattice.getNy() / (T)100;
    T nz100 = nsLattice.getNz() / (T)100;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << nx100 * ny100 * nz100 * (T)iT / totalTime << endl;
#endif

    delete nsBoundaryCondition;
    delete adBoundaryCondition;
}
