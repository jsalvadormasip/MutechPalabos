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

#include <cmath>
#include <iostream>
#include <tuple>

#include "magic_enum.hpp"
#include "palabos3D.h"
#include "palabos3D.hh"  // explicit inclusion because it contains templates
#include "setup3d.h"
#include "shapes.h"
#include "simulationParameters.cpp"  // explicit inclusion because it contains templates
#include "simulationParameters.h"

#define DESCRIPTOR D3Q27Descriptor
namespace sp = incompressible_simulation_parameters;
using namespace plb;
using namespace std;
using Real = double;
using Int = plint;

auto split_string(const std::string &text, const char delim = ' ') -> std::vector<std::string>
{
    std::string line;
    std::vector<std::string> vec;
    std::stringstream ss(text);
    while (std::getline(ss, line, delim)) {
        vec.push_back(line);
    }
    return vec;
}

// using vectorOfPairsNameValue = std::vector<std::pair<std::string,double>>;
// using vectorOfPairsNamesValue =
// std::vector<std::pair<std::vector<std::string>,double>>;

void writeVTK3d(
    MultiBlockLattice3D<Real, DESCRIPTOR> &lattice, IncomprFlowParam<Real> const &parameters,
    plint iter)
{
    Real dx = 1;  // parameters.getDeltaX();
    Real dt = 1;  // parameters.getDeltaT();
    auto slice = Box3D(lattice.getBoundingBox());
    slice.z0 = util::roundToInt(parameters.getNz() / 2) - 10;
    slice.z1 = util::roundToInt(parameters.getNz() / 2) + 10;
    VtkImageOutput3D<Real> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<3, Real>(*computeVelocity(lattice, slice), "velocity", dx / dt);
    vtkOut.writeData<Real>(*computeDensity(lattice, slice), "rho", dx / dt);
}

template <typename Real, template <typename U> class Descriptor>
auto setupLattice3d(IncomprFlowParam<Real> const &parameters, BCmodel kmin_, BSchemeInfo &xmlParse)
{
    // 1. Create a TriangleSet on the heap and dump it on a stl file
    auto location_ellipsoid = new Array<Real, 3>(
        parameters.getNx() / 3.0, parameters.getNy() / 2.0, parameters.getNz() / 2.0);
    Real radius = parameters.getResolution() / 2.;
    TriangleSet<Real> *triangle_set =
        //        readObstacle(*location_ellipsoid, radius, radius, radius,
        //        0.5*parameters.getResolution());
        generateEllipsoid(*location_ellipsoid, radius, radius, radius, 0.5);
    triangle_set->writeBinarySTL("obstacle.stl");

    // 2. Voxelize the domain
    auto [voxelized_domain, triangle_boundary] =
        voxelize_helper<Real, Descriptor>(parameters, triangle_set);

    // 3. Define the lattice
    auto lattice_ptr =
        new MultiBlockLattice3D<Real, Descriptor>(voxelized_domain->getVoxelMatrix());

    //    ImageWriter<Int> imageWriter("leeloo");

    // 4. Define the dynamics
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelized_domain->getVoxelMatrix());
    //    imageWriter.writeScaledGif("voxels",*extractSubDomain(flagMatrix,flagMatrix.getBoundingBox()));
    if (parameters.getRe() >= 2000.0) {
        defineDynamics(
            *lattice_ptr, flagMatrix, lattice_ptr->getBoundingBox(),
            new ConsistentSmagorinskyCompleteRegularizedBGKdynamics<Real, Descriptor>(
                parameters.getOmega(), 0.1),
            voxelFlag::outside);
        pcout << "Using Smagorinsky BGK dynamics." << std::endl;
    } else {
        defineDynamics(
            *lattice_ptr, lattice_ptr->getBoundingBox(),
            new BGKdynamics<Real, Descriptor>(parameters.getOmega()));
        pcout << "Using BGK dynamics." << std::endl;
    }

    // 5. Set sponge zones to reduce pressure waves reflections
    createSpongeZones(parameters, Array<plint, 6>(0, 20, 0, 0, 0, 0), lattice_ptr);

    // 6. Define outer boundary conditions
    OnLatticeBoundaryCondition3D<Real, Descriptor> *onlatt_boundary_condition =
        inject_on_lattice_bc(lattice_ptr, parameters, *location_ellipsoid);

    // 7. Define inner-offlattice boundary conditions
    OffLatticeBoundaryCondition3D<Real, Descriptor, Array<Real, 3>> *offlatt_boundary_condition =
        inject_off_lattice_bc(lattice_ptr, voxelized_domain, kmin_, xmlParse);

    // return lattice and boundary conditions
    return std::tuple {lattice_ptr, onlatt_boundary_condition, offlatt_boundary_condition};
}

int main(int argc, char *argv[])
{
    // Palabos initialization
    plbInit(&argc, &argv);
    string outdir = "tmp/";
    BCmodel kmin = BCmodel::LIgeneric;

    BSchemeInfo xmlParse("LIplusGeneric.xml");

    if (argc > 3) {
        cout << "Error! too many parameters" << endl;
        abort();
    } else if (argc == 3) {
        outdir = static_cast<string>(argv[1]);
        auto kcasts = magic_enum::enum_cast<BCmodel>(argv[2]);
        if (kcasts.has_value())
            kmin = kcasts.value();
        else {
            auto kcastf = magic_enum::enum_cast<BCmodel>(util::roundToInt(stof(argv[2])));
            if (kcastf.has_value()) {
                kmin = kcastf.value();
            } else {
                pcout << "Kmin value" << argv[2] << " wrong or not implemented, aborting..."
                      << std::endl;
                abort();
            }
        }
    } else if (argc == 2)
        outdir = static_cast<string>(argv[1]);
    else
        outdir = "tmp";
    outdir.append("/");
    if (global::mpi().getRank() == 0)
        system(("mkdir -p " + outdir).c_str());
    global::directories().setOutputDir(outdir);
    pcout << "The selected model is " << magic_enum::enum_name(kmin) << std::endl;
    // Define simulation parameters: we use two ad-hoc units helper
    const Real re = 20;  // NB: Obstacle reynolds!
    sp::Numerics<Real, Int> lu;
    sp::NonDimensional<Real> dimless;
    dimless.initReLxLyLz(re, 12, 7, 7);  // the diameter is the reference length
    dimless.printParameters();
    // initLrefluNodim initializes lu, getIncomprFlowParam() returns
    // the palabos structure IncomprFlowParam
    IncomprFlowParam<Real> parameters =
        lu.initLrefluNodim(25 /*resolution*/, &dimless, 0.04 /*u_lb*/, false /*tau, optional*/)
            .printParameters()
            .getIncomprFlowParam();

    // Define output constants
    const Real logT = (Real)0.02;
    [[maybe_unused]] const Real imSave = (Real)0.06;
    [[maybe_unused]] const Real vtkSave =
        10 * (Real)parameters.getDeltaT();  //(Real)parameters.getDeltaT();
    const Real maxT = (Real)5.1;
    writeLogFile(parameters, "Poiseuille flow");

    // Setup the simulation and allocate lattice and boundary conditions
    auto [lattice_ptr, boundaryon, boundaryoff] =
        setupLattice3d<Real, DESCRIPTOR>(parameters, kmin, xmlParse);

    auto &lattice = *lattice_ptr;

    pcout << "Saving VTK file at initialization..." << endl;
    writeVTK3d(lattice, parameters, 0);

    // Main loop over time iterations.
    for (plint iT = 0; (Real)iT * parameters.getDeltaT() < maxT; ++iT) {
        Real cd = 0., error = 0.0;
        if (iT % parameters.nStep(logT) == 0) {
            cd = 2.0 * boundaryoff->getForceOnObject()[0]
                 / (0.25 * lu.getUlb() * lu.getUlb() * lu.getLref() * lu.getLref() * M_PI);
            error = (cd - empirical_sphere_drag(parameters.getRe()))
                    / empirical_sphere_drag(parameters.getRe());
            std::string fullName = global::directories().getLogOutDir() + "Cd.dat";
            plb_ofstream ofile(fullName.c_str(), std::ostream::app);
            ofile << std::setw(15) << iT * parameters.getDeltaT() << " " << std::setw(15) << cd
                  << " " << std::setw(15) << error << std::endl;
        }
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages
        //   (getStoredAverageEnergy and getStoredAverageDensity) correspond to
        //   the previous time iT-1.

        // For the 3D case we use the VTK output
        if (iT % parameters.nStep(vtkSave) == 0 && iT > 0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK3d(lattice, parameters, iT);
        }

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "step " << iT << "; t=" << iT * parameters.getDeltaT();
        }
        global::timer("iteration").restart();

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time
        //   iT.
        if (iT % parameters.nStep(logT) == 0) {
            pcout << "; av energy = " << setprecision(10) << getStoredAverageEnergy<Real>(lattice)
                  << "; av rho = " << getStoredAverageDensity<Real>(lattice) << "; Cd = " << cd
                  << "; Drag sphere error = " << error << std::endl;
        }
        global::timer("iteration").stop();
    }

    delete boundaryoff;
    delete boundaryon;
    delete lattice_ptr;

    return 0;
}
