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

/* \file
 * External flow around a 3D obstacle.
 * This example demonstrates many features of Palabos:
 * Loading geometries from STL files.
 * Using the voxelizer.
 * Using off-lattice boundary conditions.
 * Imposing sophisticated outflow boundary conditions.
 * Using sponge zones.
 * Imposing time dependent inlet boundary conditions.
 * Integrating particles for the computation of streamlines.
 * Computing the force on objects.
 * */

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T, 3> Velocity;
#define DESCRIPTOR descriptors::D3Q19Descriptor
typedef DenseParticleField3D<T, DESCRIPTOR> ParticleFieldT;

#define PADDING 8
T angle = 7.8 * 3.1415 / 180;
static std::string outputDir("./tmp/");

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param {
    T nu;                         // Kinematic viscosity.
    T lx, ly, lz;                 // Size of computational domain, in physical units.
    T cx, cy, cz;                 // Position of the center of the obstacle, in physical units.
    plint cxLB, cyLB, czLB;       // Position of the center of the obstacle, in lattice units.
    bool freeSlipWall;            // Use free-slip condition on obstacle, as opposed to no-slip?
    bool lateralFreeSlip;         // Use free-slip lateral boundaries or periodic ones?
    T maxT, statT, imageT, vtkT;  // Time, in physical units, at which events occur.
    plint resolution;             // Number of lattice nodes along a reference length.
    T inletVelocity;              // Inlet x-velocity in physical units.
    T uLB;                        // Velocity in lattice units (numerical parameters).
    bool useSmago;                // Use a Smagorinsky LES model or not.
    T cSmago;                     // Parameter for the Smagorinsky LES model.
    plint nx, ny, nz;             // Grid resolution of bounding box.
    T omega;                      // Relaxation parameter.
    T dx, dt;                     // Discrete space and time steps.
    plint maxIter, statIter;      // Time for events in lattice units.
    plint imageIter, vtkIter;
    bool useParticles;       // Simulate particles or not.
    int particleTimeFactor;  // If the particle time factor is 2, then the integration time step
                             //   for the particles is twice that of the fluid.
    T particleProbabilityPerCell;  // Probability of injection of a particle at an injection cell at
                                   // each time step.
    T cutOffSpeedSqr;              // Criterion to eliminate particles with very small velocity.
    int maxNumParticlesToWrite;    // Maximum number of particles in the output VTK files.

    T outletSpongeZoneWidth;     // Width of the outlet sponge zone.
    plint numOutletSpongeCells;  // Number of the lattice nodes contained in the outlet sponge zone.
    int outletSpongeZoneType;    // Type of the outlet sponge zone (Viscosity or Smagorinsky).
    T targetSpongeCSmago;  // Target Smagorinsky parameter at the end of the Smagorinsky sponge
                           // Zone.
    plint initialIter;  // Number of initial iterations until the inlet velocity reaches its final
                        // value.

    Box3D inlet, outlet, lateral1;  // Outer domain boundaries in lattice units.
    Box3D lateral2, lateral3, lateral4;

    std::string geometry_fname;

    Param() { }

    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);
        document["geometry"]["filename"].read(geometry_fname);
        document["geometry"]["center"]["x"].read(cx);
        document["geometry"]["center"]["y"].read(cy);
        document["geometry"]["center"]["z"].read(cz);
        document["geometry"]["freeSlipWall"].read(freeSlipWall);
        document["geometry"]["lateralFreeSlip"].read(lateralFreeSlip);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

        document["numerics"]["nu"].read(nu);
        document["numerics"]["inletVelocity"].read(inletVelocity);
        document["numerics"]["resolution"].read(resolution);
        document["numerics"]["uLB"].read(uLB);  //velocity in lattice boltzmann units. 
        document["numerics"]["useSmago"].read(useSmago);
        if (useSmago) {
            document["numerics"]["cSmago"].read(cSmago);
        }

        document["numerics"]["useParticles"].read(useParticles);
        if (useParticles) {
            document["numerics"]["particleTimeFactor"].read(particleTimeFactor);
            document["numerics"]["particleProbabilityPerCell"].read(particleProbabilityPerCell);
            document["numerics"]["cutOffSpeedSqr"].read(cutOffSpeedSqr);
            document["numerics"]["maxNumParticlesToWrite"].read(maxNumParticlesToWrite);
        }

        document["numerics"]["outletSpongeZoneWidth"].read(outletSpongeZoneWidth);
        std::string zoneType;
        document["numerics"]["outletSpongeZoneType"].read(zoneType);
        if ((util::tolower(zoneType)).compare("viscosity") == 0) {
            outletSpongeZoneType = 0;
        } else if ((util::tolower(zoneType)).compare("smagorinsky") == 0) {
            outletSpongeZoneType = 1;
        } else {
            pcout << "The sponge zone type must be either \"Viscosity\" or \"Smagorinsky\"."
                  << std::endl;
            exit(-1);
        }
        document["numerics"]["targetSpongeCSmago"].read(targetSpongeCSmago);

        document["numerics"]["initialIter"].read(initialIter);

        document["output"]["maxT"].read(maxT);
        document["output"]["statT"].read(statT);
        document["output"]["imageT"].read(imageT);
        document["output"]["vtkT"].read(vtkT);

        computeLBparameters();
    }

    void computeLBparameters()   //maybe LB parameters need to be normalized or smthing? interesting to figure out. 
    {
        dx = lx / (resolution - 1.0);   //as said, x and y step for every particle.
        dt = (uLB / inletVelocity) * dx; //time step
        T nuLB = nu * dt / (dx * dx);
        omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nuLB + 0.5);
        if (lateralFreeSlip) {  //if false boundaries are periodic. Then, if true, they are not, thus the y and z have 1 less row for some reason
            nx = util::roundToInt(lx / dx) + 1;
            ny = util::roundToInt(ly / dx) + 1;
            nz = util::roundToInt(lz / dx) + 1; 
        } else {
            nx = util::roundToInt(lx / dx) + 1;
            ny = util::roundToInt(ly / dx);
            nz = util::roundToInt(lz / dx);  
        }
        cxLB = util::roundToInt(cx / dx);
        cyLB = util::roundToInt(cy / dx);
        czLB = util::roundToInt(cz / dx);
        maxIter = util::roundToInt(maxT / dt);
        statIter = util::roundToInt(statT / dt);
        imageIter = util::roundToInt(imageT / dt);
        vtkIter = util::roundToInt(vtkT / dt);
        numOutletSpongeCells = util::roundToInt(outletSpongeZoneWidth / dx);

        inlet = Box3D(0, 0, 0, ny - 1, 0, nz - 1);  //this means x starts at 0 and ends at 0, y starts at 0 and ends at ny-1, and z starts at 0 and ends at nz-1
        outlet = Box3D(nx - 1, nx - 1, 0, ny - 1, 0, nz - 1);
        lateral1 = Box3D(1, nx - 2, 0, 0, 0, nz - 1); //this means x starts at 1 and ends at nx-2, y starts and ends at 0 and z goes from 0 to nz-1
        lateral2 = Box3D(1, nx - 2, ny - 1, ny - 1, 0, nz - 1);
        lateral3 = Box3D(1, nx - 2, 1, ny - 2, 0, 0);
        lateral4 = Box3D(1, nx - 2, 1, ny - 2, nz - 1, nz - 1);
    }

    Box3D boundingBox() const
    {
        return Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1);  //definition of the bounding box.
    }

    T getInletVelocity(plint iIter)
    {
        static T pi = std::acos((T)-1.0);

        if (iIter >= initialIter) {
            return uLB;
        }

        if (iIter < 0) {
            iIter = 0;
        }

        return uLB;// * std::sin(pi * iIter / (2.0 * initialIter));  //This produces a sinusoidal velocity profile that smoothly ramps up from 0 to uLB over the range of iterations from 0 to initialIter.
    }  //basically this function ramps the inlet velocity to uLB following a sinusoidal distribution, reaching its peak, uLB when initalIter is achieved. Interesting...
};

Param param; //creates object param of structure Param. 

// Instantiate the boundary conditions for the outer domain.
void outerDomainBoundaries(
    MultiBlockLattice3D<T, DESCRIPTOR> *lattice, MultiScalarField3D<T> *rhoBar,
    MultiTensorField3D<T, 3> *j, OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc)
{
    Array<T, 3> uBoundary(param.getInletVelocity(0)*std::cos(angle), 0.0, param.getInletVelocity(0)*std::sin(angle));  //because the inlet is x direction, the inlet velocity you get x direction.
    
    if (param.lateralFreeSlip) {  //if there are free slip lateral boundaries, the boundaries are not periodic, and then the inlet is given dirichlet velocity and the rest free slip. 
        pcout << "Free-slip lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(false);  //they do the little shitty arrow thing because we actually called the pointers in the function for some reason, and now we want to access the data. 
        rhoBar->periodicity().toggleAll(false);
        j->periodicity().toggleAll(false);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.inlet, uBoundary);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral1, boundary::freeslip);  //this setvelocityconditiononblockboundaries function is awesome. You just give it the lattice pointer, slice, and the general velocity vector and it cooks.
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral2, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::freeslip);
        setBoundaryVelocity(*lattice, param.lateral1, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral2, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral3, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral4, uBoundary);

        // The VirtualOutlet is a sophisticated outflow boundary condition.
        Box3D globalDomain(lattice->getBoundingBox());  //Box3D is a class or structure that makes a 3d box with the coordinates that you gve it. Normally it uses x0, x1, y0, y1, z0, z1. However, now getboundingbox already gives that. 
        
        std::vector<MultiBlock3D *> bcargs; //I think bcargs stores the pointers of the different objects. 
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);
        T outsideDensity = 1.0;
        int bcType = 1;
        integrateProcessingFunctional(   // so the virtual outlet is an outflow boundary condition that works with external rhoBar from previous step.
            new VirtualOutlet<T, DESCRIPTOR>(outsideDensity, globalDomain, bcType), param.outlet,
            bcargs, 2);  //bascially the integrateProcessingFunctional puts an internal data processor
        setBoundaryVelocity(*lattice, param.outlet, uBoundary);
    } else {
        pcout << "Periodic lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(true);
        rhoBar->periodicity().toggleAll(true);
        j->periodicity().toggleAll(true);

        lattice->periodicity().toggle(0, false);  //the x direction is still not periodic. Only the laterals. 
        rhoBar->periodicity().toggle(0, false);
        j->periodicity().toggle(0, false);
        lattice->periodicity().toggle(2, false);  //the z direction is still not periodic. Only the laterals. 
        rhoBar->periodicity().toggle(2, false);
        j->periodicity().toggle(2, false);

        T outsideDensity = 1.0;
        // bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        // setBoundaryVelocity(*lattice, param.inlet, uBoundary);
        // bc->addVelocityBoundary0N(param.inlet, *lattice);  //adds velocity boundary to the inlet side, of the lattice, does not specifiy boundary type
        // setBoundaryVelocity(*lattice, param.inlet, uBoundary);  //again magic line is awesome
        // bc->addPressureBoundary0N(param.inlet,*lattice);
        // setBoundaryDensity(*lattice,param.inlet,outsideDensity);
        // bc->addVelocityBoundary2P(param.lateral4,*lattice);  
        // setBoundaryVelocity(*lattice,param.lateral4,uBoundary);
        // bc->addPressureBoundary2P(param.lateral4,*lattice);
        // setBoundaryDensity(*lattice, param.lateral4, outsideDensity);   //I ASSUMED DENSITY HERE IS 1.225. I AM NOT SURE, BECAUSE I DERIVED IT FROM 101325 Pa pressure 
        // bc->addVelocityBoundary2N(param.lateral3,*lattice);
        // setBoundaryVelocity(*lattice,param.lateral3,uBoundary);
        // bc->addPressureBoundary2N(param.lateral3,*lattice);
        // setBoundaryDensity(*lattice, param.lateral3, outsideDensity);   //I ASSUMED DENSITY HERE IS 1.225. I AM NOT SURE, BECAUSE I DERIVED IT FROM 101325 Pa pressure 
        // bc->addVelocityBoundary0P(param.outlet,*lattice);
        // setBoundaryVelocity(*lattice,param.outlet,uBoundary);
        // bc->addPressureBoundary0P(param.outlet,*lattice);
        // setBoundaryDensity(*lattice,param.outlet,outsideDensity);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.inlet, uBoundary);
        // bc->setPressureConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        // setBoundaryDensity(*lattice,param.inlet,outsideDensity);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.outlet, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.outlet, uBoundary);
        // bc->setPressureConditionOnBlockBoundaries(*lattice, param.outlet, boundary::dirichlet);
        // setBoundaryDensity(*lattice,param.outlet,outsideDensity);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.lateral3, uBoundary);
        // bc->setPressureConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::dirichlet);
        // setBoundaryDensity(*lattice,param.lateral3,outsideDensity);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.lateral4, uBoundary);
        // bc->setPressureConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::dirichlet);
        // setBoundaryDensity(*lattice,param.lateral4,outsideDensity);
        // The VirtualOutlet is a sophisticated outflow boundary condition.
        // The "globalDomain" argument for the boundary condition must be
        // bigger than the actual bounding box of the simulation for
        // the directions which are periodic.   I don't understand why.
        Box3D globalDomain(lattice->getBoundingBox());
        //pcout << "Look here chaval" << lattice->getBoundingBox().toString() << std::endl;
        globalDomain.y0 -= 2;  // y-periodicity   I don't really understand how this makes it periodic but sure. 
        globalDomain.y1 += 2;
        // globalDomain.z0 -= 2;  // z-periodicity
        // globalDomain.z1 += 2;
        // std::vector<MultiBlock3D *> bcargs;
        // bcargs.push_back(lattice);
        // bcargs.push_back(rhoBar);
        // bcargs.push_back(j);
        
        // int bcType = 1;
        // integrateProcessingFunctional(
        //     new VirtualOutlet<T, DESCRIPTOR>(outsideDensity, globalDomain, bcType), param.outlet,
        //     bcargs, 2);
        // setBoundaryVelocity(*lattice, param.outlet, uBoundary);
    }
}

// Write VTK file for the flow around the obstacle, to be viewed with Paraview.
void writeVTK(OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> &bc, plint iT)
{
    VtkImageOutput3D<T> vtkOut(createFileName("volume", iT, PADDING), param.dx);
    vtkOut.writeData<float>(
        *bc.computeVelocityNorm(param.boundingBox()), "velocityNorm", param.dx / param.dt);
    vtkOut.writeData<3, float>(
        *bc.computeVelocity(param.boundingBox()), "velocity", param.dx / param.dt);
    vtkOut.writeData<float>(
        *bc.computePressure(param.boundingBox()), "pressure",
        param.dx * param.dx / (param.dt * param.dt));   //can you read these files in python ?
}

// Write PPM images on slices.
void writePPM(OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> &bc, plint iT)
{
    Box3D xSlice(param.cxLB, param.cxLB, 0, param.ny - 1, 0, param.nz - 1);
    Box3D ySlice(0, param.nx - 1, param.cyLB, param.cyLB, 0, param.nz - 1);
    Box3D zSlice(0, param.nx - 1, 0, param.ny - 1, param.czLB, param.czLB);

    ImageWriter<T> writer("leeloo");
    writer.writeScaledPpm(
        createFileName("vnorm_xslice", iT, PADDING), *bc.computeVelocityNorm(xSlice));
    writer.writeScaledPpm(
        createFileName("vnorm_yslice", iT, PADDING), *bc.computeVelocityNorm(ySlice));
    writer.writeScaledPpm(
        createFileName("vnorm_zslice", iT, PADDING), *bc.computeVelocityNorm(zSlice));
}

void runProgram()
{
    uint32_t seed = 1;

    /*
     * Read the obstacle geometry.
     */

    pcout << std::endl << "Reading STL data for the obstacle geometry." << std::endl;
    Array<T, 3> center(param.cx, param.cy, param.cz);
    Array<T, 3> centerLB(param.cxLB, param.cyLB, param.czLB);
    // The triangle-set defines the surface of the geometry.
    TriangleSet<T> triangleSet(param.geometry_fname, DBL);  //DBL is double precision TriangleSet (std::string fname, Precision precision_=DBL, SurfaceGeometryFileFormat fformat=STL, TriangleSelector< T > *selector=0)

    // Place the obstacle in the correct place in the simulation domain.
    // Here the "geometric center" of the obstacle is computed manually,
    // by computing first its bounding cuboid. In cases that the STL
    // file with the geometry of the obstacle contains its center as
    // the point, say (0, 0, 0), then the following variable
    // "obstacleCenter" must be set to (0, 0, 0) manually.
    Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    Array<T, 3> obstacleCenter = (T)0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner); //gets center of obstacle
    
    triangleSet.translate(-obstacleCenter); //centers it at 0,0,0 i think. indeed, tested below. 
    // Cuboid<T> b1Cuboid = triangleSet.getBoundingCuboid();
    // Array<T, 3> ne1wobstacleCenter = (T)0.5 * (b1Cuboid.lowerLeftCorner + b1Cuboid.upperRightCorner);
    // pcout << "obstacle center x goes from: " << obstacleCenter[0] << "to : " << ne1wobstacleCenter[0] << std::endl;
    // pcout << "obstacle center y goes from: " << obstacleCenter[1] << "to : " << ne1wobstacleCenter[1] << std::endl;
    // pcout << "obstacle center z goes from: " << obstacleCenter[2] << "to : " << ne1wobstacleCenter[2] << std::endl;
    triangleSet.scale(1.0 / param.dx);  // In lattice units from now on...
    triangleSet.translate(centerLB);
    TriangleSet<T> triangleSet2;
    TriangleSet<T> triangleSet3;
    
    Plane<T> planeyminus(Array<T, 3>(0., 0., 0.),Array<T, 3>(0., -1., 0.) );
    Plane<T> planeyplus(Array<T, 3>(0., param.ny , 0.),Array<T, 3>(0., 1., 0.) );
    triangleSet.cutWithPlane(planeyminus, triangleSet2);
    triangleSet2.cutWithPlane(planeyplus,triangleSet3);
    // triangleSet.rotate((T) -3.1415/2.0,(T) -7.8 * 3.1415 / 180, (T) 3.1415/2 );
    triangleSet3.writeBinarySTL(outputDir + "obstacle_LB.stl");  //palabos is crazy...writeBinarySTL and done!

    // The DEFscaledMesh, and the triangle-boundary are more sophisticated data
    // structures used internally by Palabos to treat the boundary.
    plint xDirection = 0;  //reference direction idk exactly what it is. 
    plint borderWidth = 1;  // Because Guo acts in a one-cell layer.
                            // Requirement: margin>=borderWidth.
    plint margin =
        1;  // Extra margin of allocated cells around the obstacle, for the case of moving walls.
    plint blockSize = 0;  // Size of blocks in the sparse/parallel representation.
                          // Zero means: don't use sparse representation. A dense grid stores information for every node, while a sparse grid only stores information for non-empty or significant nodes. Sparse representations can save memory and improve performance, especially for large simulations with many empty or uniform regions.
    DEFscaledMesh<T> defMesh(triangleSet3, 0, xDirection, margin, Dot3D(0, 0, 0));  //DEFscaledMesh (TriangleSet< T > const &triangleSet_, plint resolution_, plint referenceDirection_, plint margin_, Dot3D location)
    TriangleBoundary3D<T> boundary(defMesh); //TriangleBoundary3D (DEFscaledMesh< T > const &defMesh, bool automaticCloseHoles=true)
    // boundary.getMesh().inflate();
    //they basically redefine the mesh and boundary to have more sophisticated data 

    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "Number of iterations in an integral time scale: " << (plint)(1.0 / param.dt) << std::endl;

    /*
     * Voxelize the domain.
     */

    // Voxelize the domain means: decide which lattice nodes are inside the obstacle and which are
    // outside. Interesting!
    pcout << std::endl << "Voxelizing the domain." << std::endl;
    plint extendedEnvelopeWidth = 2;  // Extrapolated off-lattice BCs.
    const int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain(
        boundary, flowType, param.boundingBox(), borderWidth, extendedEnvelopeWidth, blockSize); //VoxelizedDomain3D (TriangleBoundary3D< T > const &boundary_, int flowType_, Box3D const &boundingBox, plint borderWidth_, plint envelopeWidth_, plint blockSize_, plint gridLevel_=0, bool dynamicMesh_=false)
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    MultiScalarField3D<int> flagMatrix((MultiBlock3D &)voxelizedDomain.getVoxelMatrix());
    setToConstant(
        flagMatrix, voxelizedDomain.getVoxelMatrix(), voxelFlag::inside,
        flagMatrix.getBoundingBox(), 1);
    // setToConstant(
    //     flagMatrix, voxelizedDomain.getVoxelMatrix(), voxelFlag::innerBorder,
    //     flagMatrix.getBoundingBox(), 1);
    pcout << "Number of fluid cells: " << computeSum(flagMatrix) << std::endl;
    // pcout << "VoxelMatrix" << flagMatrix << std::endl;
    // pcout << voxelizedDomain.getVoxelMatrix()
    //prints this:
    // Size of the multi-block:     260-by-79-by-79
    // Number of atomic-blocks:     1
    // Smallest atomic-block:       260-by-79-by-79
    // Largest atomic-block:        260-by-79-by-79
    // Number of allocated cells:   1.62266 million
    // Fraction of allocated domain: 100 percent

    /*
     * Generate the lattice, the density and momentum blocks.
     */

    pcout << "Generating the lattice, the rhoBar and j fields." << std::endl;
    MultiBlockLattice3D<T, DESCRIPTOR> *lattice =
        new MultiBlockLattice3D<T, DESCRIPTOR>(voxelizedDomain.getVoxelMatrix());  //generates our beautiful lattice, with our voxelized domain and the descriptor.
    if (param.useSmago) {  //for now i won't look into the LES thing. 
        defineDynamics(
            *lattice, lattice->getBoundingBox(),
            new SmagorinskyBGKdynamics<T, DESCRIPTOR>(param.omega, param.cSmago));
        pcout << "Using Smagorinsky BGK dynamics." << std::endl;
    } else {
        defineDynamics(
            *lattice, lattice->getBoundingBox(), new BGKdynamics<T, DESCRIPTOR>(param.omega));  //defines the collision dynamics defineDynamics(lattice, domain, dynamics*)
        pcout << "Using BGK dynamics." << std::endl;
    }
    bool velIsJ = false;
    defineDynamics(
        *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
        new NoDynamics<T, DESCRIPTOR>(), voxelFlag::inside); //here it defines the non-existing dynamics of inside the obstacle. 
    lattice->toggleInternalStatistics(false); //This is due to the fact that a reduction acts like a synchronization barrier, a fact which is observed to have a negative
// impact on performance, especially on cluster-like parallel machines. It is therefore possible to turn off the internal
// statistics (and thus, the reduction operations), through a function call

    // The rhoBar and j fields are used at both the collision and at the implementation of the
    // outflow boundary condition.
    plint envelopeWidth = 1;
    MultiScalarField3D<T> *rhoBar =
        generateMultiScalarField<T>((MultiBlock3D &)*lattice, envelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);

    MultiTensorField3D<T, 3> *j =
        generateMultiTensorField<T, 3>((MultiBlock3D &)*lattice, envelopeWidth).release();
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock3D *> lattice_rho_bar_j_arg;
    lattice_rho_bar_j_arg.push_back(lattice);
    lattice_rho_bar_j_arg.push_back(rhoBar);
    lattice_rho_bar_j_arg.push_back(j);
    integrateProcessingFunctional(
        new ExternalRhoJcollideAndStream3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
        lattice_rho_bar_j_arg, 0);
    integrateProcessingFunctional(
        new BoxRhoBarJfunctional3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
        lattice_rho_bar_j_arg, 3);  // rhoBar and j are computed at level 3 because
                                    // the boundary conditions are on levels 0, 1 and 2.
                                    //There are different ways to control the order in which internal data processors are executed in the function call
// executeInternalProcessors(). First of all, each data processor is attributed to a processor level, and these
// processor levels are traversed in increasing order, starting with level 0. By default, all internal processors are attributed
// to level 0, but you have the possibility to put them into any other level, specified as the last, optional parameter of the
// function addInternalProcessor or integrateProcessingFunctional. Inside a processor level, the
// data processors are executed in the order in which they were added to the block. Additionally to imposing an order
// of execution, the attribution of data processors to a given level has an influence on the communication pattern inside
// multi-blocks. As a matter of fact, communication is not immediately performed after the execution of a data processor
// with write access, but only when switching from one level to the next.

    /*
     * Generate the off-lattice boundary condition on the obstacle and the outer-domain boundary
     * conditions.
     */

    pcout << "Generating boundary conditions." << std::endl;

    // OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> *boundaryCondition;  //defined above, Velocity is a size 3 array type definition.
    
    BoundaryProfiles3D<T, Velocity> profiles;
    bool useAllDirections = true;
    // OffLatticeModel3D<T, Velocity> *offLatticeModel = 0;
    if (param.freeSlipWall) {
        profiles.setWallProfile(new FreeSlipProfile3D<T>);
    } else {
        profiles.setWallProfile(new NoSlipProfile3D<T>);
    }
    BouzidiOffLatticeModel3D<T, DESCRIPTOR> *model =
        new BouzidiOffLatticeModel3D<T, DESCRIPTOR>(
            new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles),
            flowType);

    OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> *boundaryCondition =
        new OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity>(
            model->clone(), voxelizedDomain, *lattice);
    // offLatticeModel = new GuoOffLatticeModel3D<T, DESCRIPTOR>( //This class implements the Guo (GZS,2002) boundary condition on a BoundaryShape
    //     new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles), flowType, useAllDirections);
    // offLatticeModel->setVelIsJ(velIsJ);
    // boundaryCondition = new OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity>(
    //     offLatticeModel, voxelizedDomain, *lattice);

    boundaryCondition->insert();

    //this is the most difficult part of the code so far. Help!

    // The boundary condition algorithm or the outer domain.
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *outerBoundaryCondition =
        createLocalBoundaryCondition3D<T, DESCRIPTOR>();
    outerDomainBoundaries(lattice, rhoBar, j, outerBoundaryCondition);

    /*
     * Implement the outlet sponge zone.
     */

    if (param.numOutletSpongeCells > 0) {
        T bulkValue;
        Array<plint, 6> numSpongeCells;

        if (param.outletSpongeZoneType == 0) {
            pcout << "Generating an outlet viscosity sponge zone." << std::endl;
            bulkValue = param.omega;
        } else if (param.outletSpongeZoneType == 1) {
            pcout << "Generating an outlet Smagorinsky sponge zone." << std::endl;
            bulkValue = param.cSmago;
        } else {
            pcout << "Error: unknown type of sponge zone." << std::endl;
            exit(-1);
        }

        // Number of sponge zone lattice nodes at all the outer domain boundaries.
        // So: 0 means the boundary at x = 0
        //     1 means the boundary at x = nx-1
        //     2 means the boundary at y = 0
        //     and so on...
        numSpongeCells[0] = 0;
        numSpongeCells[1] = param.numOutletSpongeCells;
        numSpongeCells[2] = 0;
        numSpongeCells[3] = 0;
        numSpongeCells[4] = 0;
        numSpongeCells[5] = 0;

        std::vector<MultiBlock3D *> args;
        args.push_back(lattice);

        if (param.outletSpongeZoneType == 0) {
            applyProcessingFunctional(
                new ViscositySpongeZone3D<T, DESCRIPTOR>(
                    param.nx, param.ny, param.nz, bulkValue, numSpongeCells),
                lattice->getBoundingBox(), args);
        } else {
            applyProcessingFunctional(
                new SmagorinskySpongeZone3D<T, DESCRIPTOR>(
                    param.nx, param.ny, param.nz, bulkValue, param.targetSpongeCSmago,
                    numSpongeCells),
                lattice->getBoundingBox(), args);
        }
    }

    /*
     * Setting the initial conditions.
     */

    // Initial condition: Constant pressure and velocity-at-infinity everywhere.
    Array<T, 3> uBoundary(param.getInletVelocity(0)*std::cos(angle), (T)0.0, param.getInletVelocity(0)*std::sin(angle));
    initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), (T)1.0, uBoundary);
    applyProcessingFunctional(
        new BoxRhoBarJfunctional3D<T, DESCRIPTOR>(), lattice->getBoundingBox(),
        lattice_rho_bar_j_arg);  // Compute rhoBar and j before VirtualOutlet is executed.
    // lattice->executeInternalProcessors(1); // Execute all processors except the ones at level 0.
    // lattice->executeInternalProcessors(2);
    // lattice->executeInternalProcessors(3);

    /*
     * Particles (streamlines).
     */

    // This part of the code that relates to particles, is purely for visualization
    // purposes. Particles are used to compute streamlines essentially.

    // Definition of a particle field.
    MultiParticleField3D<ParticleFieldT> *particles = 0;

    if (param.useParticles) {
        particles = new MultiParticleField3D<ParticleFieldT>(
            lattice->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getCombinedStatistics());

        std::vector<MultiBlock3D *> particleArg;
        particleArg.push_back(particles);

        std::vector<MultiBlock3D *> particleFluidArg;
        particleFluidArg.push_back(particles);
        particleFluidArg.push_back(lattice);

        // Functional that advances the particles to their new position at each predefined time
        // step.
        integrateProcessingFunctional(
            new AdvanceParticlesEveryWhereFunctional3D<T, DESCRIPTOR>(param.cutOffSpeedSqr),
            lattice->getBoundingBox(), particleArg, 0);
        // Functional that assigns the particle velocity according to the particle's position in the
        // fluid.
        integrateProcessingFunctional(
            new FluidToParticleCoupling3D<T, DESCRIPTOR>((T)param.particleTimeFactor),
            lattice->getBoundingBox(), particleFluidArg, 1);

        // Definition of a domain from which particles will be injected in the flow field.
        Box3D injectionDomain(
            0, 0, centerLB[1] - 0.25 * param.ny, centerLB[1] + 0.25 * param.ny,
            centerLB[2] - 0.25 * param.nz, centerLB[2] + 0.25 * param.nz);

        // Definition of simple mass-less particles.
        Particle3D<T, DESCRIPTOR> *particleTemplate = 0;
        particleTemplate =
            new PointParticle3D<T, DESCRIPTOR>(0, Array<T, 3>(0., 0., 0.), Array<T, 3>(0., 0., 0.));

        // Functional which injects particles with predefined probability from the specified
        // injection domain.
        std::vector<MultiBlock3D *> particleInjectionArg;
        particleInjectionArg.push_back(particles);

        integrateProcessingFunctional(
            new InjectRandomParticlesFunctionalPPRNG3D<T, DESCRIPTOR>(
                particleTemplate, param.particleProbabilityPerCell, particles->getBoundingBox(),
                &seed),
            injectionDomain, particleInjectionArg, 0);

        // Definition of an absorbtion domain for the particles.
        Box3D absorbtionDomain(param.outlet);

        // Functional which absorbs the particles which reach the specified absorbtion domain.
        integrateProcessingFunctional(
            new AbsorbParticlesFunctional3D<T, DESCRIPTOR>, absorbtionDomain, particleArg, 0);

        particles->executeInternalProcessors();
    }

    /*
     * Starting the simulation.
     */

    plb_ofstream energyFile((outputDir + "average_energy.dat").c_str());

    pcout << std::endl;
    pcout << "Starting simulation." << std::endl;
    bool checkForErrors = true;
    for (plint i = 0; i < param.maxIter; ++i) {
        if (i <= param.initialIter) {
            Array<T, 3> uBoundary(param.getInletVelocity(i)*std::cos(angle), 0.0, param.getInletVelocity(i)*std::sin(angle));
            setBoundaryVelocity(*lattice, param.inlet, uBoundary);
            setBoundaryVelocity(*lattice, param.outlet, uBoundary);
            setBoundaryVelocity(*lattice, param.lateral3, uBoundary);
            setBoundaryVelocity(*lattice, param.lateral4, uBoundary);
            
            // pcout << "The speed in x is " << uBoundary[0] << " and in z it is " << uBoundary[1] << std::endl;
        }

        if (i % param.statIter == 0) {
            pcout << "At iteration " << i << ", t = " << i * param.dt << std::endl;
            if (i != 0) {
                Array<T, 3> force(boundaryCondition->getForceOnObject());
                T factor = util::sqr(util::sqr(param.dx)) / util::sqr(param.dt);
                pcout << "Force on object over fluid density: F[x] = " << force[0] * factor
                      << ", F[y] = " << force[1] * factor << ", F[z] = " << force[2] * factor
                      << std::endl;
                
            }
            T avEnergy = boundaryCondition->computeAverageEnergy() * util::sqr(param.dx)
                         / util::sqr(param.dt);
            pcout << "Average kinetic energy over fluid density: E = " << avEnergy << std::endl;
            energyFile << i * param.dt << "  " << avEnergy << std::endl;
            pcout << std::endl;
        }

        if (i % param.vtkIter == 0) {
            pcout << "Writing VTK at time t = " << i * param.dt << endl;
            writeVTK(*boundaryCondition, i);
            if (param.useParticles) {
                writeParticleVtk<T, DESCRIPTOR>(
                    *particles, createFileName(outputDir + "particles_", i, PADDING) + ".vtk",
                    param.dx, param.maxNumParticlesToWrite);
            }
        }

        if (i % param.imageIter == 0) {
            pcout << "Writing PPM image at time t = " << i * param.dt << endl;
            writePPM(*boundaryCondition, i);
        }

        lattice->executeInternalProcessors();
        lattice->incrementTime();
        if (param.useParticles && i % param.particleTimeFactor == 0) {
            particles->executeInternalProcessors();
        }
        seed++;

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }
    }

    energyFile.close();
    delete outerBoundaryCondition;
    delete boundaryCondition;
    if (param.useParticles) {
        delete particles;
    }
    delete j;
    delete rhoBar;
    delete lattice;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outputDir);

    // The try-catch blocks catch exceptions in case an error occurs,
    // and terminate the program properly with a nice error message.

    // 1. Read command-line parameter: the input file name.
    string xmlFileName;
    try {
        global::argv(1).read(xmlFileName);
    } catch (PlbIOException &exception) {
        pcout << "Wrong parameters; the syntax is: " << (std::string)global::argv(0)
              << " input-file.xml" << std::endl;
        return -1;
    }

    // 2. Read input parameters from the XML file.
    try {
        param = Param(xmlFileName);
    } catch (PlbIOException &exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }

    // 3. Execute the main program.
    try {
        runProgram();
    } catch (PlbIOException &exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
