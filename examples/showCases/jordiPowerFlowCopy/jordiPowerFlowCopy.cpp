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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;
typedef Array<T, 3> Velocity;

#define RESCALER   ConvectiveNoForceRescaler    
#define DESCRIPTOR descriptors::D3Q19Descriptor
T angle = 7.8 * 3.1415 / 180; //angle of attack. Note that the airfoil's chord is aligned with the domain, and it is the air that comes at an angle. 
T chordLengthPercentage = 0.2/(7.790592*2); //ratio between domain length and chord length. Useful for calculation below. 

bool teaddon = false; //later on this will be used. Basically, if the STL model has a Trailing Edge Add-on, it will still center the STL file based on the airfoil center, without considering the add-on. It is set as false bc for now we are using NACA0018
struct SimulationParameters {
    /*
     * Parameters set by the user.
     */

    // Geometry.

    std::string staticSurfaceFileName;  // Files with the static immersed surface geometries.

    // Grid.
    std::string gridDensityFunctionFile;  // File with the discrete grid density function.
    plint minLeafLevel;                   // Octree leaf levels.
    plint maxLeafLevel;
    plint nBlock;  // Block size.

    // Boundary Conditions.
    T inletVelocity;  // Inlet velocity x-component in physical units.

    int outflowBcType;  // Type of the outflow boundary condition.
                        // If 0, then a constant velocity is imposed at the outlet.
                        // If 1, then a constant pressure is imposed at the outlet.
                        // If 2, then a velocity Neumann condition is imposed at the outlet.

    Array<T, 6> spongeWidths;  // Sponge zone widths (0.0 for no sponge zone for the corresponding
                               // lattice boundary).
    // Numerics.
    Precision precision;  // Precision for geometric operations.
    T L_Ref;              // Length to define the Reynolds number.
    T u_Ref;              // Characteristic velocity.
    T u_LB;               // Lattice velocity (w.r.t the characteristic velocity).
    plint maxIter;        // Maximum number of iterations at the coarsest level.

    // Fluid.

    T rho;              // Fluid density in physical units
    T Re;               // Reynolds number.
    T nu;               // Fluid kinematic viscosity in physical units.
    T ambientPressure;  // Absolute stagnation pressure in physical units.

    // Output.

    std::string outDir;    // Output directory.
    plint statIter;        // Number of iterations for terminal output at the coarsest level.
    plint outIter;         // Number of iterations for disk output at the coarsest level.
    bool computeAverages;  // Compute average and RMS values or not?
    plint avgIter;         // Number of iterations to start averaging at the coarsest level.
    plint minOutputLevel;  // Minimum grid refinement level for output on disk.
    plint maxOutputLevel;  // Maximum grid refinement level for output on disk.

    bool outputInDomain;     // Save data on disk in a volume domain or not?
    Cuboid<T> outputCuboid;  // Volume domain for disk output.

    bool outputOnSlices;        // Save data on disk on a set of slices or not?
    std::vector<T> xPositions;  // Positions of the x-slices for output.
    std::vector<T> xyRange;     // y range of the x-slices.
    std::vector<T> xzRange;     // z range of the x-slices.
    std::vector<T> yPositions;  // Positions of the y-slices for output.
    std::vector<T> yzRange;     // z range of the y-slices.
    std::vector<T> yxRange;     // x range of the y-slices.
    std::vector<T> zPositions;  // Positions of the z-slices for output.
    std::vector<T> zxRange;     // x range of the z-slices.
    std::vector<T> zyRange;     // y range of the z-slices.

    plint cpIter;  // Number of iterations for checkpointing.
    plint abIter;  // Number of iterations for checking for user-driven program abortion.
    std::string abortFileName;        // File for signaling program abortion (inside the outDir).
    std::string xmlContinueFileName;  // XML file for restarting (inside the outDir).
    std::string baseFileName;         // Basename of the checkpoint files (inside the outDir).
    bool useParallelIO;  // For a desktop PC this should be "false", for a cluster "true".

    /*
     * Parameters NOT set by the user.
     */

    std::string surfaceName;

    OctreeGridStructure ogs;
    Cuboid<T> fullDomain;
    plint finestLevel;
    T dxCoarsest, dtCoarsest;
    T dxFinest, dtFinest;

    T rho_LB;
    Array<T, 3> inletVelocity_LB;
    std::vector<T> omega;
    Array<T, 3> physicalLocation;
    plint smallEnvelopeWidth;
    plint mediumEnvelopeWidth;
    plint largeEnvelopeWidth;
    bool incompressibleModel;
    plint fileNamePadding;
    bool saveDynamicContent;

    std::map<plint, std::vector<Box3D> > outputDomains;  // Output domains per output level.
    std::vector<std::string> outputDomainNames;
};

T toLB(T physVal, plint direction, T dx, Array<T, 3> const &location)  //this function normalizes distances to have them in lattice units, e.g. 10 voxels of distance.
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return ((physVal - location[direction]) / dx);
}

Array<T, 3> toLB(Array<T, 3> const &physVal, T dx, Array<T, 3> const &location)
{
    return ((physVal - location) / dx);
}

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters &param)
{
    XMLreader document(xmlInputFileName);

    std::vector<T> fullDomainX, fullDomainY, fullDomainZ;
    document["geometry"]["simulationDomain"]["x"].read(fullDomainX);
    PLB_ASSERT(fullDomainX.size() == 2 && fullDomainX[1] > fullDomainX[0]);
    document["geometry"]["simulationDomain"]["y"].read(fullDomainY);
    PLB_ASSERT(fullDomainY.size() == 2 && fullDomainY[1] > fullDomainY[0]);
    document["geometry"]["simulationDomain"]["z"].read(fullDomainZ);
    PLB_ASSERT(fullDomainZ.size() == 2 && fullDomainZ[1] > fullDomainZ[0]);
    param.fullDomain.lowerLeftCorner[0] = fullDomainX[0];
    param.fullDomain.lowerLeftCorner[1] = fullDomainY[0];
    param.fullDomain.lowerLeftCorner[2] = fullDomainZ[0];
    param.fullDomain.upperRightCorner[0] = fullDomainX[1];
    param.fullDomain.upperRightCorner[1] = fullDomainY[1];
    param.fullDomain.upperRightCorner[2] = fullDomainZ[1];

    document["geometry"]["staticSurfaceFileName"].read(param.staticSurfaceFileName);

    document["grid"]["gridDensityFunctionFile"].read(param.gridDensityFunctionFile);
    abortIfCannotOpenFileForReading(param.gridDensityFunctionFile);
    plint numLevels = 0;
    document["grid"]["numLevels"].read(numLevels); 
    PLB_ASSERT(numLevels >= 1);
    plint maxOctreeLevel = 0;
    document["grid"]["maxOctreeLevel"].read(maxOctreeLevel); 
    if (maxOctreeLevel - numLevels + 1 < 0) {
        maxOctreeLevel = numLevels - 1;
    }
    param.maxLeafLevel = maxOctreeLevel;
    param.minLeafLevel = maxOctreeLevel - numLevels + 1;
    pcout << "maxLeafLevel is " << param.maxLeafLevel << " and min leaf level is " << param.minLeafLevel << std::endl;
    document["grid"]["nBlock"].read(param.nBlock);
    PLB_ASSERT(param.nBlock >= 6);

    document["boundaryConditions"]["inletVelocity"].read(param.inletVelocity);

    document["boundaryConditions"]["outflowBcType"].read(param.outflowBcType);
    PLB_ASSERT(param.outflowBcType >= 0 && param.outflowBcType <= 2);

    std::vector<T> zoneWidths;
    document["boundaryConditions"]["spongeZones"]["xWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[0] = zoneWidths[0];
    param.spongeWidths[1] = zoneWidths[1];
    zoneWidths.clear();
    document["boundaryConditions"]["spongeZones"]["yWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[2] = zoneWidths[0];
    param.spongeWidths[3] = zoneWidths[1];
    zoneWidths.clear();
    document["boundaryConditions"]["spongeZones"]["zWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[4] = zoneWidths[0];
    param.spongeWidths[5] = zoneWidths[1];
    zoneWidths.clear();

    std::string precision;
    document["numerics"]["precision"].read(precision);
    PLB_ASSERT(
        precision == "FLT" || precision == "DBL" || precision == "LDBL" || precision == "INF");
    if (precision == "FLT") {
        param.precision = FLT;
    } else if (precision == "DBL") {
        param.precision = DBL;
    } else if (precision == "LDBL") {
        param.precision = LDBL;
    } else {
        param.precision = INF;
    }
    document["numerics"]["characteristicLength"].read(param.L_Ref);
    PLB_ASSERT(util::greaterThan_abs(param.L_Ref, (T)0));
    document["numerics"]["characteristicVelocity"].read(param.u_Ref);
    PLB_ASSERT(util::greaterThan_abs(param.u_Ref, (T)0));
    document["numerics"]["uLB"].read(param.u_LB);
    document["numerics"]["maxIter"].read(param.maxIter);
    PLB_ASSERT(param.maxIter > 0);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["Re"].read(param.Re);
    param.nu = param.u_Ref * param.L_Ref / param.Re;
    document["fluid"]["ambientPressure"].read(param.ambientPressure);

    std::string outDir;
    document["output"]["outDir"].read(outDir);
    if (outDir[outDir.size() - 1] != '/') {
        outDir += '/';
    }
    param.outDir = outDir;
    abortIfCannotCreateFileInDir(param.outDir, "plb-checkfile.txt");

    document["output"]["statIter"].read(param.statIter);
    PLB_ASSERT(param.statIter > 0);
    document["output"]["outIter"].read(param.outIter);
    PLB_ASSERT(param.outIter > 0);
    document["output"]["computeAverages"].read(param.computeAverages);

    if (param.computeAverages) {
        document["output"]["avgIter"].read(param.avgIter);
        PLB_ASSERT(param.avgIter >= 0);
    } else {
        param.avgIter = 0;
    }

    document["output"]["minOutputLevel"].read(param.minOutputLevel);
    document["output"]["maxOutputLevel"].read(param.maxOutputLevel);
    PLB_ASSERT(param.minOutputLevel <= param.maxOutputLevel);

    document["output"]["outputInDomain"].read(param.outputInDomain);
    if (param.outputInDomain) {
        std::vector<T> x, y, z;
        document["output"]["outputDomain"]["x"].read(x);
        PLB_ASSERT(x.size() == 2 && x[1] > x[0]);
        document["output"]["outputDomain"]["y"].read(y);
        PLB_ASSERT(y.size() == 2 && y[1] > y[0]);
        document["output"]["outputDomain"]["z"].read(z);
        PLB_ASSERT(z.size() == 2 && z[1] > z[0]);
        param.outputCuboid.lowerLeftCorner[0] = x[0];
        param.outputCuboid.lowerLeftCorner[1] = y[0];
        param.outputCuboid.lowerLeftCorner[2] = z[0];
        param.outputCuboid.upperRightCorner[0] = x[1];
        param.outputCuboid.upperRightCorner[1] = y[1];
        param.outputCuboid.upperRightCorner[2] = z[1];
    }

    document["output"]["outputOnSlices"].read(param.outputOnSlices);
    if (param.outputOnSlices) {
        document["output"]["outputSlices"]["xSlices"]["xPositions"].read(param.xPositions);
        document["output"]["outputSlices"]["xSlices"]["yRange"].read(param.xyRange);
        PLB_ASSERT(param.xyRange.size() == 2 && param.xyRange[1] > param.xyRange[0]);
        document["output"]["outputSlices"]["xSlices"]["zRange"].read(param.xzRange);
        PLB_ASSERT(param.xzRange.size() == 2 && param.xzRange[1] > param.xzRange[0]);

        document["output"]["outputSlices"]["ySlices"]["yPositions"].read(param.yPositions);
        document["output"]["outputSlices"]["ySlices"]["zRange"].read(param.yzRange);
        PLB_ASSERT(param.yzRange.size() == 2 && param.yzRange[1] > param.yzRange[0]);
        document["output"]["outputSlices"]["ySlices"]["xRange"].read(param.yxRange);
        PLB_ASSERT(param.yxRange.size() == 2 && param.yxRange[1] > param.yxRange[0]);

        document["output"]["outputSlices"]["zSlices"]["zPositions"].read(param.zPositions);
        document["output"]["outputSlices"]["zSlices"]["xRange"].read(param.zxRange);
        PLB_ASSERT(param.zxRange.size() == 2 && param.zxRange[1] > param.zxRange[0]);
        document["output"]["outputSlices"]["zSlices"]["yRange"].read(param.zyRange);
        PLB_ASSERT(param.zyRange.size() == 2 && param.zyRange[1] > param.zyRange[0]);
    }

    document["output"]["cpIter"].read(param.cpIter);
    document["output"]["abIter"].read(param.abIter);
    PLB_ASSERT(param.abIter > 0);
    document["output"]["abortFileName"].read(param.abortFileName);
    document["output"]["xmlContinueFileName"].read(param.xmlContinueFileName);
    document["output"]["baseFileName"].read(param.baseFileName);
    document["output"]["useParallelIO"].read(param.useParallelIO);

    param.abortFileName = param.outDir + param.abortFileName;
    param.xmlContinueFileName = param.outDir + param.xmlContinueFileName;
    param.baseFileName = param.outDir + param.baseFileName;
}

void createOctreeGridStructure(SimulationParameters &param)
{
    // Default parameters for the octree grid generation.
    // Octree construction parameters.
    bool useSamples = false;
    plint numSamples = -1;
    plint maxIter = 100;
    bool removeBlocks = true;
    bool fineToCoarse = true;
    // Assign processIds from the finest to the coarsest level or the other way around.
     // We choose that both overlaps and allocatable leaves are grouped in the same
     // direction (fine -> coarse, or coarse -> fine).
    int numLevelsToGroupBlocks = 0;
    int numLevelsToGroupOverlaps = -1;// The octree nodes can be grouped or not. By grouping nodes we mean
     // that we assign adjacent nodes to the same process (they have
     // the same processId).
     // If numLevelsToGroup* = 0, this means that no grouping takes
     // place.
     // If numLevelsToGroup* = 1, this means that the adjacent
     // nodes are grouped "per level", so in groups of 8 nodes at
     // most.
     // If numLevelsToGroup* = 2, this means that the groups have
     // 64 nodes at most.
     // In general, if numLevelsToGroup* = l, then the respective
     // groups have at most 8^l nodes.
    //  int numLevelsToGroupBlocks;    // For allocatable leaves.
    //  int numLevelsToGroupOverlaps;  // For overlaps.
    bool strongGrouping = false;// If strongGrouping is true, then the overlap blocks are grouped with
     // all their leaf children.
     // If strongGrouping is false, then the overlap blocks are grouped
     // only with their leaf children to which they are coupled.
    bool verbose = true;// Output parameters.
    bool stlOutput = true;
    std::string stlBaseName = "octree";
    bool xPeriodic = false;
    bool yPeriodic = false;
    bool zPeriodic = false;
    T gridDensityScaleFactor = (T)1;

    OctreeGridGenerator<T> octreeGridGenerator(
        param.fullDomain, param.gridDensityFunctionFile, param.minLeafLevel, param.maxLeafLevel,
        param.nBlock, global::mpi().getSize(), xPeriodic, yPeriodic, zPeriodic,
        gridDensityScaleFactor, useSamples, numSamples, maxIter, removeBlocks, fineToCoarse,
        numLevelsToGroupBlocks, numLevelsToGroupOverlaps, strongGrouping, param.outDir, verbose,
        stlOutput, stlBaseName);
    //Hier above is a really smart class that generates an ocree grid structure based on the density file that you input, and based on the maximum and minimum levels that you want. Check Notion user guide for more info.
    param.ogs = octreeGridGenerator.generateOctreeGridStructure();
    param.fullDomain = octreeGridGenerator.getFullDomain();
    param.dxFinest = octreeGridGenerator.getDxFinestLevel();

    param.ogs.writeXML(octreeGridGenerator.getOutDir() + "octreeGridStructure.xml");

    {
        XMLwriter gridInfoXML;
        Array<T, 6> fullDomain(
            octreeGridGenerator.getFullDomain().x0(), octreeGridGenerator.getFullDomain().x1(),
            octreeGridGenerator.getFullDomain().y0(), octreeGridGenerator.getFullDomain().y1(),
            octreeGridGenerator.getFullDomain().z0(), octreeGridGenerator.getFullDomain().z1());
        gridInfoXML["fullDomain"].set<T, 6>(fullDomain);
        gridInfoXML["dxFinestLevel"].set(octreeGridGenerator.getDxFinestLevel());
        gridInfoXML.print(octreeGridGenerator.getOutDir() + "octreeGridInfo.xml");
    }
}

Box3D computeOutputDomain( 
    SimulationParameters const &param, Cuboid<T> const &outputCuboid, plint level)
{
    T dx = param.dxFinest * (T)util::intTwoToThePower(param.finestLevel - level); //bc levels of resolution scale by factors of 2
    Box3D outputBox;
    outputBox.x0 = (plint)toLB(outputCuboid.x0(), 0, dx, param.physicalLocation); //this toLB function does (physlocation-x0)/dx
    outputBox.x1 = (plint)toLB(outputCuboid.x1(), 0, dx, param.physicalLocation) + (plint)1;
    outputBox.y0 = (plint)toLB(outputCuboid.y0(), 1, dx, param.physicalLocation);
    outputBox.y1 = (plint)toLB(outputCuboid.y1(), 1, dx, param.physicalLocation) + (plint)1;
    outputBox.z0 = (plint)toLB(outputCuboid.z0(), 2, dx, param.physicalLocation);
    outputBox.z1 = (plint)toLB(outputCuboid.z1(), 2, dx, param.physicalLocation) + (plint)1;
    return (outputBox);
}

void computeAllOutputDomains(SimulationParameters &param)
{
    std::vector<Cuboid<T> > outputCuboids; //list of all output cuboids. 
    if (param.outputInDomain) {
        outputCuboids.push_back(param.outputCuboid); //outputs the output cuboid that the user set. 
        param.outputDomainNames.push_back("domain");
    }

    if (param.outputOnSlices) { //if output onslices boolean is true
        plint numXdigits = util::val2str(param.xPositions.size()).length(); //number of digits.
        for (plint i = 0; i < (plint)param.xPositions.size(); i++) {
            Array<T, 3> llc(param.xPositions[i], param.xyRange[0], param.xzRange[0]); //first point of the slice (llc means lowerleftcorner)
            Array<T, 3> urc(param.xPositions[i], param.xyRange[1], param.xzRange[1]); //second point of the slice (upperrightcorner)
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_x_", i, numXdigits + 1));
        }

        plint numYdigits = util::val2str(param.yPositions.size()).length();
        for (plint i = 0; i < (plint)param.yPositions.size(); i++) {
            Array<T, 3> llc(param.yxRange[0], param.yPositions[i], param.yzRange[0]);
            Array<T, 3> urc(param.yxRange[1], param.yPositions[i], param.yzRange[1]);
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_y_", i, numYdigits + 1));
        }

        plint numZdigits = util::val2str(param.zPositions.size()).length();
        for (plint i = 0; i < (plint)param.zPositions.size(); i++) {
            Array<T, 3> llc(param.zxRange[0], param.zyRange[0], param.zPositions[i]);
            Array<T, 3> urc(param.zxRange[1], param.zyRange[1], param.zPositions[i]);
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_z_", i, numZdigits + 1));
        }
    }

    std::vector<plint> domainExistsInNumLevels(outputCuboids.size(), 0); //vector of 0s with size output cuboids
    for (plint iLevel = param.minOutputLevel; iLevel <= param.maxOutputLevel; iLevel++) { //loop over output levels
        std::vector<Box3D> outputDomainsAtLevel;
        for (plint iCuboid = 0; iCuboid < (plint)outputCuboids.size(); iCuboid++) {
            Box3D outputDomain =
                computeOutputDomain(param, outputCuboids[iCuboid], param.minOutputLevel);
            if (iLevel != param.minOutputLevel) {
                outputDomain =
                    outputDomain.multiply(util::intTwoToThePower(iLevel - param.minOutputLevel));
                // outputDomain.x1--;
                // outputDomain.y1--;
                // outputDomain.z1--;
            }
            if (outputDomain.getNx() <= 0 || outputDomain.getNy() <= 0 || outputDomain.getNz() <= 0) //if domain doesnt exists officially break it
            {
                outputDomain = Box3D(-1, -1, -1, -1, -1, -1);
            } else {
                domainExistsInNumLevels[iCuboid] += 1; //if its not broken set 1 to the vector of 0s
            }
            outputDomainsAtLevel.push_back(outputDomain); //add the output domain at level
        }
        PLB_ASSERT(param.outputDomainNames.size() == outputDomainsAtLevel.size());
        param.outputDomains[iLevel] = outputDomainsAtLevel; //all ouptut domains at level are put to ouputdomains vector
    }
    for (plint iCuboid = 0; iCuboid < (plint)outputCuboids.size(); iCuboid++) {
        PLB_ASSERT(domainExistsInNumLevels[iCuboid] != 0); //chekc that all domain exists in num levels.
    }
}

void calculateDerivedSimulationParameters(SimulationParameters &param)
{
    // Derived quantities.

    param.smallEnvelopeWidth = 1;
    param.mediumEnvelopeWidth = 2;
    param.largeEnvelopeWidth = 4;
    param.fileNamePadding = 8;
    param.saveDynamicContent = true;

    plint numLevels = param.ogs.getNumLevels(); //get number of levels in octree grid structure. 
    PLB_ASSERT(numLevels >= 1);
    param.finestLevel = numLevels - 1; //note that param.finestLevel is not maxOctreeLevel, but they are actually scaled from 0 to finest level. 
    if (param.minOutputLevel < 0) {
        param.minOutputLevel = 0; //if min output level is below 0, set at 0
    }
    if (param.maxOutputLevel > param.finestLevel) {
        param.maxOutputLevel = param.finestLevel; //if max output level is above finest level, set at finest level.
    }
    

    param.physicalLocation = //the physical location is the full domain 0 position.
        Array<T, 3>(param.fullDomain.x0(), param.fullDomain.y0(), param.fullDomain.z0());

    param.dxCoarsest = param.dxFinest * (T)util::intTwoToThePower(param.finestLevel); //indeed, the scaling is by 2
    param.dtFinest = (param.u_LB / param.u_Ref) * param.dxFinest;  //so dt = ulb/u*dx Note that the best way to change the time steps is by changing uLB, bc uRef is not changeable, and dx finest u want it set at a certian value
    param.dtCoarsest = param.dtFinest * (T)util::intTwoToThePower(param.finestLevel); //the time also scales by 2. 

    param.rho_LB = 1.0;
    param.inletVelocity_LB =
        Array<T, 3>(param.inletVelocity * param.dtFinest / param.dxFinest*std::cos(angle), (T)0, param.inletVelocity * param.dtFinest / param.dxFinest*std::sin(angle)); //this is the inlet velocity, I added the angle of attack.

    param.omega.resize(numLevels); //omega (relaxation factor) is different per level now as can be seen below. 
    for (plint iLevel = 0; iLevel < numLevels; iLevel++) {
        T dx = param.dxFinest * (T)util::intTwoToThePower(param.finestLevel - iLevel);
        T dt = param.dtFinest * (T)util::intTwoToThePower(param.finestLevel - iLevel);
        T nu_LB = param.nu * dt / (dx * dx);
        param.omega[iLevel] = (T)1 / (DESCRIPTOR<T>::invCs2 * nu_LB + (T)0.5); //invCs2 chat says is the inverse of the square of the speed of sound in the lattice. This is some lattice boltzmann theory outside of my scope of understanding
        
    }

    FileName fileName(param.staticSurfaceFileName);
    param.surfaceName = fileName.getName();

    computeAllOutputDomains(param);
}

void printSimulationParameters(SimulationParameters const &param)
{
    pcout << "inletVelocity = " << param.inletVelocity << std::endl;
    pcout << "outflowBcType = " << param.outflowBcType << std::endl;

    for (int iZone = 0; iZone < 6; iZone++) {
        pcout << "spongeWidths[" << iZone << "] = " << param.spongeWidths[iZone] << std::endl;
    }

    pcout << "precision = "
          << (param.precision == FLT
                  ? "FLT"
                  : (param.precision == DBL ? "DBL" : (param.precision == LDBL ? "LDBL" : "INF")))
          << std::endl;
    pcout << "L_Ref = " << param.L_Ref << std::endl;
    pcout << "u_Ref = " << param.u_Ref << std::endl;
    pcout << "u_LB = " << param.u_LB << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;

    pcout << "rho = " << param.rho << std::endl;
    pcout << "nu = " << param.nu << std::endl;
    pcout << "ambientPressure = " << param.ambientPressure << std::endl;

    pcout << "outDir = " << param.outDir << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;

    pcout << "computeAverages = " << (param.computeAverages ? "true" : "false") << std::endl;
    if (param.computeAverages) {
        pcout << "avgIter = " << param.avgIter << std::endl;
    }

    pcout << "cpIter = " << param.cpIter << std::endl;
    pcout << "abIter = " << param.abIter << std::endl;
    pcout << "abortFileName = " << param.abortFileName << std::endl;
    pcout << "xmlContinueFileName = " << param.xmlContinueFileName << std::endl;
    pcout << "baseFileName = " << param.baseFileName << std::endl;
    pcout << "useParallelIO = " << (param.useParallelIO ? "true" : "false") << std::endl;

    pcout << "finestLevel = " << param.finestLevel << std::endl;

    pcout << "inletVelocity_LB = [" << param.inletVelocity_LB[0] << ", "
          << param.inletVelocity_LB[1] << ", " << param.inletVelocity_LB[2] << "]" << std::endl;
    pcout << "Re = " << param.u_Ref * param.L_Ref / param.nu << std::endl;
    pcout << "omegaFinest = " << param.omega[param.finestLevel] << std::endl;
    pcout << "tauFinest = " << (T)1 / param.omega[param.finestLevel] << std::endl;
    pcout << "omegaCoarsest = " << param.omega[0] << std::endl;
    pcout << "tauCoarsest = " << (T)1 / param.omega[0] << std::endl;
    pcout << "dxFinest = " << param.dxFinest << std::endl;
    pcout << "dtFinest / dxFinest = " << param.dtFinest / param.dxFinest << std::endl;
    pcout << "dtFinest / (dxFinest * dxFinest) = "
          << param.dtFinest / (param.dxFinest * param.dxFinest) << std::endl;
    pcout << "physicalLocation = (" << param.physicalLocation[0] << ", "
          << param.physicalLocation[1] << ", " << param.physicalLocation[2] << ")" << std::endl;
    pcout << "incompressibleModel = " << (param.incompressibleModel ? "true" : "false")
          << std::endl;
    pcout << std::endl;
}

void createZones( //creates sponze zones. A sponge zone dampens sound waves so that they don't rebound on the walls and affect the whole simulation. This function is called for each level
    SimulationParameters const &param, MultiBlockLattice3D<T, DESCRIPTOR> &lattice, plint level)
{
    T dx = param.dxFinest * (T)util::intTwoToThePower(param.finestLevel - level);

    Array<plint, 6> numSpongeCells;
    plint totalNumSpongeCells = 0;
    for (plint iZone = 0; iZone < 6; iZone++) {
        numSpongeCells[iZone] = util::roundToInt(param.spongeWidths[iZone] / dx); //determine number of sponge cells per zone by dividing the spongeWidth from config.xml by dx, and rounding it to int.
        totalNumSpongeCells += numSpongeCells[iZone]; //keep track of total number of spongecells.
    }
    // Number of sponge zone lattice nodes at all the outer domain boundaries.
        // So: iZone 0 means the boundary at x = 0
        //     iZone 1 means the boundary at x = nx-1
        //     2 means the boundary at y = 0
        //     and so on... 
    plint nx = util::roundToInt(toLB(param.fullDomain.x1(), 0, dx, param.physicalLocation)) + 1;
    plint ny = util::roundToInt(toLB(param.fullDomain.y1(), 1, dx, param.physicalLocation)) + 1;
    plint nz = util::roundToInt(toLB(param.fullDomain.z1(), 2, dx, param.physicalLocation)) + 1;
    pcout << "nx " << nx << "ny " << ny << "nz " << nz << std::endl;
    Box3D fullBox(0, nx - 1, 0, ny - 1, 0, nz - 1); //at this current level resolution, the total domain is this big. 
    if (totalNumSpongeCells > 0) {
        pcout << "Generating viscosity sponge zone at level: " << level << std::endl;
        bool smartSpongeZone = true; //if true, my function is applied, which makes a cylindrical sponge zone. If false, just makes it cubical sponge zone. 
        T bulkValue = param.omega[level]; //bulk value is the relaxation parameter's for the sponge zone. 
        plint radius = nx*chordLengthPercentage*33; //this is the radius of my cylindrical sponge zone. 33 chord lengths. 
        plint centerX = nx / 2; //center of cylindrical sopnge zone.
        plint centerZ = nz / 2;
        if (smartSpongeZone) { //cylindrical sponge zone function:
            T x0lattice = lattice.getBoundingBox().x0;
            T x1lattice = lattice.getBoundingBox().x1;
            T y0lattice = lattice.getBoundingBox().y0;
            T y1lattice = lattice.getBoundingBox().y1;
            T z0lattice = lattice.getBoundingBox().z0;
            T z1lattice = lattice.getBoundingBox().z1;
            if (x0lattice <= centerX - radius  || x1lattice >= centerX + radius || z0lattice <= centerZ - radius || z1lattice >= centerZ + radius) { //first determine whether it is even necessary to apply the sponge zone to that specific resolution level, indeed, if the whole resolution level is constrained within the cylinder, no need to put sponge zone. 
                //below are two scalar fields, the first one is easy to modify, then will be copied into flagMatrix1 and then that will be used in the spongezone function
                ScalarField3D<int> flagMatrix(nx, ny, nz); 
                MultiScalarField3D<int> flagMatrix1(nx, ny, nz);
                pcout << (T) nx*ny*nz << "This many nodes " << std::endl;

                for (plint iX = 0; iX < nx; ++iX) {
                    for (plint iY = 0; iY < ny; ++iY) {
                        for (plint iZ = 0; iZ < nz; ++iZ) { //loop over all cells.
                            T dxi = (T)(iX - centerX);
                            T dzi = (T)(iZ - centerZ);
                            T distance = std::sqrt(dxi * dxi + dzi * dzi);
                            if (distance < radius) { //determine whether that voxel is inside the cylinder
                                flagMatrix.get(iX, iY, iZ) = 0.0;  // Inside cylinder, no sponge
                            } 
                            else {   
                                flagMatrix.get(iX, iY, iZ) = 1.0;  // Outside cylinder, apply sponge
                            }
                        }
                    }
                }
                copySerializedBlock(flagMatrix, flagMatrix1); //copy flagmatrix into flagmatrix1
                pcout << "This should be larger than 1 " << computeSum(flagMatrix1) << std::endl; //check that it works.
                
                std::vector<MultiBlock3D *> args; //arguments for our MaskedViscositySpongeZone class.
                args.push_back(&lattice);
                args.push_back(&flagMatrix1);

                applyProcessingFunctional(
                    new MaskedViscositySpongeZone3D<T, DESCRIPTOR>(
                        nx, ny, nz, bulkValue, 1, numSpongeCells),  // Note: 1 is the flag value for applying the sponge
                    lattice.getBoundingBox(), args);
                pcout << "hey5" << std::endl; // sometimes i add this prints to check where the code breaks.
            }
        }
        
        
        else { //the cubical normal sponge zone. 
            std::vector<MultiBlock3D *> args;
            args.push_back(&lattice);
            applyProcessingFunctional(
                new ViscositySpongeZone3D<T, DESCRIPTOR>(nx, ny, nz, bulkValue, numSpongeCells),
                lattice.getBoundingBox(), args);
        }
        
    }   
}

void applyOuterBoundaryConditions(
    SimulationParameters const &param, MultiLevelCoupling3D<T, DESCRIPTOR, RESCALER> &lattices,
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc)
{
    Box3D coarsestBoundingBox = lattices.getOgs().getClosedCover(0); //0 is the coarsest level, therefore this gets the coarsest bounding box. 
    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        pcout << "Generating outer domain boundary conditions at level: " << iLevel << std::endl;
        MultiBlockLattice3D<T, DESCRIPTOR> &lattice = lattices.getLevel(iLevel);

        lattice.periodicity().toggleAll(true); //this makes all boundaries periodic. 
    
        lattice.periodicity().toggle(0, false);  //the x direction is still not periodic. 
    
        lattice.periodicity().toggle(2, false);  //the z direction is still not periodic. therefore, only the y direction stays periodic. 
    

        Box3D boundingBox = coarsestBoundingBox.multiply(util::intTwoToThePower(iLevel)); //put bounding box in current level units
        Box3D box = boundingBox;
        T objecty0 = (lattices.getLevel(param.finestLevel).getBoundingBox().y0+5); //this is to implement a freeslip wall right at the tip of the blade as discussed in the presentation
        // objecty0 /= (util::intTwoToThePower(param.finestLevel));
        // objecty0 *= util::intTwoToThePower(iLevel);
        T objecty1 = (lattices.getLevel(param.finestLevel).getBoundingBox().y1-5); //same
        // objecty1 /= (util::intTwoToThePower(param.finestLevel));
        // objecty1 *= util::intTwoToThePower(iLevel);
        objecty0 = int(std::round(objecty0)); //rounding it.
        objecty1 = int(std::round(objecty1));
        pcout << "objecty0" << objecty0 << " objecty1 " << objecty1 << std::endl;
        Box3D inlet(box.x0, box.x0, box.y0, box.y1, box.z0, box.z1);
        Box3D outlet(box.x1, box.x1, box.y0 + 1, box.y1 - 1, box.z0 + 1, box.z1 - 1);
        Box3D yBottom(box.x0 + 1, box.x1, box.y0, box.y0, box.z0, box.z1);
        Box3D yTop(box.x0 + 1, box.x1, box.y1, box.y1, box.z0, box.z1);
        Box3D zBottom(box.x0 + 1, box.x1, box.y0, box.y1, box.z0, box.z0);
        Box3D zTop(box.x0 + 1, box.x1, box.y0, box.y1, box.z1, box.z1);
        Box3D yMarlon0(box.x0+1, box.x1, objecty0,objecty0, box.z0, box.z1 );  //this is where the first freeslip wall is applied. 
        Box3D yMarlon1(box.x0+1, box.x1, objecty1,objecty1, box.z0, box.z1 );
        // Box3D yMarlonedge0
        // Inlet boundary condition.

        bc->setVelocityConditionOnBlockBoundaries(lattice,  boundingBox, inlet,  boundary::dirichlet); //set inlet as dirichlet boundary.

        Array<T, 3> velocity(param.inletVelocity_LB);
        setBoundaryVelocity(lattice, inlet, velocity); //give a value for that dirichlet condition. Note that the velocity is 3 dimnesional
        // bc->setVelocityConditionOnBlockBoundaries(lattice, yMarlon0 ,boundingBox,     boundary::freeslip);
        // bc->setVelocityConditionOnBlockBoundaries(lattice, yMarlon1, boundingBox,     boundary::freeslip);
        bc->addVelocityBoundary1N(yMarlon0, lattice, boundary::freeslip); //make the first freeslip wall
        bc->addVelocityBoundary1P(yMarlon1, lattice, boundary::freeslip);
        // setBoundaryVelocity(lattice, yMarlon0, velocity);
        // setBoundaryVelocity(lattice, yMarlon1, velocity);
        Array<T, 3> zero((T)0, (T)0, (T)0);

        bc->setVelocityConditionOnBlockBoundaries(
            lattice, boundingBox, zBottom, boundary::dirichlet);
        setBoundaryVelocity(lattice, zBottom, velocity); //make zbottom dirichlet, with the same input velocity. 

        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, zTop, boundary::dirichlet);
        setBoundaryVelocity(lattice, zTop, velocity);
        
        // Outlet boundary condition.

        if (param.outflowBcType == 0) { //if 0, set as velocity boundary
            bc->setVelocityConditionOnBlockBoundaries(
                lattice, boundingBox, outlet, boundary::dirichlet);
            setBoundaryVelocity(lattice, outlet, velocity);
        } else if (param.outflowBcType == 1) {
            bc->setPressureConditionOnBlockBoundaries( // if 1 set as pressure boundary
                lattice, boundingBox, outlet, boundary::dirichlet);
            // setBoundaryVelocity(lattice, outlet, velocity);  //interestingly they set boundary velocity to set the pressure. wtf. 
        } else if (param.outflowBcType == 2) { //if 2 set as neumann velocity boundary
            bc->setVelocityConditionOnBlockBoundaries(
                lattice, boundingBox, outlet, boundary::neumann);
            setBoundaryVelocity(lattice, outlet, velocity);
        }
        setBoundaryDensity(lattice, box, param.rho_LB);  //sets boundary density (pressure) mainly for the outlet. 
            // box.y0 -= 2;  // y-periodicity  I deleted this, bc I don't think it's necessary for periodicity. This was taken from externalFlowAroundObstacle example. 
            // box.y1 += 2;     
    }
}

void initializeSimulation(
    SimulationParameters const &param, bool continueSimulation, plint &iniIter,
    MultiLevelCoupling3D<T, DESCRIPTOR, RESCALER> &lattices,
    std::vector<MultiBlock3D *> &checkpointBlocks)
{
    if (!continueSimulation) {//if we are not contnuing the simulation from a file, or from a checkpoint, then initialize normally. 
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            Array<T, 3> velocity(param.inletVelocity_LB);

            MultiBlockLattice3D<T, DESCRIPTOR> &lattice = lattices.getLevel(iLevel);
            initializeAtEquilibrium(lattice, lattice.getBoundingBox(), param.rho_LB,  velocity);//simulation is intialized at equilibrium
        } 
        lattices.initialize();
        lattices.initializeTensorFields();
    } else {
        pcout << std::endl;
        pcout << "Reading state of the simulation from file: " << param.xmlContinueFileName
              << std::endl; 
        loadState(checkpointBlocks, iniIter, param.saveDynamicContent, param.xmlContinueFileName);
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            MultiBlockLattice3D<T, DESCRIPTOR> &lattice = lattices.getLevel(iLevel);
            plint iniIterAtLevel = iniIter * util::intTwoToThePower(iLevel);
            lattice.resetTime(iniIterAtLevel);
        }

        lattices.initializeTensorFields();
        // lattices.initialize();
    }

    pcout << std::endl;
}

void writeResults( //output results.
    SimulationParameters const &param, MultiLevelCoupling3D<T, DESCRIPTOR, RESCALER> &lattices,
    std::unique_ptr<MultiLevelTensorField3D<T, 3> > &avgVel, plint iter)
{
    bool crop = true;
    for (plint iDomain = 0; iDomain < (plint)param.outputDomainNames.size(); iDomain++) {
        std::string fname =
            createFileName(param.outputDomainNames[iDomain] + "_", iter, param.fileNamePadding); //creates file name
        SparseVtkImageOutput3D sparseOut(fname); //vtk output class, it groups the vtk files.

        std::unique_ptr<MultiLevelTensorFieldForOutput3D<T, 3> > velocity = computeVelocity( //for output
            lattices, param.outputDomains.find(param.maxOutputLevel)->second[iDomain], //the second stored value gets the domain, which corresponds to the maxoutput level
            param.maxOutputLevel, crop); //computes the velocity of the whole domain.
        
        
        std::unique_ptr<MultiLevelScalarFieldForOutput3D<T> > density = computeDensity(
            lattices, param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
            param.maxOutputLevel, crop);
        
        std::unique_ptr<MultiLevelTensorFieldForOutput3D<T, 3> > outAvgVel;

        if (param.computeAverages) {
            outAvgVel = exportForOutput(
                *extractSubDomain(
                    *avgVel, param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
                    param.maxOutputLevel),
                param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
                param.maxOutputLevel, crop);
        }

        for (plint iLevel = param.minOutputLevel; iLevel <= param.maxOutputLevel; iLevel++) {
            T dx = param.dxFinest * util::intTwoToThePower(param.finestLevel - iLevel);
            T dt = param.dtFinest * util::intTwoToThePower(param.finestLevel - iLevel);

            T pressureScale = param.rho * (dx * dx) / (dt * dt) * DESCRIPTOR<T>::cs2; //scales the pressure from lattice boltzmann units to pressure. 
            T pressureOffset = param.ambientPressure - param.rho_LB * pressureScale;

            Group3D vtkGroup;
            // You can add scalar- and tensor-fields to the group with the usual "group.add()". The
            // advantage of addTransform is that is also converts the type to float and multiplies
            // by a scale factor (and adds an offset).
            addTransform<T, float, 3>(vtkGroup, velocity->getLevel(iLevel), "velocity", dx / dt); //scale factor dx/dt
            addTransform<T, float>(
                vtkGroup, density->getLevel(iLevel), "pressure", pressureScale, pressureOffset); 

            if (param.computeAverages) {
                addTransform<T, float, 3>(vtkGroup, outAvgVel->getLevel(iLevel), "avgVel", dx / dt);
            }

            // "pointData = true" is the usual VTK output. "pointData = false" inhibits
            // interpolations, and is useful for debugging.
            bool pointData = false;
            sparseOut.writeVtkBlock(vtkGroup, dx, param.physicalLocation, iLevel, pointData);
        }
    }
}

int main(int argc, char *argv[])  //the main function. 
{
    plbInit(&argc, &argv); //intialize. 

    std::cout.precision(10);//set output precision for numbers. 
    std::scientific(std::cout);

    // Command-line arguments

    if (argc != 2 && argc != 3) {
        pcerr << "Usage: " << argv[0] << " xml-input-file-name [restart]" << std::endl; //if you don't specify an input file, it asks you to do so. 
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);
    abortIfCannotOpenFileForReading(xmlInputFileName); //if it can't open the input file, abort. 

    bool continueSimulation = false;
    if (argc == 3) {
        std::string cmd(argv[2]);
        if (cmd == "restart") {
            continueSimulation = true; //if restart is written it understands that it takes the input from the checkpoint file
        }
    }

    int nproc = global::mpi().getSize(); //get number of processors. specified in teh terminal by mpirun -np 4 ./jordiPowerFlowCopy config.xml   | , 4 is the number of processors

    // global::mpi().barrier();
    global::timer("init").start();

    // Set the simulation parameters.

    SimulationParameters param;

    readUserDefinedSimulationParameters(xmlInputFileName, param);

    if (continueSimulation) {
        abortIfCannotOpenFileForReading(param.xmlContinueFileName);
    }

    createOctreeGridStructure(param);
    calculateDerivedSimulationParameters(param);

    global::directories().setOutputDir(param.outDir);
    global::IOpolicy().activateParallelIO(param.useParallelIO);

    if (nproc != param.ogs.getNumProcesses()) {
        pcerr << "The number of processes used is not the same as the one provided in the "
                 "grid-structure files."
              << std::endl;
        exit(1);
    }

    // The "order" is about how Palabos rescales the populations. For now, we work at order 0.
    // RESCALER means: we use convective scaling,
    // Palabos scales only the populations and no external scalars, and
    // this works for BGK but not for MRT.
    plint order = 0;

    pcout << std::endl;
    pcout << "Generating the lattices." << std::endl;

    order = 1;
    MultiLevelCoupling3D<T, DESCRIPTOR, RESCALER> lattices( //multilevelcoupling is a class that can couple multiple lattices. 
        param.ogs, //param.ogs is the octree grid structure. 
        new ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T, DESCRIPTOR>(  //collision dynamics for lattice boltzmann.
            param.omega[0], 0.14), //for some reason they use the first omega and the 0.14 smagorinksy factor
        order);
    pcout << "CompleteRegularizedBGKdynamics" << std::endl;
    param.incompressibleModel = false;

    for (plint iLevel = 0; iLevel < lattices.getNumLevels(); iLevel++) {
        pcout << "Info for lattice at level: " << iLevel << std::endl;
        pcout << getMultiBlockInfo(lattices.getLevel(iLevel)) << std::endl;
    }

    printSimulationParameters(param);

    // Immersed surfaces only for finest level.

    pcout << "Reading the immersed surface geometries." << std::endl;
    pcout << "Generating fluid blocks at finest level." << std::endl;

    std::vector<MultiBlock3D *> rhoBarJarg; 
    // Here we assume that the object is far from the outlet. This is why we can compute
    // the rhoBarJ field to be used with the off-lattice BC, before the FluidPressureOutlet3D
    // is executed.
    plint numScalars = 4;
    plint extendedTensorEnvelopeWidth = 2;  // Extrapolated BCs. guo = 2, generalized = 3
    MultiNTensorField3D<T> *rhoBarJfield = generateMultiNTensorField3D<T>( //multidimensional tensor field. 
        lattices.getLevel(param.finestLevel), extendedTensorEnvelopeWidth, numScalars);
    rhoBarJfield->toggleInternalStatistics(false);
    rhoBarJarg.push_back(rhoBarJfield); // add the field to the vector. 
    integrateProcessingFunctional(
        new PackedRhoBarJfunctional3D<T, DESCRIPTOR>(),
        lattices.getLevel(param.finestLevel).getBoundingBox(), lattices.getLevel(param.finestLevel),
        *rhoBarJfield, 0); // add a processing functional to calculate the rhobharJ

    pcout << "Implementing the geometry at level." << std::endl;

    // The next few lines of code are typical. They transform the surface geometry of the
    //   stl given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    TriangleSet<T> triangleSet(param.staticSurfaceFileName, DBL); //read the stl file
    triangleSet.translate(-param.physicalLocation); //translate it to the physical Location (lowerleft corner of the domain)
    triangleSet.scale((T)1 / param.dxFinest); //put it in dx finest units, bc now 1 corresponds to dx.
    
    Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    Array<T, 3> obstacleCenter = (T)0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner); //gets center of obstacle
    T scaling_factor = (lattices.getLevel(0).getBoundingBox().x1-lattices.getLevel(0).getBoundingBox().x0)/(-bCuboid.lowerLeftCorner[0] + bCuboid.upperRightCorner[0])*chordLengthPercentage*(T)util::intTwoToThePower(param.finestLevel); // this is a scaling factor to set the chord length as the desired chordLengthPercentage of the whole domain length. 
    triangleSet.scale(scaling_factor); //rescales all dimensions by the scaling factor. 
    Cuboid<T> bCuboidy = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    T yscaling_factor = (lattices.getLevel(0).getBoundingBox().y1-lattices.getLevel(0).getBoundingBox().y0)/(-bCuboidy.lowerLeftCorner[1] + bCuboidy.upperRightCorner[1])*(T)util::intTwoToThePower(param.finestLevel); //this is to elongate the stl file to cover from y0 to y1
    // T yscaling_factor = (lattices.getLevel(param.finestLevel).getBoundingBox().y1-lattices.getLevel(param.finestLevel).getBoundingBox().y0)/(-bCuboidy.lowerLeftCorner[1] + bCuboidy.upperRightCorner[1]);
    // T scaling_factor2 = (lattices.getLevel(param.finestLevel).getBoundingBox().y1-lattices.getLevel(param.finestLevel).getBoundingBox().y0)/(-bCuboid.lowerLeftCorner[1] + bCuboid.upperRightCorner[1]);
    // triangleSet.scale((T)0.0,(T) scaling_factor2*1.05 , (T) 0.0);
    // pcout << "El objetoooooooooooooo x0 " << bCuboidy.lowerLeftCorner[2] << " a " << bCuboidy.upperRightCorner[2] << std::endl;
    triangleSet.scale((T) 1.0, yscaling_factor, (T) 1.0); //scale the y dimension
    Cuboid<T> bCuboid2 = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objetoooooooooooooo x0 " << bCuboid2.lowerLeftCorner[2] << " a " << bCuboid2.upperRightCorner[2] << std::endl;
    Array<T, 3> obstacleCenter2 = (T)0.5 * (bCuboid2.lowerLeftCorner + bCuboid2.upperRightCorner); //gets center of obstacle
    if (teaddon == true) //when this is true, the obstacle center is recalculated bc we don't to center the domain at the stl's center, but at the airfoil center, (without the add on)
    {
        T cleanChordLength = (bCuboid2.upperRightCorner[0]-bCuboid2.lowerLeftCorner[0])*0.2/0.2335; //this is based on the Marlon object that we are analysing. 
        T actualCenterWithoutTail = (T)0.5 * (bCuboid2.lowerLeftCorner[0] + cleanChordLength);
        obstacleCenter2[0] = actualCenterWithoutTail;   
    }
    //fixed to center considering the tail does not count in the center calculation. 
    triangleSet.translate(-obstacleCenter2); 
    Cuboid<T> bCuboidy2 = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objetoooooooooooooo x0 " << bCuboidy2.lowerLeftCorner[2] << " a " << bCuboidy2.upperRightCorner[2] << std::endl;
    triangleSet.translate(Array<T,3>((lattices.getLevel(param.finestLevel).getBoundingBox().x1+lattices.getLevel(param.finestLevel).getBoundingBox().x0)*0.5, (lattices.getLevel(param.finestLevel).getBoundingBox().y1+lattices.getLevel(param.finestLevel).getBoundingBox().y0)*0.5, (lattices.getLevel(param.finestLevel).getBoundingBox().z1+lattices.getLevel(param.finestLevel).getBoundingBox().z0)/2)); //centers it at 0,0,0 of finest level
    Cuboid<T> bCuboidy3 = triangleSet.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy3.lowerLeftCorner[0] << " a " << bCuboidy3.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy3.lowerLeftCorner[1] << " a " << bCuboidy3.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy3.lowerLeftCorner[2] << " a " << bCuboidy3.upperRightCorner[2] << std::endl;
    
    // triangleSet.scale((T)5000);
     
    TriangleSet<T> triangleSet2;
    TriangleSet<T> triangleSetx;
    TriangleSet<T> triangleSetx2;
    TriangleSet<T> triangleSetz;
    TriangleSet<T> triangleSetz2;
    TriangleSet<T> triangleSet3;
    
    // Plane<T> planeyminus(Array<T, 3>(0., lattices.getLevel(param.finestLevel).getBoundingBox().y0+5, 0.),Array<T, 3>(0., -1., 0.) );
    Plane<T> planeyminus(Array<T, 3>(0., lattices.getLevel(param.finestLevel).getBoundingBox().y0+5, 0.),Array<T, 3>(0., -1., 0.) ); //create a cutting plane at the y0 of the finest level.
    pcout << "planeyminus " << lattices.getLevel(param.finestLevel).getBoundingBox().y0+5 << std::endl;
    // Plane<T> planeyplus(Array<T, 3>(0., lattices.getLevel(param.finestLevel).getBoundingBox().y1-5, 0.),Array<T, 3>(0., 1., 0.) );
    Plane<T> planeyplus(Array<T, 3>(0., (lattices.getLevel(param.finestLevel).getBoundingBox().y1-5), 0.),Array<T, 3>(0., 1., 0.) ); //same for y1
    pcout << "planeyplus " << (lattices.getLevel(param.finestLevel).getBoundingBox().y1-5) << std::endl;
    Plane<T> planexminus(Array<T, 3>(lattices.getLevel(param.finestLevel).getBoundingBox().x0+5, 0., 0.),Array<T, 3>(-1., 0., 0.) );
    pcout << "planexminus " << lattices.getLevel(param.finestLevel).getBoundingBox().x0+5 << std::endl;
    Plane<T> planexplus(Array<T, 3>(lattices.getLevel(param.finestLevel).getBoundingBox().x1-5, 0., 0.),Array<T, 3>(1., 0., 0.) );
    pcout << "planexplus " << lattices.getLevel(param.finestLevel).getBoundingBox().x1-5 << std::endl;
    Plane<T> planezminus(Array<T, 3>(0., 0., lattices.getLevel(param.finestLevel).getBoundingBox().z0+5),Array<T, 3>(0., 0., -1.) );
    pcout << "planezminus " << lattices.getLevel(param.finestLevel).getBoundingBox().z0+5 << std::endl;
    Plane<T> planezplus(Array<T, 3>(0., 0., lattices.getLevel(param.finestLevel).getBoundingBox().z1-5),Array<T, 3>(0., 0., 1.) );
    pcout << "planezplus " << lattices.getLevel(param.finestLevel).getBoundingBox().z1-5 << std::endl;
    // Plane<T> planeyminus(Array<T, 3>(0., lattices.getLevel(0).getBoundingBox().y0*(T)util::intTwoToThePower(param.finestLevel)+1, 0.),Array<T, 3>(0., -1., 0.) );
    // Plane<T> planeyplus(Array<T, 3>(0., lattices.getLevel(0).getBoundingBox().y1*(T)util::intTwoToThePower(param.finestLevel)-1, 0.),Array<T, 3>(0., 1., 0.) );
    if (bCuboidy3.lowerLeftCorner[1] < lattices.getLevel(param.finestLevel).getBoundingBox().y0+5) { //if the object overpasses the plane
        triangleSet.cutWithPlane(planeyminus, triangleSet2); //crop crop
    }
    else {
        triangleSet2 = triangleSet; //otherwise don't change it .
    }
    Cuboid<T> bCuboidy4 = triangleSet2.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy4.lowerLeftCorner[0] << " a " << bCuboidy4.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy4.lowerLeftCorner[1] << " a " << bCuboidy4.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy4.lowerLeftCorner[2] << " a " << bCuboidy4.upperRightCorner[2] << std::endl;
    if (bCuboidy4.upperRightCorner[1] > lattices.getLevel(param.finestLevel).getBoundingBox().y1-5) { //if object passes plane
        triangleSet2.cutWithPlane(planeyplus,triangleSetx);  //crop crop
    }
    else {
        triangleSetx = triangleSet2;
    }
    
    Cuboid<T> bCuboidy5 = triangleSetx.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy5.lowerLeftCorner[0] << " a " << bCuboidy5.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy5.lowerLeftCorner[1] << " a " << bCuboidy5.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy5.lowerLeftCorner[2] << " a " << bCuboidy5.upperRightCorner[2] << std::endl;
    // if (bCuboidy5.lowerLeftCorner[0] < lattices.getLevel(param.finestLevel).getBoundingBox().x0+5) {
    if (9>10) { //this is disabled by force (as 9 is never larger than 10) because of some bugs
        triangleSetx.cutWithPlane(planexminus,triangleSetx2);
    }
    else {
        triangleSetx2 = triangleSetx;
    }
    
    Cuboid<T> bCuboidy6 = triangleSetx2.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy6.lowerLeftCorner[0] << " a " << bCuboidy6.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy6.lowerLeftCorner[1] << " a " << bCuboidy6.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy6.lowerLeftCorner[2] << " a " << bCuboidy6.upperRightCorner[2] << std::endl;
    // if (bCuboidy6.upperRightCorner[0] > lattices.getLevel(param.finestLevel).getBoundingBox().x1-5) {
    if (9>10) { //disabled
        triangleSetx2.cutWithPlane(planexplus,triangleSetz);
    }
    else {
        triangleSetz = triangleSetx2;
    }
    
    Cuboid<T> bCuboidy7 = triangleSetz.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto below it will explode x0 " << bCuboidy7.lowerLeftCorner[0] << " a " << bCuboidy7.upperRightCorner[0] << std::endl;
    // pcout << "El objeto below it will explode y0 " << bCuboidy7.lowerLeftCorner[1] << " a " << bCuboidy7.upperRightCorner[1] << std::endl;
    // pcout << "El objeto below it will explode z0 " << bCuboidy7.lowerLeftCorner[2] << " a " << bCuboidy7.upperRightCorner[2] << std::endl;
    // if (bCuboidy7.lowerLeftCorner[2] < lattices.getLevel(param.finestLevel).getBoundingBox().z0+5) {
    if (9>10) { //disabled
        triangleSetz.cutWithPlane(planezminus,bCuboidy7,triangleSetz2);
        pcout << "ep" << std::endl;
    }
    else {
        triangleSetz2 = triangleSetz;
    }
     //here x it multiplies by a lot
    Cuboid<T> bCuboidy8 = triangleSetz2.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy8.lowerLeftCorner[0] << " a " << bCuboidy8.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy8.lowerLeftCorner[1] << " a " << bCuboidy8.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy8.lowerLeftCorner[2] << " a " << bCuboidy8.upperRightCorner[2] << std::endl;
    // if (bCuboidy7.upperRightCorner[2] > lattices.getLevel(param.finestLevel).getBoundingBox().z1-5) {
    if (9>10) { //disabled
        triangleSetz2.cutWithPlane(planezplus,bCuboidy7,triangleSet3);
    }
    else {
        triangleSet3 = triangleSetz2;
    }
    
    Cuboid<T> bCuboidy9 = triangleSet3.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy9.lowerLeftCorner[0] << " a " << bCuboidy9.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy9.lowerLeftCorner[1] << " a " << bCuboidy9.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy9.lowerLeftCorner[2] << " a " << bCuboidy9.upperRightCorner[2] << std::endl;
    Cuboid<T> bCuboid3x = triangleSet3.getBoundingCuboid();
    // pcout << "El objetoooooooooooooo x0 " << bCuboid3x.lowerLeftCorner[2] << " a " << bCuboid3x.upperRightCorner[2] << std::endl;
    Array<T, 3> obstacleCenter3 = (T)0.5 * (bCuboid3x.lowerLeftCorner + bCuboid3x.upperRightCorner); //gets center of obstacle
    if (teaddon == true)
    {
        T cleanChordLength = (bCuboid3x.upperRightCorner[0]-bCuboid3x.lowerLeftCorner[0])*0.2/0.2335; //this is based on the Marlon object that we are analysing. 
        T actualCenterWithoutTail = (T)0.5 * (bCuboid3x.lowerLeftCorner[0] + cleanChordLength);
        obstacleCenter3[0] = actualCenterWithoutTail;   
    }
    //fixed to center considering the tail does not count in the center calculation. 
    // triangleSet3.translate(-obstacleCenter3); //centers it at 0,0,0 i think. indeed, tested below. 
    Cuboid<T> bCuboidy10 = triangleSet3.getBoundingCuboid();  //gets bounding cuboid
    // pcout << "El objeto x0 " << bCuboidy10.lowerLeftCorner[0] << " a " << bCuboidy10.upperRightCorner[0] << std::endl;
    // pcout << "El objeto y0 " << bCuboidy10.lowerLeftCorner[1] << " a " << bCuboidy10.upperRightCorner[1] << std::endl;
    // pcout << "El objeto z0 " << bCuboidy10.lowerLeftCorner[2] << " a " << bCuboidy10.upperRightCorner[2] << std::endl;
    // triangleSet3.translate(Array<T,3>((lattices.getLevel(param.finestLevel).getBoundingBox().x1+lattices.getLevel(param.finestLevel).getBoundingBox().x0)*0.5, (lattices.getLevel(param.finestLevel).getBoundingBox().y1+lattices.getLevel(param.finestLevel).getBoundingBox().y0)*0.5,  (lattices.getLevel(param.finestLevel).getBoundingBox().z1+lattices.getLevel(param.finestLevel).getBoundingBox().z0)/2)); //centers it at 0,0,0 i think. indeed, tested below. 
    // triangleSet3.scale((T) 0.5);
    Cuboid<T> bCuboid3 = triangleSet3.getBoundingCuboid();
    pcout << "Hello, testing pcout " << std::endl;
    pcout << "El objetoooooooooooooo x0 " << bCuboid3.lowerLeftCorner[0] << " a " << bCuboid3.upperRightCorner[0] << std::endl;
    pcout << "Lattice x0 " << lattices.getLevel(param.finestLevel).getBoundingBox().x0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().x1 << std::endl;
    pcout << "Domain x0 " << lattices.getLevel(0).getBoundingBox().x0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().x1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    pcout << "El objetoo y0 " << bCuboid3.lowerLeftCorner[1] << " a " << bCuboid3.upperRightCorner[1] << std::endl;
    pcout << "Lattice finest -1 Y0 " << lattices.getLevel(param.finestLevel-1).getBoundingBox().y0 << "a " << lattices.getLevel(param.finestLevel-1).getBoundingBox().y1 << std::endl;
    pcout << "Lattice finest Y0 " << lattices.getLevel(param.finestLevel).getBoundingBox().y0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().y1 << std::endl;
    pcout << "Domain Y0 " << lattices.getLevel(0).getBoundingBox().y0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().y1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    pcout << "El objetoo z0 " << bCuboid3.lowerLeftCorner[2] << " a " << bCuboid3.upperRightCorner[2] << std::endl;
    pcout << "Lattice z0 " << lattices.getLevel(param.finestLevel).getBoundingBox().z0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().z1 << std::endl;
    pcout << "Domain Y0 " << lattices.getLevel(0).getBoundingBox().z0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().z1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    // triangleSet3.scale((T) 1.0 , (T) 3.0 , (T) 1.0);
    // Cuboid<T> bCuboid4 = triangleSet3.getBoundingCuboid();
    // pcout << "Hello, testing pcout " << std::endl;
    // pcout << "El objetoooooooooooooo x0 " << bCuboid4.lowerLeftCorner[0] << " a " << bCuboid4.upperRightCorner[0] << std::endl;
    // pcout << "Lattice x0 " << lattices.getLevel(param.finestLevel).getBoundingBox().x0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().x1 << std::endl;
    // pcout << "Domain x0 " << lattices.getLevel(0).getBoundingBox().x0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().x1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    // pcout << "El objetoo y0 " << bCuboid4.lowerLeftCorner[1] << " a " << bCuboid4.upperRightCorner[1] << std::endl;
    // pcout << "Lattice Y0 " << lattices.getLevel(param.finestLevel).getBoundingBox().y0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().y1 << std::endl;
    // pcout << "Domain Y0 " << lattices.getLevel(0).getBoundingBox().y0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().y1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    // pcout << "El objetoo z0 " << bCuboid4.lowerLeftCorner[2] << " a " << bCuboid4.upperRightCorner[2] << std::endl;
    // pcout << "Lattice z0 " << lattices.getLevel(param.finestLevel).getBoundingBox().z0 << "a " << lattices.getLevel(param.finestLevel).getBoundingBox().z1 << std::endl;
    // pcout << "Domain Y0 " << lattices.getLevel(0).getBoundingBox().z0*(T)util::intTwoToThePower(param.finestLevel) << "a " << lattices.getLevel(0).getBoundingBox().z1*(T)util::intTwoToThePower(param.finestLevel) << std::endl;
    // triangleSet3.scale((T) 0.5);
    // 
    DEFscaledMesh<T> defMesh(triangleSet3, 0, 0, 1, Dot3D(0, 0, 0)); //defMesh defines a mesh around the object
    defMesh.setDx(param.dxFinest);
    defMesh.setPhysicalLocation(param.physicalLocation); 
    pcout << "hey1 " << std::endl;
    TriangleBoundary3D<T> boundary(defMesh); //makes it a boundary
    pcout << "hey2 " << std::endl;
    boundary.getMesh().inflate(0.05); 
    pcout << "hey3 " << std::endl;

    //A voxelization process is used to determine which nodes are inside of the object, and which nodes are outside
    const int flowType = voxelFlag::outside;  //set the flag as outside (flow is outside)
    const int borderWidth = 1;
    const int extendedEnvelopeWidth = 2;
    const int blockSize = 0;
    VoxelizedDomain3D<T> voxelizedDomain( //voxelize the domain
        boundary, flowType, lattices.getLevel(param.finestLevel).getBoundingBox(), borderWidth,
        extendedEnvelopeWidth, blockSize);
    voxelizedDomain.reparallelize(param.ogs.getMultiBlockManagement( 
        param.finestLevel, lattices.getLevel(param.finestLevel).getBoundingBox(),
        extendedEnvelopeWidth)); //im not sure what reparallelize is for, maybe to insert the multiblockmanagement thing

    defineDynamics(
        lattices.getLevel(param.finestLevel), voxelizedDomain.getVoxelMatrix(),
        lattices.getLevel(param.finestLevel).getBoundingBox(), new NoDynamics<T, DESCRIPTOR>((T)1),
        voxelFlag::inside); //defines dynamics inside as NoDynamics (flow inside does not move )

    boundary.getMesh().writeAsciiSTL(param.outDir + "flowMesh.stl");

    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    // The boundary condition algorithm for the object.
    BoundaryProfiles3D<T, Velocity> profiles; 
    profiles.setWallProfile(new NoSlipProfile3D<T>()); //set no slip boundary condition for object

    // Filipova BC needs no dynamics inside the voxel object
    defineDynamics(
        lattices.getLevel(param.finestLevel), voxelizedDomain.getVoxelMatrix(),
        lattices.getLevel(param.finestLevel).getBoundingBox(), new NoDynamics<T, DESCRIPTOR>((T)1), //1 is the density 
        voxelFlag::innerBorder); // at the inner border also no dynamics

   
    // GuoOffLatticeModel3D< T, DESCRIPTOR > *model =
    //     new GuoOffLatticeModel3D<T, DESCRIPTOR>(
    //         new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles),
    //         flowType);
    FilippovaHaenelLocalModel3D<T, DESCRIPTOR> *model = 
        new FilippovaHaenelLocalModel3D<T, DESCRIPTOR>(
            new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles),
            flowType); //this is the filippova off lattice boundary condition model.
    // BouzidiOffLatticeModel3D<T, DESCRIPTOR> *model =
    //     new BouzidiOffLatticeModel3D<T, DESCRIPTOR>(
    //         new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles),
    //         flowType);

    OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> *boundaryCondition = //apply the boundary condition
        new OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity>(
            model->clone(), voxelizedDomain, lattices.getLevel(param.finestLevel));
    boundaryCondition->insert(rhoBarJarg); //insert the rhoBarJ field.
    PLB_ASSERT(boundaryCondition != 0);
    delete model;
    model = 0;

    if (!continueSimulation) {
        pcout << "Generating outer domain zones." << std::endl;
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            createZones(param, lattices.getLevel(iLevel), iLevel);
        }
    }

    pcout << "Generating outer domain boundary conditions." << std::endl;

    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *bc = //make on lattice boundary conditions (domain boundary conditions.)
        createInterpBoundaryCondition3D<T, DESCRIPTOR>();
    applyOuterBoundaryConditions(param, lattices, bc);
    delete bc;
    bc = 0;

    // Generation of statistics MultiLevel fields.

    std::unique_ptr<MultiLevelTensorField3D<T, 3> > vel;
    std::unique_ptr<MultiLevelTensorField3D<T, 3> > avgVelocity;  //these fields are for output of data .
    

    std::vector<MultiLevel3D *> avgVelocityArgs;

    Box3D coarsestBoundingBox = lattices.getOgs().getClosedCover(0);
    if (param.computeAverages) { 
        vel = generateMultiLevelTensorField3D<T, 3>(
            lattices.getOgs(), coarsestBoundingBox, 0, Array<T, 3>(0.0, 0.0, 0.0));
        avgVelocity = generateMultiLevelTensorField3D<T, 3>(
            lattices.getOgs(), coarsestBoundingBox, 0, Array<T, 3>(0.0, 0.0, 0.0));

        avgVelocityArgs.push_back(vel.get());
        avgVelocityArgs.push_back(avgVelocity.get());
    }

    // Initialization.

    std::vector<MultiBlock3D *> checkpointBlocks;

    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        checkpointBlocks.push_back(&lattices.getLevel(iLevel));
    }

    if (param.computeAverages) {
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            checkpointBlocks.push_back(&avgVelocity->getLevel(iLevel));
        }
    }

    plint iniIter = 0;

    initializeSimulation(param, continueSimulation, iniIter, lattices, checkpointBlocks); 

    // Use "collideAndStream" at all levels except the finest one, at
    // which "executeInternalProcessors" is used instead.
    std::map<plint, bool> useExecuteInternalProcessors;
    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        useExecuteInternalProcessors[iLevel] = false; 
    }

    // Prepare files.
    std::string fileName = param.outDir + "average_energy_finest_level.dat";
    plb_ofstream energyFile(
        fileName.c_str(), continueSimulation ? std::ofstream::app : std::ofstream::out);

    // integration of point measures and creation of output files if needed
    std::vector<std::vector<plint> > ids;
    std::vector<plb_ofstream *> probesFileNames(param.finestLevel + 1);
    std::vector<std::vector<std::vector<T> > > results(param.finestLevel + 1);
    const plint statsId = -200;

    bool computeStats = false;

    fileName = param.outDir + param.surfaceName + "_total_force.dat";
    plb_ofstream forces(
        fileName.c_str(), continueSimulation ? std::ofstream::app : std::ofstream::out);

    // Starting iterations.
    global::timer("init").stop();

    pcout << "The full initialization phase took " << global::timer("init").getTime()
          << " seconds on " << nproc << " processes." << std::endl;
    pcout << std::endl;
    pcout << "rho " << param.rho << " inlet velocity " << param.inletVelocity << " x0 " << param.fullDomain.lowerLeftCorner[0] << " x1 " << param.fullDomain.upperRightCorner[0] << " y0 " << param.fullDomain.lowerLeftCorner[1] << " y1 " << param.fullDomain.upperRightCorner[1] << std::endl;
    
    pcout << "Starting simulation." << std::endl;
    pcout << std::endl;
    bool stopExecution = false;
    bool checkForErrors = true;
    plint iter = iniIter;
    bool avgProcIntegrated = false;
    std::vector<plint> extProcFunIds;
    for (; iter < param.maxIter && !stopExecution; iter++) { //here begins the simulation. 
        if (iter % param.statIter == 0 && iter != 0) {
            pcout << "At coarsest level iteration: " << iter << ", t = " << iter * param.dtCoarsest
                  << std::endl;
            T energy = boundaryCondition->computeAverageEnergy() * param.rho //computes kinetic energy
                       * (param.dxFinest * param.dxFinest) / (param.dtFinest * param.dtFinest);

            pcout << "Average kinetic energy at the finest level: " << energy << std::endl;
            energyFile << (double)(iter * param.dtCoarsest) << " " << (double)energy << std::endl;

            // Forces on immersed surfaces.

            T forceConversion = //there needs to be a force conversion bc we go from lattice units (basically unitless) to SI units
                param.rho * (param.dxFinest * param.dxFinest * param.dxFinest * param.dxFinest)
                / (param.dtFinest * param.dtFinest);
            Array<T, 3> force = forceConversion * boundaryCondition->getForceOnObject();
            forces << (double)(iter * param.dtCoarsest) << " " << force[2]*std::cos(angle)-force[0]*std::sin(angle)/(0.5*param.rho*param.inletVelocity*param.inletVelocity*(-param.fullDomain.lowerLeftCorner[0]+param.fullDomain.upperRightCorner[0])*chordLengthPercentage*(param.fullDomain.upperRightCorner[1]-param.fullDomain.lowerLeftCorner[1])) << " " << force[0]*std::cos(angle)+force[2]*std::sin(angle)/(0.5*param.rho*param.inletVelocity*param.inletVelocity*(-param.fullDomain.lowerLeftCorner[0]+param.fullDomain.upperRightCorner[0])*chordLengthPercentage*(param.fullDomain.upperRightCorner[1]-param.fullDomain.lowerLeftCorner[1])) << " " << force[1] << std::endl; //time, cl, cd, lateral force. because lift is perpendicular to airflow, and drag parallel
            //this long line above calculated lift coefficient, drag coefficient and lateral force.  (the coordinate system is at an angle with the flow.)
            //Cl = (Fz*cos(aoa) - Fx*sin(aoa))/(1/2*rho*V**2*S) .  Cd = (Fx*cos(aoa) + Fz*sin(aoa))/(1/2*rho*V**2*S)
            if (iter > 0) {
                T totTime = global::timer("lb-iter").getTime();
                pcout << "Total coarsest level iteration: " << totTime / (T)iter << std::endl;
            }

            pcout << std::endl;
        }

        // With grid-refinement, the order of the integration of data processors is different.
        // ExternalRhoJcollideAndStream3D is integrated at the end, and not at the beginning
        // of the cycle, as it happens in codes that do not use grid-refinement. This means
        // that in these cases the boundary conditions are imposed "before" collide-and-stream
        // in the cycle. So, the results are exported before the BCs are enforced, and this
        // might cause visualization problems.
        if (iter % param.outIter == 0) {
            pcout << "Output to disk at coarsest level iteration: " << iter
                  << ", t = " << iter * param.dtCoarsest << std::endl;
            // MultiLevelTensorField3D<T,3> vels = vel;
            writeResults(param, lattices, avgVelocity, iter);
            pcout << "Output Completed" << std::endl;
        }

        if ((param.cpIter > 0 && iter % param.cpIter == 0 && iter != iniIter)
            || iter == param.maxIter - 1) {
            pcout << "Saving the state of the simulation at coarsest level iteration: " << iter
                  << std::endl;
            saveState(
                checkpointBlocks, iter, param.saveDynamicContent, param.xmlContinueFileName,
                param.baseFileName, param.fileNamePadding);
            pcout << std::endl;
        }

        if (iter % param.abIter == 0) {
            stopExecution = abortExecution(
                param.abortFileName, checkpointBlocks, iter, param.saveDynamicContent,
                param.xmlContinueFileName, param.baseFileName, param.fileNamePadding);

            if (stopExecution) {
                pcout << "Aborting execution at iteration: " << iter << std::endl;
                pcout << std::endl;
            }
        }

        global::timer("lb-iter").start();
        global::timer("solver").start();
        // for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        //     pcout << "useExecuteInteralProcessors"<<useExecuteInternalProcessors[iLevel] << std::endl;
        // }
        
        lattices.collideAndStream(
            0, useExecuteInternalProcessors, extProcFunIds, computeStats, statsId, ids, results);
        global::timer("solver").stop();
        global::timer("lb-iter").stop();

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }

        if (param.computeAverages
            && (iter == param.avgIter
                || (iter > param.avgIter && continueSimulation && !avgProcIntegrated)))
        {
            avgProcIntegrated = true;
            plint compFields = -100;
            plint compAvgs = compFields - 1;
            integrateProcessingFunctional(
                new BoxVelocityFunctional3D<T, DESCRIPTOR>(), coarsestBoundingBox, 0, lattices,
                *vel, compFields);
            integrateProcessingFunctional(
                new UpdateAveTensorTransientStatistics3D<T, 3>(iter - param.avgIter + 1),
                coarsestBoundingBox, 0, lattices, avgVelocityArgs, param.ogs.getNumLevels(),
                compAvgs);

            extProcFunIds.push_back(compFields);
            extProcFunIds.push_back(compAvgs);
        }
    }

    pcout << "The " << iter - iniIter << " iterations at the coarsest level, took "
          << global::timer("solver").getTime() << " seconds on " << nproc << " processes."
          << std::endl;

    plb_ofstream summary("execution_summary.txt");
    summary << "Summary of execution of the solver: " << argv[0] << std::endl;
    summary << "Number of processes: " << nproc << std::endl;
    summary << "Total time of the initialization phase           : "
            << (double)global::timer("init").getTime() << " s" << std::endl;
    summary << "Total time of the pure solution phase (no output): "
            << (double)global::timer("solver").getTime() << " s" << std::endl;
    summary.close();
    energyFile.close();

    delete boundaryCondition;
    delete rhoBarJfield;

    return 0;
}
