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

#include "HWLatticeModel3D.h"
#include "HWLatticeModel3D.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "setupMacroscopiFields.h"
#include "simulationParameters.h"
#include "magic_enum.hpp"
#include "exprtk/exprtk.hpp"
using namespace plb;
using namespace plb::descriptors;
namespace lu = incompressible_simulation_parameters;

#ifndef CURVED_BOUNDARIES_SETUP_H
#define CURVED_BOUNDARIES_SETUP_H

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>     expression_t;
typedef exprtk::parser<double>             parser_t;


class BSchemeInfo{
public:
    explicit BSchemeInfo(const std::string &boundaryXmlPathFromRoot);
    //    XMLreader boundaryXml("LIplusGeneric.xml");
    //    std::string main_expression_str{};
    //    symbol_table_t main_table;
    //    expression_t main_expression;
    parser_t main_parser;
    symbol_table_t symbolTable;

    std::string alphaMin_expression_str;
    expression_t alphaMin_expression;

    std::string alphaPlus_expression_str;
//    symbol_table_t alphaPlus_table;
    expression_t alphaPlus_expression;

    std::string beta_expression_str;
//    symbol_table_t beta_table;
    expression_t beta_expression;

    std::string kMin_expression_str;
//    symbol_table_t kMin_table;
    expression_t kMin_expression;

    std::string kPlus_expression_str;
//    symbol_table_t kPlus_table;
    expression_t kPlus_expression;

    double q{}, up{}, down{};
    double Lambda_min{}, Lambda_plus{},tau_plus{}, tau_min{};
    double magic{};
    void parse(symbol_table_t& table,expression_t& expression, std::string& expression_str, parser_t& parser) {
        table.add_variable("q", q);
        table.add_variable("LambdaMinus", Lambda_min);
        table.add_variable("LambdaPlus", Lambda_plus);
        table.add_variable("tauPlus", tau_plus);
        table.add_variable("tauMin", tau_min);
        expression.register_symbol_table(table);
        // Create the parsing tree (compiling) and check error
        if (not parser.compile(expression_str, expression)) {
            pcout <<
                "ERROR PARSING EXPRESSION "<< expression_str <<": Compilation error..."
                  <<std::endl;
            pcout << parser.error() << std::endl;
            {getchar(); abort();}
        }
    }
};

BSchemeInfo::BSchemeInfo(const std::string &boundaryXmlPathFromRoot)
{
    // 1. Read from xml and store as vectors of strings the information
    XMLreader boundaryXml(boundaryXmlPathFromRoot);
    std::string tmp;
    boundaryXml["generic_eli"]["alphaMin"].read(tmp);
    alphaMin_expression_str = tmp;

    boundaryXml["generic_eli"]["alphaPlus"].read(tmp);
    alphaPlus_expression_str = tmp;

    boundaryXml["generic_eli"]["beta"].read(tmp);
    beta_expression_str = tmp;

    boundaryXml["generic_eli"]["kMin"].read(tmp);
    kMin_expression_str = tmp;

    boundaryXml["generic_eli"]["kPlus"].read(tmp);
    kPlus_expression_str = tmp;

    parse(symbolTable, alphaMin_expression, alphaMin_expression_str, main_parser);
    parse(symbolTable, alphaPlus_expression, alphaPlus_expression_str, main_parser);
    parse(symbolTable, beta_expression, beta_expression_str, main_parser);
    parse(symbolTable, kMin_expression, kMin_expression_str, main_parser);
    parse(symbolTable, kPlus_expression, kPlus_expression_str, main_parser);

}

enum class BCmodel {HW, BFL, FH, MLS,
    ELIULC,
    LIgeneric, ELIUL, ELIULK1, ELIULK3, ELIULK4};
/**
 * This helper functions return a voxelized domain form a TriangleSet for an
 * external flow
 * @tparam Real
 * @tparam Descriptor
 * @param parameters
 * @param triangle_set
 * @return VoxelizedDomain3D<Real>* voxalized_domain
 */
template <typename Real, template <typename U> class Descriptor>
auto voxelize_helper(const IncomprFlowParam<Real> parameters,
                     TriangleSet<Real>* triangle_set) {
    // 1. Transform the TriangleSet in two more advanced data structures to
    // handle the boundary
    //    TriangleSet->DEFscaledMesh->TriangleBoundary3D
    plint xDirection = 0;
    plint borderWidth = 2;  // Requirement: margin>=borderWidth.
    plint margin = 2;  // Extra margin of allocated cells around the obstacle,
                       // for the case of moving walls.
    plint blockSize =
        0;  // Size of blocks in the sparse/parallel representation. Zero means:
            // don't use sparse representation.
    auto defMesh = new DEFscaledMesh<Real>(*triangle_set, 0, xDirection, margin,
                                           Dot3D(0, 0, 0));
    auto boundary = new TriangleBoundary3D<Real>(*defMesh);
    // 2. Create the voxel matrix defining the "inside" and the "outside"
    plint extendedEnvelopeWidth = 2;  // Extrapolated off-lattice BCs.
    auto bounding_box =
        Box3D(0, parameters.getNx() - 1, 0, parameters.getNy() - 1, 0,
              parameters.getNz() - 1);
    auto voxalized_domain = new VoxelizedDomain3D<Real>(
        *boundary, voxelFlag::outside, bounding_box, borderWidth,
        extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxalized_domain->getVoxelMatrix()) << std::endl;

    return std::tuple{voxalized_domain, boundary};
}


enum class ProcessorLevel{
        offLattice = 1,
        rhoBarJ = 2
    };
/**
 * This functions integrates a data-processor for the boundary condition in the
 * target_lattice. It needs as parameter also the voxalized_domain to know the
 * inside and the outside.
 * @tparam Real Template type parameter for real numbers (float or double)
 * @tparam Descriptor Template parameter for the lattice topology.
 * @param target_lattice
 * @param voxalized_domain
 * @return OffLatticeBoundaryCondition3D<Real,Descriptor,VelVector>
 */
template <typename Real, template <typename U> class Descriptor>
auto inject_off_lattice_bc(
    MultiBlockLattice3D<Real, Descriptor>* target_lattice,
    VoxelizedDomain3D<Real>* voxalized_domain, BCmodel kmin_,
    BSchemeInfo& xmlParse) {
    pcout << "Generating off lattice boundary conditions." << std::endl;
    using Array3D = Array<Real, 3>;
    OffLatticeBoundaryCondition3D<Real, Descriptor, Array3D>* boundaryCondition;
    auto profiles = new BoundaryProfiles3D<Real, Array3D>;
//    bool useAllDirections = true;
    OffLatticeModel3D<Real, Array3D>* offLatticeModel = nullptr;
    profiles->setWallProfile(new NoSlipProfile3D<Real>);

    pcout << "Setting noDynamics for voxelFlag::inside and voxelFlag::innerBorder...";
    defineDynamics(*target_lattice, voxalized_domain->getVoxelMatrix(),
                   target_lattice->getBoundingBox(),
                   new NoDynamics<Real, Descriptor>(), voxelFlag::inside);
    defineDynamics(*target_lattice, voxalized_domain->getVoxelMatrix(),
                   target_lattice->getBoundingBox(),
                   new NoDynamics<Real, Descriptor>(), voxelFlag::innerBorder);
    pcout << "done." << std::endl;

    auto coefficients = [&xmlParse](Real q, Real tauPlus,
                           Real tauMinus) -> std::array<Real, 5> {
        xmlParse.q = q;
        xmlParse.up = q > 0.5 ? 1.0 : 0.0;
        xmlParse.down = 1. - xmlParse.up;
        xmlParse.tau_plus = tauPlus;
        xmlParse.tau_min = tauMinus;
        xmlParse.Lambda_min = tauMinus-0.5;
        xmlParse.Lambda_plus = tauPlus-0.5;
        Real alphaPlus = xmlParse.alphaPlus_expression.value();
        Real alphaMinus = xmlParse.alphaMin_expression.value();
        Real beta = xmlParse.beta_expression.value();
        Real Kplus = xmlParse.kPlus_expression.value();
        Real Kmin = xmlParse.kMin_expression.value();;
        return {{alphaPlus, alphaMinus, beta, Kplus, Kmin}};
    };

    switch (kmin_) {
        case BCmodel::HW:
            offLatticeModel = new HWLatticeModel3D<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::BFL:
            offLatticeModel = new BouzidiOffLatticeModel3D<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::FH:
            offLatticeModel = new FilippovaHaenelLocalModel3D<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::MLS:
            offLatticeModel = new MeiLuoShyyModel3D<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::ELIULC:
            offLatticeModel = new ELIULC<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::ELIUL:
            offLatticeModel = new ELIUL<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::ELIULK1:
            offLatticeModel = new ELIULK1<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::ELIULK3:
            offLatticeModel = new ELIULK3<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::ELIULK4:
            offLatticeModel = new ELIULK4<Real, Descriptor>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside);
             break;
        case BCmodel::LIgeneric:
            offLatticeModel = new ELIgeneric<Real, Descriptor, decltype(coefficients)>(
                new TriangleFlowShape3D<Real, Array<Real, 3> >(
                    voxalized_domain->getBoundary(), *profiles),
                voxelFlag::outside, coefficients);
            break ;
        default:
            pcout << "Case non defined...aborting..."<<std::endl;
            abort();
    }

//    FilippovaHaenelLocalModel3D<T, DESCRIPTOR> *model =
//        new FilippovaHaenelLocalModel3D<T, DESCRIPTOR>(
//            new TriangleFlowShape3D<ST, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles),
//            flowType);
    boundaryCondition =
        new OffLatticeBoundaryCondition3D<Real, Descriptor, Array3D>(
            offLatticeModel, *voxalized_domain, *target_lattice);
    std::vector<MultiBlock3D *> rhoBarJarg;
    plint numScalars = 4;
    plint extendedEnvelopeWidth = 2;
    MultiNTensorField3D<Real> *rhoBarJfield =
        generateMultiNTensorField3D<Real>(*target_lattice, extendedEnvelopeWidth, numScalars);
    rhoBarJfield->toggleInternalStatistics(false);
    rhoBarJarg.push_back(rhoBarJfield);
    boundaryCondition->insert(rhoBarJarg,(plint)ProcessorLevel::offLattice);
//    plint processorLevel = -2;
    integrateProcessingFunctional(
        new PackedRhoBarJfunctional3D<Real, Descriptor>(), target_lattice->getBoundingBox(), *target_lattice,
        *rhoBarJfield, (plint)ProcessorLevel::rhoBarJ);

//    initializeAtEquilibrium(
//            *target_lattice, target_lattice->getBoundingBox(), 1., Array<Real, 3>(0., 0., 0.));
//    target_lattice->initialize();
    applyProcessingFunctional(
            new PackedRhoBarJfunctional3D<Real, Descriptor>(), target_lattice->getBoundingBox(), *target_lattice,
            *rhoBarJfield);
    return boundaryCondition;
}

/**
 * This function integrates the on-lattice external boundary condition in the
 * lattice using the setVelocityConditionOnBlockBoundaries helper and sets the
 * initial conditions.
 * @tparam Real real numbers type
 * @tparam Descriptor lattice topology
 * @param target_lattice target lattice for the boundary conditions
 * @param parameters IncomprFlowParam<Real>
 * @return
 */
template <typename Real, template <typename U> class Descriptor>
auto inject_on_lattice_bc(MultiBlockLattice3D<Real, Descriptor>* target_lattice,
                          IncomprFlowParam<Real> const& parameters,
                          Array<Real, 3> spherePosition) {
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D outlet(nx - 1, nx, 1, ny - 2, 0, nz - 1);
    Box3D inlet(0, 0, 1, ny - 2, 1, nz - 2);
//    Box3D (0, nx - 1, 0, 0, 1, nz - 2);
    OnLatticeBoundaryCondition3D<Real, Descriptor>* onlatt_boundary_condition =
        createLocalBoundaryCondition3D<Real, Descriptor>();

    // Sets periodicity in all directions
    target_lattice->periodicity().toggleAll(true);
    target_lattice->periodicity().toggle(0,false);

    // Create Velocity boundartoggle(y conditions everywhere. Behind the scene
    // integrates a data processors in the lattice for boundary conditions on
    // surfaces, edges and corners.
    onlatt_boundary_condition->setVelocityConditionOnBlockBoundaries(
        *target_lattice, target_lattice->getBoundingBox());
//    onlatt_boundary_condition->setVelocityConditionOnBlockBoundaries(
//        *target_lattice, Box3D(0, nx - 1, 0, 0, 1, nz - 2), boundary::neumann);
//    onlatt_boundary_condition->setVelocityConditionOnBlockBoundaries(
//        *target_lattice, Box3D(0, nx - 1, ny - 1, ny - 1, 1, nz - 2), boundary::neumann);
//    onlatt_boundary_condition->setPressureConditionOnBlockBoundaries(
//        *target_lattice, outlet, boundary::density);

    // Define the value of the imposed velocity on all nodes which have previously been
    //   defined to be velocity boundary nodes.
//    setBoundaryVelocity(*target_lattice, inlet,
//                        ConstantVelocity<Real>(parameters));
    setBoundaryVelocity(*target_lattice, target_lattice->getBoundingBox().enlarge(1),
                        AdaptiveVelocityProfileOfSphereInAChannel<Real>(parameters,spherePosition));
    setBoundaryDensity(*target_lattice, outlet,
                       ConstantDensity<Real>(1.0));
//    initializeAtEquilibrium(
//        *target_lattice, target_lattice->getBoundingBox(),
//        AdaptiveVelocityProfileOfSphereInAChannel<Real>(parameters,spherePosition));

    return onlatt_boundary_condition;
}

template <typename Real, template <typename U> class Descriptor>
void createSpongeZones(IncomprFlowParam<Real> const& parameters,
                       Array<plint, 6> const& numSpongeCells,
                       MultiBlockLattice3D<Real, Descriptor>* lattice,
                       bool useSmagorinskySponges = false,
                       Real bulkValue = 0.01, Real targetValue = 0.5) {
    if (std::accumulate(&numSpongeCells[0], &numSpongeCells[5], (Real)0.0) >
        0) {
        if (useSmagorinskySponges) {
            pcout << "Generating Smagorinsky sponge zones." << std::endl;

            std::vector<MultiBlock3D*> args;
            args.push_back(lattice);
            applyProcessingFunctional(
                new SmagorinskySpongeZone3D<Real, Descriptor>(
                    parameters.getNx(), parameters.getNy(), parameters.getNz(),
                    bulkValue, targetValue, numSpongeCells),
                lattice->getBoundingBox(), args);
        } else {
            pcout << "Generating viscosity sponge zones." << std::endl;
            bulkValue = parameters.getOmega();

            std::vector<MultiBlock3D*> args;
            args.push_back(lattice);
            applyProcessingFunctional(
                new ViscositySpongeZone3D<Real, Descriptor>(
                    parameters.getNx(), parameters.getNy(), parameters.getNz(),
                    bulkValue, numSpongeCells),
                lattice->getBoundingBox(), args);
        }
    }
}

#endif  // CURVED_BOUNDARIES_SETUP_H
