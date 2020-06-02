#include <cmath>
#include <time.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <random>

#ifdef MSVC
#define NOMINMAX // problems with defined macros
#include <windows.h>
#undef interface // problems with defined macros
#else
#include <unistd.h>
#endif

// Load Palabos Lib
#include "palabos3D.h"
#include "palabos3D.hh"

// Load npFEM Lib
#include "rbcShapeOp.h"
#include "rbcGlobal.h"
#include "rbcParticle.h"
#include "rbcCommunicate.h"
#include "visualization.h"

using namespace plb;
using namespace std;
using namespace plb::global;
using namespace plb::npfem;


// If FORCED, uncomment the 2 following lines
//#define FORCED
//#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define DESCRIPTOR descriptors::D3Q19Descriptor

typedef double T;

bool CellPacking;

// So that at least all neighbors of a bulk vertex are found in the envelope.
plint particleEnvelopeWidth;

// The shape of all walls, to be used for collisions with the bodies.
RawConnectedTriangleMesh<T>* wallMesh;
std::string wallStlFileName;
std::vector<Array<T, 3>> wallVertexNormals;

// Tag to distinguish the type of body. The key is the unique ID of the body.
std::map<pluint, pluint> bodyToType;

// To be used by PIK
T collisions_threshold_rep;
T collisions_threshold_nonRep;

// Generic solvers
// ShapeOp always in physical units
ShapeOp_Solver RBC_shapeOpSolverTemplate;
ShapeOp_Solver PLT_shapeOpSolverTemplate;

// Generic meshes, used as template to construct the other ones.
// The mesh stores all the vertices (interior + onSurface), but the connectivity
// refers only to the surface elements
// In lattice units.
RawConnectedTriangleMesh<T> *rbcTemplate, *pltTemplate;

// Indicates on which processor the ShapeOp solver of each body is located.
// This variable does not change during the time iterations, and it
// is exactly the same on each processor.
// In the GPU version, there are CPUs who do not get any body because there is
// no matched GPU.
std::map<pluint, pluint> bodyToProc;

// The ShapeOp solvers for all bodies that are on the current processor.
// This variable does not change during the time iterations, and it
// is different on each processor.
// In physical units.
std::vector<ShapeOpBody<T>> shapeOpBodies;
std::vector<ShapeOp_Solver*> shapeOpRBCs, shapeOpPLTs;
std::vector<ShapeOp_Solver*> shapeOpRBCs_empty, shapeOpPLTs_empty;
#ifdef NPFEM_CUDA
Solver_GPU sGPU;
Solver_GPU sGPU_plt;
#endif

// _p  : physical system
// _lb : lattice/ non-dimensional system
// Paramenters to be read from the XML file or to be defined manually
// Palabos params
std::string Case_Study_RBC, Case_Study_PLT;
std::string Simulation_Type, pos_filename;
bool incompressibleModel;
pluint inamuro_iT;
T nu_p, nu_lb;
T dx_p, dt_p, rho_p;
T u_p, u_lb, shear_rate;
T Re, Re_rbc, Ca, Ma;
T lx_p, ly_p, lz_p;
T tau_lb;
T PI_;
T ht, RBCsToPLT;
// ShapeOp params in Physical Units
bool applyGlobalVolumeConservation_ShapeOp_RBC, applyGlobalVolumeConservation_ShapeOp_PLT;
T rho_ShapeOp_RBC, Calpha_ShapeOp_RBC, Cbeta_ShapeOp_RBC, globalVolumeConservationWeight_ShapeOp_RBC;
T rho_ShapeOp_PLT, Calpha_ShapeOp_PLT, Cbeta_ShapeOp_PLT, globalVolumeConservationWeight_ShapeOp_PLT;
T collisions_weight_rep, collisions_weight_nonRep, beta_morse;
// ShapeOp General
pluint timestep_ShapeOp;
pluint max_iterations_ShapeOp, max_line_search_loops_ShapeOp, m_ShapeOp;
pluint max_attempts_to_solve_stagnation_ShapeOp, convergence_window_ShapeOp;
T tol_ShapeOp, gamma_ShapeOp, gamma2_ShapeOp, C_EWMA_ShapeOp;

pluint iT, numOfBlocks;

// Output & Stopping Criterion
std::string OutputDir;
pluint output_timestep;
T max_time;
T dump_RBCs_ratio_VTKs, dump_PLTs_ratio_VTKs;
T dump_RBCs_ratio_COMs, dump_PLTs_ratio_COMs;
bool withBodies;
plint fileNamePadding = 12;

// CheckPointing
bool restartFromCP, continuePLB;
std::string CP_OutputDir;
pluint checkpoint_timestep;

// Initial placing in physical units. Rotations are in degrees.
std::vector<Array<T, 3>> iniPositions, iniRotations;

// Range refering to myRank
std::pair<plint, plint> myRangeRBC;
std::pair<plint, plint> myRangePLT;

// Total number of bodies (not per task but global)
pluint numRBC, numPLT;
int loc_rank = 0, sock_rank, local_body_id = 0, total_devices_nb = 0, number_devices = 1, number_nodes = 1;
// Poiseuille Flow
T drivingForce_lb;
T pipeRadius_LB;

// Dumping of Centers of Mass for every body
MatrixXX centers_log;

// for debuging only
void flatten_mesh(RawConnectedTriangleMesh<double> *obstacle, double *obstacle_vertices, double *obstacle_normals) {

	std::vector<Array<T, 3>> local_vertices = obstacle->getVertices();
	int nb_obstacle_vertices = obstacle->getNumVertices();

	for (int i = 0; i < nb_obstacle_vertices * 3; i += 3) {
		//auto vertex = wallMesh->vertex(i/3);
		//printf("%f %f %f \n", (*vertex)[0], (*vertex)[1], (*vertex)[2]);
		obstacle_vertices[i] = local_vertices[i][0] * dx_p;
		obstacle_vertices[i + 1] = local_vertices[i][1] * dx_p;
		obstacle_vertices[i + 2] = local_vertices[i][2] * dx_p;

		auto n = obstacle->vertex(i / 3)->normal();
		obstacle_normals[i] = n[0];
		obstacle_normals[i + 1] = n[1];
		obstacle_normals[i + 2] = n[2];
	}
}

// LB units
template <typename T> class Pipe : public DomainFunctional3D {
public:
	Pipe(Array<T, 3> center_, T radius_)
		: center(center_)
		, radius(radius_)
	{ }

	virtual bool operator()(plint iX, plint iY, plint iZ) const
	{
		// True if solid
		return (std::sqrt(util::sqr((T)iX - center[0]) + util::sqr((T)iY - center[1])) >= radius);
	}

	virtual Pipe<T>* clone() const { return new Pipe<T>(*this); }
private:
	Array<T, 3> center;
	T radius;
};

// Shear Flow
void flowSetup(Group3D& group, OnLatticeBoundaryCondition3D<T, DESCRIPTOR>& boundaryCondition, T vel)
{
	// number of lattice cells in each direction
	const pluint nx = group.getBoundingBox().getNx();
	const pluint ny = group.getBoundingBox().getNy();
	const pluint nz = group.getBoundingBox().getNz();

	Box3D top = Box3D(0, nx - 1, ny - 1, ny - 1, 0, nz - 1);
	Box3D bottom = Box3D(0, nx - 1, 0, 0, 0, nz - 1);

	group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggleAll(true);
	group.getScalar<T>("rhoBar").periodicity().toggleAll(true);
	group.getTensor<T, 3>("j").periodicity().toggleAll(true);
	group.getTensor<T, 6>("piNeq").periodicity().toggleAll(true);

	// 1 corresponds to the y-direction (walls) (0->x, 2->z)
	group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(1, false);
	group.getScalar<T>("rhoBar").periodicity().toggle(1, false);
	group.getTensor<T, 3>("j").periodicity().toggle(1, false);
	group.getTensor<T, 6>("piNeq").periodicity().toggle(1, false);

    group.getDenseParticles<T, DESCRIPTOR>("vertexParticles").periodicity().toggleAll(true);

	if (Simulation_Type.compare("Shear_Flow") == 0)
    {
		boundaryCondition.addVelocityBoundary1P(
			top, group.getLattice<T, DESCRIPTOR>("lattice"));
		boundaryCondition.addVelocityBoundary1N(
			bottom, group.getLattice<T, DESCRIPTOR>("lattice"));

		// Velocity in lattice units
		setBoundaryVelocity(group.getLattice<T, DESCRIPTOR>("lattice"), top, Array<T, 3>(0., 0., vel));
		setBoundaryVelocity(group.getLattice<T, DESCRIPTOR>("lattice"), bottom, Array<T, 3>(0., 0., 0.));
		initializeAtEquilibrium(group.getLattice<T, DESCRIPTOR>("lattice"), group.getBoundingBox(), 1.0, Array<T, 3>::zero());
		initializeAtEquilibrium(group.getLattice<T, DESCRIPTOR>("lattice"), top, 1.0, Array<T, 3>(0., 0., vel));
		initializeAtEquilibrium(group.getLattice<T, DESCRIPTOR>("lattice"), bottom, 1.0, Array<T, 3>(0., 0., 0.));
	}
	else
    {
		// rho_LB = 1.
		initializeAtEquilibrium(group.getLattice<T, DESCRIPTOR>("lattice"), group.getBoundingBox(), 1., Array<T, 3>::zero());
	}

    group.getLattice<T, DESCRIPTOR>("lattice").initialize();
}

// Poiseuille Flow
void flowSetup(Group3D& group)
{
	// POISEUILLE FLOW
	group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggleAll(false);
	group.getScalar<T>("rhoBar").periodicity().toggleAll(false);
	group.getTensor<T, 3>("j").periodicity().toggleAll(false);
	group.getTensor<T, 6>("piNeq").periodicity().toggleAll(false);

	// 2 corresponds to the z-direction (flow direction)
	group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(2, true);
	group.getScalar<T>("rhoBar").periodicity().toggle(2, true);
	group.getTensor<T, 3>("j").periodicity().toggle(2, true);
	group.getTensor<T, 6>("piNeq").periodicity().toggle(2, true);

    group.getDenseParticles<T, DESCRIPTOR>("vertexParticles").periodicity().toggleAll(true);

	// rho_LB = 1.
	initializeAtEquilibrium(group.getLattice<T, DESCRIPTOR>("lattice"), group.getBoundingBox(), 1., Array<T, 3>::zero());
    group.getLattice<T, DESCRIPTOR>("lattice").initialize();
}

template <class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice, T dx, T dt, pluint iter)
{
#ifdef MSVC
	VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, fileNamePadding), dx);
#else
	ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, fileNamePadding), 1, dx);
#endif

	// From lattice to physical
	vtkOut.writeData<3, float>(*computeVelocity(lattice), "Velocity (micro-m/micro-s)", dx / dt);
	//vtkOut.writeData<float>(*computeDensity(lattice), "Density (lb units)", 1.0);
}

// Assign the ShapeOp bodies to their respective processors, and construct them.
// ShapeOp objects always in physical units, while palabos always in lattice
// units.
void iniShapeOpSolvers(std::string fname_points, std::string fname_constraints,
	std::string fname_connectivity, std::string fname_onsurface)
{
    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: Start");

	// Create the Palabos-Mesh template. Mesh in lattice units
	rbcTemplate = new RawConnectedTriangleMesh<T>(generateMesh<T>(RBC_shapeOpSolverTemplate, dx_p));
	pltTemplate = new RawConnectedTriangleMesh<T>(generateMesh<T>(PLT_shapeOpSolverTemplate, dx_p));

	if (wallStlFileName == "")
	{
		wallMesh = nullptr;
	}
	else
	{
        pcout << "I am reading the wall: " << wallStlFileName << std::endl;
		wallMesh = new RawConnectedTriangleMesh<T>(stlToConnectedTriangleMesh<T>(wallStlFileName, getEpsilon<T>(DBL)));
		wallMesh->scale(1. / dx_p); // Convert from physical to lattice units.
	}

	// Distribute the ShapeOp bodies to the processors as evenly as possible.
	pluint numProc = global::mpi().getSize();

	// First we stack all the RBCs, then the PLTs, ... (in the .pos file)
	// So the unique body IDs have a pattern.

	// Every available GPU takes a range of bodies to solve
    // In the GPU version, CPUs that are not linked to GPUs get zero bodies to solve
	std::vector<std::pair<plint, plint>> ranges_RBCs;
	std::vector<std::pair<plint, plint>> ranges_PLTs;
	// Range refering to myRank
	myRangeRBC = std::pair<plint, plint>(0, -1);
	myRangePLT = std::pair<plint, plint>(0, -1);

    if (!withBodies)
        return;

    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: Ranges start");
#ifdef NPFEM_CUDA
    //////////////////////////////////////////////////////////////////////////////////////////
    // This code works only for cases where every node has at least one GPU attached to it. //
    //////////////////////////////////////////////////////////////////////////////////////////

    // number_devices refers to per node
	total_devices_nb = number_devices*number_nodes;

	// Find local rank, meaning rank local to the cluster node.
	MPI_Comm shmcomm;
	// MPI_COMM_TYPE_SHARED: the communicator will be split into subcommunicators whose processes can create shared memory regions
	// This is highly useful when you want to have a hierarchy of communication where communicators are created based on node-level
	// memory locality, even if a shared memory region is not needed
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
	MPI_Comm_rank(shmcomm, &loc_rank); // "nodal" rank

	int *partition = new int[numProc];
	int *who_got_gpu = new int[total_devices_nb];

	MPI_Gather(&loc_rank, 1, MPI_INT, partition, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (global::mpi().getRank() == 0)
	{
		int gpu = 0;
		for (int i = 0; i < numProc; i++)
		{
            // number_devices refers to per node
			// use the nodal rank to match a cpu with a corresponding node-gpu
			if (partition[i] < number_devices)
			{
				who_got_gpu[gpu++] = i; // store the global rank that gets a gpu
			}
		}
	}

	MPI_Bcast(who_got_gpu, total_devices_nb, MPI_INT, 0, MPI_COMM_WORLD);

	// Distribute the bodies to the available devices
	if (numRBC > 0)
		util::linearRepartition(0, numRBC - 1, std::min((pluint)total_devices_nb, numRBC), ranges_RBCs);

	if (numPLT > 0)
		util::linearRepartition(numRBC, numPLT + numRBC - 1, std::min((pluint)total_devices_nb, numPLT), ranges_PLTs);

	// Traverse the available GPUs and distribute ranges of bodies
	// RBCs
	for (plint iProc = 0; iProc < (plint)ranges_RBCs.size(); ++iProc)
	{
		for (plint bodyID = ranges_RBCs[iProc].first; bodyID <= ranges_RBCs[iProc].second; ++bodyID)
			bodyToProc[bodyID] = who_got_gpu[iProc];

        if (who_got_gpu[iProc] == global::mpi().getRank())
            myRangeRBC = ranges_RBCs[iProc];
	}
	// PLTs
	for (plint iProc = 0; iProc < (plint)ranges_PLTs.size(); ++iProc)
	{
		for (plint bodyID = ranges_PLTs[iProc].first; bodyID <= ranges_PLTs[iProc].second; ++bodyID)
			bodyToProc[bodyID] = who_got_gpu[iProc];

        if (who_got_gpu[iProc] == global::mpi().getRank())
            myRangePLT = ranges_PLTs[iProc];
	}
#else
	if (numRBC > 0)
		util::linearRepartition(0, numRBC - 1, std::min((pluint)numProc, numRBC), ranges_RBCs);

	if (numPLT > 0)
		util::linearRepartition(numRBC, numPLT + numRBC - 1, std::min((pluint)numProc, numPLT), ranges_PLTs);

	for (plint iProc = 0; iProc < (plint)ranges_RBCs.size(); ++iProc)
	{
		for (plint bodyID = ranges_RBCs[iProc].first; bodyID <= ranges_RBCs[iProc].second; ++bodyID)
			bodyToProc[bodyID] = iProc;
        
        if ((int)iProc == global::mpi().getRank())
            myRangeRBC = ranges_RBCs[iProc];
	}

	for (plint iProc = 0; iProc < (plint)ranges_PLTs.size(); ++iProc)
	{
		for (plint bodyID = ranges_PLTs[iProc].first; bodyID <= ranges_PLTs[iProc].second; ++bodyID)
			bodyToProc[bodyID] = iProc;

        if ((int)iProc == global::mpi().getRank())
            myRangePLT = ranges_PLTs[iProc];
	}

	loc_rank = global::mpi().getRank();
#endif
    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: Ranges end");

	// std::pair<plint, plint> myRangeRBC(0, -1);
	// std::pair<plint, plint> myRangePLT(0, -1);
	// for CPUs with no matched GPUs: { myRangeRBC.second - myRangeRBC.first + myRangePLT.second - myRangePLT.first + 2 } == 0
	shapeOpBodies = std::vector<ShapeOpBody<T>>(myRangeRBC.second - myRangeRBC.first + myRangePLT.second - myRangePLT.first + 2);
	std::vector<int> *rbc_graph = RBC_shapeOpSolverTemplate.get_mesh_graph();
	std::vector<int> *plt_graph = PLT_shapeOpSolverTemplate.get_mesh_graph();

	auto initBody = [&](std::vector<ShapeOp_Solver*> &shapeOpP, std::string &case_name, plb::RawConnectedTriangleMesh<T> *mesh_template,
		T Calpha_ShapeOp_, T  Cbeta_ShapeOp_, T rho_, T applyGlobalVolumeConservation_ShapeOp_, T globalVolumeConservationWeight_ShapeOp_,
		std::vector<int> *mesh_graph, int  bodyID, int i, int body_type)
	{
		ShapeOp_Solver *solver = &(shapeOpBodies[i].getSolver());

		solver->bodyID_ = bodyID;
		solver->bodyType_ = body_type;

		setPointsFromCSV(*solver, "./Case_Studies/" + case_name + "/" + fname_points);

		if (!restartFromCP)
		{
			Array<T, 3> position = iniPositions[bodyID];
			Array<T, 3> rotation = iniRotations[bodyID];

			solver->rotatePoints("x", rotation[0] * std::acos((T)-1) / 180);
			solver->rotatePoints("y", rotation[1] * std::acos((T)-1) / 180);
			solver->rotatePoints("z", rotation[2] * std::acos((T)-1) / 180);

			solver->shiftPoints(Vector3(position[0], position[1], position[2]));

			setConstraintsFromCSV(*solver, "./Case_Studies/" + case_name + "/" + fname_constraints);
			setConnectivityListFromCSV(*solver, "./Case_Studies/" + case_name + "/" + fname_connectivity);
			setOnSurfaceParticle(*solver, "./Case_Studies/" + case_name + "/" + fname_onsurface);

			Matrix3X forces(3, solver->getPoints().cols());
			forces.setZero();
			addVertexForce(*solver, forces);

			shapeOpBodies[i].init(bodyID, mesh_template, mesh_graph);
			solver->initialize(Calpha_ShapeOp_, Cbeta_ShapeOp_, dt_p * (T)(timestep_ShapeOp), rho_, false, applyGlobalVolumeConservation_ShapeOp_, globalVolumeConservationWeight_ShapeOp_, dx_p, dt_p);
		}
		else
		{
			setConstraintsFromCSV(*solver, "./Case_Studies/" + case_name + "/" + fname_constraints);
			setConnectivityListFromCSV(*solver, "./Case_Studies/" + case_name + "/" + fname_connectivity);
			setOnSurfaceParticle(*solver, "./Case_Studies/" + case_name + "/" + fname_onsurface);

			Matrix3X forces(3, solver->getPoints().cols());
			forces.setZero();
			addVertexForce(*solver, forces);

			shapeOpBodies[i].init(bodyID, mesh_template, mesh_graph);
			solver->initialize(Calpha_ShapeOp_, Cbeta_ShapeOp_, dt_p * (T)(timestep_ShapeOp), rho_, false, applyGlobalVolumeConservation_ShapeOp_, globalVolumeConservationWeight_ShapeOp_, dx_p, dt_p);

			if (!continuePLB)
			{
                //setPointsFromCSV(*solver, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv");
                // This is to avoid setting the initial conditions with deformed bodies such that we avoid energy bursts
				setPointsFromCSV(*solver, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv", true); 
			}
			else
			{
				setPointsFromCSV(*solver, createFileName(CP_OutputDir + "body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
				setVelsFromCSV(*solver, createFileName(CP_OutputDir + "vels_body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
			}
		}
		shapeOpP.push_back(solver);
	};

    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: npFEM init start");
	// init bodies, if any
	if (shapeOpBodies.size())
    {
		local_body_id = 0;
		for (plint bodyID = myRangeRBC.first; bodyID <= myRangeRBC.second; ++bodyID)
        {
			initBody(shapeOpRBCs, Case_Study_RBC, rbcTemplate, Calpha_ShapeOp_RBC, Cbeta_ShapeOp_RBC, rho_ShapeOp_RBC, applyGlobalVolumeConservation_ShapeOp_RBC, globalVolumeConservationWeight_ShapeOp_RBC, rbc_graph, bodyID, local_body_id, 0);
			local_body_id++;
		}

		for (plint bodyID = myRangePLT.first; bodyID <= myRangePLT.second; ++bodyID)
        {
			initBody(shapeOpPLTs, Case_Study_PLT, pltTemplate, Calpha_ShapeOp_PLT, Cbeta_ShapeOp_PLT, rho_ShapeOp_PLT, applyGlobalVolumeConservation_ShapeOp_PLT, globalVolumeConservationWeight_ShapeOp_PLT, plt_graph, bodyID, local_body_id, 1);
			local_body_id++;
		}
	}
    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: npFEM init end");

    plb::global::logfile("progress.log").flushEntry("initializing npFEM solvers: End");
}


template <typename T, template <typename U> class Descriptor>
void iniPalabosParticles(Group3D& group, T vel)
{
    std::vector<Particle3D<T, Descriptor>*> particles;

    std::vector<pluint> verIDs(0), bdIDs(0);
    std::vector<Array<T, 3>> ps(0), vs(0);


    plb::global::logfile("progress.log").flushEntry("Injecting Particles: RBCs");
    const std::vector<bool>& RBC_onSurfaceParticle = RBC_shapeOpSolverTemplate.get_onSurfaceParticle();
    for (plint bodyID = myRangeRBC.first; bodyID <= myRangeRBC.second; ++bodyID)
    {
        ShapeOp_Solver shapeOpSolverTemplate_tmp = RBC_shapeOpSolverTemplate;

        if (!restartFromCP)
        {
            Array<T, 3> position = iniPositions[bodyID];
            Array<T, 3> rotation = iniRotations[bodyID];

            shapeOpSolverTemplate_tmp.rotatePoints("x", rotation[0] * std::acos((T)-1) / 180);
            shapeOpSolverTemplate_tmp.rotatePoints("y", rotation[1] * std::acos((T)-1) / 180);
            shapeOpSolverTemplate_tmp.rotatePoints("z", rotation[2] * std::acos((T)-1) / 180);

            shapeOpSolverTemplate_tmp.shiftPoints(Vector3(position[0], position[1], position[2]));
        }
        else
        {
            if (!continuePLB)
            {
                //setPointsFromCSV(shapeOpSolverTemplate_tmp, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv");
                // This is to avoid setting the initial conditions with deformed bodies such that we avoid energy bursts
                setPointsFromCSV(shapeOpSolverTemplate_tmp, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv", true);
            }
            else
            {
                setPointsFromCSV(shapeOpSolverTemplate_tmp, createFileName(CP_OutputDir + "body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
                setVelsFromCSV(shapeOpSolverTemplate_tmp, createFileName(CP_OutputDir + "vels_body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
            }
        }

        const Matrix3X& tmp_points = shapeOpSolverTemplate_tmp.getPoints();
        const Matrix3X& tmp_vels = shapeOpSolverTemplate_tmp.getVelocities();

        for (pluint vertexID = 0; vertexID < (pluint)tmp_points.cols(); ++vertexID)
        {
            // Palabos cares only for particles that are on the RBC surface.
            // The interiors are for the ShapeOp only (if any).
            if (!RBC_onSurfaceParticle[vertexID])
                continue;

            // From physical to lattice
            Array<T, 3> position;
            position[0] = tmp_points.col(vertexID)[0];
            position[1] = tmp_points.col(vertexID)[1];
            position[2] = tmp_points.col(vertexID)[2];
            position /= dx_p;
            Array<T, 3> velocity;
            velocity[0] = tmp_vels.col(vertexID)[0];
            velocity[1] = tmp_vels.col(vertexID)[1];
            velocity[2] = tmp_vels.col(vertexID)[2];
            velocity /= (dx_p / dt_p);

            // Check that the particles are inside the domain, else translate
            // them appropriately (lb units)
            bool x = false, y = false, z = false;
            controlPositionPeriodic<T>(position,
                group.getLattice<T, DESCRIPTOR>("lattice").getNx() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNy() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNz() - 1,
                x, y, z);

            if (global::mpi().getRank() == 0)
            {
                // particle position in lattice units
                particles.push_back(new bodyVertexParticle<T, Descriptor>(vertexID, position, bodyID, velocity));
            }
            else
            {
                verIDs.push_back(vertexID);
                ps.push_back(position);
                bdIDs.push_back(bodyID);
                vs.push_back(velocity);
            }
        }
    }


    plb::global::logfile("progress.log").flushEntry("Injecting Particles: PLTs");
    const std::vector<bool>& PLT_onSurfaceParticle = PLT_shapeOpSolverTemplate.get_onSurfaceParticle();
    for (plint bodyID = myRangePLT.first; bodyID <= myRangePLT.second; ++bodyID)
    {
        ShapeOp_Solver shapeOpSolverTemplate_tmp = PLT_shapeOpSolverTemplate;

        if (!restartFromCP)
        {
            Array<T, 3> position = iniPositions[bodyID];
            Array<T, 3> rotation = iniRotations[bodyID];

            shapeOpSolverTemplate_tmp.rotatePoints("x", rotation[0] * std::acos((T)-1) / 180);
            shapeOpSolverTemplate_tmp.rotatePoints("y", rotation[1] * std::acos((T)-1) / 180);
            shapeOpSolverTemplate_tmp.rotatePoints("z", rotation[2] * std::acos((T)-1) / 180);

            shapeOpSolverTemplate_tmp.shiftPoints(Vector3(position[0], position[1], position[2]));
        }
        else
        {
            if (!continuePLB)
            {
                //setPointsFromCSV(shapeOpSolverTemplate_tmp, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv");
                // This is to avoid setting the initial conditions with deformed bodies such that we avoid energy bursts
                setPointsFromCSV(shapeOpSolverTemplate_tmp, CP_OutputDir + "body_ID_" + util::val2str(bodyID) + ".csv", true);
            }
            else
            {
                setPointsFromCSV(shapeOpSolverTemplate_tmp, createFileName(CP_OutputDir + "body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
                setVelsFromCSV(shapeOpSolverTemplate_tmp, createFileName(CP_OutputDir + "vels_body_ID_" + util::val2str(bodyID) + "_iT_", iT, fileNamePadding) + ".csv");
            }
        }

        const Matrix3X& tmp_points = shapeOpSolverTemplate_tmp.getPoints();
        const Matrix3X& tmp_vels = shapeOpSolverTemplate_tmp.getVelocities();

        for (pluint vertexID = 0; vertexID < (pluint)tmp_points.cols(); ++vertexID)
        {
            // Palabos cares only for particles that are on the RBC surface.
            // The interiors are for the ShapeOp only (if any).
            if (!PLT_onSurfaceParticle[vertexID])
                continue;

            // From physical to lattice
            Array<T, 3> position;
            position[0] = tmp_points.col(vertexID)[0];
            position[1] = tmp_points.col(vertexID)[1];
            position[2] = tmp_points.col(vertexID)[2];
            position /= dx_p;
            Array<T, 3> velocity;
            velocity[0] = tmp_vels.col(vertexID)[0];
            velocity[1] = tmp_vels.col(vertexID)[1];
            velocity[2] = tmp_vels.col(vertexID)[2];
            velocity /= (dx_p / dt_p);

            // Check that the particles are inside the domain, else translate
            // them appropriately (lb units)
            bool x = false, y = false, z = false;
            controlPositionPeriodic<T>(position,
                group.getLattice<T, DESCRIPTOR>("lattice").getNx() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNy() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNz() - 1,
                x, y, z);

            if (global::mpi().getRank() == 0)
            {
                // particle position in lattice units
                particles.push_back(new bodyVertexParticle<T, Descriptor>(vertexID, position, bodyID, velocity));
            }
            else
            {
                verIDs.push_back(vertexID);
                ps.push_back(position);
                bdIDs.push_back(bodyID);
                vs.push_back(velocity);
            }
        }
    }


    plb::global::logfile("progress.log").flushEntry("Injecting Particles: Walls");
    if (wallMesh)
    {
        for (pluint vertexID = 0; vertexID < (pluint)wallMesh->getNumVertices(); ++vertexID)
        {
            auto vertex = wallMesh->vertex(vertexID);
            Array<T, 3> position((*vertex)[0], (*vertex)[1], (*vertex)[2]);

            // Check that the particles are inside the domain, else translate
            // them appropriately (lb units)
            bool x = false, y = false, z = false;
            controlPositionPeriodic<T>(position,
                group.getLattice<T, DESCRIPTOR>("lattice").getNx() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNy() - 1,
                group.getLattice<T, DESCRIPTOR>("lattice").getNz() - 1,
                x, y, z);

            if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
            {
                // Zero velocity on these boundaries
                particles.push_back(new bodyVertexParticle<T, Descriptor>(vertexID, position, IDforWall, Array<T, 3>::zero()));
            }
            else
            {
                if (position[1] < ((ly_p / 2.) / dx_p)) // bottom
                    particles.push_back(new bodyVertexParticle<T, Descriptor>(vertexID, position, IDforWall, Array<T, 3>(0., 0., 0.)));
                else // top
                    particles.push_back(new bodyVertexParticle<T, Descriptor>(vertexID, position, IDforWall, Array<T, 3>(0., 0., vel)));
            }

            const bool areaWeighted = true;
            wallVertexNormals.push_back(vertex->normal(areaWeighted));
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    // Communication Start: AtMainProc
    ///////////////////////////////////////////////////////////////////////////

    plb::global::logfile("progress.log").flushEntry("Injecting Particles: Communication");

    pluint containers_size;

    std::vector<MPI_Request*> mpiRequests;

    if (global::mpi().getRank() != 0) // SEND
    {
        containers_size = verIDs.size();

        mpiRequests.push_back(new MPI_Request);
        MPI_Isend(&containers_size, 1, MPI_UNSIGNED_LONG_LONG, 0,
            MPItagIniPalabosParticles, MPI_COMM_WORLD, mpiRequests.back());

        if (containers_size > 0)
        {
            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&verIDs[0], (int)containers_size, MPI_UNSIGNED_LONG_LONG, 0,
                MPItagIniPalabosParticles, MPI_COMM_WORLD, mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&ps[0][0], (int)(3 * containers_size), MPI_DOUBLE, 0,
                MPItagIniPalabosParticles, MPI_COMM_WORLD, mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&bdIDs[0], (int)containers_size, MPI_UNSIGNED_LONG_LONG, 0,
                MPItagIniPalabosParticles, MPI_COMM_WORLD, mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&vs[0][0], (int)(3 * containers_size), MPI_DOUBLE, 0,
                MPItagIniPalabosParticles, MPI_COMM_WORLD, mpiRequests.back());
        }
    }
    else // RECV
    {
        pluint numProc = global::mpi().getSize();
        pluint received = 0;

        while (received != (numProc - 1))
        {
            MPI_Status status, dummyStatus;

            MPI_Recv(&containers_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE,
                MPItagIniPalabosParticles, MPI_COMM_WORLD, &status);

            int sending_processor = status.MPI_SOURCE;

            if (containers_size > 0)
            {
                std::vector<pluint> verIDs_(containers_size), bdIDs_(containers_size);
                std::vector<Array<T, 3>> ps_(containers_size), vs_(containers_size);

                MPI_Recv(&verIDs_[0], (int)containers_size, MPI_UNSIGNED_LONG_LONG, sending_processor,
                    MPItagIniPalabosParticles, MPI_COMM_WORLD, &dummyStatus);

                MPI_Recv(&ps_[0][0], (int)(3 * containers_size), MPI_DOUBLE, sending_processor,
                    MPItagIniPalabosParticles, MPI_COMM_WORLD, &dummyStatus);

                MPI_Recv(&bdIDs_[0], (int)containers_size, MPI_UNSIGNED_LONG_LONG, sending_processor,
                    MPItagIniPalabosParticles, MPI_COMM_WORLD, &dummyStatus);

                MPI_Recv(&vs_[0][0], (int)(3 * containers_size), MPI_DOUBLE, sending_processor,
                    MPItagIniPalabosParticles, MPI_COMM_WORLD, &dummyStatus);

                for (pluint l = 0; l < containers_size; ++l)
                {
                    particles.push_back(
                        new bodyVertexParticle<T, Descriptor>(verIDs_[l], ps_[l], bdIDs_[l], vs_[l])
                    );
                }
            }

            ++received;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Wait for communication to be completed
    ///////////////////////////////////////////////////////////////////////////

    MPI_Status status;
    for (pluint i = 0; i < (pluint)mpiRequests.size(); ++i)
    {
        MPI_Wait(mpiRequests[i], &status);
        delete mpiRequests[i];
    }
    mpiRequests.clear();

    ///////////////////////////////////////////////////////////////////////////
    // Communication End: AtMainProc
    ///////////////////////////////////////////////////////////////////////////

    if (global::mpi().getRank() != 0)
    {
        for (pluint i = 0; i < particles.size(); ++i)
            delete particles[i];
        particles.clear();
    }

    plb::global::logfile("progress.log").flushEntry("Injecting Particles: Palabos");
    injectParticlesAtMainProc<T, Descriptor>(
        particles,
        group.getDenseParticles<T, Descriptor>("vertexParticles"),
        group.getBoundingBox()
    );

    if (global::mpi().getRank() == 0)
    {
        for (pluint i = 0; i < particles.size(); ++i)
            delete particles[i];
        particles.clear();
    }
}

// Reset communication data at every iteration
void resetCommunication()
{
	// Get rid of the local meshes of the previous iteration.
	std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().begin();
	for (; it != LocalMeshes<T>().end(); ++it)
		delete it->second;
	LocalMeshes<T>().clear();

	// Reset data structures for MPI communication
	for (pluint i = 0; i < shapeOpBodies.size(); ++i)
		shapeOpBodies[i].resetBeforeReceiving();
}

void forceMpiCompletion()
{
	typename std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().begin();
	for (; it != LocalMeshes<T>().end(); ++it)
		it->second->completeRequests();
}

void velocityMpiCompletion()
{
	for (pluint i = 0; i < shapeOpBodies.size(); ++i)
		shapeOpBodies[i].completeRequests();
}

void time_statistics(pluint iter, pluint delta_iter)
{
	T t_all = timer("all").getTime();
	T t_reset_com = timer("reset-communication").getTime();
	T t_mpi = timer("mpi").getTime();
	T t_shapeop = timer("shapeop").getTime();
	T t_apply_cont = timer("apply-containers").getTime();

	T t_actions1  = timer("actions1").getTime();
	T t_actions2a = timer("actions2a").getTime();
    T t_actions2b = timer("actions2b").getTime();
	T t_actions3  = timer("actions3").getTime();
	T t_actions4  = timer("actions4").getTime();
	T t_actions5  = timer("actions5").getTime();
	T t_actions6  = timer("actions6").getTime();

    T t_actionsForcingTerm = 0.;
    if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
    {
        t_actionsForcingTerm = timer("actionsForcingTerm").getTime();
    }

    T t_logging = timer("logging").getTime();
	T t_output  = timer("output").getTime();


	std::stringstream ss;
	ss << "At iteration " << std::setw(6) << std::left << iter
		<< "the statistitics over the past " << std::setw(6) << std::left
		<< delta_iter << "iterations is:\n";
	ss << "===================================================================="
		"==\n";
	ss << "Total elapsed time: " << setprecision(3) << (t_all / delta_iter) * 1.e3
		<< " ms per iteration.\n";

	ss << std::setw(44) << std::left << "reset-communication:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_reset_com / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_reset_com / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "ShapeOp-Palabos communication:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_mpi / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_mpi / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "ShapeOp:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_shapeop / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_shapeop / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "apply-containers:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_apply_cont / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_apply_cont / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions1 (BoxRhoBarJPiNeq):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions1 / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions1 / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions2a (LocMesh):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions2a / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions2a / delta_iter) * 1.e3 << " ms/iT)\n";

    ss << std::setw(44) << std::left << "actions2b (Combo):";
    ss << std::setw(6) << std::right << std::setprecision(3)
        << (t_actions2b / t_all) * 100.0 << "%";
    ss << std::setw(8) << std::right << std::setprecision(3) << "("
        << (t_actions2b / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions3 (IBM):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions3 / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions3 / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions4 (Collide and Stream):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions4 / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions4 / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions5 (LocMeshToPartVel):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions5 / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions5 / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "actions6 (AdvanceParts):";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_actions6 / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_actions6 / delta_iter) * 1.e3 << " ms/iT)\n";

    if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
    {
        ss << std::setw(44) << std::left << "action ForcingTerm:";
        ss << std::setw(6) << std::right << std::setprecision(3)
            << (t_actionsForcingTerm / t_all) * 100.0 << "%";
        ss << std::setw(8) << std::right << std::setprecision(3) << "("
            << (t_actionsForcingTerm / delta_iter) * 1.e3 << " ms/iT)\n";
    }


    ss << std::setw(44) << std::left << "Logging:";
    ss << std::setw(6) << std::right << std::setprecision(3)
        << (t_logging / t_all) * 100.0 << "%";
    ss << std::setw(8) << std::right << std::setprecision(3) << "("
        << (t_logging / delta_iter) * 1.e3 << " ms/iT)\n";

	ss << std::setw(44) << std::left << "Output:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< (t_output / t_all) * 100.0 << "%";
	ss << std::setw(8) << std::right << std::setprecision(3) << "("
		<< (t_output / delta_iter) * 1.e3 << " ms/iT)\n";


	ss << std::setw(44) << std::left << "SUM:";
	ss << std::setw(6) << std::right << std::setprecision(3)
		<< ((t_reset_com + t_mpi + t_shapeop + t_apply_cont +
			 t_actions1 + t_actions2a + t_actions2b + t_actions3 + t_actions4 + t_actions5 + t_actions6 + t_actionsForcingTerm +
			 t_output + t_logging) / t_all) * 100.0 << "%";
	ss << "\n\n";


	timer("all").restart();

	timer("reset-communication").reset();
	timer("mpi").reset();
	timer("shapeop").reset();
	timer("apply-containers").reset();

	timer("actions1").reset();
	timer("actions2a").reset();
    timer("actions2b").reset();
	timer("actions3").reset();
	timer("actions4").reset();
	timer("actions5").reset();
	timer("actions6").reset();
    if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
    {
        timer("actionsForcingTerm").reset();
    }

    timer("logging").reset();
	timer("output").reset();


	plb::global::logfile_nonparallel("profiling.log").flushEntry(ss.str());
}

void saveAllBlocks(Group3D& group, plint iteration)
{
	std::vector<MultiBlock3D*> checkpointBlocks;
	checkpointBlocks.push_back(&group.get("lattice"));
	checkpointBlocks.push_back(&group.get("rhoBar"));
	checkpointBlocks.push_back(&group.get("j"));
	checkpointBlocks.push_back(&group.get("piNeq"));
	//checkpointBlocks.push_back(&group.get("vertexParticles"));

	numOfBlocks = (pluint)checkpointBlocks.size();

	bool saveDynamicContent = true;
	std::string continueFile = "continue.xml";
	std::string checkpointFile = "checkpoint_";

	saveState(checkpointBlocks, iteration, saveDynamicContent, continueFile, checkpointFile, fileNamePadding);
}

void loadAllBlocks(Group3D& group, pluint& iteration)
{
	std::vector<MultiBlock3D*> checkpointBlocks;
	checkpointBlocks.push_back(&group.get("lattice"));
	checkpointBlocks.push_back(&group.get("rhoBar"));
	checkpointBlocks.push_back(&group.get("j"));
	checkpointBlocks.push_back(&group.get("piNeq"));
	//checkpointBlocks.push_back(&group.get("vertexParticles"));

	bool saveDynamicContent = true;
	std::string continueFile = "continue.xml";

	plint intIteration;
	loadState(checkpointBlocks, intIteration, saveDynamicContent, continueFile);
	iteration = intIteration;
}

#ifdef NPFEM_CUDA
void copy_data_to_gpu(Solver_GPU &sGPU, vector<Solver*> &cpu_solvers)
{
	// copy forces and colliding data from CPU solver to GPU solver
	int total_neighbor = 0;
	std::vector<int> neighbor(cpu_solvers.size() + 1);
	neighbor[0] = 0;
	for (pluint i = 0; i < cpu_solvers.size(); ++i)
	{
		neighbor[i + 1] = neighbor[i];
		for (auto &pieces : cpu_solvers[i]->verticesOfCollidingPieces_)
		{
			neighbor[i + 1] += pieces.cols();
		}
		sGPU.set_Palabos_Forces(cpu_solvers[i]->get_Palabos_Forces(), i);
	}
	total_neighbor = neighbor[cpu_solvers.size()];

	Matrix3X neighbor_vertices = Matrix3X::Zero(3, total_neighbor);
	Matrix3X neighbor_normals = Matrix3X::Zero(3, total_neighbor);
	int k = 0;
	for (ShapeOp_Solver *s : cpu_solvers)
	{
		for (int i = 0; i < s->verticesOfCollidingPieces_.size(); i++)
		{
			neighbor_vertices.block(0, k, 3, s->verticesOfCollidingPieces_[i].cols()) = s->verticesOfCollidingPieces_[i];
			neighbor_normals.block(0, k, 3, s->normalsOfCollidingPieces_[i].cols()) = s->normalsOfCollidingPieces_[i];
			k += s->verticesOfCollidingPieces_[i].cols();
		}
	}
	sGPU.add_collinding_points(neighbor_vertices.data(), neighbor_normals.data(), neighbor.data(), total_neighbor);
}

void copy_data_to_cpu(Solver_GPU &sGPU, vector<Solver*> &cpu_solvers)
{
	MatrixXXCuda vel = sGPU.getVelocities();
	MatrixXXCuda *points = sGPU.get_transformedPoints();
	for (pluint i = 0; i <cpu_solvers.size(); ++i)
	{
		cpu_solvers[i]->setVelocities(vel.block(3 * i, 0, 3, vel.cols()));
		cpu_solvers[i]->setPoints(points->block(3 * i, 0, 3, vel.cols()));
	}
}
#endif // GPU

bool file_exist(const std::string& name)
{
	ifstream f(name.c_str());
	return f.good();
}


int main(int argc, char* argv[])
{
	int fake_argc = 1;
	plbInit(&fake_argc, &argv);

    std::string xml_name;
    if (argc > 1) {
		xml_name = std::string(argv[1]);
	}
    else {
		if (plb::global::mpi().isMainProcessor())
            std::cout << "Please enter a parameter file" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (argc > 2) number_nodes = atoi(argv[2]);
	if (argc > 3) number_devices = atoi(argv[3]); // number of devices per node
												  // @Piz-Daint is 1 GPU per node and every hybrid node has 12 cores

	///////////////////////////////////////////////////////////////////////////

	// Make sure that there will be no data type mismatch in MPI communication
	int tsize, tsize_u;
	MPI_Type_size(MPI_LONG_LONG, &tsize);
	MPI_Type_size(MPI_UNSIGNED_LONG_LONG, &tsize_u);

	if ((sizeof(plint) != tsize) || (sizeof(pluint) != tsize_u))
	{
		pcout << "Sizes of data types differ. ABORT. \n";
		exit(EXIT_FAILURE);
	}

	///////////////////////////////////////////////////////////////////////////

	// Open the XMl file.
	XMLreader xmlFile(xml_name);

	xmlFile["System_Setup"]["OutputDir"].read(OutputDir);
	xmlFile["CheckPointing"]["CP_OutputDir"].read(CP_OutputDir);

#ifdef MSVC
	if (CreateDirectory(OutputDir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError())
	{
		;
	}
	else
	{
		// Failed to create directory.
		PLB_ASSERT(false);
	}
	if (CreateDirectory(CP_OutputDir.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError())
	{
		;
	}
	else
	{
		// Failed to create directory.
		PLB_ASSERT(false);
	}
#else
	mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir(CP_OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
	global::directories().setOutputDir(OutputDir);
	global::IOpolicy().activateParallelIO(true);

	///////////////////////////////////////////////////////////////////////////

	xmlFile["CellPacking"]["Run"].read(CellPacking);

	///////////////////////////////////////////////////////////////////////////

	// Shear or Poiseuille Flow
	xmlFile["Simulation_Type"].read(Simulation_Type);

	// Envelopes
	xmlFile["System_Setup"]["particleEnvelopeWidth"].read(particleEnvelopeWidth);

	// Walls/ Boundaries
	xmlFile["System_Setup"]["wallStlFileName"].read(wallStlFileName);

	xmlFile["Case_Study_RBC"].read(Case_Study_RBC);
	xmlFile["Case_Study_PLT"].read(Case_Study_PLT);

	// Simulation:
	// _p  : physical system
	// _lb : lattice/ non-dimensional system
	xmlFile["System_Setup"]["lx_p"].read(lx_p);
	xmlFile["System_Setup"]["ly_p"].read(ly_p);
	xmlFile["System_Setup"]["lz_p"].read(lz_p);
	xmlFile["System_Setup"]["dx_p"].read(dx_p);
	xmlFile["System_Setup"]["dt_p"].read(dt_p);
	xmlFile["System_Setup"]["tau_lb"].read(tau_lb);
	xmlFile["System_Setup"]["u_lb"].read(u_lb);
	xmlFile["System_Setup"]["u_p"].read(u_p);
	xmlFile["System_Setup"]["shear_rate"].read(shear_rate);
	xmlFile["System_Setup"]["inamuro_iT"].read(inamuro_iT);

	xmlFile["Physical_Properties"]["nu_p"].read(nu_p);
	xmlFile["Physical_Properties"]["rho_p"].read(rho_p);

	// --------------------------------------------------------------------- //
	// logfile: dumps output from the main processor
	// logfile_nonparallel: each processor dumps its own output
	// --------------------------------------------------------------------- //

	// Store Data in DATA.log
	plb::global::logfile("DATA.log").flushEntry("Simulation_Type: " + Simulation_Type);
	plb::global::logfile("DATA.log").flushEntry("Case_Study_RBC: " + Case_Study_RBC);
	plb::global::logfile("DATA.log").flushEntry("Case_Study_PLT: " + Case_Study_PLT);
	plb::global::logfile("DATA.log").flushEntry("lx_p: " + util::val2str(lx_p));
	plb::global::logfile("DATA.log").flushEntry("ly_p: " + util::val2str(ly_p));
	plb::global::logfile("DATA.log").flushEntry("lz_p: " + util::val2str(lz_p));
	plb::global::logfile("DATA.log").flushEntry("dx_p: " + util::val2str(dx_p));
	plb::global::logfile("DATA.log").flushEntry("dt_p: " + util::val2str(dt_p));
	plb::global::logfile("DATA.log").flushEntry("tau_lb: " + util::val2str(tau_lb));
	plb::global::logfile("DATA.log").flushEntry("u_lb: " + util::val2str(u_lb));
	plb::global::logfile("DATA.log").flushEntry("u_p: " + util::val2str(u_p));
	plb::global::logfile("DATA.log").flushEntry("shear_rate: " + util::val2str(shear_rate));
	plb::global::logfile("DATA.log").flushEntry("inamuro_iT: " + util::val2str(inamuro_iT));
	plb::global::logfile("DATA.log").flushEntry("nu_p: " + util::val2str(nu_p));
	plb::global::logfile("DATA.log").flushEntry("rho_p: " + util::val2str(rho_p));

	// ShapeOp
	xmlFile["ShapeOp_General"]["timestep"].read(timestep_ShapeOp);
	xmlFile["ShapeOp_General"]["m"].read(m_ShapeOp);
	xmlFile["ShapeOp_General"]["max_iterations"].read(max_iterations_ShapeOp);
	xmlFile["ShapeOp_General"]["max_line_search_loops"].read(max_line_search_loops_ShapeOp);
	xmlFile["ShapeOp_General"]["max_attempts_to_solve_stagnation"].read(max_attempts_to_solve_stagnation_ShapeOp);
	xmlFile["ShapeOp_General"]["convergence_window"].read(convergence_window_ShapeOp);
	xmlFile["ShapeOp_General"]["tol"].read(tol_ShapeOp);
	xmlFile["ShapeOp_General"]["gamma"].read(gamma_ShapeOp);
	xmlFile["ShapeOp_General"]["gamma2"].read(gamma2_ShapeOp);
    xmlFile["ShapeOp_General"]["collisions_threshold_nonRep"].read(collisions_threshold_nonRep);
    xmlFile["ShapeOp_General"]["collisions_weight_nonRep"].read(collisions_weight_nonRep);
    xmlFile["ShapeOp_General"]["collisions_threshold_rep"].read(collisions_threshold_rep);
    xmlFile["ShapeOp_General"]["collisions_weight_rep"].read(collisions_weight_rep);
    xmlFile["ShapeOp_General"]["beta_morse"].read(beta_morse);


	// ShapeOp: RBCs
	xmlFile["ShapeOp_RBC"]["rho"].read(rho_ShapeOp_RBC);
	xmlFile["ShapeOp_RBC"]["Calpha"].read(Calpha_ShapeOp_RBC);
	xmlFile["ShapeOp_RBC"]["Cbeta"].read(Cbeta_ShapeOp_RBC);
	xmlFile["ShapeOp_RBC"]["applyGlobalVolumeConservation"].read(applyGlobalVolumeConservation_ShapeOp_RBC);
    xmlFile["ShapeOp_RBC"]["globalVolumeConservationWeight"].read(globalVolumeConservationWeight_ShapeOp_RBC);
	// ShapeOp: PLTs
	xmlFile["ShapeOp_PLT"]["rho"].read(rho_ShapeOp_PLT);
	xmlFile["ShapeOp_PLT"]["Calpha"].read(Calpha_ShapeOp_PLT);
	xmlFile["ShapeOp_PLT"]["Cbeta"].read(Cbeta_ShapeOp_PLT);
	xmlFile["ShapeOp_PLT"]["applyGlobalVolumeConservation"].read(applyGlobalVolumeConservation_ShapeOp_PLT);
	xmlFile["ShapeOp_PLT"]["globalVolumeConservationWeight"].read(globalVolumeConservationWeight_ShapeOp_PLT);


    T vel_cap, vel_cap_fin;
    xmlFile["System_Setup"]["vel_cap"].read(vel_cap);
    xmlFile["System_Setup"]["vel_cap_fin"].read(vel_cap_fin);

#ifdef NPFEM_CUDA
    Mesh_info params_rbc;

    params_rbc.rho = rho_ShapeOp_RBC;
    params_rbc.dx_p = dx_p;
    params_rbc.threshold_rep = collisions_threshold_rep*dx_p;
    params_rbc.threshold_nonRep = collisions_threshold_nonRep*dx_p;
    params_rbc.weight_col_rep = collisions_weight_rep;
    params_rbc.weight_col_nonRep = collisions_weight_nonRep;
    params_rbc.beta_morse = beta_morse;
    params_rbc.volume_weight = globalVolumeConservationWeight_ShapeOp_RBC;
    params_rbc.calpha = Calpha_ShapeOp_RBC;
    params_rbc.vel_cap = vel_cap;
    params_rbc.vel_cap_fin = vel_cap_fin;

    Mesh_info params_plt;

    params_plt.rho = rho_ShapeOp_PLT;
    params_plt.dx_p = dx_p;
    params_plt.threshold_rep = collisions_threshold_rep*dx_p;
    params_plt.threshold_nonRep = collisions_threshold_nonRep*dx_p;
    params_plt.weight_col_rep = collisions_weight_rep;
    params_plt.weight_col_nonRep = collisions_weight_nonRep;
    params_plt.beta_morse = beta_morse;
    params_plt.volume_weight = globalVolumeConservationWeight_ShapeOp_PLT;
    params_plt.calpha = Calpha_ShapeOp_PLT;
    params_plt.vel_cap = vel_cap;
    params_plt.vel_cap_fin = vel_cap_fin;
#endif // GPU

	///////////////////////////////////////////////////////////////////////////

    // Convert to 1/microsec
    shear_rate = shear_rate / std::pow(10., 6.);

    if (Simulation_Type.compare("Shear_Flow") == 0)
    {
        // Velocity on top ONLY
        u_p = (shear_rate * ly_p) / 1.; // divide by 2. if top & bottom oposite vels
    }
    else
    {
        T D = ly_p; // Tube/Channel radius (just a rough approximation) in micrometers
        // Refers to the wall shear rate
        shear_rate = (4. * u_p) / D; // in 1/microsec
    }

	// See Kruger-LBM-Book p.273
	if (dt_p < 0. && u_lb < 0.) {
		PLB_ASSERT(tau_lb <= 3.);

		dt_p = DESCRIPTOR<T>::cs2 * (tau_lb - 0.5) * dx_p * dx_p / nu_p;
		plb::global::logfile("DATA.log").flushEntry("dt_p: " + util::val2str(dt_p));

		u_lb = u_p * (dt_p / dx_p);
		plb::global::logfile("DATA.log").flushEntry("u_lb: " + util::val2str(u_lb));
		PLB_ASSERT(u_lb <= 0.04);
	}
	else if (tau_lb < 0. && u_lb < 0.) {
		tau_lb = (dt_p * nu_p) / (DESCRIPTOR<T>::cs2 * dx_p * dx_p) + 0.5;
		PLB_ASSERT(tau_lb <= 3.);
		plb::global::logfile("DATA.log").flushEntry("tau_lb: " + util::val2str(tau_lb));

		u_lb = u_p * (dt_p / dx_p);
		plb::global::logfile("DATA.log").flushEntry("u_lb: " + util::val2str(u_lb));
		PLB_ASSERT(u_lb <= 0.04);
	}
	else if (dt_p < 0. && tau_lb < 0.) {
		dt_p = dx_p * (u_lb / u_p);
		plb::global::logfile("DATA.log").flushEntry("dt_p: " + util::val2str(dt_p));

		tau_lb = (dt_p * nu_p) / (DESCRIPTOR<T>::cs2 * dx_p * dx_p) + 0.5;
		PLB_ASSERT(tau_lb <= 3.);
		plb::global::logfile("DATA.log").flushEntry("tau_lb: " + util::val2str(tau_lb));
	}

    nu_lb = nu_p * (dt_p / (dx_p * dx_p));

    // Dimensionless numbers of interest

    // Mach number (smaller than 0.1)
    T lattice_speed = dx_p / dt_p;
    Ma = u_p / lattice_speed;
    
    // This is the classic Reynolds number
    Re = u_p * ly_p / nu_p;
    plb::global::logfile("DATA.log").flushEntry("Re: " + util::val2str(Re));

    // This is the Reynolds number of capsules in shear flows
    T r_RBC = 4.; // the radius of an RBC (just a rough approximation) in micrometers
    Re_rbc = (shear_rate * r_RBC * r_RBC) / nu_p;

    // Capillary number
    T heta = nu_p * rho_p;
    T ks = 5.; // RBC shear elasticity (pN/micrometers)
    Ca = (heta * shear_rate * r_RBC) / ks;

	// This constructor explicitly sets the physical domain. The
	// resolution isn't specified here. You should specify it later on,
	// for example with set_dx() or set_xResolution().
	Units3D units(Array<T, 3>(0., 0., 0.), Array<T, 3>(lx_p, ly_p, lz_p));
	if (Simulation_Type.compare("Shear_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
		units.set_periodic(true, false, true);
	else if (Simulation_Type.compare("Poiseuille_Flow") == 0)
		units.set_periodic(false, false, true);
	units.set_dx(dx_p);

	Box3D lbdomain = units.lbDomain();

	// Discretize the domain with n_ lattice sites per direction
	pluint nx = lbdomain.getNx();
	pluint ny = lbdomain.getNy();
	pluint nz = lbdomain.getNz();

	plb::global::logfile("DATA.log").flushEntry("nx: " + util::val2str(nx));
	plb::global::logfile("DATA.log").flushEntry("ny: " + util::val2str(ny));
	plb::global::logfile("DATA.log").flushEntry("nz: " + util::val2str(nz));

	pcout << "\n" << "***********************" << "\n";
	pcout << "     dx_p (micro-m)  : " << dx_p << std::endl;
	pcout << "     dt_p (micro-s)  : " << dt_p << std::endl;
	pcout << "             tau_lb  : " << tau_lb << std::endl;
	pcout << "               u_lb  : " << u_lb << std::endl;
    pcout << "***********************" << "\n";
    pcout << "    Mach number, Ma  : " << Ma << std::endl;
    pcout << "Reynolds number, Re  : " << Re << std::endl;
    pcout << "         Capsule Re  : " << Re_rbc << std::endl;
    pcout << "Capillary number, Ca : " << Ca << std::endl;
	pcout << "***********************" << "\n";

	///////////////////////////////////////////////////////////////////////////

	// Create the RBC ShapeOp template.
	RBC_shapeOpSolverTemplate.bodyType_ = 0;
	setPointsFromCSV(RBC_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_RBC + "/" + "points.csv"); // Read all points (internal + OnSurface)
	setConstraintsFromCSV(RBC_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_RBC + "/" + "constraints.csv");
	setConnectivityListFromCSV(RBC_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_RBC + "/" + "surface_connectivity.csv"); // Surface Connectivity (tetrahedra are not included)
	setOnSurfaceParticle(RBC_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_RBC + "/" + "onSurfaceParticle.csv");
	Matrix3X RBC_forces(3, RBC_shapeOpSolverTemplate.getPoints().cols());
	RBC_forces.setZero();
	addVertexForce(RBC_shapeOpSolverTemplate, RBC_forces);
	RBC_shapeOpSolverTemplate.initialize(Calpha_ShapeOp_RBC, Cbeta_ShapeOp_RBC, dt_p * (T)(timestep_ShapeOp), rho_ShapeOp_RBC, false, applyGlobalVolumeConservation_ShapeOp_RBC, globalVolumeConservationWeight_ShapeOp_RBC, dx_p, dt_p);

	// Create the PLT ShapeOp template.
	PLT_shapeOpSolverTemplate.bodyType_ = 1;
	setPointsFromCSV(PLT_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_PLT + "/" + "points.csv"); // Read all points (internal + OnSurface)
	setConstraintsFromCSV(PLT_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_PLT + "/" + "constraints.csv");
	setConnectivityListFromCSV(PLT_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_PLT + "/" + "surface_connectivity.csv"); // Surface Connectivity (tetrahedra are not included)
	setOnSurfaceParticle(PLT_shapeOpSolverTemplate, "./Case_Studies/" + Case_Study_PLT + "/" + "onSurfaceParticle.csv");
	Matrix3X PLT_forces(3, PLT_shapeOpSolverTemplate.getPoints().cols());
	PLT_forces.setZero();
	addVertexForce(PLT_shapeOpSolverTemplate, PLT_forces);
	PLT_shapeOpSolverTemplate.initialize(Calpha_ShapeOp_PLT, Cbeta_ShapeOp_PLT, dt_p * (T)(timestep_ShapeOp), rho_ShapeOp_PLT, false, applyGlobalVolumeConservation_ShapeOp_PLT, globalVolumeConservationWeight_ShapeOp_PLT, dx_p, dt_p);

#ifdef NPFEM_CUDA
    params_rbc.miu = 35.0;
    params_rbc.lambda = 5.0;
    params_rbc.kappa = 0.0;

    params_plt.miu = 35.0;
    params_plt.lambda = 5.0;
    params_plt.kappa = 0.0;
#endif // GPU


	numRBC = 0;
	numPLT = 0;
	iniPositions.clear();
	iniRotations.clear();

	xmlFile["System_Setup"]["withBodies"].read(withBodies);
	xmlFile["CheckPointing"]["restartFromCP"].read(restartFromCP);
	xmlFile["CheckPointing"]["checkpoint_timestep"].read(checkpoint_timestep);

	if (withBodies) {
		if (!CellPacking) {
			if (restartFromCP) {
				// Define number of bodies
				xmlFile["System_Setup"]["ht"].read(ht);
				xmlFile["CheckPointing"]["numRBC"].read(numRBC);
				xmlFile["CheckPointing"]["numPLT"].read(numPLT);

				pcout << "\n" << "*********************" << "\n";
				pcout << " ht (%) : " << ht << std::endl;
				pcout << " numRBC : " << numRBC << std::endl;
				pcout << " numPLT : " << numPLT << std::endl;
				pcout << "*********************" << "\n";

				for (pluint i = 0; i < (numRBC + numPLT); ++i)
                {
					if (i < numRBC)
						bodyToType[i] = 0;

					if (i >= numRBC && i < (numRBC + numPLT))
						bodyToType[i] = 1;
				}
			}
			else {
				xmlFile["System_Setup"]["pos_filename"].read(pos_filename);

				pluint iBody = 0;
				std::ifstream pos(pos_filename.c_str());
				std::string line;
				while (std::getline(pos, line)) {
					Array<T, 3> position;
					Array<T, 3> rotation;

					std::istringstream ss(line);
					std::string token;
					int field = 0;
					while (std::getline(ss, token, ',')) {
						if (field == 0) {
							pluint type_ = (pluint)std::stoi(token);
							bodyToType[iBody] = type_;
							if (type_ == 0) ++numRBC;
							if (type_ == 1) ++numPLT;

						}
						else if (field > 0 && field <= 3) {
							position[field - 1] = std::stof(token);
						}
						else {
							rotation[field - 4] = std::stof(token);
						}
						++field;
					}
					iniPositions.push_back(position);
					iniRotations.push_back(rotation);

					++iBody;
				}
			}
		}
		else // CellPacking
        { 
			// Define number of bodies
			xmlFile["System_Setup"]["ht"].read(ht);
			T Vtot;
			if (Simulation_Type.compare("Shear_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0) {
				Vtot = lx_p * ly_p * lz_p;
                pipeRadius_LB = 0.;
			}
			else if (Simulation_Type.compare("Poiseuille_Flow") == 0)
            {
				PI_ = std::acos((T)-1);

                // L = n - 1 & we put the wall 2 lu away from the box -> n - 1 - 2 - 2 
                pluint minN = std::min(nx - 5, ny - 5);
                pipeRadius_LB = (T)0.5 * (T)minN;

				pcout << " ********* \n";
				pcout << "pipeRadius (micro-m) = " << pipeRadius_LB*dx_p << "\n";
				pcout << " ********* \n";

				Vtot = PI_ * (pipeRadius_LB * dx_p) * (pipeRadius_LB * dx_p) * lz_p;
			}

			T Vrbcs = (ht / 100.) * Vtot;
			T Vrbc = RBC_shapeOpSolverTemplate.Volume0_;

			numRBC = std::floor(Vrbcs / Vrbc);
			xmlFile["System_Setup"]["RBCsToPLT"].read(RBCsToPLT);
			numPLT = std::floor(((T)numRBC) / RBCsToPLT);

			pcout << "\n" << "*********************" << "\n";
			pcout << " ht (%) : " << ht << std::endl;
			pcout << " numRBC : " << numRBC << std::endl;
			pcout << " numPLT : " << numPLT << std::endl;
			pcout << "   dt_p : " << dt_p << std::endl;
			pcout << "*********************" << "\n";


			if (!restartFromCP)
            {
				iniPositions.resize(numRBC + numPLT, Array<T, 3>(0., 0., 0.));
				iniRotations.resize(numRBC + numPLT, Array<T, 3>(0., 0., 0.));

				for (pluint i = 0; i < (numRBC + numPLT); ++i) {
					if (i < numRBC) {
						bodyToType[i] = 0;
					}
					if (i >= numRBC && i < (numRBC + numPLT)) {
						bodyToType[i] = 1;
					}
				}

				if (global::mpi().getRank() == 0) {
					// PRNG settings
					std::random_device rd; // seeding
					std::mt19937 mt(rd());

					PI_ = std::acos((T)-1);
					T offset;
					xmlFile["CellPacking"]["offset"].read(offset);
					std::uniform_real_distribution<T> dist_pos_x(0. + offset, lx_p - offset);
					std::uniform_real_distribution<T> dist_pos_y(0. + offset, ly_p - offset);
					std::uniform_real_distribution<T> dist_pos_z(0. + offset, lz_p - offset);
					std::uniform_real_distribution<T> dist_rot(0., 360.0);

                    T radial_pos;
                    T radial_bound = pipeRadius_LB*dx_p - offset;
                    Array<T, 3> pipeCenter((T)0.5 * (T)(nx - 1), (T)0.5 * (T)(ny - 1), (T)0);
                    pipeCenter *= dx_p;

					Array<T, 3> position;
					Array<T, 3> rotation;

					pluint i = 0;
					while (i < (numRBC + numPLT))
                    {
                        if (Simulation_Type.compare("Poiseuille_Flow") == 0)
                        {
                            do
                            {
                                position[0] = dist_pos_x(mt);
                                position[1] = dist_pos_y(mt);
                                radial_pos  = sqrt(pow(position[0] - pipeCenter[0], 2) + pow(position[1] - pipeCenter[1], 2));
                            }
                            while (radial_pos >= radial_bound);
                        }
                        else
                        {
                            position[0] = dist_pos_x(mt);
                            position[1] = dist_pos_y(mt);
                        }

                        position[2] = dist_pos_z(mt);

						rotation[0] = dist_rot(mt);
						rotation[1] = dist_rot(mt);
						rotation[2] = dist_rot(mt);

						iniPositions[i] = position;
						iniRotations[i] = rotation;

						++i;
					}
				}

				MPI_Bcast(&iniPositions[0][0], (int)((numRBC + numPLT) * 3), MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&iniRotations[0][0], (int)((numRBC + numPLT) * 3), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
            else
            {
                pcout << "There is no option to restart from CP in the CellPacking ... ABORT!" << std::endl;
                exit(EXIT_FAILURE);
            }
		}
	}

	///////////////////////////////////////////////////////////////////////////

	// Check if velocity files exist in CPs
	// If yes, then set continuePLB = true

	continuePLB = false;

	xmlFile["System_Setup"]["max_time"].read(max_time);
	xmlFile["System_Setup"]["output_timestep"].read(output_timestep);

	pluint max_iTs = (max_time*std::pow(10., 6.)) / dt_p;
	for (pluint jmp = 0; jmp <= max_iTs; jmp += checkpoint_timestep)
	{
		if (file_exist(createFileName(CP_OutputDir + "vels_body_ID_0_iT_", jmp, fileNamePadding) + ".csv"))
		{
			continuePLB = true;
			break;
		}
	}

	///////////////////////////////////////////////////////////////////////////

	iT = 0;

	if (continuePLB)
    {
		// Check the latest iT
		std::string continueFile = "continue.xml";
		XMLreader cont(continueFile);
		cont["continue"]["iteration"].read(iT);

        if (global::mpi().getRank() == 0)
        {
            std::cout << "\n WARNING: \n"
                      << " You need at least two checkpoints in order to restart. \n"
                      << " The latest checkpoint (CP) is always discarded. \n"
                      << " This is happening to avoid any corruption in the latest CP. \n"
                      << " For example, the simulation times out while you write the CP. \n" 
                      << std::endl;
        }

		iT -= checkpoint_timestep;
	}

    plb::global::logfile("progress.log").flushEntry("initialize ShapeOp Solvers");
	iniShapeOpSolvers("points.csv", "constraints.csv", "surface_connectivity.csv", "onSurfaceParticle.csv");

	///////////////////////////////////////////////////////////////////////////

#ifdef NPFEM_CUDA
	// GPU init
	auto init_gpu = [&](Solver_GPU &sGPU, Solver &s_template, std::vector<Solver*> &solvers, plb::Array<double, 3>*iniPositions, vector<int> *mesh_graph, Mesh_info params) {
		sGPU = Solver_GPU(s_template.getPoints(),
		s_template.getConnectivityList(),
		s_template.get_onSurfaceParticle(),
		s_template.getConstraints(), params, s_template.Cbeta_,
		solvers.size(), dt_p, (T)(timestep_ShapeOp), solvers[0]->bodyID_, mesh_graph);


		for (int i = 0; i< solvers.size(); i++) {
			sGPU.set_gpu_starting_position(solvers[i]->getPoints(), i);
			sGPU.set_gpu_starting_velocities(solvers[i]->getVelocities(), i);
		}
		sGPU.compute_first_centroid();
	};

	if (loc_rank < number_devices)
    {
		cudaError_t err = cudaSetDevice(loc_rank);

		if (shapeOpRBCs.size() > 0)
			init_gpu(sGPU, RBC_shapeOpSolverTemplate, shapeOpRBCs, iniPositions.data(), shapeOpBodies[getSolverID(shapeOpRBCs[0]->bodyID_, shapeOpBodies)].mesh_graph, params_rbc);
		
		if (shapeOpPLTs.size() > 0)
			init_gpu(sGPU_plt, PLT_shapeOpSolverTemplate, shapeOpPLTs, iniPositions.data(), shapeOpBodies[getSolverID(shapeOpPLTs[0]->bodyID_, shapeOpBodies)].mesh_graph, params_plt);
	}
#endif

	///////////////////////////////////////////////////////////////////////////
	T omega = 1. / tau_lb;

	Group3D group;

	if (Simulation_Type.compare("Shear_Flow") == 0) {
		group.add(new MultiBlockLattice3D<T, DESCRIPTOR>(nx, ny, nz,
			new CompleteRegularizedBGKdynamics<T, DESCRIPTOR>(omega)),
			"lattice");
		incompressibleModel = false;
	}
	else if (Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
    {
		// rho_LB = 1.
		Dynamics<T, DESCRIPTOR>* bbDynamics = new BounceBack<T, DESCRIPTOR>(1.);
		group.add(new MultiBlockLattice3D<T, DESCRIPTOR>(
			nx, ny, nz, bbDynamics->clone()),
			"lattice");

		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			0, true);
		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			1, false); // y-direction Walls
		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			2, true);
		defineDynamics(group.getLattice<T, DESCRIPTOR>("lattice"),
			group.getBoundingBox(), bbDynamics->clone());

		Box3D initDomain(group.getBoundingBox());
		initDomain.y0++;
		initDomain.y1--;
        initDomain.y0++;
        initDomain.y1--;

#ifdef FORCED
		Dynamics<T, DESCRIPTOR>* dynamics
			= new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega);
		incompressibleModel = false;
#else
		Dynamics<T, DESCRIPTOR>* dynamics
			= new CompleteRegularizedBGKdynamics<T, DESCRIPTOR>(omega);
		incompressibleModel = false;
#endif
		defineDynamics(group.getLattice<T, DESCRIPTOR>("lattice"), initDomain,
			dynamics->clone());

		delete dynamics;
		dynamics = nullptr;

		delete bbDynamics;
		bbDynamics = nullptr;
	}
	else if (Simulation_Type.compare("Poiseuille_Flow") == 0)
    {
		// rho_LB = 1.
		Dynamics<T, DESCRIPTOR>* bbDynamics = new BounceBack<T, DESCRIPTOR>(1.);
		group.add(new MultiBlockLattice3D<T, DESCRIPTOR>(
			nx, ny, nz, bbDynamics->clone()),
			"lattice");

		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			0, false); // x-direction Walls
		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			1, false); // y-direction Walls
		group.getLattice<T, DESCRIPTOR>("lattice").periodicity().toggle(
			2, true); // z-direction
		defineDynamics(group.getLattice<T, DESCRIPTOR>("lattice"),
			group.getBoundingBox(), bbDynamics->clone());

		Box3D initDomain(group.getBoundingBox());
        // This is to make sure that the outer shell of layers will be
        // for sure bounce back.
		initDomain.x0++;
		initDomain.x1--;
		initDomain.y0++;
		initDomain.y1--;

#ifdef FORCED
		Dynamics<T, DESCRIPTOR>* dynamics
			= new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(omega);
		incompressibleModel = false;
#else
		Dynamics<T, DESCRIPTOR>* dynamics
			= new CompleteRegularizedBGKdynamics<T, DESCRIPTOR>(omega);
		incompressibleModel = false;
#endif
		defineDynamics(group.getLattice<T, DESCRIPTOR>("lattice"), initDomain,
			dynamics->clone());
		delete dynamics;
		dynamics = nullptr;

		// Bounce-Back walls for the pipe.
        // L = n - 1 & we put the wall 2 lu away from the box -> n - 1 - 2 - 2 
		pluint minN = std::min(nx - 5, ny - 5);
		pipeRadius_LB = (T)0.5 * (T)minN;

		pcout << " ********* \n";
		pcout << "pipeRadius (micro-m) = " << pipeRadius_LB*dx_p << "\n";
		pcout << " ********* \n";

		Array<T, 3> pipeCenter_LB((T)0.5 * (T)(nx - 1), (T)0.5 * (T)(ny - 1), (T)0);
		defineDynamics(group.getLattice<T, DESCRIPTOR>("lattice"),
			initDomain, new Pipe<T>(pipeCenter_LB, pipeRadius_LB),
			bbDynamics->clone());

		delete bbDynamics;
		bbDynamics = nullptr;

#ifdef FORCED
		drivingForce_lb = 8. * nu_lb * (u_lb * 0.5) / pipeRadius_LB / pipeRadius_LB;

		setExternalVector(group.getLattice<T, DESCRIPTOR>("lattice"), group.getLattice<T, DESCRIPTOR>("lattice").getBoundingBox(),
			DESCRIPTOR<T>::ExternalField::forceBeginsAt,
			Array<T, 3>((T)0, (T)0, drivingForce_lb));
#endif
	}

	group.generateScalar<T>("rhoBar", 4);
	group.generateTensor<T, 3>("j", 4);
	group.generateTensor<T, 6>("piNeq", 4);
	group.generateDenseParticles<T, DESCRIPTOR>("vertexParticles", particleEnvelopeWidth);

	// Info on the discretized domain
	plb::global::logfile("progress.log").flushEntry(getMultiBlockInfo(group.getLattice<T, DESCRIPTOR>("lattice")));

	plb::global::logfile("progress.log").flushEntry("Injecting Particles");
	iniPalabosParticles<T, DESCRIPTOR>(group, u_lb);

	plb::global::logfile("progress.log").flushEntry("Setting up the flow");
	OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition3D<T, DESCRIPTOR>();
	if (Simulation_Type.compare("Shear_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
	{
		flowSetup(group, *boundaryCondition, u_lb);
	}
	else if (Simulation_Type.compare("Poiseuille_Flow") == 0)
	{
		flowSetup(group);
	}
	plb::global::logfile("progress.log").flushEntry("Flow setup complete");

	///////////////////////////////////////////////////////////////////////////

	if (continuePLB)
    {
		// Re-write continueFile.xml with the one before iT
		std::string xmlFileName = "continue.xml";
		std::string baseFileName = "checkpoint_";
		std::string fname_base = createFileName(baseFileName, iT, fileNamePadding);

		std::vector<MultiBlock3D*> checkpointBlocks;
		checkpointBlocks.push_back(&group.get("lattice"));
		checkpointBlocks.push_back(&group.get("rhoBar"));
		checkpointBlocks.push_back(&group.get("j"));
		checkpointBlocks.push_back(&group.get("piNeq"));
		//checkpointBlocks.push_back(&group.get("vertexParticles"));

		XMLwriter restart;
		XMLwriter& entry = restart["continue"];
		entry["name"].setString(FileName(fname_base).defaultPath(global::directories().getOutputDir()));
		entry["num_blocks"].set(checkpointBlocks.size());
		entry["iteration"].set(iT);
		restart.print(xmlFileName);
	}

	if (continuePLB)
		loadAllBlocks(group, iT);

	///////////////////////////////////////////////////////////////////////////

	// At every restart of the simulation, dump the COMs at different files.

	pluint CP_ID = 0;

	if (continuePLB)
	{
		do
		{
			++CP_ID;
		} while (file_exist(OutputDir + util::val2str(CP_ID) + "_body_" + util::val2str(numRBC) + "_ComPos.log"));
	}

	///////////////////////////////////////////////////////////////////////////

	plb::global::logfile("progress.log").flushEntry("Creating Actions 1 : BoxRhoBarJPiNeqfunctional3D");
	Actions3D actions1;
	plint latticeID_act1 = actions1.addBlock(group.get("lattice"));
	plint rhoBarID_act1 = actions1.addBlock(group.get("rhoBar"));
	plint jID_act1 = actions1.addBlock(group.get("j"));
	plint piNeqID_act1 = actions1.addBlock(group.get("piNeq"));
	actions1.addBlock(group.get("vertexParticles"));
	// 1. Compute density, velocity, and stress tensor.
	plint idRhoBarJ_1 = actions1.addProcessor(new BoxRhoBarJPiNeqfunctional3D<T, DESCRIPTOR>(), latticeID_act1,
		rhoBarID_act1, jID_act1, piNeqID_act1, group.getBoundingBox());
	plint idRhoBarJ_2 = actions1.addCommunication(rhoBarID_act1, modif::staticVariables);
	plint idRhoBarJ_3 = actions1.addCommunication(jID_act1, modif::staticVariables);
	plint idRhoBarJ_4 = actions1.addCommunication(piNeqID_act1, modif::staticVariables);
	actions1.execute(idRhoBarJ_1);
	actions1.execute(idRhoBarJ_2);
	actions1.execute(idRhoBarJ_3);
	actions1.execute(idRhoBarJ_4);

#ifndef FORCED
	Actions3D actionsForcingTerm;
	if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
	{
		plb::global::logfile("progress.log").flushEntry("Creating Actions ForcingTerm");

		if (Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
			pipeRadius_LB = ny - 5;

        // Presure drop (dP/dz = force density)
        drivingForce_lb = 8. * nu_lb * (u_lb * 0.5) / pipeRadius_LB / pipeRadius_LB;
        
        // For pressure driven flow between infinite plates (see Wikipedia)
        if (Simulation_Type.compare("Box_Poiseuille_Flow") == 0)
            drivingForce_lb *= 2.;
		
        Array<T, 3> drivingForce_lb_vec((T)0, (T)0, drivingForce_lb);

		plint latticeID_actFT = actionsForcingTerm.addBlock(group.get("lattice"));
		plint rhoBarID_actFT = actionsForcingTerm.addBlock(group.get("rhoBar"));
		plint jID_actFT = actionsForcingTerm.addBlock(group.get("j"));
		actionsForcingTerm.addProcessor(new AddConstForceToMomentum3D<T, DESCRIPTOR>(drivingForce_lb_vec),
			latticeID_actFT, rhoBarID_actFT, jID_actFT, group.getBoundingBox());
		actionsForcingTerm.addCommunication(jID_actFT, modif::staticVariables);
	}
#endif

    // 2. Compute the force acting on each surface particle, and store it into
	// the corresponding local mesh. This also constructs the local meshes on the
	// fly. Collisions are handled from here.
	plb::global::logfile("progress.log")
		.flushEntry("Creating Actions 2 : ConstructLocalMeshesFromParticles, Combo");
	Actions3D actions2a;
	plint particleID_act2a = actions2a.addBlock(group.get("vertexParticles"));
	actions2a.addProcessor(
		new ConstructLocalMeshesFromParticles<T, DESCRIPTOR>(rbcTemplate, pltTemplate, bodyToType, particleEnvelopeWidth),
		particleID_act2a, group.getBoundingBox());

    Actions3D actions2b;
    plint rhoBarID_act2b = actions2b.addBlock(group.get("rhoBar"));
    plint piNeqID_act2b = actions2b.addBlock(group.get("piNeq"));
    plint particleID_act2b = actions2b.addBlock(group.get("vertexParticles"));
	// From lattice to physical units
	T Cf = rho_p * (dx_p * dx_p * dx_p * dx_p) / (dt_p * dt_p); // force
	T Cp = rho_p * (dx_p * dx_p) / (dt_p * dt_p); // pressure
	T Ca = dx_p * dx_p; // area
    T densityOffset = 1.0;
	actions2b.addProcessor(new CollisionsForcesCombo<T, DESCRIPTOR>(dx_p, omega, densityOffset, Cf, Cp, Ca, collisions_threshold_rep, collisions_threshold_nonRep, wallVertexNormals, bodyToType, CellPacking),
		particleID_act2b, rhoBarID_act2b, piNeqID_act2b, group.getBoundingBox());

	// 3. Immersed boundary method.
	plb::global::logfile("progress.log").flushEntry("Creating Actions 3 : IBM");
	Actions3D actions3;
	plint rhoBarID_act3 = actions3.addBlock(group.get("rhoBar"));
	plint jID_act3 = actions3.addBlock(group.get("j"));
	plint particleID_act3 = actions3.addBlock(group.get("vertexParticles"));
	for (pluint i = 0; i < inamuro_iT; i++)
	{
		actions3.addProcessor(
			new MultiDirectForcingImmersedBoundaryIteration3D<T,
			DESCRIPTOR>(1. / omega, incompressibleModel, wallMesh),
			particleID_act3, rhoBarID_act3, jID_act3,
			group.getBoundingBox());
		actions3.addCommunication(jID_act3, modif::staticVariables);
	}

	// 4. Collision and streaming.
	plb::global::logfile("progress.log").flushEntry("Creating Actions 4 : Collide and Stream");
	Actions3D actions4;
	plint latticeID_act4 = actions4.addBlock(group.get("lattice"));
	plint rhoBarID_act4 = actions4.addBlock(group.get("rhoBar"));
	plint jID_act4 = actions4.addBlock(group.get("j"));
	actions4.addProcessor(new ExternalRhoJcollideAndStream3D<T, DESCRIPTOR>(),
		latticeID_act4, rhoBarID_act4, jID_act4, group.getBoundingBox());
	plint level = 0;
	actions4.addInternalProcessors(latticeID_act4, level); // Add the boundary conditions that were added to the lattices.
	actions4.addCommunication(latticeID_act4, modif::staticVariables);

	// 5. Give new velocities to the particles
	plb::global::logfile("progress.log").flushEntry("Creating Actions 5 : LocalMeshToParticleVelocity3D");
	Actions3D actions5;
	plint particleID_act5 = actions5.addBlock(group.get("vertexParticles"));
	actions5.addProcessor(new LocalMeshToParticleVelocity3D<T, DESCRIPTOR>(), particleID_act5, group.getBoundingBox());
	actions5.addCommunication(particleID_act5, modif::dynamicVariables);

	// 6. Advance the particles
	plb::global::logfile("progress.log").flushEntry("Creating Actions 6 : AdvanceParticlesEveryWhereFunctional3D");
	Actions3D actions6;
	plint particleID_act6 = actions6.addBlock(group.get("vertexParticles"));
	T cutoffSpeedSqr = -1.0;
	actions6.addProcessor(
		new AdvanceParticlesEveryWhereFunctional3D<T, DESCRIPTOR>(
			cutoffSpeedSqr),
		particleID_act6, group.getBoundingBox());
	actions6.addCommunication(particleID_act6, modif::dynamicVariables);

	///////////////////////////////////////////////////////////////////////////

	// allocate memory for output dump
	centers_log = MatrixXX(output_timestep / timestep_ShapeOp, 3 * shapeOpBodies.size());

	// if restart from a CP then the first output is discarded
    // because it is the corresponding zero time and thus it is a void one.
	bool discardFlag = true;
	if (continuePLB)
		discardFlag = false;

	//////////////////////////////////////////////////////////////////////

	// Output related
	bool printVTKs;
	xmlFile["System_Setup"]["printVTKs"].read(printVTKs);

    xmlFile["System_Setup"]["dump_RBCs_ratio_VTKs"].read(dump_RBCs_ratio_VTKs);
    xmlFile["System_Setup"]["dump_PLTs_ratio_VTKs"].read(dump_PLTs_ratio_VTKs);
    xmlFile["System_Setup"]["dump_RBCs_ratio_COMs"].read(dump_RBCs_ratio_COMs);
    xmlFile["System_Setup"]["dump_PLTs_ratio_COMs"].read(dump_PLTs_ratio_COMs);

    pluint numRBCs_toPrint_COMs = (pluint)((T)shapeOpRBCs.size() * dump_RBCs_ratio_COMs);
    pluint numPLTs_toPrint_COMs = (pluint)((T)shapeOpPLTs.size() * dump_PLTs_ratio_COMs);

	// how many CPs so far
	pluint CPs = 0;

	// Loop over main time iteration.
	plb::global::logfile("progress.log").flushEntry("\n\nStarting Main Time Loop");

	// At every iteration we want to advance both the bodies and the fluid from t to t+1
	// Thus we choose the sequence below:
	// 1. Compute density, velocity and stress fields (t)
	// 2. Fext(t) on solid bodies (& colliding neighbors)
	// 3. Solve the bodies with Fext(t) - But do not advance yet at t+1
	// 4. Change the momentum of the fluid using IBM & any other external force (Shan-Chen forcing scheme).
	//    The fluid needs to know the position & velocity of the immersed body at time t and based on this
	//    we are going to advance to t+1. Our schemes here are explicit and so we compute t+1 with the info
	//    that we have at t. This is why we did not advance at the previous step
	// 5. Collide & Stream (FLUID at t+1)
	// 6. Advance particles and copy new vels (SOLID at t+1)

	timer("all").start();
	while (true)
    {
		if (iT>0 && iT%output_timestep == 0 && discardFlag)
			time_statistics(iT, output_timestep);

        timer("logging").start();
		plb::global::logfile("progress.log").flushEntry("\n\nStarting Time Loop " + util::val2str(iT));
        timer("logging").stop();

#ifdef ENABLE_LOGS
		plb::global::logfile_nonparallel("localProgress.log").flushEntry("\n\nStarting Time Loop " + util::val2str(iT));
		plb::global::logfile_nonparallel("particles.log").flushEntry("\n\nStarting Time Loop " + util::val2str(iT));
#endif

        timer("logging").start();
		plb::global::logfile("progress.log").flushEntry("Reset communication data");
        timer("logging").stop();
		timer("reset-communication").start();
		resetCommunication(); // Reset all communication data.
		timer("reset-communication").stop();

		if (!CellPacking) 
        {
            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Execute Actions 1");
            timer("logging").stop();
			// Execute a time iteration.
			timer("actions1").start();
			actions1.execute();
			// - Compute density, velocity, and stresses (all fields refer at time t)
			timer("actions1").stop();
		}

        timer("logging").start();
		plb::global::logfile("progress.log").flushEntry("Execute Actions 2a");
        timer("logging").stop();
		timer("actions2a").start();
		actions2a.execute();
		// - Construct local meshes
        // need this for Area Computation (IBM)
		timer("actions2a").stop();
        
		// ShapeOp is executed every timestep_ShapeOp steps
		// Solved but not yet fully advanced at t+1 (SO side at t+1 but not Palabos side which is still at t)
		if (iT%timestep_ShapeOp == 0)
        {
            timer("logging").start();
            plb::global::logfile("progress.log").flushEntry("Execute Actions 2b");
            timer("logging").stop();
            timer("actions2b").start();
            actions2b.execute();
            // - Collisions detection
            // - Copy forces to local meshes (Fext(t))
            timer("actions2b").stop();

            // If temporal avg THEN communicate forces at every step
            timer("logging").start();
            plb::global::logfile("progress.log").flushEntry("Send forces and Collision Containers (Palabos -> ShapeOp)");
            timer("logging").stop();
            timer("mpi").start();
            sendBodyForcesAndCollisionData<T>(bodyToProc, shapeOpBodies); // Send forces from Palabos local meshes to ShapeOp
            timer("mpi").stop();

            timer("logging").start();
            plb::global::logfile("progress.log").flushEntry("Receive forces and Collision Containers (Palabos -> ShapeOp)");
            timer("logging").stop();
            timer("mpi").start();
            receiveBodyForcesAndCollisionData<T>(shapeOpBodies); // Receive forces on the ShapeOp end.
            timer("mpi").stop();

            // After mpi-communication with Isend, we need to wait and be sure
            // that all data we want are in their place.
            timer("logging").start();
            plb::global::logfile("progress.log").flushEntry("MPI completion for Containers");
            timer("logging").stop();
            timer("mpi").start();
            forceMpiCompletion();
            timer("mpi").stop();

            timer("logging").start();
            plb::global::logfile("progress.log").flushEntry("Apply Containers");
            timer("logging").stop();
            timer("apply-containers").start();
            for (pluint i = 0; i < shapeOpBodies.size(); ++i)
            {
                shapeOpBodies[i].applyForces(CellPacking);
                shapeOpBodies[i].applyCollidingNeighbors();
            }
            timer("apply-containers").stop();

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Run the ShapeOp solver");
            timer("logging").stop();
			timer("shapeop").start();
#ifdef NPFEM_CUDA
			if (shapeOpRBCs.size()) {
				copy_data_to_gpu(sGPU, shapeOpRBCs);
				sGPU.Palabos_iT_ = iT;
				sGPU.solve(max_iterations_ShapeOp, tol_ShapeOp, true, gamma_ShapeOp, max_line_search_loops_ShapeOp, m_ShapeOp, gamma2_ShapeOp);
				// GPU counterpart of imposePlbPeriodicityToSO
				sGPU.make_periodic(group.getLattice<T, DESCRIPTOR>("lattice").getNx() - 1,
					group.getLattice<T, DESCRIPTOR>("lattice").getNy() - 1,
					group.getLattice<T, DESCRIPTOR>("lattice").getNz() - 1, dx_p, iT);
				copy_data_to_cpu(sGPU, shapeOpRBCs);
				sGPU.free_collinding_points();

				double *centers = sGPU.get_centers();
				int n_cells = sGPU.get_nb_cells();
				int i = ((iT - timestep_ShapeOp) % output_timestep) / timestep_ShapeOp;
				for (int j = 0; j < n_cells * 3; j += 3) {
					centers_log(i, j) = centers[j];
					centers_log(i, j + 1) = centers[j + 1];
					centers_log(i, j + 2) = centers[j + 2];
				}
			}

			if (shapeOpPLTs.size()) {

				copy_data_to_gpu(sGPU_plt, shapeOpPLTs);
				sGPU_plt.Palabos_iT_ = iT;
				sGPU_plt.solve(max_iterations_ShapeOp, tol_ShapeOp, true, gamma_ShapeOp, max_line_search_loops_ShapeOp, m_ShapeOp, gamma2_ShapeOp);
				// GPU counterpart of imposePlbPeriodicityToSO
				sGPU_plt.make_periodic(group.getLattice<T, DESCRIPTOR>("lattice").getNx() - 1,
					group.getLattice<T, DESCRIPTOR>("lattice").getNy() - 1,
					group.getLattice<T, DESCRIPTOR>("lattice").getNz() - 1, dx_p, iT);
				copy_data_to_cpu(sGPU_plt, shapeOpPLTs);
				sGPU_plt.free_collinding_points();

				double *centers = sGPU_plt.get_centers();
				int last_cells = sGPU_plt.get_nb_cells() + shapeOpRBCs.size();
				int i = ((iT - timestep_ShapeOp) % output_timestep) / timestep_ShapeOp;
				int k = 0;
				for (int j = 3 * shapeOpRBCs.size(); j < 3 * last_cells; j += 3) {
					centers_log(i, j) = centers[k];
					centers_log(i, j + 1) = centers[k + 1];
					centers_log(i, j + 2) = centers[k + 2];
					k += 3;
				}
			}
#else
			Vector3 COM;
			pluint i = ((iT - timestep_ShapeOp) % output_timestep) / timestep_ShapeOp;
			pluint j = 0;
			for (auto& so : shapeOpRBCs)
			{
				so->Palabos_iT_ = iT;
				so->solve(m_ShapeOp
					, max_iterations_ShapeOp
					, max_line_search_loops_ShapeOp
					, max_attempts_to_solve_stagnation_ShapeOp
					, convergence_window_ShapeOp
					, tol_ShapeOp
					, gamma_ShapeOp
					, gamma2_ShapeOp
					, collisions_threshold_rep*dx_p
					, collisions_weight_rep
                    , collisions_threshold_nonRep*dx_p
                    , collisions_weight_nonRep
                    , beta_morse
					, vel_cap
                    , vel_cap_fin);

				COM = so->getPoints().rowwise().mean();

				centers_log(i, j) = COM[0];
				centers_log(i, j + 1) = COM[1];
				centers_log(i, j + 2) = COM[2];

				j += 3;
			}

			for (auto& so : shapeOpPLTs)
			{
				so->Palabos_iT_ = iT;
				so->solve(m_ShapeOp
					, max_iterations_ShapeOp
					, max_line_search_loops_ShapeOp
					, max_attempts_to_solve_stagnation_ShapeOp
					, convergence_window_ShapeOp
					, tol_ShapeOp
					, gamma_ShapeOp
					, gamma2_ShapeOp
                    , collisions_threshold_rep*dx_p
                    , collisions_weight_rep
                    , collisions_threshold_nonRep*dx_p
                    , collisions_weight_nonRep
                    , beta_morse
					, vel_cap
                    , vel_cap_fin);

				COM = so->getPoints().rowwise().mean();

				centers_log(i, j) = COM[0];
				centers_log(i, j + 1) = COM[1];
				centers_log(i, j + 2) = COM[2];

				j += 3;
			}
#endif
			timer("shapeop").stop();

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Send velocities (ShapeOp -> PLB)");
            timer("logging").stop();
			timer("mpi").start();
			sendBodiesVelocities<T>(shapeOpBodies);
			timer("mpi").stop();

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Receive velocities (ShapeOp -> PLB)");
            timer("logging").stop();
			timer("mpi").start();
			receiveBodiesVelocities<T>(bodyToProc, shapeOpBodies, dx_p, dt_p, rho_p);
			timer("mpi").stop();

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("MPI completion for velocities");
            timer("logging").stop();
			timer("mpi").start();
			velocityMpiCompletion();
			timer("mpi").stop();
		}

		// Add forces (everything refers at time t)
		if (!CellPacking)
        {
#ifndef FORCED
			if (Simulation_Type.compare("Poiseuille_Flow") == 0 || Simulation_Type.compare("Box_Poiseuille_Flow") == 0) {
				timer("actionsForcingTerm").start();
				actionsForcingTerm.execute();
				timer("actionsForcingTerm").stop();
			}
#else
            timer("actionsForcingTerm").start();
			setExternalVector(group.getLattice<T, DESCRIPTOR>("lattice"), group.getLattice<T, DESCRIPTOR>("lattice").getBoundingBox(),
				DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>((T)0, (T)0, drivingForce_lb));
            timer("actionsForcingTerm").stop();
#endif

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Execute Actions 3");
            timer("logging").stop();
			timer("actions3").start();
			actions3.execute();
			timer("actions3").stop();
			// - Implement the immersed boundary method

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Execute Actions 4");
            timer("logging").stop();
			timer("actions4").start();
			actions4.execute();
			timer("actions4").stop();
			// - Collide And Stream (FLUID at t+1)
		}

		if (iT%timestep_ShapeOp == 0)
		{
            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Execute Actions 5");
            timer("logging").stop();
			timer("actions5").start();
			actions5.execute();
			timer("actions5").stop();
			// - Copy velocities from local meshes to particle velocities
		}

        timer("logging").start();
		plb::global::logfile("progress.log").flushEntry("Execute Actions 6");
        timer("logging").stop();
		timer("actions6").start();
		actions6.execute();
		timer("actions6").stop();
		// - Advance the particles (SOLID at t+1)

#ifndef NPFEM_CUDA
		if (iT%timestep_ShapeOp == 0)
		{
            timer("shapeop").start();
			for (pluint i = 0; i < shapeOpBodies.size(); ++i)
			{
				shapeOpBodies[i].imposePlbPeriodicityToSO(dx_p,
					group.getLattice<T, DESCRIPTOR>("lattice").getNx(),
					group.getLattice<T, DESCRIPTOR>("lattice").getNy(),
					group.getLattice<T, DESCRIPTOR>("lattice").getNz());
			}
            timer("shapeop").stop();
		}
#endif

        timer("logging").start();
		plb::global::logfile("progress.log").flushEntry("\nEnd of Cycle");
        timer("logging").stop();
#ifdef ENABLE_LOGS
		plb::global::logfile_nonparallel("localProgress.log").flushEntry("\nEnd of Cycle");
#endif

		///////////////////////////////////////////////////////////////////////

		if (iT%output_timestep == 0 && discardFlag)
		{
            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Writing VTK files & COMs");
            timer("logging").stop();

			timer("output").start();

            //if (printVTKs)
            //    writeVTK(group.getLattice<T, DESCRIPTOR>("lattice"), dx_p, dt_p, iT);

			if (shapeOpBodies.size())
			{
				// Write the Centers of Mass (COMs)
				if (iT > 0)
				{
					pluint bodyID_tmp;
					plint cnt = 0;
                   
                    bool print = true, RBC_printing = false, PLT_printing = false;
                    pluint numRBCs_toPrint_tmp = 0;
                    pluint numPLTs_toPrint_tmp = 0;

					for (pluint j = 0; j < centers_log.cols(); j += 3)
					{
						if ((j / 3) < shapeOpRBCs.size())
						{
                            RBC_printing = true;
                            PLT_printing = false;
                            ++numRBCs_toPrint_tmp;

							bodyID_tmp = shapeOpRBCs[cnt]->bodyID_;

							if ((j / 3) == (shapeOpRBCs.size() - 1))
								cnt = -1;
						}
						else
						{
                            RBC_printing = false;
                            PLT_printing = true;
                            ++numPLTs_toPrint_tmp;

							bodyID_tmp = shapeOpPLTs[cnt]->bodyID_;
						}

                        if (RBC_printing)
                            if (numRBCs_toPrint_tmp > numRBCs_toPrint_COMs)
                                print = false;
                            else
                                print = true;
                        
                        if (PLT_printing)
                            if (numPLTs_toPrint_tmp > numPLTs_toPrint_COMs)
                                print = false;
                            else
                                print = true;

                        if (print)
                        {
                            std::ofstream ofile;
                            ofile.open(OutputDir + util::val2str(CP_ID) + "_body_" + plb::util::val2str(bodyID_tmp) + "_ComPos.log", std::ofstream::out | std::ofstream::app);
                            for (int i = 0; i < centers_log.rows(); i++)
                            {
                                ofile << plb::util::val2str(iT - output_timestep + i*timestep_ShapeOp + timestep_ShapeOp) + "," +
                                    plb::util::val2str(centers_log(i, j)) + "," +
                                    plb::util::val2str(centers_log(i, j + 1)) + "," +
                                    plb::util::val2str(centers_log(i, j + 2))
                                    << std::endl;
                            }
                            ofile.close();
                        }

						++cnt;
					}
				}

				if (printVTKs)
				{
					// For visualization purposes, it is better to have them separately
					// Write Meshes in lattice units
					RawConnectedTriangleMesh<T> RBCs_mesh = shapeOpBodies[0].mergeMultipleMeshes(shapeOpRBCs, shapeOpPLTs_empty, dx_p, dump_RBCs_ratio_VTKs, dump_PLTs_ratio_VTKs);
					RawConnectedTriangleMesh<T> PLTs_mesh = shapeOpBodies[0].mergeMultipleMeshes(shapeOpRBCs_empty, shapeOpPLTs, dx_p, dump_RBCs_ratio_VTKs, dump_PLTs_ratio_VTKs);

					// RBCs
					multiProcWriteVTK(RBCs_mesh,
						createFileName(OutputDir + util::val2str(global::mpi().getRank()) + "_RBCs_iT_", iT, fileNamePadding) + ".vtk",
						group.getLattice<T, DESCRIPTOR>("lattice").getNx(),
						group.getLattice<T, DESCRIPTOR>("lattice").getNy(),
						group.getLattice<T, DESCRIPTOR>("lattice").getNz(),
						dx_p, dt_p, rho_p);
					// PLTs
					multiProcWriteVTK(PLTs_mesh,
						createFileName(OutputDir + util::val2str(global::mpi().getRank()) + "_PLTs_iT_", iT, fileNamePadding) + ".vtk",
						group.getLattice<T, DESCRIPTOR>("lattice").getNx(),
						group.getLattice<T, DESCRIPTOR>("lattice").getNy(),
						group.getLattice<T, DESCRIPTOR>("lattice").getNz(),
						dx_p, dt_p, rho_p);
				}
			}

			timer("output").stop();
		}

		// Access averages from internal statistics ( their value is defined only after the call to lattice.collideAndStream() )
		if (iT%output_timestep == 0)
		{
			timer("output").start();

			pcout << "step " << iT << "; t(micro-s)= " << iT * dt_p;
			pcout << "; av energy=" << setprecision(10)
				<< computeAverageEnergy<T>(
					group.getLattice<T, DESCRIPTOR>("lattice"))
				<< "; av rho=" << setprecision(10)
				<< computeAverageDensity<T>(
					group.getLattice<T, DESCRIPTOR>("lattice"))
				<< "; max uLB=" << setprecision(10)
				<< computeMax(*computeVelocityNorm(
					group.getLattice<T, DESCRIPTOR>("lattice")))
				<< "; max u_p(micro-m/micro-s)="
				<< computeMax(*computeVelocityNorm(
					group.getLattice<T, DESCRIPTOR>("lattice"))) * (dx_p / dt_p)
				<< endl;

            timer("logging").start();
			plb::global::logfile("energy.log").flushEntry(util::val2str(computeAverageEnergy<T>(group.getLattice<T, DESCRIPTOR>("lattice"))));
            timer("logging").stop();

			timer("output").stop();
		}

		// Checkpointing
		if (iT%checkpoint_timestep == 0 && iT>0 && discardFlag)
		{
			timer("output").start();

            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Checkpointing files");
            timer("logging").stop();

			// Palabos related ones
            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("Palabos related ones");
            timer("logging").stop();
            if (!CellPacking)
			    saveAllBlocks(group, iT);

			// ShapeOp related ones
            timer("logging").start();
			plb::global::logfile("progress.log").flushEntry("ShapeOp related ones");
            timer("logging").stop();
            if (!CellPacking)
            {
                for (pluint i = 0; i < shapeOpBodies.size(); ++i)
                {
                    savePointsToCSV(shapeOpBodies[i].getSolver(), createFileName(CP_OutputDir + "body_ID_" + util::val2str(shapeOpBodies[i].getID()) + "_iT_", iT, fileNamePadding) + ".csv", iT, timestep_ShapeOp);
                    saveVelsToCSV(shapeOpBodies[i].getSolver(), createFileName(CP_OutputDir + "vels_body_ID_" + util::val2str(shapeOpBodies[i].getID()) + "_iT_", iT, fileNamePadding) + ".csv");
                }
            }
            else
            {
                for (pluint i = 0; i < shapeOpBodies.size(); ++i)
                {
                    savePointsToCSV(shapeOpBodies[i].getSolver(), CP_OutputDir + "body_ID_" + util::val2str(shapeOpBodies[i].getID()) + ".csv", iT, timestep_ShapeOp);
                }
            }

            ++CPs;

			// Do a cleaning
			if (CPs > 2 && !CellPacking)
			{
				pluint tmpID;
                pluint iter = iT - 2 * checkpoint_timestep;

				for (auto& bd : shapeOpBodies)
				{
					tmpID = bd.getID();
					remove((createFileName(CP_OutputDir + "body_ID_" + util::val2str(tmpID) + "_iT_", iter, fileNamePadding) + ".csv").c_str());
					remove((createFileName(CP_OutputDir + "vels_body_ID_" + util::val2str(tmpID) + "_iT_", iter, fileNamePadding) + ".csv").c_str());
				}

                if (global::mpi().getRank() == 0)
                {
                    for (pluint i = 0; i < numOfBlocks; ++i)
                    {
                        remove((createFileName(OutputDir + "checkpoint_", iter, fileNamePadding) + "_" + util::val2str(i) + ".dat").c_str());
                        remove((createFileName(OutputDir + "checkpoint_", iter, fileNamePadding) + "_" + util::val2str(i) + ".plb").c_str());
                    }
                }
			}

            timer("output").stop();
		}

		if ((iT*dt_p) >= max_time*std::pow(10., 6.))
			break; // Stopping Criterion

		discardFlag = true;

		++iT;
	} // End of time loop

	delete boundaryCondition;
}
