#ifndef SHAPEOPWRAPPER_H
#define SHAPEOPWRAPPER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "npFEM/Solver.h"
#ifndef NPFEM_SA
#ifdef NPFEM_CUDA
#include "npFEM/Solver_GPU.h"
#endif
#endif // !NPFEM_SA
#include "npFEM/Constraint.h"
#include "npFEM/Force.h"

namespace plb {

typedef plb::npfem::Solver ShapeOp_Solver;

void setPointsFromCSV(ShapeOp_Solver& s, std::string filename, bool best_fit_T = false);
void setVelsFromCSV(ShapeOp_Solver& s, std::string filename);

void savePointsToCSV(ShapeOp_Solver& s, std::string filename, size_t iT = 1, size_t dt_ShapeOp = 1);
void saveVelsToCSV(ShapeOp_Solver& s, std::string filename);

void setConstraintsFromCSV(ShapeOp_Solver& s, std::string filename);

void setConnectivityListFromCSV(ShapeOp_Solver& s, std::string filename);

void setForcesFromCSV(ShapeOp_Solver& s, std::string filename);

void saveForcesToCSV(ShapeOp_Solver& s, std::string filename);

void setOnSurfaceParticle(ShapeOp_Solver& s, std::string filename);

// This runs only once as an initialization. Afterwards, modify forces with
// editVertexForce
void addVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id = 0);

void editVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id = 0);

}

///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "shapeOpWrapper.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif
