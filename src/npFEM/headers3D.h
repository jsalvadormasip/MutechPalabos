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
 * Groups all the include files for npFEM.
 */
#include "npFEM/Constraint.h"
#include "npFEM/Force.h"
#include "npFEM/LSSolver.h"
#include "npFEM/Solver.h"
#include "npFEM/Types.h"
#include "npFEM/rbcCommunicate.h"
#include "npFEM/rbcFiltering.h"
#include "npFEM/rbcGlobal.h"
#include "npFEM/rbcParticle.h"
#include "npFEM/rbcShapeOp.h"
#include "npFEM/shapeOpWrapper.h"
#include "npFEM/visualization.h"
#include "npFEM/common.h"

#ifdef NPFEM_CUDA
#include "npFEM/Constraint_Flattening.h"
#include "npFEM/GPU_data.h"
#include "npFEM/Solver_GPU.h"
#include "npFEM/device_utilities.h"
#include "npFEM/projections_GPU_MATH.h"
#include "npFEM/projections_GPU_soa.h"
#include "npFEM/quasy_newton.h"
#include "npFEM/quasy_newton3.h"
#include "npFEM/sparse_matrix.h"
#include "npFEM/sum_cuda.h"
#include "npFEM/svd3_cuda.h"
#endif //NPFEM_CUDA