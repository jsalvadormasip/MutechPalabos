///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifndef DEVICE_UTILITIES_H
#define DEVICE_UTILITIES_H
///////////////////////////////////////////////////////////////////////////////
//#include "Types.h"
#include "npFEM/common.h"
#include "npFEM/GPU_data.h"
#include "npFEM/Constraint_Flattening.h"
#include <iostream>
#include <cstdio>
#include <string.h>
#include <vector>

#include "npFEM/sparse_matrix.h"

namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
inline void CUDA_HandleError( cudaError_t err, const char *file, int line );
#define CUDA_HANDLE_ERROR( err ) (CUDA_HandleError( err, __FILE__, __LINE__ ))
///////////////////////////////////////////////////////////////////////////////
void GPU_Mem_check();
///////////////////////////////////////////////////////////////////////////////
void GPU_Init(Mesh_info        *mesh_info_,
			  Mesh_data        *mesh_data_d,        Mesh_data *mesh_data_h,
              Simulation_input *simulation_input_d, Simulation_input *simulation_input_h,
              Simulation_data  *simulation_data_d ,
              Collision_data   *collision_data_d  , Collision_data   *collision_data_h, cuda_scalar **matrices_d);

///////////////////////////////////////////////////////////////////////////////
void compute_next_frame_rbc(Mesh_info  *info, Mesh_data *mesh, Simulation_input *input,
	Simulation_data *sim, Collision_data *coll, ShapeOpScalar h, ShapeOpScalar h2, int it_max, cudaStream_t stream);
///////////////////////////////////////////////////////////////////////////////
struct gradient_functor;
///////////////////////////////////////////////////////////////////////////////
void median_filter(graph_data graph, ShapeOpScalar *force, int nb_cells, int nb);
///////////////////////////////////////////////////////////////////////////////
void send_graph(graph_data *graph_h, graph_data *graph_d, int n);
///////////////////////////////////////////////////////////////////////////////
void first_center(Mesh_info info, Simulation_input input, Simulation_data sim, double *center_d);
///////////////////////////////////////////////////////////////////////////////
void set_cells_initial_position(const double *center_h, double *center_d, const int nb_cells, Mesh_info info, Mesh_data mesh, Simulation_input input, Simulation_data sim);
///////////////////////////////////////////////////////////////////////////////
void reset_position_d(const double *center_h, double *center_d, ShapeOpScalar *points, Mesh_info info, Mesh_data mesh, Simulation_input input_d, Simulation_data sim);
///////////////////////////////////////////////////////////////////////////////
void shift_cells_points_g(ShapeOpScalar *points, ShapeOpScalar *center_h, ShapeOpScalar *center_d, int n, int nb_cells);
///////////////////////////////////////////////////////////////////////////////
void send_GPU_collinding_points(ShapeOpScalar *points, ShapeOpScalar **colid_points_d, ShapeOpScalar *normals, ShapeOpScalar **colid_normals_d, int *n_per_cell, int *n_per_cell_d, int nb, int nb_cell);
///////////////////////////////////////////////////////////////////////////////
void free_GPU_pointer(void *pointer);
///////////////////////////////////////////////////////////////////////////////
void external_forces_from_Host_to_Device( int n_points, ShapeOpScalar *Palabos_Forces_h, ShapeOpScalar *Palabos_Forces_d, cudaStream_t stream);
///////////////////////////////////////////////////////////////////////////////
void points_from_Host_to_Device(int n_points, ShapeOpScalar *points_d, ShapeOpScalar *points_h, int cell);
///////////////////////////////////////////////////////////////////////////////
void points_from_Device_to_Host( int n_points, ShapeOpScalar *points_d, ShapeOpScalar *points_h, cudaStream_t stream);
///////////////////////////////////////////////////////////////////////////////
void debug_matrix_from_gpu(double *mat_d, int n);
///////////////////////////////////////////////////////////////////////////////
void points_time_mat_3X3(cuda_scalar *mat_d, const cuda_scalar *mat_h, double *points, const int n_cells, const int n);
///////////////////////////////////////////////////////////////////////////////
void device_make_periodic(ShapeOpScalar *points_d, ShapeOpScalar *center, float nx, float ny, float nz, int ncells, int npoints,int it);
///////////////////////////////////////////////////////////////////////////////
}
}
#endif // DEVICE_UTILITIES_H
///////////////////////////////////////////////////////////////////////////////