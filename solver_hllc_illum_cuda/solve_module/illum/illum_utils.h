#ifndef ILLUM_UTILS_H
#define ILLUM_UTILS_H

#include "../../global_def.h"
#include "../../global_headers.h"
#include "../solve_config.h"
#include "../solve_global_struct.h"
#if defined ILLUM 

int GetDirectionIllumFromFace(const int size_grid, const int num_dir, const Type* illum_on_face, std::vector<Type>& illum_in_cell);
#endif

#if defined ILLUM && defined SOLVE
#ifdef USE_CUDA
void CalculateParamOnCuda(const grid_directions_t& grid_dir, grid_t& grid);
#endif

int InitIllum(file_name main_dir, grid_t& grid);

int CalculateIllum(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x, const std::vector < std::vector<int>>& sorted_id_cell,
	//const std::vector<Type>& res_inner_bound, 
	grid_t& grid, std::vector<Type>& Illum, std::vector<Type>& int_scattering);

int SolveIllumAndHLLC(const Type tau, grid_t& grid);
int CalculateIllumParam(const grid_directions_t& grid_direction, grid_t& grid);

Type BoundaryConditions(const int type_bound, Vector3& inter_coef);

#ifdef USE_MPI
int InitSendDispIllumArray(const int myid, const int np, const int count_directions, const int count_cells);
void MPI_INIT(const int myid, const int np, const int count_directions, const grid_t& grid);
int InitPhysOmpMpi(const int count_cells);

int MPI_CalculateIllum(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x, const std::vector < std::vector<int>>& sorted_id_cell,
	//const std::vector<Type>& res_inner_bound, 
	grid_t& grid);

int MPI_CalculateIllumAsync(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell, grid_t& grid);
#endif

#endif //ILLUM
#endif //ILLUM_UTILS_H