#ifndef ILLUM_UTILS_H
#define ILLUM_UTILS_H

#include "../../global_def.h"
#include "../../global_headers.h"
#include "../solve_config.h"
#include "../solve_global_struct.h"
#if defined ILLUM && defined SOLVE

int GetDirectionIllumFromFace(const int size_grid, const int num_dir, const std::vector<Type>& illum_on_face, std::vector<Type>& illum_in_cell);

int InitIllum(file_name main_dir, grid_t& grid);

int CalculateIllum(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x, const std::vector < std::vector<int>>& sorted_id_cell,
	//const std::vector<Type>& res_inner_bound, 
	grid_t& grid, std::vector<Type>& Illum, std::vector<Type>& int_scattering);

int SolveIllumAndHLLC(const Type tau, std::vector<elem_t>& cells);
int CalculateIllumParam(const grid_directions_t& grid_direction, grid_t& grid);

#endif //ILLUM
#endif //ILLUM_UTILS_H