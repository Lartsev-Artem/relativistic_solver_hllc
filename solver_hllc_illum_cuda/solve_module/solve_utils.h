#ifndef SOLVE_UTILS_H
#define SOLVE_UTILS_H
#include "solve_config.h"

#ifdef SOLVE

#include "solve_global_struct.h"


void ReBuildNeighStruct(
	std::vector<int>& neighbours_id_faces,
	std::vector<Normals>& normals,
	std::vector<Type>& squares_faces,
	std::vector<Type>& volume,
	std::vector<Vector3>& centers,
	std::vector<face_t>& faces, std::vector<elem_t>& cells);

int ReWriteGeoFiles(file_name name_file_geometry_faces, file_name name_file_geometry_cells);

int HLLC_STEP(const Type tau, grid_t& grid);
int GetTimeStep(hllc_value_t& hllc_cfg, const grid_t& grid);
int HLLC_INIT(file_name file_settings_hllc, hllc_value_t& hllc_set, file_name file_init_value, std::vector<elem_t>& cells);


int StartLowDimensionTask(file_name main_dir);

int TestDivStream(const std::vector<Vector3>& centers_face, grid_t& grid);
int TestDivStream(file_name BASE_ADRESS);
#endif //SOLVE
#endif //SOLVE_UTILS_H