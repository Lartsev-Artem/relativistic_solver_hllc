#pragma once
#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "solve_short_characteristics_global_structure.h"
#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_calculations.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_direction, const int num_in_face, const Vector3& x, const std::vector<Vector2>& X0,
	 std::vector<cell>& grid, const std::vector<int>& neighbours_id_face,
	int& id_try_pos, int& pos_in_res, int& posX0,
	const uint64_t ShiftRes, const uint64_t ShiftX0, const int ShiftTry);

Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	const vector<Type>& int_scattering);
Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	const vector<Type>& int_scattering, const vector<VectorX>& U);
#ifdef  USE_VTK
Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares);
#endif //  USE_VTK
#endif
