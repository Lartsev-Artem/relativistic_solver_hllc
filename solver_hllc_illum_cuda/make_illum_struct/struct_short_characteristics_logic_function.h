#pragma once
#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_calculations.h"


#ifndef ReadNow

int GetNodes(const int num_cur_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<cell>& nodes_value,
	std::vector<Type>& vec_res_bound, std::vector<Type>& vec_s, std::vector<Vector3>& vec_x, std::vector<Vector2>& vec_x0_local, std::vector<ShortId>& vec_in_id);
int CalculateNodeValue(const int num_cur_cell, const std::vector<Normals>& normals, const std::vector<Face>& grid,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	std::vector<Type>& vec_res_bound, std::vector<Type>& vec_s, std::vector<Vector2>& vec_x0_local, std::vector<ShortId>& vec_in_id);
Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value,
	std::vector<Type>& vec_res_bound, std::vector<Vector2>& vec_x0_local);
#else
int GetNodes(const int num_cur_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<cell>& nodes_value,
	FILE* file_res_bound, FILE* file_s, FILE* file_x, FILE* file_x0_local, FILE* file_in_id);

int CalculateNodeValue(const int num_cur_cell, const std::vector<Normals>& normals, const std::vector<Face>& grid,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	FILE* file_res_bound, FILE* file_s, FILE* file_x0_local, FILE* file_in_id);

int CalculateIllumeOnInnerFace(const int num_cur_cell, const ShortId num_in_face, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	Vector3& x0, std::vector<cell>& nodes_value, FILE* file_res_bound, FILE* file_x0_local);
#endif

#endif
