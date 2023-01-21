#ifndef SHORT_CHARACTERISTICS_LOGIC_H
#define SHORT_CHARACTERISTICS_LOGIC_H

#include "../global_headers.h"

#ifndef ONLY_GEO_DATA

int GetNodes(const int num_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int face_state, const Vector3& direction, const std::vector<Normals>& normals,
	const std::vector<int>& all_pairs_face, const BasePointTetra& x,
	std::vector<Type>& vec_res_bound, std::vector<cell_local>& vec_x0);

int MakeArrayX(const std::vector<Eigen::Matrix4d>& vertexs, std::vector<BasePointTetra>& vec_x);

#endif //ONLY_GEO_DATA
#endif //SHORT_CHARACTERISTICS_LOGIC_H