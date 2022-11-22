#pragma once
#ifndef SHORT_CHARACTERISTICS_HLLC_H
#define SHORT_CHARACTERISTICS_HLLC_H

#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_global_structure.h"

int ReBuildDataForHLLC(const int N, std::vector<VectorX>& data);
inline void MakeRotationMatrix(const Vector3& n, Eigen::MatrixXd& T);
int HLLC_step(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Eigen::VectorXd>& U_full, std::vector<Eigen::VectorXd>& U_full_prev);

VectorX HLLC_stepToOMP(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev);

int HLLC(const int N, const Type tau, const std::vector<int>& bound_cells, const std::vector<int>& neighbours_id_faces,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	const std::vector<VectorX>& U_full_prev);
#ifdef RHLLC
void HLLC_Rel(const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, std::vector<VectorX>& U_full_prev);

Type FormTimeStepToRHLLC(const int n, const Type h, const Type k);
int ReBuildDataForHLLCRel(const int N, std::vector<VectorX>& data);
#endif

#ifdef HLLC_1D
int HLLC_1d(std::string& main_dir);
#endif

#ifdef RHLLC_1D
int RHLLC_1d(std::string& main_dir);
#endif


#endif