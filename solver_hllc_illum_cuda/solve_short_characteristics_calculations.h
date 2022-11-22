#pragma once
#ifndef SHORT_CHARACTERISTICS_CALCULATIONS_H
#define SHORT_CHARACTERISTICS_CALCULATIONS_H

#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_global_structure.h"

#ifdef USE_VTK
size_t ReWriteDataArray(const size_t class_file_vtk, const std::string name_file_vtk, const std::string& main_dir,
	int& size_grid,	const bool is_print = false);

size_t ReadDataArray(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, int& size_grid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print = false);

size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print = false);

int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);
int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell);
int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra);
size_t CenterOfTetra(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra);
size_t FindAllCenterOfTetra(const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, std::vector<Vector3>& point_in_tetra);
bool InTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX);
size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result);
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	vtkSmartPointer<vtkUnstructuredGrid>& u_grid);

size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls,
	vtkSmartPointer<vtkUnstructuredGrid>& u_grid);
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::string& file_vtk);
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::string& file_vtk);

size_t ReBuildDataArray(const size_t class_file_vtk, const std::string name_file_vtk, const std::string name_file_out,
	const std::string& main_dir, const bool is_print=false);

//size_t ReBuildDataArrayFull(const size_t class_file_vtk, const std::string name_file_vtk, const std::string name_file_out,
	//const std::string& main_dir, const bool is_print=false);

#else
size_t ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	int& size_grid, const bool is_print = false);

size_t ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	std::vector<Vector3>& velocity, std::vector<Type>& pressure,
	int& size_grid, const bool is_print=false);

size_t WriteFileSolution(const std::string& main_dir, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::vector<Type>& vector_div_stream);
size_t WriteFileSolution(const std::string& main_dir, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::vector<Type>& vector_div_stream,
	const std::vector<Vector3>& vector_div_impuls, const std::vector<VectorX>& vector_U);


#endif

size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Vector3>& directions_all, 
	vector<Type>& squares, Type& square_surface);

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face, Matrix3& inclined_face_inverse, Matrix3& straight_face_inverse);


int InitNodesValue(const std::vector<int>& all_pairs_face, std::vector<cell>& nodes_value);

int ResetNodesValue(std::vector<cell>& nodes_value);
int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals);

Type NormIllum(const std::vector<Type>& Illum, const std::vector<Type>& Illum2);
Type NormIllumOmp(const std::vector<Type>& Illum, const std::vector<Type>& Illum2);
int ReadGraph(const std::string name_file_graph, std::vector<int>& sorted_id_cell);
int ReadGraphBin(const std::string name_file_graph, std::vector<int>& sorted_id_cell);
size_t SetBasis(const Type* start_point, Vector3& normal, Matrix3& basis);
size_t Make2dPoint(const Type* start, const Matrix3& local_basis, const Type* point, Vector3& new_point);

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state);

int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord);
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord);

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord);
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord);

Type BoundaryFunction(const int id_cell, const Vector3& x, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);

 Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);
 Vector3 GetInterpolationCoefInverse(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);

Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares);
int CalculateInt(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering);
int CalculateIntOmp(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering);
int Min(const int a, const int b);

int MakeStream(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface, vector<Vector3>& stream);
int MakeImpuls(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface, vector<Matrix3>& impuls);
size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy);
size_t MakeEnergy(const Type* Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy);
int MakeDivStream(const vector<Vector3>& stream, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, std::vector<int>& neighbours_id_face, vector<Type>& div_stream);

Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface);

int ReadSizes(const std::string& name_file_size, int& countX, int& countX0, int& countOutC,
	int& countOut, int& countIn, int& countS, int& countRes, int& countTry);

size_t ReadCompactFastGridData(const int count_dir, const int N, const Str_Type& file_logs,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	Str_Type& name_file_shift_out, Str_Type& name_file_shift_res, Str_Type& name_file_shift_x0, Str_Type& name_file_shift_try,
	std::vector<cell>& grid, std::vector<int>& neighbours_id_face, std::vector<ShortId>& OutC,
	std::vector<ShortId>& Out, std::vector<ShortId>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<int>& sorted_id_cell,
	vector<uint64_t>& ShiftOut, vector<uint64_t>& ShiftRes, vector<uint64_t>& ShiftX0, vector<int>& ShiftTry);

size_t ReadCompactFastGridDataOptMemory(const int count_dir, const int N, const Str_Type& main_dir,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	Str_Type& name_file_shift_out, Str_Type& name_file_shift_res, Str_Type& name_file_shift_x0, Str_Type& name_file_shift_try,
	std::vector<cell>& grid, std::vector<int>& neighbours_id_face, std::vector<ShortId>& OutC,
	std::vector<ShortId>& Out, std::vector<ShortId>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<int>& sorted_id_cell,
	vector<int>& ShiftOut, vector<int>& ShiftRes, vector<int>& ShiftX0, vector<int>& ShiftTry);

int ReadCentersOfTetra(const std::string name_file_centers, std::vector<Vector3>& centers);

#endif