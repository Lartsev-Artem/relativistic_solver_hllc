#include "../file_module/writer_bin.h"
#include "../utils/grid_geometry/geometry_data.h"

#include "../global_value.h"

#ifdef MAKE

#ifdef USE_VTK

int WriteCellFaces(const std::string name_file_cells, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	FILE* f;
	f = fopen(name_file_cells.c_str(), "wb");

	//Type A[3];
	//Type B[3];
	//Type C[3];
	vtkPoints* points_face;
	std::vector<Type>pp(9);

	int n = unstructured_grid->GetNumberOfCells();
	fwrite_unlocked(&n, sizeof(int), 1, f);

	//std::vector<Vector3>centers_face(n * 4);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < 4; j++) {
			points_face = unstructured_grid->GetCell(i)->GetFace(j)->GetPoints();

			points_face->GetPoint(0, pp.data());
			points_face->GetPoint(1, pp.data() + 3);
			points_face->GetPoint(2, pp.data() + 6);

			//centers_face[i*4+j] = Vector3((pp[0]+pp[3]+pp[6])/3, (pp[1] + pp[4] + pp[7]) / 3, (pp[2] + pp[5] + pp[8]) / 3 );

			fwrite_unlocked(pp.data(), sizeof(Type), 9, f);

		}
	}

	fclose(f);
	return 0;
}
int WriteVertex(const std::string name_file_vertex, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	const int n = unstructured_grid->GetNumberOfCells();
	std::vector<Matrix4> vertexs(n);
			
	for (size_t i = 0; i < n; i++) 
	{
		SetVertexMatrix(i, unstructured_grid, vertexs[i]);		
	}

	WriteSimpleFileBin(name_file_vertex, vertexs);

	return 0;
}

#endif//USE_VTK

#ifndef ONLY_GEO_DATA

int ReadCellFaces(const std::string name_file_cells, std::vector<Face>& grid) {

	FILE* f;
	f = fopen(name_file_cells.c_str(), "rb");

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	grid.resize(n * 4);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < 4; j++) {
			fread_unlocked(grid[i * 4 + j].A.data(), sizeof(Type), 3, f);
			fread_unlocked(grid[i * 4 + j].B.data(), sizeof(Type), 3, f);
			fread_unlocked(grid[i * 4 + j].C.data(), sizeof(Type), 3, f);
		}
	}

	fclose(f);

	return 0;
}

int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord)
{
	// vertex_tetra -> [X;Y;Z;1]
	// возможно надо будет использовать transpose из-за инициализации матрицы перехода

	Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
	// сразу с транспонированием
	for (int i = 0; i < 3; ++i)
	{
		local_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			local_coord[i] += vertex_tetra_inverse(i, j) * global_coord[j];
		local_coord[i] += vertex_tetra_inverse(i, 3);
	}

	//local_coord = vertex_tetra * global_coord;
	return 0;
}
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord) {
	// vertex_tetra -> [X,Y,Z,1]
	// написать в ручную т.к. преобразования известны, и 4я строка постоянна и меняться не должна

	Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
	for (int i = 0; i < 3; ++i)
	{
		global_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			global_coord[i] += vertex_tetra(i, j) * local_coord[j];
		global_coord[i] += vertex_tetra(i, 3) * eta4;
	}

	//global_coord = vertex_tetra.inverse() * local_coord;

	return 0;
}

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}


Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})

	return interpolation_nodes.partialPivLu().solve(function_value);

}
Vector3 GetInterpolationCoefInverse(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})

	return interpolation_nodes * function_value;
}

#endif //ONLY_GEO_DATA

#include "struct_short_characteristics_calculations.h"
BaseTetra_t::BaseTetra_t()
{
	// 3 узла интерполяции
	{
		straight_face << 1. / 6, 1. / 6, 1,
			2. / 3, 1. / 6, 1,
			1. / 6, 2. / 3, 1;
	}

	// 3 узла интерполяции на наклонной плоскости
	{
		inclined_face <<
			0, sqrt(2. / 3), 1,
			sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
			-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
	}

	//Матрицы перехода из стандартного тетраэдра в координаты наклонной плоскости 
	{ transform_matrix <<
		-1. / sqrt(2), 1. / sqrt(2), 0,
		-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
		1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
	}

	//Матрицы перехода из наклонной плоскости в  координаты стандартного тетраэдра
	{
		inverse_transform_matrix <<
			-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
			1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
			0, sqrt(2. / 3), 1. / sqrt(3);
	}

	// Начало координата плоскости
	start_point_plane_coord << 0.5, 0.5, 0;

	inclined_face_inverse = inclined_face.inverse();  // в решении
	straight_face_inverse = straight_face.inverse(); // в решении	
}

#endif //MAKE