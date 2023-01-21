#ifndef SHORT_CHARACTERISTICS_CALCULATIONS_H
#define SHORT_CHARACTERISTICS_CALCULATIONS_H
#include "struct_short_characteristics_global_structure.h"

#ifdef MAKE

#ifdef  USE_VTK
int WriteCellFaces(const std::string name_file_cells, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);
int WriteVertex(const std::string name_file_vertex, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);
#endif //  USE_VTK

int ReadCellFaces(const std::string name_file_cells, std::vector<Face>& grid);
int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord);
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord);

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord);
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord);

inline Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);
inline Vector3 GetInterpolationCoefInverse(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);

#endif
#endif //MAKE