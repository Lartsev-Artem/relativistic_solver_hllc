#ifndef GEOMETRY_DATA
#define GEOMETRY_DATA

#include "../../prj_config.h"
#include "../../global_headers.h"

#if (defined USE_VTK)

int GetNeighborFace3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);

int GetNeighborFace2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);

int GetNormalAndSquares3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Normals>& normals, std::vector<Type>& area);
int GetNormalAndSquares2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Normals>& normals, std::vector<Type>& area);

int GetVolume3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Type>& volume);
int GetVolume2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Type>& volume);


int GetCentersOfTetra3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers);
int GetCentersOfTetra2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers);

int GetCentersOfFaces3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers);

int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra);

#endif //USE_VTK

int SetTypeOfBound(const std::vector<Vector3>& centers, const std::vector<Normals>& normals, std::vector<int>& all_pairs_face);
int ReWriteGeoFiles(file_name name_file_geometry_faces, file_name name_file_geometry_cells);

#endif // !GEOMETRY_DATA
