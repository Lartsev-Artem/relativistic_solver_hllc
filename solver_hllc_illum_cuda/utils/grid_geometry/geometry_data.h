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

#endif //USE_VTK

#if 1 //defined BUILD
int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result);
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Vector3& XX);
#endif

int SetTypeOfBound(const std::vector<Vector3>& centers, const std::vector<Normals>& normals, std::vector<int>& all_pairs_face);

#endif // !GEOMETRY_DATA
