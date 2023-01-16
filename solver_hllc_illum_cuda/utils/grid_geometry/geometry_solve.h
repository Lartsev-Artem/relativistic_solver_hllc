#ifndef GEOMETRY_SOLVE
#define GEOMETRY_SOLVE

#include "../../prj_config.h"
#include "../../global_headers.h"

size_t IntersectionWithPlaneDisk(const Vector3& X0, const Vector3& n, Vector3& res);
int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result);
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Vector3& XX);

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix, Matrix3& straight_face, Matrix3& inclined_face);

int GetInAndOutFaces(const Vector3& direction, const Normals& normals, int& face_state);
#endif // !GEOMETRY_SOLVE
