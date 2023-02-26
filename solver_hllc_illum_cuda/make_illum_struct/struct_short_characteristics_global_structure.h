#if !defined SHORT_CHARACTERISTICS_GLOBAL_H && defined MAKE
#define SHORT_CHARACTERISTICS_GLOBAL_H


#include "../global_def.h"
#include "../global_headers.h"
#include "../global_value.h"

struct BaseTetra_t
{
	Vector3 start_point_plane_coord;   // ������ ��������� ���������
	Matrix3 transform_matrix;          // ������� �������� �� �������� ��������� � ���������
	Matrix3 inverse_transform_matrix;  // ������� �������� �� ��������� � ������� ��������

	Matrix3	straight_face;  // 3 ���� ������������
	Matrix3 inclined_face;  // 3 ���� ������������ �� ��������� ���������

	Matrix3	straight_face_inverse;  // 3 ���� ������������
	Matrix3 inclined_face_inverse;  // 3 ���� ������������ �� ��������� ���������

	//	Vector3 center_local_sphere;  // ����� ��������� ����� ����� ������������ ���������

	BaseTetra_t();
};
const BaseTetra_t base_tetra_geo;

// ��������� ����� � ���������� �����:
const Type Rsphere = 0.001;
const Type R1disk = 0.001;
const Type R2disk = 0.09;

extern std::vector<Vector3> x_try_surface;
extern std::vector<int> id_try_surface;

#endif