#ifndef GLOBAL_VALUE
#define GLOBAL_VALUE

#include <string>
#include "prj_config.h"
#include "global_def.h"

extern std::string BASE_ADRESS;


#define PI 3.14159265358979323846


extern std::vector<Vector3> X;
extern std::vector<Vector2> X0;
extern std::vector<Type> S;

const double eps = 1e-10;

extern Vector3 start_point_plane_coord;   // ������ ��������� ���������
extern Matrix3 transform_matrix;          // ������� �������� �� �������� ��������� � ���������
extern Matrix3 inverse_transform_matrix;  // ������� �������� �� ��������� � ������� ��������

extern Matrix3	straight_face;  // 3 ���� ������������
extern Matrix3 inclined_face;  // 3 ���� ������������ �� ��������� ���������

//std::vector<Type> Illum2;


#ifdef Sphere
const Vector3 center_point(0, 0, 0);
const Type inner_radius = 0.51; // ������ ���������� ����� (� �������)
#else 
const Vector3 center_point(10, 0, 0);
const Type inner_radius = 0.12; // ������ ���������� ����� (� �������)
#endif

#endif //GLOBAL_VALUE
