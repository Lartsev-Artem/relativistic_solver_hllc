#include "../../prj_config.h"

#include "../../global_headers.h"
#include "../../global_def.h"
#include "../../global_value.h"


static size_t SetBasis(const Type* start_point, Vector3& normal, Matrix3& basis) {
	/*�� ��������� ����� � ������� ������ ��������� ����� ��������� ��������� (vec1, vec2).
	  ������� ����. ������ ���� ������ �����������(������������ normal). ������ ������ �� ���������� ������������*/
	Vector3 vec_1;
	Vector3 vec_2;

	if (abs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //��-�� ���������� ������������ (N, vec1)==0
	}

	// ���������� ���������� ������ ���������
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// ������� ��������� ���������. Eigen �������� � �� �����!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
static size_t Make2dPoint(const Type* start, const Matrix3& local_basis, const Type* point, Vector3& new_point) {

	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//������� 3d ����� � 2d (� ��������� ������ {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}

size_t IntersectionWithPlaneDisk(const Vector3& X0, const Vector3& n, Vector3& res) {

	//  ----------������ ������. �.�. ���� �������� ���������� ����������, ��������� ����� ������ ����--------------
	/*
	 {
	std::vector<Vector3> curface(3);		 // ����� �������� ��������� �����
			curface[0][0] = 1;
			curface[0][1] = 0;
			curface[0][2] = 0;

			curface[1][0] = 0;//0;
			curface[1][1] = 0.9928768384869221;//
			curface[1][2] = 0.11914522061843064;//;

			curface[2][0] = 2;//;
			curface[2][1] = 0;
			curface[2][2] = 0;// ;   // Wolfram
			}

	* const std::vector<Vector3>& face,
	Vector3 A = face[0];
	Vector3 B = face[1];
	Vector3 C = face[2];

	Type a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	Type b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	Type c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	Type d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	Type t = -(a * X0[0] + b * X0[1] + c * X0[2] + d) / (a * n[0] + b * n[1] + c * n[2]);
	*/

	/*
	a= 0
	b= 0.1191452206184306
	c= -0.9928768384869221
	d= 0
	*/

	const Type b = 0.1191452206184306;
	const Type c = -0.9928768384869221;

	const Type t = -(b * X0[1] + c * X0[2]) / (b * n[1] + c * n[2]);

	res = (t * n + X0);

	/*for (size_t i = 0; i < 3; i++)
		res[i] = (n[i] * t + X0[i]);*/

	return 0;
}
int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//������� ������������

	Type a, b, c, d;  // ��������� ��������� ���������
	Type t;

	a = face.A[1] * (face.B[2] - face.C[2]) + face.B[1] * (face.C[2] - face.A[2]) + face.C[1] * (face.A[2] - face.B[2]);
	b = face.A[0] * (face.C[2] - face.B[2]) + face.B[0] * (face.A[2] - face.C[2]) + face.C[0] * (face.B[2] - face.A[2]);
	c = face.A[0] * (face.B[1] - face.C[1]) + face.B[0] * (face.C[1] - face.A[1]) + face.C[0] * (face.A[1] - face.B[1]);
	d = face.A[0] * (face.C[1] * face.B[2] - face.B[1] * face.C[2]) + face.B[0] * (face.A[1] * face.C[2] -
		face.C[1] * face.A[2]) + face.C[0] * (face.B[1] * face.A[2] - face.A[1] * face.B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // ����� ����������� ����  (start->direction) � ����������!!! face

	return 0;
}
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Vector3& XX)
{
	/*face --- �����������, X --- ����� ��� ��������*/

	// ������� ������������
	const Type* AA = cell_face.A.data();
	const Type* BB = cell_face.B.data();
	const Type* CC = cell_face.C.data();

	Vector3 A, B, C, X;  // ����� ����� �� ���������
	{
		Eigen::Matrix3d basis;
		Vector3 n = normals_cell.n[number_face % 4];
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, XX.data(), X);
	}

	// �������� �������
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face) {
	// 3 ���� ������������
		{
			straight_face << 1. / 6, 1. / 6, 1,
				2. / 3, 1. / 6, 1,
				1. / 6, 2. / 3, 1;
		}

		// 3 ���� ������������ �� ��������� ���������
		{
			inclined_face <<
				0, sqrt(2. / 3), 1,
				sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
				-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
		}

		//������� �������� �� ������������ ��������� � ���������� ��������� ��������� 
		{ transform_matrix <<
			-1. / sqrt(2), 1. / sqrt(2), 0,
			-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
			1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
		}

		//������� �������� �� ��������� ��������� �  ���������� ������������ ���������
		{
			inverse_transform_matrix <<
				-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				0, sqrt(2. / 3), 1. / sqrt(3);
		}

		// ������ ���������� ���������
		start_point_plane_coord << 0.5, 0.5, 0;
		return 0;
}

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face, Matrix3& inclined_face_inverse) 
{

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

	inclined_face_inverse = inclined_face.inverse();  // � �������

	//straight_face_inverse = straight_face.inverse(); // � �������	

	return 0;
}


int GetInAndOutFaces(const Vector3& direction, const Normals& normals, int& face_state) 
{
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	for (size_t i = 0; i < base; ++i) 
	{
		Type sca = normals.n[i].dot(direction);

		if (normals.n[i].dot(direction) < -eps)
			face_state = set_bit(face_state, i);
		else
			face_state = clear_bit(face_state, i); // ���� ����� �����������, �������, ��� ��� �� �������� ������������
	}

	return 0;
}


void MakeRotationMatrix(const Vector3& n, MatrixX& T) {

	T = MatrixX::Zero(5, 5);
	T(0, 0) = T(4, 4) = 1;

	if (fabs(n[2] * n[2] - 1) > eps)
	{

		T(1, 1) = n[0];
		T(1, 2) = n[1];
		T(1, 3) = n[2];

		Type sqr = sqrt(1 - n[2] * n[2]);

		if (sqr < eps * eps - eps / 10)
			printf("Err T\n");

		T(2, 1) = -n[1] / sqr;
		T(2, 2) = n[0] / sqr;

		T(3, 1) = -n[0] * n[2] / sqr;
		T(3, 2) = -n[1] * n[2] / sqr;
		T(3, 3) = sqr;
	}
	else if (n[2] > 0)  // n_z == 1
	{
		T(1, 3) = 1;
		T(2, 2) = 1;
		T(3, 1) = -1;
	}
	else  // n_z == -1
	{
		T(1, 3) = -1;
		T(2, 2) = -1;
		T(3, 1) = 1;
	}

}
void MakeRotationMatrix(const Vector3& n, Matrix4& T) {

	// n=(x,y,0)!!!!

	T = Matrix4::Zero();
	T(0, 0) = T(3, 3) = 1;

	//T(1, 1) = n[0];
	//T(1, 2) = n[1];

	//T(2, 1) = -n[1];
	//T(2, 2) = n[0];

	T(1, 1) = -n[0];
	T(1, 2) = -n[1];

	T(2, 1) = n[1];
	T(2, 2) = -n[0];

}
int MakeRotationMatrix(const Vector3& n, Matrix3& T)
{
	T = Matrix3::Zero();

	if (fabs(n[2] * n[2] - 1) > eps)
	{

		T(0, 0) = n[0];
		T(0, 1) = n[1];
		T(0, 2) = n[2];

		Type sqr = sqrt(1 - n[2] * n[2]);

		if (sqr < eps * eps - eps / 10)
			printf("Err T\n");

		T(1, 0) = -n[1] / sqr;
		T(1, 1) = n[0] / sqr;

		T(2, 0) = -n[0] * n[2] / sqr;
		T(2, 1) = -n[1] * n[2] / sqr;
		T(2, 2) = sqr;
	}
	else if (n[2] > 0)  // n_z == 1
	{
		T(0, 2) = 1;
		T(1, 1) = 1;
		T(2, 0) = -1;
	}
	else  // n_z == -1
	{
		T(0, 2) = -1;
		T(1, 1) = -1;
		T(2, 0) = 1;
	}

	return 0;
}

#if 0 // �������� �� ������ ������ �������
static inline void MakeRotationMatrix(const Vector3& n, Matrix3& T)
{
	const Type x = n[0];
	const Type y = n[1];
	const Type theta = atan2(y, x);
	const Type phi = atan2(sqrt(x * x + y * y), n[2]);

	T(0, 0) = cos(theta) * sin(phi);
	T(0, 1) = sin(theta) * sin(phi);
	T(0, 2) = cos(phi);

	T(1, 0) = -sin(theta);
	T(1, 1) = cos(theta);
	T(1, 2) = 0;

	T(2, 0) = -cos(theta) * cos(phi);
	T(2, 1) = -sin(theta) * cos(phi);
	T(2, 2) = sin(phi);

	return;
}
#endif