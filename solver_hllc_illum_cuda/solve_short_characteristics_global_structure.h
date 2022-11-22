#pragma once
#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "solve_short_characteristics_headers.h"

typedef uint8_t ShortId;
typedef double Type;
typedef std::string Str_Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;



#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

#define PRINTF(a) {printf(a); return 1;}

using namespace std::chrono;
using namespace std;

struct Normals {
	std::vector<Vector3> n;
	Normals() {
	}

	Normals(const int size) {
		n.resize(size);
	}
};

struct cell {
	//int id;  - номер в массиве
	std::vector<Vector3> nodes_value;
	//std::vector<int> neighbours_id_face;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3::Zero());//Vector3(-666, -666, -666));	
	}
};

enum eBoundaryTypes
{
	eBound_FreeBound = -1,
	eBound_LockBound = -2,
	eBound_OutSource = -3,
	eBound_InnerSource = -4
};


#define PI 3.14159265358979323846
const double eps = 1e-14;
const int count_threads = 8;

extern Vector3 start_point_plane_coord;   // начало координат плоскости
extern Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
extern Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

extern Matrix3	straight_face;  // 3 узла интерпол€ции
extern Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

extern Matrix3	straight_face_inverse;  // 3 узла интерпол€ции
extern Matrix3 inclined_face_inverse;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

extern size_t class_file_vtk;
extern int size_grid;

// скал€рные данные сетки (unstructured_grid)
#ifdef USE_VTK
extern vtkDataArray* density;
extern vtkDataArray* pressure;
extern vtkDataArray* velocity;
extern vtkDataArray* absorp_coef;
extern vtkDataArray* rad_en_loose_rate;

extern vtkDataArray* e_substance;
#else
extern std::vector<Type> density;
extern std::vector<Type> absorp_coef;
extern std::vector<Type> rad_en_loose_rate;
extern std::vector<Type> pressure;
extern std::vector<Vector3> velocity;
extern std::vector<Type> e_substance;  // не хранить, а пересчитывать?

extern std::vector<VectorX> U_full;
#endif
extern Type square_surface;  // площадь поверхности дискретной 

extern Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

extern std::vector<Type> res_inner_bound;  // значение на внутренней границе
extern std::vector<int>id_try_surface;

const Type res_on_inner_bound = 5;

const double gamma1 = 5. / 3; // 5. / 3;  // показатель адиабаты
const double gamma_g = gamma1 / (gamma1 - 1);


#ifdef SAVE_DUMP_HLLC
extern uint32_t size_dump;
extern uint32_t w_ptr;
extern std::vector<double> dump;
extern int start_section;
extern bool bad_hllc_flag;
#endif

#endif