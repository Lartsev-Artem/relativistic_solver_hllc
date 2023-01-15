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

	~Normals() 
	{
		n.clear();
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

#ifdef USE_MPI
#include "mpi.h"
#define RETURN(a) return a;
#define MPI_END MPI_Finalize();
#define MPI_RETURN(a) { MPI_END RETURN(a) }
#else
#define MPI_RETURN(a) return a;
#endif //USE_MPI


#if 1
//#ifdef NEW_CLASS

#define base 4
struct illum_value
{
	std::vector<Type> int_scattering; // (count_cells* count_directions, 0);

	Type energy;
	Vector3 stream;

	//Type prev_energy;
	//Vector3 prev_stream;

	Matrix3 impuls;
	Type div_stream;
	Vector3 div_impuls;

	illum_value(const int num_dir = 0);
};
struct flux
{
	Type d;
	Vector3 v;
	Type p;

	flux();
	flux(const Type a, const Type b, const Type c, const Type dd, const Type e);
	
	flux operator+ (const flux& x); 
	flux operator+= (const flux& x);
	flux operator-= (const flux& x);
	flux operator* (const Type x); 
	flux operator- (const flux& x);
	flux operator/ (const Type x); 

	// это временно дл€ свз€и со старым кодом
	Type operator[](const int i);	
	Type operator()(const int i);

	//private:
	//	flux(const flux& f) {};
};
struct face
{
	flux f;
#ifndef ONLY_HLLC	
	std::vector<Type> Illum;	
	Vector3 x[base][base - 1];	   //[номер грани][номер точки интерпол€ции]	
	Vector2 x0[base][base - 1];	   // [номер грани][номер точки интерпол€ции]	
	
#endif

	int id_l;
	int id_r;

	Vector3 n;
	Type S;

	face(const int num_dir = 0);	
};

struct elem
{
#ifndef ONLY_ILLUM
	flux val;
	flux phys_val;
#endif

#ifndef ONLY_HLLC	
	illum_value illum_val;	

	std::vector<uint8_t> enter_face; //вход€ща€ или выход€ща€ грань {0-вход€ща€ наклонна€, 1-вход€ща€ пр€ма€, 2- выход€ща€ наклонна€,3-выход€ща€ пр€ма€}

#if 1
	std::vector<uint8_t> inner_face; //номер св€занной вход€щей грани
	std::vector<Type> s_line; //(size= вычисл€етс€ в make, разный по всем направлени€м [])
#else
	std::vector<std::pair<uint8_t, Type>> face_and_s;
#endif

#endif

	int id_faces[base];
	Type V;
	bool sign_n[base];		
	Vector3 center;

	elem(const int num_dir = 0);
private:
	elem(const elem& el) {};
};
#endif

#endif