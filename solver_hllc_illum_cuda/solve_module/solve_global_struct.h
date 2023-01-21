#ifndef SOVLE_GLOBAL_STRUCT_H
#define SOVLE_GLOBAL_STRUCT_H

#include "solve_config.h"

struct hllc_value_t
{		
	Type T;
	Type CFL;	
	Type h;
	Type print_timer;
	
	Type tau;
};

struct solve_mode_t
{
	int size_grid;
	int class_vtk;
	int max_number_of_iter;
	Type accuracy;
	bool use_cuda;
	int cuda_mod; //1 - все массивы, 0 - min
	std::string name_file_solve;
	
	solve_mode_t()
	{
		size_grid = 0;
		class_vtk = 0;
		max_number_of_iter = 1;
		accuracy = 1e-5;
		use_cuda = true;
		cuda_mod = 1; //1 - все массивы, 0 - min
		name_file_solve = "Solve";
	}
};

struct flux_t
{
public:
	Type d;
	Vector3 v;
	Type p;

	flux_t();
	flux_t(const Type a, const Type b, const Type c, const Type dd, const Type e);
	
	flux_t operator+ (const flux_t& x);
	flux_t operator+= (const flux_t& x);
	flux_t operator-= (const flux_t& x);
	flux_t operator* (const Type x);
	flux_t operator- (const flux_t& x);
	flux_t operator/ (const Type x);

	// это временно дл€ свз€и со старым кодом
	Type operator[](const int i);	
	Type operator()(const int i);

//private:
	flux_t(const flux_t& f);
};

struct geo_face_t
{
	int id_l;
	int id_r;

	//int bound; // признак границы(наверное удобнее хранить его отдельно, но id_r включает его)

	Vector3 n;
	Type S;

	geo_face_t();
};
struct face_t
{
public:
	flux_t  f;
	geo_face_t geo; // геометри€ €чейки
	face_t() {};
};

struct geo_cell_t
{
	int id_faces[base];
	Type V;
	bool sign_n[base];
	Vector3 center;

	geo_cell_t();
};
struct illum_value_t
{
	//std::vector<Type> int_scattering; // (count_cells* count_directions, 0);

	std::vector<Type> illum; //num_dir*base
	Type energy;
	Vector3 stream;

	//Type prev_energy;
	//Vector3 prev_stream;

	Matrix3 impuls;
	Type div_stream;
	Vector3 div_impuls;

	// сейчас эти переменные лежат отдельно(они нужны дл€ статического расчЄта излучени€)

	//Type density;
	Type absorp_coef;
	Type rad_en_loose_rate;

	Vector3 coef_inter[base]; // коэффициенты интерпол€ции 

	illum_value_t(const int num_dir = 0);
};
struct elem_t
{
	flux_t  phys_val;

#if defined HLLC || defined RHLLC
	flux_t  conv_val;	
#endif

#ifdef ILLUM	
	illum_value_t illum_val;	
#endif

	geo_cell_t geo; //геометри€ элемента
	
//private:
//	elem_t(const elem_t& el) {};
};

struct grid_t
{
	int size;
	std::vector<elem_t> cells;
	std::vector<face_t> faces;

	grid_t() { size = 0; }
};

extern solve_mode_t solve_mode;
#endif //SOVLE_GLOBAL_STRUCT_H