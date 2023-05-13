#ifndef GLOBAL_DEF
#define GLOBAL_DEF

#include "prj_config.h"
#include "global_headers.h"

#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__

#define PRINT_POS  printf("Calls from: %s: %s, %d c.\n",__FILE__, __FUNCTION__,__LINE__);

#ifdef USE_CUDA
#define CUDA_ERR(str){ WRITE_LOG_ERR(str); solve_mode.use_cuda = false;}
#else
#define CUDA_ERR(str){}
#endif //USE_CUDA

#define EXIT_ERR(str){printf(str); EXIT(1);}
#define EXIT_ERRS(str,val){printf(str,val); EXIT(1);}

#define RETURN_ERRS(str, val) { printf(str,val); PRINT_POS return 1; }
#define RETURN_ERR(str) { printf(str); PRINT_POS return 1; }

#define base (NUMBER_OF_MEASUREMENTS + 1)

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Vector4d Vector4;
typedef Eigen::Matrix4d Matrix4;
typedef Eigen::MatrixXd MatrixX;

typedef double Type;
typedef int IntId;
typedef uint8_t State;

typedef uint8_t ShortId;
typedef std::string Str_Type;

typedef const std::string& file_name;

using namespace std;
using namespace std::chrono;

enum eBoundaryTypes
{
	eBound_FreeBound = -1,
	eBound_LockBound = -2,
	eBound_OutSource = -3,
	eBound_InnerSource = -4
};

struct Normals {
	std::vector<Vector3> n;
	Normals(const int size = base) { n.resize(size); }
	~Normals() { n.clear(); }
};

struct Face {
	Vector3 A;
	Vector3 B;
	Vector3 C;
	Face& operator=(const Face& face) 
	{
		A = face.A;
		B = face.B;
		C = face.C;
		return *this;
	}
};

struct FaceCell {
	int face_id;
	Face face;
	FaceCell(const int id = 0, const Face& face_init = Face()) 
	: face_id(id), face(face_init) {}
};

struct direction_s
{
	Vector3 dir;
	Type area;
};

struct grid_directions_t
{
	int size;
	std::vector<direction_s> directions;
	Type full_area;
};

struct BasePointTetra //узлы интерпол€ции всех тетраэдров // ¬ перспективе можно уйти к гран€м
{
	Vector3 x[base][NUMBER_OF_MEASUREMENTS];

	Vector3 operator()(const int i, const int j) const
	{
		return x[i][j];
	}
};
struct cell_local // дл€ каждой €чейки и каждого направлени€
{
	Type s;    //рассто€ние x0-x
	Vector2 x0; //локальна€ координата входного узла дл€ интерпол€ции
	ShortId in_face_id; //id выходной грани

	friend std::ostream& operator<< (std::ostream& out, const cell_local& point);
};

enum class_vtk_t
{
	e_class_grid_default = 0,
	e_class_grid_static_illum = 10,
};

#define OPEN_FSTREAM(file, namefile) \
file.open(namefile); \
if (!file.is_open()) RETURN_ERRS("Error : file %s is not open\n", namefile);

#ifdef CLASTER
#define Files_log F_LOG
#else
#define Files_log std::string(glb_files.base_adress + F_LOG).c_str()
#endif

#define WRITE_LOG_ERR(str){  \
std::ofstream ofile; \
ofile.open(Files_log, std::ios::app); \
ofile << str; \
ofile.close(); }


#define D_LD \
/*__pragma("omp critical")*/ \
{\
std::ofstream out(Files_log, std::ios::app); \
out << "DIE: " << __FILE__ << " " << __FUNCTION__ << ", " <<__LINE__<<"c.\n"; \
out.close(); \
/*MPI_END  MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COMM);*/ \
exit(1); \
}

#define DIE_IF(_cond) if(_cond) D_LD;

#ifndef CLASTER
#define CREATE_DIR(dir) \
struct stat st; \
if (stat(std::string(dir).c_str(), &st)) \
{ \
	if(mkdir(std::string(dir).c_str())) D_LD; \
}
#else
#define CREATE_DIR(dir){}
#endif

#ifdef USE_MPI
#define MPI_GET_INF(_number_of_nodes, _num_cur_id){ \
MPI_Comm_size(MPI_COMM_WORLD, &_number_of_nodes); \
MPI_Comm_rank(MPI_COMM_WORLD, &_num_cur_id);}
#else
#define MPI_GET_INF(_number_of_nodes, _num_cur_id){_number_of_nodes = 1; _num_cur_id =0;}
#endif

#ifdef WRITE_GLOBAL_LOG	

#ifdef WRITE_LOG_ON_SCREAN
#define WRITE_LOG(str){std::cout<<str;}
#else

#define WRITE_LOG(str) WRITE_LOG_ERR(str)

#ifdef WRITE_MPI_LOG
#define WRITE_LOG_MPI(str, id){  \
std::ofstream ofile; \
ofile.open(Files_log + std::to_string(id)+".txt", std::ios::app); \
ofile << str; \
ofile.close(); }
#else
#define WRITE_LOG_MPI(str, id) {CONVERT_TO_STRING(str);CONVERT_TO_STRING(id); }
#endif


#endif  //WRITE_LOG_ON_SCREAN
#else
#define WRITE_LOG(str) {CONVERT_TO_STRING(str);}
#define WRITE_LOG_MPI(str, id) {CONVERT_TO_STRING(str);CONVERT_TO_STRING(id); }
#endif //WRITE_GLOBAL_LOG

#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

#define check_bit(word, idx) ((word >> (idx)) & 0x1)  // проверка i-го бита
#define clear_bit(word, idx) (word & (~(1 << (idx)))) // выключение i-го бита
#define set_bit(word, idx) (word | (1 << (idx)))      // установка i-го бита

#define WRITE_FILE_VECTOR(name_file, data, value) \
{ \
FILE* f;\
int n = data.size(); \
f = fopen(name_file, "wb"); \
if(!f) RETURN_ERRS("file %s not open\n",name_file); \
fwrite(&n, sizeof(int), 1, f); \
for (auto& el : data) \
{	\
	fwrite(&el.value, sizeof(el.value), 1, f);	\
}	\
fclose(f);\
}

#define WRITE_FILE(name_file, data, n) \
{ \
FILE* f;\
f = fopen(name_file, "wb"); \
if(!f) RETURN_ERRS("file %s not open\n",name_file); \
fwrite(&n, sizeof(int), 1, f); \
fwrite(data, sizeof(data[0]), n, f);\
fclose(f);\
}

#define WRITE_FILE_PHYS(name_file, data, value, type, param) \
{ \
FILE* f;\
int n = data.size(); \
f = fopen(name_file, "wb"); \
if(!f) RETURN_ERRS("file %s not open\n",name_file); \
fwrite(&n, sizeof(int), 1, f); \
for (auto& el : data) \
{	\
type x = el.value*param; \
fwrite(&x, sizeof(x), 1, f);	\
}	\
fclose(f); \
}


#ifdef USE_MPI
#define MPI_START(argc, argv) MPI_Init(&argc, &argv);
#define MPI_END MPI_Finalize();
#define EXIT(a) {WRITE_LOG_ERR("Err calls from"<<__FILE__<<": " <<__FUNCTION__<<", "<<__LINE__<<" c.\n"); PRINT_POS /*MPI_END*/ exit(a); }
#else
#define MPI_START(argc, argv) {}
#define MPI_END {}
#define EXIT(a) {PRINT_POS exit(a);}
#endif //USE_MPI

#define SIGN(a) (a < 0.0 ? -1.0 : 1.0) 

#ifdef  USE_MPI
extern MPI_Datatype MPI_flux_t;
extern MPI_Datatype MPI_hllc_value_t;
extern MPI_Datatype MPI_flux_all_t;
extern MPI_Datatype MPI_flux_elem_t;
#endif

#define FOR_EACH_1(what,x)		what(x);
#define FOR_EACH_2(what,x, ...) what(x); FOR_EACH_1(what,__VA_ARGS__);
#define FOR_EACH_3(what,x, ...) what(x); FOR_EACH_2(what,__VA_ARGS__);
#define FOR_EACH_4(what,x, ...) what(x); FOR_EACH_3(what,__VA_ARGS__);
#define FOR_EACH_5(what,x, ...) what(x); FOR_EACH_4(what,__VA_ARGS__);
#define FOR_EACH_6(what,x, ...) what(x); FOR_EACH_4(what,__VA_ARGS__);
#define FOR_EACH_7(what,x, ...) what(x); FOR_EACH_4(what,__VA_ARGS__);
#define FOR_EACH_8(what,x, ...) what(x); FOR_EACH_4(what,__VA_ARGS__);
#define FOR_EACH_ID(_1,_2,_3,_4,_5,_6,_7,_8,NAME,...) NAME

//#define CREATE_TYPE(type, ...) type __VA_ARGS__
//
//#define CREATE_STRUCT_(name, ...) struct name{ FOR_EACH_ID(__VA_ARGS__,\
//FOR_EACH_8, FOR_EACH_7, FOR_EACH_6, FOR_EACH_5, FOR_EACH_4, FOR_EACH_3, FOR_EACH_2, FOR_EACH_1)\
//(CREATE_TYPE, __VA_ARGS__) }
//#define CREATE_STRUCT(name, type, ...) CREATE_STRUCT_(name, CREATE_TYPE(type, __VA_ARGS__))

#define CREATE_STRUCT(name, type, ...) struct name {type __VA_ARGS__ ;}

#define PRINT_STRUCT(st, type) \
type* str = (type*)&st; \
while (str < (type*)&st + sizeof(st) / sizeof(type)) \
{\
	std::cout << *str++ << '\n'; \
}

#define FILL_STRUCT(st, type, val) {\
type* str = (type*)&st; \
while (str < (type*)&st + sizeof(st) / sizeof(type)) \
{\
	*str++ = val; \
}}

#if 0 // под intel22 собираетс€
#define STRINGIZE(arg) #arg
#define CONCATENATE(arg1, arg2)  arg1##arg2

#define FOR_EACH_1(what, x, ...) what(x);
#define FOR_EACH_2(what, x, ...) what(x); FOR_EACH_1(what,  __VA_ARGS__);
#define FOR_EACH_3(what, x, ...) what(x); FOR_EACH_2(what, __VA_ARGS__);
#define FOR_EACH_4(what, x, ...) what(x); FOR_EACH_3(what,  __VA_ARGS__);
#define FOR_EACH_5(what, x, ...) what(x); FOR_EACH_4(what,  __VA_ARGS__);
#define FOR_EACH_6(what, x, ...) what(x); FOR_EACH_5(what,  __VA_ARGS__); 
#define FOR_EACH_7(what, x, ...) what(x); FOR_EACH_6(what,  __VA_ARGS__);
#define FOR_EACH_8(what, x, ...) what(x); FOR_EACH_7(what,  __VA_ARGS__);

#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, N, ...) N 
#define FOR_EACH_RSEQ_N() 8, 7, 6, 5, 4, 3, 2, 1, 0
#define FOR_EACH_NARG_(...) FOR_EACH_ARG_N(__VA_ARGS__) 
#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())


#define FOR_EACH_(N, what, x, ...) CONCATENATE(FOR_EACH_, N)(what, x, __VA_ARGS__)
#define FOR_EACH(what, x, ...) FOR_EACH_(FOR_EACH_NARG(x, __VA_ARGS__), what, x, __VA_ARGS__)

//#define CREATE_STRUCT(name, type, ...) struct name { FOR_EACH(type, __VA_ARGS__) }
#endif

#endif //GLOBAL_DEF
