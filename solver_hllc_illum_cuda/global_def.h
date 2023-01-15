#ifndef GLOBAL_DEF
#define GLOBAL_DEF

#include "prj_config.h"
#include "global_headers.h"

#ifdef USE_MPI
#define MPI_START(argc, argv) MPI_Init(&argc, &argv);
#define MPI_END MPI_Finalize();
#define EXIT(a) { MPI_END exit(a); }
#else
#define MPI_START(argc, argv) {}
#define EXIT(a) exit(a);
#endif //USE_MPI

#define RETURN_ERR(str) { printf(str); return 1; }

#define base (NUMBER_OF_MEASUREMENTS + 1)

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;

typedef double Type;
typedef int IntId;
typedef uint8_t State;

typedef uint8_t ShortId;
typedef std::string Str_Type;

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
	Face& operator=(const Face& face) {
		A = face.A;
		B = face.B;
		C = face.C;
	}
};

struct FaceCell {
	int face_id;
	Face face;
	FaceCell(const int id = 0, const Face& face_init = Face()) {
		face_id = id;
		face = face_init;
	}
};

//сравнить с solve
struct cell {
	//int id;  - номер в массиве
	std::vector<Vector3> nodes_value;
	std::vector<int> neighbours_id_face;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3(-666, -666, -666));
		neighbours_id_face.resize(4, -1);
	}
};


#define OPEN_FSTREAM(file, namefile) \
file.open(namefile); \
if (!file.is_open()) { RETURN_ERR("Error : file %s is not open\n", namefile) }


#ifdef WRITE_GLOBAL_LOG	
#define WRITE_LOG(str){  \
ofstream ofile; \
ofile.open(BASE_ADRESS + "File_Logs.txt", std::ios::app); \
ofile << str; \
ofile.close(); }
#else
#define WRITE_LOG(ofile, str) {}
#endif

#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

#endif //GLOBAL_DEF
