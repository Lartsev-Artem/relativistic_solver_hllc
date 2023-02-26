#if !defined BUILD_GRAPH_STRUCT && defined BUILD
#define BUILD_GRAPH_STRUCT
#include "../global_def.h"

extern std::vector<int> id_try_surface;		 // id граней, определяющих внутренюю границу
extern std::vector<double> dist_try_surface; // расстояния между точками (через полость внутри) 
extern std::vector<Vector3> x_try_surface;   // x точка выхода

extern uint64_t id_try_size;
extern uint64_t dist_try_size;
extern uint64_t x_try_size;

struct TrySolve {
	int id_1;
	int id_2;
	int id_3;

	double s_1;
	double s_2;
	double s_3;

	Vector3 x1;
	Vector3 x2;
	Vector3 x3;
};
extern TrySolve buf_try;

#endif 
