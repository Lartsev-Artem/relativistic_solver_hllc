#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "../prj_config.h"
#ifdef MAKE

#include "../global_def.h"
#include "../global_headers.h"

//#define ReadNow  // при выполнение чиатать данные из файла по мере использования данных. (Для экономии ОЗУ)

extern Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра
extern int num_cur_direction; // номер текущего направления
extern Vector3 cur_direction;

// параметры диска и внутренней сферы:
const Type Rsphere = 0.001;
const Type R1disk = 0.001;
const Type R2disk = 0.09;

extern std::vector<Vector3> x_try_surface;
extern std::vector<int> id_try_surface;

extern int pos_x_try;
extern int posX;
extern int posX0;
extern int posOutC;
extern int posOut;
extern int posIn;
extern int posS;
extern int posRes;

extern ShortId n_out;

#endif //MAKE

#endif