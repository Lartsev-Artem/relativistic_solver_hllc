#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "../prj_config.h"
#ifdef MAKE

#include "../global_def.h"
#include "../global_headers.h"
#include "../global_value.h"

// параметры диска и внутренней сферы:
const Type Rsphere = 0.001;
const Type R1disk = 0.001;
const Type R2disk = 0.09;

extern std::vector<Vector3> x_try_surface;
extern std::vector<int> id_try_surface;

#endif //MAKE

#endif