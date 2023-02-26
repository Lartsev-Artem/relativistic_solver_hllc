#if !defined GET_GRID_DATA_H && defined UTILS
#define GET_GRID_DATA_H

#include "../../prj_config.h"
#include "../../global_headers.h"

#if defined USE_VTK
int Trace2D(int argc, char* argv[]);
int GetAverageArea(int argc, char* argv[]);
#endif

#endif //GET_GRID_DATA_H