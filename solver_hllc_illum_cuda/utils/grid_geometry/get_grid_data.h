#ifndef GET_GRID_DATA_H
#define GET_GRID_DATA_H

#include "../../prj_config.h"
#include "../../global_headers.h"

#if defined USE_VTK && defined UTILS
int Trace2D(int argc, char* argv[]);
int GetAverageArea(int argc, char* argv[]);
#endif

#endif //GET_GRID_DATA_H