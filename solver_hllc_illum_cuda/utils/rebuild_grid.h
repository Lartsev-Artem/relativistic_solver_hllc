#ifndef REBUILD_GRID_H
#define REBUILD_GRID_H
#include "../prj_config.h"
#ifdef UTILS

#include "../global_headers.h"

int ReBuildNetgenToVTK2d(int argc, char* argv[]);
int ReNumberingGrid(int argc, char* argv[]);
int ReBuildNetgenToMetis(int argc, char* argv[]);

#ifdef USE_VTK
int SetScalarDataVtkFromFile(int argc, char* argv[]);
#endif //USE_VTK

#endif // UTILS
#endif //REBUILD_GRID_H