#if !defined REBUILD_GRID_H && defined UTILS
#define REBUILD_GRID_H

#include "../prj_config.h"
#include "../global_headers.h"

int ReBuildNetgenToVTK2d(int argc, char* argv[]);
int ReNumberingGrid(int argc, char* argv[]);
int ReBuildNetgenToMetis(int argc, char* argv[]);
void GenerateMeshSizeFileNetgen();

#ifdef USE_VTK
int SetScalarDataVtkFromFile(int argc, char* argv[]);
#endif //USE_VTK

#endif //REBUILD_GRID_H