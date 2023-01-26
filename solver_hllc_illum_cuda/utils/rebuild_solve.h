#ifndef REBUILD_SOLVE_H
#define REBUILD_SOLVE_H

#include "../prj_config.h"
#include "../global_headers.h"
#include "../global_def.h"

#if defined USE_VTK && defined UTILS
int rebuild_solve(int argc, char* argv[], file_name name_file_settings);
int rewrite_vtk_array(int argc, char* argv[], file_name name_file_settings);
int  RunMake1d(int argc, char* argv[]);
int BuildHLLC_1dTime(int argc, char* argv[]);

int GetPhysScale();
#endif

#endif //REBUILD_SOLVE_H