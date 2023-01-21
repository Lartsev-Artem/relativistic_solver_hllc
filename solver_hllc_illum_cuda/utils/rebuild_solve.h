#ifndef REBUILD_SOLVE_H
#define REBUILD_SOLVE_H

#include "../prj_config.h"
#include "../global_headers.h"

#if defined USE_VTK && defined UTILS
int rebuild_solve(file_name name_file_settings);
int rewrite_vtk_array(file_name name_file_settings);
int  RunMake1d(int argc, char* argv[]);
int BuildHLLC_1dTime(int argc, char* argv[]);
#endif

#endif //REBUILD_SOLVE_H