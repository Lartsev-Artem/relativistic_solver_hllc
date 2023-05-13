#if !defined REBUILD_SOLVE_H && defined UTILS
#define REBUILD_SOLVE_H

#include "../prj_config.h"
#include "../global_headers.h"
#include "../global_def.h"

#if defined USE_VTK
int rebuild_solve(int argc, char* argv[], file_name name_file_settings);
int rebuild_solve(int argc, char* argv[]);
int recalc_grid_data(int argc, char* argv[]);

int rewrite_vtk_array(int argc, char* argv[], file_name name_file_settings);
int  RunMake1d(int argc, char* argv[]);
int BuildHLLC_1dTime(int argc, char* argv[]);

int GetPhysScale();
int MakeHllcInitFile(file_name base_adress);

int MakeHllcInitFromGrid(int argc, char* argv[]);

int CompareFiles(int argc, char* argv[]);
#endif

#endif //REBUILD_SOLVE_H