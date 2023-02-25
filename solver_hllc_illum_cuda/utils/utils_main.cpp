#include "../prj_config.h"
#ifdef UTILS
#include "rebuild_solve.h"
#include "rebuild_grid.h"
#include "grid_geometry/get_grid_data.h"
#include "../global_value.h"

int RunUtilsModule(int argc, char* argv[], const std::string& settings)
{		
	CompareFiles(argc, argv);

//	ReNumberingGrid(argc, argv);
	
//	ReBuildNetgenToMetis(argc, argv);

//	GetPhysScale();

//	MakeHllcInitFile(BASE_ADRESS);

#ifdef USE_VTK

//  MakeHllcInitFromGrid(argc, argv, settings);

//	rebuild_solve(argc, argv, settings);

//	BuildHLLC_1dTime(argc, argv);

//	RunMake1d(argc, argv);

//	rewrite_vtk_array(argc, argv, settings);

//	SetScalarDataVtkFromFile(argc, argv);

#if NUMBER_OF_MEASUREMENTS == 2
	// Trace2D(argc, argv);

	// GetAverageArea(argc, argv);

	// ReBuildNetgenToVTK2d(argc, argv);
#endif //2d

#endif //USE_VTK
	return 0;
}
#endif //UTILS
