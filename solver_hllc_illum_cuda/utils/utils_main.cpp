#ifdef UTILS
#include "../prj_config.h"
#include "rebuild_solve.h"
#include "rebuild_grid.h"
#include "grid_geometry/get_grid_data.h"
#include "../global_value.h"
#include "files_utils.h"

int RunUtilsModule(int argc, char* argv[], const std::string& settings)
{		
	//MakeSimpleTxtFromData(argc, argv);
	//WriteGeoFiles(argc, argv);

	//DeleteSolveFiles(argc, argv);
	//ReduceNameSolveFiles(argc, argv);
	//CopySolveFiles(argc, argv);

	//CompareFiles(argc, argv);

	//ReNumberingGrid(argc, argv);

	//GetReNumberGraphMetisToMPI(argc, argv);
	
	//RenumberNodesMetis(argc, argv);
	//GetMpiConfToHLLC(argc, argv);

	//GetReNumberGraphMetis(argc, argv);
	
//	ReBuildNetgenToMetis(argc, argv);

//	GetPhysScale();

//	MakeHllcInitFile(glb_files.base_adress);

#ifdef USE_VTK

//	recalc_grid_data(argc, argv);
//  MakeHllcInitFromGrid(argc, argv);

//	rebuild_solve(argc, argv, settings);
//	rebuild_solve(argc, argv);

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
