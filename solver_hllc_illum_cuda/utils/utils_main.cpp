#include "../prj_config.h"
#ifdef UTILS
#include "rebuild_solve.h"
#include "rebuild_grid.h"
#include "grid_geometry/get_grid_data.h"
#include "../global_value.h"

int RunUtilsModule(int argc, char* argv[], const std::string& settings)
{			
//	ReNumberingGrid(argc, argv);
	
//	ReBuildNetgenToMetis(argc, argv);

//	GetPhysScale();

//	MakeHllcInitFile(BASE_ADRESS);

#ifdef USE_VTK

	//delete files
	{
		int a = 0;
		int b = 0;
		int step = 4;
		std::string file_solve;
		std::string foo; int foo2;
		ReadStartSettings(settings, foo2, foo, foo, foo, BASE_ADRESS, file_solve, foo2);

		for (int i = a; i < b; i+=step)
		{
			//2,6,10
			std::remove((file_solve + std::to_string(i) + ".vtk").c_str());
		}
		return 0;
	}

 MakeHllcInitFromGrid(argc, argv, settings);

	//rebuild_solve(argc, argv, settings);

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
