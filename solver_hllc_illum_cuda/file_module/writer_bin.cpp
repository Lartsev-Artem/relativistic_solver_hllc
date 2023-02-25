#include "../prj_config.h"

#include "../global_headers.h"
#include "../global_def.h"

#include "reader_vtk.h"
#include "writer_bin.h"

#if 1 //def SOLVE
#include "../solve_module/solve_config.h"
#include "../solve_module/solve_global_struct.h"
#include "../solve_module/illum/illum_utils.h"

size_t WriteFileSolution(const std::string& main_dir, const grid_t& grid) 
{

#ifdef ILLUM 
	std::vector<Type> illum;
	GetDirectionIllumFromFace(grid.size, 0, grid.Illum, illum);
	WriteSimpleFileBin(main_dir + "Illum.bin", illum);

#if !defined USE_CUDA
	WRITE_FILE_VECTOR((main_dir + "energy.bin").c_str(), grid.cells, illum_val.energy);

	WRITE_FILE_VECTOR((main_dir + "stream.bin").c_str(), grid.cells, illum_val.stream);

	WRITE_FILE_VECTOR((main_dir + "impuls.bin").c_str(), grid.cells, illum_val.impuls);

	WRITE_FILE_VECTOR((main_dir + "divstream.bin").c_str(), grid.cells, illum_val.div_stream);

	WRITE_FILE_VECTOR((main_dir + "divimpuls.bin").c_str(), grid.cells, illum_val.div_impuls);
#else 
	WRITE_FILE((main_dir + "energy.bin").c_str(), grid.energy, grid.size);

	WRITE_FILE((main_dir + "stream.bin").c_str(), grid.stream, grid.size);	

	WRITE_FILE((main_dir + "impuls.bin").c_str(), grid.impuls, grid.size);

	WRITE_FILE((main_dir + "divstream.bin").c_str(), grid.divstream, grid.size);

	WRITE_FILE((main_dir + "divimpuls.bin").c_str(), grid.divimpuls, grid.size);
#endif
#endif

#if defined HLLC || defined RHLLC

#ifdef Cone	
	WRITE_FILE_PHYS((main_dir + "density.bin").c_str(), grid.cells, phys_val.d, Type, DENSITY);

	WRITE_FILE_PHYS((main_dir + "pressure.bin").c_str(), grid.cells, phys_val.p, Type, PRESSURE);

	WRITE_FILE_VECTOR((main_dir + "velocity.bin").c_str(), grid.cells, phys_val.v);

	//WRITE_FILE_PHYS((main_dir + "velocity.bin").c_str(), cells, phys_val.v, Vector3, VELOCITY);
#else
	WRITE_FILE((main_dir + "density.bin").c_str(), cells, phys_val.d);

	WRITE_FILE((main_dir + "pressure.bin").c_str(), cells, phys_val.p);

	WRITE_FILE((main_dir + "velocity.bin").c_str(), cells, phys_val.v);
#endif
#endif

	return 0;
}


int WriteGeometryGrid(const std::string& file_cells, const std::string& file_faces, grid_t& grid)
{		
	WRITE_FILE_VECTOR(file_faces.c_str(), grid.faces, geo);
	WRITE_FILE_VECTOR(file_cells.c_str(), grid.cells, geo);
	return 0;
}

#endif //SOLVE

