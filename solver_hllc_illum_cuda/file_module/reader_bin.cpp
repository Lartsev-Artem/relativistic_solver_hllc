#include "../prj_config.h"

#include "reader_bin.h"
//#if defined BUILD

#include "../global_headers.h"
#include "../global_def.h"

#ifdef USE_VTK
#include "reader_vtk.h"

#endif // USE_VTK

#ifdef  USE_MPI
#include "../solve_module/solve_utils.h"
#endif


int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals) {

	FILE* f;
	f = fopen(name_file_normals.c_str(), "rb");

	int n;
	fread(&n, sizeof(int), 1, f);
	normals.resize(n);

	Normals norm(base);
	for (size_t i = 0; i < n; i++)
	{
		for (int j = 0; j < base; j++)
			fread(&norm.n[j], sizeof(Vector3), 1, f);

		normals[i] = norm;
	}
	fclose(f);

	return 0;
}




int ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	std::vector<Vector3>& velocity, std::vector<Type>& pressure, const bool is_print/*=false*/) {

	switch (class_file_vtk)
	{
	case 1: //test grid
		ReadSimpleFileBin(main_dir + "alpha.bin", density);
		ReadSimpleFileBin(main_dir + "alpha.bin", absorp_coef);
		ReadSimpleFileBin(main_dir + "Q.bin", rad_en_loose_rate);
		break;
	case 2: //main grid to illum
		ReadSimpleFileBin(main_dir + "density.bin", density);
		ReadSimpleFileBin(main_dir + "AbsorpCoef.bin", absorp_coef);
		ReadSimpleFileBin(main_dir + "radEnLooseRate.bin", rad_en_loose_rate);
		break;

	case 3: //full_test_grid
		ReadSimpleFileBin(main_dir + "density.bin", density);
		ReadSimpleFileBin(main_dir + "alpha.bin", absorp_coef);
		ReadSimpleFileBin(main_dir + "Q.bin", rad_en_loose_rate);
		ReadSimpleFileBin(main_dir + "velocity.bin", velocity);
		ReadSimpleFileBin(main_dir + "pressure.bin", pressure);
		break;
	default:
		printf("Grid without data\n");
		return 0;
	}

	if (is_print)
	{
		std::cout << "density_Size: " << density.size() << '\n';
		std::cout << "absorp_coef_Size: " << absorp_coef.size() << '\n';
		std::cout << "Q_Size: " << rad_en_loose_rate.size() << '\n';
	}

	const int a = density.size();
	const int b = rad_en_loose_rate.size();
	const int c = absorp_coef.size();
	if (a == b && a == c && b == c)
	{
		return 0;
	}

	RETURN_ERR("Error size data\n");
	return 1;
}


#ifdef SOLVE
#ifdef ILLUM
int ReadDataArray(const solve_mode_t& mode, file_name main_dir, grid_t& grid)
{
	switch (mode.class_vtk)
	{
	case 1: //test grid
		READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells, phys_val.d);
		READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells, illum_val.absorp_coef);
		READ_FILE((main_dir + "Q.bin").c_str(), grid.cells, illum_val.rad_en_loose_rate);
		break;

	case 2: //main grid to illum
		READ_FILE((main_dir + "density.bin").c_str(), grid.cells, phys_val.d);
		READ_FILE((main_dir + "AbsorpCoef.bin").c_str(), grid.cells, illum_val.absorp_coef);
		READ_FILE((main_dir + "radEnLooseRate.bin").c_str(), grid.cells, illum_val.rad_en_loose_rate);
		break;

	case 3: //full_test_grid
		READ_FILE((main_dir + "density.bin").c_str(), grid.cells, phys_val.d);
		READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells, illum_val.absorp_coef);
		READ_FILE((main_dir + "Q.bin").c_str(), grid.cells, illum_val.rad_en_loose_rate);
		READ_FILE((main_dir + "velocity.bin").c_str(), grid.cells, phys_val.v);
		READ_FILE((main_dir + "pressure.bin").c_str(), grid.cells, phys_val.p);
		break;
	default:
		printf("Grid without data\n");
		return 0;
	}

	return 0;
}
#endif //ILLUM
int ReadGeometryGrid(const std::string& file_cells, const std::string& file_faces, grid_t& grid)
{
	READ_FILE(file_faces.c_str(), grid.faces, geo);
	READ_FILE(file_cells.c_str(), grid.cells, geo);	
	grid.size = grid.cells.size();
	return 0;
}

int ReadIllumGeometry(const int count_dir,
	file_name file_x, file_name file_state, file_name file_x0, file_name file_graph, file_name file_res,
	std::vector<BasePointTetra>& vec_x,
	std::vector <std::vector<int>>& face_states,
	std::vector <std::vector<cell_local>>& vec_x0,
	std::vector < std::vector<int>>& sorted_id_cell,
	std::vector<Type>& vec_res_bound)
{
	if (ReadSimpleFileBin(file_x, vec_x)) return 1;
	if (ReadSimpleFileBin(file_res, vec_res_bound)) return 1;

	std::vector<int> disp;
	std::vector<int> send;

#ifdef  USE_MPI
	int np, myid;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	GetDisp(np, count_dir, disp);
	GetSend(np, count_dir, send);		
#else
	disp.resize(1, 0);
	send.resize(1, count_dir);
	int myid = 0;
#endif

	vec_x0.resize(send[myid]);
	sorted_id_cell.resize(send[myid]);
	face_states.resize(send[myid]);

	for (int i = 0; i < send[myid]; i++)	
	{
		if (ReadSimpleFileBin(file_x0 + std::to_string(disp[myid]+i) + ".bin", vec_x0[i])) return 1;
		if (ReadSimpleFileBin(file_graph + std::to_string(disp[myid]+i) + ".bin", sorted_id_cell[i])) return 1;
		if (ReadSimpleFileBin(file_state + std::to_string(disp[myid]+i) + ".bin", face_states[i])) return 1;
	}

	return 0;
}

int ReadValueGrid(const std::string& main_dir, grid_t& grid)
{	
	READ_FILE((main_dir + "phys_val.bin").c_str(), grid.cells, phys_val);

#ifdef ILLUM

#endif

	return 0;
}

#if defined HLLC || defined RHLLC
int ReadHllcInit(file_name file_init_value, std::vector<elem_t>& cells)
{
	READ_FILE((file_init_value).c_str(), cells, phys_val);	
	return 0;
}
#endif

#endif //SOLVE
//#endif //BUILD