#ifdef SOLVE
#include "solve_utils.h"
#include "../global_def.h"


#include "hllc/hllc_utils.h"
#include "rhllc/rhllc_utils.h"

int HLLC_INIT(file_name file_settings_hllc, hllc_value_t& hllc_set,
	file_name file_init_value, std::vector<elem_t>& cells)
{
#if !defined RHLLC_MPI
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 3 
	return InitHLLC(file_settings_hllc, hllc_set, file_init_value, cells);  //Начальные данные для HLLC
#elif defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 
	return InitRHLLC(file_settings_hllc, hllc_set, file_init_value, cells);
#endif

#else //!USE_MPI
#if defined RHLLC&& NUMBER_OF_MEASUREMENTS == 3
		
	if (InitRHLLC(file_settings_hllc, hllc_set, file_init_value, cells))
	{
		D_LD;
	}

	InitMPI_RHllc(cells);

	return 0;
#endif
	return 0;
#endif//USE_MPI
	return 0;
}

int HLLC_STEP(const Type tau, grid_t& grid)
{
#if !defined RHLLC_MPI
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 3 
	return HLLC_3d(tau, grid);
#elif defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 
	return RHLLC_3d(tau, grid);
#endif

#else //!USE_MPI
	return 0;
#endif//USE_MPI
	return 0;
}

int GetTimeStep(hllc_value_t& hllc_cfg, const grid_t& grid)
{
#if !defined RHLLC_MPI
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 3 
	return HllcGetTimeStep(hllc_cfg, grid.cells);
#elif defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 
	return RHllcGetTimeStep(hllc_cfg, grid.cells);
#endif	

#else //!USE_MPI
	return 0;
#endif//USE_MPI
	return 0;
}


int StartLowDimensionTask(file_name main_dir)
{
	const std::string name_file_id_neighbors = glb_files.base_adress + "pairs.bin";
	const std::string name_file_normals = glb_files.base_adress + "normals.bin";
	const std::string name_file_centers = glb_files.base_adress + "centers.bin";
	const std::string name_file_squares = glb_files.base_adress + "squares.bin";
	const std::string name_file_volume = glb_files.base_adress + "volume.bin";

#if NUMBER_OF_MEASUREMENTS == 3
	return 0;
#endif

#if defined HLLC_1D
	HLLC_1d(main_dir);
	printf("HLLC 1d end\n");
	return 0;

#elif defined RHLLC_1D
	RHLLC_1d(main_dir);
	printf("RHLLC 1d end\n");
	return 0;
#endif

#if defined HLLC_2D

	HLLC2d(main_dir, centers, neighbours_id_faces, normals, squares_cell, volume);
	printf("End 2d hllc\n");
	return 0;

#elif defined RHLLC_2D
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	_clock = -omp_get_wtime();

	MPI_RHLLC(main_dir, centers, neighbours_id_faces, normals, squares_cell, volume);

	_clock += omp_get_wtime();
	printf("Full mpi time= %lf\n", _clock);
	printf("End 2d hllc\n");

	return 0;
#else
	_clock = -omp_get_wtime();

	RHLLC2d(main_dir, centers, neighbours_id_faces, normals, squares_cell, volume);

	_clock += omp_get_wtime();
	printf("Full omp time= %lf\n", _clock);
	printf("End 2d hllc\n");
	return 0;
#endif
#endif
}

#ifdef RUN_TEST
int TestDivStream(file_name glb_files.base_adress)
{
	const std::string name_file_geometry_faces = glb_files.base_adress + "geo_faces.bin";
	const std::string name_file_geometry_cells = glb_files.base_adress + "geo_cells.bin";
	const std::string name_file_centers_faces = glb_files.base_adress + "center_face.bin";

	grid_t grid;
	if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid))
	{
		WRITE_LOG("Error reading grid, try read parts geo\n");

		ReWriteGeoFiles(name_file_geometry_faces, name_file_geometry_cells);

		if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid)) RETURN_ERR("Error reading grid\n");
	}
	WRITE_LOG("Reading geometry grid\n");


	std::vector<Vector3> center_cells;
	ReadSimpleFileBin(name_file_centers_faces, center_cells);
	TestDivStream(center_cells, grid);
	WriteFileSolution(glb_files.base_adress+"Solve\\Solve0", std::vector<Type>(), grid.cells); //печать начальной сетки

	return 0;
}
#endif //RUN_TEST
#endif //SOLVE

#ifdef USE_MPI
// аналог классического getSend, но с разгрузкой управляющего узла на coef
void GetDispSend(const int np, const int n, const int coef, std::vector<int>& send_count,  std::vector<int>& disp)
{
	int first = std::max(0, n / np - coef);

	const int N = n - first;
	const int NP = np - 1;

	send_count.resize(np, N / NP);
	send_count[0] = first;

	if (N % NP)  // если число процессов не кратно размерности задачи 
	{
		for (int i = 0; i < N % NP; i++) // первые процессы берут на единицу больше тел
			++send_count[i + 1];
	}

	disp.resize(np, 0);
	for (int i = 0; i < np - 1; i++)
	{
		disp[i + 1] = disp[i] + send_count[i];
	}
	return;
}

void GetSend(const int np, const int n, std::vector<int>& send_count)
{
	// вычисление  числа тел на узел
	send_count.resize(np, n / np);

	if (n % np)  // если число процессов не кратно размерности задачи 
	{
		for (int i = 0; i < n % np; i++) // первые процессы берут на единицу больше тел
			++send_count[i];
	}
}

void GetDisp(const int np, const int n, std::vector<int>& disp)
{
	// вычисление сдвигов  на узел
	disp.resize(np, 0);
	int a = n % np;
	for (int i = 1; i < np; i++)
		disp[i] = i * (n / np) + a;

	if (n % np)  // если число процессов не кратно размерности задачи 
	{
		for (int i = 1; i < a; i++) // смещения для процессов за ними увеличивается начиная со второго
			disp[a - i] -= i;
	}
}


void Init_MPI()
{
	// структура потоков
	{
		int len[3 + 1] = { 1,3,1,  1 };
		MPI_Aint pos[4] = { offsetof(flux_t, d), offsetof(flux_t, v),offsetof(flux_t, p) ,sizeof(flux_t) };
		MPI_Datatype typ[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(4, len, pos, typ, &MPI_flux_t);
		MPI_Type_commit(&MPI_flux_t);
	}

	// перессылка потоков из ячейки
	{
		int len[2 + 1] = { 1,1,1 };
		MPI_Aint pos[3] = { offsetof(elem_t, phys_val), offsetof(elem_t, conv_val) ,sizeof(elem_t) };
		MPI_Datatype typ[3] = { MPI_flux_t,MPI_flux_t , MPI_UB };
		MPI_Type_struct(3, len, pos, typ, &MPI_flux_elem_t);
		MPI_Type_commit(&MPI_flux_elem_t);
	}

	// физические и консервативные потоки
	{
		int len[2 + 1] = { 1,1,1 };
		MPI_Aint pos[3] = { offsetof(flux_all_t, phys_val), offsetof(flux_all_t, conv_val) ,sizeof(flux_all_t) };
		MPI_Datatype typ[3] = { MPI_flux_t,MPI_flux_t , MPI_UB };
		MPI_Type_struct(3, len, pos, typ, &MPI_flux_all_t);
		MPI_Type_commit(&MPI_flux_all_t);
	}

	//настройки динамического расчета
	{
		int len[5 + 1] = { 1,1,1,1,1,  1 };
		MPI_Aint pos[6] = { offsetof(hllc_value_t,T),offsetof(hllc_value_t,CFL),offsetof(hllc_value_t,h)
			,offsetof(hllc_value_t,print_timer) ,offsetof(hllc_value_t,tau) ,sizeof(hllc_value_t) };
		MPI_Datatype typ[6] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(6, len, pos, typ, &MPI_hllc_value_t);
		MPI_Type_commit(&MPI_hllc_value_t);
	}

	int np, myid;
	MPI_GET_INF(np, myid);
	ReadMpiConf(glb_files.base_adress + F_MPI_CONF, myid, np);

	return;
}

#endif //USE_MP