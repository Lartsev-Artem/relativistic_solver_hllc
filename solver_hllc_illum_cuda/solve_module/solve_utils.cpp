#include "solve_utils.h"
#ifdef SOLVE

#if NUMBER_OF_MEASUREMENTS == 3 
void ReBuildNeighStruct(
	std::vector<int>& neighbours_id_faces, 
	std::vector<Normals>& normals, 
	std::vector<Type>& squares_faces,
	std::vector<Type>& volume,
	std::vector<Vector3>& centers,
	std::vector<face_t>& faces, std::vector<elem_t>& cells)
{	
	const int N = centers.size();
	cells.resize(N);	

	int cc = 0;
	for (int i = 0; i < N * base; i++)
	{
		int idx = neighbours_id_faces[i];

		if (idx != -10)
		{
			face_t f;
			f.geo.id_l = i / base; //ячейка

			f.geo.n = normals[i / base].n[i % base];
			f.geo.S = squares_faces[i];

			cells[i / base].geo.sign_n[i % base] = true;
			cells[i / base].geo.id_faces[i % base] = cc;

			neighbours_id_faces[i] = -10;
			if (idx >= 0)
			{
				f.geo.id_r = idx / base; //сосед

				cells[idx / base].geo.sign_n[idx % base] = false;
				cells[idx / base].geo.id_faces[idx % base] = cc;

				neighbours_id_faces[idx] = -10;
			}
			else
			{
				f.geo.id_r = idx; // код границы
			}

			faces.push_back(f); // как потом искать с ячейками?
			cc++;
		}
	}

	neighbours_id_faces.clear();
	normals.clear();
	squares_faces.clear();

	for (int i = 0; i < N; i++)
	{
		cells[i].geo.center = centers[i];
		cells[i].geo.V = volume[i];

	}
	
	volume.clear();
	centers.clear();

	return;
}
#endif //3d

#include "../file_module/reader_bin.h"
#include "../file_module/writer_bin.h"
int ReWriteGeoFiles(file_name name_file_geometry_faces, file_name name_file_geometry_cells)
{
	grid_t grid;
	WRITE_LOG("Error reading grid, try read parts geo\n");
	const std::string name_file_id_neighbors = BASE_ADRESS + "pairs.bin";
	const std::string name_file_normals = BASE_ADRESS + "normals.bin";
	const std::string name_file_centers = BASE_ADRESS + "centers.bin";
	const std::string name_file_squares = BASE_ADRESS + "squares.bin";
	const std::string name_file_volume = BASE_ADRESS + "volume.bin";

	std::vector<int> neighbours_id_faces;
	std::vector<Normals> normals;
	std::vector<Type> squares_faces;
	std::vector<Type> volume;
	std::vector<Vector3> centers;

	if (ReadSimpleFileBin(name_file_id_neighbors, neighbours_id_faces)) RETURN_ERR("Error reading file neighbours\n");
	if (ReadNormalFile(name_file_normals, normals)) RETURN_ERR("Error reading file normals\n");
	if (ReadSimpleFileBin(name_file_squares, squares_faces)) RETURN_ERR("Error reading file squares_faces\n");
	if (ReadSimpleFileBin(name_file_volume, volume)) RETURN_ERR("Error reading file volume\n");
	if (ReadSimpleFileBin(name_file_centers, centers)) RETURN_ERR("Error reading file centers\n");

	ReBuildNeighStruct(neighbours_id_faces, normals, squares_faces, volume, centers, grid.faces, grid.cells);
	grid.size = grid.cells.size();

	return WriteGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid);
}


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
	const std::string name_file_id_neighbors = BASE_ADRESS + "pairs.bin";
	const std::string name_file_normals = BASE_ADRESS + "normals.bin";
	const std::string name_file_centers = BASE_ADRESS + "centers.bin";
	const std::string name_file_squares = BASE_ADRESS + "squares.bin";
	const std::string name_file_volume = BASE_ADRESS + "volume.bin";

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
int TestDivStream(file_name BASE_ADRESS)
{
	const std::string name_file_geometry_faces = BASE_ADRESS + "geo_faces.bin";
	const std::string name_file_geometry_cells = BASE_ADRESS + "geo_cells.bin";
	const std::string name_file_centers_faces = BASE_ADRESS + "center_face.bin";

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
	WriteFileSolution(BASE_ADRESS+"Solve\\Solve0", std::vector<Type>(), grid.cells); //печать начальной сетки

	return 0;
}
#endif //RUN_TEST
#endif //SOLVE

#ifdef USE_MPI
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

#endif //USE_MP