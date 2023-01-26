#include "solve_global_struct.h"
#include "../global_def.h"

#ifdef SOLVE
#include "../global_headers.h"
#include "../global_value.h"

#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_bin.h"

#include "../utils/grid_geometry/geometry_solve.h"

#include "solve_utils.h"

#include"illum/illum_utils.h"
#include "hllc/hllc_utils.h"  //сделать общий для rhllc и hllc

solve_mode_t solve_mode;

/*
	ReadData (vtk or not)

	Rebuild/ rewrite section

	BuildHLLC struct

	Solve HLLC

	Solve illum -> divW, U, divT

	Solve HLLC + Illum

	next time step

*/

int RunSolveModule(const std::string& name_file_settings)
{

#ifdef  USE_MPI
	int np, myid;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#else
	int myid = 0;
#endif //  USE_MPI

	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string adress_graph_file;
	std::string adress_solve;
	std::string name_file_value_init = BASE_ADRESS + "hllc_init_value.bin";

	if (ReadStartSettings(name_file_settings, solve_mode.class_vtk, name_file_vtk, name_file_sphere_direction, adress_graph_file, BASE_ADRESS, adress_solve, solve_mode.max_number_of_iter,
		name_file_value_init))
	{
		RETURN_ERR("Error reading solve settings\n");
	}
	WRITE_LOG("start solve module\n");	

#if defined HLLC || defined RHLLC
	StartLowDimensionTask(BASE_ADRESS);
#endif

#ifdef RUN_TEST
	TestDivStream(BASE_ADRESS);
	return 0;
#endif // RUN_TEST
	
	const std::string name_file_hllc_set = BASE_ADRESS + "hllc_settings.txt";

	//--------------------------Файлы расчётных данных(отдельные файлы по направлениям)-----------------------------------//
	const std::string name_file_state_face = BASE_ADRESS + "state_face";
	const std::string name_file_x0_loc = BASE_ADRESS + "LocX0";
	const std::string name_file_x = BASE_ADRESS + "X.bin";

	//--------------------------Файлы трассировки сквозь внутреннюю границу----------------------------------------//
	const std::string name_file_dist_try = BASE_ADRESS + "dist_defining_faces";
	const std::string name_file_id_try = BASE_ADRESS + "id_defining_faces";
	const std::string name_file_res = BASE_ADRESS + "ResBound";
	const std::string name_file_neib = BASE_ADRESS + "pairs.bin";
	
#if 0
	//--------------------------Файлы сдвигов по направлениям в расчётных файлах-----------------------------------//
	const std::string name_file_shift_out = BASE_ADRESS + "ShiftOut";
	const std::string name_file_shift_res = BASE_ADRESS + "ShiftRes";
	const std::string name_file_shift_x0 = BASE_ADRESS + "ShiftX0";
	const std::string name_file_shift_try = BASE_ADRESS + "ShiftTry";
#endif

	//--------------------------Файлы геометрии--------------------------------------------------------------------//
	const std::string name_file_geometry_faces = BASE_ADRESS + "geo_faces.bin";
	const std::string name_file_geometry_cells = BASE_ADRESS + "geo_cells.bin";

	grid_t grid;
	if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid))		
	{
		WRITE_LOG("Error reading grid, try read parts geo\n");
		
		ReWriteGeoFiles(name_file_geometry_faces, name_file_geometry_cells);

		if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid)) RETURN_ERR("Error reading grid\n");
	}
	WRITE_LOG("Reading geometry grid\n");

	double _clock;	

	//------------------------------ Illume section-----------------------------
	std::vector<Type> Illum; // отдельно от сетки из-за расчета на cuda
	std::vector<Type> int_scattering; // отдельно от сетки из-за расчета на cuda

#ifdef ILLUM
	grid_directions_t grid_direction;
	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, grid_direction)) RETURN_ERR("Error reading the sphere direction\n");

	WRITE_LOG("reading sphere_direction\n");

	std::vector<BasePointTetra> vec_x;
	std::vector<Type> vec_res_bound;
	std::vector<int> pairs;

	const int cont_dir = grid_direction.size;
	std::vector < std::vector<int>> face_states(cont_dir); //битовое поле: 0=> выходящая грань,  1=> входящая   
	std::vector < std::vector<cell_local>> vec_x0(cont_dir);
	std::vector < std::vector<int>> sorted_id_cell(cont_dir); 	// Упорядоченные индексы ячеек по данному направлению	

	_clock = -omp_get_wtime();
	if (ReadIllumGeometry(cont_dir, name_file_x, name_file_state_face, name_file_x0_loc,
		adress_graph_file, name_file_res,
		vec_x, face_states, vec_x0, sorted_id_cell, vec_res_bound))
		RETURN_ERR("Error fast reading the data grid vtk\n");

	_clock += omp_get_wtime();

	WRITE_LOG("Reading time of the Illum Geometry file : " << _clock << "\n");

	Illum.resize(base * grid.size * grid_direction.size, 0);
	int_scattering.resize(grid.size * grid_direction.size, 0);
	
#ifdef USE_CUDA
	if (solve_mode.use_cuda)
	{
		if (CheckDevice()) 
		{			
			CUDA_ERR("Cuda wasn't found\n");
		}
		else
		{
			if (InitDevice(grid_direction.size, grid.size, solve_mode.cuda_mod)) CUDA_ERR("Cuda Error init\n");
			if (HostToDevice(grid_direction, Illum, 1))  CUDA_ERR("Error cuda hostTodevice\n");
		}
	}
#endif

	ReadSimpleFileBin(name_file_neib, pairs);
	for (auto& el : grid.cells)
	{
		el.illum_val.illum.resize(grid_direction.size * base);
	}
	InitIllum(BASE_ADRESS, grid);
#endif //ILLUM
	//------------------------------------------------------------------------------------------------------------

	//---------------------------------------Solve section--------------------------------------------------------------
			
	int res_count = 0; // счётчик решений	
	Type full_time = -omp_get_wtime();

#if defined HLLC || defined RHLLC
	hllc_value_t hllc_cfg;
	Type t = 0.0;
	Type cur_timer = 0;
	
	if (HLLC_INIT(name_file_hllc_set, hllc_cfg, name_file_value_init, grid.cells)) RETURN_ERR("Bad init hllc\n");  //Начальные данные для HLLC
	
	WriteFileSolution(adress_solve + std::to_string(res_count++), Illum, grid.cells); //печать начальной сетки

	WRITE_LOG("Start main task\n");

	while (t < hllc_cfg.T)
	{	
		Type time_step = -omp_get_wtime();
		HLLC_STEP(hllc_cfg.tau, grid);

#ifdef ILLUM

		CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
		
		CalculateIllumParam(grid_direction, grid);

		SolveIllumAndHLLC(hllc_cfg.tau, grid.cells);
		
		//energy.swap(prev_energy);
		//stream.swap(prev_stream);
#endif
		t += hllc_cfg.tau;
		cur_timer += hllc_cfg.tau;
		
		GetTimeStep(hllc_cfg, grid); //формирование шага по времени

		if (cur_timer >= hllc_cfg.print_timer)
		{
			WriteFileSolution(adress_solve + std::to_string(res_count++), Illum, grid.cells);						

			//WRITE_LOG("\nt= " << t << "; tau= " << hllc_cfg.tau << "; step= " << res_count << '\n');
			printf("\n t= %f,  tau= %lf,  res_step= %d\n", t, hllc_cfg.tau, res_count);
			cur_timer = 0;
		}

		time_step += omp_get_wtime();
		WRITE_LOG("\nt= " << t << "; tau= " << hllc_cfg.tau << "; step= " << res_count << "time= " << time_step << " c\n");
	}// while(t < T)

#elif defined ILLUM
	CalculateIllum(grid_direction, face_states,pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
	CalculateIllumParam(grid_direction, grid);
#else
	WRITE_LOG("No solve. Bad config\n");
	RETURN_ERR("No solve. Bad config\n");
#endif

	full_time += omp_get_wtime();
	printf("Time while(): %f\n", full_time);

#ifdef USE_CUDA
	if (solve_mode.use_cuda)
	{
		ClearDevice(solve_mode.cuda_mod);
	}
#endif

	WriteFileSolution(adress_solve + std::to_string(res_count), Illum, grid.cells);
	
	WRITE_LOG("End solve module\n");

	return EXIT_SUCCESS;
}

#endif //SOLVE