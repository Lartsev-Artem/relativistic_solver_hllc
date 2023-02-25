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

#include "../cuda/cuda_solve.h"

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
	const std::string name_file_state_face = BASE_ADRESS + "illum_geo/state_face";
	const std::string name_file_x0_loc = BASE_ADRESS + "illum_geo/LocX0";
	const std::string name_file_x = BASE_ADRESS + "illum_geo/X.bin";

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
	if (myid == 0)
	{
		if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid))
		{
			WRITE_LOG("Error reading grid, try read parts geo\n");

			ReWriteGeoFiles(name_file_geometry_faces, name_file_geometry_cells);

			if (ReadGeometryGrid(name_file_geometry_cells, name_file_geometry_faces, grid)) RETURN_ERR("Error reading grid\n");
		}
		WRITE_LOG("Reading geometry grid\n");
	}

#ifdef USE_MPI
	MPI_Bcast(&grid.size, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif // USE_MPI
	
	double _clock;

	//------------------------------ Illume section-----------------------------
	//std::vector<Type> Illum; // отдельно от сетки из-за расчета на cuda
	//std::vector<Type> int_scattering; // отдельно от сетки из-за расчета на cuda

#ifdef ILLUM
	grid_directions_t grid_direction;
	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, grid_direction)) RETURN_ERR("Error reading the sphere direction\n");

	WRITE_LOG("reading sphere_direction\n");

	std::vector<BasePointTetra> vec_x;
	std::vector<Type> vec_res_bound;
	std::vector<int> pairs;

	const int count_dir = grid_direction.size;
	std::vector < std::vector<int>> face_states; //битовое поле: 0=> выходящая грань,  1=> входящая   
	std::vector < std::vector<cell_local>> vec_x0;
	std::vector < std::vector<int>> sorted_id_cell; 	// Упорядоченные индексы ячеек по данному направлению	

	_clock = -omp_get_wtime();
	if (ReadIllumGeometry(count_dir, name_file_x, name_file_state_face, name_file_x0_loc,
		adress_graph_file, name_file_res,
		vec_x, face_states, vec_x0, sorted_id_cell, vec_res_bound))
		RETURN_ERR("Error fast reading the data grid vtk\n");

	_clock += omp_get_wtime();

	WRITE_LOG("Reading time of the Illum Geometry file : " << _clock << "\n");

	ReadSimpleFileBin(name_file_neib, pairs);

	if (myid == 0)
	{
		/*grid.Illum;

		Illum.resize(base * grid.size * grid_direction.size, 0);
		int_scattering.resize(grid.size * grid_direction.size, 0);*/

#ifdef USE_CUDA //пока только на основном узле		
		SetDevice(0);
		InitDevice(grid_direction, grid);
#else
		grid.InitMemory(grid_direction.size);				
#endif			
		InitIllum(BASE_ADRESS, grid);
	} //myid==0

#ifdef USE_MPI
	InitSendDispIllumArray(myid, np, grid_direction.size, grid.size);
#endif
#endif //ILLUM
	//------------------------------------------------------------------------------------------------------------

	//---------------------------------------Solve section--------------------------------------------------------------
	int res_count = 0; // счётчик решений	
	Type full_time = -omp_get_wtime();

#if defined HLLC || defined RHLLC
	hllc_value_t hllc_cfg;
	Type t = 0.0;
	Type cur_timer = 0;

	if (myid == 0)
	{
		if (HLLC_INIT(name_file_hllc_set, hllc_cfg, name_file_value_init, grid.cells)) RETURN_ERR("Bad init hllc\n");  //Начальные данные для HLLC

		WriteFileSolution(adress_solve + std::to_string(res_count++), grid); //печать начальной сетки

		WRITE_LOG("Start main task\n");
	}

#ifdef USE_MPI
	{		
		int len[5 + 1] = { 1,1,1,1,1,  1 };
		MPI_Aint pos[6] = { offsetof(hllc_value_t,T),offsetof(hllc_value_t,CFL),offsetof(hllc_value_t,h)
			,offsetof(hllc_value_t,print_timer) ,offsetof(hllc_value_t,tau) ,sizeof(hllc_value_t)};
		MPI_Datatype typ[6] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(6, len, pos, typ, &MPI_hllc_value_t);
		MPI_Type_commit(&MPI_hllc_value_t);
	}
	MPI_Bcast(&hllc_cfg, 1, MPI_hllc_value_t, 0, MPI_COMM_WORLD);

#endif // USE_MPI
	while (t < hllc_cfg.T)
	{
		Type time_step = -omp_get_wtime();
		
		if (myid == 0)
		{
			HLLC_STEP(hllc_cfg.tau, grid);

			WRITE_LOG("\n hllc time= " << time_step + omp_get_wtime()<< " c\n");
		}

#ifdef ILLUM
		{
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD); //ждем газодинамического расчёта
			MPI_CalculateIllumAsync(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid);
			MPI_Barrier(MPI_COMM_WORLD); //ждем  расчёта излучения
#else
			CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
#endif // USE_MPI
		}

			if (myid == 0)
			{
				Type tt = -omp_get_wtime();

				CalculateIllumParam(grid_direction, grid);
				
				if (SolveIllumAndHLLC(hllc_cfg.tau, grid) == 1)
				{
					WRITE_LOG("bad illum solve\n");
					//non phys conv value -> continue(надо учесть mpi)
				}

				WRITE_LOG("\n paramIllum time= " << tt + omp_get_wtime() << " c\n");
				//energy.swap(prev_energy);
				//stream.swap(prev_stream);
			}		
#endif //iLLUM

		t += hllc_cfg.tau;
		
		if (myid == 0)
		{
			cur_timer += hllc_cfg.tau;

			GetTimeStep(hllc_cfg, grid); //формирование шага по времени

			if (cur_timer >= hllc_cfg.print_timer)
			{
				//calculate full array
				WriteFileSolution(adress_solve + std::to_string(res_count++), grid);
				
				printf("\n t= %f,  tau= %lf,  res_step= %d\n", t, hllc_cfg.tau, res_count);
				cur_timer = 0;
			}

			time_step += omp_get_wtime();
			WRITE_LOG("\nt= " << t << "; tau= " << hllc_cfg.tau << "; step= " << res_count << " time= " << time_step << " c\n");
		}

#ifdef USE_MPI
		MPI_Bcast(&hllc_cfg, 1, MPI_hllc_value_t, 0, MPI_COMM_WORLD);
#endif

		EXIT_ERR("debug exit\n");
		
	}// while(t < T)

#elif defined ILLUM

#ifdef USE_MPI
	MPI_CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
#else
	CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
#endif
	

	if (myid == 0)
	{
		CalculateIllumParam(grid_direction, grid);
	}
#else
	WRITE_LOG("No solve. Bad config\n");
	RETURN_ERR("No solve. Bad config\n");
#endif
	if (myid == 0)
	{
		full_time += omp_get_wtime();
		printf("Time while(): %f\n", full_time);

#ifdef USE_CUDA // только на главном узле		
		ClearDevice();		
#endif

		WriteFileSolution(adress_solve + std::to_string(res_count), grid);

		WRITE_LOG("End solve module\n");

	}

	return EXIT_SUCCESS;
}

#endif //SOLVE