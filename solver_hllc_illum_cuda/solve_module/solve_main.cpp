#ifdef SOLVE
#include "solve_global_struct.h"
#include "../global_def.h"


#include "../global_headers.h"
#include "../global_value.h"

#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_bin.h"

#include "../utils/grid_geometry/geometry_solve.h"
#include "../utils/grid_geometry/geometry_data.h"

#include "solve_utils.h"

#include"illum/illum_utils.h"
#include "hllc/hllc_utils.h"  //сделать общий для rhllc и hllc

#include "../cuda/cuda_solve.h"

#include "MPI_utils.h"

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

int RunSolveModule(int argc, char* argv[], const std::string& name_file_settings)
{
	int np, myid;
	MPI_GET_INF(np, myid);

	int res_count = 0; // счётчик решений	
	if (argc > 1)
	{
		glb_files.name_file_settings = argv[1];
		if (argc > 2) res_count = std::stoi(argv[2]);
	}
	else
	{
		glb_files.name_file_settings = name_file_settings;
	}

	if (ReadStartSettings(glb_files, solve_mode))
	{
		RETURN_ERR("Error reading solve settings\n");
	}

	WRITE_LOG_ERR("start solve module " << myid << "\n");

#if (defined HLLC || defined RHLLC) &&  NUMBER_OF_MEASUREMENTS < 3
	StartLowDimensionTask(glb_files.glb_files.base_adress);
	WRITE_LOG_ERR("end low dimension task\n");
	return 0;
#endif

#ifdef RUN_TEST
	TestDivStream(glb_files.base_adress);
	return 0;
#endif // RUN_TEST

	grid_t grid;

	if (ReadGeometryGrid(glb_files.name_file_geometry_cells, glb_files.name_file_geometry_faces, grid))
	{
#pragma warning "grid on each nodes!"
		WRITE_LOG_ERR("Error reading grid, try read parts geo\n");
	}
	WRITE_LOG("Reading geometry grid\n");

#ifdef USE_MPI
	Init_MPI();
	//MPI_Bcast(&grid.size, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif // USE_MPI

	double _clock;

	//------------------------------ Illume section-----------------------------
	//std::vector<Type> Illum; // отдельно от сетки из-за расчета на cuda
	//std::vector<Type> int_scattering; // отдельно от сетки из-за расчета на cuda

#ifdef ILLUM
	grid_directions_t grid_direction;
	if (ReadSphereDirectionDecartToSpherical(glb_files.name_file_sphere_direction, grid_direction)) RETURN_ERR("Error reading the sphere direction\n");

	WRITE_LOG("reading sphere_direction\n");

	std::vector<BasePointTetra> vec_x;
	std::vector<Type> vec_res_bound;
	std::vector<int> pairs;

	const int count_dir = grid_direction.size;
	std::vector < std::vector<int>> face_states; //битовое поле: 0=> выходящая грань,  1=> входящая   
	std::vector < std::vector<cell_local>> vec_x0;
	std::vector < std::vector<int>> sorted_id_cell; 	// Упорядоченные индексы ячеек по данному направлению	

	_clock = -omp_get_wtime();
	if (ReadIllumGeometry(count_dir, glb_files, vec_x, face_states, vec_x0, sorted_id_cell, vec_res_bound))
	{
		RETURN_ERR("Error fast reading the data grid vtk\n");
	}
	_clock += omp_get_wtime();

	WRITE_LOG("Reading time of the Illum Geometry file : " << _clock << "\n");

	ReadSimpleFileBin(glb_files.name_file_neib, pairs);

	WRITE_LOG_MPI("start Init illum\n", myid);

#ifdef USE_MPI
	InitSendDispIllumArray(myid, np, grid_direction.size, grid.size);	
	InitPhysOmpMpi(grid.size);
#endif

	{
#ifdef USE_CUDA 	
		SetDevice(0);
		InitDevice(grid_direction, grid, GetDisp(myid), GetDisp(myid) + GetSend(myid));
#else
		grid.InitMemory(grid_direction.size);		
#endif			
		if (myid == 0) InitIllum(glb_files.base_adress, grid);
	}

	MPI_INIT(myid, np, grid_direction.size, grid);
	WRITE_LOG_MPI("Init illum\n", myid);
#endif //ILLUM
	//------------------------------------------------------------------------------------------------------------

	//---------------------------------------Solve section--------------------------------------------------------------	


#if defined HLLC || defined RHLLC
	hllc_value_t hllc_cfg;

	if (HLLC_INIT(glb_files.name_file_hllc_set, hllc_cfg, glb_files.hllc_init_value, grid.cells)) RETURN_ERR("Bad init hllc\n");  //Начальные данные для HLLC

#ifdef USE_MPI	
	MPI_Bcast(&hllc_cfg, 1, MPI_hllc_value_t, 0, MPI_COMM_WORLD);
#endif // USE_MPI

	if (myid == 0)
	{
		WriteFileSolution(glb_files.solve_adress + std::to_string(res_count++), grid); //печать начальной сетки
		WRITE_LOG_ERR("Start main task\n");
	}

#ifdef DEBUG	
	CREATE_STRUCT(timer_t, Type, full_time, hllc_time, illum_time, illum_param_time, solve_hllc_time, step) timer;
	timer.full_time = -omp_get_wtime();
#endif

	MPI_Barrier(MPI_COMM_WORLD); // ждем пока все считается с дисков

	int skip_count = 0;
	const int skip_size = solve_mode.class_vtk != e_class_grid_static_illum ? 10 : 1e10; //10-стационарное излучение (пропускаем всегда)

	Type t = 0.0;
	Type cur_timer = 0;

	while (t < hllc_cfg.T)
	{
#ifdef DEBUG		
		timer.step = -omp_get_wtime();
		timer.hllc_time = -omp_get_wtime();
#endif

#ifdef RHLLC_MPI
		RHLLC_3d_MPI(hllc_cfg.tau, grid);
#else
#ifdef ILLUM
		MPI_RHLLC_3d(myid, hllc_cfg.tau, grid);
#else
		if (myid == 0)
		{
			HLLC_STEP(hllc_cfg.tau, grid);
		}
#endif
#endif //RHLLC_MPI

#ifdef DEBUG
		if (myid == 0)
		{
			WRITE_LOG("\n hllc time= " << timer.step + omp_get_wtime() << " c\n");
		}
		timer.hllc_time += omp_get_wtime();
#endif

#ifdef ILLUM
#ifdef DEBUG
		timer.illum_time = -omp_get_wtime();
#endif

		if (skip_count++ % skip_size == 0)
		{
#ifdef USE_MPI			
			MPI_CalculateIllumAsync(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid);
#ifndef USE_CUDA

			WRITE_LOG_MPI("CalculateIllumParam 1\n", myid);
			CalculateIllumParam(grid_direction, grid);
			WRITE_LOG_MPI("CalculateIllumParam 2\n", myid);
#endif		
#else
			CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid, Illum, int_scattering);
#endif // USE_MPI
		}
#ifdef DEBUG
		timer.illum_time += omp_get_wtime();
		timer.solve_hllc_time = -omp_get_wtime();
#endif

		
#ifndef RHLLC_MPI
		if (myid == 0)
#else
		{			
			if (SolveIllumAndHLLC(hllc_cfg.tau, grid) == 1)
			{
				WRITE_LOG_ERR("bad illum solve\n");
				D_LD;
				//non phys conv value -> continue(надо учесть mpi)
			}			

			//energy.swap(prev_energy);
			//stream.swap(prev_stream);
		}
#endif
#endif //iLLUM
#ifdef DEBUG
		timer.solve_hllc_time += omp_get_wtime();
#endif

		t += hllc_cfg.tau;
		cur_timer += hllc_cfg.tau;

		//GetTimeStep(hllc_cfg, grid); //формирование шага по времени

		if (cur_timer >= hllc_cfg.print_timer)
		{
#ifdef RHLLC_MPI

			if (myid != 0)
			{
				MPI_Send(grid.cells.data() + mpi_conf[myid].left, mpi_conf[myid].right - mpi_conf[myid].left, MPI_flux_elem_t, 0, myid, MPI_COMM_WORLD);
			}
			else
			{
				for (int tag = 1; tag < np; tag++)
				{
					MPI_Recv(grid.cells.data() + mpi_conf[tag].left, mpi_conf[tag].right - mpi_conf[tag].left,
						MPI_flux_elem_t, tag, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
#endif
			if (myid == 0)
			{
				GetTimeStep(hllc_cfg, grid); //формирование шага по времени

				WriteFileSolution(glb_files.solve_adress + std::to_string(res_count++), grid);

				WRITE_LOG_ERR("\nt= " << t << "; tau= " << hllc_cfg.tau << "; step= " << res_count << " time_step= " << timer.step + omp_get_wtime() << " c, time" << timer.full_time + omp_get_wtime() << " c\n");
				printf("\n t= %f,  tau= %lf,  res_step= %d\n", t, hllc_cfg.tau, res_count);

			}
			cur_timer = 0;						
		}

		timer.step += omp_get_wtime();

#ifdef USE_MPI		
		MPI_Bcast(&hllc_cfg, 1, MPI_hllc_value_t, 0, MPI_COMM_WORLD);
#endif

		WRITE_LOG_MPI("t= " << t << ", hllc= " << timer.hllc_time << ", illum " << timer.illum_time
			<< ", illum_param " << timer.illum_param_time
			<< ", solve " << timer.solve_hllc_time
			<< ", step " << timer.step << "\n\n", myid);

	}// while(t < T)

#elif defined ILLUM

#ifdef USE_MPI
	MPI_CalculateIllumAsync(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid);
	//MPI_CalculateIllum(grid_direction, face_states, pairs, vec_x0, vec_x, sorted_id_cell, grid);
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

#ifdef DEBUG
	timer.full_time += omp_get_wtime();
	WRITE_LOG_MPI("full time= " << timer.full_time << '\n', myid);
	printf("Time while(): %f\n", timer.full_time);
#endif

#ifdef USE_CUDA // только на главном узле		
	ClearDevice();
#endif

	if (myid == 0)
	{		
		WriteFileSolution(glb_files.solve_adress + std::to_string(res_count), grid);
		WRITE_LOG_ERR("End solve module\n");
	}
	return EXIT_SUCCESS;
}

#endif //SOLVE