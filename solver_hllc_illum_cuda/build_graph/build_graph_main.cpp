#if defined BUILD

#include "../global_headers.h"
#include "../global_def.h"
#include "../global_value.h"

#include "build_graph_prj_config.h"
#include "build_graph_structures.h"

#include "build_graph_read_write.h"
#include "build_graph_calculation.h"

#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"
#include "../file_module/writer_bin.h"


#ifdef USE_MPI 
#define GetTime MPI_Wtime()
#else
#define GetTime omp_get_wtime()
#endif

std::vector<int> id_try_surface;
std::vector<Type> dist_try_surface;
std::vector<Vector3> x_try_surface;
TrySolve buf_try;

uint64_t id_try_size;
uint64_t dist_try_size;
uint64_t x_try_size;

int RunBuildModule(const std::string& name_file_settings)
{
#ifdef  USE_MPI
	int np, myid;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#else
	int myid = 0;
#endif //  USE_MPI

	int class_vtk = 0;
	std::string name_file_vtk;
	std::string name_file_sphere_direction; // = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph;
	std::string str;
	int iter;

	if (ReadStartSettings(name_file_settings, class_vtk, name_file_vtk, name_file_sphere_direction, name_file_graph, str, BASE_ADRESS, str, iter))
	{
		RETURN_ERR("Error reading build settings\n");
	}
	std::remove((BASE_ADRESS + "File_Logs.txt").c_str());

	std::string name_file_normals = BASE_ADRESS + "normals.bin";
	std::string name_file_pairs = BASE_ADRESS + "pairs.bin";
	std::string name_file_inner_boundary = BASE_ADRESS + "inner_bound.txt";
	std::string name_file_init_boundary = BASE_ADRESS + "init_boundary.txt";
	std::string name_file_face_and_id = BASE_ADRESS + "faceId.txt";

	std::string name_file_squares = BASE_ADRESS + "squares.bin";
	std::string name_file_volume = BASE_ADRESS + "volume.bin";
	std::string name_file_centers = BASE_ADRESS + "centers.bin";

	std::string name_file_centers_face = BASE_ADRESS + "center_face.bin";

	remove((std::string(BASE_ADRESS) + "File_Logs.txt").c_str());

	double t = 0;
#ifdef WriteFiles
	if (myid == 0)
	{
		t = -GetTime;

		BuildSetForClaster(name_file_vtk, name_file_pairs, name_file_init_boundary, name_file_normals,
			name_file_inner_boundary, name_file_face_and_id,
			name_file_squares, name_file_volume, name_file_centers, name_file_centers_face);

		t += GetTime;

		printf("Time for build start set: %lf\n", t);
	}

#endif //WriteFiles

#ifdef  ReadFiles
	{
		if (myid == 0)  printf("OMP num_threads: %d\n", omp_get_num_threads());

		t = -GetTime;

		std::vector<IntId> all_pairs_id;
		std::set<IntId> set_inner_boundary_cells;
		std::vector<Normals> normals;
		std::map<IntId, FaceCell> inner_faces;
		std::vector<Vector3> directions;

		if (ReadSimpleFileBin(name_file_pairs, all_pairs_id)) RETURN_ERR("Error reading file neighbours\n");

		if (ReadInnerCellBoundary(name_file_inner_boundary, set_inner_boundary_cells)) RETURN_ERR("Error reading file bound_cell");
		std::cout << "Inner boundary has " << set_inner_boundary_cells.size() << "faces\n";

		if (ReadNormalFile(name_file_normals, normals)) RETURN_ERR("Error reading file normals\n");

		if (ReadInnerCellOfSphereAndId(name_file_face_and_id, inner_faces)) RETURN_ERR("Error reading file inner_faces\n");

		if (ReadSphereDirectionDecart(name_file_sphere_direction, directions)) RETURN_ERR("Error reading file directions\n");

		if (myid == 0)
		{
			t += GetTime;
			printf("Time reading in main procces: %lf\n", t);
			t = -GetTime;
		}


		const int num_cells = all_pairs_id.size() / 4;
		const int count_direction = directions.size();

		std::vector<IntId> count_in_face(num_cells, 0);
		std::vector<IntId> count_knew_face(num_cells, 0);
		std::vector<IntId> graph(num_cells, 0);

		std::set<IntId> inner_part;
		std::set<IntId> outter_part;

		std::vector<State> faces_state;
		Vector3 direction;

		bool flag = true;
		int count = 0;

#ifdef USE_OMP
		std::vector<IntId> next_step_el_OMP;
		std::list<IntId> cur_el_OMP;
#else
		std::vector<IntId> cur_el;
		std::set<IntId> next_step_el;
#endif // USE_OMP


		//	std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces" + ".bin").c_str(), "wb"), fclose);
		//	if (!file_id) { printf("file_id is not opened for writing\n"); return 1; }

		//	std::unique_ptr<FILE, int(*)(FILE*)> file_dist(fopen((std::string(BASE_ADRESS) + "dist_defining_faces" + ".bin").c_str(), "wb"), fclose);
		//	if (!file_dist) { printf("file_dist is not opened for writing\n"); return 1; }


		//	std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((std::string(BASE_ADRESS) + "x_defining_faces" + ".bin").c_str(), "wb"), fclose);
		//	if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }

		//	std::unique_ptr<FILE, int(*)(FILE*)> file_shift_try(fopen((std::string(BASE_ADRESS) + "ShiftTry" + ".bin").c_str(), "wb"), fclose);
		//	if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }


//========================================================================================//

#if defined ONLY_ONE_DIRECTION
		for (int cur_direction = 0; cur_direction < 1; cur_direction++)
#else
#ifdef USE_MPI		
		for (int cur_direction = myid; cur_direction < count_direction; cur_direction += np)
#else
		for (int cur_direction = 0; cur_direction < count_direction; ++cur_direction)
#endif
#endif //ONLY_ONE_DIRECTION

		{
			if (myid == 0)
			{
				WRITE_LOG("Diretion #: " << cur_direction << '\n');
			}

			flag = true;
			InitFacesState(all_pairs_id, faces_state, inner_faces);
			direction = directions[cur_direction];  // 13 -- error direction for test.vtk grid
			
			FractionInnerBoundary(direction, normals, inner_faces, set_inner_boundary_cells, inner_part, outter_part);

			//id_try_surface.clear();
			//dist_try_surface.clear();
			//x_try_surface.clear();

			//x_try_surface.reserve(outter_part.size() * 3);
			//id_try_surface.reserve(outter_part.size() * 3);
			//dist_try_surface.reserve(outter_part.size() * 3);

			//-------------------------------------
#ifdef USE_OMP
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el_OMP);
#else
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el);
#endif // USE_OMP
			//-------------------------------------

			int count_graph = 0; // число €чеек вошедших в граф 
			graph.assign(num_cells, -1); // дл€ отлавливани€ ошибочных направлений
			bool try_restart = true;

#ifdef USE_OMP
			while (next_step_el_OMP.size() && flag)
#else
			while (next_step_el.size() && flag)
#endif // USE_OMP
			{

#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY	
				IntId id_cell = FindCurCellWithHole(next_step_el_OMP, count_in_face, count_knew_face, cur_el_OMP, inner_part, outter_part, \
					inner_faces, all_pairs_id, direction, normals);
#else
				IntId id_cell = FindCurCell(next_step_el_OMP, count_in_face, count_knew_face, cur_el_OMP);
#endif //GRID_WITH_INNER_BOUNDARY

#else  //no use omp

#ifdef GRID_WITH_INNER_BOUNDARY	
				IntId id_cell = FindCurCellWithHole(next_step_el, count_in_face, count_knew_face, cur_el, inner_part, outter_part, \
					inner_faces, all_pairs_id, direction, normals);
#else
				IntId id_cell = FindCurCell(next_step_el, count_in_face, count_knew_face, cur_el);
#endif //GRID_WITH_INNER_BOUNDARY
#endif // USE_OMP

				if (id_cell == -1)
				{
#ifdef USE_MPI
					std::cout << "Error proc: " << myid << '\n';
#endif // USE_MPI
					//WriteFileGraph(cur_direction, name_file_graph, graph);
					std::cout << "Error num_direction: " << cur_direction << '\n';
					std::cout << "Error.Complete " << count_graph << " cells\n";
					WRITE_LOG("Error num_direction: " << cur_direction << "\nError.Complete " << count_graph << " cells\n");

					flag = false; // break;
				}

#ifdef USE_OMP

				NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el_OMP, next_step_el_OMP);

				for (auto el : cur_el_OMP) {
					graph[count_graph] = el;
					count_graph++;
				}
#else
				NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el, next_step_el);
				for (auto el : cur_el)
				{
					graph[count_graph] = el;
					count_graph++;
				}
#endif // USE_OMP

				if (count_graph < graph.size() && next_step_el.size() == 0 && try_restart)
				{
					try_restart = !try_restart;

					flag = true;

#ifdef GRID_WITH_INNER_BOUNDARY
					if (outter_part.size() != 0)
					{
						cur_el.clear();
						next_step_el.clear();
						next_step_el.insert(outter_part.begin(), outter_part.end());
						printf("Error. try_short_restart %d\n", cur_direction);

						WRITE_LOG("Error. try_short_restart" << cur_direction << "\n");

						//пытатьс€ определ€ть внутреннюю границу до последнего
						static int cc = 0;
						if (cc++ < 10)
							try_restart = !try_restart;

						continue;
					}
#endif //GRID_WITH_INNER_BOUNDARY

					printf("\n        Error. try_restart %d            \n", cur_direction);
					//ofile << "----------   Error. try_restart   ----------" << cur_direction << "\n";

					cur_el.clear();
					next_step_el.clear();

					std::cout << outter_part.size() << '\n';
					for (size_t i = 0; i < num_cells; i++)
					{
						if (count_in_face[i] > count_knew_face[i])
						{
							next_step_el.emplace(i);  // €чейка была изменена, проверить ее готовность на следующем шаге
						}

					}
				}

			}//while

			try_restart = !try_restart;
			if (count_graph < graph.size())
			{
				printf("-------------------------Error size graph-------------------------------\n");
				WRITE_LOG("----------  Error size graph   ----------(" << cur_direction << ") Size: " << count_graph << "\n");
			}

			//fwrite_unlocked(&id_try_size, sizeof(int), 1, file_shift_try.get());

			if (WriteSimpleFileBin(name_file_graph + std::to_string(cur_direction) + ".bin", graph))
			{
				RETURN_ERR("file_graph is not opened for writing\n");
			}

			printf("Id_proc: %d. Graph construction in the direction %d is completed\n", myid, cur_direction);
		}

		//	fclose(file_id.get());
		//  fclose(file_dist.get());
		//	fclose(file_x.get());
		//	fclose(file_shift_try.get());

		std::ofstream ofile2(std::string(BASE_ADRESS) + "Size.txt");
		ofile2 << id_try_size << '\n';
		// "id_try_size: " << id_try_size;
		/*ofile2 << "\ndist_try_size: " << dist_try_size;
		ofile2 << "\nx_try_size: " << x_try_size << '\n';*/
		ofile2.close();

		if (myid == 0)
		{
			t += GetTime;
			printf("Full time: %lf\n", t);
		}
	}
#endif //  ReadFiles

	return EXIT_SUCCESS;
}

#endif //BUILD