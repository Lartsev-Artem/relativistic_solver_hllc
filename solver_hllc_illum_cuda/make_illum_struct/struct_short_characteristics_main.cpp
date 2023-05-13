#ifdef MAKE
#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_calculations.h"
#include "struct_short_characteristics_logic_function.h"

#include "../file_module/reader_vtk.h"
#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_bin.h"

#include "../utils/grid_geometry/geometry_solve.h"

int RunMakeModule(std::string name_file_settings, int a, int b)
{
	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string name_file_graph;
	std::string name_file_Illum_geo;
	std::string out_file_grid_vtk;
	std::string str; int class_vtk; int iter;

	if (ReadStartSettings(name_file_settings, class_vtk, name_file_vtk, name_file_sphere_direction, name_file_graph, name_file_Illum_geo, glb_files.base_adress, str, iter))
	{
		RETURN_ERR("Error reading the start settings\n");
	}

	//-----------файлы с данными сетки. ѕостроены здесь на метки USE_VTK-----------------------
	const std::string name_file_normals = glb_files.base_adress + "normals.bin";
	const std::string name_file_cells = glb_files.base_adress + "grid.bin";
	const std::string name_file_vertex = glb_files.base_adress + "vertex.bin";
	const std::string name_file_pairs = glb_files.base_adress + "pairs.bin";
	//-------------------читающиес€ файлы, построенные в build_graph---------------------------	
	const std::string name_file_id_defining_faces = glb_files.base_adress + "pairs.bin";
	const std::string name_file_x_defining_faces = glb_files.base_adress + "x_defining_faces.bin";
	const std::string name_file_size = glb_files.base_adress + "Size.txt";  // на ƒќ«јѕ»—№
	//--------------------------------создающиес€ файлы----------------------------------------
	const std::string name_file_state_face = name_file_Illum_geo + "state_face";
	const std::string name_file_x = name_file_Illum_geo + "X.bin";
	const std::string name_file_x0_loc = name_file_Illum_geo + "LocX0";
	const std::string name_file_res_bound = name_file_Illum_geo + "ResBound";
	
#ifdef USE_VTK

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_vtk, unstructured_grid))
	{
		RETURN_ERR("Error reading the file vtk\n");
	}

	WriteCellFaces(name_file_cells, unstructured_grid);
	WriteVertex(name_file_vertex, unstructured_grid);
#endif

#ifdef ONLY_GEO_DATA
	return 0;
#else
	Type _clock = 0;
	// make
	std::vector<Face> grid;
	std::vector<Matrix4> vertexs;
	std::vector<Normals> normals;
	std::vector<int> all_pairs_face;	

	grid_directions_t grid_direction;
	
	_clock = -omp_get_wtime();
	{		
		if (ReadSimpleFileBin(name_file_pairs, all_pairs_face)) return 1;	
		if (ReadCellFaces(name_file_cells, grid)) return 1;
		if (ReadSimpleFileBin(name_file_vertex, vertexs)) return 1;
		if (ReadNormalFile(name_file_normals, normals)) return 1;

		if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, grid_direction)) return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\nReading time of the sphere_direction file: " << _clock << "\n";

	const int count_cells = vertexs.size();
	const int count_directions = grid_direction.size;
	if (b == 0) b = count_directions;

	std::vector < std::vector<int>> sorted_id_cell(count_directions);

	for (size_t i = a; i < b; i++)
	{
		ReadSimpleFileBin(name_file_graph + std::to_string(i) + ".bin", sorted_id_cell[i]);
	}

	{
		/*	int n;
			std::ifstream ifile(name_file_size);
			if (!ifile.is_open()) RETURN("Error read files size\n")
			ifile >> n;
			ifile.close();

			FILE* f;
			f = fopen(name_file_id_defining_faces.c_str(), "rb");
			if (!f) RETURN("name_file_id_defining_faces not open\n")
				id_try_surface.resize(n);
			fread_unlocked(id_try_surface.data(), sizeof(int), n, f);
			fclose(f);

			f = fopen(name_file_x_defining_faces.c_str(), "rb");
			if (!f) RETURN("name_file_x_defining_faces not open\n")
				x_try_surface.resize(n);
			fread_unlocked(x_try_surface.data(), sizeof(Vector3), n, f);
			fclose(f);*/

	}
	
	std::vector<BasePointTetra> vec_x(count_cells);
	MakeArrayX(vertexs, vec_x);

	_clock = -omp_get_wtime();

	register int num_cell;
	Vector3 direction;

	// массивы записи в файл:
	std::vector<int> face_states(count_cells, 0); //битовое поле: 0=> выход€ща€ грань,  1=> вход€ща€   
	std::vector<cell_local> vec_x0;
	std::vector<Type> vec_res_bound;

	/*---------------------------------- далее FOR по направлени€м----------------------------------*/
	for (int num_direction = a; num_direction < b; ++num_direction)
	{
		direction = grid_direction.directions[num_direction].dir;

		vec_x0.clear();
		vec_x0.reserve(base * count_cells);

		/*---------------------------------- далее FOR по €чейкам----------------------------------*/
		for (int h = 0; h < count_cells; ++h)
		{
			num_cell = sorted_id_cell[num_direction][h];
			
			int face_state = 0;
			GetInAndOutFaces(direction, normals[num_cell], face_state);
			face_states[num_cell] = face_state;


			for (ShortId num_out_face = 0; num_out_face < base; ++num_out_face)
			{
				if (!check_bit(face_state, num_out_face)) // выход€щие грани
				{
					GetNodes(num_cell, grid, num_out_face, vertexs[num_cell],
						face_state, direction, normals, all_pairs_face,
						vec_x[num_cell], vec_res_bound, vec_x0);
				}
			}
		}
		/*---------------------------------- конец FOR по €чейкам----------------------------------*/

		if (WriteSimpleFileBin(name_file_state_face + std::to_string(num_direction) + ".bin", face_states)) RETURN_ERR("Error vec_x");
		if (WriteSimpleFileBin(name_file_x0_loc + std::to_string(num_direction) + ".bin", vec_x0)) RETURN_ERR("Error vec_x");
		
		WRITE_LOG("End direction number #" << num_direction << '\n');
	}
	/*---------------------------------- конец FOR по направлени€м----------------------------------*/

	if (WriteSimpleFileBin(name_file_x, vec_x)) RETURN_ERR("Error vec_x");
	if (WriteSimpleFileBin(name_file_res_bound, vec_res_bound)) RETURN_ERR("Error vec_res_bound");

	_clock += omp_get_wtime();
	std::cout << "Time of make: " << _clock << '\n';
	WRITE_LOG("Time of make: " << _clock << '\n');

#endif // ONLY_BUILD_DATA

	printf("End Programm\n");

	return 0;
}

#endif //MAKE