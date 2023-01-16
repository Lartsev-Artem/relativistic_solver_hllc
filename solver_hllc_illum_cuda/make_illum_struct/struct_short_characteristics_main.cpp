#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_calculations.h"
#include "struct_short_characteristics_logic_function.h"

#include "../file_module/reader_vtk.h"
#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_bin.h"

#include "../utils/grid_geometry/geometry_solve.h"
#if 0
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int num_cur_direction; // номер текущего направлени€
Vector3 cur_direction;

std::vector<Vector3> x_try_surface;
std::vector<int> id_try_surface;

struct make_illum_val pos_cnt;
#endif

int RunMakeBuild(std::string name_file_settings, int a, int b)
{
	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string out_file_grid_vtk;
	std::string str; int class_vtk;

	if (ReadStartSettings(name_file_settings, class_vtk, name_file_vtk, name_file_sphere_direction, str, BASE_ADRESS, str))
	{
		RETURN_ERR("Error reading the start settings\n");
	}

	//-----------файлы с данными сетки. ѕостроены здесь на метки USE_VTK-----------------------
	const std::string name_file_normals = BASE_ADRESS + "normals.bin";
	const std::string name_file_cells = BASE_ADRESS + "grid.bin";
	const std::string name_file_vertex = BASE_ADRESS + "vertex.bin";
	const std::string name_file_pairs = BASE_ADRESS + "pairs.bin";
	//-------------------читающиес€ файлы, построенные в build_graph---------------------------
	const std::string name_file_graph = BASE_ADRESS + "graph";
	const std::string name_file_id_defining_faces = BASE_ADRESS + "pairs.bin";
	const std::string name_file_x_defining_faces = BASE_ADRESS + "x_defining_faces.bin";
	const std::string name_file_size = BASE_ADRESS + "Size.txt";  // на ƒќ«јѕ»—№
	//--------------------------------создающиес€ файлы----------------------------------------
	const std::string name_file_state_face = BASE_ADRESS + "state_face";
	const std::string name_file_x = BASE_ADRESS + "X";
	const std::string name_file_x0_loc = BASE_ADRESS + "LocX0";
	const std::string name_file_res_bound = BASE_ADRESS + "ResBound";

	Type _clock;
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

#ifdef ONLY_BUILD_DATA
	return 0;
#else
	// make
	std::vector<Face> grid;
	std::vector<Eigen::Matrix4d> vertexs;
	std::vector<Normals> normals;
	std::vector<cell> nodes_value;

	std::vector<direction_s> directions;
	Type square_surface;

	_clock = -omp_get_wtime();
	{
		std::vector<int> all_pairs_face;
		if (ReadSimpleFileBin(name_file_pairs, all_pairs_face)) return 1;
		InitNodesValue(all_pairs_face, nodes_value);
		all_pairs_face.clear();

		if (ReadCellFaces(name_file_cells, grid)) return 1;
		if (ReadSimpleFileBin(name_file_vertex, vertexs)) return 1;
		if (ReadNormalFile(name_file_normals, normals)) return 1;

		if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, square_surface)) return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\nReading time of the sphere_direction file: " << _clock << "\n";

	const int count_cells = vertexs.size();
	const int count_directions = directions.size();
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

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

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
		direction = directions[num_direction].dir;

		vec_x0.clear();
		vec_x0.reserve(base * count_cells);

		/*---------------------------------- далее FOR по €чейкам----------------------------------*/
		for (int h = 0; h < count_cells; ++h)
		{
			num_cell = sorted_id_cell[num_direction][h];

			int face_state = 0;
			GetInAndOutFaces(direction, normals[num_cell], face_state);
			face_states[num_cell] = face_state;

			for (ShortId num_out_face = 0; num_out_face < 4; ++num_out_face)
			{
				if (!check_bit(face_state, num_out_face)) // выход€щие грани
				{
					GetNodes(num_cell, grid, num_out_face, vertexs[num_cell],
						face_state, direction, normals, nodes_value,
						vec_x[num_cell], vec_res_bound, vec_x0);
				}
			}
		}
		/*---------------------------------- конец FOR по €чейкам----------------------------------*/

		if (WriteSimpleFileBin(name_file_state_face + std::to_string(num_direction) + ".bin", face_states)) RETURN_ERR("Error vec_x");
		if (WriteSimpleFileBin(name_file_x0_loc + std::to_string(num_direction) + ".bin", vec_x0)) RETURN_ERR("Error vec_x");

		printf("End direction number #%d\n", num_direction);
		WRITE_LOG("End direction number #" << num_direction << '\n');
	}
	/*---------------------------------- конец FOR по направлени€м----------------------------------*/

	if (WriteSimpleFileBin(name_file_x, vec_x)) RETURN_ERR("Error vec_x");

	_clock += omp_get_wtime();
	std::cout << "Time of make: " << _clock << '\n';
	WRITE_LOG("Time of make: " << _clock << '\n');

#endif // ONLY_BUILD_DATA

	printf("End Programm\n");

	return 0;
}


