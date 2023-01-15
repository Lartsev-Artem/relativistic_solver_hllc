
#include "struct_short_characteristics_global_structure.h"
#include "struct_short_characteristics_calculations.h"
#include "struct_short_characteristics_logic_function.h"

#include "../global_value.h"
#include "../file_module/reader_vtk.h"
#include "../file_module/reader_bin.h"
#include "../file_module/reader_txt.h"


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

ShortId n_out;

int pos_x_try;
int posX;
int posX0;
int posOutC;
int posOut;
int posIn;
int posS;
int posRes;

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
	const std::string name_file_id_defining_faces = BASE_ADRESS + "id_defining_faces.bin";
	const std::string name_file_x_defining_faces = BASE_ADRESS + "x_defining_faces.bin";
	const std::string name_file_size = BASE_ADRESS + "Size.txt";  // на ƒќ«јѕ»—№
//--------------------------------создающиес€ файлы----------------------------------------
	const std::string name_file_out_id = BASE_ADRESS + "OutId.bin";
	const std::string name_file_in_id = BASE_ADRESS + "InId.bin";
	const std::string name_file_x = BASE_ADRESS + "X.bin";
	const std::string name_file_x0_loc = BASE_ADRESS + "LocX0.bin";

	const std::string name_file_s = BASE_ADRESS + "S.bin";
	const std::string name_file_count_out_id = BASE_ADRESS + "CountOutId.bin";
	const std::string name_file_res_bound = BASE_ADRESS + "ResBound.bin";

	const std::string name_file_shift_out = BASE_ADRESS + "ShiftOut.bin";
	const std::string name_file_shift_res = BASE_ADRESS + "ShiftRes.bin";
	const std::string name_file_shift_x0 = BASE_ADRESS + "ShiftX0.bin";

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

	vector<Vector3> directions;
	vector<Type> squares;
	
	_clock = -omp_get_wtime();
	{
		std::vector<int> all_pairs_face;
		if (ReadSimpleFileBin(name_file_pairs, all_pairs_face)) return 1;
		InitNodesValue(all_pairs_face, nodes_value);
		all_pairs_face.clear();

		if (ReadCellFaces(name_file_cells, grid)) return 1;
		if (ReadVertex(name_file_vertex, vertexs)) return 1;
		if (ReadNormalFile(name_file_normals, normals)) return 1;	

		if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface)) return 1;

		
	}
	_clock += omp_get_wtime();
	std::cout << "\nReading time of the sphere_direction file: " << _clock << "\n";

	const int count_cells = vertexs.size();
	const int count_directions = directions.size();
	if (b == 0) b = count_directions;

	//vector<int> sorted_id_cell(count_cells * count_directions); // ”пор€доченные индексы €чеек по данному направлению

	vector < vector<int>> sorted_id_cell(count_directions);

	FILE* f;
	

		for (size_t i = a; i < b; i++)
		{
			f = fopen((name_file_graph + std::to_string(i) + ".bin").c_str(), "rb");
			if (!f) RETURN("name_file_graph not open\n")
			sorted_id_cell[i].resize(count_cells);
			fread_unlocked(sorted_id_cell[i].data(), sizeof(int), count_cells, f);
			fclose(f);
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

	std::ofstream ofile;
	ofile.open(BASE_ADRESS + "File_Make_Logs.txt");

	posX = 0;
	posX0 = 0;
	posOutC = 0;
	posOut = 0;
	posIn = 0;
	posS = 0;
	posRes = 0;
	pos_x_try = 0;

	uint64_t shift_out = 0;
	uint64_t shift_x0 = 0;
	uint64_t shift_res = 0;

	_clock = -omp_get_wtime();
	{		
		FILE* file_out_id = fopen(name_file_out_id.c_str(), "wb");
		FILE* file_in_id = fopen(name_file_in_id.c_str(), "wb");
		FILE* file_x = fopen(name_file_x.c_str(), "wb");
		FILE* file_x0_local = fopen(name_file_x0_loc.c_str(), "wb");
		FILE* file_s = fopen(name_file_s.c_str(), "wb");
		FILE* file_count_out_id = fopen(name_file_count_out_id.c_str(), "wb");

		FILE* file_res_bound = fopen(name_file_res_bound.c_str(), "wb");

		FILE* file_shift_out = fopen(name_file_shift_out.c_str(), "wb");
		FILE* file_shift_res = fopen(name_file_shift_res.c_str(), "wb");
		FILE* file_shift_x0 = fopen(name_file_shift_x0.c_str(), "wb");

		if (!file_out_id) RETURN("Outid not open\n")
		if (!file_in_id) RETURN("Inid not open\n")
		if (!file_x) RETURN("Xid not open\n")
		if (!file_x0_local) RETURN("X0id not open\n")
		if (!file_s) RETURN("Sid not open\n")
		if (!file_count_out_id) RETURN("file_count_out_id not open\n")
		
		if (!file_res_bound) RETURN("file_res_bound not open\n")

		if (!file_shift_out) RETURN("file_shift_out not open\n");
		if (!file_shift_res) RETURN("file_shift_res not open\n");
		if (!file_shift_x0) RETURN("file_shift_x0 not open\n");

		Vector3 x;
		Vector3 x0;
		register int num_cell;		
		Vector3 direction;

// массивы записи в файл:
		std::vector<ShortId> vec_n_out(count_cells, 0);
#ifndef ReadNow		
		std::vector<ShortId> vec_out;
		std::vector<ShortId> vec_in;

		std::vector<Type> vec_s;
		std::vector<Vector3> vec_x;
		std::vector<Vector2> vec_x0_loc;

		std::vector<Type> vec_res_bound;

		
#endif // !ReadNow

		/*---------------------------------- далее FOR по направлени€м----------------------------------*/
		for (int num_direction = a; num_direction < b; ++num_direction)
		{
			direction = directions[num_direction];

			num_cur_direction = num_direction;
			cur_direction = direction;
			
#ifndef ReadNow
			vec_out.clear();
			vec_in.clear();

			vec_s.clear();
			vec_x.clear();
			vec_x0_loc.clear();

			const int N = 4 * count_cells;
			vec_out.reserve(N);
			vec_in.reserve(N);

			vec_s.reserve(N);
			vec_x.reserve(N);
			vec_x0_loc.reserve(N);

#endif // ReadNow


			int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   
		
			/*---------------------------------- далее FOR по €чейкам----------------------------------*/
			for (int h = 0; h < count_cells; ++h) {
				//if (count_cells * num_direction + h > sorted_id_cell.size()) RETURN("err size graph\n")
					//num_cell = sorted_id_cell[count_cells * num_direction + h];
					num_cell = sorted_id_cell[num_direction][h];
				
				FindInAndOutFaces(direction, num_cell, normals, face_state);

				vec_n_out[h] = 0;
				n_out = 0;

				for (ShortId num_out_face = 0; num_out_face < 4; ++num_out_face) 
				{
					if (!face_state[num_out_face]) // выход€щие грани
					{  
#ifndef ReadNow
						vec_out.push_back(num_out_face);
						GetNodes(num_cell, grid, num_out_face, vertexs[num_cell],
							face_state, direction, normals, nodes_value,
							vec_res_bound, vec_s, vec_x, vec_x0_loc, vec_in);
#else

						GetNodes(num_cell, grid, num_out_face, vertexs[num_cell],
							face_state, direction, normals, nodes_value,
							file_res_bound, file_s, file_x, file_x0_local, file_in_id);

						fwrite_unlocked(&num_out_face, sizeof(ShortId), 1, file_out_id);
#endif // ReadNow
												
						posOut++;
						vec_n_out[h]++;//n_out++;
					}
				}

				//fwrite_unlocked(&n_out, sizeof(ShortId), 1, file_count_out_id); posOutC++;
			}
			/*---------------------------------- конец FOR по €чейкам----------------------------------*/
#ifndef ReadNow

			fwrite_unlocked(vec_out.data(), sizeof(ShortId), vec_out.size(), file_out_id);
			fwrite_unlocked(vec_in.data(), sizeof(ShortId), vec_in.size(), file_in_id);

			fwrite_unlocked(vec_x0_loc.data(), sizeof(Vector2), vec_x0_loc.size(), file_x0_local);
			fwrite_unlocked(vec_x.data(), sizeof(Vector3), vec_x.size(), file_x);
			fwrite_unlocked(vec_s.data(), sizeof(Type), vec_s.size(), file_s);
			fwrite_unlocked(vec_res_bound.data(), sizeof(Type), vec_res_bound.size(), file_res_bound);

//==========================================================
			
#endif // ReadNow
			fwrite_unlocked(vec_n_out.data(), sizeof(ShortId), count_cells, file_count_out_id);
			fwrite_unlocked(&shift_out, sizeof(uint64_t), 1, file_shift_out);  shift_out = posOut;
			fwrite_unlocked(&shift_res, sizeof(uint64_t), 1, file_shift_res); shift_res = posRes;
			fwrite_unlocked(&shift_x0, sizeof(uint64_t), 1, file_shift_x0); shift_x0 = posX0;

			printf("End direction number #%d\n", num_direction);
			ofile << "End direction number #" << num_direction << '\n';
		}
		/*---------------------------------- конец FOR по направлени€м----------------------------------*/

		
		fclose(file_res_bound);

		fclose(file_out_id);
		fclose(file_in_id);

		fclose(file_x);
		fclose(file_x0_local);
		fclose(file_s);
		fclose(file_count_out_id);

		fclose(file_shift_out);
		fclose(file_shift_res);
		fclose(file_shift_x0);


		posOutC = count_cells * count_directions;
		WriteSize(name_file_size);		
	}


	_clock += omp_get_wtime();
	std::cout << "Time of iter: " << _clock << '\n';

	ofile << "Time of iter: " << _clock << '\n';
	ofile.close();
#endif // USE_VTK

	printf("End Programm\n");

	return 0;
}


