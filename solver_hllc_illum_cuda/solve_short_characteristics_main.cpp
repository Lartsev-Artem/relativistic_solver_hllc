#include "solve_short_characteristics_main.h"
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Matrix3	straight_face_inverse;  // 3 узла интерпол€ции
Matrix3 inclined_face_inverse;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
#ifdef USE_VTK
vtkDataArray* density;
vtkDataArray* absorp_coef;
vtkDataArray* rad_en_loose_rate;
#else
std::vector<Type> density;
std::vector<Type> absorp_coef;
std::vector<Type> rad_en_loose_rate;
std::vector<Type> pressure;
std::vector<Vector3> velocity;
std::vector<Type> e_substance;

#endif
std::vector<VectorX> U_full;

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра
size_t class_file_vtk;

int size_grid;
std::vector<Type> res_inner_bound;  // значение на внутренней границе
std::vector<int> id_try_surface;    // id определ€ющих €чеек внутренней границе

//
//int main2(int argc, char* argv[])
//{
//	std::string name_file_settings = "";
//	int max_number_of_iter = 1;
//	Type accuracy = 1e-5;
//	bool use_cuda = true;
//	const int cuda_mod = 1; //1 - все массивы, 0 - min
//
//	if (argc <= 1)
//		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_solve.txt";
//	else
//		name_file_settings = argv[1];
//	if (argc > 2)
//		max_number_of_iter = std::stoi(argv[2]);
//
//	std::cout << "Max_number_of_iter= " << max_number_of_iter << '\n';
//
//
//	std::string name_file_vtk;
//	std::string name_file_sphere_direction;
//	std::string out_file_grid_vtk;
//
//std:string main_dir;
//
//	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk,
//		main_dir)) {
//		std::cout << "Error reading the start settings\n";
//		return 1;
//	}
//	std::string name_file_graph = main_dir + "graph";
//	std::string name_file_in_faces = main_dir + "InId";
//	std::string name_file_out_faces = main_dir + "OutId";
//	std::string name_file_count_out_faces = main_dir + "CountOutId";
//	std::string name_file_local_x0 = main_dir + "LocX0";
//	std::string name_file_x = main_dir + "X";
//	std::string name_file_s = main_dir + "S";
//
//	std::string name_file_id_neighbors = main_dir + "pairs";
//	std::string name_file_centers = main_dir + "centers";
//
//	std::string name_file_dist_try = main_dir + "dist_defining_faces";
//	std::string name_file_id_try = main_dir + "id_defining_faces";
//	std::string name_file_res = main_dir + "ResBound";
//	std::string name_file_sizes = main_dir + "Size";
//
//	std::string name_file_shift_out = main_dir + "ShiftOut";
//	std::string name_file_shift_res = main_dir + "ShiftRes";
//	std::string name_file_shift_x0 = main_dir + "ShiftX0";
//	std::string name_file_shift_try = main_dir + "ShiftTry";
//
//	std::string name_file_normals = main_dir + "normals.bin";
//	std::string name_file_squares = main_dir + "squares.bin";
//	std::string name_file_volume = main_dir + "volume.bin";
//
//	Type _clock;
//
//#ifdef USE_VTK
//
//#ifdef ReBuildSolve
//	if (ReBuildDataArray(class_file_vtk, name_file_vtk, out_file_grid_vtk, main_dir, true)) printf("Error ReBuild\n");
//	return 0;
//#endif
//#ifdef ReWriteSolve
//	/*vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
//		vtkSmartPointer<vtkUnstructuredGrid>::New();
//	 _clock = -omp_get_wtime();
//	if (ReadDataArray(class_file_vtk, name_file_vtk, unstructured_grid, size_grid, density, absorp_coef, rad_en_loose_rate, true)) {
//		std::cout << "Error reading the data array vtk\n";
//		return 1;
//	}
//	_clock += omp_get_wtime();
//	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";*/
//
//	if (ReWriteDataArray(class_file_vtk, name_file_vtk, main_dir, size_grid, true)) {
//		std::cout << "Error rewrite the data array vtk\n";
//		return 1;
//	}
//	return 0;
//#endif
//
//
//#else
//	 _clock = -omp_get_wtime();
//	if (ReadDataArray(class_file_vtk, main_dir, density, absorp_coef, rad_en_loose_rate, size_grid, true)) {
//		std::cout << "Error reading the data array vtk\n";
//		return 1;
//	}
//	_clock += omp_get_wtime();
//	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";
//#endif // USE_VTK
//
//	vector<Vector3> directions;
//	vector<Type> squares;
//
//	_clock = -omp_get_wtime();
//	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface)) {
//		std::cout << "Error reading the sphere direction\n";
//		return 1;
//	}
//	_clock += omp_get_wtime();
//	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";
//
//	const int count_directions = directions.size();
//	const int count_cells = size_grid;
//
//	std::vector<int> neighbours_id_faces;
//	std::vector<Vector3> nodes_value;
//	std::vector<cell> grid;
//	std::vector<ShortId> OutC;
//	std::vector<ShortId> Out;
//	std::vector<ShortId> In;
//	std::vector<Type>S;
//	std::vector<Vector3> X;
//	std::vector<Vector2> X0;
//
//	std::vector<int> ShiftOut;
//	std::vector<int> ShiftRes;
//	std::vector<int> ShiftX0;
//	std::vector<int> ShiftTry;
//
//	vector<int> sorted_id_cell(count_cells * count_directions); 	// ”пор€доченные индексы €чеек по данному направлению
//
//	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face, inclined_face_inverse,
//		straight_face_inverse);
//
//	std::vector<Type> Illum(4 * count_cells * count_directions, 0);
//	//std::vector<Type> illum(count_cells * count_directions);
//	std::vector<Type> int_scattering(count_cells * count_directions, 0);
//
//	int count = 0;
//	Type norm = 0;
//
//	ofstream ofile;
//	ofile.open(main_dir + "File_with_Logs_solve.txt");
//
//	_clock = -omp_get_wtime();
//	if (ReadCompactFastGridData(count_directions, size_grid, name_file_in_faces, name_file_out_faces, name_file_count_out_faces,
//		name_file_local_x0, name_file_x, name_file_s, name_file_id_neighbors, name_file_centers,
//		name_file_dist_try, name_file_id_try, name_file_res, name_file_sizes, name_file_graph, name_file_shift_out, name_file_shift_res,
//		name_file_shift_x0, name_file_shift_try, grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
//		res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry)) {
//		std::cout << "Error fast reading the data grid vtk\n";
//		return 1;
//	}
//	_clock += omp_get_wtime();
//	std::cout << "\n Reading time of the data_grid file: " << _clock << "\n";
//
//	if (use_cuda)
//	{
//		if (CheckDevice()) use_cuda = false;
//		else
//		{
//			if (InitDevice(count_directions, count_cells, cuda_mod)) use_cuda = false; //CUDA_RETURN("Error cuda init\n")
//			if (HostToDevice(directions, squares, Illum, 1))use_cuda = false;// CUDA_RETURN("Error cuda hostTodevice\n")
//		}
//	}
//
//
//	do {
//		Type _clock = -omp_get_wtime();
//
//		norm = -1;
//		/*---------------------------------- далее FOR по направлени€м----------------------------------*/
//
//		for (register int num_direction = 0; num_direction < count_directions; ++num_direction)
//		{
//			int posX = 0;
//			int posX0 = 0;
//			int posOut = 0;
//			int posIn = 0;
//			int posS = 0;
//
//			int pos_in_res = 0;
//			int id_try_pos = 0;
//
//
//			register int num_cell;
//			ShortId n_out;
//			Vector3 I;
//			Type sumI = 0;
//
//
//			/*---------------------------------- далее FOR по €чейкам----------------------------------*/
//			for (int h = 0; h < count_cells; ++h) {
//				num_cell = sorted_id_cell[num_direction * count_cells + h];
//
//				n_out = OutC[num_direction * count_cells + h];
//
//				sumI = 0;
//
//				for (ShortId i = 0; i < n_out; ++i) {
//					ShortId num_out_face = Out[ShiftOut[num_direction] + posOut++];
//
//					//GetNodes
//					for (int num_node = 0; num_node < 3; ++num_node) {
//						Vector3 x = X[3 * ShiftOut[num_direction] + posX++];
//
//						ShortId num_in_face = In[3 * ShiftOut[num_direction] + posIn++];
//
//						Type I_x0 = CalculateIllumeOnInnerFace(num_cell, num_direction, num_in_face, x, X0, grid, neighbours_id_faces,
//							id_try_pos, pos_in_res, posX0, ShiftRes[num_direction], ShiftX0[num_direction], ShiftTry[num_direction]);
//
//						Type s = S[3 * ShiftOut[num_direction] + posS++];
//
//						I[num_node] = CurGetIllum(num_cell, num_direction, x, s, I_x0, int_scattering);
//						sumI += I[num_node];
//
//					}//num_node
//
//					// если хранить не знгачени€ у коэфф. интерпол€ции
//					{
//						//Vector3 coef;
//						//if (num_out_face == 3)
//						//	coef = inclined_face_inverse * I;// GetInterpolationCoefInverse(inclined_face_inverse, I);
//						//else
//						//	coef = straight_face_inverse * I;// GetInterpolationCoefInverse(straight_face_inverse, I);
//					}
//
//					grid[num_cell].nodes_value[num_out_face] = I; //coef;  // храним значение узлов а не коэф. интерпол€ции
//
//					//Illum[num_direction * (4*count_cells) + num_cell * 4 + num_out_face] = (I[0] + I[1] + I[2]) / 3;					
//
//					int neighbor_id_face = neighbours_id_faces[num_cell * 4 + num_out_face];
//					if (neighbor_id_face >= 0)
//						grid[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] = I;// coef;
//
//
//				} //num_out_face
//
//				for (size_t i = 0; i < 4; i++)
//				{
//					Vector3 Il = grid[num_cell].nodes_value[i];
//					const Type curI = (Il[0] + Il[1] + Il[2]) / 3;
//
//					{
//						Type buf_norm = fabs(Illum[num_direction * (4 * count_cells) + num_cell * 4 + i] - curI) / curI;
//						if (buf_norm > norm) norm = buf_norm;
//					}
//					Illum[num_direction * (4 * count_cells) + num_cell * 4 + i] = curI;
//				}
//				//illum[num_direction * count_cells + num_cell] = sumI / (3 * n_out);				
//			}
//			/*---------------------------------- конец FOR по €чейкам----------------------------------*/
//
//			//std::cout << "End direction number: " << num_direction << '\n';
//		}
//		/*---------------------------------- конец FOR по направлени€м----------------------------------*/
//
//
//
//		if (max_number_of_iter > 1)  // пропуск первой итерации
//		{
//			if (use_cuda)
//			{
//				if (CalculateIntScattering(32, count_cells, count_directions, Illum, int_scattering)) { // если была ошибка пересчитать на CPU 
//					CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
//					ClearDevice(cuda_mod);
//					use_cuda = false;
//				}
//			}
//			else
//				CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
//		}
//
//		//Illum.swap(Illum2);
//
//		_clock += omp_get_wtime();
//
//		//if (count % 2 == 0) norm = NormIllumOmp(Illum, Illum2);  // вычисл€ть норму раз в 2 итерации
//
//		std::cout << "Error:= " << norm << '\n';
//		std::cout << "Time of iter: " << _clock << '\n';
//		std::cout << "End iter_count number: " << count << '\n';
//
//		ofile << "Error:= " << norm << '\n';
//		ofile << "Time of iter: " << _clock << '\n';
//		ofile << "End iter_count number: " << count << '\n';
//
//		count++;
//	} while (norm > accuracy && count < max_number_of_iter);
//
//	//Illum.swap(Illum2);
//	ofile.close();
//
//	//Illum2.clear();	
//
//	std::vector<Type> illum(count_cells * count_directions);
//
//	for (size_t i = 0; i < count_directions; i++)
//		for (size_t j = 0; j < count_cells; j++)
//		{
//			const int N = i * 4 * count_cells + j * 4;
//			Type I = (Illum[N] + Illum[N + 1] + Illum[N + 2] + Illum[N + 3]) / 4;
//			illum[i * count_cells + j] = I;
//		}
//
//
//
//	std::vector<Vector3> Stream(size_grid * 4);
//	MakeStream(Illum, directions, squares, square_surface, Stream);
//	Illum.clear();
//
//
//	vector<Type> energy(size_grid);
//	vector<Vector3> stream(count_cells);
//	vector<Matrix3> impuls(size_grid);
//	std::vector<Type> div_stream(size_grid, 0);
//
//	std::vector<Normals> normals;
//	std::vector<Type> squares_cell;
//	std::vector<Type> volume;
//
//	{
//		ReadNormalFile(name_file_normals, normals);
//		ReadSimpleFileBin(name_file_squares, squares_cell);
//		ReadSimpleFileBin(name_file_volume, volume);
//	
//	/*	std::vector<Vector3> points;
//		ReadSimpleFileBin(main_dir + "face_center.bin", points);
//		Vector3 Stream;
//		for (size_t i = 0; i < size_grid; i++) {
//			for (size_t j = 0; j < 4; j++) {
//				Vector3 x = points[i * 4 + j];
//				Stream = Vector3(sin(x[0]),0,0);
//				div_stream[i] += normals[i].n[j].dot(Stream) * squares_cell[i * 4 + j];
//			}
//			div_stream[i] /= volume[i];
//		}*/
//	}
//
//	for (size_t i = 0; i < size_grid; i++) 
//	{
//		div_stream[i] = 0;
//		for (size_t j = 0; j < 4; j++)
//			div_stream[i] += normals[i].n[j].dot(Stream[i * 4 + j]) * squares_cell[i * 4 + j];
//		div_stream[i] /= volume[i];
//	}
//
//	//debug output:
//	{
//		ofstream outf(main_dir + "divTest.txt");
//
//		for (size_t i = 0; i < size_grid; i++)
//		{
//			for (size_t j = 0; j < 4; j++)
//				if (neighbours_id_faces[4 * i + j] == -2)
//				{
//				
//					outf << "Cell:" << i << '\n';
//					
//					outf << "Neig: ";
//					for (size_t k = 0; k < 4; k++)
//						outf << neighbours_id_faces[4 * i + k] << " ";
//					outf << '\n';
//
//					outf << "Area: \n";
//					for (size_t k = 0; k < 4; k++)
//						outf << setprecision(16) << squares_cell[4 * i + k] << "\n";					
//
//					outf << "Volume:" << setprecision(16) << volume[i] << '\n';
//
//					outf << "Normals: \n";
//					for (size_t k = 0; k < 4; k++)
//						outf << setprecision(16) << normals[i].n[k](0) << " :: " << normals[i].n[k](1) << " ::  " << normals[i].n[k](2) << "\n";
//					
//					outf << "Stream:\n";
//					for (size_t k = 0; k < 4; k++) {
//						outf << setprecision(16) << Stream[i * 4 + k](0) << " :: " << Stream[i * 4 + k](1) << " :: " << Stream[i * 4 + k](2) << " \n";
//					}
//					outf << "divStream= " << setprecision(16) << div_stream[i] << '\n';		
//					outf << "======================================================================\n";
//					break;
//				}
//		}
//
//		outf.close();
//	}
//
//	for (size_t i = 0; i < count_cells; i++)		
//		{
//		const int N = i * 4;
//			Vector3  S = (Stream[N] + Stream[N + 1] + Stream[N + 2] + Stream[N + 3]) / 4;
//			stream[i] = S;
//		}
//
//	for (size_t i = 0; i < size_grid; i++)
//	{
//		div_stream[i] = 0;
//		for (size_t j = 0; j < 4; j++)
//		{
//			int neig = neighbours_id_faces[i * 4 + j];
//			if (neig >= 0)
//				div_stream[i] += normals[i].n[j].dot(stream[neig / 4]) * squares_cell[i * 4 + j];
//			else
//				div_stream[i] += normals[i].n[j].dot(stream[i]) * squares_cell[i * 4 + j];
//		}
//		div_stream[i] /= volume[i];
//	}
//
//	use_cuda = false;
//	if (use_cuda && cuda_mod == 1) {
//		CalculateEnergy(32, count_cells, count_directions, energy);
//		CalculateStream(32, count_cells, count_directions, stream);
//		CalculateImpuls(32, count_cells, count_directions, impuls);
//
//	}
//	else {
//		MakeEnergy(illum, squares, square_surface, energy);
//		//MakeStream(illum, directions, squares, square_surface, stream);
//		MakeImpuls(illum, directions, squares, square_surface, impuls);
//
//		//MakeDivStream(stream, normals, squares_cell, volume, neighbours_id_faces, div_stream);
//	}
//	if (use_cuda) {
//		ClearDevice(cuda_mod);
//	}
//
//#ifdef USE_VTK
//	WriteFileSolution(out_file_grid_vtk, illum, energy, stream, impuls, name_file_vtk);
//#else
//	WriteFileSolution(main_dir, illum, energy, stream, impuls, div_stream);
//#endif
//	//WriteFileSolution(out_file_grid_vtk, Illum, energy, name_file_vtk);
//
//	std::cout << "End programm\n";
//	return 0;
//}
//
//
//
