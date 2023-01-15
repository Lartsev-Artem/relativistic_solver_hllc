#include "solve_short_characteristics_main.h"
#include "solve_short_characteristics_hllc.h"

#ifdef WRITE_LOG	
#define WRITE(ofile, str)  \
ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app); \
ofile << str; \
ofile.close();
#else
#define WRITE(ofile, str) {}
#endif

#define WRITE_FILE(name_file, data, value) \
{ \
FILE* f;\
int n = cells.size(); \
f = fopen(name_file, "wb"); \
fwrite(&n, sizeof(int), 1, f); \
for (auto& el : data) \
{	\
	fwrite(&el.value, sizeof(el.value), 1, f);	\
}	\
fclose(f);\
}


static int RebuildConvToPhysHLLC(std::vector<elem>&cells)
{	
	for (auto& el : cells)
	{
		const Type d = el.val.d;
		el.phys_val.d = d;
		el.phys_val.v = el.val.v / d;		
		const Type vv = el.phys_val.v.dot(el.phys_val.v);
		el.phys_val.p = (el.val.p - vv * d / 2.) * (gamma1 - 1);	
	}
	return 0;
}
size_t WriteFileSolution(const std::string& main_dir, const std::vector<Type>& vector_illum, const std::vector<elem>& cells) {
	
	int n;
	FILE* f;
#ifndef ONLY_HLLC

	f = fopen((main_dir + "Illum.bin").c_str(), "wb");
	n = vector_illum.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_illum.data(), sizeof(Type), n, f);
	fclose(f);

	WRITE_FILE((main_dir + "energy.bin").c_str(), cells, illum_val.energy);

	WRITE_FILE((main_dir + "stream.bin").c_str(), cells, illum_val.stream);

	WRITE_FILE((main_dir + "impuls.bin").c_str(), cells, illum_val.impuls);

	WRITE_FILE((main_dir + "divstream.bin").c_str(), cells, illum_val.div_stream);

	WRITE_FILE((main_dir + "divimpuls.bin").c_str(), cells, illum_val.div_impuls);
#endif

#ifndef ONLY_ILLUM

	WRITE_FILE((main_dir + "density.bin").c_str(), cells, phys_val.d);

	WRITE_FILE((main_dir + "pressure.bin").c_str(), cells, phys_val.v);

	WRITE_FILE((main_dir + "velocity.bin").c_str(), cells, phys_val.p);
#endif

	return 0;
}

int WriteStepTimeSolve(const std::string& main_dir, const int num_dir_to_illum, const std::vector<elem>& cells,
	const std::vector<Type>& Illum) {

	std::vector<Type> illum;

#ifndef ONLY_HLLC
	const int size_grid = cells.size();
	illum.resize(size_grid, 0);
	const int i = num_dir_to_illum;
	for (size_t j = 0; j < size_grid; j++)
	{
		const int N = i * 4 * size_grid + j * 4;
		Type I = (Illum[N] + Illum[N + 1] + Illum[N + 2] + Illum[N + 3]) / 4;
		illum[j] = I;
	}
#endif

	WriteFileSolution(main_dir, illum, cells);

	return 0;
}


static Type FormTimeStepToHLLC(const Type h, const Type k, const std::vector<elem>& cells) {
	
	Type c_max = -1;
	for(auto& el : cells)	
	{		
		const Type d = el.val.d;		
		const Type vv = el.val.v.dot(el.val.v) / (d * d);
		const Type p = (el.val.p - vv * d / 2.) * (gamma1 - 1);
		const Type a = sqrt(gamma1 * p / d);
		const Type c = sqrt(vv) + a;
		if (c > c_max) c_max = c;
	}
	return k * h / c_max;
}

static int ReBuildHllcPhysToConv(std::vector<elem>& cells)
{
	for(auto& el : cells)	
	{
		const Type v = el.val.v.dot(el.val.v);
		const Type d = el.val.d;

		//el.val.d = d;
		el.val.v *= d;
		el.val.p = el.val.p / (gamma1 - 1) + d * v / 2;		
	}

	return 0;
}
static int ReBuildRHllcPhysToConv(std::vector<elem>& cells)
{
	for (auto& el : cells)
	{
		const Type v = el.phys_val.v.dot(el.phys_val.v);
		const Type d = el.phys_val.d;
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * el.phys_val.p / d;
		const Type dhGG = d * h * Gamma * Gamma;

		el.val.d = Gamma * d;
		el.val.v *= dhGG;
		el.val.p = dhGG - el.phys_val.p;
	}
	return 0;
}
int InitHLLC(std::vector<elem>& cells)
{	
	size_grid = cells.size();
	const int n = size_grid;
	
#if defined Cylinder 
	for (size_t i = 0; i < size_grid; i++)
	{
		Vector3 x = centers[i];
		if (Vector2(x[1], x[2]).norm() < 0.2 && x[0] < 0.5)
		{
			density[i] = 0.1;
			pressure[i] = 0.01;
			velocity[i] = Vector3(0.99, 0, 0);
		}
		else
		{
			density[i] = 10;
			pressure[i] = 0.01;
			velocity[i] = Vector3(0, 0, 0);
		}
	}
#endif
#if defined RHLLC && defined Cube
	for (size_t i = 0; i < size_grid; i++)
	{
		Type x = centers[i][0];
		if (x < 0.5)
		{
			density[i] = 1;
			pressure[i] = 1;
			velocity[i] = Vector3(0.9, 0, 0);
		}
		else
		{
			density[i] = 1;
			pressure[i] = 10;
			velocity[i] = Vector3(0, 0, 0);
		}
		/*		Type x = centers[i][0];
				if (x < 0.5)
				{
					density[i] = 1;
					pressure[i] = 10;
					velocity[i] = Vector3(-0.6, 0, 0);
				}
				else
				{
					density[i] = 10;
					pressure[i] = 20;
					velocity[i] = Vector3(0.5, 0, 0);
				}*/
				/*if (centers[i][0] < 0.5) {
					density[i] = 1;
					pressure[i] = 1;
					velocity[i] = Vector3(0, 0, 0);
				}
				else {
					density[i] = 0.125;
					pressure[i] = 0.1;
					velocity[i] = Vector3(0, 0, 0);
				}*/
	}

#endif

#if defined Jet
	for (size_t i = 0; i < size_grid; i++)
	{
		const Type betta = 0.1;
		const Type a = 1;
		const Type b = 0.001;
		Type x = centers[i][0];
		density[i] = a * exp(-x * x / betta) + b;
		pressure[i] = a * exp(-x * x / betta) + (1e-5);
		velocity[i] = Vector3(1e-4, 0, 0);
	}
#endif
#if (defined Cone || defined Sphere) && (!defined Jet)
	for (size_t i = 0; i < size_grid; i++)
	{
		{
			// a~0.4
			density[i] = 0.1;
			pressure[i] = 0.01;
			velocity[i] = Vector3(0, 0, 0);
		}
	}
#endif
#if defined Step
	for (size_t i = 0; i < size_grid; i++)
	{
		density[i] = 0.1;
		pressure[i] = 0.01;
		velocity[i] = Vector3(0, 0, 0);
	}
#endif
#if Plane
	for (size_t i = 0; i < size_grid; i++)
	{
		if (centers[i][0] > 0 && centers[i][2] > 0)
		{
			density[i] = 0.1;
			pressure[i] = 0.01;
			velocity[i] = Vector3(0, 0, 0);
		}
		else if (centers[i][0] < 0 && centers[i][2] > 0)
		{
			density[i] = 0.1;
			pressure[i] = 1;
			velocity[i] = Vector3(0.99, 0, 0);
		}
		else if (centers[i][0] < 0 && centers[i][2] < 0)
		{
			density[i] = 0.5;
			pressure[i] = 1;
			velocity[i] = Vector3(0, 0, 0);
		}
		else if (centers[i][0] > 0 && centers[i][2] < 0)
		{
			density[i] = 0.1;
			pressure[i] = 1;
			velocity[i] = Vector3(0, 0, 0.99);
		}
	}
#endif
#if defined Cube  && !defined RHLLC//SODA

	for(auto& el : cells)	
	{
		if (el.center[0] < 0.5) 
		{
			el.val.d = 1;
			el.val.p = 1;
			el.val.v = Vector3(0, 0, 0);
		}
		else {
			el.val.d = 0.125;
			el.val.p = 0.1;
			el.val.v = Vector3(0, 0, 0);
		}		
	}
#endif	

#if !defined RHLLC
	ReBuildHllcPhysToConv(cells);
#else
	ReBuildRHllcPhysToConv(cells);
#endif

	return 0;
}

int SetHLLCValue(Type& tau, Type& CFL, Type& print_timer, Type& h, Type& T) {
	Type cur_timer = 10;  // больше print_timer для вывода первого шага
#ifdef ONLY_HLLC
	 tau = 1e-5;
	 CFL = 0.1;
	 print_timer = 0.01;
#else
	 tau = 1e-8;
	 CFL = 0.001;
	 print_timer = 1e-5;
#endif

	Type t = 0.0;
	 T = 0.4;
#if defined Sphere
	 h = 0.0128079221422811;
#elif defined Cone && !defined Jet
	//const Type h = 0.0032396496530313; // test grid
	 h = 0.0022530803678209; //grid: 112023 cells
#elif defined Jet
	 h = 0.0012548169651948;
#elif defined Cube
	//const Type h = 0.0007123669658939; // Soda1d_2
	//const Type h  = 0.0010828369115320; // Soda1d
	 h = 0.0010307259619874; // Soda1d_3
#elif defined Step
	 h = 0.0018751819368151;
	T = 5;
	CFL = 0.5;
	print_timer = 0.1;
#elif defined Cylinder
	 h = 0.0064085233935784;
	T = 10;
	CFL = 0.05;
	print_timer = 0.0000001;
#else
	 h = 1;
	printf("Bad geometry define\n");
#endif	

}

/*
	ReadData (vtk or not)

	Rebuild/ rewrite section

	BuildHLLC struct

	Solve HLLC

	Solve illum -> divW, U, divT

	Solve HLLC + Illum

	next time step

*/

int main(int argc, char* argv[])
{

#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	RHLLC_MPI(main_dir, centers, neighbours_id_faces, normals, squares_cell, volume);
	MPI_RETURN(0);
#endif


	std::string name_file_settings = "";
	int max_number_of_iter = 1;
	Type accuracy = 1e-5;
	bool use_cuda = true;
	const int cuda_mod = 1; //1 - все массивы, 0 - min

	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_solve.txt";
	else
		name_file_settings = argv[1];
	if (argc > 2)
		max_number_of_iter = std::stoi(argv[2]);

	std::cout << "Max_number_of_iter= " << max_number_of_iter << '\n';

	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string out_file_grid_vtk;
	std::string solve_direction;

	std:string main_dir;

	ofstream ofile;
#ifdef WRITE_LOG	
	std::remove((main_dir + "File_with_Logs_solve.txt").c_str());	
#endif

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, main_dir, solve_direction))
	{
		ERR_RETURN("Error reading the start settings\n")
	}

	WRITE(ofile, "start settings\n");

#if defined HLLC_1D
	HLLC_1d(main_dir);
	printf("HLLC 1d end\n");
	return 0;

#elif defined RHLLC_1D
	RHLLC_1d(main_dir);
	printf("RHLLC 1d end\n");
	return 0;
#endif

	std::string name_file_graph = main_dir + "graph";
	std::string name_file_solve = main_dir + "Solve";

	//--------------------------Файлы расчётных данных-----------------------------------//
	std::string name_file_in_faces = main_dir + "InId";
	std::string name_file_out_faces = main_dir + "OutId";
	std::string name_file_count_out_faces = main_dir + "CountOutId";
	std::string name_file_local_x0 = main_dir + "LocX0";
	std::string name_file_x = main_dir + "X";
	std::string name_file_s = main_dir + "S";

	//--------------------------Файлы трассировки сквозь внутреннюю границу----------------------------------------//
	std::string name_file_dist_try = main_dir + "dist_defining_faces";
	std::string name_file_id_try = main_dir + "id_defining_faces";
	std::string name_file_res = main_dir + "ResBound";
	std::string name_file_sizes = main_dir + "Size";

	//--------------------------Файлы сдвигов по направлениям в расчётных файлах-----------------------------------//
	std::string name_file_shift_out = main_dir + "ShiftOut";
	std::string name_file_shift_res = main_dir + "ShiftRes";
	std::string name_file_shift_x0 = main_dir + "ShiftX0";
	std::string name_file_shift_try = main_dir + "ShiftTry";

	//--------------------------Файлы геометрии--------------------------------------------------------------------//
	std::string name_file_id_neighbors = main_dir + "pairs";
	std::string name_file_normals = main_dir + "normals.bin";
	std::string name_file_centers = main_dir + "centers.bin";
	std::string name_file_squares = main_dir + "squares.bin";
	std::string name_file_volume = main_dir + "volume.bin";
	
#ifdef USE_VTK
#else
	double _clock = -omp_get_wtime();
	if (ReadDataArray(class_file_vtk, main_dir, density, absorp_coef, rad_en_loose_rate, velocity, pressure, size_grid, true))
	{
		ERR_RETURN("Error reading the data array vtk\n")
	}

	WRITE(ofile, "reading the data array\n");

	
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";

	WRITE(ofile, "Init_hllc\n");
	
	//------------------------------ Illume section-----------------------------
	std::vector<Vector3> directions;
	std::vector<Type> squares;

	_clock = -omp_get_wtime();
	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface))
	{
		ERR_RETURN("Error reading the sphere direction\n")
	}
		_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

	WRITE(ofile, "reading sphere_direction\n");

	const int count_directions = directions.size();
	const int count_cells = size_grid;
	printf("\nReal size task = %d\n", size_grid);

	//std::vector<int> neighbours_id_faces;
	std::vector<Vector3> nodes_value;
	std::vector<cell> grid;
	std::vector<ShortId> OutC;
	std::vector<ShortId> Out;
	std::vector<ShortId> In;
	std::vector<Type>S;
	std::vector<Vector3> X;
	std::vector<Vector2> X0;

	std::vector<uint64_t> ShiftOut;
	std::vector<uint64_t> ShiftRes;
	std::vector<uint64_t> ShiftX0;
	std::vector<int> ShiftTry;

	//std::vector<Normals> normals;
	//std::vector<Type> squares_cell;
	//std::vector<Type> volume;
	//std::vector<Vector3> centers;

	std::vector<face> faces;
	std::vector<elem> cells;

	if (ReadSimpleFileBin(name_file_centers, faces)) ERR_RETURN("Error reading file centers\n");
	if (ReadSimpleFileBin(name_file_squares, cells)) ERR_RETURN("Error reading file squares_cell\n");

	std::vector<int> inner_bound;
	std::vector<int> sorted_id_cell(count_cells * count_directions); 	// Упорядоченные индексы ячеек по данному направлению
	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face, inclined_face_inverse,
		straight_face_inverse);

	/*if (ReadNormalFile(name_file_normals, normals)) ERR_RETURN("Error reading file normals\n");
	if (ReadCentersOfTetra(name_file_centers, centers)) ERR_RETURN("Error reading file centers\n");
	if (ReadSimpleFileBin(name_file_squares, squares_cell)) ERR_RETURN("Error reading file squares_cell\n");
	if (ReadSimpleFileBin(name_file_volume, volume)) ERR_RETURN("Error reading file volume\n");

	if (ReadSimpleFileBin(name_file_id_neighbors + ".bin", neighbours_id_faces)) ERR_RETURN("Error reading file neighbours\n");*/
	//if (ReadSimpleFile("D:\\Desktop\\FilesCourse\\graphSet\\inner_bound.txt", inner_bound)) ERR_RETURN("Error reading file inner_bound\n");

	WRITE(ofile, "reading simple files\n");



#ifndef ONLY_HLLC  
	_clock = -omp_get_wtime();
	if (ReadCompactFastGridData(count_directions, size_grid, main_dir + "File_with_Logs_solve.txt", name_file_in_faces, name_file_out_faces, name_file_count_out_faces,
		name_file_local_x0, name_file_x, name_file_s, name_file_id_neighbors, name_file_centers,
		name_file_dist_try, name_file_id_try, name_file_res, name_file_sizes, name_file_graph, name_file_shift_out, name_file_shift_res,
		name_file_shift_x0, name_file_shift_try, grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
		res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry))
		ERR_RETURN("Error fast reading the data grid vtk\n")

		_clock += omp_get_wtime();
	std::cout << "\n Reading time of the data_grid file: " << _clock << "\n";
#endif

	WRITE(ofile, "reading  data_grid file\n");


	std::vector<Type> Illum; // отдельно от сетки из-за расчета на cuda
	std::vector<Type> int_scattering; // отдельно от сетки из-за расчета на cuda
//
//	std::vector<Type> energy; // (size_grid, 0);
//	std::vector<Vector3> stream; //(size_grid, Vector3(0,0,0));
//
//	//vector<Type> prev_energy(size_grid, 0);
//	//vector<Vector3> prev_stream(size_grid, Vector3(0, 0, 0));
//
//	std::vector<Matrix3> impuls; //(size_grid, Matrix3::Zero());
//	std::vector<Type> div_stream; //(size_grid, 0);
//	std::vector<Vector3> div_impuls; //(size_grid, Vector3(0,0,0));
#ifndef ONLY_HLLC  
	Illum.resize(4 * count_cells * count_directions, 0);
	int_scattering.resize(count_cells * count_directions, 0);
//
//	energy.resize(size_grid, 0);
//	stream.resize(size_grid, Vector3(0, 0, 0));
//
//	impuls.resize(size_grid, Matrix3::Zero());
//	div_stream.resize(size_grid, 0);
//	div_impuls.resize(size_grid, Vector3(0, 0, 0));
#endif

	if (use_cuda)
	{
		if (CheckDevice()) use_cuda = false;
		else
		{
			if (InitDevice(count_directions, count_cells, cuda_mod)) use_cuda = false; //CUDA_RETURN("Error cuda init\n")
			if (HostToDevice(directions, squares, Illum, 1))use_cuda = false;// CUDA_RETURN("Error cuda hostTodevice\n")
		}
	}
	//------------------------------------------------------------------------------------------------------------


	//---------------------------------------HLLC section--------------------------------------------------------------
	Type cur_timer = 10;  // больше print_timer для вывода первого шага
	Type tau, CFL, print_timer, h, T;
	Type t = 0.0;

#ifndef ONLY_ILLUM
	//---------------Начальные данные для HLLC---------------------------------//
	InitHLLC(cells);
	SetHLLCValue(tau, CFL, print_timer, h, T);
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
	
//------------------------------------------------------------------------------------------------------------

#if 0 //ONLY_HLLC  // инициализация составляющий для излучения (не надо для квазилинейного случая)
		CalculateIllumOptMemory(main_dir, use_cuda, cuda_mod, accuracy, max_number_of_iter, directions, Illum, int_scattering, squares,
			grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
			res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry, U_full_prev);

	CalculateIllumParamOptMemory(size_grid, count_directions,
		directions, squares, normals, squares_cell, volume,
		Illum, energy, stream, impuls, div_stream, div_impuls);

	printf("Init Illum completed\n\n");
#endif

	//------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------
	int res_count = 0; // счётчик решений	
	Type full_time = -omp_get_wtime();
	while (t < T)
	{	
		if (cur_timer >= print_timer)
		{
			std::string file_solve = solve_direction + "Solve" + std::to_string(res_count);
			RebuildConvToPhysHLLC(cells);
			WriteStepTimeSolve(file_solve, 0, cells, Illum);
			cur_timer = 0;
			res_count++;

			printf("\n t= %f,  tau= %lf,  res_step= %d\n", t, tau, res_count);
		}
#ifdef RHLLC
		HLLC_Rel(tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);
#elif defined NEW_CLASS
		HLLC(tau, faces, cells);
#else
		HLLC(size_grid, tau, bound_cells, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);
#endif
#endif

#ifndef ONLY_HLLC

		CalculateIllum(main_dir, use_cuda, cuda_mod, accuracy, max_number_of_iter, directions, Illum, int_scattering, squares,
			grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
			res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry, U_full);

		CalculateIllumParam(size_grid, count_directions,
			directions, squares, normals, squares_cell, volume,
			Illum, energy, stream, impuls, div_stream, div_impuls);

		SolveIllumAndHLLC(size_grid, tau, U_full, energy, energy, stream, stream, div_stream, div_impuls);

		//energy.swap(prev_energy);
		//stream.swap(prev_stream);
#endif
		t += tau;
		cur_timer += tau;

		// ---------------формирование шага по времени--------------------//
#ifdef RHLLC
		tau = FormTimeStepToRHLLC(size_grid, h, CFL);
#else
		tau = FormTimeStepToHLLC(h, CFL, cells);
#endif
		if (tau < 0)
		{
			printf("Error tau = %lf\n", tau);
			break;
		}

		WRITE(ofile, "\nt= " << t << "; tau= " << tau << "; step= " << res_count << '\n');		

	}// while(t < T)

	full_time += omp_get_wtime();
	printf("Time while(): %f\n", full_time);

	if (use_cuda)
	{
		ClearDevice(cuda_mod);
	}

	RebuildConvToPhysHLLC(cells);
	WriteStepTimeSolve(solve_direction + "Solve" + std::to_string(res_count), 0, cells, Illum);
#endif // USE_VTK

	std::cout << "End programm\n";

	return 0;
}


#define READ_FILE(name_file, arr)\
{\
	FILE* f;\
	f = fopen(name_file, "rb");\
	if (!f) { printf("Error reading file: %s\n", name_file); return 1; }	\
	fread_unlocked(arr.data(), sizeof(arr[0]), arr.size(), f);\
	fclose(f);\
}

size_t ReadCompactFastGridData(const int count_dir, const int N, const Str_Type& main_dir,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	Str_Type& name_file_shift_out, Str_Type& name_file_shift_res, Str_Type& name_file_shift_x0, Str_Type& name_file_shift_try,
	std::vector<cell>& grid, std::vector<int>& neighbours_id_face, std::vector<ShortId>& OutC,
	std::vector<ShortId>& Out, std::vector<ShortId>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<int>& sorted_id_cell,
	vector<uint64_t>& ShiftOut, vector<uint64_t>& ShiftRes, vector<uint64_t>& ShiftX0, vector<int>& ShiftTry,
	std::vector<face>& faces, std::vector<elem>& cells) 
{
	ofstream ofile;
	FILE* f;
	  
	ShiftOut.resize(count_dir);
	READ_FILE((name_file_shift_out + ".bin").c_str(), ShiftOut);
	WRITE(ofile, "reading  ShiftOut file\n");

	ShiftRes.resize(count_dir);
	READ_FILE((name_file_shift_res + ".bin").c_str(), ShiftRes);
	WRITE(ofile, "reading  ShiftRes file\n");

	ShiftX0.resize(count_dir);
	READ_FILE((name_file_shift_x0 + ".bin").c_str(), ShiftX0);
	WRITE(ofile, "reading  ShiftX0 file\n");

#ifdef TRY_SHIFT
	ShiftTry.resize(count_dir);
	READ_FILE((name_file_shift_try + ".bin").c_str(), ShiftTry);
	WRITE(ofile, "reading  ShiftTry file\n");
#endif
	
#if ONE_GRAPH
	READ_FILE((name_file_graph + ".bin").c_str(), sorted_id_cell);	
#else
	for (size_t i = 0; i < count_dir; i++)
	{		
		f = fopen((name_file_graph + std::to_string(i) + ".bin").c_str(), "rb");
		fread_unlocked(sorted_id_cell.data() + i * N, sizeof(int), N, f);
		fclose(f);
	}
#endif

	WRITE(ofile, "reading  sorted_id_cell file\n");

	grid.resize(N);
	int countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry;
	ReadSizes(name_file_sizes, countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry);

	WRITE(ofile, "reading  size file\n");


	res_inner_bound.resize(countRes);
	READ_FILE((name_file_res + ".bin").c_str(), res_inner_bound);
	WRITE(ofile, "reading  res_inner_bound file\n");


#ifdef TRY_SHIFT
	id_try_surface.resize(countTry);
	READ_FILE((name_file_id_try + ".bin").c_str(), id_try_surface);
	WRITE(ofile, "reading  id_try_surface file\n");
#endif
	
	S.resize(countS);
	READ_FILE((name_file_s + ".bin").c_str(), S);
	WRITE(ofile, "reading  S file\n");

	
	OutC.resize(countOutC);
	READ_FILE((name_file_count_out_faces + ".bin").c_str(), OutC);
	WRITE(ofile, "reading  OutC file\n");

	Out.resize(countOut);
	READ_FILE((name_file_out_faces + ".bin").c_str(), Out);
	WRITE(ofile, "reading  Out file\n");


	In.resize(countIn);
	READ_FILE((name_file_in_faces + ".bin").c_str(), In);
	WRITE(ofile, "reading  In file\n");

	X.resize(countX);
	READ_FILE((name_file_x + ".bin").c_str(), X);
	WRITE(ofile, "reading  X file\n");
	
	X0.resize(countX0);
	READ_FILE((name_file_local_x0 + ".bin").c_str(), X0);
	WRITE(ofile, "reading  X0 file\n");

	return 0;
}