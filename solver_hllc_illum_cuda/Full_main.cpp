#include "solve_short_characteristics_main.h"
#include "solve_short_characteristics_hllc.h"


#ifdef SAVE_DUMP_HLLC
uint32_t size_dump = 0;
uint32_t w_ptr = 0;
std::vector<double> dump;
int start_section = -666;
bool bad_hllc_flag = false;

void WriteDump(std::string name_file_dump, std::vector<double>& buffer)
{

	FILE* f;
	f = fopen(name_file_dump.c_str(), "wb");	
	fwrite(&size_dump, sizeof(int), 1, f);
	fwrite(buffer.data(), sizeof(double), size_dump, f);
	fclose(f);
}
#endif

int ReBuild(std::string name_file_sphere_direction, std::string main_dir);

Type FormTimeStepToHLLC(const int n, const Type h, const Type k, const std::vector<VectorX>& U_array) {

	// формирование шага по времени	
	Eigen::VectorXd U(5);

	Type c_max = -1;
	for (size_t i = 0; i < n; i++)
	{
		U = U_array[i];
		const Type d = U[0];
		const Vector3 vel(U(1) / d, U(2) / d, U(3) / d);
		const Type v = vel.norm();
		const Type p = (U(4) - v * v * d / 2.) * (gamma1 - 1);
		const Type a = sqrt(gamma1 * p / d);
		const Type c = v + a;
		if (c > c_max) c_max = c;
	}
	return k * h / c_max;
}
int BoundDataToHLLC(const std::vector<Normals>& normals, const std::vector<int>& neighbours_id_faces,
	const std::vector<Vector3> centers, std::vector<VectorX>& U_array)
 {
#if defined RHLLC
#if !defined Cylinder // jet 2d

	VectorX buf(5);
	for (size_t num_cell = 0; num_cell < size_grid; num_cell++)
	{
		for (size_t i = 0; i < 4; i++) // по граням
		{
			const int neig = neighbours_id_faces[4 * num_cell + i];

			if (neig == eBound_OutSource)
			{
				const Vector3 velocity(0.99, 0, 0);
				const Type pressure = 0.01;
				const Type d = 0.1;

				const Type v = velocity.dot(velocity);
				const Type Gamma = 1. / sqrt(1 - v);
				const Type h = 1 + gamma_g * pressure / d;
				const Type dhGG = d * h * Gamma * Gamma;

				buf[0] = Gamma * d;

				buf[1] = dhGG * velocity[0];
				buf[2] = dhGG * velocity[1];
				buf[3] = dhGG * velocity[2];

				buf[4] = dhGG - pressure;

				U_array[num_cell] = buf;
				break;
			}

		}
	}

	return 0;
#endif
	return 0;
#endif

#if defined Sphere

	VectorX buf(5);
	for (int i = 0; i < size_grid; i++)
	{
		Type d = U_array[i][0];
		Vector3 vel(U_array[i][1] / d, U_array[i][2] / d, U_array[i][3] / d);
		Type v = vel.norm();
		Type p = (U_array[i](4) - v * v * d / 2.) * (gamma1 - 1);
		//if (p < 1e-5) U_array[i][4] = (1e-5) / (gamma1 - 1) + d * v * v / 2;

		for (int j = 0; j < 4; j++)
			if (neighbours_id_faces[i * 4 + j] == -2)
			{				
				{
					// звук - sqrt(gamma*p/rho)
					//pressure[i] = (vector_U[i](4) - v * v * d / 2.)* (gamma1 - 1);

					d = 1;
					v = 1.5; //~a-0.4 
					Vector3 VV(d * v * normals[i].n[j][0], d * v * normals[i].n[j][1], d * v * normals[i].n[j][2]);
					Type  V = VV.norm();
					p = 0.1;
					double e = p / (gamma1 - 1) + d * V * V / 2;
					//buf << d, d* v, 0, 0, e;// 6.286;
					buf << d, -d*VV[0], -d* VV[1], -d* VV[2], e;
					U_array[i] = buf;
					break;
				}
			}
	}
#endif

#if defined Jet
	return 0;
#endif

#if defined Cone   //// задание константы на левой границе для конуса
	VectorX buf(5);

	for (int i = 0; i < size_grid; i++)
	{
		Type d = U_array[i][0];
		Vector3 vel(U_array[i][1] / d, U_array[i][2] / d, U_array[i][3] / d);
		Type v = vel.norm();
		Type p = (U_array[i](4) - v * v * d / 2.)* (gamma1 - 1);
		//if (p < 1e-5) U_array[i][4] = (1e-5) / (gamma1 - 1) + d * v * v / 2;

		for (int j = 0; j < 4; j++)
			if ((normals[i].n[j] - Vector3(-1, 0, 0)).norm() < eps && neighbours_id_faces[i * 4 + j] < 0)
			{
				/*входной поток не по всей площади*/ //if (sqrt(centers[i][2] * centers[i][2] + centers[i][1] * centers[i][1]) < 0.1)
				{
					// звук - sqrt(gamma*p/rho)
					//pressure[i] = (vector_U[i](4) - v * v * d / 2.)* (gamma1 - 1);
					d = 1;  
					v = 1.5; //~a-0.4 
					p = 0.1;
					double e = p / (gamma1 - 1) + d * v * v / 2;
					buf << d, d*v , 0, 0, e;// 6.286;
					U_array[i] = buf;
					break;
				}
			}
	}
#endif


#if 0 //// задание константы на левой границе 2d
	VectorX buf(5);

	for (auto& el : left)
	{
		buf << 1, 0, 3, 0, 6.286;
		U_full_prev[el] = buf;
	}
#endif
#if defined Step
	VectorX buf(5);
	for (size_t num_cell = 0; num_cell < size_grid; num_cell++)
	{
		for (size_t i = 0; i < 4; i++) // по граням
		{
			const int neig = neighbours_id_faces[4 * num_cell + i];
			{
				if (neig ==  eBound_FreeBound)
					if ((normals[num_cell].n[i] - Vector3(-1, 0, 0)).norm() < 1e-5)
					{
						// звук - sqrt(gamma*p/rho)
						//pressure[i] = (vector_U[i](4) - v * v * d / 2.)* (gamma1 - 1);
						Type  d = 1;
						Type  v = 1.2; //~a-0.36 -> 3M
						Type  p = 0.1;
						double e = p / (gamma1 - 1) + d * v * v / 2;
						buf << d, d* v, 0, 0, e;// 6.286;
						U_array[num_cell] = buf;
						break;
					}
			}
		}	
	}
#endif

#if defined Const_1d	// задание константы на левой границе 1d
	VectorX buf(5);
	for (size_t num_cell = 0; num_cell < size_grid; num_cell++)
	{
		for (size_t i = 0; i < 4; i++) // по граням
		{
			const int neig = neighbours_id_faces[4 * num_cell + i];
			{
				if (neig < 0)
					if ((normals[num_cell].n[i] - Vector3(-1, 0, 0)).norm() < 0.001)
					{
						//U_R << d, d* velocity[neig][0], d* velocity[neig][1], d* velocity[neig][2], e_substance[neig];
						buf << 1, 0, 0, 0, 1 / (gamma1 - 1);
						U_array[num_cell] = buf;
						break;
					}
			}
		}
	}
#endif
	return 0;
}


int WriteStepTimeSolve(const std::string& main_dir, const int num_dir_to_illum,
	const std::vector<Type>& Illum, const std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	const std::vector<Type>& div_stream, const vector<Vector3>& div_impuls, const std::vector<VectorX>& U_full);
int SolveIllumAndHLLC(const int N, const Type tau, std::vector<VectorX>& U_full,
	const vector<Type>& prev_energy, const vector<Type>& energy,
	const vector<Vector3>& prev_stream, const vector<Vector3>& stream,
	const std::vector<Type>& div_stream, const  vector<Vector3>& div_impuls);
int CalculateIllumParam(const int count_cells, const int count_directions,
	const std::vector<Vector3>& directions, const std::vector<Type>& squares,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Type>& Illum, std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	std::vector<Type>& div_stream, vector<Vector3>& div_impuls);
int CalculateIllumParamOptMemory(const int count_cells, const int count_directions,
	const std::vector<Vector3>& directions, const std::vector<Type>& squares,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Type>& Illum, std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	std::vector<Type>& div_stream, vector<Vector3>& div_impuls);

int MakeIllum(const int count_cells, const int count_directions, const std::vector<Type>& Illum, std::vector<Type>& illum);
int MakeDivStream(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume, vector<Vector3>& stream,
	vector<Type>& div_stream);
int MakeDivImpuls(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume, vector<Matrix3>& impuls,
	vector<Vector3>& div_impuls);

#ifdef USE_VTK
size_t ReBuildDataArrayFull(const size_t class_file_vtk, const std::string name_file_vtk, const std::string name_file_out,
	const std::string& main_dir, const bool is_print=false) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	if (ReadFileVtk(0, name_file_vtk, unstructured_grid, foo, foo, foo, is_print)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> EnergyArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> DivStreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> DensityArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> PressureArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> StreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	StreamArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> DivImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	DivImpulsArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> VelocityArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	VelocityArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> ImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	ImpulsArray->SetNumberOfComponents(9);

	std::vector<Type> vector_energy;
	std::vector<Type> vector_illum;
	std::vector<Type> vector_pressure;
	std::vector<Type> vector_density;
	std::vector<Type> vector_div_stream;
	std::vector<Vector3> vector_stream;
	std::vector<Vector3> vector_velocity;
	std::vector<Vector3> vector_div_impuls;
	std::vector<Matrix3> vector_impuls;
	
	int size;

	std::unique_ptr<FILE, int(*)(FILE*)> file_illum(fopen((main_dir + "Illum.bin").c_str(), "rb"), fclose);
	if (!file_illum) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_illum.get());
	//if (n != size) { printf("Error size data\n"); return 1; }
	vector_illum.resize(size);
	fread_unlocked(vector_illum.data(), sizeof(Type), size, file_illum.get());
	fclose(file_illum.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_divstream(fopen((main_dir + "divstream.bin").c_str(), "rb"), fclose);
	if (!file_divstream) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_divstream.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_div_stream.resize(size);
	fread_unlocked(vector_div_stream.data(), sizeof(Type), size, file_divstream.get());
	fclose(file_divstream.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_energy(fopen((main_dir + "energy.bin").c_str(), "rb"), fclose);
	if (!file_energy) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_energy.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_energy.resize(size);
	fread_unlocked(vector_energy.data(), sizeof(Type), size, file_energy.get());
	fclose(file_energy.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_density(fopen((main_dir + "density.bin").c_str(), "rb"), fclose);
	if (!file_density) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_density.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_density.resize(size);
	fread_unlocked(vector_density.data(), sizeof(Type), size, file_density.get());
	fclose(file_density.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_pressure(fopen((main_dir + "pressure.bin").c_str(), "rb"), fclose);
	if (!file_pressure) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_pressure.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_pressure.resize(size);
	fread_unlocked(vector_pressure.data(), sizeof(Type), size, file_pressure.get());
	fclose(file_pressure.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_stream(fopen((main_dir + "stream.bin").c_str(), "rb"), fclose);
	if (!file_stream) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_stream.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_stream.resize(size);
	fread_unlocked(vector_stream.data(), sizeof(Vector3), size, file_stream.get());
	fclose(file_stream.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_velocity(fopen((main_dir + "velocity.bin").c_str(), "rb"), fclose);
	if (!file_velocity) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_velocity.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_velocity.resize(size);
	fread_unlocked(vector_velocity.data(), sizeof(Vector3), size, file_velocity.get());
	fclose(file_velocity.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_impuls(fopen((main_dir + "impuls.bin").c_str(), "rb"), fclose);
	if (!file_impuls) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_impuls.get());
	if (n != size) { printf("Error size data (%d - %d)\n", n, size); return 1; }
	vector_impuls.resize(size);
	fread_unlocked(vector_impuls.data(), sizeof(Matrix3), size, file_impuls.get());
	fclose(file_impuls.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_div_impuls(fopen((main_dir + "divimpuls.bin").c_str(), "rb"), fclose);
	if (!file_div_impuls) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_div_impuls.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_div_impuls.resize(size);
	fread_unlocked(vector_div_impuls.data(), sizeof(Vector3), size, file_div_impuls.get());
	fclose(file_div_impuls.get());


	for (size_t i = 0; i < n; i++) {
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		DivStreamArray->InsertNextTuple1(vector_div_stream[i]);
		DensityArray->InsertNextTuple1(vector_density[i]);
		PressureArray->InsertNextTuple1(vector_pressure[i]);
		IllumArray->InsertNextTuple1(vector_illum[i]);  // по первому направлению		
		StreamArray->InsertNextTuple3(vector_stream[i](0), vector_stream[i](1), vector_stream[i](2));
		VelocityArray->InsertNextTuple3(vector_velocity[i](0), vector_velocity[i](1), vector_velocity[i](2));
		DivImpulsArray->InsertNextTuple3(vector_div_impuls[i](0), vector_div_impuls[i](1), vector_div_impuls[i](2));
		ImpulsArray->InsertNextTuple9(vector_impuls[i](0, 0), vector_impuls[i](0, 1), vector_impuls[i](0, 2),
			vector_impuls[i](1, 0), vector_impuls[i](1, 1), vector_impuls[i](1, 2),
			vector_impuls[i](2, 0), vector_impuls[i](2, 1), vector_impuls[i](2, 2));
	}

	EnergyArray->SetName("energy");
	unstructured_grid->GetCellData()->AddArray(EnergyArray);

	DivStreamArray->SetName("divStream");
	unstructured_grid->GetCellData()->AddArray(DivStreamArray);

	IllumArray->SetName("illum");
	unstructured_grid->GetCellData()->AddArray(IllumArray);

	PressureArray->SetName("pressure");
	unstructured_grid->GetCellData()->AddArray(PressureArray);

	DensityArray->SetName("density");
	unstructured_grid->GetCellData()->AddArray(DensityArray);


	StreamArray->SetName("stream");
	unstructured_grid->GetCellData()->SetActiveVectors("stream");
	unstructured_grid->GetCellData()->SetVectors(StreamArray);

	VelocityArray->SetName("velocity");
	unstructured_grid->GetCellData()->SetActiveVectors("velocity");
	unstructured_grid->GetCellData()->SetVectors(VelocityArray);

	DivImpulsArray->SetName("divimpuls");
	unstructured_grid->GetCellData()->SetActiveVectors("divimpuls");
	unstructured_grid->GetCellData()->SetVectors(DivImpulsArray);

	ImpulsArray->SetName("impuls");
	unstructured_grid->GetCellData()->SetActiveTensors("impuls");
	unstructured_grid->GetCellData()->SetTensors(ImpulsArray);

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(unstructured_grid);
	writer->Write();

	return 0;
}

size_t ReBuildDataArrayFullSave (const size_t class_file_vtk, const std::string name_file_vtk, const std::string name_file_out,
	const std::string& main_dir, const bool is_print = false) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	if (ReadFileVtk(0, name_file_vtk, unstructured_grid, foo, foo, foo, is_print)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> EnergyArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> DivStreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> DensityArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> PressureArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> StreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	StreamArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> DivImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	DivImpulsArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> VelocityArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	VelocityArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> ImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	ImpulsArray->SetNumberOfComponents(9);

	std::vector<Type> vector_energy;
	std::vector<Type> vector_illum;
	std::vector<Type> vector_pressure;
	std::vector<Type> vector_density;
	std::vector<Type> vector_div_stream;
	std::vector<Vector3> vector_stream;
	std::vector<Vector3> vector_velocity;
	std::vector<Vector3> vector_div_impuls;
	std::vector<Matrix3> vector_impuls;

	int size;

#ifndef ONLY_HLLC
	std::unique_ptr<FILE, int(*)(FILE*)> file_illum(fopen((main_dir + "Illum.bin").c_str(), "rb"), fclose);
	if (!file_illum) { 
		printf("file_data is not opened for reading\n"); 
		vector_illum.resize(n, 0);
		//return 1; 
	}
	else 
	{
		fread_unlocked(&size, sizeof(int), 1, file_illum.get());
		//if (n != size) { printf("Error size data\n"); return 1; }
		vector_illum.resize(size);
		fread_unlocked(vector_illum.data(), sizeof(Type), size, file_illum.get());
		fclose(file_illum.get());
	}

	std::unique_ptr<FILE, int(*)(FILE*)> file_divstream(fopen((main_dir + "divstream.bin").c_str(), "rb"), fclose);
	if (!file_divstream)
	{
		printf("file_data is not opened for reading\n"); 
		vector_div_stream.resize(n, 0);
		//return 1;
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_divstream.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_div_stream.resize(size);
		fread_unlocked(vector_div_stream.data(), sizeof(Type), size, file_divstream.get());
		fclose(file_divstream.get());
	}

	std::unique_ptr<FILE, int(*)(FILE*)> file_energy(fopen((main_dir + "energy.bin").c_str(), "rb"), fclose);
	if (!file_energy) 
	{
		printf("file_data is not opened for reading\n");
		vector_energy.resize(n, 0);
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_energy.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_energy.resize(size);
		fread_unlocked(vector_energy.data(), sizeof(Type), size, file_energy.get());
		fclose(file_energy.get());
	}
#endif

	std::unique_ptr<FILE, int(*)(FILE*)> file_density(fopen((main_dir + "density.bin").c_str(), "rb"), fclose);
	if (!file_density) 
	{ 
		printf("file_data is not opened for reading\n"); 
		vector_density.resize(n,0);
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_density.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_density.resize(size);
		fread_unlocked(vector_density.data(), sizeof(Type), size, file_density.get());
		fclose(file_density.get());
	}

	std::unique_ptr<FILE, int(*)(FILE*)> file_pressure(fopen((main_dir + "pressure.bin").c_str(), "rb"), fclose);
	if (!file_pressure) 
	{
		printf("file_data is not opened for reading\n"); 
		vector_pressure.resize(n,0);
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_pressure.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_pressure.resize(size);
		fread_unlocked(vector_pressure.data(), sizeof(Type), size, file_pressure.get());
		fclose(file_pressure.get());
	}

#ifndef ONLY_HLLC
	std::unique_ptr<FILE, int(*)(FILE*)> file_stream(fopen((main_dir + "stream.bin").c_str(), "rb"), fclose);
	if (!file_stream)
	{
		printf("file_data is not opened for reading\n"); 
		vector_stream.resize(n,Vector3::Zero());
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_stream.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_stream.resize(size);
		fread_unlocked(vector_stream.data(), sizeof(Vector3), size, file_stream.get());
		fclose(file_stream.get());
	}
#endif

	std::unique_ptr<FILE, int(*)(FILE*)> file_velocity(fopen((main_dir + "velocity.bin").c_str(), "rb"), fclose);
	if (!file_velocity) 
	{
		printf("file_data is not opened for reading\n"); 
		vector_velocity.resize(n, Vector3::Zero());
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_velocity.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_velocity.resize(size);
		fread_unlocked(vector_velocity.data(), sizeof(Vector3), size, file_velocity.get());
		fclose(file_velocity.get());
	}

#ifndef ONLY_HLLC
	std::unique_ptr<FILE, int(*)(FILE*)> file_impuls(fopen((main_dir + "impuls.bin").c_str(), "rb"), fclose);
	if (!file_impuls) 
	{ 
		printf("file_data is not opened for reading\n");
		vector_impuls.resize(n, Matrix3::Zero());
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_impuls.get());
		if (n != size) { printf("Error size data (%d - %d)\n", n, size); return 1; }
		vector_impuls.resize(size);
		fread_unlocked(vector_impuls.data(), sizeof(Matrix3), size, file_impuls.get());
		fclose(file_impuls.get());
	}

	std::unique_ptr<FILE, int(*)(FILE*)> file_div_impuls(fopen((main_dir + "divimpuls.bin").c_str(), "rb"), fclose);
	if (!file_div_impuls) 
	{
		printf("file_data is not opened for reading\n"); 
		vector_div_impuls.resize(n, Vector3::Zero());
		//return 1; 
	}
	else
	{
		fread_unlocked(&size, sizeof(int), 1, file_div_impuls.get());
		if (n != size) { printf("Error size data\n"); return 1; }
		vector_div_impuls.resize(size);
		fread_unlocked(vector_div_impuls.data(), sizeof(Vector3), size, file_div_impuls.get());
		fclose(file_div_impuls.get());
	}
#endif

	for (size_t i = 0; i < n; i++) {
		DensityArray->InsertNextTuple1(vector_density[i]);
		PressureArray->InsertNextTuple1(vector_pressure[i]);
		VelocityArray->InsertNextTuple3(vector_velocity[i](0), vector_velocity[i](1), vector_velocity[i](2));

#ifndef ONLY_HLLC
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		DivStreamArray->InsertNextTuple1(vector_div_stream[i]);

		IllumArray->InsertNextTuple1(vector_illum[i]);  // по первому направлению		
		StreamArray->InsertNextTuple3(vector_stream[i](0), vector_stream[i](1), vector_stream[i](2));		
		DivImpulsArray->InsertNextTuple3(vector_div_impuls[i](0), vector_div_impuls[i](1), vector_div_impuls[i](2));
		ImpulsArray->InsertNextTuple9(vector_impuls[i](0, 0), vector_impuls[i](0, 1), vector_impuls[i](0, 2),
			vector_impuls[i](1, 0), vector_impuls[i](1, 1), vector_impuls[i](1, 2),
			vector_impuls[i](2, 0), vector_impuls[i](2, 1), vector_impuls[i](2, 2));
#endif
	}
#ifndef ONLY_HLLC
	EnergyArray->SetName("energy");
	unstructured_grid->GetCellData()->AddArray(EnergyArray);

	DivStreamArray->SetName("divStream");
	unstructured_grid->GetCellData()->AddArray(DivStreamArray);

	IllumArray->SetName("illum");
	unstructured_grid->GetCellData()->AddArray(IllumArray);

	StreamArray->SetName("stream");
	unstructured_grid->GetCellData()->SetActiveVectors("stream");
	unstructured_grid->GetCellData()->SetVectors(StreamArray);

	DivImpulsArray->SetName("divimpuls");
	unstructured_grid->GetCellData()->SetActiveVectors("divimpuls");
	unstructured_grid->GetCellData()->SetVectors(DivImpulsArray);

	ImpulsArray->SetName("impuls");
	unstructured_grid->GetCellData()->SetActiveTensors("impuls");
	unstructured_grid->GetCellData()->SetTensors(ImpulsArray);

#endif

	PressureArray->SetName("pressure");
	unstructured_grid->GetCellData()->AddArray(PressureArray);

	DensityArray->SetName("density");
	unstructured_grid->GetCellData()->AddArray(DensityArray);

	VelocityArray->SetName("velocity");
	unstructured_grid->GetCellData()->SetActiveVectors("velocity");
	unstructured_grid->GetCellData()->SetVectors(VelocityArray);

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(unstructured_grid);
	writer->Write();

	return 0;
}
#endif

int main(int argc, char* argv[])
{

	/*
		ReadData (vtk or not)

		Rebuild/ rewrite section

		BuildHLLC struct

		Solve HLLC

		Solve illum -> divW, U, divT

		Solve HLLC + Illum

		next time step

	*/

	std::string name_file_settings = "";
	int max_number_of_iter = 1;
	Type accuracy = 1e-5;
	bool use_cuda =  true;
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
#ifdef WRITE_LOG
	ofstream ofile;
	ofile.open(main_dir + "File_with_Logs_solve.txt");
	ofile.close();
#endif

	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, main_dir, solve_direction))
		ERR_RETURN("Error reading the start settings\n")

#ifdef WRITE_LOG		
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "start settings\n";
	ofile.close();
#endif


#if defined HLLC_1D && !defined RHLLC_1D

	HLLC_1d(main_dir);
	printf("HLLC 1d end\n");
	return 0;
#endif

#ifdef RHLLC_1D

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

	std::remove((main_dir + "File_with_Logs_solve.txt").c_str()); // удаление существующего файла лога
	
#ifdef USE_VTK

#ifdef ReBuildSolve
	/*if (ReBuildDataArrayFull(0, name_file_vtk, out_file_grid_vtk, main_dir, true))
		printf("Error ReBuild\n");	
	return 0;	*/

	for (size_t i = 0; i < max_number_of_iter; i++)
	{
		std::string file_solve = solve_direction + "Solve" + std::to_string(i);
		std::string vtk_file_solve = solve_direction + "Solve" + std::to_string(i) + ".vtk";
		if (ReBuildDataArrayFullSave(0, name_file_vtk, vtk_file_solve, file_solve, true))
		{
			printf("Error ReBuild\n");
			continue;
			//ERR_RETURN("Error ReBuild\n");
		}
			

		std::remove((file_solve + "Illum.bin").c_str());
		std::remove((file_solve + "energy.bin").c_str());
		std::remove((file_solve + "stream.bin").c_str());
		std::remove((file_solve + "impuls.bin").c_str());
		std::remove((file_solve + "divstream.bin").c_str());
		std::remove((file_solve + "divimpuls.bin").c_str());

		std::remove((file_solve + "density.bin").c_str());
		std::remove((file_solve + "pressure.bin").c_str());
		std::remove((file_solve + "velocity.bin").c_str());

		printf("grid %d\n", i);
	}
	return 0;
	//------------------------------------------------------------------------------------

	std::string file_start = "D:\\Desktop\\FilesCourse\\HLLC\\Build1d.txt";

	std::string make_1d_programm = "D:\\MiniProgramm\\Make1dSolveAndTrace\\Run.exe";
	std::string main_solve_file = "D:\\Desktop\\FilesCourse\\HLLC\\Solve";
	std::string main_res_file = "D:\\Desktop\\FilesCourse\\HLLC\\SolveD";
	Vector3 direction(1, 0, 0);
	Vector3 start_p(0, 0.05, 0.05);
	std::string name_data = "density";
	bool delete_vtk = false;


	std::ifstream ifile(file_start);
	if (!ifile.is_open()) ERR_RETURN("File build wasn't opened\n");

	getline(ifile, make_1d_programm);
	getline(ifile, main_solve_file);
	getline(ifile, main_res_file);
	getline(ifile, name_data);
	ifile >> direction[0] >> direction[1] >> direction[2];
	ifile >> start_p[0] >> start_p[1] >> start_p[2];

	ifile >> delete_vtk;
	ifile.close();
	
	for (size_t i = 0; i < max_number_of_iter; i++)
	{
		std::string vtk_file_solve = main_solve_file + std::to_string(i) + ".vtk";
		std::string file_out = main_res_file + std::to_string(i);

		char buf[512];

		sprintf_s(buf, "%s %s %s %s %lf %lf %lf %lf %lf %lf", make_1d_programm.c_str(), vtk_file_solve.c_str(), file_out.c_str(), name_data.c_str(),
			direction[0], direction[1], direction[2], start_p[0], start_p[1], start_p[2]);
		//printf(buf);

		system(buf);

		//std::remove((file_out + "trace.txt").c_str());
		std::remove((file_out + "trace.vtk").c_str());
		std::remove((file_out + "value.txt").c_str());

		//if (delete_vtk) std::remove(vtk_file_solve.c_str());
	}
	//------------------------------------------------------------------------------------
	
	std::vector<int> trace;
	{
		std::string file_trace = main_res_file + std::to_string(0) + "trace.txt";

		std::ifstream ifile(file_trace);
		if (!ifile.is_open()) ERR_RETURN("File file_trace(i) wasn't opened\n");

		std::string str;

		int a;
		while (!ifile.eof()) {			
			ifile >> a;
			trace.push_back(a);
		}
		ifile.close();

		std::remove((file_trace).c_str());
		std::vector<Vector3> centers;
		ReadCentersOfTetra(name_file_centers, centers);

		vtkDataArray* density;
		vtkDataArray* pressure;
		vtkDataArray* velocity;	

		vtkSmartPointer<vtkUnstructuredGrid> u_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
		
		std::ofstream ofileP(main_solve_file + "P.txt");
		std::ofstream ofileD(main_solve_file + "D.txt");
		std::ofstream ofileV(main_solve_file + "V.txt");
		if (!ofileP.is_open()) ERR_RETURN("File result wasn't opened\n");

		for (size_t i = 0; i < max_number_of_iter; i++)
		{
			std::string vtk_file_solve = main_solve_file + std::to_string(i) + ".vtk";
			ReadFileVtk(5, vtk_file_solve, u_grid, density, pressure, velocity, false);
			for(auto id: trace)		
			{
				const Type x = centers[id][0];
				ofileP << x << " " << pressure->GetTuple1(id) << '\n';
				ofileD << x << " " << density->GetTuple1(id) << '\n';				
				Vector3 v(velocity->GetTuple3(id));
				ofileV << x << " " << v[0] << '\n';
			}

			std::remove((main_res_file + std::to_string(i) + "trace.txt").c_str());
			
		}

		ofileP.close();
		ofileD.close();
		ofileV.close();
	}

	return 0;
	ofile.open(main_res_file+".txt");
	if (!ofile.is_open()) ERR_RETURN("File result wasn't opened\n");

	for (size_t i = 0; i < max_number_of_iter; i++)
	{	
		std::string file_out = main_res_file + std::to_string(i) + "value.txt";
		
		std::ifstream ifile(file_out);
		if (!ifile.is_open()) ERR_RETURN("File solve(i) wasn't opened\n");

		std::string str;

		while (!ifile.eof())
		{
			getline(ifile, str);
			ofile << str << '\n';
		}

		ifile.close();

		std::remove(file_out.c_str());
	}

	ofile.close();
	
	return 0;
#else 
#ifdef ReWriteSolve

	if (ReWriteDataArray(class_file_vtk, name_file_vtk, main_dir, size_grid, true))
	{
		std::cout << "Error rewrite the data array vtk\n";
		return 1;
	}
	return 0;
#endif
#endif

#else
	Type _clock = -omp_get_wtime();
	if (ReadDataArray(class_file_vtk, main_dir, density, absorp_coef, rad_en_loose_rate, velocity, pressure, size_grid, true))
		ERR_RETURN("Error reading the data array vtk\n")

#ifdef WRITE_LOG		
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "reading the data array\n";
	ofile.close();
#endif
	//---------------Начальные данные для HLLC---------------------------------//
	std::vector<Vector3> centers;
#ifndef ONLY_ILLUM
	{	
		//std::vector<Vector3> centers;
		ReadCentersOfTetra(name_file_centers, centers);
		size_grid = centers.size();
		int n = size_grid;
		density.resize(n);
		pressure.resize(n);
		velocity.resize(n);

#if defined Cylinder 
		for (size_t i = 0; i < size_grid; i++)
		{
			Vector3 x = centers[i];
			if (Vector2(x[1], x[2]).norm() < 0.3 && x[0] < 0.5)
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
			} else if (centers[i][0] < 0 && centers[i][2] > 0)
			{
				density[i] = 0.1;
				pressure[i] = 1;
				velocity[i] = Vector3(0.99, 0, 0);
			}else if (centers[i][0] < 0 && centers[i][2] < 0)
			{
				density[i] = 0.5;
				pressure[i] = 1;
				velocity[i] = Vector3(0, 0, 0);
			}else if (centers[i][0] > 0 && centers[i][2] < 0)
			{
				density[i] = 0.1;
				pressure[i] = 1;
				velocity[i] = Vector3(0, 0, 0.99);
			}
		}
#endif
#if defined Cube  && !defined RHLLC//SODA
		for (size_t i = 0; i < size_grid; i++)
		{
			if (centers[i][0] < 0.5) {
				density[i] = 1;
				pressure[i] = 1;
				velocity[i] = Vector3(0, 0, 0);
			}
			else {
				density[i] = 0.125;
				pressure[i] = 0.1;
				velocity[i] = Vector3(0, 0, 0);
			}

			//Type x = centers[i][0];
			//if (x < 0.5)
			//{
			//	density[i] = 2;
			//	pressure[i] = 25;
			//	velocity[i] = Vector3(0, 0, 0);
			//}
			//else
			//{
			//	density[i] = 1;
			//	pressure[i] = 1;
			//	velocity[i] = Vector3(0, 0, 0);
			//}

		}
#endif
	}
#endif 
		_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";

#ifdef WRITE_LOG		
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "Init_hllc\n";
	ofile.close();
#endif

	//------------------------------ Illume section-----------------------------
	vector<Vector3> directions;
	vector<Type> squares;

	_clock = -omp_get_wtime();
	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface))
		ERR_RETURN("Error reading the sphere direction\n")
		_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

#ifdef WRITE_LOG		
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "reading sphere_direction\n";
	ofile.close();
#endif

	const int count_directions = directions.size();
	const int count_cells = size_grid;
	printf("\nReal size task = %d\n", size_grid);

	std::vector<int> neighbours_id_faces;
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

	std::vector<Normals> normals;
	std::vector<Type> squares_cell;
	std::vector<Type> volume;

	std::vector<int> inner_bound;

	vector<int> sorted_id_cell(count_cells* count_directions); 	// Упорядоченные индексы ячеек по данному направлению

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face, inclined_face_inverse,
		straight_face_inverse);

	if (ReadNormalFile(name_file_normals, normals)) ERR_RETURN("Error reading file normals\n")
	if (ReadSimpleFileBin(name_file_squares, squares_cell)) ERR_RETURN("Error reading file squares_cell\n")
	if (ReadSimpleFileBin(name_file_volume, volume)) ERR_RETURN("Error reading file volume\n")

	if (ReadSimpleFileBin(name_file_id_neighbors + ".bin", neighbours_id_faces)) ERR_RETURN("Error reading file neighbours\n")
	//if (ReadSimpleFile("D:\\Desktop\\FilesCourse\\graphSet\\inner_bound.txt", inner_bound)) ERR_RETURN("Error reading file inner_bound\n")
		
#ifdef WRITE_LOG		
		ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "reading simple files\n";
	ofile.close();
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

#ifdef WRITE_LOG		
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "reading  data_grid file\n";
	ofile.close();
#endif

#ifdef SAVE_DUMP_HLLC
	size_dump = (size_grid * 5 + 1) * 5;
	dump.resize(size_dump);
	w_ptr = 0;
	start_section = -666;
	bad_hllc_flag = false;
#endif

	std::vector<Type> Illum; // (4 * count_cells * count_directions, 0);
	std::vector<Type> int_scattering; // (count_cells* count_directions, 0);

	vector<Type> energy; // (size_grid, 0);
	vector<Vector3> stream; //(size_grid, Vector3(0,0,0));

	//vector<Type> prev_energy(size_grid, 0);
	//vector<Vector3> prev_stream(size_grid, Vector3(0, 0, 0));

	vector<Matrix3> impuls; //(size_grid, Matrix3::Zero());
	std::vector<Type> div_stream; //(size_grid, 0);
	vector<Vector3> div_impuls; //(size_grid, Vector3(0,0,0));
#ifndef ONLY_HLLC  
	Illum.resize(4 * count_cells * count_directions, 0);
	int_scattering.resize(count_cells* count_directions, 0);

	energy.resize(size_grid, 0);
	stream.resize(size_grid, Vector3(0, 0, 0));
	
	impuls.resize(size_grid, Matrix3::Zero());
	div_stream.resize(size_grid, 0);
	div_impuls.resize(size_grid, Vector3(0, 0, 0));
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
#ifdef ONLY_HLLC
	Type tau = 1e-5;
	Type CFL = 0.1;
	Type print_timer = 0.015;
#else
	Type tau = 1e-8;
	Type CFL = 0.001;
	Type print_timer = 1e-5;
#endif

	Type t = 0.0;	
	Type T =  0.7;
#if defined Sphere
	const Type h = 0.0128079221422811;
#elif defined Cone && !defined Jet
	//const Type h = 0.0032396496530313; // test grid
	const Type h = 0.0022530803678209; //grid: 112023 cells
#elif defined Jet
	const Type h = 0.0012548169651948;
#elif defined Cube
	const Type h = 0.0007123669658939; // Soda1d_2
	//const Type h  = 0.0010828369115320; // Soda1d
	//const Type h = 0.0010307259619874; // Soda1d_3
#elif defined Step
	const Type h = 0.0018751819368151;
	T = 5;
	CFL = 0.5;
	print_timer = 0.1;
#elif defined Cylinder
	const Type h = 0.0066864401941437;
	T = 10;
	CFL = 0.05;
	print_timer = 0.1;
#else
	const Type h = 1;
	printf("Bad geometry define\n");
#endif	
	

	std::vector<VectorX> U_full_prev;
#ifndef ONLY_ILLUM
	U_full.resize(size_grid, VectorX(5));	//  может не пройти инициализация VectorX!!
#if defined RHLLC
	ReBuildDataForHLLCRel(size_grid, U_full_prev);
#else
	if (ReBuildDataForHLLC(size_grid, U_full_prev)) ERR_RETURN("Error rebuild data to HLLC\n")
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
	std::vector<int> bound_cells;
#if 0 // ONLY_ILLUM (отдельный массив по границам. Должно обрабатываться в потоках)
	for (int i = 0; i < neighbours_id_faces.size(); i++)
	{
		if (neighbours_id_faces[i] < 0) // все границы. 
		{
			if (  (normals[i / 4].n[i % 4] - Vector3(1, 0, 0)).norm() > 1e-2  
				 && (normals[i / 4].n[i % 4] -Vector3(-1, 0, 0)).norm() > 1e-2
			   )  // возможно это условие надо использовать далее
				bound_cells.push_back(i / 4);
		}
	}
	std::unique(bound_cells.begin(), bound_cells.end());
#endif

	//------------------------------------------------------------------------------------------------------------
	int res_count = 0; // счётчик решений	
	Type full_time = -omp_get_wtime();
	while (t < T)
	{		
		BoundDataToHLLC(normals, neighbours_id_faces, centers, U_full_prev);  //ГУ для HLLC

		if (cur_timer >= print_timer)
		{
			std::string file_solve = solve_direction +"Solve" + std::to_string(res_count); //(count / step);
			WriteStepTimeSolve(file_solve, 0, Illum, energy, stream, impuls, div_stream, div_impuls, U_full_prev);
			cur_timer = 0;
			res_count++;

			printf("\n t= %f,  tau= %lf,  res_step= %d\n", t, tau, res_count);			
		}
#ifdef RHLLC
		HLLC_Rel(tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);
#else
		HLLC(size_grid, tau, bound_cells, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);
#endif

#ifdef SAVE_DUMP_HLLC
		dump[w_ptr] = start_section--; w_ptr = (w_ptr + 1) % size_dump;
		for (size_t i = 0; i < size_grid; i++)
		{
			dump[w_ptr] = U_full[i][0]; w_ptr = (w_ptr + 1) % size_dump;
			dump[w_ptr] = U_full[i][1]; w_ptr = (w_ptr + 1) % size_dump;
			dump[w_ptr] = U_full[i][2]; w_ptr = (w_ptr + 1) % size_dump;
			dump[w_ptr] = U_full[i][3]; w_ptr = (w_ptr + 1) % size_dump;
			dump[w_ptr] = U_full[i][4]; w_ptr = (w_ptr + 1) % size_dump;
		}
		dump[w_ptr] = tau; w_ptr = (w_ptr + 1) % size_dump;

		if (bad_hllc_flag)
		{
			WriteDump(main_dir + "dump.bin", dump);

			printf("start_section= %d\n", start_section+1);
			printf("Error hllc\n");
			return 1;
		}
		bad_hllc_flag = false;
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
#ifndef RHLLC
		U_full_prev.swap(U_full);
#endif
		t += tau;
		cur_timer += tau;
		
		// ---------------формирование шага по времени--------------------//
#ifdef RHLLC
		tau = FormTimeStepToRHLLC(size_grid, h, CFL);
#else
		tau = FormTimeStepToHLLC(size_grid, h, CFL, U_full_prev); 			
#endif
		if (tau < 0)
		{
			printf("Error tau = %lf\n", tau);
			break;
		}

#ifdef WRITE_LOG
		ofstream ofile;
		ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
		ofile << "\nt= " << t << "; tau= " << tau << "; step= " << res_count << '\n';
		ofile.close();
#endif
	
	}// while(t < T)
	
	full_time += omp_get_wtime();
	printf("Time while(): %f\n", full_time);
	
	if (use_cuda) 
	{
		ClearDevice(cuda_mod);
	}

	WriteStepTimeSolve(solve_direction + "Solve" + std::to_string(res_count), 0, Illum, energy, stream, impuls, div_stream, div_impuls, U_full_prev);
#endif // USE_VTK

	std::cout << "End programm\n";


	//ReBuild(name_file_sphere_direction, main_dir);

	return 0;
}

int WriteStepTimeSolve(const std::string& main_dir, const int num_dir_to_illum,
	const std::vector<Type>& Illum, const std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	const std::vector<Type>& div_stream, const vector<Vector3>& div_impuls, const std::vector<VectorX>& U_full) {

	std::vector<Type> illum;
	
#ifndef ONLY_HLLC
	illum.resize(size_grid, 0);
	const int i = num_dir_to_illum;
		for (size_t j = 0; j < size_grid; j++)
		{
			const int N = i * 4 * size_grid + j * 4;
			Type I = (Illum[N] + Illum[N + 1] + Illum[N + 2] + Illum[N + 3]) / 4;
			illum[j] = I;
		}
#endif

#ifdef USE_VTK
	//WriteFileSolution(out_file_grid_vtk, illum, energy, stream, impuls, name_file_vtk);
#else
	WriteFileSolution(main_dir, illum, energy, stream, impuls, div_stream, div_impuls, U_full);
#endif

	return 0;
}

int SolveIllumAndHLLC(const int N,const Type tau, std::vector<VectorX>& U_full,
	const vector<Type>& prev_energy, const vector<Type>& energy, 
	const vector<Vector3>& prev_stream, const vector<Vector3>& stream,
	const std::vector<Type>& div_stream, const  vector<Vector3>& div_impuls) {
	
	const Type c = 3 * 1e+8;  // скорость света
	// U[i] = tau*BU + U  // нет зависимости от соседей. Делаем U+= 
	for (size_t i = 0; i < N; i++)
	{
		//U_full[i][0] += 0;

		U_full[i][1] += tau * (-div_impuls[i][0]/c /*- (stream[i][0] - prev_stream[i][0]) / tau / c / c*/);  // +F_g //vx
		U_full[i][2] += tau * (-div_impuls[i][1]/c /*- (stream[i][1] - prev_stream[i][1]) / tau / c / c*/);  //vy
		U_full[i][3] += tau * (-div_impuls[i][2] /c /*- (stream[i][2] - prev_stream[i][2]) / tau / c / c*/);  //vz

		U_full[i][4] += tau * (/*-(energy[i] - prev_energy[i]) / tau / c*/ - div_stream[i]); //+F_g.dot(vel)  //e
	}
	
	return 0;
}

int CalculateIllumParam(const int count_cells, const int count_directions,
	const std::vector<Vector3>& directions, const std::vector<Type>& squares,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Type>& Illum, std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	std::vector<Type>& div_stream, vector<Vector3>& div_impuls) {

	std::vector<Type> illum; 
	MakeIllum(count_cells, count_directions, Illum, illum);
	MakeEnergy(illum, squares, square_surface, energy);

	MakeDivStream(Illum, directions, squares, square_surface, normals, squares_cell, volume, stream, div_stream);
	MakeDivImpuls(Illum, directions, squares, square_surface, normals, squares_cell, volume, impuls, div_impuls);

	return 0;
}
int CalculateIllumParamOptMemory(const int count_cells, const int count_directions,
	const std::vector<Vector3>& directions, const std::vector<Type>& squares,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Type>& Illum, std::vector<Type>& energy, std::vector<Vector3>& stream, vector<Matrix3>& impuls,
	std::vector<Type>& div_stream, vector<Vector3>& div_impuls) {

	//std::vector<Type> illum;
	//MakeIllum(count_cells, count_directions, Illum, illum);

	Type* illum = new Type[count_cells * count_directions];

	FILE* f;
	f = fopen("D:\\Desktop\\FilesCourse\\Test126\\FormIllum.bin", "rb");
	fread(illum, sizeof(Type), count_cells * count_directions, f);
	fclose(f);

	MakeEnergy(illum, squares, square_surface, energy);

	//MakeDivStream(Illum, directions, squares, square_surface, normals, squares_cell, volume, stream, div_stream);
	//MakeDivImpuls(Illum, directions, squares, square_surface, normals, squares_cell, volume, impuls, div_impuls);

	return 0;
}

int MakeIllum(const int count_cells, const int count_directions, const std::vector<Type>& Illum, std::vector<Type>& illum) {

	if(illum.size() != (count_cells * count_directions))
		illum.resize(count_cells * count_directions);

	for (size_t i = 0; i < count_directions; i++)
		for (size_t j = 0; j < count_cells; j++)
		{
			const int N = i * 4 * count_cells + j * 4;
			Type I = (Illum[N] + Illum[N + 1] + Illum[N + 2] + Illum[N + 3]) / 4;
			illum[i * count_cells + j] = I;
		}
	return 0;
}

int MakeDivStream(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume, vector<Vector3>& stream,
	vector<Type>& div_stream) {

	
	std::vector<Vector3> Stream(size_grid * 4);
	MakeStream(Illum, directions, squares, square_surface, Stream);

	if(stream.size() != size_grid)
		stream.resize(size_grid);
	for (size_t i = 0; i < size_grid; i++)
	{
		stream[i] = (Stream[i * 4] + Stream[i * 4 + 1] + Stream[i * 4 + 2] + Stream[i * 4 + 3]) / 4;
	}

	{
			/*std::vector<Vector3> points;
			ReadSimpleFileBin("D:\\Desktop\\FilesCourse\\IllumGrid\\center_face.bin", points);
			std::vector<Vector3> centers;
			ReadCentersOfTetra("D:\\Desktop\\FilesCourse\\IllumGrid\\centers.bin", centers);
			Vector3 Stream;
			div_stream.resize(size_grid, 0);
			for (size_t i = 0; i < size_grid; i++) 
			{
				for (size_t j = 0; j < 4; j++) {
					Vector3 x = points[i * 4 + j];
					Stream = Vector3(sin(x[0]),0,0);
					div_stream[i] += normals[i].n[j].dot(Stream) * squares_cell[i * 4 + j];
				}
				div_stream[i] /= volume[i];
			}

			Type norm = 0;
			for (size_t i = 0; i < size_grid; i++) 
			{
				norm += (cos(centers[i][0]) - div_stream[i]) * (cos(centers[i][0]) - div_stream[i]);
			}
			norm = sqrt(norm) / size_grid;
			printf("Error div= %f\n", norm);

			return 0;*/
	}

	if (div_stream.size() != size_grid)
		div_stream.resize(size_grid);

	for (size_t i = 0; i < size_grid; i++)
	{
		div_stream[i] = 0;
		for (size_t j = 0; j < 4; j++)
			div_stream[i] += normals[i].n[j].dot(Stream[i * 4 + j]) * squares_cell[i * 4 + j];
		div_stream[i] /= volume[i];
	}

	Stream.clear();
	return 0;
}

int MakeDivImpuls(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume, vector<Matrix3>& impuls,
	vector<Vector3>& div_impuls) {

	std::vector<Matrix3> Impuls(size_grid * 4);
	MakeImpuls(Illum, directions, squares, square_surface, Impuls);  // проверить для расчета по граням

	if (impuls.size() != size_grid)
		impuls.resize(size_grid);
	for (size_t i = 0; i < size_grid; i++)
	{
		impuls[i] = (Impuls[i * 4] + Impuls[i * 4 + 1] + Impuls[i * 4 + 2] + Impuls[i * 4 + 3]) / 4;
	}
	
	/*for (size_t i = 0; i < size_grid; i++)
	{
		div_stream[i] = 0;
		for (size_t j = 0; j < 4; j++)
			div_stream[i] += normals[i].n[j].dot(Stream[i * 4 + j]) * squares_cell[i * 4 + j];
		div_stream[i] /= volume[i];
	}*/

	if (div_impuls.size() != size_grid)
		div_impuls.resize(size_grid);

	for (size_t i = 0; i < size_grid; i++)
	{
		div_impuls[i] = Vector3::Zero();
		for (size_t j = 0; j < 4; j++) 
		{
			div_impuls[i] += (Impuls[i * 4 + j] * (normals[i].n[j])) * squares_cell[i * 4 + j];

			//div_stream[i] += normals[i].n[j].dot(Stream[i * 4 + j]) * squares_cell[i * 4 + j];
		//	div_impuls[i][0] += Impuls[i * 4 + j].row(0).dot(normals[i].n[j]) * squares_cell[i * 4 + j];
		//	div_impuls[i][1] += Impuls[i * 4 + j].row(1).dot(normals[i].n[j]) * squares_cell[i * 4 + j];
		//	div_impuls[i][2] += Impuls[i * 4 + j].row(2).dot(normals[i].n[j]) * squares_cell[i * 4 + j];
		}
		div_impuls[i] /= volume[i];
	}

//	impuls.clear();
	return 0;
}









Type IntegarteDirection(const int num_cell, const std::vector<std::vector<Type>>& Illum, const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();  // number of directions
	int m = size_grid; // Illum.size() / n;  // number of cells

	for (size_t i = 0; i < n; i++) {
		res += Illum[i][num_cell] * squares[i];
	}

	return res / scuare_surface;
}
size_t MakeEnergy(const std::vector<std::vector<Type>>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i)
	{
		energy[i] = IntegarteDirection(i, Illum, squares, scuare_surface);
	}

	return 0;
}

int ReBuild(std::string name_file_sphere_direction, std::string main_dir) {

	//size_grid = 8392463;
	size_grid = 20475;
	vector<Vector3> directions;
	vector<Type> squares;
	vector<Type> energy(size_grid);

	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface))
		ERR_RETURN("Error reading the sphere direction\n")

		printf("Here\n");
	
	FILE* f;
	f = fopen((main_dir + "FormIllum.bin").c_str(), "rb");
	std::vector < std::vector < Type>>illum(directions.size());
	for (size_t i = 0; i < directions.size(); i++)
	{
		illum[i].resize(size_grid);
		fread(illum[i].data(), sizeof(Type), size_grid, f);
	}
	fclose(f);


	cout << "EndRead" << '\n';

	MakeEnergy(illum, squares, square_surface, energy);

	cout << energy.size() << '\n';

	FILE* fe;
	fe = fopen((main_dir + "Solve0Illum.bin").c_str(), "wb");
	fwrite(&size_grid, sizeof(int), 1, fe);
	fwrite(illum[0].data(), sizeof(Type), illum[0].size(), fe);
	fclose(fe);

	fe = fopen((main_dir + "Solve0pressure.bin").c_str(), "wb");
	fwrite(&size_grid, sizeof(int), 1, fe);
	fwrite(illum[70].data(), sizeof(Type), illum[70].size(), fe);
	fclose(fe);

	fe = fopen((main_dir + "Solve0density.bin").c_str(), "wb");
	fwrite(&size_grid, sizeof(int), 1, fe);
	fwrite(illum[125].data(), sizeof(Type), illum[125].size(), fe);
	fclose(fe);

	fe = fopen((main_dir + "Solve0energy.bin").c_str(), "wb");
	int n = energy.size();
	fwrite(&n, sizeof(int), 1, fe);
	fwrite(energy.data(), sizeof(Type), size_grid, fe);
	fclose(fe);

	return 0;
}