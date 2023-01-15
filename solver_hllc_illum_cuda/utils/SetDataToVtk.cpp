#include "Header.h"

int FindBoundCells(std::string name_file_data, std::string name_file_output) {

	std::vector<int> vector_data;
	int size;

	std::unique_ptr<FILE, int(*)(FILE*)> file_data(fopen(name_file_data.c_str(), "rb"), fclose);
	if (!file_data) {
		printf("file_data is not opened for reading\n");
		return 1;
	}
	else
	{
		fread(&size, sizeof(int), 1, file_data.get());
		//if (n != size) { printf("Error size data\n"); return 1; }
		vector_data.resize(size);
		fread(vector_data.data(), sizeof(int), size, file_data.get());
		fclose(file_data.get());
	}

	std::vector<int> vector_out(size / 4, 0);

	for (size_t i = 0; i < size / 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			if (vector_data[i * 4 + j] < 0)
			{
				vector_out[i] += vector_data[i * 4 + j];
			}
		}
	}


	std::unique_ptr<FILE, int(*)(FILE*)> file_out(fopen(name_file_output.c_str(), "wb"), fclose);
	if (!file_out) {
		printf("file_output is not opened for reading\n");
		return 1;
	}
	else
	{
		size = vector_out.size();
		fwrite(&size, sizeof(int), 1, file_out.get());
		//if (n != size) { printf("Error size data\n"); return 1; }
		vector_out.resize(size);
		fwrite(vector_out.data(), sizeof(int), size, file_out.get());
		fclose(file_out.get());
	}

	return 0;
}


int GetAverageArea(int argc, char* argv[])
{
	std::string name_file_grid = "";
	if (argc <= 1)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input: path\\file.vtk\n");
		return 1;
#else
		name_file_grid = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.vtk";
#endif		
	}
	else
	{
		name_file_grid = argv[1];
	}

	vtkSmartPointer<vtkUnstructuredGrid>  ugrid =
		vtkSmartPointer<vtkUnstructuredGrid> ::New();

	vtkDataArray* foo;
	if (ReadFileVtk(0, name_file_grid, ugrid, foo, foo, foo, true))
	{
		printf("File vtk not opened\n");
		return 1;
	}

	Type h_min = 1000;
	Type h_mid = 0;
	const int N = ugrid->GetNumberOfCells();

	Type p[3], p1[3];
	vtkPoints* point_ptr;
	
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			point_ptr = ugrid->GetCell(i)->GetEdge(j)->GetPoints();
			point_ptr->GetPoint(0, p);
			point_ptr->GetPoint(1, p1);
						
			const Type dx = p[0] - p1[0];
			const Type dy = p[1] - p1[1];
			const Type h = sqrt(dx * dx + dy * dy);

			h_mid += h;
			if (h_min > h)
			{
				h_min = h;
			}
		}
	}

	std::cout << setprecision(16) << "Min size cell: " << h_min << '\n';
	std::cout << setprecision(16) << "Average size cell: " << h_mid/N/3 << '\n';

	return 0;
}
int ReBuildNetgen2dToVTK(int argc, char* argv[])
{
	std::string name_file_in = "";
	std::string name_file_out = "";

	if (argc <= 2)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input: path\\input_file, path\\output_file.vtk\n");
		return 1;
#else
		name_file_in = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.txt";
		name_file_out = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.vtk";
#endif
	}
	else
	{
		name_file_in = argv[1];
		name_file_out = argv[2];
	}


	std::ifstream ifile(name_file_in);

	if (!ifile.is_open())
	{
		printf("Error. Input file not opened\n");
		return 1;
	}

	int n;
	ifile >> n;
	std::vector<Eigen::Vector2d> point(n);
	for (size_t i = 0; i < n; i++)
	{
		ifile >> point[i][0] >> point[i][1];
	}

	ifile >> n;
	std::vector<Eigen::Vector3i> cell(n);
	for (size_t i = 0; i < n; i++)
	{
		ifile >> cell[i][0] >> cell[i][0] >> cell[i][1] >> cell[i][2];
	}


	std::ofstream ofile(name_file_out);

	if (!ofile.is_open())
	{
		printf("Error. Output file not opened\n\n");
		return 1;
	}

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS " << point.size() << " double\n";

	for (size_t i = 0; i < point.size(); i++)
	{
		ofile << setprecision(16) << point[i][0] << ' ' << point[i][1] <<" " << 2 << '\n';
	}

	ofile << "CELLS " << cell.size() << " " << cell.size() * 4 << '\n';

	for (size_t i = 0; i < cell.size(); i++)
	{
		ofile << 3 << " " << cell[i][0] - 1 << ' ' << cell[i][1] - 1 << ' ' << cell[i][2] - 1 << '\n';
	}

	ofile << "CELL_TYPES " << cell.size() << "\n";
	for (size_t i = 0; i < cell.size(); i++)
	{
		ofile << 5 /*VTK_TRIANGLE*/ << '\n';
	}
	
	ofile.close();
	return 0;
}


static int WriteDataToGrid(std::string name_file_grid, std::string name_file_data, std::string name_file_output,
	std::string name_data) {
	
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	if (ReadFileVtk(0, name_file_grid, unstructured_grid, foo, foo, foo, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	const int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> DataArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	std::vector<int> vector_data;
	int size;

	std::string formats = name_file_data.substr(name_file_data.find(".") + 1, 3);
	if (formats == "txt")
	{
		std::ifstream ifile(name_file_data);

		if (!ifile.is_open())
		{
			printf("file_data is not opened for reading\n");
			vector_data.resize(n, 0);
			return 1;
		}
		else
		{
			ifile >> size;
			
			vector_data.resize(size);
			for (int i = 0; i < size; i++)
			{
				ifile >> vector_data[i];
			}

			ifile.close();
		}
	}
	else if (formats == "bin")
	{

		std::unique_ptr<FILE, int(*)(FILE*)> file_data(fopen(name_file_data.c_str(), "rb"), fclose);
		if (!file_data)
		{
			printf("file_data is not opened for reading\n");
			vector_data.resize(n, 0);
			return 1;
		}
		else
		{
			fread(&size, sizeof(int), 1, file_data.get());
			vector_data.resize(size);
			fread(vector_data.data(), sizeof(int), size, file_data.get());
			fclose(file_data.get());
		}
	}
	else
	{
		printf("Error formats data file. use .bin or .txt\n");
		return 1;
	}

	if (size != n)
	{
		printf("Error size grid (%d) and data (%d)\n", n, size);
		return 1;
	}

	for (size_t i = 0; i < n; i++)
	{
		DataArray->InsertNextTuple1(vector_data[i]);
	}

	DataArray->SetName(name_data.c_str());
	unstructured_grid->GetCellData()->AddArray(DataArray);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_output.c_str());
	writer->SetInputData(unstructured_grid);
	writer->Write();

	return 0;

}

int SetScalarDataVtkFromFile(int argc, char* argv[])
{
	std::string file_grid = "";
	std::string file_data = "";
	std::string file_out = "";
	std::string name_data = "data";

	if (argc < 4)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input Format: path\\file.vtk,  path\\data.bin,  path\\outfile.vtk\n");		
		printf("Additional parameters: name_field_data");
		printf("Data format: size a b c ... \n");
		return 1;
#else
		file_grid = "D:\\Desktop\\FilesCourse\\Grids\\step.vtk";
		file_data = "D:\\Desktop\\FilesCourse\\Step\\data.bin";
		file_out = "D:\\Desktop\\FilesCourse\\Step\\out.vtk";		
#endif
	}
	else
	{
		file_grid = argv[1];
		file_data = argv[2];
		file_out = argv[3];
		
		if (argc > 4)
		{
			name_data = argv[4];
		}
	}

	if (WriteDataToGrid(file_grid, file_data, file_out, name_data))
	{
		printf("Error write data to vtk grid\n");
		return 1;
	}

	return 0;
}


int ReNumberingGrid(int argc, char* argv[])
{
	std::string file_grid = "";
	std::string file_new_grid= "";
	std::string file_graph = "";

	if (argc < 4)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input Format: path\\file.vtk,  path\\new_file.vtk,  path\\graph.bin\n");
		return 1;
#else
		file_grid = "D:\\Desktop\\FilesCourse\\Grids\\small_test_grid_54.vtk";
		file_new_grid = "D:\\Desktop\\FilesCourse\\Test\\new_small_test_grid_54.vtk";
		file_graph = "D:\\Desktop\\FilesCourse\\Test\\graph0.bin";
#endif
	}
	else
	{
		file_grid = argv[1];
		file_new_grid = argv[2];
		file_graph = argv[3];
	}

	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkDataArray* foo;
	ReadFileVtk(0, file_grid, ugrid, foo, foo, foo, true);

	const int N = ugrid->GetNumberOfCells();
	std::vector<int> graph(N, 0);
	
#if 1
	
	//try
	{
		FILE* f = fopen((file_graph).c_str(), "rb");
		fread(graph.data(), sizeof(int), N, f);
		fclose(f);
	}
//	catch(int num_ex)
	{
	//	printf("Error reading graph. Maybe false size. ()\n");
	//	return 1;
	}
#endif

	std::ofstream ofile(file_new_grid);

	if (!ofile.is_open())
	{
		printf("Error. Output file not opened\n\n");
		return 1;
	}

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS " << ugrid->GetNumberOfPoints() << " double\n";

	double p[3];
	for (size_t i = 0; i < ugrid->GetNumberOfPoints(); i++)
	{
		ugrid->GetPoint(i, p);
		ofile << setprecision(16) << p[0] << ' ' << p[1] << " " << p[2] << '\n';
	}
	
	int n = ugrid->GetCell(0)->GetPointIds()->GetNumberOfIds();
	ofile << "CELLS " << N << " " << N * (n+1) << '\n';

	for (size_t i = 0; i < N; i++)
	{
		int sort_id = graph[i];
		vtkIdList* idp = ugrid->GetCell(sort_id)->GetPointIds();

		ofile << idp->GetNumberOfIds() << " ";
		for (size_t id = 0; id < idp->GetNumberOfIds(); id++)
		{
			ofile << idp->GetId(id) << ' ';			
		}
		ofile << '\n';		
	}

	ofile << "CELL_TYPES " << N << "\n";
	for (size_t i = 0; i < N; i++)
	{
		ofile << ugrid->GetCell(graph[i])->GetCellType()  /*VTK_TRIANGLE*/ << '\n';
	}

	ofile.close();
		
	return 0;
}
