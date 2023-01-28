#include "rebuild_grid.h"

#if defined UTILS

#include "../file_module/reader_bin.h"
#include "../file_module/reader_vtk.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_vtk.h"
#include "../file_module/writer_bin.h"

static int ReadNetgenGrid(file_name name_file_in, std::vector<Vector3>& point, std::vector<Eigen::Vector4i>& cell, const int size = 3)
{
	std::ifstream ifile;
	OPEN_FSTREAM(ifile, name_file_in);

	int n;
	ifile >> n;
	point.resize(n);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			ifile >> point[i][j];
		}
	}

	ifile >> n;
	cell.resize(n);
	for (size_t i = 0; i < n; i++)
	{
		ifile >> cell[i][0];
		for (size_t j = 0; j < size + 1; j++)
		{
			ifile >> cell[i][j];
		}
	}
	ifile.close();

	return 0;
}
static int WriteVtkFile2d(file_name name_file_out, const std::vector<Vector3>& point, const std::vector<Eigen::Vector4i>& cell)
{
	std::ofstream ofile;
	OPEN_FSTREAM(ofile, name_file_out.c_str());

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS " << point.size() << " double\n";

	for (size_t i = 0; i < point.size(); i++)
	{
		ofile << setprecision(16) << point[i][0] << ' ' << point[i][1] << " " << 2 << '\n';
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

int ReBuildNetgenToMetis(int argc, char* argv[])
{
	std::string name_file_in = "";
	std::string name_file_out = "";
	int size = 3; // по умолчанию 3d сетка

	if (argc <= 2)
	{
		printf("Error input data!\n");
		printf("Input: path\\input_file, path\\output_file.vtk\n");
		RETURN_ERR("");

	}
	else
	{
		name_file_in = argv[1];
		name_file_out = argv[2];

		if (argc == 4)
		{
			size = std::stoi(argv[3]);
		}
	}

	std::vector<Eigen::Vector3d> point;
	std::vector<Eigen::Vector4i> cell;

	ReadNetgenGrid(name_file_in, point, cell, size);

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, name_file_out);

	ofile << cell.size() << "\n";
	for (size_t i = 0; i < cell.size(); i++)
	{
		for (size_t j = 0; j < size + 1; j++)
		{
			ofile << cell[i][j] << ' ';
		}
		ofile << '\n';
	}

	ofile.close();
	return 0;
}

int ReBuildNetgenToVTK2d(int argc, char* argv[])
{
	std::string name_file_in = "";
	std::string name_file_out = "";

	if (argc <= 2)
	{
		RETURN_ERR("Error input data!\n Input: path\\input_file, path\\output_file.vtk\n");		
	}
	else
	{
		name_file_in = argv[1];
		name_file_out = argv[2];
	}

	std::vector<Vector3> point;
	std::vector<Eigen::Vector4i> cell;

	ReadNetgenGrid(name_file_in, point, cell, 2);

	WriteVtkFile2d(name_file_out, point, cell);
	
	return 0;
}

int ReNumberingGrid(int argc, char* argv[])
{
	std::string file_grid = "";
	std::string file_new_grid = "";
	std::string file_graph = "";

	if (argc < 4)
	{		
		printf("Error input data!\n");
		RETURN_ERR("Input Format: path\\file.vtk,  path\\new_file.vtk,  path\\graph.bin\n");	
	}
	else
	{
		file_grid = argv[1];
		file_new_grid = argv[2];
		file_graph = argv[3];
	}

	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();	
	if (ReadFileVtk(file_grid, ugrid)) RETURN_ERR("Error reading grid\n");

	const int N = ugrid->GetNumberOfCells();
	std::vector<int> graph(N, 0);

	if (ReadSimpleFileBin(file_graph, graph)) RETURN_ERR("Error reading graph\n");

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, file_new_grid.c_str());	

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
	ofile << "CELLS " << N << " " << N * (n + 1) << '\n';

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

#ifdef USE_VTK
static int WriteDataToGrid(file_name name_file_grid, file_name name_file_data, file_name name_file_output, file_name name_data) 
{
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
	if (ReadFileVtk(name_file_grid, unstructured_grid))  RETURN_ERR("Error reading the file vtk\n");
	
	const int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> DataArray = vtkSmartPointer<vtkDoubleArray>::New();

	std::vector<int> vector_data;

	std::string formats = name_file_data.substr(name_file_data.find(".") + 1, 3);
	if (formats == "txt")
	{
		if (ReadSimpleFileTxt(name_file_data, vector_data)) RETURN_ERR("");
	}
	else if (formats == "bin")
	{
		if (ReadSimpleFileBin(name_file_data, vector_data))RETURN_ERR("");
	}
	else
	{
		RETURN_ERR("Error formats data file. use .bin or .txt\n");
	}

	if (vector_data.size() != n)
	{
		RETURN_ERR("Error size grid and data\n");
	}

	for (size_t i = 0; i < n; i++)
	{
		DataArray->InsertNextTuple1(vector_data[i]);
	}
	DataArray->SetName(name_data.c_str());
	unstructured_grid->GetCellData()->AddArray(DataArray);


	WriteVtkGrid(name_file_output, unstructured_grid);	
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

#endif//USE_VTK

#endif //UTILS
