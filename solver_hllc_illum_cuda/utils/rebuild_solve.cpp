#include "rebuild_solve.h"

#ifdef UTILS
#include "../file_module/reader_bin.h"
#include "../file_module/writer_bin.h"

#if defined USE_VTK
#include "../file_module/reader_vtk.h"
#include "../file_module/writer_vtk.h"


static int SetVtkData1(file_name main_dir, file_name name_data, vtkSmartPointer<vtkUnstructuredGrid>& ugrid, const bool delete_file = true)
{
	vtkSmartPointer<vtkDoubleArray> data_vtk = vtkSmartPointer<vtkDoubleArray>::New();
	data_vtk->SetNumberOfComponents(1);

	std::vector<Type> data_bin;	

	if (ReadSimpleFileBin(main_dir + name_data + ".bin", data_bin)) return 1;
	if (data_bin.size() != ugrid->GetNumberOfCells()) RETURN_ERR("Error size data and grid\n");
	if (delete_file) std::remove((main_dir + name_data + ".bin").c_str());

	for (auto el : data_bin)
	{
		data_vtk->InsertNextTuple1(el);
	}

	data_vtk->SetName(name_data.c_str());
	
	ugrid->GetCellData()->AddArray(data_vtk);
	
	return 0;
}
static int SetVtkData3(file_name main_dir, file_name name_data, vtkSmartPointer<vtkUnstructuredGrid>& ugrid, const bool delete_file = true)
{
	vtkSmartPointer<vtkDoubleArray> data_vtk = vtkSmartPointer<vtkDoubleArray>::New();
	data_vtk->SetNumberOfComponents(3);

	std::vector<Vector3> data_bin;

	if (ReadSimpleFileBin(main_dir + name_data + ".bin", data_bin)) return 1;
	if (data_bin.size() != ugrid->GetNumberOfCells()) RETURN_ERR("Error size data and grid\n");
	if (delete_file) std::remove((main_dir + name_data + ".bin").c_str());

	for (auto el : data_bin)
	{
		data_vtk->InsertNextTuple(el.data());
	}

	data_vtk->SetName(name_data.c_str());
	ugrid->GetCellData()->SetActiveVectors(name_data.c_str());
	ugrid->GetCellData()->SetVectors(data_vtk);

	return 0;
}
static int SetVtkData9(file_name main_dir, file_name name_data, vtkSmartPointer<vtkUnstructuredGrid>& ugrid, const bool delete_file = true)
{
	vtkSmartPointer<vtkDoubleArray> data_vtk = vtkSmartPointer<vtkDoubleArray>::New();
	data_vtk->SetNumberOfComponents(9);

	std::vector<Matrix3> data_bin;

	if (ReadSimpleFileBin(main_dir + name_data + ".bin", data_bin)) return 1;
	if (data_bin.size() != ugrid->GetNumberOfCells()) RETURN_ERR("Error size data and grid\n");
	if(delete_file) std::remove((main_dir + name_data + ".bin").c_str());

	for (auto el : data_bin)
	{
		data_vtk->InsertNextTuple(el.data());
	}

	data_vtk->SetName(name_data.c_str());
	
	ugrid->GetCellData()->SetActiveTensors(name_data.c_str());
	ugrid->GetCellData()->SetTensors(data_vtk);

	return 0;
}

static size_t ReBuildDataVtkArray(file_name name_file_vtk, file_name name_file_out, file_name main_dir) 
{
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_vtk, ugrid)) RETURN_ERR("Error reading the file vtk\n");
		
#ifdef ILLUM
	SetVtkData1(main_dir, "Illum", ugrid);
	SetVtkData1(main_dir, "divstream", ugrid);
	SetVtkData1(main_dir, "energy", ugrid);

	SetVtkData3(main_dir, "stream", ugrid);
	SetVtkData3(main_dir, "divimpuls", ugrid);

	SetVtkData9(main_dir, "impuls", ugrid);
#endif

#ifdef HLLC
	SetVtkData1(main_dir, "density", ugrid);
	SetVtkData1(main_dir, "pressure", ugrid);
	SetVtkData3(main_dir, "velocity", ugrid);
#endif

	WriteVtkGrid(name_file_out, ugrid, true);
	
	return 0;
}

int rebuild_solve(file_name name_file_settings)
{
	std::string name_file_vtk;	
	std::string foo;
	std::string adress_solve;
	int max_number_of_iter;
	int class_vtk;

	if (ReadStartSettings(name_file_settings, class_vtk, name_file_vtk, foo, foo, foo, adress_solve, max_number_of_iter))
	{
		RETURN_ERR("Error reading solve settings\n");
	}

	for (int i = 0; i < max_number_of_iter; i++)
	{
		std::string file_solve = adress_solve + "Solve" + std::to_string(i);
		std::string vtk_file_solve = adress_solve + "Solve" + std::to_string(i) + ".vtk";

		if (ReBuildDataVtkArray(name_file_vtk, vtk_file_solve, file_solve))
		{
			printf("Error ReBuild\n");
			continue;
			//ERR_RETURN("Error ReBuild\n");
		}

		printf("grid %d\n", i);
	}
	return 0;
}

int rewrite_vtk_array(file_name name_file_settings)
{
	std::string name_file_vtk;
	std::string foo;
	std::string base_adress;
	int max_number_of_iter;
	int class_vtk;

	if (ReadStartSettings(name_file_settings, class_vtk, name_file_vtk, foo, foo, base_adress, foo, max_number_of_iter))
	{
		RETURN_ERR("Error reading solve settings\n");
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_vtk, unstructured_grid)) RETURN_ERR("Error reading the file vtk\n");

	const int size_grid = unstructured_grid->GetNumberOfCells();

	if (unstructured_grid->GetCellData()->GetNumberOfArrays() == 0) RETURN_ERR("grid hasn't data\n");

	std::vector<Type> v_data;
	int size;
	for (size_t i = 0; i < unstructured_grid->GetCellData()->GetNumberOfArrays(); i++)
	{
		std::string name_data(unstructured_grid->GetCellData()->GetArrayName(i));
		vtkDataArray* data = unstructured_grid->GetCellData()->GetScalars(name_data.c_str());
		size = data->GetSize() - 1;
		v_data.resize(size);		
		for (size_t i = 0; i < size; i++)
		{
			v_data[i] = data->GetTuple1(i);
		}

		WriteSimpleFileBin(base_adress + name_data + ".bin", v_data);

	}

	return 0;
}

static int Make1dFromTrace(file_name file_trace, file_name file_centers, file_name file_solve, const int max_iter)
{
	std::vector<int> trace;
	{
		std::ifstream ifile;
		OPEN_FSTREAM(ifile, file_trace.c_str());			
		int a;
		while (!ifile.eof()) 
		{
			ifile >> a;
			trace.push_back(a);
		}
		ifile.close();

		printf("trace has %d elem\n", trace.size());
		
		std::vector<Vector3> centers;
		ReadSimpleFileBin(file_centers, centers);

		vtkDataArray* density;
		vtkDataArray* pressure;
		vtkDataArray* velocity;

		vtkSmartPointer<vtkUnstructuredGrid> u_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

		std::ofstream ofileP;
		std::ofstream ofileD;
		std::ofstream ofileV;

		OPEN_FSTREAM(ofileP, (file_solve + "P.txt").c_str());
		OPEN_FSTREAM(ofileD, (file_solve + "D.txt").c_str());
		OPEN_FSTREAM(ofileV, (file_solve + "V.txt").c_str());		

		for (size_t i = 0; i < max_iter; i++)
		{
			std::string vtk_file_solve = file_solve + std::to_string(i) + ".vtk";
			ReadFileVtk(5, vtk_file_solve, u_grid, density, pressure, velocity, false);
			for (auto id : trace)
			{
				const Type x = centers[id][0];
				ofileP << x << " " << pressure->GetTuple1(id) << '\n';
				ofileD << x << " " << density->GetTuple1(id) << '\n';
				Vector3 v(velocity->GetTuple3(id));
				ofileV << x << " " << v[0] << '\n';
			}

			//std::remove((main_res_file + std::to_string(i) + "trace.txt").c_str());

		}

		ofileP.close();
		ofileD.close();
		ofileV.close();
	}

	printf("Succsess time_trace\n");
	return 0;
}

int  RunMake1d(int argc, char* argv[])
{
	std::string file_trace = "";
	std::string file_centers = "";
	std::string file_solve = "";
	int max_iter = 0;


	if (argc <= 4)
	{
		printf("Error input data!\n");
		RETURN_ERR("Input Format: path\\trace.txt, path\\centers.bin, path\\Solve_i, number_of_iter \n");		
	}
	else
	{
		file_trace = argv[1];
		file_centers = argv[2];
		file_solve = argv[3];
		max_iter = std::stoi(argv[4]);
	}

	printf("You input:\n");
	printf("file_trace:  %s\n", file_trace.c_str());
	printf("file_centers:  %s\n", file_centers.c_str());
	printf("file_solve:  %s\n", file_solve.c_str());
	printf("iter:  %d\n", max_iter);

	return Make1dFromTrace(file_trace, file_centers, file_solve, max_iter);
}

int BuildHLLC_1dTime(int argc, char* argv[])
{
	std::string file_start = "";
	int max_number_of_iter = 1;

	if (argc < 2)
	{
		printf("Error input data!\n");
		printf("Input Format: path\\file_settings.txt, number_of_data_files \n");
		printf("format file_settings: \n");

		printf("path_make_1d_programm\\Run.exe");
		printf("main_solve_files:( like: path\\Solve. In direction: Solve0,vtk, Solve1.vtk, ...)\n");
		printf("main_result_file_name\n");
		printf("name_data_vtk\n");
		printf("direction: x y z\n");
		printf("start_point x y z \n");
		printf("name_file_centers\n");
		RETURN_ERR("");
	}
	else
	{
		file_start = argv[1];
		max_number_of_iter = std::stoi(argv[2]);
	}

	//std::cout << "Max_number_of_iter= " << max_number_of_iter << '\n';
	std::string make_1d_programm; // = "D:\\MiniProgramm\\Make1dSolveAndTrace\\Run.exe";
	std::string main_solve_file;// = "D:\\Desktop\\FilesCourse\\Soda\\Solve\\Solve";
	std::string main_res_file;// = "D:\\Desktop\\FilesCourse\\Soda\\Solve\\Solve1d";
	Vector3 direction(1, 0, 0);
	Vector3 start_p(0, 0, 0);
	std::string name_data;// = "density";	
	std::string  name_file_centers;// = "D:\\Desktop\\FilesCourse\\Soda\\centers.bin";

	std::ifstream ifile;
	OPEN_FSTREAM(ifile, file_start.c_str());

	getline(ifile, make_1d_programm);
	getline(ifile, main_solve_file);
	getline(ifile, main_res_file);
	getline(ifile, name_data);
	ifile >> direction[0] >> direction[1] >> direction[2];
	ifile >> start_p[0] >> start_p[1] >> start_p[2];
	getline(ifile, name_file_centers); // ������ ������� ������
	getline(ifile, name_file_centers);
	ifile.close();

	//------------------------------------------------------------------------------------
	{
		std::string vtk_file_solve = main_solve_file + "0.vtk"; // �������� � ���� ��������� �� ����� ����� ����������
		std::string file_out = main_res_file + "0";

		char buf[512];

		// ������������ �������
		sprintf_s(buf, "%s %s %s %s %lf %lf %lf %lf %lf %lf", make_1d_programm.c_str(), vtk_file_solve.c_str(), file_out.c_str(), name_data.c_str(),
			direction[0], direction[1], direction[2], start_p[0], start_p[1], start_p[2]);

		system(buf); // ������ ������������

		//std::remove((file_out + "trace.txt").c_str());
		std::remove((file_out + "trace.vtk").c_str());
		std::remove((file_out + "value.txt").c_str());
	}
	//------------------------------------------------------------------------------------

	std::vector<int> trace;
	{
		std::string file_trace = main_res_file + "0trace.txt";

		if (Make1dFromTrace(file_trace, name_file_centers, main_solve_file, max_number_of_iter)) RETURN_ERR("Error make1dFrom trace\n");

		std::remove(file_trace.c_str());
	}
	return 0;
}

#endif // USE_VTK
#endif //UTILS