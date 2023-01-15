#include "Header.h"

int ReadCentersOfTetra(const std::string name_file_centers, std::vector<Vector3>& centers) {

	FILE* f;
	f = fopen(name_file_centers.c_str(), "rb");

	if (!f)
	{
		printf("File centers wasn't opened\n %s\n", name_file_centers.c_str());
		return 1;
	}

	int n;
	fread(&n, sizeof(int), 1, f);

	centers.resize(n);
	for (size_t i = 0; i < n; i++) {

		fread(centers[i].data(), sizeof(Type), 3, f);
	}
	fclose(f);
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
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input Format: path\\trace.txt, path\\centers.bin, path\\Solve_i, number_of_iter \n");		
		return 1;
#else
		 file_trace = "D:\\Desktop\\FilesCourse\\Plane2D\\trace0.txt";
		 file_centers = "D:\\Desktop\\FilesCourse\\Plane2D\\centers.bin";
		 file_solve = "D:\\Desktop\\FilesCourse\\Plane2D\\Solve\\Solve";
		 max_iter = 14;
#endif
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

int Make1dFromTrace(const std::string& file_trace, const std::string& file_centers, const std::string& file_solve,
	const int max_iter)
{
	/*std::string file_trace = "";
	std::string file_centers = "";
	std::string file_solve = "";
	int max_iter = 1;*/

	std::vector<int> trace;
	{
		//std::string file_trace = file_trace;  //trace.txt

		std::ifstream ifile(file_trace);
		if (!ifile.is_open()) ERR_RETURN("File file_trace(i) wasn't opened\n");

		std::string str;

		int a;
		while (!ifile.eof()) {
			ifile >> a;
			trace.push_back(a);
		}
		ifile.close();

		printf("trace has %d elem\n", trace.size());

		//std::remove((file_trace).c_str());
		std::vector<Vector3> centers;
		ReadCentersOfTetra(file_centers, centers);

		vtkDataArray* density;
		vtkDataArray* pressure;
		vtkDataArray* velocity;

		vtkSmartPointer<vtkUnstructuredGrid> u_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();


		std::ofstream ofileP(file_solve + "P.txt");
		std::ofstream ofileD(file_solve + "D.txt");
		std::ofstream ofileV(file_solve + "V.txt");
		if (!ofileP.is_open()) ERR_RETURN("File result wasn't opened\n");

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

int BuildHLLC_1dTime(int argc, char* argv[])
{
	std::string file_start = "";
	int max_number_of_iter = 1;

	if (argc < 2)
	{
#ifdef RELEASE_VER
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
		return 1;
#else
		file_start = "D:\\Desktop\\FilesCourse\\settings_file_build1d_hllc.txt";
#endif
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
	Vector3 direction (1, 0, 0);
	Vector3 start_p(0, 0, 0);
	std::string name_data;// = "density";	
	std::string  name_file_centers;// = "D:\\Desktop\\FilesCourse\\Soda\\centers.bin";


	std::ifstream ifile(file_start);
	if (!ifile.is_open()) ERR_RETURN("File settings wasn't opened\n");

	getline(ifile, make_1d_programm);
	getline(ifile, main_solve_file);
	getline(ifile, main_res_file);
	getline(ifile, name_data);
	ifile >> direction[0] >> direction[1] >> direction[2];
	ifile >> start_p[0] >> start_p[1] >> start_p[2];
	getline(ifile, name_file_centers); // просто перевод строки
	getline(ifile, name_file_centers);
	ifile.close();

	for (size_t i = 0; i < 1/*max_number_of_iter*/; i++)  // проекция у всех состояний на одной сетке одинаковая
	{
		std::string vtk_file_solve = main_solve_file + std::to_string(i) + ".vtk";
		std::string file_out = main_res_file + std::to_string(i);

		char buf[512];

		// формирование команды
		sprintf_s(buf, "%s %s %s %s %lf %lf %lf %lf %lf %lf", make_1d_programm.c_str(), vtk_file_solve.c_str(), file_out.c_str(), name_data.c_str(),
			direction[0], direction[1], direction[2], start_p[0], start_p[1], start_p[2]);		

		system(buf); // запуск подпрограммы

		//std::remove((file_out + "trace.txt").c_str());
		std::remove((file_out + "trace.vtk").c_str());
		std::remove((file_out + "value.txt").c_str());			
	}
	//------------------------------------------------------------------------------------
	
	std::vector<int> trace;
	{
		std::string file_trace = main_res_file + std::to_string(0) + "trace.txt";

		if (Make1dFromTrace(file_trace, name_file_centers, main_solve_file, max_number_of_iter))
		{
			printf("Error make1dFrom trace\n");
			return 1;
		}

		for (size_t i = 0; i < max_number_of_iter; i++)
		{
			std::remove((main_res_file + std::to_string(i) + "trace.txt").c_str());
		}
		
		return 0;

		/*std::ifstream ifile(file_trace);
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
			for (auto id : trace)
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
		ofileV.close();*/
	}

	return 0;
}

int Trace2D(int argc, char* argv[])
{
	std::string file_vtk = "";
	std::string file_trace = "";

	Eigen::Vector2d ray(1, 0);
	Eigen::Vector2d start_p(-1, 0.25);
	
	if (argc <= 2)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input:  path\\file_grid.vtk, path\\output_file,\n");
		printf("Additional param:\n");
		printf("Input:  direction x y\n");
		printf("Input:  start point x y\n");
		return 1;
#else
		file_vtk = "D:\\Desktop\\FilesCourse\\Plane2D\\Solve\\Solve0.vtk";
		file_trace = "D:\\Desktop\\FilesCourse\\Plane2D\\trace0.txt";
#endif
	}
	else
	{
		if (argc >= 3)
		{
			file_vtk = argv[1];
			file_trace = argv[2];
		}

		if (argc >= 5)
		{
			ray(0) = std::stod(argv[3]);
			ray(1) = std::stod(argv[4]);
		
			ray.normalize();
		}
		if (argc == 7)
		{

			start_p(0) = std::stod(argv[5]);
			start_p(1) = std::stod(argv[6]);
		}
	}

	printf("You input:\n");
	printf("file_grid.vtk:  %s\n", file_vtk.c_str());
	printf("file_trace.txt:  %s\n", file_trace.c_str());
	printf("direction x= %lf y= %lf\n", ray(0), ray(1));
	printf("start point x= %lf y= %lf\n", start_p(0), start_p(1));

	vtkSmartPointer <vtkUnstructuredGrid> ugrid = vtkSmartPointer < vtkUnstructuredGrid>::New();
	vtkDataArray* foo;
	if (ReadFileVtk(0, file_vtk, ugrid, foo, foo, foo, true))
	{
		ERR_RETURN("Error read file vtk\n");
	}

	const int N = ugrid->GetNumberOfCells();

	std::vector<int> id_trace; 
	id_trace.reserve(N);

	
	for (size_t i = 0; i < N; i++)
	{
		vtkCell* cell = ugrid->GetCell(i);

		Type a[3];
		Type b[3];
		Type c[3];
			
		cell->GetPoints()->GetPoint(0, a);
		cell->GetPoints()->GetPoint(1, b);
		cell->GetPoints()->GetPoint(2, c);
		
		auto Line{ [&ray, &start_p](Type x) {
			return (ray(1) * (x - start_p(0))) / ray(0) + start_p(1); // y из уравнения прямой
		} };
		
		bool flag_less = false;
		bool flag_more = false;

		
		if (a[1] < Line(a[0]) || b[1] < Line(b[0]) || c[1] < Line(c[0])) // точка ниже прямой
		{
			flag_less = true;
		}

		if (a[1] > Line(a[0]) || b[1] > Line(b[0]) || c[1] > Line(c[0])) // точка выше прямой
		{
			flag_more = true;
		}

		if (flag_more && flag_less) // есть точки по разные стороны от прямой
		{
			id_trace.push_back(i);
		}
	}

	std::ofstream ofile(file_trace);

	if (!ofile.is_open())
	{
		printf("File trace not opend\n");
		return 1;
	}
	
	for (auto el: id_trace)
	{
		ofile << el << '\n';
	}
	ofile.close();

	printf("Trace has %d cells\n", id_trace.size());
	
	return 0;
}