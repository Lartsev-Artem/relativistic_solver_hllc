#include "get_grid_data.h"

#ifdef UTILS
#if defined USE_VTK
#include "../../file_module/reader_vtk.h"

int GetAverageArea(int argc, char* argv[])
{
	std::string name_file_grid = "";
	if (argc <= 1)
	{
		RETURN_ERR("Error input data!\n Input: path\\file.vtk\n");
	}
	else
	{
		name_file_grid = argv[1];
	}

	vtkSmartPointer<vtkUnstructuredGrid>  ugrid = vtkSmartPointer<vtkUnstructuredGrid> ::New();

	if (ReadFileVtk(name_file_grid, ugrid)) RETURN_ERR("File vtk not opened\n");

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
	std::cout << setprecision(16) << "Average size cell: " << h_mid / N / 3 << '\n';

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
		printf("Error input data!\n");
		printf("Input:  path\\file_grid.vtk, path\\output_file,\n");
		printf("Additional param:\n");
		printf("Input:  direction x y\n");
		printf("Input:  start point x y\n");
		RETURN_ERR("");
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
	if (ReadFileVtk(file_vtk, ugrid)) RETURN_ERR("Error read file vtk\n");
	
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

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, file_trace.c_str());

	for (auto el : id_trace)
	{
		ofile << el << '\n';
	}
	ofile.close();

	printf("Trace has %d cells\n", id_trace.size());

	return 0;
}

#endif //USE_VTK
#endif //UTILS
