#if defined UTILS
#include "rebuild_grid.h"
#include <algorithm>
#include <numeric> //iota

#include "../file_module/reader_bin.h"
#include "../file_module/reader_vtk.h"
#include "../file_module/reader_txt.h"

#include "../file_module/writer_vtk.h"
#include "../file_module/writer_bin.h"

static int ReadNetgenGrid(file_name name_file_in, std::vector<Vector3>& point, std::vector<Eigen::Vector4i>& cell, const int size = 3)
{
	std::ifstream ifile;
	OPEN_FSTREAM(ifile, name_file_in.c_str());

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
	OPEN_FSTREAM(ofile, name_file_out.c_str());

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

void GenerateMeshSizeFileNetgen()
{
	Type left = 0;
	Type right = 0.175;
	Type R1 = 0.05 / sqrt(2);
	Type R2 = 0.1;

	Type h = 0.0025;
	Type h_x = h;
	Type h_y = h;
	Type h_z = h;

	int Nx = (right - left) / h_x;
	int Ny = 2 * (R1 / h_y - 1) + 1;
	int Nz = 2 * (R1 / h_z - 1) + 1;

	std::ofstream out("D:\\Desktop\\FilesCourse\\Grids\\size.msz");

	out << Nx * Ny * Nz << '\n';
	int N = 0;
	Vector3 x(left + h_x / 2, -R1 + h_y, -R1 + h_z);
	for (int i = 0; i < Nx; i++)
	{
		x[1] = -R1 + h_y;
		for (int j = 0; j < Ny; j++)
		{
			x[2] = -R1 + h_z;
			for (int k = 0; k < Nz; k++)
			{
				out << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << h_x << '\n';
				x[2] += h_z;
				N++;
			}
			x[1] += h_y;
		}
		x[0] += h_x;
	}

	out << 0 << '\n';
	out.close();
	printf("N=%d\n", N);
	return;
}

int GetReNumberGraphMetis(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\metis_file.txt (simple_format)\n");
		printf("path\\graph.bin\n");		
		return 1;
	}

	const std::string metis_file = argv[1];
	const std::string graph_file = argv[2];

	std::vector<int> metis_nodes;
	std::vector<int> new_graph; // на i-ой позиции новый номер i-й ячейки

	if (ReadSimpleFileTxt(metis_file, metis_nodes)) RETURN_ERR("Metis file not open\n");

	const int N = metis_nodes.size();
	new_graph.resize(N, -1);

	int max = *std::max_element(metis_nodes.begin(), metis_nodes.end()) + 1;
	printf("Number of nodes: %d\n", max);
	
	int pos = 0;
	for (int id = 0; id < max; id++)
	{
		for (int i = 0; i < N; i++)
		{
			if (metis_nodes[i] == id)
			{
				new_graph[pos++] = i;
			}
		}
	}

	if (WriteSimpleFileBin(graph_file, new_graph)) RETURN_ERR("graph not written");
	if (WriteSimpleFileTxt(graph_file+".txt", new_graph)) RETURN_ERR("graph not written");
	return 0;
}

int GetReNumberGraphMetisToMPI(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\metis_file.txt (simple_format)\n");
		printf("path\\pairs.bin (after first rebuild!!!) \n");
		printf("path\\graph.bin\n");
		return 1;
	}

	const std::string metis_file = argv[1];
	const std::string pairs_file = argv[2];
	const std::string graph_file = argv[3];

	/*std::string metis_file = "D:\\Desktop\\FilesCourse\\Test\\metis\\cone_metis.txt.epart.2.txt";
	std::string pairs_file = "D:\\Desktop\\FilesCourse\\Test\\metis\\pairs.bin";
	std::string graph_file = "D:\\Desktop\\FilesCourse\\Test\\metis\\graph_mpi.bin";*/

	std::vector<int> metis_nodes;
	std::vector<int> pairs;	

	if (ReadSimpleFileTxt(metis_file, metis_nodes)) RETURN_ERR("Metis file not open\n");
	if (ReadSimpleFileBin(pairs_file, pairs)) RETURN_ERR("pairs file not open\n");

	const int N = metis_nodes.size();	

	int max = *std::max_element(metis_nodes.begin(), metis_nodes.end()) + 1;
	printf("Number of nodes: %d\n", max);

	int pos = 0;

	std::vector<int> size_claster(max, 0);

	for (int i = 0; i < N; i++)
	{
		size_claster[metis_nodes[i]]++;
	}

	
	std::vector<int> graph(N, 0);
	for (int i = 0; i < N; i++) graph[i] = i; // инициализация

	int left = 0;
	int right = 0;
	for (int id = 0; id < max; id++)
	{
		left = right;
		right += size_claster[id];
		
		int end_id = right - 1;
		int start_id = left;

		for (int i = left; i < right; i++)
		{
			for (int j = 0; j < base; j++)
			{
				int neig = pairs[graph[i] * base + j];
				if (neig / base > N) D_LD;

				if (neig >= 0)
				{
					int cell = neig / base;
					if (cell >= right)
					{
						if (end_id< 0) 
							D_LD;
						if (i < 0) D_LD;
						std::swap(graph[end_id--], graph[i--]);						
						break; // ячейка уже перестроена, пропускаем её
					}

					if (cell < left)
					{
						if (start_id >= N) D_LD;
						if (i < 0)
						{
							D_LD;
						}
						std::swap(graph[start_id++], graph[i--]);
						break; // ячейка уже перестроена, пропускаем её
					}
				}
			}
		}		
	}	

	if (WriteSimpleFileBin(graph_file, graph)) RETURN_ERR("graph not written");	

	printf("SUCCESS\n");
	return 0;
}

int RenumberNodesMetis(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\metis_file.txt (simple_format)\n");
		printf("path\\centers.bin \n");		
		return 1;
	}
	std::string  file_metis = argv[1];
	std::string  file_centers = argv[2];

	/*std::string  file_metis = "D:\\Desktop\\FilesCourse\\ConeCyl\\conecyl_m_epart_3.txt";
	std::string  file_centers = "D:\\Desktop\\FilesCourse\\ConeCyl\\centers.bin";*/

	std::vector<int> metis_id; // i-cells, a[i] - id
	std::vector<Vector3> centers;

	if (ReadSimpleFileTxt(file_metis, metis_id)) RETURN_ERR("error reading metis\n");
	if (ReadSimpleFileBin(file_centers, centers)) RETURN_ERR("error reading centers\n");
	
	int max = *std::max_element(metis_id.begin(), metis_id.end()) + 1;
	printf("Number of nodes: %d\n", max);

	std::vector<std::pair<Type, Type>> pos(max, std::make_pair(metis_id.size() + 1, -1));
	Vector3 base_point(0, 0, 0);
	
	for (int i = 0; i < metis_id.size(); i++)
	{
		int id = metis_id[i];
		Type r = centers[i][0] - base_point[0];

		pos[id].first = std::min(r, pos[id].first);
		pos[id].second = std::max(r, pos[id].second);
	}
	centers.clear();

	std::vector<std::pair<int, Type>> node_center_point(max);
	std::pair<int, Type> n(-1, 0);
	std::generate(node_center_point.begin(), node_center_point.end(), 
		[&n, &pos] { n.first++; n.second= (pos[n.first].first + pos[n.first].second)/2; return n; });

	std::sort(node_center_point.begin(), node_center_point.end(),
		[](std::pair<int, Type> l, std::pair<int, Type> r) { return l.second < r.second; });

	
	std::vector<int> metis_id_new(metis_id.size(), -1);
	for (int i = 0; i < metis_id.size(); i++)
	{
		metis_id_new[i] = node_center_point[metis_id[i]].first;		
	}
	
	WriteSimpleFileTxt(file_metis, metis_id_new);

	for (int i = 0; i < max; i++)
	{
		printf("renum nodes %d -> %d\n",i,node_center_point[i].first);
	}

	return 0;
}


// узлы должны иметь одного соседа. При соседстве  с разными узлами результат не гарантируется
int GetMpiConfToHLLC(int argc, char* argv[])
{
	/*
	Пока делаем так
	for(0,N){if face in node -> do, else continue}

	wait mpi rcv

	for(0,N){if cell in node, face not -> do}, else continue

	Да, грани гоняем дважды и по всему диапазону. Но зато это не требует пересортировки граней, а ячейки уже отсортированы

	0 и N можно сократить до min_max on node

	*/
	if (argc != 4)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\file_geometry_faces\n");
		printf("path\\metis.txt file\n");
		printf("path\\out_file.txt\n");

		printf("path\\FORMAT:\nnp \nid_i \nleft_i \nright_i \nleft_min_i \nleft_max_i \nright_min_i \nright_max_i\n");
	//	return 1;
	}

	//std::string file_geometry_faces = argv[1];
	//std::string  file_metis = argv[2]; // "D:\\Desktop\\FilesCourse\\Test\\metis\\cone_metis.txt.epart.2.txt";
	//std::string	file_output = argv[3];

	std::string file_geometry_faces = "D:\\Desktop\\FilesCourse\\ConeCyl\\geo_faces.bin";
	std::string  file_metis = "D:\\Desktop\\FilesCourse\\ConeCyl\\conecyl_m_epart_3.txt";
	std::string	file_output = "D:\\Desktop\\FilesCourse\\ConeCyl\\mpi_conf.txt";

	grid_t grid;
	READ_FILE(file_geometry_faces.c_str(), grid.faces, geo);

	std::vector<int> metis_id; // i-cells, a[i] - id
	if (ReadSimpleFileTxt(file_metis, metis_id)) RETURN_ERR("error reading metis\n");

	int np = *std::max_element(metis_id.begin(), metis_id.end()) + 1;

	std::vector<int> size_claster(np, 0);
	for (int i = 0; i < metis_id.size(); i++)
	{
		size_claster[metis_id[i]]++;
	}

	metis_id.assign(metis_id.size(), -1);
	int k = 0;
	for (int id = 0; id < np; id++)
	{
		for (int j = 0; j < size_claster[id]; j++)
		{
			metis_id[k++] = id; // new id-cell chain
		}
	}
#ifdef DEBUG
	{
		auto res = std::min_element(metis_id.begin(), metis_id.end());
		if (*res < 0)
		{
			RETURN_ERR("Bad new config\n");
		}
	}
#endif
//	metis_id.clear();

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, file_output.c_str());

	ofile << np << '\n';

	int disp = 0;
	for (int id = 0; id < np; id++)
	{
		std::vector<int> neig_cells_left;
		std::vector<int> neig_cells_right;
		int neig_r = -1;
		int neig_l = -1;

		int left = disp;
		int right = left + size_claster[id];

		for (auto& f : grid.faces)
		{
			if (f.geo.id_l >= left && f.geo.id_l < right) //текущий узел
			{
				if ((f.geo.id_r < left || f.geo.id_r >= right) && f.geo.id_r >= 0) //сосед на другом узле
				{
					if(neig_r==-1) neig_r = metis_id[f.geo.id_r];
					else if (neig_r != metis_id[f.geo.id_r])
					{
						RETURN_ERR("bad config neig right\n");
					}
											
					if (metis_id[f.geo.id_r] != id + 1)
					{
					//	RETURN_ERR("bad config neig right\n");
					}
					neig_cells_right.push_back(f.geo.id_r); //id соседей в глобальной нумерации					
				}
			}

			if (f.geo.id_r >= left && f.geo.id_r < right) //текущий узел
			{
				if ((f.geo.id_l < left || f.geo.id_l >= right) && f.geo.id_l >= 0) //сосед на другом узле
				{
					int buf = metis_id[f.geo.id_l];
					if (neig_l == -1) neig_l = metis_id[f.geo.id_l];
					else if (neig_l != metis_id[f.geo.id_l])
					{
						RETURN_ERR("bad config neig left\n");
					}

					
					if (metis_id[f.geo.id_l] != id - 1)
					{
					//	RETURN_ERR("bad config neig left\n");
					}
					neig_cells_left.push_back(f.geo.id_l); //id соседей в глобальной нумерации					
				}
			}
		}

		std::sort(neig_cells_left.begin(), neig_cells_left.end());
		std::sort(neig_cells_right.begin(), neig_cells_right.end());

		auto last = std::unique(neig_cells_left.begin(), neig_cells_left.end());
		neig_cells_left.erase(last, neig_cells_left.end());

		last = std::unique(neig_cells_right.begin(), neig_cells_right.end());
		neig_cells_right.erase(last, neig_cells_right.end());

		ofile << id << '\n';

		ofile << left << '\n';
		ofile << right << '\n';

		{
			if (neig_cells_left.size())
			{
				ofile << *std::min_element(neig_cells_left.begin(), neig_cells_left.end()) << '\n';
				ofile << *std::max_element(neig_cells_left.begin(), neig_cells_left.end()) << '\n';
			}
			else
			{
				ofile << -1 << '\n' << -1 << '\n';
			}

			if (neig_cells_right.size())
			{
				ofile << *std::min_element(neig_cells_right.begin(), neig_cells_right.end()) << '\n';
				ofile << *std::max_element(neig_cells_right.begin(), neig_cells_right.end()) << '\n';
			}
			else
			{
				ofile << -1 << '\n' << -1 << '\n';
			}
		}

		//ofile << neig_cells_left.size() << '\n';
		//ofile << neig_cells_right.size() << '\n';

		disp += size_claster[id];
	}

	ofile.close();
	return 0;
}

int GetMpiConfToHLLC_debug_manual(int argc, char* argv[])
{
	/*
	Пока делаем так
	for(0,N){if face in node -> do, else continue}

	wait mpi rcv

	for(0,N){if cell in node, face not -> do}, else continue

	Да, грани гоняем дважды и по всему диапазону. Но зато это не требует пересортировки граней, а ячейки уже отсортированы

	0 и N можно сократить до min_max on node

	*/
	if (argc != 4)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\file_geometry_faces\n");
		printf("path\\metis.txt file\n");
		printf("path\\out_file.txt\n");

		printf("path\\FORMAT:\nnp \nid_i \nleft_i \nright_i \nleft_min_i \nleft_max_i \nright_min_i \nright_max_i\n");
		//	return 1;
	}

	//std::string file_geometry_faces = argv[1];
	//std::string  file_metis = argv[2]; // "D:\\Desktop\\FilesCourse\\Test\\metis\\cone_metis.txt.epart.2.txt";
	//std::string	file_output = argv[3];

	std::string file_geometry_faces = "D:\\Desktop\\FilesCourse\\ConeCyl\\geo_faces.bin";
	std::string  file_metis = "D:\\Desktop\\FilesCourse\\ConeCyl\\conecyl_m_epart_3.txt";
	std::string	file_output = "D:\\Desktop\\FilesCourse\\ConeCyl\\mpi_conf.txt";

	grid_t grid;
	READ_FILE(file_geometry_faces.c_str(), grid.faces, geo);

	std::vector<int> metis_id; // i-cells, a[i] - id
	if (ReadSimpleFileTxt(file_metis, metis_id)) RETURN_ERR("error reading metis\n");

	int np = *std::max_element(metis_id.begin(), metis_id.end()) + 1;

	std::vector<int> size_claster(np, 0);
	for (int i = 0; i < metis_id.size(); i++)
	{
		size_claster[metis_id[i]]++;
	}

	metis_id.assign(metis_id.size(), -1);
	int k = 0;
	for (int id = 0; id < np; id++)
	{
		for (int j = 0; j < size_claster[id]; j++)
		{
			metis_id[k++] = id; // new id-cell chain
		}
	}
#ifdef DEBUG
	{
		auto res = std::min_element(metis_id.begin(), metis_id.end());
		if (*res < 0)
		{
			RETURN_ERR("Bad new config\n");
		}
	}
#endif
	//	metis_id.clear();

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, file_output.c_str());

	ofile << np << '\n';

	int disp = 0;
	for (int id = 0; id < np; id++)
	{
		/*std::vector<int> neig_cells_left;
		std::vector<int> neig_cells_right;*/
		std::vector < std::vector<int>> neig_cells_left(np);
		std::vector < std::vector<int>> neig_cells_right(np);
		int neig_r = -1;
		int neig_l = -1;

		int left = disp;
		int right = left + size_claster[id];

		for (auto& f : grid.faces)
		{
			if (f.geo.id_l >= left && f.geo.id_l < right) //текущий узел
			{
				if ((f.geo.id_r < left || f.geo.id_r >= right) && f.geo.id_r >= 0) //сосед на другом узле
				{
					if (neig_r == -1) neig_r = metis_id[f.geo.id_r];
					else if (neig_r != metis_id[f.geo.id_r])
					{
						//RETURN_ERR("bad config neig right\n");
					}

					if (metis_id[f.geo.id_r] != id + 1)
					{
						//	RETURN_ERR("bad config neig right\n");
					}
					neig_cells_right[metis_id[f.geo.id_r]].push_back(f.geo.id_r); //id соседей в глобальной нумерации					
				}
			}

			if (f.geo.id_r >= left && f.geo.id_r < right) //текущий узел
			{
				if ((f.geo.id_l < left || f.geo.id_l >= right) && f.geo.id_l >= 0) //сосед на другом узле
				{
					int buf = metis_id[f.geo.id_l];
					if (neig_l == -1) neig_l = metis_id[f.geo.id_l];
					else if (neig_l != metis_id[f.geo.id_l])
					{
						//RETURN_ERR("bad config neig left\n");
					}


					if (metis_id[f.geo.id_l] != id - 1)
					{
						//	RETURN_ERR("bad config neig left\n");
					}
					neig_cells_left[metis_id[f.geo.id_l]].push_back(f.geo.id_l); //id соседей в глобальной нумерации					
				}
			}
		}

		std::sort(neig_cells_left.begin(), neig_cells_left.end());
		std::sort(neig_cells_right.begin(), neig_cells_right.end());

		auto last = std::unique(neig_cells_left.begin(), neig_cells_left.end());
		neig_cells_left.erase(last, neig_cells_left.end());

		last = std::unique(neig_cells_right.begin(), neig_cells_right.end());
		neig_cells_right.erase(last, neig_cells_right.end());

		ofile << id << '\n';

		ofile << left << '\n';
		ofile << right << '\n';

		{
			if (neig_cells_left[id].size())
			{
				ofile << *std::min_element(neig_cells_left[id].begin(), neig_cells_left[id].end()) << '\n';
				ofile << *std::max_element(neig_cells_left[id].begin(), neig_cells_left[id].end()) << '\n';
			}
			else
			{
				ofile << -1 << '\n' << -1 << '\n';
			}

			if (neig_cells_right[id].size())
			{
				ofile << *std::min_element(neig_cells_right[id].begin(), neig_cells_right[id].end()) << '\n';
				ofile << *std::max_element(neig_cells_right[id].begin(), neig_cells_right[id].end()) << '\n';
			}
			else
			{
				ofile << -1 << '\n' << -1 << '\n';
			}
		}

		//ofile << neig_cells_left.size() << '\n';
		//ofile << neig_cells_right.size() << '\n';

		disp += size_claster[id];
	}

	ofile.close();
	return 0;
}

#include "grid_geometry/geometry_data.h"
int WriteGeoFiles(int argc, char* argv[])
{
	if (argc != 2)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\base_adress\n");		
		return 1;
	}
	glb_files.base_adress = argv[1];
	std::string name_file_geometry_faces = glb_files.base_adress + F_GEO_FACES;
	std::string name_file_geometry_cells = glb_files.base_adress + F_GEO_CELLS;

	ReWriteGeoFiles(name_file_geometry_faces, name_file_geometry_cells);

	return 0;
}

#ifdef USE_VTK
static int WriteDataToGrid(file_name name_file_grid, file_name name_file_data, file_name name_file_output, file_name name_data) 
{
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	
	if (ReadFileVtk(name_file_grid, unstructured_grid))  RETURN_ERR("Error reading the file vtk\n");
	
	const int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> DataArray = vtkSmartPointer<vtkDoubleArray>::New();

	std::vector<Type> vector_data;
	//std::vector<int> vector_data;

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
	if (argc < 4)
	{
		printf("Error input data!\n");
		printf("Input Format: path\\file.vtk,  path\\data.bin,  path\\outfile.vtk\n");
		printf("Additional parameters: name_field_data");
		printf("Data format: size a b c ... \n");
		return 1;
	}

	std::string file_grid = argv[1];
	std::string file_data = argv[2];
	std::string file_out = argv[3];
	std::string name_data = "data";
	
	if (argc > 4)
	{
		name_data = argv[4];
	}
	
	if (WriteDataToGrid(file_grid, file_data, file_out, name_data))
	{
		RETURN_ERR("Error write data to vtk grid\n");		
	}

	return 0;
}

#endif//USE_VTK

#endif //UTILS
