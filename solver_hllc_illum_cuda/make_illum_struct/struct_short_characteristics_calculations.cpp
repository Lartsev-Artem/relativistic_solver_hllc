#include "struct_short_characteristics_calculations.h"

#ifdef USE_VTK

int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell) {

	vtkIdList* idc;

	int x, y, z;
	for (size_t i = 0; i < 4; i++)
	{
		idc = neighbor_cell->GetFace(i)->GetPointIds();
		x = idc->GetId(0);
		y = idc->GetId(1);
		z = idc->GetId(2);

		if (a == x && b == y && c == z) return i;
		else if (a == x && b == z && c == y) return i;
		else if (a == y && b == x && c == z) return i;
		else if (a == y && b == z && c == x) return i;
		else if (a == z && b == x && c == y) return i;
		else if (a == z && b == y && c == x) return i;

	}
	return -2;
}
int WriteCellFaces(const std::string name_file_cells, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	FILE* f;
	f = fopen(name_file_cells.c_str(), "wb");

	//Type A[3];
	//Type B[3];
	//Type C[3];
	vtkPoints* points_face;
	std::vector<Type>pp(9);

	int n = unstructured_grid->GetNumberOfCells();
	fwrite_unlocked(&n, sizeof(int), 1, f);

	//std::vector<Vector3>centers_face(n * 4);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < 4; j++) {
			points_face = unstructured_grid->GetCell(i)->GetFace(j)->GetPoints();

			points_face->GetPoint(0, pp.data());
			points_face->GetPoint(1, pp.data() + 3);
			points_face->GetPoint(2, pp.data() + 6);

			//centers_face[i*4+j] = Vector3((pp[0]+pp[3]+pp[6])/3, (pp[1] + pp[4] + pp[7]) / 3, (pp[2] + pp[5] + pp[8]) / 3 );

			fwrite_unlocked(pp.data(), sizeof(Type), 9, f);

		}

	}

	fclose(f);

	/*f = fopen((name_file_cells + "face_center.bin").c_str(), "wb"); 
	n *= 4;
	fwrite_unlocked(&n, sizeof(int), 1, f);
	fwrite_unlocked(centers_face.data(), sizeof(Vector3), centers_face.size(), f);
	fclose(f);*/

	return 0;
}

int WriteVertex(const std::string name_file_vertex, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	FILE* f;
	f = fopen(name_file_vertex.c_str(), "wb");

	const int n = unstructured_grid->GetNumberOfCells();
	fwrite_unlocked(&n, sizeof(int), 1, f);
	Eigen::Matrix4d vertex_tetra;

	for (size_t i = 0; i < n; i++) {
		SetVertexMatrix(i, unstructured_grid, vertex_tetra);
		fwrite_unlocked(vertex_tetra.data(), sizeof(Eigen::Matrix4d), 1, f);
	}
	fclose(f);

	return 0;
}

int WriteIdPairs(const std::string name_file_pairs, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::vector<int> all_pairs_face;

	FindNeighborsPairFaceAndBoundaries(unstructured_grid, all_pairs_face);

	FILE* f;
	f = fopen(name_file_pairs.c_str(), "wb");
	if (!f) printf("id_neighbors not open\n");

	int n = all_pairs_face.size();
	fwrite_unlocked(&n, sizeof(int), 1, f);
	fwrite_unlocked(all_pairs_face.data(), sizeof(int), all_pairs_face.size(), f);

	fclose(f);
	return 0;
}
int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra) {

	// 4 вершины треугольника(по столбцам и единицы в нижний строке)

	vtkPoints* points = unstructured_grid->GetCell(number_cell)->GetPoints();
	for (size_t j = 0; j < 3; j++) {
		vertex_tetra(j, 2) = points->GetPoint(2)[j];
		vertex_tetra(j, 0) = points->GetPoint(0)[j];
		vertex_tetra(j, 1) = points->GetPoint(1)[j];
		vertex_tetra(j, 3) = points->GetPoint(3)[j];
	}

	for (size_t i = 0; i < 4; i++)
		vertex_tetra(3, i) = 1;

	return 0;
}

#endif

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face) {
	// 3 узла интерполяции
		{
			straight_face << 1. / 6, 1. / 6, 1,
				2. / 3, 1. / 6, 1,
				1. / 6, 2. / 3, 1;
		}

		// 3 узла интерполяции на наклонной плоскости
		{
			inclined_face <<
				0, sqrt(2. / 3), 1,
				sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
				-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
		}

		//Матрицы перехода из стандартного тетраэдра в координаты наклонной плоскости 
		{ transform_matrix <<
			-1. / sqrt(2), 1. / sqrt(2), 0,
			-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
			1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
		}

		//Матрицы перехода из наклонной плоскости в  координаты стандартного тетраэдра
		{
			inverse_transform_matrix <<
				-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				0, sqrt(2. / 3), 1. / sqrt(3);
		}

		// Начало координата плоскости
		start_point_plane_coord << 0.5, 0.5, 0;
		return 0;
}


int InitNodesValue(const std::vector<int>& all_pairs_face, std::vector<cell>& nodes_value) {


	const int n = all_pairs_face.size()/4;
	nodes_value.resize(n);

	for (size_t i = 0; i < n; ++i)
	{
		//nodes_value[i].id = i;	
		for (int j = 0; j < 4; j++) {
			nodes_value[i].neighbours_id_face[j] = all_pairs_face[i * 4 + j];
			nodes_value[i].nodes_value[j] = Vector3(-666, -666, -666);
		}
	}
	return 0;
}

int ResetNodesValue(std::vector<cell>& nodes_value) {


	const int n = nodes_value.size();
	
	for (size_t i = 0; i < n; ++i)
	{
		//nodes_value[i].id = i;	
		for (int j = 0; j < 4; j++) {
			nodes_value[i].nodes_value[j] = Vector3(-666, -666, -666);
		}
	}
	return 0;
}


Type NormIllum(const std::vector<Type>& Illum, const std::vector<Type>& Illum2) {
	Type max = -1;
	Type buf;
	for (size_t i = 0; i < Illum.size(); i++)
	{
		buf = fabs(Illum[i] - Illum2[i]);
		if (buf > max)
			max = buf;
	}
	return max;
}


int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state) {
	//face_state  -0=> выходящая грань,  1=> входящая  face_state.size=4!!!  

	Vector3 normal;

	for (size_t i = 0; i < 4; ++i) {

		normal = normals[number_cell].n[i];

		if (normal.dot(direction) < -eps)
			face_state[i] = 1;
		else
			face_state[i] = 0;  // если грань параллельна, считаем, что она не является определяющей

	}

	return 0;
}


int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord)
{
	// vertex_tetra -> [X;Y;Z;1]
	// возможно надо будет использовать transpose из-за инициализации матрицы перехода

	Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
	// сразу с транспонированием
	for (int i = 0; i < 3; ++i)
	{
		local_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			local_coord[i] += vertex_tetra_inverse(i, j) * global_coord[j];
		local_coord[i] += vertex_tetra_inverse(i, 3);
	}

	//local_coord = vertex_tetra * global_coord;
	return 0;
}
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord) {
	// vertex_tetra -> [X,Y,Z,1]
	// написать в ручную т.к. преобразования известны, и 4я строка постоянна и меняться не должна

	Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
	for (int i = 0; i < 3; ++i)
	{
		global_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			global_coord[i] += vertex_tetra(i, j) * local_coord[j];
		global_coord[i] += vertex_tetra(i, 3) * eta4;
	}

	//global_coord = vertex_tetra.inverse() * local_coord;

	return 0;
}

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}



Type BoundaryFunction(const int id_cell, const Vector3& x, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares) {
	//Type betta = 0.5;

	//Type I0 = betta * GetSPart(id_cell, direction, illum_old, directions, squares, square_surface);

	return 0;
}
Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
	return (3. * (1 + pow(direction.dot(direction2),2))) / 4;
	} };


	Vector3 cur_direction;
	Type S = 0;
	const int N_dir = directions.size();
	const int N_cell = illum_old.size() / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		S += Gamma(directions[num_direction], direction) * illum_old[num_direction * N_cell + num_cell] * squares[num_direction];
	}
	return S / square_surface;     // было *4PI, но из-за нормировки Gamma разделили на 4PI
}


Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})
	return interpolation_nodes.partialPivLu().solve(function_value);
}

int Min(const int a, const int b) {
	if (a < b) return a;
	return b;
}

size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i) {
		energy[i] = IntegarteDirection(i, Illum, squares, scuare_surface);
	}

	return 0;
}
Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();  // number of directions
	int m = Illum.size() / n;  // number of cells

	std::vector<pair<Type, Type>> I(n);


	for (size_t i = 0; i < n; i++) {
		I[i] = make_pair(Illum[m * i + num_cell], squares[i]);
	}

	auto cmp{ [](const pair<Type,Type> left, const pair<Type,Type> right) {

		return left.first < right.first;
} };
	std::sort(I.begin(), I.end(), cmp);

	for (size_t i = 0; i < n; i++) {
		res += I[i].first * I[i].second;
	}

	/*for (size_t i = 0; i < n; i++) {
		res += Illum[m * i + num_cell] * squares[i];
		}*/

	return res / scuare_surface;
}


int ReadCellFaces(const std::string name_file_cells, std::vector<Face>& grid) {

	FILE* f;
	f = fopen(name_file_cells.c_str(), "rb");

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	grid.resize(n * 4);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < 4; j++) {
			fread_unlocked(grid[i * 4 + j].A.data(), sizeof(Type), 3, f);
			fread_unlocked(grid[i * 4 + j].B.data(), sizeof(Type), 3, f);
			fread_unlocked(grid[i * 4 + j].C.data(), sizeof(Type), 3, f);
		}
	}

	fclose(f);

	return 0;
}

int WriteSize(const std::string& name_file_size) {

	std::ifstream ofile(name_file_size);
	uint64_t tr;
	ofile >> tr;
	ofile.close();

	std::ofstream ofile2(name_file_size);

	ofile2 << tr << '\n';
	ofile2 << posX << '\n';
	ofile2 << posX0 << '\n';
	ofile2 << posOutC << '\n';
	ofile2 << posOut << '\n';
	ofile2 << posIn << '\n';
	ofile2 << posS << '\n';
	ofile2 << posRes << '\n';
	/*
		X=3a
		x0
		OutC=N*M
		out=3a
		In=3a
		S=3a
		Res??
	*/

	ofile2.close();
}
int ReadVertex(const std::string name_file_vertex, std::vector<Eigen::Matrix4d>& vertexs) {

	FILE* f;
	f = fopen(name_file_vertex.c_str(), "rb");

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	vertexs.resize(n);
	Eigen::Matrix4d vertex_tetra;

	for (size_t i = 0; i < n; i++) {
		fread_unlocked(vertexs[i].data(), sizeof(Eigen::Matrix4d), 1, f);
	}
	fclose(f);

	return 0;
}


size_t IntersectionWithPlaneDisk(const Vector3& X0, const Vector3& n, Vector3& res) {

	//  ----------полный расчет. Т.к. диск задается постоянной плоскостью, параметры можно задатб явно--------------
	/*
	 {
	std::vector<Vector3> curface(3);		 // точки задающие плоскость диска
			curface[0][0] = 1;
			curface[0][1] = 0;
			curface[0][2] = 0;

			curface[1][0] = 0;//0;
			curface[1][1] = 0.9928768384869221;//
			curface[1][2] = 0.11914522061843064;//;

			curface[2][0] = 2;//;
			curface[2][1] = 0;
			curface[2][2] = 0;// ;   // Wolfram
			}

	* const std::vector<Vector3>& face,
	Vector3 A = face[0];
	Vector3 B = face[1];
	Vector3 C = face[2];

	Type a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	Type b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	Type c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	Type d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	Type t = -(a * X0[0] + b * X0[1] + c * X0[2] + d) / (a * n[0] + b * n[1] + c * n[2]);
	*/

	/*
	a= 0
	b= 0.1191452206184306
	c= -0.9928768384869221
	d= 0
	*/

	const Type b = 0.1191452206184306;
	const Type c = -0.9928768384869221;

	const Type t = -(b * X0[1] + c * X0[2]) / (b * n[1] + c * n[2]);

	res = (t * n + X0);

	/*for (size_t i = 0; i < 3; i++)
		res[i] = (n[i] * t + X0[i]);*/

	return 0;
}