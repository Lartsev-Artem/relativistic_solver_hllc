#include "solve_short_characteristics_calculations.h"

#ifdef USE_VTK

size_t ReadDataArray(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	int& size_grid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, is_print)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	size_grid = unstructured_grid->GetNumberOfCells();

	return 0;
}
size_t ReWriteDataArray(const size_t class_file_vtk, const std::string name_file_vtk, const std::string& main_dir,
	int& size_grid, const bool is_print/*=false*/) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	if (ReadFileVtk(0, name_file_vtk, unstructured_grid, foo, foo, foo, is_print)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}

	size_grid = unstructured_grid->GetNumberOfCells();

	if (unstructured_grid->GetCellData()->GetNumberOfArrays() == 0) { printf("grid hasn't data\n"); return 0; }

	std::vector<Type> v_data;
	int size;
	for (size_t i = 0; i < unstructured_grid->GetCellData()->GetNumberOfArrays(); i++)
	{
		std::string name_data(unstructured_grid->GetCellData()->GetArrayName(i));
		vtkDataArray* data = unstructured_grid->GetCellData()->GetScalars(name_data.c_str());
		size = data->GetSize()-1;
		v_data.resize(size);
		for (size_t i = 0; i < size; i++)
			v_data[i] = data->GetTuple1(i);

		std::unique_ptr<FILE, int(*)(FILE*)> file_data(fopen((main_dir + name_data + ".bin").c_str(), "wb"), fclose);

		if (!file_data) { printf("file_data is not opened for writing\n"); return 1; }

		fwrite_unlocked(&size, sizeof(int), 1, file_data.get());
		fwrite_unlocked(v_data.data(), sizeof(Type), data->GetSize(), file_data.get());

		fclose(file_data.get());

	}

	return 0;
}
size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	/*Чтение исходного файла и запись в vtkUnstructuredGrid*/

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();


	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	switch (class_file_vtk) {
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	case 5:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("pressure");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("velocity");  // GetScalars("radEnLooseRate");
		break;
	}

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}

int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face) {

	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	all_pairs_face.resize(N * 4);
	for (int i = 0; i < N * 4; i++)
		all_pairs_face[i] = -2;

	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b, id_c;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

		for (int num_face = 0; num_face < 4; ++num_face) {
			if (all_pairs_face[num_cell * 4 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);
			id_c = idp->GetId(2);

			/*Может быть проблема с указателями на списки!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));
				all_pairs_face[num_cell * 4 + num_face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = num_cell * 4 + num_face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[num_cell * 4 + num_face] = -1; // граничная ячейка
			else
				std::cout << "More than 1 neighbor????\n";
		}

	}

	return count_unique_face;
}
size_t NormalToFace(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, int number_face, Vector3& n) {

	vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(num_cell)->GetFace(number_face)->GetPointIds();

	Type P0[3], P1[3], P2[3];

	unstructured_grid->GetPoints()->GetPoint(idp->GetId(0), P0);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(1), P1);
	unstructured_grid->GetPoints()->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];


	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructured_grid->GetCell(num_cell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++) {
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructured_grid->GetCell(num_cell)->GetPoints()->GetPoint(idp2->GetId(id), P3);


	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}
bool InTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type AA[3];
	Type BB[3];
	Type CC[3];
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(0, AA);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(1, BB);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(2, CC);


	Vector3 A, B, C, X;
	{
		Type Xx[3] = { XX[0],XX[1],XX[2] };
		Vector3 n;
		Matrix3 basis;
		NormalToFace(num_cell, unstructuredgrid, number_face, n);
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, Xx, X);
	}

	// линейная алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}
size_t CenterOfTetra(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra) {

	auto MakeS{ [](Type* P0,Type* P1,Type* P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };

	Type P0[3], P1[3], P2[3], P3[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(number_cell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };


	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		point_in_tetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}
	return 0;
}
size_t FindAllCenterOfTetra(const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, std::vector<Vector3>& point_in_tetra) {

	const int n = unstructuredgrid->GetNumberOfCells();
	point_in_tetra.resize(n);

	for (size_t i = 0; i < n; ++i)
		CenterOfTetra(i, unstructuredgrid, point_in_tetra[i]);

	return 0;
}
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
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	vtkSmartPointer<vtkUnstructuredGrid>& u_grid) {

	int n = u_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> EnergyArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	for (size_t i = 0; i < n; i++) {
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		IllumArray->InsertNextTuple1(vector_illum[i]);  // по первому направлению
	}

	EnergyArray->SetName("energy");
	u_grid->GetCellData()->AddArray(EnergyArray);

	IllumArray->SetName("illum");
	u_grid->GetCellData()->AddArray(IllumArray);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(u_grid);
	writer->Write();
	return 0;
}
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls,
	vtkSmartPointer<vtkUnstructuredGrid>& u_grid) {

	int n = u_grid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> EnergyArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkDoubleArray> StreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	StreamArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> ImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	ImpulsArray->SetNumberOfComponents(9);


	for (size_t i = 0; i < n; i++) {
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		IllumArray->InsertNextTuple1(vector_illum[i]);  // по первому направлению		
		StreamArray->InsertNextTuple3(vector_stream[i](0), vector_stream[i](1), vector_stream[i](2));
		ImpulsArray->InsertNextTuple9(vector_impuls[i](0, 0), vector_impuls[i](0, 1), vector_impuls[i](0, 2),
			vector_impuls[i](1, 0), vector_impuls[i](1, 1), vector_impuls[i](1, 2),
			vector_impuls[i](2, 0), vector_impuls[i](2, 1), vector_impuls[i](2, 2));
	}

	EnergyArray->SetName("energy");
	u_grid->GetCellData()->AddArray(EnergyArray);

	IllumArray->SetName("illum");
	u_grid->GetCellData()->AddArray(IllumArray);


	StreamArray->SetName("stream");
	u_grid->GetCellData()->SetActiveVectors("stream");
	u_grid->GetCellData()->SetVectors(StreamArray);


	ImpulsArray->SetName("impuls");
	u_grid->GetCellData()->SetActiveTensors("impuls");
	u_grid->GetCellData()->SetTensors(ImpulsArray);

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(u_grid);
	writer->Write();
	return 0;
}

size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::string& file_vtk) {

	vtkSmartPointer<vtkUnstructuredGrid> ugrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	ReadFileVtk(0, file_vtk, ugrid, foo, foo, foo, false);

	WriteFileSolution(name_file_out, vector_illum, vector_energy, vector_stream, vector_impuls, ugrid);

	return 0;
}
size_t WriteFileSolution(const std::string name_file_out, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::string& file_vtk) {

	vtkSmartPointer<vtkUnstructuredGrid> ugrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* foo;
	ReadFileVtk(0, file_vtk, ugrid, foo, foo, foo, false);

	WriteFileSolution(name_file_out, vector_illum, vector_energy, ugrid);

	return 0;
}
size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//вершины треугольника

	Type A[3];
	Type B[3];
	Type C[3];
	face->GetPoints()->GetPoint(0, A);
	face->GetPoints()->GetPoint(1, B);
	face->GetPoints()->GetPoint(2, C);


	Type a, b, c, d;  // параметры уравнения плоскости
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}

size_t ReBuildDataArray(const size_t class_file_vtk, const std::string name_file_vtk, const std::string name_file_out,
	const std::string& main_dir, const bool is_print/*=false*/) {

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

	vtkSmartPointer<vtkDoubleArray> StreamArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	StreamArray->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> ImpulsArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	ImpulsArray->SetNumberOfComponents(9);

	std::vector<Type> vector_energy;
	std::vector<Type> vector_illum;
	std::vector<Type> vector_div_stream;
	std::vector<Vector3> vector_stream;
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

	std::unique_ptr<FILE, int(*)(FILE*)> file_stream(fopen((main_dir + "stream.bin").c_str(), "rb"), fclose);
	if (!file_stream) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_stream.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_stream.resize(size);
	fread_unlocked(vector_stream.data(), sizeof(Vector3), size, file_stream.get());
	fclose(file_stream.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_impuls(fopen((main_dir + "impuls.bin").c_str(), "rb"), fclose);
	if (!file_impuls) { printf("file_data is not opened for reading\n"); return 1; }
	fread_unlocked(&size, sizeof(int), 1, file_impuls.get());
	if (n != size) { printf("Error size data\n"); return 1; }
	vector_impuls.resize(size);
	fread_unlocked(vector_impuls.data(), sizeof(Matrix3), size, file_impuls.get());
	fclose(file_impuls.get());


	for (size_t i = 0; i < n; i++) {
		EnergyArray->InsertNextTuple1(vector_energy[i]);
		DivStreamArray->InsertNextTuple1(vector_div_stream[i]);
		IllumArray->InsertNextTuple1(vector_illum[i]);  // по первому направлению		
		StreamArray->InsertNextTuple3(vector_stream[i](0), vector_stream[i](1), vector_stream[i](2));
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


	StreamArray->SetName("stream");
	unstructured_grid->GetCellData()->SetActiveVectors("stream");
	unstructured_grid->GetCellData()->SetVectors(StreamArray);


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





#else

size_t ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	std::vector<Vector3>& velocity, std::vector<Type>& pressure,
	int& size_grid, const bool is_print/*=false*/) {

	if (class_file_vtk == 0) {
		size_grid = 51167;// 19888;// 20475;
		printf("No grid data\n");
	}

	else if (class_file_vtk == 1) {
		FILE* f;
		f = fopen((main_dir + "alpha.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		density.resize(size_grid);
		fread(density.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "alpha.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		absorp_coef.resize(size_grid);
		fread(absorp_coef.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "Q.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		rad_en_loose_rate.resize(size_grid);
		fread(rad_en_loose_rate.data(), sizeof(Type), size_grid, f);
		fclose(f);


		//	std::unique_ptr<FILE, int(*)(FILE*)> file_data(fopen((main_dir + "alpha.bin").c_str(), "rb"), fclose);
		//	if (!file_data) { printf("file_data is not opened for reading\n"); return 1; }
		//	fread_unlocked(&size_grid, sizeof(int), 1, file_data.get());
		//	density.resize(size_grid);
		//	fread_unlocked(density.data(), sizeof(Type), size_grid, file_data.get());
		//	fclose(file_data.get());   

		 //	ifstream ifile(main_dir + "alpha.bin", std::ios::binary);
			//ifile >> size_grid;
			//density.resize(size_grid);
			//for (size_t i = 0; i < size_grid; i++)		
				//ifile >> density[i];		
			//ifile.close();


	//		std::unique_ptr<FILE, int(*)(FILE*)> file_data1(fopen((main_dir + "alpha.bin").c_str(), "rb"), fclose);
		//	if (!file_data1) { printf("file_data is not opened for reading\n"); return 1; }
			//fread_unlocked(&size_grid, sizeof(int), 1, file_data1.get());
			//absorp_coef.resize(size_grid);
		  //for(int i=0; i<size_grid;i++)
			//fread_unlocked(absorp_coef.data()+i, sizeof(Type), 1, file_data1.get());
			//fclose(file_data1.get());

		//	ifstream ifile3(main_dir + "alpha.bin", std::ios::binary);
		//	ifile3 >> size_grid;
	//		absorp_coef.resize(size_grid);
	//		for (size_t i = 0; i < size_grid; i++)
	//		{
	//			ifile3 >> absorp_coef[i];
	//		}
	//		ifile3.close();

		//	std::unique_ptr<FILE, int(*)(FILE*)> file_data2(fopen((main_dir + "Q.bin").c_str(), "rb"), fclose);
		//	if (!file_data2) { printf("file_data is not opened for reading\n"); return 1; }
		//	fread_unlocked(&size_grid, sizeof(int), 1, file_data2.get());
		//	rad_en_loose_rate.resize(size_grid);
	 //     for(int i=0; i<size_grid;i++)
		//	fread_unlocked(rad_en_loose_rate.data()+i, sizeof(Type), 1, file_data2.get());
		//	fclose(file_data2.get());

		//	ifstream ifile2(main_dir + "Q.bin", std::ios::binary);
		//	ifile2 >> size_grid;
	//		rad_en_loose_rate.resize(size_grid);
		//	for (size_t i = 0; i < size_grid; i++)
	//		{
	//			ifile2 >> rad_en_loose_rate[i];
	//		}
	//		ifile2.close();

		if (rad_en_loose_rate.size() != density.size() ||
			rad_en_loose_rate.size() != absorp_coef.size() ||
			absorp_coef.size() != density.size()) {

			printf("Error size data\n");
			return 1;
		}

	}
	else if (class_file_vtk == 2)
	{
		FILE* f;
		f = fopen((main_dir + "density.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		density.resize(size_grid);
		fread(density.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "AbsorpCoef.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		absorp_coef.resize(size_grid);
		fread(absorp_coef.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "radEnLooseRate.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		rad_en_loose_rate.resize(size_grid);
		fread(rad_en_loose_rate.data(), sizeof(Type), size_grid, f);
		fclose(f);
		//printf("There isn't settings for fight grid data\n");
	}
	
	else if (class_file_vtk == 3) {
		FILE* f;
		f = fopen((main_dir + "density.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		density.resize(size_grid);
		fread(density.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "alpha.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		absorp_coef.resize(size_grid);
		fread(absorp_coef.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "Q.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		rad_en_loose_rate.resize(size_grid);
		fread(rad_en_loose_rate.data(), sizeof(Type), size_grid, f);
		fclose(f);

		f = fopen((main_dir + "velocity.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		velocity.resize(size_grid);
		fread(velocity.data(), sizeof(Vector3), size_grid, f);
		fclose(f);


		f = fopen((main_dir + "pressure.bin").c_str(), "rb");
		fread(&size_grid, sizeof(int), 1, f);
		pressure.resize(size_grid);
		fread(pressure.data(), sizeof(Type), size_grid, f);
		fclose(f);

	}


	if (is_print) {
		std::cout << "Grid has " << size_grid << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density.size() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef.size() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate.size() << '\n';
		}
	}
	return 0;
}

size_t ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	int& size_grid, const bool is_print/*=false*/) {




	if (class_file_vtk == 0) {
		size_grid = 19888; //20475;
		printf("No grid data\n");
	}

	else if (class_file_vtk == 1) {
 		FILE* f;
		f = fopen((main_dir + "alpha.bin").c_str(), "rb");		
		fread(&size_grid, sizeof(int), 1, f); 
		density.resize(size_grid);
		fread(density.data(), sizeof(Type), size_grid, f);
		fclose(f); 
   
   f = fopen((main_dir + "alpha.bin").c_str(), "rb");		
		fread(&size_grid, sizeof(int), 1, f);
		absorp_coef.resize(size_grid);
		fread(absorp_coef.data(), sizeof(Type), size_grid, f);
		fclose(f); 

   f = fopen((main_dir + "Q.bin").c_str(), "rb");		
		fread(&size_grid, sizeof(int), 1, f);
		rad_en_loose_rate.resize(size_grid);
		fread(rad_en_loose_rate.data(), sizeof(Type), size_grid, f);
		fclose(f); 


	//	std::unique_ptr<FILE, int(*)(FILE*)> file_data(fopen((main_dir + "alpha.bin").c_str(), "rb"), fclose);
	//	if (!file_data) { printf("file_data is not opened for reading\n"); return 1; }
	//	fread_unlocked(&size_grid, sizeof(int), 1, file_data.get());
	//	density.resize(size_grid);
	//	fread_unlocked(density.data(), sizeof(Type), size_grid, file_data.get());
	//	fclose(file_data.get());   
   	
     //	ifstream ifile(main_dir + "alpha.bin", std::ios::binary);
		//ifile >> size_grid;
		//density.resize(size_grid);
		//for (size_t i = 0; i < size_grid; i++)		
			//ifile >> density[i];		
		//ifile.close();
  

//		std::unique_ptr<FILE, int(*)(FILE*)> file_data1(fopen((main_dir + "alpha.bin").c_str(), "rb"), fclose);
	//	if (!file_data1) { printf("file_data is not opened for reading\n"); return 1; }
		//fread_unlocked(&size_grid, sizeof(int), 1, file_data1.get());
		//absorp_coef.resize(size_grid);
      //for(int i=0; i<size_grid;i++)
		//fread_unlocked(absorp_coef.data()+i, sizeof(Type), 1, file_data1.get());
		//fclose(file_data1.get());
   
	//	ifstream ifile3(main_dir + "alpha.bin", std::ios::binary);
	//	ifile3 >> size_grid;
//		absorp_coef.resize(size_grid);
//		for (size_t i = 0; i < size_grid; i++)
//		{
//			ifile3 >> absorp_coef[i];
//		}
//		ifile3.close();
 
	//	std::unique_ptr<FILE, int(*)(FILE*)> file_data2(fopen((main_dir + "Q.bin").c_str(), "rb"), fclose);
	//	if (!file_data2) { printf("file_data is not opened for reading\n"); return 1; }
	//	fread_unlocked(&size_grid, sizeof(int), 1, file_data2.get());
	//	rad_en_loose_rate.resize(size_grid);
 //     for(int i=0; i<size_grid;i++)
	//	fread_unlocked(rad_en_loose_rate.data()+i, sizeof(Type), 1, file_data2.get());
	//	fclose(file_data2.get());
 
 	//	ifstream ifile2(main_dir + "Q.bin", std::ios::binary);
	//	ifile2 >> size_grid;
//		rad_en_loose_rate.resize(size_grid);
	//	for (size_t i = 0; i < size_grid; i++)
//		{
//			ifile2 >> rad_en_loose_rate[i];
//		}
//		ifile2.close();

		if (rad_en_loose_rate.size() != density.size() ||
        rad_en_loose_rate.size() != absorp_coef.size() ||
        absorp_coef.size() != density.size()) {
  
			printf("Error size data\n");
			return 1;
		}

	}
	else if (class_file_vtk == 2)
		printf("There isn't settings for fight grid data\n");



	if (is_print) {
		std::cout << "Grid has " << size_grid << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density.size() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef.size() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate.size() << '\n';
		}
	}
	return 0;
}

size_t WriteFileSolution(const std::string& main_dir, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::vector<Type>& vector_div_stream) {
	int n;
 
	FILE* f;
	f = fopen((main_dir + "Illum.bin").c_str(), "wb");	
 	n = vector_illum.size();
    fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_illum.data(), sizeof(Type), n, f);	
	fclose(f);
 
	//std::unique_ptr<FILE, int(*)(FILE*)> file_Illum(fopen((main_dir + "Illum.bin").c_str(), "wb"), fclose);
//	if (!file_Illum) { printf("file_data is not opened for reading\n"); return 1; }
//	n = vector_illum.size();
//	fwrite_unlocked(&n, sizeof(int), 1, file_Illum.get());
//	fwrite_unlocked(vector_illum.data(), sizeof(Type), n, file_Illum.get());
//	fclose(file_Illum.get());

	f = fopen((main_dir + "energy.bin").c_str(), "wb");	
 	n = vector_energy.size();
  fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_energy.data(), sizeof(Type), n, f);	
	fclose(f);
 
//	std::unique_ptr<FILE, int(*)(FILE*)> file_energy(fopen((main_dir + "energy.bin").c_str(), "wb"), fclose);
//	if (!file_energy) { printf("file_data is not opened for reading\n"); return 1; }
//	n = vector_energy.size();
//	fwrite_unlocked(&n, sizeof(int), 1, file_energy.get());
//	fwrite_unlocked(vector_energy.data(), sizeof(Type), n, file_energy.get());
//	fclose(file_energy.get());

	f = fopen((main_dir + "stream.bin").c_str(), "wb");	
 	n = vector_stream.size();
  fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_stream.data(), sizeof(Vector3), n, f);	
	fclose(f);

//	std::unique_ptr<FILE, int(*)(FILE*)> file_stream(fopen((main_dir + "stream.bin").c_str(), "wb"), fclose);
//	if (!file_stream) { printf("file_data is not opened for reading\n"); return 1; }
//	n = vector_stream.size();
//	fwrite_unlocked(&n, sizeof(int), 1, file_stream.get());
//	fwrite_unlocked(vector_stream.data(), sizeof(Vector3), n, file_stream.get());
//	fclose(file_stream.get());

	f = fopen((main_dir + "impuls.bin").c_str(), "wb");	
 	n = vector_impuls.size();
  fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_impuls.data(), sizeof(Matrix3), n, f);	
	fclose(f);
 
//	std::unique_ptr<FILE, int(*)(FILE*)> file_impuls(fopen((main_dir + "impuls.bin").c_str(), "wb"), fclose);
//	if (!file_impuls) { printf("file_data is not opened for reading\n"); return 1; }
//	n = vector_impuls.size();
//	fwrite_unlocked(&n, sizeof(int), 1, file_impuls.get());
//	fwrite_unlocked(vector_impuls.data(), sizeof(Matrix3), n, file_impuls.get());
//	fclose(file_impuls.get());

	f = fopen((main_dir + "divstream.bin").c_str(), "wb");
	n = vector_div_stream.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_div_stream.data(), sizeof(Type), n, f);
	fclose(f);

	return 0;
}

size_t WriteFileSolution(const std::string& main_dir, const std::vector<Type>& vector_illum, const std::vector<Type>& vector_energy,
	const std::vector<Vector3>& vector_stream, const std::vector<Matrix3>& vector_impuls, const std::vector<Type>& vector_div_stream,
	const std::vector<Vector3>& vector_div_impuls, const std::vector<VectorX>& vector_U) {
	int n;

	FILE* f;
	f = fopen((main_dir + "Illum.bin").c_str(), "wb");
	n = vector_illum.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_illum.data(), sizeof(Type), n, f);
	fclose(f);


	f = fopen((main_dir + "energy.bin").c_str(), "wb");
	n = vector_energy.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_energy.data(), sizeof(Type), n, f);
	fclose(f);

	
	f = fopen((main_dir + "stream.bin").c_str(), "wb");
	n = vector_stream.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_stream.data(), sizeof(Vector3), n, f);
	fclose(f);

	
	f = fopen((main_dir + "impuls.bin").c_str(), "wb");
	n = vector_impuls.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_impuls.data(), sizeof(Matrix3), n, f);
	fclose(f);
	

	f = fopen((main_dir + "divstream.bin").c_str(), "wb");
	n = vector_div_stream.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_div_stream.data(), sizeof(Type), n, f);
	fclose(f);

	f = fopen((main_dir + "divimpuls.bin").c_str(), "wb");
	n = vector_div_impuls.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(vector_div_impuls.data(), sizeof(Vector3), n, f);
	fclose(f);

#ifndef RHLLC
	n = vector_U.size();
	std::vector<Type> density(n);
	std::vector<Type> pressure(n);
	std::vector<Vector3> velocity (n);

	for (size_t i = 0; i < n; i++)
	{
		const Type d = vector_U[i](0);
		density[i] = d;
		velocity[i](0) = vector_U[i](1) / d;
		velocity[i](1) = vector_U[i](2) / d;
		velocity[i](2) = vector_U[i](3) / d;
		const Type v = velocity[i].norm();
		pressure[i] = (vector_U[i](4) - v * v * d / 2.)* (gamma1 - 1);
	}
#endif

	f = fopen((main_dir + "density.bin").c_str(), "wb");
	n = density.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(density.data(), sizeof(Type), n, f);
	fclose(f);

	f = fopen((main_dir + "pressure.bin").c_str(), "wb");
	n = pressure.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(pressure.data(), sizeof(Type), n, f);
	fclose(f);

	f = fopen((main_dir + "velocity.bin").c_str(), "wb");
	n = velocity.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(velocity.data(), sizeof(Vector3), n, f);
	fclose(f);

	return 0;
}
#endif

size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Vector3>& directions_all, 
	vector<Type>& squares, Type& square_surface) {

	ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read(not open) file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	directions_all.resize(N);
	squares.resize(N);

	for (int i = 0; i < N; i++) {
		ifile >> squares[i];
		ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
	}
	ifile >> square_surface;
	ifile.close();

	return 0;
}

int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face, Matrix3& inclined_face_inverse, Matrix3& straight_face_inverse) {
	// 3 узла интерполяции
		{
			straight_face << 1. / 6, 1. / 6, 1,
				2. / 3, 1. / 6, 1,
				1. / 6, 2. / 3, 1;

			straight_face_inverse = straight_face.inverse(); // в решении
		}

		// 3 узла интерполяции на наклонной плоскости
		{
			inclined_face <<
				0, sqrt(2. / 3), 1,
				sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
				-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;

			inclined_face_inverse = inclined_face.inverse();  // в решении
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
			//nodes_value[i].neighbours_id_face[j] = all_pairs_face[i * 4 + j];
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

int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals) {

	FILE* f;
	f = fopen(name_file_normals.c_str(), "rb");

	int n;	
	fread_unlocked(&n, sizeof(int), 1, f);
	normals.resize(n);

	Normals norm(4);
	for (size_t i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++)
			fread_unlocked(norm.n[j].data(), sizeof(Type), 3, f);
		normals[i] = norm;
	}	
	fclose(f);




	//std::ifstream ifile;

	//ifile.open(name_file_normals);
	//if (!ifile.is_open()) {
	//	std::cout << "Error read file normals\n";
	//	return 1;
	//}

	//int N;
	//ifile >> N;
	//normals.resize(N);

	//Normals norm(4);
	//for (size_t i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//		ifile >> norm.n[j][0] >> norm.n[j][1] >> norm.n[j][2];
	//	normals[i] = norm;
	//}

	//ifile.close();
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
Type NormIllumOmp(const std::vector<Type>& Illum, const std::vector<Type>& Illum2) {
	Type max = -1;
	
	const int n = Illum.size();
	
#pragma omp parallel default(none) shared(n, max, Illum, Illum2)
	{
		Type buf;
#pragma omp for
		for (int i = 0; i < n; i++)
		{
			buf = fabs(Illum[i] - Illum2[i]);// / fabs(Illum[i]);

			if (buf > max)
			{
#pragma omp critical 
				{
					if (buf > max)
						max = buf;
				}
			}
		}
	}
	return max;
}

int ReadGraph(const std::string name_file_graph, std::vector<int>& sorted_id_cell) {
	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	int i;
	int count = 0;
	while (ifile >> i) {
		sorted_id_cell[count++] = i;
	}
	ifile.close();
	return 0;
}

int ReadGraphBin(const std::string name_file_graph, std::vector<int>& sorted_id_cell) {

	std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen(name_file_graph.c_str(), "rb"), fclose);
	if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

	const int n = sorted_id_cell.size();
	fread_unlocked(sorted_id_cell.data(), sizeof(int), n, file_graph.get());

	fclose(file_graph.get());
	return 0;
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

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord,
	Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord,
	Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}



size_t SetBasis(const Type* start_point, Vector3& normal, Matrix3& basis) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/
	Vector3 vec_1;
	Vector3 vec_2;

	if (abs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Matrix3& local_basis, const Type* point, Vector3& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
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
		Type I = (illum_old[num_direction * N_cell + num_cell]     + illum_old[num_direction * N_cell + num_cell + 1] +
			      illum_old[num_direction * N_cell + num_cell + 2] + illum_old[num_direction * N_cell + num_cell + 3]) / 4;

		S += Gamma(directions[num_direction], direction) * I * squares[num_direction];
	}
	return S / square_surface;     // было *4PI, но из-за нормировки Gamma разделили на 4PI
}

int CalculateInt(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering) {

	Vector3 direction;
	for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
		direction = directions[num_direction];
		for (size_t cell = 0; cell < num_cells; cell++)
		{
			int_scattering[num_direction * num_cells + cell] = GetS(4 * cell, direction, illum, directions, squares);;
		}
	}

	return 0;
}

int CalculateIntOmp(const int num_cells, const int num_directions, const std::vector<Type>& illum,
	const vector<Vector3>& directions, const vector<Type>& squares, vector<Type>& int_scattering) {
	
#pragma omp parallel default(none) shared(num_cells, num_directions, illum, directions, squares, int_scattering, square_surface)
		{
			Vector3 direction;

#pragma omp for
			for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
				direction = directions[num_direction];
				for (size_t cell = 0; cell < num_cells; cell++)
				{
					int_scattering[num_direction * num_cells + cell] = GetS(4*cell, direction, illum, directions, squares);
				}
			}
		}
	return 0;
}

Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})
	
	return interpolation_nodes.partialPivLu().solve(function_value);
	
}
Vector3 GetInterpolationCoefInverse(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})

	return interpolation_nodes * function_value;
}

int Min(const int a, const int b) {
	if (a < b) return a;
	return b;
}



Matrix3 IntegarteDirection9(const int num_cell, const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface) {
	Matrix3 res = Matrix3::Zero();
	int n = squares.size();  // number of directions
	int m = Illum.size() / n;  // number of cells

	for (size_t i = 0; i < 3; i++)
		for (size_t k = 0; k < 3; k++)
	
	for (size_t j = 0; j < n; j++) {
		res(i,k) += directions[j][i] * directions[j][k] * (Illum[m * j + num_cell] * squares[j]);
	}

	return res / scuare_surface;
}
int MakeImpuls(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface, vector<Matrix3>& impuls) {

	const int n = impuls.size();

	for (size_t i = 0; i < n; ++i) {
		impuls[i] = IntegarteDirection9(i, Illum, directions, squares, scuare_surface);
	}

	return 0;
}

Vector3 IntegarteDirection3(const int num_cell, const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface) {
	Vector3 res = Vector3::Zero();
	int n = squares.size();  // number of directions
	int m = Illum.size() / n;  // number of cells
	
	for (size_t i = 0; i < n; i++) {
		res += directions[i] * (Illum[m * i + num_cell] * squares[i]);
		}

	return res / scuare_surface;
}
int MakeStream(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface, vector<Vector3>& stream) {

	const int n = stream.size();

	for (size_t i = 0; i < n; ++i) {
		stream[i] = IntegarteDirection3(i, Illum, directions, squares, scuare_surface);
	}

	return 0;
}
int MakeDivStream(const vector<Vector3>& stream, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, std::vector<int>& neighbours_id_face,vector<Type>& div_stream) {
	//	std::vector<Vector3> Stream(size_grid * 4);
	//	
	//	MakeStream(Illum, directions, squares, square_surface, Stream);

	//	/*for (size_t i = 0; i < size_grid; i++) {
	//		for (size_t j = 0; j < 4; j++)
	//			div_stream[i] += normals[i].n[j].dot(Stream[i * 4 + j]) * squares_cell[i * 4 + j];
	//		div_stream[i] /= volume[i];
	//	}*/
	const int n = stream.size();

	for (size_t i = 0; i < n; i++) {
		div_stream[i] = 0;
		for (size_t j = 0; j < 4; j++) {

			Vector3 S_neib;
			const int neighbor_id_face = neighbours_id_face[i * 4 + j];

			if (neighbor_id_face >= 0)
				S_neib = (stream[neighbor_id_face / 4] + stream[i]) / 2;
			else
				S_neib = stream[i];

			div_stream[i] += ((normals[i].n[j].dot(S_neib)) * squares_cell[i * 4 + j]);
		}
		div_stream[i] /= volume[i];
	}

	return 0;


}
Type IntegarteDirection(const int num_cell, const Type* Illum, const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();  // number of directions
	int m = size_grid; // Illum.size() / n;  // number of cells

//	std::vector<pair<Type, Type>> I(n);
//
//
//	for (size_t i = 0; i < n; i++) {
//		I[i] = make_pair(Illum[m * i + num_cell], squares[i]);
//	}
//
//	auto cmp{ [](const pair<Type,Type> left, const pair<Type,Type> right) {
//
//		return left.first < right.first;
//} };
//	std::sort(I.begin(), I.end(), cmp);
//
//	for (size_t i = 0; i < n; i++) {
//		res += I[i].first * I[i].second;
//	}

	for (size_t i = 0; i < n; i++) {
		res += Illum[m * i + num_cell] * squares[i];
	}

	return res / scuare_surface;
}


size_t MakeEnergy(const Type* Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i)
	{
		energy[i] = IntegarteDirection(i, Illum, squares, scuare_surface);
	}

	return 0;
}

size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i) 
	{
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



int ReadSizes(const std::string& name_file_size, int& countX, int& countX0, int& countOutC,
				int& countOut, int& countIn, int& countS, int& countRes, int& countTry) {

	std::ifstream ifile(name_file_size+".txt");

	if (!ifile.is_open()) {
		printf("Error read files size\n");
		return 1;
	}

	ifile >> countTry;
	ifile >> countX;
	ifile >> countX0;
	ifile >> countOutC;
	ifile >> countOut;
	ifile >> countIn;
	ifile >> countS;
	ifile >> countRes;

	ifile.close();
	return 0;
}

size_t ReadCompactFastGridDataOptMemory(const int count_dir, const int N, const Str_Type& main_dir,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	Str_Type& name_file_shift_out, Str_Type& name_file_shift_res, Str_Type& name_file_shift_x0, Str_Type& name_file_shift_try,
	std::vector<cell>& grid, std::vector<int>& neighbours_id_face, std::vector<ShortId>& OutC,
	std::vector<ShortId>& Out, std::vector<ShortId>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<int>& sorted_id_cell,
	vector<int>& ShiftOut, vector<int>& ShiftRes, vector<int>& ShiftX0, vector<int>& ShiftTry) {

	ShiftOut.resize(count_dir); ShiftRes.resize(count_dir); ShiftX0.resize(count_dir); ShiftTry.resize(count_dir);
	FILE* f;
	f = fopen((name_file_shift_out + ".bin").c_str(), "rb");
	fread_unlocked(ShiftOut.data(), sizeof(int), ShiftOut.size(), f);
	fclose(f);

	f = fopen((name_file_shift_res + ".bin").c_str(), "rb");
	fread_unlocked(ShiftRes.data(), sizeof(int), ShiftRes.size(), f);
	fclose(f);

	f = fopen((name_file_shift_x0 + ".bin").c_str(), "rb");
	fread_unlocked(ShiftX0.data(), sizeof(int), ShiftX0.size(), f);
	fclose(f);

	f = fopen((name_file_shift_try + ".bin").c_str(), "rb");
	fread_unlocked(ShiftTry.data(), sizeof(int), ShiftTry.size(), f);
	fclose(f);

	printf("Read Shift files\n");


	f = fopen((name_file_graph + ".bin").c_str(), "rb");
	fread_unlocked(sorted_id_cell.data(), sizeof(int), sorted_id_cell.size(), f);
	fclose(f);

#ifdef WRITE_LOG
	ofstream ofile;
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "Read graph files\n";
	ofile.close();
#endif
	printf("Read graph files\n");

	grid.resize(N);
	int countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry;
	ReadSizes(name_file_sizes, countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry);


	//name_file_id_neighbors
	f = fopen((name_file_id_neighbors + ".bin").c_str(), "rb");
	int size;
	fread_unlocked(&size, sizeof(int), 1, f);
	neighbours_id_face.resize(size, -5);
	fread_unlocked(neighbours_id_face.data(), sizeof(int), size, f);
	fclose(f);

	printf("Read neighbors files\n");


	//res_inner_bound.resize(countRes);
	//id_try_surface.resize(countTry);
	OutC.resize(countOutC);
	Out.resize(countOut);
	In.resize(countIn);
	//S.resize(countS);
	//X.resize(countX);
	//X0.resize(countX0);

	/*f = fopen((name_file_res + ".bin").c_str(), "rb");
	fread_unlocked(res_inner_bound.data(), sizeof(Type), countRes, f);
	fclose(f);

	f = fopen((name_file_id_try + ".bin").c_str(), "rb");
	fread(id_try_surface.data(), sizeof(int), countTry, f);
	fclose(f);*/

	printf("Read try, res files\n");

	//f = fopen((name_file_s + ".bin").c_str(), "rb");
	//fread_unlocked(S.data(), sizeof(Type), S.size(), f);
	//fclose(f);

	printf("Read s files\n");

	f = fopen((name_file_count_out_faces + ".bin").c_str(), "rb");
	fread(OutC.data(), sizeof(ShortId), OutC.size(), f);
	fclose(f);


	f = fopen((name_file_out_faces + ".bin").c_str(), "rb");
	fread(Out.data(), sizeof(ShortId), Out.size(), f);
	fclose(f);


	f = fopen((name_file_in_faces + ".bin").c_str(), "rb");
	fread(In.data(), sizeof(ShortId), In.size(), f);
	fclose(f);

#ifdef WRITE_LOG
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
	ofile << "Read  files\n";
	ofile.close();
#endif

	return 0;
}


size_t ReadCompactFastGridData(const int count_dir, const int N, const Str_Type& file_logs,
	Str_Type& name_file_in_faces, Str_Type& name_file_out_faces, Str_Type& name_file_count_out_faces, Str_Type& name_file_local_x0, Str_Type& name_file_x,
	Str_Type& name_file_s, Str_Type& name_file_id_neighbors, Str_Type& name_file_centers,
	Str_Type& name_file_dist_try, Str_Type& name_file_id_try, Str_Type& name_file_res, Str_Type& name_file_sizes, Str_Type& name_file_graph,
	Str_Type& name_file_shift_out, Str_Type& name_file_shift_res, Str_Type& name_file_shift_x0, Str_Type& name_file_shift_try,
	std::vector<cell>& grid, std::vector<int>& neighbours_id_face, std::vector<ShortId>& OutC,
	std::vector<ShortId>& Out, std::vector<ShortId>& In, std::vector<Type>& S, std::vector<Vector3>& X, std::vector<Vector2>& X0,
	std::vector<Type>& res_inner_bound, std::vector<int>& id_try_surface, vector<int>& sorted_id_cell,
	vector<uint64_t>& ShiftOut, vector<uint64_t>& ShiftRes, vector<uint64_t>& ShiftX0, vector<int>& ShiftTry) {

#ifdef WRITE_LOG
	ofstream ofile;
#endif

	ShiftOut.resize(count_dir); ShiftRes.resize(count_dir); ShiftX0.resize(count_dir); ShiftTry.resize(count_dir);
 	FILE* f;
	f = fopen((name_file_shift_out + ".bin").c_str(), "rb");	
	fread_unlocked(ShiftOut.data(), sizeof(uint64_t), ShiftOut.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  ShiftOut file\n";
	ofile.close();
#endif

	f = fopen((name_file_shift_res + ".bin").c_str(), "rb");
	fread_unlocked(ShiftRes.data(), sizeof(uint64_t), ShiftRes.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  ShiftRes file\n";
	ofile.close();
#endif

	f = fopen((name_file_shift_x0 + ".bin").c_str(), "rb");
	fread_unlocked(ShiftX0.data(), sizeof(uint64_t), ShiftX0.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  ShiftX0 file\n";
	ofile.close();
#endif

#ifdef TRY_SHIFT
	f = fopen((name_file_shift_try + ".bin").c_str(), "rb");
	fread_unlocked(ShiftTry.data(), sizeof(int), ShiftTry.size(), f);
	fclose(f);
#endif

	printf("Read Shift files\n");

#if ONE_GRAPH
 	f = fopen((name_file_graph + ".bin").c_str(), "rb");
	fread_unlocked(sorted_id_cell.data(), sizeof(int), sorted_id_cell.size(), f);
	fclose(f);
#else

	for (size_t i = 0; i < count_dir; i++)
	{
		f = fopen((name_file_graph+std::to_string(i) + ".bin").c_str(), "rb");
		fread_unlocked(sorted_id_cell.data() + i*size_grid, sizeof(int), size_grid, f);
		fclose(f);
	}
#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  graph file\n";
	ofile.close();
#endif
#endif

	printf("Read graph files\n");

	grid.resize(N);
	int countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry;
	ReadSizes(name_file_sizes, countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  size file\n";
	ofile.close();
#endif

	//name_file_id_neighbors
  	f = fopen((name_file_id_neighbors + ".bin").c_str(), "rb");
	int size;
	fread_unlocked(&size, sizeof(int), 1, f);
	neighbours_id_face.resize(size, -5);
	fread_unlocked(neighbours_id_face.data(), sizeof(int), size, f);	
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  neighbors file\n";
	ofile.close();
#endif

	printf("Read neighbors files\n");
 

	res_inner_bound.resize(countRes);
	id_try_surface.resize(countTry);
	OutC.resize(countOutC);
	Out.resize(countOut);
	In.resize(countIn);
	S.resize(countS);
	X.resize(countX);
	X0.resize(countX0);
 
#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  resize file\n";
	ofile.close();
#endif

	f = fopen((name_file_res + ".bin").c_str(), "rb");
	fread_unlocked(res_inner_bound.data(), sizeof(Type), countRes, f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  res_inner_bound file\n";
	ofile.close();
#endif
 
#ifdef TRY_SHIFT
 	f = fopen((name_file_id_try + ".bin").c_str(), "rb");
	fread(id_try_surface.data(), sizeof(int), countTry, f);
	fclose(f);
#endif
	printf("Read try, res files\n");

 	f = fopen((name_file_s + ".bin").c_str(), "rb");
	fread_unlocked(S.data(), sizeof(Type), S.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  fileS file\n";
	ofile.close();
#endif
 
	printf("Read s files\n");

  	f = fopen((name_file_count_out_faces + ".bin").c_str(), "rb");
	fread(OutC.data(), sizeof(ShortId), OutC.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  count_out_faces file\n";
	ofile.close();
#endif


  f = fopen((name_file_out_faces + ".bin").c_str(), "rb");
  fread(Out.data(), sizeof(ShortId), Out.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  out_faces file\n";
	ofile.close();
#endif

  f = fopen((name_file_in_faces + ".bin").c_str(), "rb");
  fread(In.data(), sizeof(ShortId), In.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  in_faces file\n";
	ofile.close();
#endif

	printf("Read in out files\n");

  f = fopen((name_file_x + ".bin").c_str(), "rb");
  fread(X.data(), sizeof(Vector3), X.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  X file\n";
	ofile.close();
#endif

	printf("Read x files\n");

  
  f = fopen((name_file_local_x0 + ".bin").c_str(), "rb");
  if (!f) { printf("Err x0\n"); return 1; }
  fread(X0.data(), sizeof(Vector2), X0.size(), f);
	fclose(f);

#ifdef WRITE_LOG		
	ofile.open(file_logs, std::ios::app);
	ofile << "reading  X0 file\n";
	ofile.close();
#endif

	printf("Read x0 files\n");

	return 0;
}

int ReadCentersOfTetra(const std::string name_file_centers, std::vector<Vector3>& centers) {

	FILE* f;
	f = fopen(name_file_centers.c_str(), "rb");

	if (!f)
	{
		printf("File centers wasn't opened\n %s\n", name_file_centers.c_str());
		return 1;
	}

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);

	centers.resize(n);
	for (size_t i = 0; i < n; i++) {

		fread_unlocked(centers[i].data(), sizeof(Type), 3, f);
	}
	fclose(f);
	return 0;


	//std::unique_ptr<FILE, int(*)(FILE*)> file_centers(fopen(name_file_centers.c_str(), "rb"), fclose);
	//if (!file_centers) { printf("file_centers is not opened for readings\n"); return 1; }
	//
	//int n;
	//fread_unlocked(&n, sizeof(int), 1, file_centers.get());

	//centers.resize(n);
	//for (size_t i = 0; i < n; i++) {
	//	
	//	fread_unlocked(centers[i].data(), sizeof(Type), 3, file_centers.get());
	//}
	//fclose(file_centers.get());
}