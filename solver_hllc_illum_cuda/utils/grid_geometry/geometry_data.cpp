#include "../../prj_config.h"

#include "../../global_headers.h"
#include "../../global_def.h"
#include "../../global_value.h"

#ifdef USE_VTK
int GetNeighborFace3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face) {

	auto GetNumberNeighborFace{ [](const int a, const int b, const int c, vtkCell* neighbor_cell)
		{

			vtkIdList* idc;

			int x, y, z;
			for (int i = 0; i < 4; i++)
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
		} };

	int count_unique_face = 0;
	
	const int N = unstructured_grid->GetNumberOfCells();
	all_pairs_face.resize(N * 4, -2);
	
	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b, id_c;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) 
	{

		for (int num_face = 0; num_face < 4; ++num_face) {
			if (all_pairs_face[num_cell * 4 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);
			id_c = idp->GetId(2);

			/*ћожет быть проблема с указател€ми на списки!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);
			int face = num_cell * 4 + num_face;

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));

				all_pairs_face[face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[face] = -1; // гранична€ €чейка
			else
			{
				RETURN_ERR("More than 1 neighbor????\n");
			}
		}

	}

	return 0;	
}
int GetNeighborFace2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face)
{
	auto GetNumberNeighborFace{ [](const int a, const int b, vtkCell* neighbor_cell)
			{

				vtkIdList* idc;

				int x, y;
				for (int i = 0; i < 3; i++)
				{
					idc = neighbor_cell->GetEdge(i)->GetPointIds();
					x = idc->GetId(0);
					y = idc->GetId(1);

					if (a == x && b == y) return i;
					else if (a == y && b == x) return i;
				}
				return -2;
			} };

	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	all_pairs_face.resize(N * 3, -2);
	
	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell)
	{

		for (int num_face = 0; num_face < 3; ++num_face) {
			if (all_pairs_face[num_cell * 3 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetEdge(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);

			/*ћожет быть проблема с указател€ми на списки!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);
			int face = num_cell * 3 + num_face;

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, unstructured_grid->GetCell(id_neighbor_cell));

				all_pairs_face[face] = id_neighbor_cell * 3 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 3 + id_neighbor_face] = face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[face] = -1; // гранична€ €чейка
			else
			{
				RETURN_ERR("More than 1 neighbor????\n");
			}				
		}

	}

	return 0;
}

static int NormalAndSquareFace3D(const int NumberCell, const int NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type& S, Vector3& n) {

	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	S = n.norm() / 2;
	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

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
	unstructuredgrid->GetPoint(idp2->GetId(id), P3);
	/*for (size_t i = 0; i < 3; i++){
		sum += n[i] * (P3[i] - P0[i]);
	}*/

	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] = -n[i];
	return 0;
}
int GetNormalAndSquares3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Normals>& normals, std::vector<Type>& area) {
	
	const int n = unstructured_grid->GetNumberOfCells();
	normals.resize(n, Normals(base));
	area.resize(n * (base));
	
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < base; j++)
		{
			NormalAndSquareFace3D(i, j, unstructured_grid, area[i* base +j], normals[i].n[j]);
		}
	}

	return 0;
}


static int NormalAndSquareFace2D(size_t NumberCell, size_t NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type& S, Vector3& n) {

	vtkIdList* idp = unstructuredgrid->GetCell(NumberCell)->GetEdge(NumberFace)->GetPointIds();

	Type P0[3], P1[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);

	Type a[2];
	for (size_t i = 0; i < 2; i++) {
		a[i] = P1[i] - P0[i];
	}
	n[0] = P0[1] - P1[1];
	n[1] = P1[0] - P0[0];
	n[2] = 0;


	S = sqrt(a[0] * a[0] + a[1] * a[1]); // длина отрезка

	n[0] /= S;
	n[1] /= S;
	//n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 3; i++) {
		int count = 0;
		for (size_t j = 0; j < 2; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 2) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructuredgrid->GetPoint(idp2->GetId(id), P3);
	for (size_t i = 0; i < 2; i++) {
		sum += n[i] * (P3[i] - P0[i]);
	}

	if (sum < 0)
		for (size_t i = 0; i < 2; i++)
			n[i] *= -1;
	return 0;
}
int GetNormalAndSquares2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Normals>& normals, std::vector<Type>& area)
{
	const int n = unstructured_grid->GetNumberOfCells();
	normals.resize(n, Normals(base));
	area.resize(n * (base));

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < base; j++)
		{
			NormalAndSquareFace2D(i, j, unstructured_grid, area[i*base+j], normals[i].n[j]);			
		}
	}

	return 0;
}


int GetVolume3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Type>& volume)
{
	auto GetVolumeCell{ [](size_t NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
		Type V = 0;
		Type P0[3], P1[3], P2[3], P3[3];

		vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();
		unstructuredgrid->GetPoint(idp->GetId(0), P0);
		unstructuredgrid->GetPoint(idp->GetId(1), P1);
		unstructuredgrid->GetPoint(idp->GetId(2), P2);
		unstructuredgrid->GetPoint(idp->GetId(3), P3);

		Vector3 a, b, c;
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
			c[i] = P3[i] - P0[i];
		}
		return -(a.dot(b.cross(c))) / 6.0;


		/*Type a[3], b[3], c[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
			c[i] = P3[i] - P0[i];
		}

		V = a[0] * (b[1] * c[2] - c[1] * b[2]) - a[1] * (b[0] * c[2] - c[0] * b[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
		return fabs(V) / 6;*/
	} };
	
	const int n = unstructured_grid->GetNumberOfCells();
	volume.resize(n);	
	
	for (size_t i = 0; i < n; i++) 
	{
		volume[i] = GetVolumeCell(i, unstructured_grid);		
	}

	return 0;
}

int GetVolume2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Type>& volume)
{

	auto GetVolumeCell{ [](size_t NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
		Type V = 0;
		Type P0[3], P1[3], P2[3];

		vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();
		unstructuredgrid->GetPoint(idp->GetId(0), P0);
		unstructuredgrid->GetPoint(idp->GetId(1), P1);
		unstructuredgrid->GetPoint(idp->GetId(2), P2);

		Eigen::Vector2d a, b;
		for (size_t i = 0; i < 2; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}
		return (Vector3(a[0],a[1],0).cross(Vector3(b[0], b[1], 0))).norm() / 2;
	} };

	const int n = unstructured_grid->GetNumberOfCells();
	volume.resize(n);

	for (size_t i = 0; i < n; i++) {
		volume[i] = GetVolumeCell(i, unstructured_grid);
	}

	return 0;
}

static size_t CenterOfTetra3D(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra) {

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
int GetCentersOfTetra3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers) 
{	
	const int n = unstructured_grid->GetNumberOfCells();
	centers.resize(n);
	
	for (int i = 0; i < n; i++) 
	{
		CenterOfTetra3D(i, unstructured_grid, centers[i]);	
	}	
	return 0;
}

static size_t CenterOfTetra2D(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& point_in_tetra) {

	Type P0[3], P1[3], P2[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(number_cell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	for (int k = 0; k < 2; k++)
	{
		point_in_tetra[k] = (P0[k] + P1[k] + P2[k]) / 3;
	}
	point_in_tetra[2] = 0;
	return 0;
}

int GetCentersOfTetra2D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers) {

	const int n = unstructured_grid->GetNumberOfCells();
	centers.resize(n);

	for (int i = 0; i < n; i++)
	{
		CenterOfTetra2D(i, unstructured_grid, centers[i]);
	}

	return 0;
}

int GetCentersOfFaces3D(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<Vector3>& centers)
{	
	int n = unstructured_grid->GetNumberOfCells();
	centers.resize(n * base);
	
	for (size_t i = 0; i < n ; i++)
	{
		for (size_t j = 0; j < base; j++)
		{
			vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds();

			Type P0[3], P1[3], P2[3];
			unstructured_grid->GetPoint(idp->GetId(0), P0);
			unstructured_grid->GetPoint(idp->GetId(1), P1);
			unstructured_grid->GetPoint(idp->GetId(2), P2);
			
			for (int k = 0; k < 3; k++)
			{
				centers[i * base + j][k] = (P0[k] + P1[k] + P2[k]) / 3;
			}
		}

	}
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
#endif // USE_VTK

int SetTypeOfBound(const std::vector<Vector3>& centers, const std::vector<Normals>& normals, std::vector<int>& all_pairs_face)
{
	const size_t n = centers.size();
	for (size_t num_cell = 0; num_cell < n; num_cell++)
	{
		Vector3 P = centers[num_cell];
		for (size_t num_face = 0; num_face < base; num_face++)
		{
			const int id = num_cell * base + num_face;
			if (all_pairs_face[id] < 0) // граница
			{

#ifdef Cone
				//if (P[0] > inner_radius)  //x>R
				if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-5)
				{
					all_pairs_face[id] = eBound_OutSource; // излучающее дно			
				}
				else if ((normals[num_cell].n[num_face] - Vector3(1, 0, 0)).norm() < 1e-5)
				{
					all_pairs_face[id] = eBound_FreeBound; // свободна€ поверхность
				}
				else
				{
					all_pairs_face[id] = eBound_LockBound; // бокова€ поверхность	
				}
#endif
#ifdef Sphere
				if ((P - center_point).norm() > inner_radius)
					all_pairs_face[id] = eBound_FreeBound; // внешн€€ сфера
				else
					all_pairs_face[id] = eBound_InnerSource; // внутренн€€ сфера		
#endif
#ifdef Cube
				if (P[0] > 0.5)
					all_pairs_face[id] = eBound_FreeBound; // внешн€€ сфера	
				else
					all_pairs_face[id] = eBound_FreeBound; // внешн€€ сфера	
#endif
#ifdef Step
				if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-5)
				{
					all_pairs_face[id] = eBound_FreeBound; // лева€ стенка	
				}
				else if ((normals[num_cell].n[num_face] - Vector3(1, 0, 0)).norm() < 1e-5 && P[0] > 0.7)
				{
					all_pairs_face[id] = eBound_FreeBound; // горизонтальные стенки
				}
				else
				{
					all_pairs_face[id] = eBound_LockBound; //ступенька	
				}


#endif
#ifdef Cylinder
				if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-5 &&
					Vector2(P[1], P[2]).norm() < 0.2)
				{
					all_pairs_face[id] = eBound_OutSource; // лева€ стенка	
				}
				else
				{
					all_pairs_face[id] = eBound_FreeBound; // лева€ стенка	
				}

#endif
} //if
		} //face
	} //cell

	return 0;
}


//

#if 0
int FromDecartToSphere(const Type* decart, Type& fi, Type& theta) {
	Type x = decart[0];
	Type y = decart[1];

	theta = atan(sqrt(x * x + y * y) / decart[2]) + PI / 2;

	if (x <= 0)
		fi = atan(y / x) + PI;
	else if (x > 0 && y < 0)
		fi = atan(y / x) + 2 * PI;
	else if (x > 0 && y >= 0)
		fi = atan(y / x);

	return 0;
}
int FromSphericalToDecart(const int number_cur, const std::vector<Type>& all_directions, Vector3& direction) {
	Type theta = all_directions[number_cur];
	Type fi = all_directions[all_directions.size() / 2 + number_cur];

	direction[0] = sin(theta) * cos(fi);
	direction[1] = sin(theta) * sin(fi);
	direction[2] = cos(theta);
	return 0;
}
#endif