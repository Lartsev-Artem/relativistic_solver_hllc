#include "build_graph_prj_config.h"

#ifdef BUILD
#include "../global_value.h"

#ifdef WriteFiles
#include "../file_module/reader_vtk.h"

#include "../utils/grid_geometry/geometry_data.h"
#include "../file_module/writer_bin.h"
#include "../file_module/reader_bin.h"


static int ReWritePairsByType(const std::string& name_file_pairs, const std::string& name_file_normals, const std::string& name_file_centers)
{
	std::vector<int>all_pairs_face;
	std::vector<Normals>normals;
	std::vector<Vector3>centers;

	ReadSimpleFileBin(name_file_centers, centers);
	ReadSimpleFileBin(name_file_pairs, all_pairs_face);
	ReadNormalFile(name_file_normals, normals);

	SetTypeOfBound(centers, normals, all_pairs_face);

	WriteSimpleFileBin(name_file_pairs, all_pairs_face);

	return 0;
}
static int WritePairsId(const std::string& name_file_pairs, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid)
{
	std::vector<int>all_pairs_face;

#ifdef TASK_3D
	if (GetNeighborFace3D(unstructured_grid, all_pairs_face)) RETURN_ERR("Neighbor wasn't found\n");
#elif defined TASK_2D
	if (GetNeighborFace2D(unstructured_grid, all_pairs_face)) RETURN_ERR("Neighbor wasn't found\n");
#endif
	
	WriteSimpleFileBin(name_file_pairs, all_pairs_face);
	return 0;
}
static int WriteNormalAndSquaresFile(const std::string& name_file_normals, const std::string name_file_squares, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::vector<Normals> normals;
	std::vector<Type> area;
#ifdef TASK_3D
	GetNormalAndSquares3D(unstructured_grid, normals, area);
#else
	GetNormalAndSquares2D(unstructured_grid, normals, area);
#endif

	WriteSimpleFileBin(name_file_squares, area);

	std::unique_ptr<FILE, int(*)(FILE*)> file_norm(fopen(name_file_normals.c_str(), "wb"), fclose);
	if (!file_norm) { printf("file normals is not opened for writing\n"); return 1; }


	int n = unstructured_grid->GetNumberOfCells();
	fwrite(&n, sizeof(int), 1, file_norm.get());
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < base; j++)
		{			
			//fwrite(normals.data(), sizeof(Normals), n, file_norm.get());
			fwrite(&normals[i].n[j], sizeof(Vector3), 1, file_norm.get());
		}
	}
	fclose(file_norm.get());
	return 0;
}
static int WriteVolume(const std::string name_file_volume, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) 
{

	std::vector<Type>volume;

#ifdef TASK_3D
	if (GetVolume3D(unstructured_grid, volume)) RETURN_ERR("volume wasn't found\n");
#elif defined TASK_2D
	if (GetVolume2D(unstructured_grid, volume)) RETURN_ERR("volume wasn't found\n");
#endif

	WriteSimpleFileBin(name_file_volume, volume);	
	return 0;
}
static int WriteCentersOfTetra(const std::string name_file_centers, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::vector<Vector3>centers;

#ifdef TASK_3D
	if (GetCentersOfTetra3D(unstructured_grid, centers)) RETURN_ERR("centers wasn't found\n");
#elif defined TASK_2D
	if (GetCentersOfTetra2D(unstructured_grid, centers)) RETURN_ERR("centers wasn't found\n");
#endif

	WriteSimpleFileBin(name_file_centers, centers);
	
	return 0;
}
#ifdef TASK_3D

static int WriteInitBoundarySetAndInnerBoundaryFace(const std::string name_file_boundary, const std::string name_file_boundary_inner,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::set<IntId> boundary_cells;
	std::set<IntId> inner_boundary_faces;

	boundary_cells.clear();
	inner_boundary_faces.clear();

	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {

				Vector3 P(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(0));
				if ((P - center_point).norm() > inner_radius) { // внешняя сфера
					boundary_cells.emplace(i);
				}
				else {
					inner_boundary_faces.emplace(i * 4 + j);
					boundary_cells.emplace(i);
				}
				break;
			}
		}
	}


	std::ofstream ofile;

	ofile.open(name_file_boundary);
	if (!ofile.is_open()) {
		std::cout << "Error write file boundary\n";
		return 1;
	}

	ofile << boundary_cells.size() << '\n';
	for (auto el : boundary_cells)
		ofile << el << '\n';
	ofile.close();


	ofile.open(name_file_boundary_inner);
	if (!ofile.is_open()) {
		std::cout << "Error write file boundary_inner\n";
		return 1;
	}

	ofile << inner_boundary_faces.size() << '\n';
	for (auto el : inner_boundary_faces)
		ofile << el << '\n';
	ofile.close();

	return 0;
}


static int WriteCentersOfSquares(const std::string name_file_centers_squares, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::vector<Vector3>centers_squares;

	if (GetCentersOfFaces3D(unstructured_grid, centers_squares)) RETURN_ERR("centers_squares wasn't found\n");

	WriteSimpleFileBin(name_file_centers_squares, centers_squares);	
	return 0;
}

static int ReadInnerBoundaryForWrite(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face) {

	std::ifstream ifile;

	ifile.open(name_file_boundary_inner);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary_inner\n";
		return 1;
	}

	id_inner_boundary_face.clear();

	IntId buf;
	ifile >> buf;
	std::cout << "Inner boundary has " << buf << "faces\n";
	while (ifile >> buf)
	{
		id_inner_boundary_face.emplace(buf);
	}

	ifile.close();
	return 0;
}

static int WriteInnerCellOfSphere(const std::string name_file_inner_sphere, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const std::string name_file_boundary_inner) {



	std::set<IntId> id_inner_boundary_face; // хранит номера граней в формате num_cell*4+num_face

	if (ReadInnerBoundaryForWrite(name_file_boundary_inner, id_inner_boundary_face)) return 1;


	std::ofstream ofile;

	ofile.open(name_file_inner_sphere);
	if (!ofile.is_open()) {
		std::cout << "Error write inner_sphere\n";
		return 1;
	}

	Type A[3];
	Type B[3];
	Type C[3];
	vtkPoints* points_face;

	ofile << id_inner_boundary_face.size() << '\n';

	for (auto el : id_inner_boundary_face) {
		points_face = unstructured_grid->GetCell(el / 4)->GetFace(el % 4)->GetPoints();

		points_face->GetPoint(0, A);
		points_face->GetPoint(0, B);
		points_face->GetPoint(0, C);

		ofile << setprecision(16) << A[0] << ' ' << A[1] << ' ' << A[2]
			<< B[0] << ' ' << B[1] << ' ' << B[2]
			<< C[0] << ' ' << C[1] << ' ' << C[2] << '\n';
	}

	ofile.close();
	return 0;
}

static int WriteInnerCellAndIdFaces(const std::string name_file_boundary_inner, const std::string name_file_face_and_id,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::set<IntId> id_inner_boundary_face; // хранит номера граней в формате num_cell*4+num_face
	if (ReadInnerBoundaryForWrite(name_file_boundary_inner, id_inner_boundary_face)) return 1;

	std::ofstream ofile;
	ofile.open(name_file_face_and_id);
	if (!ofile.is_open()) {
		std::cout << "Error write inner_sphere\n";
		return 1;
	}

	Type A[3];
	Type B[3];
	Type C[3];
	vtkPoints* points_face;

	ofile << id_inner_boundary_face.size() << '\n';

	for (auto el : id_inner_boundary_face) {
		points_face = unstructured_grid->GetCell(el / 4)->GetFace(el % 4)->GetPoints();

		points_face->GetPoint(0, A);
		points_face->GetPoint(1, B);
		points_face->GetPoint(2, C);

		ofile << setprecision(16) << el << ' ' << A[0] << ' ' << A[1] << ' ' << A[2] << ' '
			<< B[0] << ' ' << B[1] << ' ' << B[2] << ' '
			<< C[0] << ' ' << C[1] << ' ' << C[2] << '\n';
	}

	ofile.close();

	return 0;
}

#endif //TASK_3D

int BuildSetForClaster(const std::string name_file_vtk, const std::string name_file_pairs,
	const std::string name_file_boundary, const std::string name_file_normals, const std::string name_file_boundary_inner,
	const std::string name_file_face_and_id, const std::string name_file_squares, const std::string name_file_volume,
	const std::string name_file_centers, const std::string name_file_centers_faces)
{

	std::cout << "Start build \n";

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_vtk.c_str(), unstructured_grid)) 
	{
		RETURN_ERR("Error reading file vtk\n");		
	}

	if (WritePairsId(name_file_pairs, unstructured_grid)) return 1;

	if (WriteNormalAndSquaresFile(name_file_normals, name_file_squares, unstructured_grid)) return 1;

	if (WriteVolume(name_file_volume, unstructured_grid)) return 1;

	if (WriteCentersOfTetra(name_file_centers, unstructured_grid)) return 1;

#ifdef TASK_3D

	if (WriteInitBoundarySetAndInnerBoundaryFace(name_file_boundary, name_file_boundary_inner, unstructured_grid)) return 1;
	
	if (WriteCentersOfSquares(name_file_centers_faces, unstructured_grid)) return 1;		

	//if (WriteInnerCellOfSphere(name_file_inner_sphere, unstructured_grid, name_file_boundary_inner)) return 1;

	if (WriteInnerCellAndIdFaces(name_file_boundary_inner, name_file_face_and_id, unstructured_grid)) return 1;

#endif //TASK_3D

	ReWritePairsByType(name_file_pairs, name_file_normals, name_file_centers);

	std::cout << "Build end\n";
	return 0;
}

#endif //WriteFiles


#ifdef ReadFiles
#include <map>
#include "build_graph_structures.h"

int ReadInnerCellBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_cell) {

	std::ifstream ifile;

	OPEN_FSTREAM((ifile), name_file_boundary_inner.c_str());

	id_inner_boundary_cell.clear();

	IntId buf;
	ifile >> buf;
	//	std::cout << "Inner boundary has " << buf << "faces\n";
	while (ifile >> buf)
	{
		id_inner_boundary_cell.emplace(buf / 4);
	}

	ifile.close();
	return 0;
}
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, FaceCell>& inner_cells) {

	std::ifstream ifile;

	ifile.open(name_file_face_and_id);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	Face face;
	IntId id;
	for (int i = 0; i < N; ++i) {
		ifile >> id;
		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_cells.emplace(id / 4, FaceCell(id, face));
	}

	ifile.close();
	return 0;
}

#ifdef USE_STRANGE_FUNCTION
int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId>& boundary_cells) {

	std::ifstream ifile;

	ifile.open(name_file_boundary);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary\n";
		return 1;
	}

	IntId N;
	ifile >> N;
	//std::cout << "Boundary has " << N << "cells\n";
	while (ifile >> N) {
		boundary_cells.emplace(N);
	}

	ifile.close();
	return 0;
}
int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face) {

	std::ifstream ifile;

	ifile.open(name_file_boundary_inner);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary_inner\n";
		return 1;
	}

	id_inner_boundary_face.clear();

	IntId buf;
	ifile >> buf;
//	std::cout << "Inner boundary has " << buf << "faces\n";
	while (ifile >> buf)
	{
		id_inner_boundary_face.emplace(buf);
	}

	ifile.close();
	return 0;
}
int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face>& inner_faces) {

	std::ifstream ifile;

	ifile.open(name_file_inner_sphere);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	inner_faces.resize(N);

	Face face;
	for (int i = 0; i < N; ++i){

		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_faces[i] = face;
	}

	ifile.close();
	return 0;
}
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId,Face>& inner_faces) {

	std::ifstream ifile;

	ifile.open(name_file_face_and_id);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	Face face;
	IntId id;
	for (int i = 0; i < N; ++i) {
		ifile >> id;
		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_faces.emplace(id, face);
	}

	ifile.close();
	return 0;
}
#endif


int WriteFileGraph(const int i, const std::string& name_file_graph, const std::vector<IntId>& graph) 
{

	std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen((name_file_graph + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

	const int n = graph.size();
	fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

	fclose(file_graph.get());


	std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_id) { printf("file_id is not opened for writing\n"); return 1; }

	int size = id_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_id.get());
	fwrite_unlocked(id_try_surface.data(), sizeof(IntId), size, file_id.get());

	fclose(file_id.get());


	std::unique_ptr<FILE, int(*)(FILE*)> file_dist(fopen((std::string(BASE_ADRESS) + "dist_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_dist) { printf("file_dist is not opened for writing\n"); return 1; }

	size = dist_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_dist.get());
	fwrite_unlocked(dist_try_surface.data(), sizeof(Type), size, file_dist.get());

	fclose(file_dist.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((std::string(BASE_ADRESS) + "x_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }

	size = x_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_x.get());
	fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), size, file_x.get());

	fclose(file_x.get());


	return 0;
}

int WriteFileGraph(std::unique_ptr<FILE, int(*)(FILE*)>& file_graph, std::unique_ptr<FILE, int(*)(FILE*)>& file_id,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_dist, std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	const int i, const int n, const std::vector<IntId>& graph) {

	fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

	id_try_size += id_try_surface.size();	
	fwrite_unlocked(id_try_surface.data(), sizeof(IntId), id_try_surface.size(), file_id.get());


	dist_try_size += dist_try_surface.size();	
	fwrite_unlocked(dist_try_surface.data(), sizeof(Type), dist_try_surface.size(), file_dist.get());


	x_try_size += x_try_surface.size();	
	fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), x_try_surface.size(), file_x.get());

	return 0;
}

#endif //ReadFiles
#endif // BUILD