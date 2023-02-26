#if !defined BUILD_GRAPH_READ_WRITE && defined BUILD
#define BUILD_GRAPH_READ_WRITE
#include "../global_def.h"

#include "build_graph_prj_config.h"
#include "build_graph_structures.h"
#include <map>

#ifdef WriteFiles

int BuildSetForClaster(const std::string name_file_vtk, const std::string name_file_pairs,
	const std::string name_file_boundary, const std::string name_file_normals, const std::string name_file_boundary_inner,
	const std::string name_file_face_and_id, const std::string name_file_squares, const std::string name_file_volume, 
	const std::string name_file_centers, const std::string name_file_centers_faces);
#endif

#ifdef ReadFiles

int ReadInnerCellBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_cell);
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, FaceCell>& inner_cells);

int WriteFileGraph(const int i, const std::string& name_file_graph, const std::vector<IntId>& graph);

int WriteFileGraph(std::unique_ptr<FILE, int(*)(FILE*)>& file_graph, std::unique_ptr<FILE, int(*)(FILE*)>& file_id,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_dist, std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	const int i, const int n, const std::vector<IntId>& graph);


#ifdef USE_STRANGE_FUNCTION
int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId>& boundary_cells);
int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face);
int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face>& inner_faces);
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, Face>& inner_faces);
#endif

#endif

#endif //BUILD_GRAPH_READ_WRITE