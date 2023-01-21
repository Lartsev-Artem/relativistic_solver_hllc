#ifndef BUILD_GRAPH_CALCULATION
#define BUILD_GRAPH_CALCULATION

#include "build_graph_prj_config.h"

#ifdef BUILD

#include "build_graph_structures.h"

#include <map>
#include<set>
#include <bitset>
#include<list>

#ifdef  USE_OMP
#include <omp.h>
#endif 

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<Normals>& normals, const int cur_cell, int* id);

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, std::bitset<4>& face_state);
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const  std::map<IntId, FaceCell>& inner_boundary_faces);

int FractionInnerBoundary(const Vector3& direction, const std::vector<Normals>& normals, const std::map<IntId, FaceCell>& inner_cells,
	const std::set<IntId>& full_boundary, std::set<IntId>& inner_part, std::set<IntId>& outter_part);

//=======================================OMP=======================================
#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY
int FindCurCellWithHole(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el, const std::set<IntId>& inner_part, std::set<IntId>& outter_part,
	const std::map<IntId, FaceCell>& inner_cells, const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals);
#else
int FindCurCell(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el);
#endif

int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& next_step_el);

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::list<IntId>& cur_el,
	std::vector<IntId>& next_step_el);

#else

#ifdef GRID_WITH_INNER_BOUNDARY
int FindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el,
	const std::set<IntId>& inner_part, std::set<IntId>& outter_part, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals);
#else
int FindCurCell(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el);
#endif

int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::set<IntId>& next_step_el);

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, const std::vector<IntId>& cur_el,
	std::set<IntId>& next_step_el);

#endif //USE_OMP

#endif //BUILD

#endif

