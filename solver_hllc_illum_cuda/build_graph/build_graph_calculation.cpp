#ifdef BUILD
#include "build_graph_calculation.h"
#include "../utils/grid_geometry/geometry_data.h"
#include "../utils/grid_geometry/geometry_solve.h"

//-------------------------------------------------------------------------------------

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<Normals>& normals, const int cur_cell, int* id) {

	// вершины треугольника.
	Face face = inner_cells.find(cur_cell)->second.face;
	Vector3 P1(face.A.data());
	Vector3 P2(face.B.data());
	Vector3 P3(face.C.data());

	// середины сторон (противолежащий точке P1,P2,P3 соответственно)
	Vector3 P11 = (P3 + P2) / 2;
	Vector3 P22 = (P3 + P1) / 2;
	Vector3 P33 = (P2 + P1) / 2;

	// точки на медианах
	Vector3 vertex1 = P1 + (P11 - P1) / 3;
	Vector3 vertex2 = P2 + (P22 - P2) / 3;
	Vector3 vertex3 = P3 + (P33 - P3) / 3;

	// ищем пересечения vectrex->direction c гранями внутренней границы
	//Vector3 intersect_point1;  //buf_try.x1
	//Vector3 intersect_point2;  //buf_try.x2
	//Vector3 intersect_point3;  //buf_try.x3

	FaceCell plane_face;
	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex1, direction, buf_try.x1);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x1))
			//if (in_id != cur_cell && ((intersect_point - vertex1).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/

			buf_try.id_1 = plane_face.face_id;
			id[0] = in_id;
			break;
		}
	}

	if (id[0] == -1) return 1;

	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex2, direction, buf_try.x2);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x2))
			//if (in_id != cur_cell && ((intersect_point - vertex2).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/
			id[1] = in_id;
			buf_try.id_2 = plane_face.face_id;
			break;
		}
	}

	if (id[1] == -1) return 1;

	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex3, direction, buf_try.x3);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x3))
			//if (in_id != cur_cell && ((intersect_point - vertex3).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/
			id[2] = in_id;
			buf_try.id_3 = plane_face.face_id;
			break;
		}
	}

	if (id[2] == -1) return 1;

	buf_try.s_1 = (vertex1 - buf_try.x1).norm();
	buf_try.s_2 = (vertex2 - buf_try.x2).norm();
	buf_try.s_3 = (vertex3 - buf_try.x3).norm();


	return 0;
}
int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, std::bitset<4>& face_state) {
	//face_state  -0=> выходящая грань,  1=> входящая  face_state.size=4!!!  

	Normals cell = normals[number_cell];
	const Type eps = 1e-10;
	for (size_t i = 0; i < 4; ++i) {

		if (cell.n[i].dot(direction) < -eps)
			face_state.set(i);
		else
			face_state.set(i, false);
	}

	return 0;
}
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const std::map<IntId,FaceCell>& inner_boundary_faces) {

	const int n = all_pairs_id.size();
	faces_state.assign(n, 0);

	for (size_t i = 0; i < n; i++) {
		if (all_pairs_id[i] < 0)
			faces_state[i] = 1;
	}

	for (auto& el : inner_boundary_faces)
		faces_state[el.second.face_id] = 0;

	return 0;
}

int FractionInnerBoundary(const Vector3& direction, const std::vector<Normals>& normals, const std::map<IntId, FaceCell>& inner_cells,
	const std::set<IntId>& full_boundary, std::set<IntId>& inner_part, std::set<IntId>& outter_part) 
{
	inner_part.clear();
	outter_part.clear();

	std::bitset<4> state;

	int bound_face_id = -1;
	for (auto el : full_boundary)
	{
		FindInAndOutFaces(direction, el, normals, state);

		bound_face_id = inner_cells.find(el)->second.face_id;

		if (state[bound_face_id % 4] != 0)  // выходящая грань
			outter_part.emplace(el);
		else
			inner_part.emplace(el);
	}

	return 0;
}

//=======================================OMP=======================================
#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY	
int FindCurCellWithHole(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el,
	const std::set<IntId>& inner_part, std::set<IntId>& outter_part, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();

	const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face,cur_el,inner_part,outter_part, all_pairs_id, direction,normals,inner_cells, N)  
	{
		int cell;
#pragma omp for
		for (int i = 0; i < N; ++i)
		{
			cell = next_step_el[i];
			if (outter_part.count(cell) != 0) {
				if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
					IntId try_id[3] = { -1,-1,-1 };

					FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);

					if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1)
						int a;
					else if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
						count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
						count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена
#pragma omp critical 
							{
								if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
									count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
									count_in_face[try_id[2]] == count_knew_face[try_id[2]]) {
									cur_el.push_back(cell);
									outter_part.erase(cell);
									count_knew_face[cell]++;
								}
							}
					}
				}
				else if (count_in_face[cell] == count_knew_face[cell]) {  // граница определена
#pragma omp critical 
					{
						if (count_in_face[cell] == count_knew_face[cell]) {
							cur_el.push_back(cell);
							outter_part.erase(cell);
						}
					}
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical 
				{
					if (count_in_face[cell] == count_knew_face[cell])
						cur_el.push_back(cell);
				}
			}
		}


	}

	if (cur_el.size() == 0) {

		// плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell);//outter_part.erase(cell);
				cur_el.push_back(cell);   //добавляет все найденные ранее границы!!!
			}
		}

		for (auto el : buf_erase)
			outter_part.erase(el);

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}
#else
int FindCurCell(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el) {

	cur_el.clear();

	const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face,cur_el,  N)  
	{
#pragma omp for
		for (int i = 0; i < N; i++)
		{
			int cell = next_step_el[i];
			if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical 
				{
					if (count_in_face[cell] == count_knew_face[cell])
						cur_el.push_back(cell);
				}
			}
		}
	}

	if (cur_el.size() == 0) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	return 0;
}
#endif

// число входящих граней для каждой ячейки + число известных из них + начальная граница
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& next_step_el) {

	const int N = normals.size();  // число ячеек

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// для i-ой ячейки добавляем:
				count_in_face[i]++;  // входную грань
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //определнную входную грань (т.е. граничную)
					next_step_el.push_back(i);  // начальный набор кандидатов
				}
			}
		}
	}
	return 0;
}

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::list<IntId>& cur_el,
	std::vector<IntId>& next_step_el) {

	int buf_size = next_step_el.size();
	next_step_el.clear();
	next_step_el.reserve(buf_size); // резервируем память на основе предыдущего шага(предпологая, что порядок величины будет тот же)

	int N = cur_el.size();

	std::set<IntId> next_step;

	for (auto cell : cur_el)
	{
		// по всем соседям
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour  < 0) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				count_knew_face[neighbour]++;  // всегда ли эта грань будет входящей (проверка по нормалям??)
				if (next_step.count(neighbour) == 0) {
					next_step.emplace(neighbour);
					next_step_el.push_back(neighbour);  // ячейка была изменена, проверить ее готовность на следующем шаге
				}
			}
		}

	}

	return 0;
}

#else

// число входящих граней для каждой ячейки + число известных из них + начальная граница
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::set<IntId>& next_step_el) {

	const int N = normals.size();  // число ячеек

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// для i-ой ячейки добавляем:
				count_in_face[i]++;  // входную грань
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //определнную входную грань (т.е. граничную)
					next_step_el.emplace(i);  // начальный набор кандидатов
				}
			}
		}
	}
	return 0;
}

int FindCurCell(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el) {

	cur_el.clear();

	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (count_in_face[cell] == count_knew_face[cell])
			cur_el.push_back(cell);
	}

	if (cur_el.size() == 0) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	return 0;
}

int FindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el, const std::set<IntId>& inner_part, std::set<IntId>& outter_part,
	const std::map<IntId, FaceCell>& inner_cells, const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();


	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (outter_part.count(cell) != 0) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };

				
				if (FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id)) continue;

				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					//id_try_surface.push_back(buf_try.id_1);
					//id_try_surface.push_back(buf_try.id_2);
					//id_try_surface.push_back(buf_try.id_3);

					//x_try_surface.push_back(buf_try.x1);
					//x_try_surface.push_back(buf_try.x2);
					//x_try_surface.push_back(buf_try.x3);

					//dist_try_surface.push_back(buf_try.s_1);
					//dist_try_surface.push_back(buf_try.s_2);
					//dist_try_surface.push_back(buf_try.s_3);

					cur_el.push_back(cell);
					outter_part.erase(cell);
					count_knew_face[cell]++;
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {  // граница определена
				cur_el.push_back(cell);
				outter_part.erase(cell);
			}
		}
		else if (count_in_face[cell] == count_knew_face[cell]) {
			cur_el.push_back(cell);
		}
	}

	if (cur_el.size() == 0) {

		// плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell); //outter_part.erase(cell);
				cur_el.push_back(cell);
			}
		}

		for (auto el : buf_erase)
			outter_part.erase(el);

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	const std::vector<IntId>& cur_el,
	std::set<IntId>& next_step_el) {

	next_step_el.clear();
	int N = cur_el.size();

	for (size_t i = 0; i < N; i++)
	{
		int cell = cur_el[i];
		//	count_knew_face[cell]++;
			// по всем соседям
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour < 0) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				next_step_el.emplace(neighbour);  // ячейка была изменена, проверить ее готовность на следующем шаге
				count_knew_face[neighbour]++;  // всегда ли эта грань будет входящей (проверка по нормалям??)
			}
			//else if(count_in_face[neighbour] == count_knew_face[neighbour]) count_knew_face[neighbour]++; //попытка настроить restart
		}

	}

	return 0;
}

#endif // USE_OMP

#endif //BUILD