#include "struct_short_characteristics_logic_function.h"

#if 0
int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const std::vector<Type>& illum_old,
	const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			X.push_back(x);
			//fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get());
			//ofile_x << x[0] << ' ' << x[1] << ' ' << x[2] << ' ';

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани		
			X.push_back(x);
			/*fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get());
			ofile_x << x[0] << ' ' << x[1] << ' ' << x[2] << ' ';*/

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			X.push_back(x);
			/*fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get());
			ofile_x << x[0] << ' ' << x[1] << ' ' << x[2] << ' ';*/

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}// x->координата узла на выходящей грани		}
	//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			X.push_back(x);
			/*fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get());
			ofile_x << x[0] << ' ' << x[1] << ' ' << x[2] << ' ';
		*/
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, illum_old, directions, squares);
		}
		//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}




	// дублирование на соседнюю ячейку
	/*int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];

	if (neighbor_id_face != -1)
		nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
			nodes_value[num_cur_cell].nodes_value[num_cur_out_face];*/

	return 0;
}

int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_s,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_in_id) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

//			X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get()); posX++;

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани		
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get()); posX++;

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get()); posX++;

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}// x->координата узла на выходящей грани		}	
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x.get()); posX++;

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x,
				file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}




	// дублирование на соседнюю ячейку
	/*int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];

	if (neighbor_id_face != -1)
		nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
			nodes_value[num_cur_cell].nodes_value[num_cur_out_face];*/

	return 0;
}


int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_res_bound,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_s,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_x0_local,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_in_id) {

	Vector3 x0;

	for (ShortId num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0)) {

			Type s = (x - x0).norm();

			fwrite_unlocked(&s, sizeof(Type), 1, file_s.get());  posS++;
			fwrite_unlocked(&num_in_face, sizeof(ShortId), 1, file_in_id.get()); posIn++;

			/*S.push_back(s);
			in_id.push_back(num_in_face);*/

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value,
				file_res_bound, file_x0_local);

			break;
		}

	}//for num_in_face


	return 0;
}

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0)) {

			Type s = (x - x0).norm();

			//S.push_back(s);

			//in_id.push_back(num_in_face);
			//ofile_in_id << num_in_face << ' '; 
		/*	ofile_s << s << ' ';
			fwrite_unlocked(&s, sizeof(Type), 1, file_s.get());*/

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value);

			//Type I = CurGetIllum(num_cur_cell, x0, s, I_x0, direction, illum_old, directions, squares);

			//nodes_value[num_cur_cell].nodes_value[num_cur_out_face][num_node] = I;

			break;
		}

	}//for num_in_face


	return 0;
}

#endif // USE_VTK


#ifndef ReadNow

int GetNodes(const int num_cur_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<cell>& nodes_value,
	std::vector<Type>& vec_res_bound, std::vector<Type>& vec_s, std::vector<Vector3>& vec_x, std::vector<Vector2>& vec_x0_local, std::vector<ShortId>& vec_in_id) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			vec_x.push_back(x);
			
			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, vec_res_bound, vec_s, vec_x0_local, vec_in_id);
		}
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани		
			
			vec_x.push_back(x);
			
			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, vec_res_bound, vec_s, vec_x0_local, vec_in_id);
		}
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			
			vec_x.push_back(x);			

			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, vec_res_bound, vec_s, vec_x0_local, vec_in_id);
		}// x->координата узла на выходящей грани		}	
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			vec_x.push_back(x);
			
			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x,
				vec_res_bound, vec_s, vec_x0_local, vec_in_id);
		}
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}

	posX += 3;


	// дублирование на соседнюю ячейку
	/*int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];

	if (neighbor_id_face != -1)
		nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
			nodes_value[num_cur_cell].nodes_value[num_cur_out_face];*/

	return 0;
}

int CalculateNodeValue(const int num_cur_cell, const std::vector<Normals>& normals, const std::vector<Face>& grid,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	std::vector<Type>& vec_res_bound, std::vector<Type>& vec_s, std::vector<Vector2>& vec_x0_local, std::vector<ShortId>& vec_in_id) {

	Vector3 x0;

	for (ShortId num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани

		IntersectionWithPlane(grid[num_cur_cell * 4 + num_in_face], x, direction, x0);

		if (InTriangle(num_in_face, grid[num_cur_cell * 4 + num_in_face], normals[num_cur_cell], x0)) {

			Type s = (x - x0).norm(); 

			/*fwrite_unlocked(&s, sizeof(Type), 1, file_s);  posS++;
			fwrite_unlocked(&num_in_face, sizeof(ShortId), 1, file_in_id); posIn++;*/

			vec_s.push_back(s); posS++;
			vec_in_id.push_back(num_in_face);  posIn++;

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value,
				vec_res_bound, vec_x0_local);

			break;
		}

	}//for num_in_face


	return 0;
}

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::vector<cell>& nodes_value,
	std::vector<Type>& vec_res_bound, std::vector<Vector2>& vec_x0_local) {

	Type I_x0 = 0;
	const int neib_id = nodes_value[num_cell].neighbours_id_face[num_in_face];
	
	if (neib_id == eBound_FreeBound || neib_id == eBound_LockBound)
	{
		/*Граничные условия*/
		return I_x0;
	}
	else if (neib_id == eBound_OutSource) // граница основания конуса
	{
		return I_x0;
	}
	else if (neib_id == eBound_InnerSource) {

		I_x0 = 10;
		vec_res_bound.push_back(I_x0);
		//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound); 
		posRes++;
		return I_x0;

		// внутренняя граница (//пересечение с диском / сферой)		
		pos_x_try++;
		static int pos_id_try = 0;
		pos_id_try++;

		Vector3 res;
		// пересечние луча с плоскостью диска
		IntersectionWithPlaneDisk(x, cur_direction, res);

		const Vector3 v1(1, 0, 0);
		const Vector3 v2(0, -0.992877, -0.119145); // Wolfram

		//в плоскости
		Vector3 LocRes(res.dot(v1), res.dot(v2), 0);  // точка пересечения в локальных координатах плоскости
		LocRes -= center_point;
		const Type dist = LocRes.dot(LocRes);

		const Type A = cur_direction.dot(cur_direction);

		const Type buf = (x - center_point).dot(cur_direction);

		const Type radical = 3 * buf * buf - 4 * A * (1 - 2 * x[0] + x.dot(x) - Rsphere * Rsphere);

		// есть пересечение со сферой
		if (radical >= 0)
		{
			const Type t = (cur_direction[0] - cur_direction.dot(x) - sqrt(radical) / 2) / A;

			const Vector3 inSphere = cur_direction * t + x;

			// не пересекает плоскость			
			if (dist <= R1disk * R1disk || dist >= R2disk * R2disk) {
				I_x0 = 50;
				//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound); 
				vec_res_bound.push_back(I_x0);
				posRes++;
				return I_x0; // 2;
			}

			// с чем луч встречается раньше?
			{
				const Vector3 Xsphere = x - inSphere;
				const Vector3 Xres = x - res;
				const Type LenSpehere = Xsphere.dot(Xsphere);
				const Type LenPlane = Xres.dot(Xres);

				if (LenSpehere > LenPlane) { //!!!!!!!!!!!!!!>
					I_x0 = 20;
					//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound);
					posRes++;
					vec_res_bound.push_back(I_x0);
					return I_x0; // 1;
				}
				else {
					I_x0 = 50;
					//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound); 
					vec_res_bound.push_back(I_x0);
					posRes++;
					return I_x0; // 2;
				}
			}
		}
		else if ((dist < R2disk * R2disk) && (dist > R1disk * R1disk)) {
			I_x0 = 20;
			//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound); 
			vec_res_bound.push_back(I_x0);
			posRes++;
			return I_x0; // 1;
		}
		else // внутренняя граница не пересекла ни диск ни сферу 
		{

			Vector3 try_x0;
			if (pos_x_try - 1 >= x_try_surface.size()) printf("err size try_id\n");
			FromGlobalToLocalTetra(vertex_tetra, x_try_surface[pos_x_try - 1], try_x0);


			switch (id_try_surface[pos_id_try - 1] % 4) {
			case 3:
				Vector3 local_plane_x0;
				FromTetraToPlane(transform_matrix, start_point_plane_coord, try_x0, local_plane_x0);
				//fwrite_unlocked(Vector2(local_plane_x0[0], local_plane_x0[1]).data(), sizeof(Type), 2, file_x0_local);
				vec_x0_local.push_back(Vector2(local_plane_x0[0], local_plane_x0[1]));
				posX0++;
				break;

			case 1:
				vec_x0_local.push_back(Vector2(try_x0[1], try_x0[2]));
				//fwrite_unlocked(Vector2(try_x0[1], try_x0[2]).data(), sizeof(Type), 2, file_x0_local);
				posX0++;
				break;
			case 2:
				vec_x0_local.push_back(Vector2(try_x0[0], try_x0[2]));
				//fwrite_unlocked(Vector2(try_x0[0], try_x0[2]).data(), sizeof(Type), 2, file_x0_local);
				posX0++;
				break;
			case 0:
				vec_x0_local.push_back(Vector2(try_x0[0], try_x0[1]));
				//fwrite_unlocked(Vector2(try_x0[0], try_x0[1]).data(), sizeof(Type), 2, file_x0_local);
				posX0++;
				break;
			}

			I_x0 = -10;
			vec_res_bound.push_back(I_x0);
			//fwrite_unlocked(&I_x0, sizeof(Type), 1, file_res_bound); 
			posRes++;
			return I_x0; // 2;
		}

	}
	else
	{

		Vector3 x0_local;
		FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);

		switch (num_in_face) 
		{
		case 3:
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);
			vec_x0_local.push_back(Vector2(local_plane_x0[0], local_plane_x0[1]));
			//fwrite_unlocked(Vector2(local_plane_x0[0], local_plane_x0[1]).data(), sizeof(Type), 2, file_x0_local);
			posX0++;

			/*I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];
			Vector3 coef = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);*/
			break;
		case 1:
			vec_x0_local.push_back(Vector2(x0_local[1], x0_local[2]));
			//fwrite_unlocked(Vector2(x0_local[1], x0_local[2]).data(), sizeof(Type), 2, file_x0_local);
			posX0++;
			break;
		case 2:
			vec_x0_local.push_back(Vector2(x0_local[0], x0_local[2]));
			//fwrite_unlocked(Vector2(x0_local[0], x0_local[2]).data(), sizeof(Type), 2, file_x0_local);
			posX0++;
			break;
		case 0:
			vec_x0_local.push_back(Vector2(x0_local[0], x0_local[1]));
			//fwrite_unlocked(Vector2(x0_local[0], x0_local[1]).data(), sizeof(Type), 2, file_x0_local);
			posX0++;
			break;
		}

		return I_x0;
	}
}


#else
int GetNodes(const int num_cur_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<cell>& nodes_value,
	FILE* file_res_bound, FILE* file_s, FILE* file_x, FILE* file_x0_local, FILE* file_in_id) {

	Vector3 x;
	Vector3 node;

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

//			X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x); posX++;

			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани		
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x); posX++;

			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x); posX++;

			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x, file_res_bound, file_s, file_x0_local, file_in_id);
		}// x->координата узла на выходящей грани		}	
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			//X.push_back(x);
			fwrite_unlocked(x.data(), sizeof(Type), 3, file_x); posX++;

			CalculateNodeValue(num_cur_cell, normals, grid, num_cur_out_face, face_state, direction, nodes_value,
				num_node, vertex_tetra, x,
				file_res_bound, file_s, file_x0_local, file_in_id);
		}
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}




	// дублирование на соседнюю ячейку
	/*int neighbor_id_face = nodes_value[num_cur_cell].neighbours_id_face[num_cur_out_face];

	if (neighbor_id_face != -1)
		nodes_value[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] =
			nodes_value[num_cur_cell].nodes_value[num_cur_out_face];*/

	return 0;
}


int CalculateNodeValue(const int num_cur_cell, const std::vector<Normals>& normals, const std::vector<Face>& grid,
	const int num_cur_out_face, const int* face_state, const Vector3& direction,
	std::vector<cell>& nodes_value, const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	FILE* file_res_bound, FILE* file_s, FILE* file_x0_local, FILE* file_in_id) {

	Vector3 x0;

	for (ShortId num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани

		IntersectionWithPlane(grid[num_cur_cell * 4 + num_in_face], x, direction, x0);

		if (InTriangle(num_in_face, grid[num_cur_cell * 4 + num_in_face], normals[num_cur_cell], x0)) {

			Type s = (x - x0).norm();

			fwrite_unlocked(&s, sizeof(Type), 1, file_s);  posS++;
			fwrite_unlocked(&num_in_face, sizeof(ShortId), 1, file_in_id); posIn++;

			/*S.push_back(s);
			in_id.push_back(num_in_face);*/

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_cur_cell, num_in_face, vertex_tetra, x, x0, nodes_value,
				file_res_bound, file_x0_local);

			break;
		}

	}//for num_in_face


	return 0;
}

int CalculateIllumeOnInnerFace(const int num_cur_cell, const ShortId num_in_face, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	Vector3& x0, std::vector<cell>& nodes_value,
	FILE*file_res_bound, FILE* file_x0_local) {
	printf("Foo CalculateIllumeOnInnerFace \n");
	return 0;
}



#endif
