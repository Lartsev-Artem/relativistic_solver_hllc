#include "struct_short_characteristics_global_structure.h"
#ifdef MAKE
#include "struct_short_characteristics_calculations.h"
#include "../global_value.h"
#include "../global_def.h"
#include "../utils/grid_geometry/geometry_solve.h"

#ifndef ONLY_GEO_DATA
int MakeArrayX(const std::vector<Eigen::Matrix4d>& vertexs, std::vector<BasePointTetra>& vec_x)
{

	const int count_cells = vertexs.size();
	BasePointTetra x;
	Vector3 node;
	Eigen::Matrix4d vertex_tetra;

	vec_x.reserve(count_cells);

	for (int num_cell = 0; num_cell < count_cells; ++num_cell)
	{
		vertex_tetra = vertexs[num_cell];
		for (int num_cur_out_face = 0; num_cur_out_face < base; num_cur_out_face++)
		{
			switch (num_cur_out_face)
			{

			case 1:// 1->2
				for (size_t num_node = 0; num_node < 3; ++num_node) {
					node[0] = 0;
					node[1] = base_tetra_geo.straight_face.row(num_node)[0];
					node[2] = base_tetra_geo.straight_face.row(num_node)[1];
					FromLocalToGlobalTetra(vertex_tetra, node, x.x[num_cur_out_face][num_node]);  // x->координата узла на выходящей грани					
				}
				break;
			case 2://2->0
				for (size_t num_node = 0; num_node < 3; ++num_node) {
					node[0] = base_tetra_geo.straight_face.row(num_node)[0];
					node[1] = 0;
					node[2] = base_tetra_geo.straight_face.row(num_node)[1];
					FromLocalToGlobalTetra(vertex_tetra, node, x.x[num_cur_out_face][num_node]);  // x->координата узла на выходящей грани		

				}
				break;
			case 0: //0->3
				for (size_t num_node = 0; num_node < 3; ++num_node) {
					node[0] = base_tetra_geo.straight_face.row(num_node)[0];
					node[1] = base_tetra_geo.straight_face.row(num_node)[1];
					node[2] = 0;
					FromLocalToGlobalTetra(vertex_tetra, node, x.x[num_cur_out_face][num_node]);

				}// x->координата узла на выходящей грани		}	
				break;
			case 3: //3->1
				for (size_t num_node = 0; num_node < 3; ++num_node) {
					node[0] = base_tetra_geo.inclined_face.row(num_node)[0];
					node[1] = base_tetra_geo.inclined_face.row(num_node)[1];
					node[2] = 0;
					FromPlaneToTetra(base_tetra_geo.inverse_transform_matrix, base_tetra_geo.start_point_plane_coord, node, node);
					FromLocalToGlobalTetra(vertex_tetra, node, x.x[num_cur_out_face][num_node]);

				}
				break;
			default:
				RETURN_ERR("Number face is not {0,1,2,3}????\n");
			}
		}

		vec_x[num_cell] = x;
	}

	return 0;
}

static Vector2 CalculateIllumeOnInnerFace(const int num_in_face, const int neib_id, const Matrix4& vertex_tetra,
	const Vector3& x, const Vector3& x0,  std::vector<Type>& vec_res_bound) 
{
	Vector2 x0_local(-100, -100);

	Type I_x0 = 0;
	switch (neib_id)
	{
		case eBound_FreeBound:
			return x0_local;

		case eBound_LockBound:
			return x0_local;

		case eBound_OutSource:
			return x0_local;

		case eBound_InnerSource:
		{
#if 0
			I_x0 = 10;
			vec_res_bound.push_back(I_x0);			
			posRes++;
			return x0_local;

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
					vec_res_bound.push_back(I_x0);
					posRes++;
					return x0_local; // 2;
				}

				// с чем луч встречается раньше?
				{
					const Vector3 Xsphere = x - inSphere;
					const Vector3 Xres = x - res;
					const Type LenSpehere = Xsphere.dot(Xsphere);
					const Type LenPlane = Xres.dot(Xres);

					if (LenSpehere > LenPlane) {
						I_x0 = 20;						
						posRes++;
						vec_res_bound.push_back(I_x0);
						return x0_local; // 1;
					}
					else {
						I_x0 = 50;						
						vec_res_bound.push_back(I_x0);
						posRes++;
						return x0_local; // 2;
					}
				}
			}
			else if ((dist < R2disk * R2disk) && (dist > R1disk * R1disk))
			{
				I_x0 = 20;				
				vec_res_bound.push_back(I_x0);
				posRes++;
				return x0_local; // 1
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
					x0_local = (Vector2(local_plane_x0[0], local_plane_x0[1]));					
					break;

				case 1:
					x0_local = (Vector2(try_x0[1], try_x0[2]));										
					break;
				case 2:
					x0_local = (Vector2(try_x0[0], try_x0[2]));										
					break;
				case 0:
					x0_local = (Vector2(try_x0[0], try_x0[1]));										
					break;
				default:
					EXIT_ERR("Error\n");					
				}				
				I_x0 = -10;
				vec_res_bound.push_back(I_x0);				
				posRes++;
			}
#endif							
			return x0_local; // 2;
		}
		
		default:
			Vector3 x0_loc;
			FromGlobalToLocalTetra(vertex_tetra, x0, x0_loc);

			Vector3 local_plane_x0;
			switch (num_in_face)
			{
			case 3:				
				FromTetraToPlane(base_tetra_geo.transform_matrix, base_tetra_geo.start_point_plane_coord, x0_loc, local_plane_x0);
				return Vector2(local_plane_x0[0], local_plane_x0[1]);
			case 1:
				return Vector2(x0_loc[1], x0_loc[2]);								
			case 2:
				return Vector2(x0_loc[0], x0_loc[2]);
			case 0:
				return Vector2(x0_loc[0], x0_loc[1]);
			default:
				EXIT_ERRS("Error num_in_face = %d", num_in_face);
			}		
	}		
	return Vector2(0, 0);
}

static int CalculateNodeValue(const int num_cell, const Normals& normal, const std::vector<Face>& grid,
 const int face_state, const Vector3& direction, const Matrix4& vertex_tetra, const Vector3& x,
	const std::vector<int>& all_pairs_face,  std::vector<Type>& vec_res_bound, std::vector<cell_local>& vec_x0)
{
	Vector3 x0;

	for (ShortId num_in_face = 0; num_in_face < base; ++num_in_face) 
	{
		const int face_id = num_cell * base + num_in_face;
		if (!check_bit(face_state,num_in_face)) continue;  // обрабатываем только входные грани

		IntersectionWithPlane(grid[face_id], x, direction, x0);

		if (InTriangle(num_in_face, grid[face_id], normal, x0))
		{			
			cell_local x0_local;
			x0_local.s = (x - x0).norm();
			x0_local.in_face_id = num_in_face;			
			x0_local.x0 = CalculateIllumeOnInnerFace(num_in_face, all_pairs_face[face_id], vertex_tetra, x, x0, vec_res_bound);

			vec_x0.push_back(x0_local);
			break;
		}

	}//for num_in_face

	return 0;
}

int GetNodes(const int num_cell, const std::vector<Face>& grid, const ShortId num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int face_state, const Vector3& direction, const std::vector<Normals>& normals,
	const std::vector<int>& all_pairs_face, const BasePointTetra& x,
	std::vector<Type>& vec_res_bound, std::vector<cell_local>& vec_x0) 
{

	switch (num_cur_out_face)
	{

	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) 
		{	
			CalculateNodeValue(num_cell, normals[num_cell], grid, face_state, direction, vertex_tetra, x(num_cur_out_face, num_node), 
				all_pairs_face, vec_res_bound, vec_x0);
		}
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			CalculateNodeValue(num_cell, normals[num_cell], grid, face_state, direction, vertex_tetra, x(num_cur_out_face, num_node),
				all_pairs_face, vec_res_bound, vec_x0);
		}
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {

			CalculateNodeValue(num_cell, normals[num_cell], grid, face_state, direction, vertex_tetra, x(num_cur_out_face, num_node),
				all_pairs_face, vec_res_bound, vec_x0);
		}// x->координата узла на выходящей грани		}	
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {

			CalculateNodeValue(num_cell, normals[num_cell], grid, face_state, direction, vertex_tetra, x(num_cur_out_face, num_node),
				all_pairs_face, vec_res_bound, vec_x0);
		}
		break;
	default:
		RETURN_ERR("Number face is not {0,1,2,3}????\n");				
	}

	return 0;
}
#endif // ONLY_GEO_DATA
#endif //MAKE
