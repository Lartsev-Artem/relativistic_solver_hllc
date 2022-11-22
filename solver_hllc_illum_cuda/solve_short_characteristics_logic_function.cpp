#include "solve_short_characteristics_logic_function.h"

Type CalculateIllumeOnInnerFace(const int num_cell, const int num_direction, const int num_in_face ,const Vector3& x, 
	const std::vector<Vector2>&X0, 
	 std::vector<cell>& grid, const std::vector<int>& neighbours_id_face,
	int& id_try_pos, int& pos_in_res, int& posX0, 
	const uint64_t ShiftRes, const uint64_t ShiftX0, const int ShiftTry) {
	Type I_x0 = 0;
	
	if (neighbours_id_face[num_cell * 4 + num_in_face] ==  eBound_OutSource) // дно конуса
	{
		grid[num_cell].nodes_value[num_in_face] = Vector3(30, 30, 30);
		return 30;
	}
	else if (neighbours_id_face[num_cell * 4 +  num_in_face] == eBound_FreeBound ||
		neighbours_id_face[num_cell * 4 + num_in_face] == eBound_LockBound )
	{
		grid[num_cell].nodes_value[num_in_face] = Vector3(0, 0, 0);
		/*Граничные условия*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;
	}
	else if (neighbours_id_face[num_cell * 4 + num_in_face] == eBound_InnerSource) // внутренняя граница
	{		
		id_try_pos++;
		grid[num_cell].nodes_value[num_in_face] = Vector3(res_on_inner_bound, res_on_inner_bound, res_on_inner_bound);
		return res_on_inner_bound;

		Type data = res_inner_bound[ShiftRes + pos_in_res++]; // защита на выход из диапазона??
		if (data >= 0) //данные от пересечения с диском или шаром
		{
		/*
			 Проверить порядок. Данные в массиве должны лежать в соответствии с упорядоченным графом 
			 от направления к направлению
		*/
			return res_on_inner_bound; // I_x0;
			return data;  
		}
		else // определяющимм являются противолежаащие грани (возможен расчет с учетом s=(x-x0).norm())
		{

			// результат не вполне понятен. Пока лучше использовать константу или другие параметры области (шар == граница)

			//+dist_try_surface			
			int id = id_try_surface[ShiftTry + id_try_pos - 1];  // будет лежать id грани			
			const int cell = id / 4;
			const int face = id % 4;		
			
			Vector3 coef = grid[cell].nodes_value[face];
			Vector2	x0_local = X0[ShiftX0 + posX0++];//grid[num_cell].x0_loc[num_in_face_dir];

			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
			
			//if (I_x0 < 0) I_x0 = 0;
			return  res_on_inner_bound; // I_x0;

		}
	}
	else 
	{

		Vector3 coef = grid[num_cell].nodes_value[num_in_face];
	
		//сейчас храним значения а не коэффициента интерполяции
				
		//Vector2	x0_local = X0[ShiftX0 + posX0++]; // grid[num_cell].x0_loc[num_in_face_dir];
		//I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
		
		I_x0 = (coef[0] + coef[1] + coef[2]) / 3;

		if (I_x0 < 0) 
		{			
			return 0;
		}

		return I_x0;
	}
}


Type CurGetIllum(const int cur_id, const int cur_direction, const Vector3 x, const Type s, const Type I_node_prev,
	const vector<Type>& int_scattering, const vector<VectorX>& U) {
	
	
	if (class_file_vtk == 10) // для конуса (считаем, что излучающая часть не изменяется в зависимости от газораспределения)
	{
		Type Q = 0;
		Type alpha = 0.5;
		Type betta = 0.5;
		Type S = int_scattering[cur_direction * size_grid + cur_id];

		
		//if (x[0] < 0.06) // излучающий слой
		//{
		//	Q = 10; alpha = 1;  betta = 2; 
		//}

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_node_prev + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	// HLLC + Illum
	if (class_file_vtk == 11) // для конуса
	{

		Type S = int_scattering[cur_direction * size_grid + cur_id];
		Type Ie = 1;		
		Type alpha = 0.5;
		Type betta = 0.5;

		{
			const Type R = 8.314;
			const Type c = 3 * 1e8;
			const Type h = 6.62*1e-34;
			const Type k = 1.38*1e-23;
			const Type sigma = 6.652 * 1e-29;
			const Type m = 1.6735575 * 1e-27;

			Type d = U[cur_id](0);
			Vector3 vel(
			U[cur_id](1) / d,
			U[cur_id](2) / d,
			U[cur_id](3) / d);
			 Type v = vel.norm();
			 Type p = (U[cur_id](4) - v * v * d / 2.) * (gamma1 - 1);
			 Type T = p / (d * R);

			if (x[0] < 0.05) 
			{
				d = 0.1;
				p = 0.01;
				T = p / (d * R);		
			}

			Ie = 2 * pow(k * PI * T, 4) / (15 * h * h * h * c * c);
			betta = sigma * d / m;
			alpha = betta;  			//alpha = ???			
		}
		
		Type Q = alpha * Ie;				

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_node_prev + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}


	// без интеграла рассеивания  (излучающий шар)
	if (class_file_vtk == 0)
	{
		Type Q = 0;
		Type alpha = 2;
		Type betta = 1;
		Type S = 0;// int_scattering[cur_direction * size_grid + cur_id];

		if ((x - Vector3(1, 0, 0)).norm() > 0.09) { Q = 0; alpha = 0.5;  betta = 0.5; }

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_node_prev + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;

	//-------------------------------------------
			/*Type Ie = 10;
			Type k = 10;
			if ((x - Vector3(1, 0, 0)).norm() > 0.09) { Ie = 0; k = 1; }

			Type I;

			if (k > 1e-10)
				I = Ie * (1 - exp(-s * k)) + I_node_prev * exp(-s * k);
			else
				I = I_node_prev * (1 - s * k) + Ie * s * k;

			if (I < 0)
				I = 0;
			return I;*/
	}

	// test task
	if (class_file_vtk == 1)
	{

#ifdef USE_VTK
		Type S = int_scattering[cur_direction * size_grid + cur_id]; //GetS(cur_id, cur_direction, illum_old, directions, squares);
		Type Q = 0; // rad_en_loose_rate->GetTuple1(cur_id);  //Q=alpha*Ie
		Type alpha = 0; // absorp_coef->GetTuple1(cur_id);
#else
		//U_Full[cur_id][0] -> density;
		Type S = int_scattering[cur_direction * size_grid + cur_id]; //GetS(cur_id, cur_direction, illum_old, directions, squares);
		Type Q = rad_en_loose_rate[cur_id];  //Q=alpha*Ie
		Type alpha = absorp_coef[cur_id];
#endif


		Type betta = alpha / 2;  // просто из головы
		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_node_prev + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	if (class_file_vtk == 2) //
	{
		const Type ss = s * 388189 * 1e5;  // числа --- переход к размерным параметрам
#ifndef USE_VTK
		Type Q = rad_en_loose_rate[cur_id];  //Q=alpha*Ie
		Type alpha = density[cur_id] * absorp_coef[cur_id];
#else
		Type Q = 0;
		Type alpha = 1;
#endif


		Type I;
		if (alpha > 1e-15)
			I = I_node_prev * exp(-alpha * ss) + Q * (1 - exp(-alpha * ss)) / alpha;
		else
			I = I_node_prev * (1 - alpha * ss) + Q * ss;

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	if (class_file_vtk == 3)
	{

		const Type d = U[cur_id][0];
		Type S = int_scattering[cur_direction * size_grid + cur_id]; //GetS(cur_id, cur_direction, illum_old, directions, squares);

		//alpha = d*XXX...
		//betta = d*sigma*XXX...
		Type alpha = 0;// absorp_coef[cur_id];
		Type betta = alpha / 2;  // просто из головы
		Type Q = 0;//  rad_en_loose_rate[cur_id];  //Q=alpha*Ie

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_node_prev * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_node_prev + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}
}

#ifdef  USE_VTK
Type GetValueInCenterCell(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const Vector3 center,
	const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::vector<cell>& nodes_value,
	const std::vector<Type>& illum_old, const vector<Vector3>& directions, const vector<Type>& squares) {
	/*Все грани должно быть определены*/
	Type value = -666;
	Vector3 x0;

	for (size_t i = 0; i < 4; i++) {

		IntersectionWithPlane(cur_cell->GetFace(i), center, direction, x0);
		if (InTriangle(num_cell, unstructuredgrid, cur_cell, i, x0)) {
			if ((center - x0).dot(direction) <= 0) continue;
			
			Type s = (center - x0).norm();
		//	Type I_x0 = CalculateIllumeOnInnerFace(num_cell, i, vertex_tetra, center, x0, nodes_value);

			//value = CurGetIllum(num_cell, x0, s, I_x0, direction, illum_old, directions, squares);
			break;
		}
	}
	if (value < 0)
		return 0;
	return value;
}
#endif //  USE_VTK