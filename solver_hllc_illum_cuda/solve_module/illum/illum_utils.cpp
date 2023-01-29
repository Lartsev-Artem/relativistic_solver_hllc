#include "illum_utils.h"
#if defined ILLUM

#include "../../file_module/reader_bin.h"

/*! Пересчет излучения с грани в ячейку в заданном направлении
num_dir - номер направления
illum_on_face - излучение на гранях по всем направлениям
*/
int GetDirectionIllumFromFace(const int size_grid, const int num_dir, const std::vector<Type>& illum_on_face, std::vector<Type>& illum_in_cell)
{	
	if (illum_on_face.size() < size_grid* base +size_grid * num_dir * base) RETURN_ERR("illum_on_face hasn't enough data\n");

	illum_in_cell.resize(size_grid, 0);	
	for (size_t j = 0; j < size_grid; j++)
	{
		const int N = num_dir * base * size_grid + j * base;
		for (size_t k = 0; k < base; k++)
		{
			illum_in_cell[j] += illum_on_face[N + k];
		}
		illum_in_cell[j] /=  base;		
	}
	return 0;
}
#endif
#if defined ILLUM && defined SOLVE
int SolveIllumAndHLLC(const Type tau, std::vector<elem_t>& cells)
{
	for (auto& el : cells)
	{
		//el.phys_val.d += 0;
		//el.phys_val.v += (C_LIGHT*DIST/TIME/PRESSURE) * tau * (-el.illum_val.div_impuls * C_LIGHT_INV);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
		//el.phys_val.p += (1./TIME/RADIATION)*tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e

		el.conv_val.v +=   tau * (-el.illum_val.div_impuls * 1);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
		el.conv_val.p +=   tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e
	}	
	return 0;
}

int InitIllum(file_name main_dir, grid_t& grid)
{
	double _clock = -omp_get_wtime();

	if (ReadDataArray(solve_mode, main_dir, grid))
	{
		RETURN_ERR("Error reading the data array vtk\n");
	}
	_clock += omp_get_wtime();

	WRITE_LOG("Reading the data array " << _clock << " c\n");	

#if 0 // инициализация составляющий для излучения (не надо для квазилинейного случая)
	CalculateIllumOptMemory(main_dir, use_cuda, cuda_mod, accuracy, max_number_of_iter, directions, Illum, int_scattering, squares,
		grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
		res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry, U_full_prev);

	CalculateIllumParamOptMemory(size_grid, count_directions,
		directions, squares, normals, squares_cell, volume,
		Illum, energy, stream, impuls, div_stream, div_impuls);

	printf("Init Illum completed\n\n");
#endif

	return 0;
}


Type BoundaryConditions(const int type_bound, Vector3& inter_coef)
{
	Type I0 = 0;
	switch (type_bound)
	{
	case eBound_OutSource: // дно конуса
		I0 = 30;
		break;
	case eBound_FreeBound:
		I0 = 0;
		break;

	case eBound_LockBound:
		I0 = 0;
		break;

	case eBound_InnerSource:  // внутренняя граница	
	{
#if 0
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
#endif
		I0 = 30;
		break;
	}
	default:

		WRITE_LOG("unknown bound type in illum\n");
		I0 = 0;
		break;
	}

	inter_coef = Vector3(I0, I0, I0);
	return I0;
}

#endif //ILLUM
