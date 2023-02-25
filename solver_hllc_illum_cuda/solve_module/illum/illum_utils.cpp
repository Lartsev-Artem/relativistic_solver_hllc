#include "illum_utils.h"
#if defined ILLUM

#include "../../file_module/reader_bin.h"

/*! Пересчет излучения с грани в ячейку в заданном направлении
num_dir - номер направления
illum_on_face - излучение на гранях по всем направлениям
*/

int GetDirectionIllumFromFace(const int size_grid, const int num_dir, const Type* illum_on_face, 
	std::vector<Type>& illum_in_cell)
{	
	//if (illum_on_face.size() < size_grid* base +size_grid * num_dir * base) RETURN_ERR("illum_on_face hasn't enough data\n");
	if (illum_on_face == nullptr) RETURN_ERR("illum_on_face hasn't enough data\n");

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

#ifdef USE_CUDA
#include "../../cuda/cuda_solve.h"
void CalculateParamOnCuda(const grid_directions_t& grid_dir, grid_t& grid)
{
	/* //это для нестационарного излучения (todo: не считать stream и impuls отдельно
	CalculateEnergy(32, grid.size, grid_direction.size, energy);
	CalculateStream(32, grid.size, grid_direction.size, stream);
	CalculateImpuls(32, grid.size, grid_direction.size, impuls);
	*/
									
	CalculateEnergy(grid_dir, grid);
	CalculateDivStream(grid_dir, grid);
	CalculateDivImpuls(grid_dir, grid);

	return;
}
#endif

static int rhllc_get_phys_value3(const flux_t& U, flux_t& W)
{
	Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));
	Type p = W.p;
	
	const Type h = 1 + gamma_g * p / W.d;

	Type W0 = W.d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m = U.v;
	Type mm = m.dot(m);


	Vector3 v = W.v;

	Type D = U.d;
	Type E = U.p;

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		//W.d * h * Gamma0 * Gamma0 - p - E;
		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
		W0 -= (fW / dFdW);

		Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v = m / W0;

		err -= W0;
		cc++;

		if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
		{
			p = sqrt(mm) - E;
			v = m / (E + p);
			Gamma0 = 1. / sqrt(1 - v.dot(v));
			//printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
			break;
		}

	} while (fabs(err / W0) > 1e-14);

	if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
	{
		printf("W= %lf,(%lf), %lf", W.d, W.v.norm(), W.p);
		printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
		EXIT(1);
	}

	W.d = D / Gamma0;
	W.v = v;
	W.p = p;

	return 0;
}

static int rhllc_get_phys_value_debug2(const flux_t& U, flux_t& W)
{
	Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));	
	Type p = W.p;
	const Type h = 1 + gamma_g * p / W.d;

	Type W0 = W.d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m = U.v;
	Type mm = m.dot(m);


	Vector3 v = W.v;

	Type D = U.d;
	Type E = U.p;

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		//W.d * h * Gamma0 * Gamma0 - p - E;
		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
		W0 -= (fW / dFdW);

		//if(mm < W0*W0 )
			Gamma0 = /*1. /*/ sqrt(1 - mm / (W0 * W0));
		//else
		//{
		//	Gamma0 = 1. / 1e-50;
		//	const Type h = 1 + gamma_g * p / (D/Gamma0);
		//	Type W0 = D * h * Gamma0; //U[0] * Gamma0 * h;
		//}


		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v = m / W0;

		err -= W0;
		cc++;

		if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
		{
			printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
			//	EXIT(1);
		}

	} while (fabs(err / W0) > 1e-14);


	W.d = D / Gamma0;
	W.v = v;
	W.p = p;

	return 0;
}
static int rhllc_get_phys_value_debug(const flux_t& U, flux_t& W)
{
	Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));
	Type p = W.p;
	const Type h = 1 + gamma_g * p / W.d;

	Type W0 = W.d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m = U.v;
	Type mm = m.dot(m);


	Vector3 v = W.v;

	Type D = U.d;
	Type E = U.p;

	W0 = (2 * E + sqrt(4 * E * E - 3 * mm) / 3 - D); // условие обеспечивающее положительность давления при условии физически правилных консервативных переенных

	p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

	Gamma0 = /*1. /*/ sqrt(1 - mm / (W0 * W0));
	p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		//W.d * h * Gamma0 * Gamma0 - p - E;
		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));

		Type W2 = W0 * W0;
		dFdW = (D * (1 - gamma1) * mm * sqrt(1 - mm / W2) * W0 - 3 * mm * W2 + W2 * W2 + gamma1 * mm * (mm + W2)) / (gamma1 * (mm - W2) * (mm - W2));

		W0 -= (fW / dFdW);

		Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v = m / W0;

		err -= W0;
		cc++;

		if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
		{
			printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
			//	EXIT(1);
		}

	} while (fabs(err / W0) > 1e-14);

	if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
	{
		printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
		//EXIT(1);
	}

	W.d = D / Gamma0;
	W.v = v;
	W.p = p;

	return 0;
}
int rhllc_get_conv_value(const flux_t& W, flux_t& U)
{
	const Type d = W.d;
	const Type Gamma = 1. / sqrt(1 - W.v.dot(W.v));
	const Type h = 1 + gamma_g * W.p / d;
	const Type dhGG = d * h * Gamma * Gamma;

	U.d = Gamma * W.d;
	U.v = dhGG * W.v;
	U.p = dhGG - W.p;

	return 0;
}

int SolveIllumAndHLLC(const Type tau, grid_t& grid)
{
#if 0
	Type min = 1;
	int num = 0;
	int ret_flag = 0;
	for (auto& el : cells)
	{
		if (el.conv_val.p - sqrt(el.conv_val.d * el.conv_val.d + el.conv_val.v.dot(el.conv_val.v)) < 0)
		{
			ret_flag = 1;// return 1;
		}
		num++;
		//if (el.conv_val.p < tau * el.illum_val.div_stream) //слишком большое приращение
		//{		
		//}
		//else
		{
			//el.phys_val.d += 0;
			//el.phys_val.v += (C_LIGHT*DIST/TIME/PRESSURE) * tau * (-el.illum_val.div_impuls * C_LIGHT_INV);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
			//el.phys_val.p += (1./TIME/RADIATION)*tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e

			el.conv_val.v += tau * (-el.illum_val.div_impuls * 1);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
			el.conv_val.p += tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e			
#if 0
			rhllc_get_phys_value3(el.conv_val, el.phys_val);

			if (el.conv_val.p - sqrt(el.conv_val.d * el.conv_val.d + el.conv_val.v.dot(el.conv_val.v)) < 0)
			{

				el.phys_val.p = 1e-10; //E=W-p
				/*Type W0 = el.conv_val.p - el.phys_val.p;
				Type Gamma0 = 1. / sqrt(1 - el.conv_val.v.dot(el.conv_val.v) / (W0 * W0));				
				el.phys_val.v = el.conv_val.v / W0;
				el.phys_val.d = el.conv_val.d / Gamma0;*/

				// что со скоростью?
				rhllc_get_conv_value(el.phys_val, el.conv_val);
				printf("Err cond %d, %d\n", (&el - &cells[0])/sizeof(elem_t), num);

			}
#endif
			//rhllc_get_phys_value_debug(el.conv_val, el.phys_val);
		}				
	}	
#endif
	Type min = 1;	
	int ret_flag = 0;

#pragma omp parallel  default(none) shared(ret_flag, grid, solve_mode) 
	{
		const int n = grid.size;
#pragma omp for
		for (int i = 0; i < n; i++)
		{
			elem_t &el = grid.cells[i];

#ifdef USE_CUDA
			el.conv_val.v += tau * (-grid.divimpuls[i] * 1);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
			el.conv_val.p += tau * (-grid.divstream[i]); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e			

#else
			el.conv_val.v += tau * (-el.illum_val.div_impuls * 1);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
			el.conv_val.p += tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e			
#endif
			
			if (el.conv_val.p - sqrt(el.conv_val.d * el.conv_val.d + el.conv_val.v.dot(el.conv_val.v)) < 0)
			{
				ret_flag = 1;// return 1;
			}
		}
	}
	return ret_flag;
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
		I0 = 1e4;
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
		I0 = 100;// 1e5;
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
