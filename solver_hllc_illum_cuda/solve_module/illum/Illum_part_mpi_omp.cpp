#include "../solve_config.h"
#if defined ILLUM && defined SOLVE && defined USE_MPI

#include "../../file_module/reader_bin.h"
#include "../../file_module/reader_txt.h"
#include "../../file_module/writer_bin.h"

#include "illum_utils.h"

#include "../solve_global_struct.h"
#include "../../global_value.h"
#include "../../global_def.h"
#include "../solve_utils.h"

#include "../../cuda/cuda_solve.h"

static std::vector<int> send_count;
static std::vector<int> disp;

static std::vector<int> send_count_illum;
static std::vector<int> disp_illum;

static std::vector<int> send_count_scattering;
static std::vector<int> disp_scattering;

static std::vector<flux_t> phys_local; 
//static std::vector<Vector3> inter_coef; 
static std::vector<std::vector<Vector3>> inter_coef_all;
static std::vector<Type> loc_illum; 	
static std::vector<Type>int_scattering_local;

int InitSendDispIllumArray(const int myid, const int np, const int count_directions, const int count_cells)
{
	GetSend(np, count_directions, send_count);
	GetDisp(np, count_directions, disp);

	send_count_illum.resize(np);
	send_count_scattering.resize(np);

	disp_illum.resize(np);
	disp_scattering.resize(np);

	for (int i = 0; i < np; i++)
	{
		send_count_illum[i] = send_count[i] * base * count_cells;
		send_count_scattering[i] = send_count[i] * count_cells;

		disp_illum[i] = disp[i] * base * count_cells;
		disp_scattering[i] = disp[i] * count_cells;
	}

	{
		int len[3 + 1] = { 1,3,1,  1 };
		MPI_Aint pos[4] = { offsetof(flux_t, d), offsetof(flux_t, v),offsetof(flux_t, p) ,sizeof(flux_t) };
		MPI_Datatype typ[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(4, len, pos, typ, &MPI_flux_t);
		MPI_Type_commit(&MPI_flux_t);
	}

	const int local_size = send_count[myid];

	phys_local.resize(count_cells);

	inter_coef_all.resize(omp_get_max_threads());
	for (int i = 0; i < inter_coef_all.size(); i++)
	{
		inter_coef_all[i].resize(count_cells * base);
	}

	WRITE_LOG("Omp threads in illum= " << inter_coef_all.size() << "\n");
	//inter_coef.resize(count_cells * base, Vector3::Zero());

	loc_illum.resize(count_cells * base * local_size, 0); //кроме id 0	
	int_scattering_local.resize(count_cells * local_size, 0);

	return 0;
}

static Type CalculateIllumeOnInnerFace(const int neigh_id, Vector3& inter_coef)
{		
	Type I_x0 = 0;

	if (neigh_id < 0)
	{
		I_x0 = BoundaryConditions(neigh_id, inter_coef);
	}
	else	
	{
		Vector3 coef = inter_coef;

		//сейчас храним значения а не коэффициента интерполяции
		I_x0 = (coef[0] + coef[1] + coef[2]) / 3;

		if (I_x0 < 0)
		{
			WRITE_LOG("Error illum value\n");
			return 0;
		}			
	}

	return I_x0;
}

static Type GetS(const int num_cell, const Vector3& direction, const grid_t&grid,
	const grid_directions_t& grid_direction) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
		Type dot = direction.dot(direction2);
	return (3. * (1 + dot* dot)) / 4;
	} };

	Type S = 0;
	const int N_dir = grid_direction.size;
	const int N_cell = grid.size * base; // / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		Type I = 0;
		for (int i = 0; i < base; i++)
		{
			I+= grid.Illum[num_direction * N_cell + num_cell + i];
		}
		I /= base;
		
		S += Gamma(grid_direction.directions[num_direction].dir, direction) * I * grid_direction.directions[num_direction].area;
	}
	return S / grid_direction.full_area;     // было *4PI, но из-за нормировки Gamma разделили на 4PI
}

static int CalculateIntCPU(const int num_cells, const grid_directions_t& grid_direction, grid_t& grid)
{

#pragma omp parallel default(none) shared(num_cells, grid_direction, grid)
	{
		Vector3 direction;
		const int num_directions = grid_direction.size;
#pragma omp for
		for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
			direction = grid_direction.directions[num_direction].dir;
			for (int cell = 0; cell < num_cells; cell++)
			{
				grid.scattering[num_direction * num_cells + cell] = GetS(base * cell, direction, grid, grid_direction);
			}
		}
	}
	return 0;
}


static Type GetCurIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, const flux_t& cell)
{
	switch (solve_mode.class_vtk)
	{
	case 0: // без интеграла рассеивания  (излучающий шар)
	{
		Type Q = 0;
		Type alpha = 2;
		Type betta = 1;
		Type S = int_scattering;// 0;

		//if ((x - Vector3(1, 0, 0)).norm() > 0.09) { Q = 0; alpha = 0.5;  betta = 0.5; }

		if ((x - Vector3(0.5, 0, 0)).norm() < 0.05) { Q = 0; alpha = 0.5;  betta = 0.5; }

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_0 * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_0 + s * (Q + S * betta));

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

	case 1: // test task
	{

		//U_Full[cur_id][0] -> density;
		Type S = int_scattering;
		Type Q = 1;// cell.illum_val.rad_en_loose_rate;  //Q=alpha*Ie
		Type alpha = 1;// cell.illum_val.absorp_coef;

		Type betta = alpha / 2;  // просто из головы
		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_0 * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_0 + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	case 2:
	{
		const Type ss = s * 388189 * 1e5;  // числа --- переход к размерным параметрам

		//Type Q = cell.illum_val.rad_en_loose_rate;  //Q=alpha*Ie
		//Type alpha = cell.d * cell.illum_val.absorp_coef;
		Type alpha = 1;
		Type Q = 1;

		Type I;
		if (alpha > 1e-15)
			I = I_0 * exp(-alpha * ss) + Q * (1 - exp(-alpha * ss)) / alpha;
		else
			I = I_0 * (1 - alpha * ss) + Q * ss;

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	case 3:
	{

		const Type d = cell.d;
		Type S = int_scattering;

		//alpha = d*XXX...
		//betta = d*sigma*XXX...
		Type alpha = 0;// absorp_coef[cur_id];
		Type betta = alpha / 2;  // просто из головы
		Type Q = 0;//  rad_en_loose_rate[cur_id];  //Q=alpha*Ie

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_0 * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_0 + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	case 10: // для конуса (считаем, что излучающая часть не изменяется в зависимости от газораспределения)
	{
		Type Q = 0;
		Type alpha = 0.5;
		Type betta = 0.5;
		Type S = int_scattering;


		//if (x[0] < 0.06) // излучающий слой
		//{
		//	Q = 10; alpha = 1;  betta = 2; 
		//}

		Type k = alpha + betta;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * s) * (I_0 * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - s * k) * (I_0 + s * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		return I;
	}

	case 11: // HLLC + Illum для конуса
	{

		Type S = int_scattering * RADIATION;

		Type d = cell.d * DENSITY;		
		Type p = cell.p * PRESSURE;

		Type T = 1e-5 * p / (d * R_gas);

		Type 	Ie = 2 * pow(k_boltzmann * PI * T, 4) / (15 * h_plank * h_plank * h_plank * c_light * c_light);
		Type 	betta = sigma_thomson * d / m_hydrogen;
		Type	alpha = betta; 			//alpha = ???			


		Type Q = alpha * Ie;

		Type k = alpha + betta;

		Type ss = s * DIST;
		Type I0 = I_0 * RADIATION;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * ss) * (I0 * k + (exp(k * ss) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - ss * k) * (I0 + ss * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

		/*static Type max = 0;
		if (I / RADIATION > max) {
			max = I / RADIATION;
			printf("I=%lf\n", I / RADIATION);
		}*/

		return  I / RADIATION;
	}


	default:
		EXIT_ERR("unknow class_vtk in get illum\n");
	}
}

//static Type ReCalcIllum(const int num_dir, std::vector<elem_t>& cells, std::vector<Type>& Illum)
static Type ReCalcIllum(const int num_dir, const std::vector<Vector3>& inter_coef, std::vector<Type>& Illum)
{	
	Type norm = -1;
//#pragma omp parallel default(none) shared(num_dir, cells, Illum, norm)
	{		
		Type norm_loc = -1;
		const int size = inter_coef.size()/base;
		const int shift_dir = num_dir * size;
	
//#pragma omp for
		for (int num_cell = 0; num_cell < size; num_cell++)
		{
			for (int i = 0; i < base; i++)
			{
				Vector3 Il = inter_coef[num_cell * base + i]; //cells[num_cell].illum_val.coef_inter[i];
				const Type curI = (Il[0] + Il[1] + Il[2]) / 3;
				const int id = base * (shift_dir + num_cell) + i;
				{					
					Type buf_norm = fabs(Illum[id] - curI) / curI;
					//if (curI < 1e-250) buf_norm = 1;
					if (buf_norm > norm_loc) norm_loc = buf_norm;
				}
				Illum[id] = curI;	
			}
			//cells[num_cell].illum_val.illum[num_dir] /= base; //пересчёт по направлению для расчёта энергий и т.д.
		}

//#pragma omp critical
		{
			if (norm_loc > norm)
			{
				norm = norm_loc;
			}		
		}
	} //parallel

	return norm;
}

#ifndef USE_CUDA
static Type ReCalcIllumGlobal(const int dir_size, grid_t& grid)
{
	
//#pragma omp parallel default(none) shared(dir_size, cells,Illum)
	{
		const int size = grid.size;
//#pragma omp for
		for (int num_dir = 0; num_dir < dir_size; num_dir++)
		{
			const int shift_dir = num_dir * size;
			for (int num_cell = 0; num_cell < size; num_cell++)
			{
				for (int i = 0; i < base; i++)
				{
					const int id = base * (shift_dir + num_cell) + i;
					grid.cells[num_cell].illum_val.illum[num_dir * base + i] = grid.Illum[id]; //на каждой грани по направлениям
				}
			}
		}
	}
	return 0;
}
#endif

//#include "../../file_module/writer_bin.h"
static int GetIntScattering(const int count_cells, const grid_directions_t& grid_direction, grid_t& grid)
{
	if (solve_mode.max_number_of_iter > 1)  // пропуск первой итерации
	{
#ifdef USE_CUDA		
		if (CalculateIntScattering(grid_direction, grid))  // если была ошибка пересчитать на CPU 
		{
			ClearDevice();
			EXIT_ERR("");
		}
#else
		CalculateIntCPU(count_cells, grid_direction, grid);
#endif
	}
	else
	{
#ifdef USE_CUDA
		CopyIllumOnDevice(count_cells * base * grid_direction.size, grid.Illum);
#endif
	}
		
	return 0;
}

int MPI_CalculateIllum(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell,
	//const std::vector<Type>& res_inner_bound, 
	grid_t& grid)
{	
	const int count_directions = grid_direction.size;	
	const int count_cells = grid.size;
	
	int np = 1, myid = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
		
	const int local_size = send_count[myid];
	const int local_disp = disp[myid];
		
	int count = 0;
	Type norm = -1;
	std::vector<Type> norms;
	
	if (myid == 0)
	{
		int i = 0;
		for (auto& el : grid.cells)
		{
			phys_local[i++] = el.phys_val;
			// сюда же absorp_coef, если понадобиться
		}

		norms.resize(np, -10);
	}
	MPI_Bcast(phys_local.data(), count_cells , MPI_flux_t, 0, MPI_COMM_WORLD);
		
	do {
		
		Type _clock = -omp_get_wtime();
		norm = -1;		


#pragma omp parallel default(none) shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
		norm, inter_coef_all,local_size, int_scattering_local, phys_local, loc_illum, disp, myid, np)
		{
			Type loc_norm = -1;
			const int num = omp_get_thread_num();
			std::vector<Vector3>* inter_coef = &inter_coef_all[num];

			/*---------------------------------- далее FOR по направлениям----------------------------------*/
#pragma omp for
			for (register int num_direction = 0; num_direction < local_size; ++num_direction)
			{
				int posX0 = 0;
				Vector3 I;

				/*---------------------------------- далее FOR по ячейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h)
				{

					const int num_cell = sorted_id_cell[num_direction][h];
					const int face_block_id = num_cell * base;

					for (ShortId num_out_face = 0; num_out_face < base; ++num_out_face)
					{
						const int out_id_face = pairs[face_block_id + num_out_face];
						if (check_bit(face_states[num_direction][num_cell], (int)num_out_face))
						{
							//формально обрабатываем мы только выходящие грани, т.к. входящие к этому моменту определены,
							// но границу надо инициализировать. Не для расчета, но для конечного интегрирования
							if (out_id_face < 0)
							{
								BoundaryConditions(out_id_face, (*inter_coef)[face_block_id + num_out_face]);
							}
							continue;
						}

						//GetNodes
						for (int num_node = 0; num_node < 3; ++num_node)
						{
							const Vector3 x = vec_x[num_cell].x[num_out_face][num_node];

							const cell_local x0 = vec_x0[num_direction][posX0++];

							const ShortId num_in_face = x0.in_face_id;
							const Type s = x0.s;
							const Vector2 X0 = x0.x0;

							const Type I_x0 = CalculateIllumeOnInnerFace(pairs[face_block_id + num_in_face], (*inter_coef)[face_block_id + num_in_face]);

							I[num_node] = GetCurIllum(x, s, I_x0, int_scattering_local[num_direction * count_cells + num_cell], phys_local[num_cell]);

						}//num_node

						(*inter_coef)[face_block_id + num_out_face] = I;
						if (out_id_face >= 0)
						{
							(*inter_coef)[out_id_face] = I;
						}

						//сюда сразу в loc_illum. сходимость ез нормы задать числом итераций

					} //num_out_face	

				}

				/*---------------------------------- конец FOR по ячейкам----------------------------------*/
				loc_norm = ReCalcIllum(num_direction, *inter_coef, loc_illum);
			}

			/*---------------------------------- конец FOR по направлениям----------------------------------*/
			
			if (loc_norm > norm)
			{
#pragma omp critical
				{
					if (loc_norm > norm)
					{
						norm = loc_norm;
					}
				}
			}
		} // end omp

		//MPI_Waitall();
		MPI_Barrier(MPI_COMM_WORLD);

		//WRITE_LOG_MPI("iter time without send = " << _clock + omp_get_wtime() << '\n', myid);

		MPI_Gatherv(loc_illum.data(), loc_illum.size(), MPI_DOUBLE, grid.Illum, send_count_illum.data(), disp_illum.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);				
		MPI_Gather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//WRITE_LOG_MPI("time after MPI_Gatherv = " << _clock + omp_get_wtime() << '\n', myid);

		if (myid == 0)		
		{			
			Type int_time = -omp_get_wtime();			
			GetIntScattering(count_cells, grid_direction, grid);
				
		//	WRITE_LOG_MPI("time int= " << int_time + omp_get_wtime() << '\n', myid);
		}

		MPI_Barrier(MPI_COMM_WORLD); // ждем расчёт интеграла рассеяния		
		MPI_Scatterv(grid.scattering, send_count_scattering.data(), disp_scattering.data(), MPI_DOUBLE,
					 int_scattering_local.data(), int_scattering_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//WRITE_LOG_MPI("time after MPI_Scatterv = " << _clock + omp_get_wtime() << '\n', myid);

		_clock += omp_get_wtime();
		if (myid == 0)
		{
			for (auto n : norms) if (n > norm) norm = n; 		
						
			WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');
		}
		MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		count++;
	} while (norm > solve_mode.accuracy && count < solve_mode.max_number_of_iter);
	
#ifndef USE_CUDA
	if (myid == 0) // это уйдет, когда все интегралы перейдут в cuda
	{		
		ReCalcIllumGlobal(grid_direction.size, grid);		
		//WriteFileSolution(BASE_ADRESS + "Illum_result.txt", Illum);
	}	
#endif
		
	return 0;
}

//-----------------------------------------------------//


int MPI_CalculateIllumAsync(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell, grid_t& grid)
{
	const int count_directions = grid_direction.size;
	const int count_cells = grid.size;

	int np = 1, myid = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	const int local_size = send_count[myid];
	const int local_disp = disp[myid];

	int count = 0;
	Type norm = -1;
	std::vector<Type> norms;

	if (myid == 0)
	{
		int i = 0;
		for (auto& el : grid.cells)
		{
			phys_local[i++] = el.phys_val;
			// сюда же absorp_coef, если понадобиться
		}

		norms.resize(np, -10);
	}
	MPI_Bcast(phys_local.data(), count_cells, MPI_flux_t, 0, MPI_COMM_WORLD);

	std::vector<MPI_Request> requests; //все запросы сообщений отпраки и принятия	
	std::vector<MPI_Status> status; //статусы всех обменов
	std::vector<int> flags_send_to_gpu;  //флаги указывающие на отправку пакета на gpu
	if (myid == 0)
	{
		const int N = count_directions - local_size;
		requests.resize(N);
		status.resize(N);
		flags_send_to_gpu.resize(N, 0);
	}
	else
	{
		requests.resize(local_size);
		status.resize(local_size);
	}

	do {

		Type _clock = -omp_get_wtime();
		norm = -1;
	
		// инициализируем весь прием сразу (от первого процесса ничего не принимаем)
		if (myid == 0)
		{
			int size_msg = grid.size * base;
			int src = 1;
			for (size_t src = 1; src < np; src++)
			{
				for (size_t j = 0; j < send_count[src]; j++)
				{
					int tag = disp[src] + j;					

					MPI_Irecv(grid.Illum + size_msg * tag, size_msg, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &requests[tag - local_size]);

					//WRITE_LOG("id= " << myid << "recv_teg= " << tag << " src= " << src << "\n");
				}
			}

			flags_send_to_gpu.assign(flags_send_to_gpu.size(), 0);			
		}
		
#pragma omp parallel default(none) shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
		norm, inter_coef_all,local_size, int_scattering_local, phys_local, loc_illum, disp, myid, np, \
requests, status,flags_send_to_gpu, BASE_ADRESS)
		{
			Type loc_norm = -1;
			const int num = omp_get_thread_num();
			std::vector<Vector3>* inter_coef = &inter_coef_all[num];

			/*---------------------------------- далее FOR по направлениям----------------------------------*/
#pragma omp for
			for (register int num_direction = 0; num_direction < local_size; ++num_direction)
			{
				int posX0 = 0;
				Vector3 I;

				/*---------------------------------- далее FOR по ячейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h)
				{
					const int num_cell = sorted_id_cell[num_direction][h];
					const int face_block_id = num_cell * base;

					for (ShortId num_out_face = 0; num_out_face < base; ++num_out_face)
					{
						const int out_id_face = pairs[face_block_id + num_out_face];
						if (check_bit(face_states[num_direction][num_cell], (int)num_out_face))
						{
							//формально обрабатываем мы только выходящие грани, т.к. входящие к этому моменту определены,
							// но границу надо инициализировать. Не для расчета, но для конечного интегрирования
							if (out_id_face < 0)
							{
								BoundaryConditions(out_id_face, (*inter_coef)[face_block_id + num_out_face]);
							}
							continue;
						}

						//GetNodes
						for (int num_node = 0; num_node < 3; ++num_node)
						{
							const Vector3 x = vec_x[num_cell].x[num_out_face][num_node];

							const cell_local x0 = vec_x0[num_direction][posX0++];

							const ShortId num_in_face = x0.in_face_id;
							const Type s = x0.s;
							const Vector2 X0 = x0.x0;

							const Type I_x0 = CalculateIllumeOnInnerFace(pairs[face_block_id + num_in_face], (*inter_coef)[face_block_id + num_in_face]);

							I[num_node] = GetCurIllum(x, s, I_x0, int_scattering_local[num_direction * count_cells + num_cell], phys_local[num_cell]);

						}//num_node

						(*inter_coef)[face_block_id + num_out_face] = I;
						if (out_id_face >= 0)
						{
							(*inter_coef)[out_id_face] = I;
						}

						//сюда сразу в loc_illum. сходимость ез нормы задать числом итераций

					} //num_out_face	
				}

				/*---------------------------------- конец FOR по ячейкам----------------------------------*/
				loc_norm = ReCalcIllum(num_direction, *inter_coef, loc_illum);

				const int n = count_cells * base;
				if (myid != 0)
				{
#pragma omp critical
					{
						const int teg = disp[myid] + num_direction;  //teg соответствует номеру направления
						MPI_Isend(loc_illum.data() + num_direction * n, n, MPI_DOUBLE, 0, teg, MPI_COMM_WORLD, &requests[num_direction]);

						//WRITE_LOG_MPI("id= " << myid <<"send_teg= " << teg<< "\n", myid);
					}
				}
				else
				{
#pragma omp critical//master
					{
						for (int i = 0; i < requests.size(); i++)
						{
							if (!flags_send_to_gpu[i])
							{								
								MPI_Request_get_status(requests[i], &flags_send_to_gpu[i], &status[i]); //проверяем все запросы принятия сообщения
								if (flags_send_to_gpu[i]) // если обмен завершён, но отправки не было
								{
#ifdef USE_CUDA
									CudaSendIllumAsync(n, n * (i + local_size), grid.Illum);  //переслать данные на gpu									
#endif
								}
							}
						}
					}
#pragma omp critical
					{
#ifdef USE_CUDA
						CudaSendIllumAsync(n, (num_direction * n), loc_illum.data());
#endif
					}				
				}
			}

			/*---------------------------------- конец FOR по направлениям----------------------------------*/

			if (loc_norm > norm)
			{
#pragma omp critical
				{
					if (loc_norm > norm)
					{
						norm = loc_norm;
					}
				}
			}
		} // end omp	
		
		MPI_Waitall(requests.size(), requests.data(), status.data());		

		if ((myid == 0) && (solve_mode.max_number_of_iter > 1))
		{
#ifdef USE_CUDA
			CalculateIntScatteringAsync(grid_direction, grid);

			for (int i = 0; i < loc_illum.size(); i++)
			{
				grid.Illum[i] = loc_illum[i];
			}
#else
			for (int i = 0; i < loc_illum.size(); i++)
			{
				grid.Illum[i] = loc_illum[i];
			}
			CalculateIntCPU(count_cells, grid_direction, grid);
#endif


		}

		MPI_Gather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (myid == 0)
		{
			for (auto n : norms) if (n > norm) norm = n;
		}
		MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef USE_CUDA
		if (myid == 0)
		{
			CudaWait();
		}
#endif

		// в теории эту отправку тоже можно разбить на выбранные направления
		MPI_Scatterv(grid.scattering, send_count_scattering.data(), disp_scattering.data(), MPI_DOUBLE,
			int_scattering_local.data(), int_scattering_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		_clock += omp_get_wtime();
		if (myid == 0)
		{
			WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');
		}
		count++;
	} while (norm > solve_mode.accuracy && count < solve_mode.max_number_of_iter);


#ifndef USE_CUDA
	if (myid == 0) // это уйдет, когда все интегралы перейдут в cuda
	{
		ReCalcIllumGlobal(grid_direction.size, grid);
		//WriteFileSolution(BASE_ADRESS + "Illum_result.txt", Illum);
	}
#endif

	return 0;
}

#ifndef USE_CUDA

#define GET_FACE_TO_CELL(val, data, init){ \
val = init; \
for (int pos = 0; pos < base; pos++) \
{ \
	val += data[pos]; \
} \
val /= base; \
} 

static Type IntegarteDirection(const vector<Type>& Illum, const grid_directions_t& grid_direction)
{		
	Type res = 0;
	Type illum_cell;
	int i = 0;
	for (auto& dir :grid_direction.directions)
	{		
		GET_FACE_TO_CELL(illum_cell, (&Illum[i * base]), 0);		
		res += illum_cell * dir.area;
		i++; // = base;
	}	
	return res / grid_direction.full_area;
}
static int MakeEnergy(const grid_directions_t& grid_direction, std::vector<elem_t>& cells) {
		
	//for (auto &el : cells)
#pragma omp parallel  default(none) shared(grid_direction, cells) 
	{
#pragma omp for
		for (int i = 0; i < cells.size(); i++)
		{
			elem_t& el = cells[i];
			el.illum_val.energy = IntegarteDirection(el.illum_val.illum, grid_direction);
		}
	}
	return 0;
}

static int IntegarteDirection3(const vector<Type>& Illum, const grid_directions_t& grid_direction, Vector3* stream_face)
{
	int i = 0;
	for (int f = 0; f < base; f++)
	{
		stream_face[f] = Vector3::Zero();
	}

	for (auto& dir : grid_direction.directions)
	{
		for (int f = 0; f < base; f++)
		{
			stream_face[f] += Illum[i++] * dir.area * dir.dir;
		}
	}

	for (int f = 0; f < base; f++)
	{
		stream_face[f] /= grid_direction.full_area;
	}

	return 0;
}
static int MakeStream(const grid_directions_t& grid_direction, grid_t& grid) 
{
#pragma omp parallel  default(none) shared(grid_direction, grid) 
	{
#pragma omp for
		for (int i = 0; i < grid.cells.size(); i++)
		{
			elem_t& el = grid.cells[i];

			Vector3 Stream[base];
			IntegarteDirection3(el.illum_val.illum, grid_direction, Stream);

			GET_FACE_TO_CELL(el.illum_val.stream, Stream, Vector3::Zero());

			el.illum_val.div_stream = 0;
			for (int j = 0; j < base; j++)
			{
				geo_face_t* geo_f = &grid.faces[el.geo.id_faces[j]].geo;
				if (el.geo.sign_n[j])
				{
					el.illum_val.div_stream += Stream[j].dot(geo_f->n) * geo_f->S;
				}
				else
				{
					el.illum_val.div_stream -= Stream[j].dot(geo_f->n) * geo_f->S;
				}
			}
			el.illum_val.div_stream /= el.geo.V;
		}
	}
	return 0;
}


static int IntegarteDirection9(const vector<Type>& Illum, const grid_directions_t& grid_direction, Matrix3* impuls_face) {
	
	int i = 0;
	for (int f = 0; f < base; f++)
	{
		impuls_face[f] = Matrix3::Zero();
	}

	for (auto& dir : grid_direction.directions)
	{
		for (int f = 0; f < base; f++)
		{			
			for (size_t h = 0; h < 3; h++)
				for (size_t k = 0; k < 3; k++)
					{
						impuls_face[f](h, k) += dir.dir[h] * dir.dir[k] * (Illum[i] * dir.area);
					}
			i++;			
		}
	}

	for (int f = 0; f < base; f++)
	{
		impuls_face[f] /= grid_direction.full_area;
	}
	
	return 0;	
}
static int MakeDivImpuls(const grid_directions_t& grid_direction, grid_t& grid)
{
#pragma omp parallel  default(none) shared(grid_direction, grid) 
	{
#pragma omp for
		for (int i = 0; i < grid.cells.size(); i++)
		{
			elem_t& el = grid.cells[i];

			Matrix3 Impuls[base];
			IntegarteDirection9(el.illum_val.illum, grid_direction, Impuls);

			GET_FACE_TO_CELL(el.illum_val.impuls, Impuls, Matrix3::Zero());

			el.illum_val.div_impuls = Vector3::Zero();
			for (int j = 0; j < base; j++)
			{
				geo_face_t* geo_f = &grid.faces[el.geo.id_faces[j]].geo;
				if (el.geo.sign_n[j])
				{
					el.illum_val.div_impuls += Impuls[j] * (geo_f->n) * geo_f->S;
				}
				else
				{
					el.illum_val.div_impuls += Impuls[j] * (-geo_f->n) * geo_f->S;
				}
			}
			el.illum_val.div_impuls /= el.geo.V;
		}

	}
	return 0;
}

int TestDivStream(const std::vector<Vector3>& centers_face, grid_t& grid)
{	
#pragma omp parallel  default(none) shared(centers_face, grid) 
	{
#pragma omp for
		for (int i = 0; i < grid.cells.size(); i++)
		{
			elem_t& el = grid.cells[i];

			Vector3 Stream[base];
			for (int j = 0; j < base; j++)
			{
				Vector3 x = centers_face[i * base + j];
				Stream[j] = Vector3(3 * x[0], 0, 0);
			}

			GET_FACE_TO_CELL(el.illum_val.stream, Stream, Vector3::Zero());

			el.illum_val.div_stream = 0;
			for (int j = 0; j < base; j++)
			{
				geo_face_t* geo_f = &grid.faces[el.geo.id_faces[j]].geo;
				if (el.geo.sign_n[j])
				{
					el.illum_val.div_stream += Stream[j].dot(geo_f->n) * geo_f->S;
				}
				else
				{
					el.illum_val.div_stream -= Stream[j].dot(geo_f->n) * geo_f->S;
				}
			}

			el.illum_val.div_stream /= el.geo.V;
		}
	}
	return 0;
}
#endif //!CUDA

int CalculateIllumParam(const grid_directions_t& grid_direction, grid_t& grid)
{
#ifdef USE_CUDA
	CalculateParamOnCuda(grid_direction, grid);
#else
	{
		MakeEnergy(grid_direction, grid.cells);
		MakeStream(grid_direction, grid);
		MakeDivImpuls(grid_direction, grid);
	}
#endif
	return 0;
}


#endif //ILLUM