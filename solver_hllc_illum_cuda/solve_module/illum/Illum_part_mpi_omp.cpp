#include "../solve_config.h"
#if 0// defined ILLUM && defined SOLVE && defined USE_MPI

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

std::vector<flux_t> phys_local; 
//static std::vector<Vector3> inter_coef; 
static std::vector<std::vector<Vector3>> inter_coef_all;
static std::vector<Type> loc_illum; 	
static std::vector<Type>int_scattering_local;

static std::vector<MPI_Request> requests; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status; //статусы всех обменов
static std::vector<int> flags_send_to_gpu;  //флаги указывающие на отправку пакета на gpu

static std::vector<MPI_Request> requests_send; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_send; //статусы всех обменов

static std::vector<MPI_Request> requests_hllc; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_hllc; //статусы всех обменов
std::vector<int> send_hllc;
std::vector<int> disp_hllc;

static int size_first_section;
static std::vector<MPI_Request> requests_first_section; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_first_section; //статусы всех обменов
static std::vector<int> flags_send_to_gpu_first_section;  //флаги указывающие на отправку пакета на gpu
static std::vector<MPI_Request> requests_send_first_section; //все запросы сообщений отпраки и принятия	
static std::vector<int> sizes_first_section;


static std::vector<MPI_Request> requests_second_section; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_second_section; //статусы всех обменов
static std::vector<int> flags_send_to_gpu_second_section;  //флаги указывающие на отправку пакета на gpu
static std::vector<MPI_Request> requests_send_second_section; //все запросы сообщений отпраки и принятия	

int GetSend(const int myid) { return send_count[myid]; }
int GetDisp(const int myid) { return disp[myid]; }

int InitPhysOmpMpi(const int count_cells)
{
	const int nthr = 6; //разделим иходный массив на 4 равных куска
	
	GetSend(nthr, count_cells, send_hllc);
	GetDisp(nthr, count_cells, disp_hllc);
		
	requests_hllc.resize(nthr);
	status_hllc.resize(nthr);

	return 0;

}
void SendPhysValue(flux_t* phys, const int size, const int msg_id)
{	
	MPI_Ibcast(phys, size, MPI_flux_t, 0, MPI_COMM_WORLD, &requests_hllc[msg_id]);	
}

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

	//loc_illum.resize(count_cells * base * local_size, 0); //кроме id 0	
	//int_scattering_local.resize(count_cells * local_size, 0);
	
	return 0;
}
void MPI_INIT(const int myid, const int np, const int count_directions, const grid_t& grid)
{
	size_first_section = 3 * send_count[0] / 4; // первый узел будем потенциально разгружать(на всех узла должно хватать направлений)
	
	sizes_first_section.resize(np, grid.size * base * size_first_section);
	//{
	//	//==================first rcv ==============================
	//	requests_first_section.resize(np - 1, MPI_REQUEST_NULL);
	//	status_first_section.resize(np - 1);
	//	flags_send_to_gpu_first_section.resize(np - 1, 0);
	//	const int size_msg = grid.size * base * size_first_section;

	//	int cc = 0;
	//	for (int src = 0; src < np; src++)
	//	{
	//		if (src == myid) continue;
	//		int tag = src;
	//		MPI_Recv_init(grid.Illum + size_msg * tag, size_msg, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &requests_first_section[cc++]);
	//	}

	//	//==================first send ==============================

	//	cc = 0;
	//	requests_send_first_section.resize(np - 1, MPI_REQUEST_NULL);
	//	for (int id = 0; id < np; id++)
	//	{
	//		if (id == myid) continue;
	//		MPI_Send_init(grid.Illum + size_msg * myid, size_msg, MPI_DOUBLE, id, myid, MPI_COMM_WORLD, &requests_send_first_section[cc++]);
	//	}
	//}
	//======================MPI_INIT=========================

	const int local_size = send_count[myid];
	const int size_second_section = local_size - size_first_section;
	{
		const int N = count_directions - (size_first_section * np) - size_second_section;
		requests_second_section.resize(N, MPI_REQUEST_NULL);
		status_second_section.resize(N);
		flags_send_to_gpu_second_section.resize(N, 0);

		int size_msg = grid.size * base;
		int cc = 0;
		for (int src = 0; src < np; src++)
		{
			if (src == myid) continue;
			for (int j = size_first_section; j < send_count[src]; j++)
			{
				int tag = disp[src] + j;
				MPI_Recv_init(grid.Illum + size_msg * tag, size_msg, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &requests_second_section[cc++/*tag - local_size*/]);
			}
		}

		requests_send_second_section.resize(size_second_section * (np - 1), MPI_REQUEST_NULL);

		for (size_t num_direction = size_first_section; num_direction < local_size; num_direction++)
		{
			const int tag = disp[myid] + num_direction;  //teg соответствует номеру направления
			cc = 0;
			for (int id = 0; id < np; id++)
			{
				if (id == myid) continue;

				MPI_Send_init(grid.Illum + (disp[myid] + num_direction) * size_msg, size_msg, MPI_DOUBLE, id, tag, MPI_COMM_WORLD, &requests_send_second_section[(np - 1) * (num_direction- size_first_section)+cc++]);
			}
		}
	}

	return;
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
			WRITE_LOG_ERR("I_x0= " << I_x0 << '\n');
			D_LD;
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
	int np, myid;
	MPI_GET_INF(np, myid);
	const int disp_loc = disp[myid];
	const int size_loc = send_count[myid];


#pragma omp parallel default(none) shared(num_cells, grid_direction, grid, disp_loc, size_loc)
	{
		Vector3 direction;
		const int num_directions = grid_direction.size;
#pragma omp for
		for (int num_direction = 0; num_direction < size_loc; ++num_direction) {
			direction = grid_direction.directions[disp_loc + num_direction].dir;
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
static Type ReCalcIllum(const int num_dir, const std::vector<Vector3>& inter_coef, Type*/*std::vector<Type>&*/ Illum)
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
	
#pragma omp parallel default(none) shared(dir_size, grid)
	{
		const int size = grid.size;
#pragma omp for
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
				loc_norm = ReCalcIllum(num_direction, *inter_coef, loc_illum.data());
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

	int np, myid;
	MPI_GET_INF(np, myid);

	const int local_size = send_count[myid];
	const int local_disp = disp[myid];

	int count = 0;
	Type norm = -1;
	std::vector<Type> norms(np, -10);

	struct {
		Type phys_time;
		Type rcv_init_time;

		Type send_wait_time;
		Type dir_time;
		Type wait_all_time;
		Type cuda_time;
		Type cuda_wait_time;
		Type time_sender;
		Type wait_scat_time;
		Type norm_cast_time;
	}timer;


	static MPI_Request rq_first_section = MPI_REQUEST_NULL;
	int first_flag = 0;
	do {

		Type _clock = -omp_get_wtime();
		norm = -1;
					
		//timer.rcv_init_time = -omp_get_wtime();
		//// инициализируем весь прием сразу (от первого процесса ничего не принимаем)		
		//if (np > 1)
		//{
		//	flags_send_to_gpu_first_section.assign(flags_send_to_gpu_first_section.size(), 0);
		//	if (requests_first_section[0] != MPI_REQUEST_NULL)
		//	{				
		//		MPI_Waitall(requests_first_section.size(), requests_first_section.data(), MPI_STATUSES_IGNORE); // status.data());
		//	}
		//	MPI_Startall(requests_first_section.size(), requests_first_section.data());
		//}
		//timer.rcv_init_time += omp_get_wtime();


		timer.phys_time = -omp_get_wtime();
		if (count == 0)
		{
			MPI_Waitall(requests_hllc.size(), requests_hllc.data(), MPI_STATUSES_IGNORE); //ждём сообщнения с газовым расчётом		
		}
		timer.phys_time += omp_get_wtime();

	/*	timer.send_wait_time = -omp_get_wtime();
		if (requests_send_first_section[0] != MPI_REQUEST_NULL)
		{
			MPI_Waitall(requests_send_first_section.size(), requests_send_first_section.data(), MPI_STATUSES_IGNORE);
		}
		timer.send_wait_time += omp_get_wtime();*/

#ifdef USE_CUDA
		// это перед вызовами теперь, т.к. мы не ждём для рассылки
		timer.cuda_wait_time = -omp_get_wtime();
		CudaWait();
		timer.cuda_wait_time += omp_get_wtime();
#endif
		
		timer.dir_time = -omp_get_wtime();
		timer.time_sender = 0;

		WRITE_LOG_MPI("start dir\n", myid);
		/*---------------------------------- далее FOR по направлениям----------------------------------*/

#pragma omp parallel default(none) shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
		norm, inter_coef_all,local_size, int_scattering_local, phys_local, loc_illum, disp, myid, np, \
requests, status,flags_send_to_gpu, BASE_ADRESS,send_count,requests_send,disp_illum,timer, \
		size_first_section, requests_first_section, status_first_section,flags_send_to_gpu_first_section,\
 requests_second_section, status_second_section,flags_send_to_gpu_second_section,\
		requests_send_first_section, requests_send_second_section, sizes_first_section,rq_first_section,first_flag)
		{
			const int n = count_cells * base;
			const int num_of_th = omp_get_num_threads();

			Type loc_norm = -1;
			const int num = omp_get_thread_num();
			std::vector<Vector3>* inter_coef = &inter_coef_all[num];

#pragma omp for
			for (register int num_direction = 0; num_direction < size_first_section; num_direction++)
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

							I[num_node] = GetCurIllum(x, s, I_x0, grid.scattering[num_direction * count_cells + num_cell], phys_local[num_cell]);

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
				loc_norm = ReCalcIllum(num_direction, *inter_coef, grid.Illum + disp_illum[myid]);

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
			}
			/*--------------------------------конец первой секции---------------------------------*/
		//////}
				
#pragma omp flush

#pragma omp single
			{
				first_flag = false;
				MPI_Iallgatherv(MPI_IN_PLACE/*grid.Illum + disp_illum[myid]*/, size_first_section*base*count_cells, MPI_DOUBLE, grid.Illum, sizes_first_section.data(), disp_illum.data(),
					MPI_DOUBLE, MPI_COMM_WORLD, &rq_first_section);

				//allgatherv?
		//#pragma omp single //nowait
				{
					WRITE_LOG_MPI("first end dir\n", myid);

					if (np > 1)
					{
						//		MPI_Startall(requests_send_first_section.size(), requests_send_first_section.data());
					}
					WRITE_LOG_MPI("1\n", myid);
				}
#ifdef USE_CUDA
				//#pragma omp single //nowait
				{
					const int n = count_cells * base;
					CudaSendIllumAsync(size_first_section * n, disp_illum[myid], grid.Illum);
				}
#endif
				WRITE_LOG_MPI("2\n", myid);
				if (np > 1)
				{
					//#pragma omp single //nowait
					{
						flags_send_to_gpu_second_section.assign(flags_send_to_gpu_second_section.size(), 0);
					}
					//#pragma omp single //nowait
					{
						if (requests_second_section[0] != MPI_REQUEST_NULL)
						{
							MPI_Waitall(requests_second_section.size(), requests_second_section.data(), MPI_STATUSES_IGNORE); // status.data());
						}
					}
					WRITE_LOG_MPI("3\n", myid);
					//#pragma omp single //nowait
					{
						MPI_Startall(requests_second_section.size(), requests_second_section.data());
					}
					//#pragma omp single //nowait
					WRITE_LOG_MPI("4\n", myid);
					{
						if (requests_send_second_section[0] != MPI_REQUEST_NULL)
						{
							MPI_Waitall(requests_send_second_section.size(), requests_send_second_section.data(), MPI_STATUSES_IGNORE);
						}
					}
					WRITE_LOG_MPI("5\n", myid);

				}

				WRITE_LOG_MPI("start second part \n", myid);
			}

//#pragma omp barrier
////#pragma omp parallel default(none) shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
////		norm, inter_coef_all,local_size, int_scattering_local, phys_local, loc_illum, disp, myid, np, \
////requests, status,flags_send_to_gpu, BASE_ADRESS,send_count,requests_send,disp_illum,timer, \
////		size_first_section, requests_first_section, status_first_section,flags_send_to_gpu_first_section,\
//// requests_second_section, status_second_section,flags_send_to_gpu_second_section,\
////		requests_send_first_section, requests_send_second_section)
//			{
//				const int n = count_cells * base;
//				const int num_of_th = omp_get_num_threads();
//
//				Type loc_norm =norm;
//				const int num = omp_get_thread_num();
//				std::vector<Vector3>* inter_coef = &inter_coef_all[num];
#pragma omp for
			for (register int num_direction = size_first_section; num_direction < local_size; num_direction++)
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

							I[num_node] = GetCurIllum(x, s, I_x0, grid.scattering[num_direction * count_cells + num_cell], phys_local[num_cell]);

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
				loc_norm = ReCalcIllum(num_direction, *inter_coef, grid.Illum + disp_illum[myid]);

#pragma omp critical
				{
					//todo: dynamic shadule?
					if (num == 0 && np > 1) //вместо критической секции пусть всегда первая нить управляет отправкой
					{
						// пересылаем первую пачку сообщений
#ifdef USE_CUDA			
						if (!first_flag)
						{
							MPI_Status st;
							MPI_Request_get_status(rq_first_section, &first_flag, &st);

							if (first_flag)
							{
								for (int src = 0; src < np; src++)
								{
									if (src == myid) continue;
									CudaSendIllumAsync(n * size_first_section, n * disp[src], grid.Illum);  //переслать данные на gpu									
								}
							}
						}
						//for (int i = 0; i < requests_first_section.size(); i++)
						//{
						//	if (!flags_send_to_gpu_first_section[i])
						//	{
						//		MPI_Request_get_status(requests_first_section[i], &flags_send_to_gpu_first_section[i], &status_first_section[i]); //проверяем все запросы принятия сообщения
						//		if (flags_send_to_gpu_first_section[i]) // если обмен завершён, но отправки не было
						//		{
						//			const int src = status_first_section[i].MPI_TAG;
						//			CudaSendIllumAsync(n * size_first_section, n * disp[src], grid.Illum);  //переслать данные на gpu									
						//		}
						//	}
						//}
#endif
					}

					if (np > 1)
					{
						//#pragma omp critical
						{
							MPI_Startall(np - 1, requests_send_second_section.data() + ((num_direction - size_first_section) * (np - 1)));
						}
					}
#ifdef USE_CUDA	
					if ((num == (num_of_th - 1)) && np > 1) // отличный от нулевого поток
					{
						for (int i = 0; i < requests_second_section.size(); i++)
						{
							if (!flags_send_to_gpu_second_section[i])
							{
								MPI_Request_get_status(requests_second_section[i], &flags_send_to_gpu_second_section[i], &status_second_section[i]); //проверяем все запросы принятия сообщения
								if (flags_send_to_gpu_second_section[i]) // если обмен завершён, но отправки не было
								{
									CudaSendIllumAsync(n, n * (status_second_section[i].MPI_TAG), grid.Illum);  //переслать данные на gpu									
								}
							}
						}
					}

					//#pragma omp critical
					{
						CudaSendIllumAsync(n, ((disp[myid] + num_direction) * n), grid.Illum);
					}
#endif		
				}
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
			}
			/*---------------------------------- конец FOR по направлениям----------------------------------*/
			

		} // end task	
		timer.dir_time += omp_get_wtime();

		WRITE_LOG_MPI("end second part \n", myid);
	
		timer.norm_cast_time = -omp_get_wtime();
		//todo: проверить можно ли здесь с блокировкой
		MPI_Request rq_norm;		
		MPI_Iallgather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD, &rq_norm);	
		MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
		for (auto n : norms) if (n > norm) norm = n;
		timer.norm_cast_time += omp_get_wtime();

		WRITE_LOG_MPI("wait \n", myid);
		timer.wait_all_time = -omp_get_wtime();
#ifdef USE_CUDA
		bool ready = true;
		do
		{
			const int n = count_cells * base;
			ready = true;

			if (!first_flag)
			{
				ready = false;				
				MPI_Request_get_status(rq_first_section, &first_flag, MPI_STATUS_IGNORE);
				//MPI_Test();

				if (first_flag)
				{
					for (int src = 0; src < np; src++)
					{
						if (src == myid) continue;
						CudaSendIllumAsync(n * size_first_section, n * disp[src], grid.Illum);  //переслать данные на gpu									
					}
				}
			}

			//for (int i = 0; i < requests_first_section.size(); i++)
			//{
			//	if (!flags_send_to_gpu_first_section[i])
			//	{
			//		WRITE_LOG_MPI("first "<<i<<" \n", myid);
			//		ready = false;
			//		MPI_Request_get_status(requests_first_section[i], &flags_send_to_gpu_first_section[i], &status_first_section[i]); //проверяем все запросы принятия сообщения
			//		if (flags_send_to_gpu_first_section[i]) // если обмен завершён, но отправки не было
			//		{
			//			const int src = status_first_section[i].MPI_TAG;
			//			CudaSendIllumAsync(n* size_first_section, n* disp[src], grid.Illum);
			//		}
			//	}
			//}
			
			for (int i = 0; i < requests_second_section.size(); i++)
			{
				if (!flags_send_to_gpu_second_section[i])
				{
					WRITE_LOG_MPI("second " << i << "\n", myid);
					ready = false;
					MPI_Request_get_status(requests_second_section[i], &flags_send_to_gpu_second_section[i], &status_second_section[i]); //проверяем все запросы принятия сообщения
					if (flags_send_to_gpu_second_section[i]) // если обмен завершён, но отправки не было
					{						
						CudaSendIllumAsync(n, n * (status_second_section[i].MPI_TAG), grid.Illum);  //переслать данные на gpu									
					}
				}
			}

		} while (!ready);
#else		
	
		MPI_Wait(&rq_first_section, MPI_STATUS_IGNORE);
		//MPI_Waitall(requests_first_section.size(), requests_first_section.data(), MPI_STATUSES_IGNORE);		
		MPI_Waitall(requests_second_section.size(), requests_second_section.data(), MPI_STATUSES_IGNORE);
#endif
		timer.wait_all_time += omp_get_wtime();
		WRITE_LOG_MPI("end wait \n", myid);

		timer.cuda_time = -omp_get_wtime();
		if (solve_mode.max_number_of_iter > 1)
		{
#ifdef USE_CUDA
			CalculateIntScatteringAsyncMPI(grid_direction, grid, local_size);
#else

			CalculateIntCPU(count_cells, grid_direction, grid);
#endif
		}
		timer.cuda_time += omp_get_wtime();

		_clock += omp_get_wtime();

		if (myid == 0)
		{
			WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');
		}

		WRITE_LOG_MPI(count << ": phys= " << timer.phys_time
			<< ", rcv " << timer.rcv_init_time
			<< ", send " << timer.send_wait_time
			<< ", dir " << timer.dir_time
			<< ", wait_all " << timer.wait_all_time
			<< ", cuda " << timer.cuda_time
			<< ", cuda_wait " << timer.cuda_wait_time
			<< ", norm_cast_time " << timer.norm_cast_time
			<< ", sender " << timer.time_sender
			<< ", all_step " << _clock << "\n\n", myid);

		count++;

		//MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); когда не будет лога

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
	CalculateAllParam(grid_direction, grid);
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

#if 0
int MPI_CalculateIllumAsyncTask(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell, grid_t& grid)
{
	const int count_directions = grid_direction.size;
	const int count_cells = grid.size;

	int np, myid;
	MPI_GET_INF(np, myid);

	const int local_size = send_count[myid];
	const int local_disp = disp[myid];

	int count = 0;
	Type norm = -1;
	std::vector<Type> norms;

	struct {
		Type phys_time;
		Type rcv_init_time;
		Type dir_time;
		Type wait_all_time;
		Type cuda_time;
		Type cuda_wait_time;
		Type scat_time;
	}timer;

	timer.phys_time = -omp_get_wtime();
	if (myid == 0)
	{
		//int i = 0;
		//for (auto& el : grid.cells)
		//{
		//	phys_local[i++] = el.phys_val;
		//	// сюда же absorp_coef, если понадобиться
		//}

		norms.resize(np, -10);
	}

	//MPI_Bcast(phys_local.data(), count_cells, MPI_flux_t, 0, MPI_COMM_WORLD);
	timer.phys_time += omp_get_wtime();

	//MPIWAITPHYS
#pragma omp parallel default(none) shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
		norm, inter_coef_all,local_size, int_scattering_local, phys_local, loc_illum, disp, myid, np, \
requests, status,flags_send_to_gpu, BASE_ADRESS, count,timer,grid_direction,status_hllc,requests_hllc, \
send_count_scattering,disp_scattering, norms, solve_mode)
	{
#pragma omp single
		{
			do {

				Type _clock = -omp_get_wtime();
				norm = -1;

				timer.rcv_init_time = -omp_get_wtime();

				// инициализируем весь прием сразу (от первого процесса ничего не принимаем)
				if (myid == 0)
				{
					MPI_Startall(requests.size(), requests.data());
					flags_send_to_gpu.assign(flags_send_to_gpu.size(), 0);
				}
				timer.rcv_init_time += omp_get_wtime();

				timer.dir_time = -omp_get_wtime();

				if (count == 0)
				{
					MPI_Waitall(requests_hllc.size(), requests_hllc.data(), status_hllc.data()); //ждём сообщнения с газовым расчётом
				}

				/*---------------------------------- далее FOR по направлениям----------------------------------*/

				const int NN = omp_get_num_threads();
				for (int id_thread = 0; id_thread < NN; id_thread++)
				{
#pragma omp task //firstprivate(id_thread)
					{
						Type loc_norm = -1;
						const int num = id_thread; // omp_get_thread_num();
						std::vector<Vector3>* inter_coef = &inter_coef_all[num];

						//#pragma omp for
						for (register int num_direction = id_thread; num_direction < local_size; num_direction += NN)
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
									MPI_Start(requests.data() + num_direction);
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

					} // end task	
				} //end threads

				timer.dir_time += omp_get_wtime();
				timer.wait_all_time = -omp_get_wtime();
				MPI_Waitall(requests.size(), requests.data(), status.data());
				timer.wait_all_time += omp_get_wtime();

				timer.cuda_time = -omp_get_wtime();
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

				timer.cuda_time += omp_get_wtime();

				MPI_Gather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (myid == 0)
				{
					for (auto n : norms) if (n > norm) norm = n;
				}
				MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				timer.cuda_wait_time = -omp_get_wtime();
#ifdef USE_CUDA
				if (myid == 0)
				{
					CudaWait();
				}
#endif
				timer.cuda_wait_time += omp_get_wtime();

				timer.scat_time = -omp_get_wtime();
				// в теории эту отправку тоже можно разбить на выбранные направления

				MPI_Scatterv(grid.scattering, send_count_scattering.data(), disp_scattering.data(), MPI_DOUBLE,
					int_scattering_local.data(), int_scattering_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				timer.scat_time += omp_get_wtime();
				_clock += omp_get_wtime();
				if (myid == 0)
				{
					WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');
				}

				WRITE_LOG_MPI(count << ": phys= " << timer.phys_time << ", rcv " << timer.rcv_init_time
					<< ", dir " << timer.dir_time
					<< ", wait_all " << timer.wait_all_time
					<< ", cuda " << timer.cuda_time
					<< ", cuda_wait " << timer.cuda_wait_time
					<< ", scater " << timer.scat_time
					<< ", all_step " << _clock << "\n\n", myid);

				count++;
			} while (norm > solve_mode.accuracy && count < solve_mode.max_number_of_iter);

		} //single
	}//parallel

#ifndef USE_CUDA
	if (myid == 0) // это уйдет, когда все интегралы перейдут в cuda
	{
		ReCalcIllumGlobal(grid_direction.size, grid);
		//WriteFileSolution(BASE_ADRESS + "Illum_result.txt", Illum);
	}
#endif

	return 0;
}
#endif