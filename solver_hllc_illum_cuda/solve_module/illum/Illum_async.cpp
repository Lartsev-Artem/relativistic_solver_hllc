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

static std::vector<int> disp_illum;

std::vector<flux_t> phys_local;
static std::vector<std::vector<Vector3>> inter_coef_all;

static std::vector<MPI_Request> requests_hllc; //все запросы сообщений отпраки и принятия	
std::vector<int> send_hllc;
std::vector<int> disp_hllc;

static int size_1_section;
static std::vector<MPI_Request> requests_rcv_1_section; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_rcv_1_section; //статусы всех обменов
static std::vector<int> flags_send_to_gpu_1_section;  //флаги указывающие на отправку пакета на gpu
static std::vector<MPI_Request> requests_send_1_section; //все запросы сообщений отпраки и принятия	

static std::vector<MPI_Request> requests_rcv_2_section; //все запросы сообщений отпраки и принятия	
static std::vector<MPI_Status> status_rcv_2_section; //статусы всех обменов
static std::vector<int> flags_send_to_gpu_2_section;  //флаги указывающие на отправку пакета на gpu
static std::vector<MPI_Request> requests_send_2_section; //все запросы сообщений отпраки и принятия	

int GetSend(const int myid) { return send_count[myid]; }
int GetDisp(const int myid) { return disp[myid]; }

int InitPhysOmpMpi(const int count_cells)
{
	int nthr = 1;   //не делим массив, если только один узел
	int np, id;
	MPI_GET_INF(np, id);
	if (np != 1)
	{
		nthr = 6;  //разделим иходный массив на 6 равных куска
		requests_hllc.resize(nthr);
	}

	GetSend(nthr, count_cells, send_hllc);
	GetDisp(nthr, count_cells, disp_hllc);
	
	return 0;
}
void SendPhysValue(flux_t* phys, const int size, const int msg_id)
{
	if (requests_hllc.size())
	{
		MPI_Ibcast(phys, size, MPI_flux_t, 0, MPI_COMM_WORLD, &requests_hllc[msg_id]);
	}
}

int InitSendDispIllumArray(const int myid, const int np, const int count_directions, const int count_cells)
{
	GetSend(np, count_directions, send_count);
	GetDisp(np, count_directions, disp);

	disp_illum.resize(np);
	for (int i = 0; i < np; i++)
	{	
		disp_illum[i] = disp[i] * base * count_cells;	
	}

	{
		int len[3 + 1] = { 1,3,1,  1 };
		MPI_Aint pos[4] = { offsetof(flux_t, d), offsetof(flux_t, v),offsetof(flux_t, p) ,sizeof(flux_t) };
		MPI_Datatype typ[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(4, len, pos, typ, &MPI_flux_t);
		MPI_Type_commit(&MPI_flux_t);
	}

	
	phys_local.resize(count_cells);

	inter_coef_all.resize(omp_get_max_threads());
	for (int i = 0; i < inter_coef_all.size(); i++)
	{
		inter_coef_all[i].resize(count_cells * base);
	}

	WRITE_LOG("Omp threads in illum= " << inter_coef_all.size() << "\n");	
	return 0;
}
void MPI_INIT(const int myid, const int np, const int count_directions, const grid_t& grid)
{
	const int size_first_section = 1 * send_count[0] / 3; // первый узел будем потенциально разгружать(на всех узла должно хватать направлений)
	size_1_section = size_first_section;
	{
		//==================first rcv ==============================
		requests_rcv_1_section.resize(np - 1, MPI_REQUEST_NULL);
		status_rcv_1_section.resize(np - 1);
		flags_send_to_gpu_1_section.resize(np - 1, 0);
		const int size_msg = grid.size * base * size_first_section;

		int cc = 0;
		for (int src = 0; src < np; src++)
		{
			if (src == myid) continue;
			int tag = src;
			MPI_Recv_init(grid.Illum + disp_illum[tag] /*size_msg * tag*/, size_msg, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &requests_rcv_1_section[cc++]);
		}

		//==================first send ==============================

		cc = 0;
		requests_send_1_section.resize(np - 1, MPI_REQUEST_NULL);
		for (int id = 0; id < np; id++)
		{
			if (id == myid) continue;
			MPI_Send_init(grid.Illum + disp_illum[myid]/*size_msg * myid*/, size_msg, MPI_DOUBLE, id, myid, MPI_COMM_WORLD, &requests_send_1_section[cc++]);
		}
	}
	
	//======================MPI_INIT=========================

	const int local_size = send_count[myid];
	const int size_second_section = local_size - size_first_section;
	{
		const int N = count_directions - (size_first_section * np) - size_second_section;
		requests_rcv_2_section.resize(N, MPI_REQUEST_NULL);
		status_rcv_2_section.resize(N);
		flags_send_to_gpu_2_section.resize(N, 0);

		int size_msg = grid.size * base;
		int cc = 0;
		for (int src = 0; src < np; src++)
		{
			if (src == myid) continue;
			for (int j = size_first_section; j < send_count[src]; j++)
			{
				int tag = disp[src] + j;
				MPI_Recv_init(grid.Illum + size_msg * tag, size_msg, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &requests_rcv_2_section[cc++/*tag - local_size*/]);
			}
		}

		requests_send_2_section.resize(size_second_section * (np - 1), MPI_REQUEST_NULL);

		for (size_t num_direction = size_first_section; num_direction < local_size; num_direction++)
		{
			const int tag = disp[myid] + num_direction;  //teg соответствует номеру направления
			cc = 0;
			for (int id = 0; id < np; id++)
			{
				if (id == myid) continue;

				MPI_Send_init(grid.Illum + (disp[myid] + num_direction) * size_msg, size_msg, MPI_DOUBLE, id, tag, MPI_COMM_WORLD, &requests_send_2_section[(np - 1) * (num_direction - size_first_section) + cc++]);
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

#ifdef DEBUG
		if (I_x0 < 0)
		{
			WRITE_LOG("Error illum value\n");
			WRITE_LOG_ERR("I_x0= " << I_x0 << '\n');
			D_LD;
			return 0;
		}
#endif
	}

	return I_x0;
}

static Type GetS(const int num_cell, const Vector3& direction, const grid_t& grid, const grid_directions_t& grid_direction) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
		Type dot = direction.dot(direction2);
	return (3. * (1 + dot * dot)) / 4;
	} };

	Type S = 0;
	const int N_dir = grid_direction.size;
	const int N_cell = grid.size * base; // / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		Type I = 0;
		for (int i = 0; i < base; i++)
		{
			I += grid.Illum[num_direction * N_cell + num_cell + i];
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


#pragma omp parallel default(none) firstprivate(num_cells, disp_loc, size_loc) shared(grid_direction, grid)
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
		Type Q = 0.01;
		Type alpha = 0.5;
		Type betta = 150; ////0.5;
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
#ifdef DEBUG
		static Type max = 0;
		if ((I / RADIATION) > max)
		{
			max = I / RADIATION;
			WRITE_LOG("I= " << max << ", S= " << S <<  '\n');
		}
#endif
		return I;
	}

	case 11: // HLLC + Illum для конуса
	{

		const Type S = int_scattering * RADIATION;

		const Type d = cell.d * DENSITY;
		const Type p = cell.p * PRESSURE;

		const Type T = std::min( 1e-6 * p / (d * R_gas), 5000.);

		const Type T1 = k_boltzmann * PI * T;
		const Type T2 = T1 * T1;
		const Type T4 = T2 * T2;

		const Type 	Ie = T4 * 2 / (15 * h_plank * h_plank * h_plank * c_light * c_light);
		const Type 	betta = d * (sigma_thomson / m_hydrogen);
		const Type	alpha = betta; 			//alpha = ???			


		const Type Q = alpha * Ie;

		const Type k = alpha + betta;

		const Type ss = s * DIST;
		const Type I0 = I_0 * RADIATION;

		Type I;
		if (k > 1e-10)
			I = (exp(-k * ss) * (I0 * k + (exp(k * ss) - 1) * (Q + S * betta))) / k;
		else
			I = (1 - ss * k) * (I0 + ss * (Q + S * betta));

		if (I < 0)
		{
			return 0;
		}

#ifdef DEBUG
		static Type max = 0;
		if ((I / RADIATION) > max) 
		{
			max = I / RADIATION;
			WRITE_LOG("I= " << max<<", T= "<<T<<", ie= " << Ie <<", I0="<< I0<< ", S="<< S << '\n');
		}
#endif

		return  I / RADIATION;
	}

	default:
		D_LD; // ("unknow class_vtk in get illum\n");
	}
}

static Type ReCalcIllum(const int num_dir, const std::vector<Vector3>& inter_coef, Type* Illum)
{
	Type norm = -1;
	//#pragma omp parallel default(none) shared(num_dir, cells, Illum, norm)
	{
		Type norm_loc = -1;
		const int size = inter_coef.size() / base;
		const int shift_dir = num_dir * size;
		Vector3 Il;

		//#pragma omp for
		for (int num_cell = 0; num_cell < size; num_cell++)
		{
			for (int i = 0; i < base; i++)
			{
				Il = inter_coef[num_cell * base + i]; 
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

static Type BoundaryConditionsNoOff(const int type_bound, const Vector3& x, Vector3& inter_coef)
{
	Type I0 = 0;

	const Type p = x.norm();
	const Type a = 1e3;
	const Type b = 0.01;
	const Type c = 10;

	switch (type_bound)
	{
	case eBound_OutSource: // дно конуса
		I0 = a * exp(-(p / b) * (p / b)) + c;
		//I0 = 1e4;
		break;
	case eBound_FreeBound:
		I0 = 0;
		break;

	case eBound_LockBound:
		I0 = 0;
		break;

	case eBound_InnerSource:  // внутренняя граница	
	{
		//I0 = a * exp(-(p / b) * (p / b)) + c;
		//I0 = 100;// 1e5;
		I0 = 0;
		break;
	}
	case 100:// дно конуса(апроксимация экспонеты: a Exp[- (x/b)^2] + c)  eBound_
				
#if 0 //через ряд до 2 и 3 члена
		if (p > 0.5)
		{
			// задать константой... без abc
		}
		else //ряд в нуле: (a + c) - (a x ^ 2) / b ^ 2 + (a x ^ 4) / (2 b ^ 4)
		{
			const Type x2 = p * p;
			const Type ax2 = a * x2;
			const Type b2 = 1./(b * b);
			I0 = (a + c) - ax2 * b2 + (ax2 * x2) * (0.5 * b2 * b2);
		}
#endif
		break;
	default:
		D_LD;
	}

	inter_coef = Vector3(I0, I0, I0);
	return I0;
}

int MPI_CalculateIllumAsync(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell, grid_t& grid)
{
	const int count_directions = grid_direction.size;
	const int count_cells = grid.size;
	const int n_illum = count_cells * base;

	int np, myid;
	MPI_GET_INF(np, myid);

	const int local_size = send_count[myid];
	const int local_disp = disp[myid];

	int count = 0;
	Type norm = -1;
	std::vector<Type> norms(np, -10);

	struct
	{
		Type phys_time;


		Type send1_wait_time;
		Type rcv1_wait_time;
		Type rcv2_wait_time;

		Type rcv2_init_time;
		Type send2_wait_time;

		Type dir_time;

		Type cuda_time;
		Type cuda_wait_time_1;
		Type cuda_wait_time_2;

		Type norm_gather;
	}timer;

	timer.phys_time = -omp_get_wtime();
	MPI_Waitall(requests_hllc.size(), requests_hllc.data(), MPI_STATUSES_IGNORE); //ждём сообщнения с газовым расчётом			
	timer.phys_time += omp_get_wtime();


	do {

		Type _clock = -omp_get_wtime();
		norm = -1;
		// эта проверка скрывает расчёт видеокарты
		timer.send1_wait_time = -omp_get_wtime();
		{
			if (np > 1 && requests_send_1_section[0] != MPI_REQUEST_NULL)
			{
				MPI_Waitall(requests_send_1_section.size(), requests_send_1_section.data(), MPI_STATUSES_IGNORE);
			}
		}
		timer.send1_wait_time += omp_get_wtime();
		timer.send2_wait_time = -omp_get_wtime();
		{
			if (np > 1 && requests_send_2_section[0] != MPI_REQUEST_NULL)
			{
				MPI_Waitall(requests_send_2_section.size(), requests_send_2_section.data(), MPI_STATUSES_IGNORE);
			}
		}
		timer.send2_wait_time += omp_get_wtime();

#ifdef USE_CUDA			
		timer.cuda_wait_time_1 = -omp_get_wtime();
		//CudaWait(); // это перед вызовами теперь, т.к. мы не ждём для рассылки
		CudaSyncStream(eCuda_scattering_1);
		if (myid == 0)
		{
			//самый высоких приоретет, т.к. надо расчитать, до конфликта с асинхронной отправкой
			CalculateAllParamStream(grid_direction, grid, eCuda_params); //запустим расчёт параметров здесь
		}
		// на выходе получим ответ за 1 шаг до сходимости, но зато без ожидания на выходе

		timer.cuda_wait_time_1 += omp_get_wtime();
#endif	

		timer.dir_time = -omp_get_wtime();				
		/*---------------------------------- далее FOR по направлениям----------------------------------*/

#pragma omp parallel default(none) firstprivate(local_size, np, n_illum, myid, size_1_section, local_disp) \
		shared(sorted_id_cell, pairs, face_states, vec_x0, vec_x, grid, \
		norm, inter_coef_all, phys_local, disp_illum,\
 BASE_ADRESS,timer, \
flags_send_to_gpu_1_section, requests_rcv_1_section, status_rcv_1_section, requests_send_1_section,\
flags_send_to_gpu_2_section, requests_rcv_2_section, status_rcv_2_section, requests_send_2_section)
		{

#pragma omp single
			{
				flags_send_to_gpu_1_section.assign(flags_send_to_gpu_1_section.size(), 0);
				if (np > 1)
				{
					MPI_Startall(requests_rcv_1_section.size(), requests_rcv_1_section.data());
				}
			}
			
			const int num_of_th = omp_get_num_threads();
			const int num = omp_get_thread_num();
						
			Type loc_norm = -1;			
			std::vector<Vector3>* inter_coef = &inter_coef_all[num];

#pragma omp for
			for (register int num_direction = 0; num_direction < size_1_section; num_direction++)
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
				Type buf_norm = ReCalcIllum(num_direction, *inter_coef, grid.Illum + disp_illum[myid]);
				if (buf_norm > loc_norm)
				{
					loc_norm = buf_norm;
				}
			}
			/*--------------------------------конец первой секции---------------------------------*/
		//////}			

#pragma omp flush

#pragma omp single
			{
				if (np > 1)
				{					
					MPI_Startall(requests_send_1_section.size(), requests_send_1_section.data());
				}
		
#ifdef USE_CUDA
				//#pragma omp single //nowait
				{					
					CudaSendIllumAsync(size_1_section * n_illum, disp_illum[myid],  grid.Illum);
				}
#endif			
			}			

#pragma omp single
			{
				if (np > 1)
				{
					//#pragma omp single //nowait
					{
						flags_send_to_gpu_2_section.assign(flags_send_to_gpu_2_section.size(), 0);
					}				
					{
						MPI_Startall(requests_rcv_2_section.size(), requests_rcv_2_section.data());
					}							
				}	

#ifdef USE_CUDA			
				timer.cuda_wait_time_2 = -omp_get_wtime();				
				CudaSyncStream(eCuda_scattering_2);
				timer.cuda_wait_time_2 += omp_get_wtime();
#endif
			}
			
			
#pragma omp for
			for (register int num_direction = size_1_section; num_direction < local_size; num_direction++)
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
				Type buf_norm = ReCalcIllum(num_direction, *inter_coef, grid.Illum + disp_illum[myid]);
				if (buf_norm > loc_norm)
				{
					loc_norm = buf_norm;
				}

#pragma omp critical
				{
					//todo: dynamic shadule?
					if (num == 0 && np > 1) //вместо критической секции пусть всегда первая нить управляет отправкой
					{
						// пересылаем первую пачку сообщений
#ifdef USE_CUDA			
						for (int i = 0; i < requests_rcv_1_section.size(); i++)
						{
							if (!flags_send_to_gpu_1_section[i])
							{
								MPI_Test(&requests_rcv_1_section[i], &flags_send_to_gpu_1_section[i], &status_rcv_1_section[i]); //проверяем все запросы принятия сообщения
								if (flags_send_to_gpu_1_section[i]) // если обмен завершён, но отправки не было
								{
									const int src = status_rcv_1_section[i].MPI_TAG;
									CudaSendIllumAsync(n_illum* size_1_section, disp_illum[src], grid.Illum);								
								}
							}
						}
#endif
					}

					if (np > 1)
					{
						//#pragma omp critical
						{
							MPI_Startall(np - 1, requests_send_2_section.data() + ((num_direction - size_1_section) * (np - 1)));
						}
					}
#ifdef USE_CUDA	
					if ((num == (num_of_th - 1)) && np > 1) // отличный от нулевого поток
					{
						for (int i = 0; i < requests_rcv_2_section.size(); i++)
						{
							if (!flags_send_to_gpu_2_section[i])
							{
								MPI_Test(&requests_rcv_2_section[i], &flags_send_to_gpu_2_section[i], &status_rcv_2_section[i]); //проверяем все запросы принятия сообщения
								if (flags_send_to_gpu_2_section[i]) // если обмен завершён, но отправки не было
								{
									CudaSendIllumAsync(n_illum, n_illum * (status_rcv_2_section[i].MPI_TAG), grid.Illum);  //переслать данные на gpu									
								}
							}
						}
					}

					//#pragma omp critical
					{
						CudaSendIllumAsync(n_illum, ((local_disp + num_direction) * n_illum), grid.Illum);
					}
#endif		
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
		timer.dir_time += omp_get_wtime();
				
		//todo: проверить можно ли здесь с блокировкой
		MPI_Request rq_norm;
		MPI_Iallgather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD, &rq_norm);
		

#ifdef USE_CUDA
		
		bool ready = true;
		timer.rcv1_wait_time = -omp_get_wtime();		
		do
		{			
			ready = true;			
			for (int i = 0; i < requests_rcv_1_section.size(); i++)
			{
				if (!flags_send_to_gpu_1_section[i])
				{
					ready = false;
					MPI_Test(&requests_rcv_1_section[i], &flags_send_to_gpu_1_section[i], &status_rcv_1_section[i]); //проверяем все запросы принятия сообщения
					if (flags_send_to_gpu_1_section[i]) // если обмен завершён, но отправки не было
					{
						int src = (status_rcv_1_section[i].MPI_TAG);
						CudaSendIllumAsync(n_illum* size_1_section, disp_illum[src], grid.Illum);
					}
				}
			}
		} while (!ready);
		timer.rcv1_wait_time += omp_get_wtime();
		
		timer.rcv2_wait_time = -omp_get_wtime();
		do
		{
			ready = true;
			for (int i = 0; i < requests_rcv_2_section.size(); i++)
			{
				if (!flags_send_to_gpu_2_section[i])
				{					
					ready = false;
					MPI_Test(&requests_rcv_2_section[i], &flags_send_to_gpu_2_section[i], &status_rcv_2_section[i]); //проверяем все запросы принятия сообщения
					if (flags_send_to_gpu_2_section[i]) // если обмен завершён, но отправки не было
					{
						CudaSendIllumAsync(n_illum, n_illum* (status_rcv_2_section[i].MPI_TAG), grid.Illum);  //переслать данные на gpu									
					}
				}
			}

		} while (!ready);
		timer.rcv2_wait_time += omp_get_wtime();
#else		
		timer.rcv1_wait_time = -omp_get_wtime();
		MPI_Waitall(requests_rcv_1_section.size(), requests_rcv_1_section.data(), MPI_STATUSES_IGNORE);
		timer.rcv1_wait_time += omp_get_wtime();
		
		timer.rcv2_wait_time = -omp_get_wtime();
		MPI_Waitall(requests_rcv_2_section.size(), requests_rcv_2_section.data(), MPI_STATUSES_IGNORE);
		timer.rcv2_wait_time += omp_get_wtime();
#endif
			

		timer.cuda_time = -omp_get_wtime();
		if (solve_mode.max_number_of_iter > 1)
		{
#ifdef USE_CUDA
			CalculateIntScatteringAsyncMPIStream(grid_direction, grid, 0, size_1_section, eCuda_scattering_1);
			CalculateIntScatteringAsyncMPIStream(grid_direction, grid, size_1_section, local_size, eCuda_scattering_2);						
			//CalculateIntScatteringAsyncMPI(grid_direction, grid, local_size);
#else

			CalculateIntCPU(count_cells, grid_direction, grid);
#endif
		}
		timer.cuda_time += omp_get_wtime();

		timer.norm_gather = -omp_get_wtime();
		MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
		for (auto n : norms) if (n > norm) norm = n;
		timer.norm_gather += omp_get_wtime();

		_clock += omp_get_wtime();

		if (myid == 0)
		{
			WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');
		}

		WRITE_LOG_MPI(count 
			<< ": phys= " << timer.phys_time
			<< ", rcv1 " << timer.rcv1_wait_time
			<< ", rcv2 " << timer.rcv2_wait_time
			<< ", send1 " << timer.send1_wait_time
			<< ", send2 " << timer.send2_wait_time
			<< ", dir " << timer.dir_time			
			//<< ", cuda " << timer.cuda_time
			<< ", cuda_wait1 " << timer.cuda_wait_time_1
			<< ", cuda_wait2 " << timer.cuda_wait_time_2
			<< ", norm_cast " << timer.norm_gather			
			<< ", all_step " << _clock << "\n\n", myid);

		count++;		

	} while (norm > solve_mode.accuracy && count < solve_mode.max_number_of_iter);


#ifndef USE_CUDA
	if (myid == 0) // это уйдет, когда все интегралы перейдут в cuda
	{
		ReCalcIllumGlobal(grid_direction.size, grid);
		//WriteFileSolution(BASE_ADRESS + "Illum_result.txt", Illum);
	}
#else
	CudaSyncStream(eCuda_params);
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
	for (auto& dir : grid_direction.directions)
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