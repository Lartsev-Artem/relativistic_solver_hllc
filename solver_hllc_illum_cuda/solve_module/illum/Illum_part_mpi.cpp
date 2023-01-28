#include "../solve_config.h"
#if defined ILLUM && defined SOLVE && defined USE_MPI

#include "illum_utils.h"

#include "../solve_global_struct.h"
#include "../../global_value.h"
#include "../../global_def.h"
#include "../solve_utils.h"


static Type CalculateIllumeOnInnerFace(const int neigh_id, Vector3& inter_coef)
{	
	//const int id_face_ = cell->geo.id_faces[num_in_face]; //����� �����
	//const int neigh_id = faces[id_face_].geo.id_r;  // ������� �� + ����� � ���������� ����������
	Type I_x0 = 0;

	switch (neigh_id)
	{
	case eBound_OutSource: // ��� ������
		inter_coef = Vector3(30, 30, 30);//cell->illum_val.coef_inter[num_in_face] = Vector3(30, 30, 30);
		return 30;
	case eBound_FreeBound:
		//cell->illum_val.coef_inter[num_in_face] = Vector3(0, 0, 0);
		inter_coef = Vector3(0, 0, 0);
		/*��������� �������*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;

	case eBound_LockBound:
		inter_coef = Vector3(0, 0, 0); //cell->illum_val.coef_inter[num_in_face] = Vector3(0, 0, 0);
		/*��������� �������*/
		//I_x0 = BoundaryFunction(num_cell, x, direction, illum_old, directions, squares);
		return I_x0;

	case eBound_InnerSource:  // ���������� �������	
	{
#if 0
		id_try_pos++;
		grid[num_cell].nodes_value[num_in_face] = Vector3(res_on_inner_bound, res_on_inner_bound, res_on_inner_bound);
		return res_on_inner_bound;

		Type data = res_inner_bound[ShiftRes + pos_in_res++]; // ������ �� ����� �� ���������??
		if (data >= 0) //������ �� ����������� � ������ ��� �����
		{
			/*
				 ��������� �������. ������ � ������� ������ ������ � ������������ � ������������� ������
				 �� ����������� � �����������
			*/
			return res_on_inner_bound; // I_x0;
			return data;
		}
		else // ������������� �������� ��������������� ����� (�������� ������ � ������ s=(x-x0).norm())
		{

			// ��������� �� ������ �������. ���� ����� ������������ ��������� ��� ������ ��������� ������� (��� == �������)

			//+dist_try_surface			
			int id = id_try_surface[ShiftTry + id_try_pos - 1];  // ����� ������ id �����			
			const int cell = id / 4;
			const int face = id % 4;

			Vector3 coef = grid[cell].nodes_value[face];
			Vector2	x0_local = X0[ShiftX0 + posX0++];//grid[num_cell].x0_loc[num_in_face_dir];

			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

			//if (I_x0 < 0) I_x0 = 0;
			return  res_on_inner_bound; // I_x0;

		}
#endif
		inter_coef = Vector3(30, 30, 30);//cell->illum_val.coef_inter[num_in_face] = Vector3(30, 30, 30);
		return 30;
	}

	default:
	{

		//Vector3 coef = grid[num_cell].nodes_value[num_in_face];
		Vector3 coef = inter_coef; // cell->illum_val.coef_inter[num_in_face];

		//������ ������ �������� � �� ������������ ������������

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
}

static Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const grid_directions_t& grid_direction) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
		Type dot = direction.dot(direction2);
	return (3. * (1 + dot* dot)) / 4;
	} };


	Vector3 cur_direction;
	Type S = 0;
	const int N_dir = grid_direction.size;
	const int N_cell = illum_old.size() / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		Type I = 0;
		for (int i = 0; i < base; i++)
		{
			I+=illum_old[num_direction * N_cell + num_cell + i];
		}
		I /= base;
		
		S += Gamma(grid_direction.directions[num_direction].dir, direction) * I * grid_direction.directions[num_direction].area;
	}
	return S / grid_direction.full_area;     // ���� *4PI, �� ��-�� ���������� Gamma ��������� �� 4PI
}

static int CalculateIntCPU(const int num_cells, const std::vector<Type>& illum, const grid_directions_t& grid_direction,
	vector<Type>& int_scattering)
{

#pragma omp parallel default(none) shared(num_cells, illum, grid_direction, int_scattering)
	{
		Vector3 direction;
		const int num_directions = grid_direction.size;
#pragma omp for
		for (int num_direction = 0; num_direction < num_directions; ++num_direction) {
			direction = grid_direction.directions[num_direction].dir;
			for (int cell = 0; cell < num_cells; cell++)
			{
				int_scattering[num_direction * num_cells + cell] = GetS(base * cell, direction, illum, grid_direction);
			}
		}
	}
	return 0;
}


static Type GetCurIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, const flux_t& cell)
{
	switch (solve_mode.class_vtk)
	{
	case 0: // ��� ��������� �����������  (���������� ���)
	{
		Type Q = 0;
		Type alpha = 2;
		Type betta = 1;
		Type S = int_scattering;// 0;

		if ((x - Vector3(1, 0, 0)).norm() > 0.09) { Q = 0; alpha = 0.5;  betta = 0.5; }

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

		Type betta = alpha / 2;  // ������ �� ������
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
		const Type ss = s * 388189 * 1e5;  // ����� --- ������� � ��������� ����������

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
		Type betta = alpha / 2;  // ������ �� ������
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

	case 10: // ��� ������ (�������, ��� ���������� ����� �� ���������� � ����������� �� �����������������)
	{
		Type Q = 0;
		Type alpha = 0.5;
		Type betta = 0.5;
		Type S = int_scattering;


		//if (x[0] < 0.06) // ���������� ����
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

	case 11: // HLLC + Illum ��� ������
	{

		Type S = int_scattering * RADIATION;

		Type d = cell.d * DENSITY;
		Type v = cell.v.norm() * VELOCITY;
		Type p = cell.p * PRESSURE;

		Type T = p / (d * R_gas);

		if (x[0] < 0.05)
		{
			//d = 0.1;
			//p = 0.01;
			//T = p / (d * R);
		}

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

		return I / RADIATION;
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
					if (buf_norm > norm_loc) norm_loc = buf_norm;
				}
				Illum[id] = curI;			
			}
			//cells[num_cell].illum_val.illum[num_dir] /= base; //�������� �� ����������� ��� ������� ������� � �.�.
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

static Type ReCalcIllumGlobal(const int dir_size, std::vector<elem_t>& cells, std::vector<Type>& Illum)
{
	const int size = cells.size();
	for (size_t num_dir = 0; num_dir < dir_size; num_dir++)
	{
		const int shift_dir = num_dir * size;
		for (int num_cell = 0; num_cell < size; num_cell++)
		{
			for (int i = 0; i < base; i++)
			{
				const int id = base * (shift_dir + num_cell) + i;
				cells[num_cell].illum_val.illum[num_dir * base + i] = Illum[id]; //�� ������ ����� �� ������������
			}
		}
	}
	return 0;
}
static int GetIntScattering(const int count_cells, const grid_directions_t& grid_direction,  std::vector<Type>& Illum, std::vector<Type>& int_scattering)
{
#ifdef USE_CUDA
	if (solve_mode.use_cuda)
	{
		if (CalculateIntScattering(32, count_cells, grid_direction.size, Illum, int_scattering))  // ���� ���� ������ ����������� �� CPU 
		{
			CalculateIntCPU(count_cells, Illum, grid_direction, int_scattering);
			ClearDevice(solve_mode.cuda_mod);
			solve_mode.use_cuda = false;
		}
	}
	else
#endif
	{
		CalculateIntCPU(count_cells, Illum, grid_direction, int_scattering);
	}

	return 0;
}

int CalculateIllum(const grid_directions_t& grid_direction, const std::vector< std::vector<int>>& face_states, const std::vector<int>& pairs,
	const std::vector < std::vector<cell_local>>& vec_x0, std::vector<BasePointTetra>& vec_x,
	const std::vector < std::vector<int>>& sorted_id_cell,
	//const std::vector<Type>& res_inner_bound, 
	grid_t& grid, std::vector<Type>& Illum, std::vector<Type>& int_scattering)
{
	// ����� ���� ����� ����� �� ���� �����
	const int count_directions = grid_direction.size;

	
	/*const*/ int count_cells = grid.size;

	MPI_Bcast(&count_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int np = 1, myid = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	std::vector<int> send_count;
	std::vector<int> disp;

	std::vector<int> send_count_illum(np,0);
	std::vector<int> disp_illum(np,0);

	std::vector<int> send_count_scattering(np,0);
	std::vector<int> disp_scattering(np,0);

	GetSend(np, count_directions, send_count);
	GetDisp(np, count_directions, disp);

	for (int i = 0; i < np; i++)
	{
		send_count_illum[i] = send_count[i] * base * count_cells;
		send_count_scattering[i] = send_count[i] * count_cells;

		disp_illum[i] = disp[i] * base * count_cells;
		disp_scattering[i] = disp[i] * count_cells;
	}

	const int local_size = send_count[myid];
	const int local_disp = disp[myid];

	/*
	�����:
	����� ��������� ��������� (��������� ����)
	illum after recalc -> �������� ����
	---
	����� �� ����� �������, �� �� ���������� ����� ���������� �� ������(phys_val + illum_val)
	� ��������� ����� ������� (illum_val)

	struct flux_t
	{
		Type d;
		Vector3 v;
		Type p;
	}

	struct illum_value_t
	{
		std::vector<Type> illum; //num_dir*base
		Type energy;
		Vector3 stream;

		Matrix3 impuls;
		Type div_stream;
		Vector3 div_impuls;
		Type absorp_coef;
		Type rad_en_loose_rate;
	}

	struct elem_t
	{
		flux_t  phys_val;
		flux_t  conv_val;
		illum_value_t illum_val;
	}
	*/

	MPI_Datatype MPI_flux_t;
	{
		int len[3 + 1] = { 1,3,1,  1 };
		MPI_Aint pos[4] = { 0,sizeof(Type), sizeof(Vector3) ,sizeof(flux_t) };
		MPI_Datatype typ[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_UB };
		MPI_Type_struct(4, len, pos, typ, &MPI_flux_t);
		MPI_Type_commit(&MPI_flux_t);
	}

	MPI_Datatype MPI_illum_t;
	{
		int len[2 + 1] = { 1,1, 1 };
		MPI_Aint pos[3] = { 0, sizeof(Type),sizeof(illum_value_t) };
		MPI_Datatype typ[3] = { MPI_DOUBLE, MPI_DOUBLE,MPI_UB };
		MPI_Type_struct(3, len, pos, typ, &MPI_illum_t);
		MPI_Type_commit(&MPI_illum_t);
	}
	MPI_Datatype MPI_elem_t;
	{
		int len[2 + 1] = { 1,1, 1 };
		MPI_Aint pos[3] = { 0, 2 * sizeof(flux_t),sizeof(elem_t) };
		MPI_Datatype typ[3] = { MPI_flux_t,MPI_illum_t, MPI_UB };
		MPI_Type_struct(3, len, pos, typ, &MPI_elem_t);
		MPI_Type_commit(&MPI_elem_t);
	}
	MPI_Datatype MPI_elem_t_phys;
	{
		int len[1 + 1] = { 1, 1 };
		MPI_Aint pos[2] = { 0, sizeof(elem_t) };
		MPI_Datatype typ[2] = { MPI_flux_t, MPI_UB };
		MPI_Type_struct(2, len, pos, typ, &MPI_elem_t_phys);
		MPI_Type_commit(&MPI_elem_t_phys);
	}


	printf("Run\n");

	int count = 0;
	Type norm = -1;
	
	std::vector<flux_t> phys_local(count_cells);

	static std::vector<Vector3> inter_coef(count_cells * base);

	static std::vector<Type> loc_illum(count_cells * base * local_size, 0); //����� id 0

	static std::vector<Type>int_scattering_local(count_cells * local_size);

	do {
		Type _clock = -omp_get_wtime();

		norm = -1;
		/*---------------------------------- ����� FOR �� ������������----------------------------------*/

		for (register int num_direction = 0; num_direction < local_size; ++num_direction)
		{
			int posX0 = 0;
			Vector3 I;
			/*---------------------------------- ����� FOR �� �������----------------------------------*/
			for (int h = 0; h < count_cells; ++h)
			{
				const int num_cell = sorted_id_cell[num_direction][h];

				//elem_t* cell = &grid.cells[num_cell];

				//sumI = 0;

				for (ShortId num_out_face = 0; num_out_face < base; ++num_out_face)
				{
					if (check_bit(face_states[num_direction][num_cell], num_out_face)) continue;

					//GetNodes
					for (int num_node = 0; num_node < 3; ++num_node)
					{
						Vector3 x = vec_x[num_cell].x[num_out_face][num_node];

						cell_local x0 = vec_x0[num_direction][posX0++];

						ShortId num_in_face = x0.in_face_id;
						Type s = x0.s;
						Vector2 X0 = x0.x0;

						Type I_x0 = CalculateIllumeOnInnerFace(pairs[num_cell * base + num_in_face], inter_coef[num_cell * base + num_in_face]);

						I[num_node] = GetCurIllum(x, s, I_x0, int_scattering_local[num_direction * count_cells + num_cell], phys_local[num_cell]);

					}//num_node

					inter_coef[num_cell * base + num_out_face] = I;
					const int id_face = pairs[num_cell * base + num_out_face];
					if (id_face >= 0)
					{
						inter_coef[id_face] = I;
					}

				} //num_out_face	

			}

			/*---------------------------------- ����� FOR �� �������----------------------------------*/
			norm = ReCalcIllum(num_direction, inter_coef, loc_illum);

		}
		/*---------------------------------- ����� FOR �� ������������----------------------------------*/
				

		MPI_Gatherv(loc_illum.data(), loc_illum.size(), MPI_DOUBLE, Illum.data(), send_count_illum.data(), disp_illum.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (myid == 0)		
		{
			if (solve_mode.max_number_of_iter > 1)  // ������� ������ ��������
			{
				GetIntScattering(count_cells, grid_direction, Illum, int_scattering);
			}

			_clock += omp_get_wtime();
			WRITE_LOG("Error:= " << norm << '\n' << "End iter_count number: " << count << " time= " << _clock << '\n');									
		}

		MPI_Barrier(MPI_COMM_WORLD); // ���� ������ ��������� ���������

		MPI_Scatterv(int_scattering.data(), send_count_scattering.data(), disp_scattering.data(), MPI_DOUBLE,
			int_scattering_local.data(), int_scattering_local.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		

		count++;
	} while (norm > solve_mode.accuracy && count < solve_mode.max_number_of_iter);

	// wait ��� ����
	if (myid == 0)
	{
		ReCalcIllumGlobal(grid_direction.size, grid.cells, Illum);
	}

	return 0;
}

//-----------------------------------------------------//


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

int CalculateIllumParam(const grid_directions_t& grid_direction, grid_t& grid) 
{
	MakeEnergy(grid_direction, grid.cells);
	MakeStream(grid_direction, grid);
	MakeDivImpuls(grid_direction, grid);

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


#endif //ILLUM