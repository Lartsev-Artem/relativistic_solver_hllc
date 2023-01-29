#include "../solve_config.h"
#if defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 && defined RHLLC_MPI
#include "../solve_global_struct.h"
#include "../../file_module/reader_bin.h"
#include "../../utils/grid_geometry/geometry_solve.h"

#define MPI_RETURN(a) EXIT(a)

class VectorVal
{
private:
	Type _val[base + 1];
public:
	VectorVal() {};
	VectorVal(const Type val) { for (int i = 0; i < base + 1; i++) _val[i] = val; }
	VectorVal(const Type a, const Type b, const Type c, const Type d, const Type e)
	{
		_val[0] = a; _val[1] = b; _val[2] = c; _val[3] = d; _val[4] = e;
	}
	VectorVal(const Type* data) { for (int i = 0; i < base + 1; i++) _val[i] = data[i]; }
	VectorVal(const VectorVal& v) { for (int i = 0; i < base + 1; i++) _val[i] = v[i]; };

	VectorVal& operator=(const VectorVal& x)
	{
		for (int i = 0; i < base + 1; i++) _val[i] = x._val[i];
		return *this;
	}
	Type operator [] (const int i) const { return _val[i]; }
	Type operator () (const int i) const { return _val[i]; }
	Type& operator [] (const int i) { return *(_val + i); }
	Type& operator () (const int i) { return *(_val + i); }

	VectorVal operator+(const VectorVal& val) const {
		return VectorVal(
			_val[0] + val[0],
			_val[1] + val[1],
			_val[2] + val[2],
			_val[3] + val[3],
			_val[4] + val[4]);
	}

	VectorVal operator-(const VectorVal& val) const {
		return VectorVal(
			_val[0] - val[0],
			_val[1] - val[1],
			_val[2] - val[2],
			_val[3] - val[3],
			_val[4] - val[4]);
	}

	VectorVal operator* (const Type x) const
	{
		return VectorVal(
			_val[0] * x,
			_val[1] * x,
			_val[2] * x,
			_val[3] * x,
			_val[4] * x);
	}

	VectorVal operator/ (const Type x) const
	{
		return VectorVal(
			_val[0] / x,
			_val[1] / x,
			_val[2] / x,
			_val[3] / x,
			_val[4] / x);
	}

	void operator+= (const Type x) { for (size_t i = 0; i < base + 1; i++) _val[i] += x; }
	
	void operator+= (const VectorVal& val) { for (size_t i = 0; i < base + 1; i++) _val[i] += val[i];	}

	inline int size() { return base + 1; }
	inline Type* data() { return _val; }
	inline void zero() { for (size_t i = 0; i < base + 1; i++) _val[i] = 0; }

	friend VectorVal operator*(const Type x, const VectorVal& val) 
	{
		return VectorVal(val[0] * x,
			val[1] * x,
			val[2] * x,
			val[3] * x,
			val[4] * x);
	}

	/// Да. Это жестко
	friend VectorVal operator*(const Eigen::MatrixXd& T, const VectorVal& val)
	{
		VectorX loc(5); loc << val[0], val[1], val[2], val[3], val[4];
		loc = T * loc;
		return VectorVal(loc[0],
			loc[1],
			loc[2],
			loc[3],
			loc[4]);
	}
};

void RHLLC_Flux(const VectorVal& W_R, const VectorVal& U_R, const VectorVal& W_L, const VectorVal& U_L, VectorVal& F) {

	/*
	An HLLC Riemann solver for relativistic flows – I. Hydrodynamics
	A. Mignone and G. Bodo
	INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
	Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
	Mon. Not. R. Astron. Soc. 364, 126–136 (2005)

	https://github.com/PrincetonUniversity/Athena-Cversion/blob/master/src/rsolvers/hllc_sr.c

	\note: комбинация кода из mignone 2005 и 2006. в части hllc
	*/

	//==================== Кэшируем физические переменные слева и справа============================//
	// нормальная сокорость
	const Vector3 Vel_L(W_L[1], W_L[2], W_L[3]);  //T * velocity[num_cell];
	const Vector3 Vel_R(W_R[1], W_R[2], W_R[3]);  //T * velocity[neig / 4];

	const Type d_L = W_L(0);
	const Type d_R = W_R(0);

	const Type p_L = W_L(4);
	const Type p_R = W_R(4);

	const Type VV_L = Vel_L.dot(Vel_L);
	const Type VV_R = Vel_R.dot(Vel_R);

	//========================================================================================//

	//=========================Вычисляем релятивистикие параметры============================//				
	const Type g_L = 1. / sqrt(1 - VV_L);	// фактор Лоренца
	const Type g_R = 1. / sqrt(1 - VV_R);

	const Type h_L = 1 + gamma_g * p_L / d_L; // энтальпия
	const Type h_R = 1 + gamma_g * p_R / d_R;

	const Type cs_L = sqrt((gamma1 * p_L) / (d_L * h_L)); // скорость звука
	const Type cs_R = sqrt((gamma1 * p_R) / (d_R * h_R));

	const Type sigmaS_L = (cs_L * cs_L) / (g_L * g_L * (1 - cs_L * cs_L)); // что-то для расчета собственных чисел HHL
	const Type sigmaS_R = (cs_R * cs_R) / (g_R * g_R * (1 - cs_R * cs_R));

	//========================================================================================//

	const Type sqr_L = sqrt(sigmaS_L * (1 - Vel_L[0] * Vel_L[0] + sigmaS_L));
	const Type sqr_R = sqrt(sigmaS_R * (1 - Vel_R[0] * Vel_R[0] + sigmaS_R));

	// здесь встречалась альтернатива сравнения с нулем min(0,L), max(0,R)
	const Type lambda_L = min((Vel_L[0] - sqr_L) / (1 + sigmaS_L), (Vel_R[0] - sqr_R) / (1 + sigmaS_R));
	const Type lambda_R = max((Vel_L[0] + sqr_L) / (1 + sigmaS_L), (Vel_R[0] + sqr_R) / (1 + sigmaS_R));

	if (lambda_R <= 0) // если верно выполнить всегда
	{
		F(0) = U_R[0] * Vel_R[0]; //D*v_x
		F(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
		F(2) = U_R[2] * Vel_R[0];
		F(3) = U_R[3] * Vel_R[0];
		F(4) = U_R[1];
		//continue;			
	}
	else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
	{
		F(0) = U_L[0] * Vel_L[0]; //D*v_x
		F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
		F(2) = U_L[2] * Vel_L[0];
		F(3) = U_L[3] * Vel_L[0];
		F(4) = U_L[1];
		//continue;			
	}
	else
	{
		//====================Расчёт потоков и приближений hll=========================================//
		VectorVal F_L(5);
		VectorVal F_R(5);
		VectorVal U_hll(5);
		VectorVal F_hll(5);

		F_R(0) = U_R[0] * Vel_R[0]; //D*v_x
		F_R(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
		F_R(2) = U_R[2] * Vel_R[0];
		F_R(3) = U_R[3] * Vel_R[0];
		F_R(4) = U_R[1];

		F_L(0) = U_L[0] * Vel_L[0]; //D*v_x
		F_L(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
		F_L(2) = U_L[2] * Vel_L[0];
		F_L(3) = U_L[3] * Vel_L[0];
		F_L(4) = U_L[1];

		//	cout << "F_L\n" << F_L << "\nF_R\n" << F_R << '\n';
		//	cout << "U_L\n" << U_L << "\nU_R\n" << U_R << '\n';			
		F_hll = (lambda_R * F_L - lambda_L * F_R + (lambda_R * lambda_L * (U_R - U_L))) / (lambda_R - lambda_L);
		U_hll = (lambda_R * U_R - lambda_L * U_L + (F_L - F_R)) / (lambda_R - lambda_L);

		//	cout << "F_hll\n" << F_hll << "\nU_hll\n" << U_hll << '\n';
#ifdef ONLY_RHLL
		F = F_hll;
#endif

		//============================================================================================//
#ifndef ONLY_RHLL		
//=========================Поиск скорости промежуточной волны===============================//
		const Type a = F_hll[4];			//F_E^hll
		const Type b = -U_hll[4] - F_hll[1]; // (E_hll + F_mx^hll)
		const Type c = U_hll[1];			//mx_hll

#if 1 // как описано в Mignone...
		Type quad = -0.5 * (b + SIGN(b) * sqrt(b * b - 4 * a * c));
		Type _lambda = c / quad;

#endif		

		if (_lambda >= 0.0)
		{
			//============================Поиск промежуточного давления ===================================//
			const Type _p = -F_hll[4] * _lambda + F_hll[1];
			//============================================================================================//

			//==========================Финальный поток HLLC=============================================//
			VectorVal _U_L;
			const Type dif_L = 1.0 / (lambda_L - _lambda);

			_U_L[0] = (U_L[0] * (lambda_L - Vel_L[0])) * dif_L;
			_U_L[1] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
			_U_L[2] = (U_L[2] * (lambda_L - Vel_L[0])) * dif_L;
			_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0])) * dif_L;
			_U_L[4] = (U_L[4] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

			F = F_L + lambda_L * (_U_L - U_L);

			//============================================================================================//
		}
		else //(_S <= 0)
		{
			//============================Поиск промежуточного давления ===================================//
			const Type _p = -F_hll[4] * _lambda + F_hll[1];
			//============================================================================================//
			VectorVal _U_R;
			const Type dif_R = 1.0 / (lambda_R - _lambda);

			_U_R[0] = (U_R[0] * (lambda_R - Vel_R[0])) * dif_R;
			_U_R[1] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
			_U_R[2] = (U_R[2] * (lambda_R - Vel_R[0])) * dif_R;
			_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0])) * dif_R;
			_U_R[4] = (U_R[4] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

			F = F_R + lambda_R * (_U_R - U_R);
		}
	}
#endif
	return;
}


inline int GetLocalId(const int neighb_glob, const int shift_node)
{
	return neighb_glob >= 0 ? (((neighb_glob / base) - shift_node) * base + neighb_glob % base) : neighb_glob;

	//const int neigh = neighbours_id_faces_local[base * i + k];	
	//if (neighb_glob >= 0)
	//{
	//	/*const int neigh_cell = neighb_glob / base;
	//	const int new_idx_cell = neigh_cell - shift_node;
	//	return new_idx_cell * base + neighb_glob % base;*/		
	//	return ((neighb_glob / base) - shift_node) * base + neighb_glob % base;
	//}
	//else
	//{
	//	return neighb_glob;
	//}
}

int  RHLLC_MPI(std::string& main_dir,
	std::vector<Vector3>& centerts, std::vector<int>& neighbours_id_faces,
	std::vector<Normals>& normals, std::vector<Type>& squares_cell, std::vector<Type>& volume)
{
	int np, myid;

	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef WRITE_LOG
	{
		remove((main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt").c_str());
	}
#endif

	MPI_Datatype MPI_VectorVal;
	int len[2] = { base + 1, 1 };
	MPI_Aint pos[2] = { 0 ,sizeof(VectorVal) };
	MPI_Datatype typ[2] = { MPI_DOUBLE, MPI_UB };
	MPI_Type_struct(2, len, pos, typ, &MPI_VectorVal);
	MPI_Type_commit(&MPI_VectorVal);
	//
	printf("Run\n");

	const int size_grid = centerts.size();

	std::vector<VectorVal> W_full_3d;
	std::vector<VectorVal> U_full_3d;
	std::vector<VectorVal> U_full_3d_prev;

	std::vector<VectorVal> W_full_3d_local;
	std::vector<VectorVal> U_full_3d_local;
	std::vector<VectorVal> U_full_3d_prev_local;

	std::string name_file_bound_cur = main_dir + "bound_cells_mpi" + to_string(myid) + "_cur.bin";
	std::string name_file_bound_next = main_dir + "bound_cells_mpi" + to_string(myid) + "_next.bin";
	std::string name_file_cells = main_dir + "cells_node_mpi" + to_string(myid) + ".bin";

	std::vector<int> id_cells_local;		 // регулярные ячейки на узле
	std::vector<int> id_bound_cur_local;  // граница на этом узле
	std::vector<int> id_bound_next_local; // граница c другого

	ReadSimpleFileBin(name_file_cells, id_cells_local);
	ReadSimpleFileBin(name_file_bound_cur, id_bound_cur_local);
	ReadSimpleFileBin(name_file_bound_next, id_bound_next_local);
	 int local_size = id_cells_local.size();
	 int local_size_cur = id_bound_cur_local.size();

	const int size = local_size + local_size_cur;
	std::vector<int> neighbours_id_faces_local(base * size);
	std::vector<Normals> normals_local(size);
	std::vector<Type> squares_cell_local(base * size);
	std::vector<Type> volume_local(size);

#if 0
	for (size_t i = 0; i < local_size; i++)
	{
		int id = id_cells_local[i];
		volume_local[i] = volume[id];
		normals_local[i] = normals[id];
		for (size_t j = 0; j < base; j++)
		{
			squares_cell_local[i * base + j] = squares_cell[id * base + j];
			neighbours_id_faces_local[i * base + j] = neighbours_id_faces[id * base + j];  // здесь пересчёт
		}
	}
	for (size_t k = 0; k < local_size_cur; k++)
	{
		int i = k + local_size;// границу в конец

		int id = id_bound_cur_local[i];
		volume_local[i] = volume[id];
		normals_local[i] = normals[id];
		for (size_t j = 0; j < base; j++)
		{
			squares_cell_local[i * base + j] = squares_cell[id * base + j];
			neighbours_id_faces_local[i * base + j] = neighbours_id_faces[id * base + j];  // здесь пересчёт
		}
	}
	neighbours_id_faces.clear();
	normals.clear();
	squares_cell.clear();
	volume.clear();

	//re index
	{
		//*Здесь мы теряем связь с нумерацийе глобальной границы!!!!*//
		int loc_id_cell = 0;
		for (int i = 0; i < local_size/*+local_size_cur*/; i++)  //в рамках узла
		{
			auto min_el = std::min_element(id_bound_cur_local.begin(), id_bound_cur_local.end());
			int shift_node = min(id_cells_local[0], (int)std::distance(id_bound_cur_local.begin(), min_el));//может не работать на с++11. Если что написать свой find_min
			for (int k = 0; k < base; k++)
			{
				const int neigh = neighbours_id_faces_local[base * i + k];
				const int neigh_cell = neigh / base;
				if (neigh >= 0)
				{
					const int new_idx_cell = neigh_cell - shift_node;
					const int new_idx_face = new_idx_cell * base + neigh % base;
					neighbours_id_faces_local[base * i + k] = new_idx_face;
				}
				else
				{
					neighbours_id_faces_local[base * i + k] = neigh;
				}
			}
		}
	}


	if (myid == 0)
	{
		// формируем массивы на управляющем узле

		RHLLC_Init_3d(size_grid, centerts, W_full_3d);

		U_full_3d.resize(size_grid);
		U_full_3d_prev.resize(size_grid);

		ReBuildConvValue_3d(W_full_3d, U_full_3d_prev);
	}

#endif

	std::vector<int> send_count;
	std::vector<int> disp;
	GetSend(np, size_grid, send_count);
	GetDisp(np, size_grid, disp);

	local_size = send_count[myid];
	const int local_disp = disp[myid];

	std::vector<int> pairs_local(local_size * base);  // соседи на узле в глобальной нумерации

	std::vector<int> pairs_local_reg(local_size * base);  // соседи на узле регулярные в лок
	std::vector<int> pairs_local_irr_left;  // соседи на узле с другого узла в глоб
	std::vector<int> pairs_local_irr_right;  // соседи на узле с другого узла в глоб

	std::vector<int> irr_left_id;
	std::vector<int> irr_right_id;

	int left_s = 0;
	int right_s = 0;
	int reg_s = 0;

	bool reg_cell = true;
	bool irr_l_cell = false;
	bool irr_r_cell = false;

	for (size_t i = 0; i < local_size; i++)
	{
		reg_cell = true; irr_l_cell = false; irr_r_cell = false;

		for (size_t j = 0; j < base; j++)
		{
			int neig = neighbours_id_faces[(local_disp + i) * base + j];
			pairs_local[i * base + j] = neig;
			if (neig < 0)
			{
				pairs_local_reg[i * base + j] = neig; // ГУ
			}
			else
			{
				int cell_id = neig / base;
				if (cell_id >= local_disp && cell_id < local_disp + local_size)  // <= or <??
				{
					pairs_local_reg[i * base + j] = GetLocalId(neig, local_disp);
				}
				else
				{
					reg_cell = false;
					if (cell_id < local_disp)
					{
						pairs_local_irr_left.push_back(neig);
						pairs_local_reg[i * base + j] = -10; // признак на соседний узел
						irr_l_cell = true;
					}
					else
					{
						pairs_local_irr_right.push_back(neig);
						pairs_local_reg[i * base + j] = -20; // признак на соседний узел
						irr_r_cell = true;
					}
				}
			}
		}

		if (irr_l_cell)
		{
			irr_left_id.push_back(i);
		}
		if (irr_r_cell)
		{
			irr_right_id.push_back(i);
		}
	}

	std::vector<VectorVal> U_prev_local_irr_left;
	std::vector<VectorVal> W_local_irr_left;

	std::vector<VectorVal> U_prev_local_irr_right;
	std::vector<VectorVal> W_local_irr_right;


	std::vector<VectorVal> W_local(local_size);  // в нем есть регулярные левые и правые
	std::vector<VectorVal> U_local(local_size);  // в нем есть регулярные левые и правые
	std::vector<VectorVal> U_prev_local(local_size);  // в нем есть регулярные левые и правые

	Eigen::MatrixXd T(5, 5);
	VectorVal tU;

	VectorVal F;
	VectorVal SumF;

	VectorVal U, U_L, U_R;
	VectorVal W, W_L, W_R;

	int cnt_left = 0;
	int cnt_right = 0;

	Type tau = 0.1;
	for (int num_cell = 0; num_cell < local_size; num_cell++)
	{
		SumF.zero();
		U = U_prev_local[num_cell];
		W = W_local[num_cell];

		for (int i = 0; i < base; i++)
		{
			const int neig = pairs_local_reg[base * num_cell + i];

			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;
				W_R = W;
				break;

			case -10:
				U_R = U_prev_local_irr_left[cnt_left];
				W_R = W_local_irr_left[cnt_left++];
				break;

			case -20:
				U_R = U_prev_local_irr_right[cnt_right];
				W_R = W_local_irr_right[cnt_right++];
				break;

			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				U_R = U_prev_local[neig / base];
				W_R = W_local[neig / base];
				break;
			}

			MakeRotationMatrix(normals_local[num_cell].n[i], T);

			U_L = T * U;   U_R = T * U_R;
			W_L = T * W;   W_R = T * W_R;

			RHLLC_Flux(W_R, U_R, W_L, U_L, F);

			F = (T.transpose()) * F;
			SumF += (F * squares_cell_local[base * num_cell + i]);

		}

		U_local[num_cell] = (U - SumF * tau / volume_local[num_cell]);
	}


	if (myid == 0)
	{
		if (np > 1)
		{
			int i = 0;
			W_local_irr_right.resize(irr_right_id.size());
			U_prev_local_irr_right.resize(irr_right_id.size());
			for (auto id : irr_right_id)
			{
				W_local_irr_right[i++] = W_local[id];
				U_prev_local_irr_right[i++] = U_prev_local[id];
			}

			MPI_Send(W_local_irr_right.data(), W_local_irr_right.size(), MPI_VectorVal, myid + 1, 0, MPI_COMM_WORLD);
		}
	}
	else if (myid == (np - 1))
	{
		int i = 0;
		W_local_irr_left.resize(irr_left_id.size());
		U_prev_local_irr_left.resize(irr_right_id.size());
		for (auto id : irr_left_id)
		{
			W_local_irr_left[i++] = W_local[id];
			U_prev_local_irr_left[i++] = U_prev_local[id];
		}

		MPI_Status st;
		MPI_Recv(W_local_irr_left.data(), W_local_irr_left.size(), MPI_VectorVal, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
	}
	else
	{
		MPI_Status st;
		MPI_Sendrecv(W_local_irr_right.data(), W_local_irr_right.size(), MPI_VectorVal, myid + 1, 0,
			W_local_irr_left.data(), W_local_irr_left.size(), MPI_VectorVal, myid - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
	}

	U_local.swap(U_prev_local);
	MPI_RETURN(0);

	//--------------------------------------------до while(T)_--------------------
}

#endif