#include "../solve_config.h"
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 3 && !defined USE_MPI && defined SOLVE 
#include "../solve_global_struct.h"
#include "../../file_module/reader_bin.h"
#include "../../utils/grid_geometry/geometry_solve.h"

#ifdef OLD_CLASS
static VectorX HLLC_stepToOMP(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_face_ts, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev) {

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
	Eigen::VectorXd tU;

	bool flag_bound = 1; // граничная ячейка или нет	

	for (size_t i = 0; i < 4; i++) // по граням
	{
		const int neig = neighbours_id_face_ts[4 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;

#if defined Cone && defined Jet // условие разрежения. См как заданы границы в make_struct.prj

				d = U_R[0];
				vel << U_R[1] / d, U_R[2] / d, U_R[3] / d;
				v = vel.norm();
				pressure = (U_R[4] - v * v * d / 2.) * (gamma1 - 1);

				d *= 0.5;
				vel[0] *= 1.1;

				v = vel.dot(vel);
				U_R[0] = d;
				U_R[1] = d * vel[0];
				U_R[2] = d * vel[1];
				U_R[3] = d * vel[2];
				U_R[4] = pressure / (gamma1 - 1) + d * v / 2;
#endif
				break;
			case eBound_InnerSource:
				U_R = U;
				break;
			case eBound_OutSource:
				MakeRotationMatrix(normals[num_cell].n[i], T);
				tU = T * U;
				tU[1] = -tU[1];
				U_R = (T.transpose()) * tU;

				break;
			case eBound_LockBound:

				MakeRotationMatrix(normals[num_cell].n[i], T);
				tU = T * U;
				tU[1] = -tU[1];
				U_R = (T.transpose()) * tU;

				break;
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				flag_bound = 0;
				U_R = U_full_prev[neig / 4];
				break;
			}
		}
#if 0
		if (neig >= 0)
		{
			flag_bound = 0;
			U_R = U_full_prev[neig / 4];
		}
		else
		{
			flag_bound = 1;
			if (neig == -1) // внешняя  граница
			{
				if (normals[num_cell].n[i].dot(Vector3(1, 0, 0)) > 0)
				{
					Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
					MakeRotationMatrix(normals[num_cell].n[i], T);
					Eigen::VectorXd tU = T * U;
					tU[1] = -tU[1];
					U_R = (T.transpose()) * tU;
				}
				else
					U_R = U;

			}
			else
				U_R = U;
		}
#endif			
		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????
		U_L = T * U;
		U_R = T * U_R;

#if 0 // 2d
		auto it_r = find(b_right.begin(), b_right.end(), num_cell);
		auto it_l = find(b_left.begin(), b_left.end(), num_cell);
		if (neig < 0) // если "какая-то" граница
			if (fabs((normals[num_cell].n[i] - Vector3(0, 1, 0)).norm()) > eps)
				if (it_r == b_right.end() && it_l == b_left.end()) U_R[1] = 0;  // если не правая и не левая граница, то условие стенки
#endif

		Type d_L = U_L[0]; //density[num_cell];
		Type d_R = U_R[0]; // density[neig];

		Vector3 vel(U_L(1) / d_L, U_L(2) / d_L, U_L(3) / d_L);
		Type v = vel.dot(vel);
		Type p_L = (U_L(4) - v * d_L / 2.) * (gamma1 - 1);

		vel << U_R(1) / d_R, U_R(2) / d_R, U_R(3) / d_R;
		v = vel.dot(vel);
		Type p_R = (U_R(4) - v * d_R / 2.) * (gamma1 - 1);

#ifdef SAVE_DUMP_HLLC
		if (p_L < 0 || p_R < 0 || d_L < 0 || d_R < 0)
		{
			printf("cell %d\n", num_cell);
			printf("%lf  %lf  %lf  %lf\n", p_L, p_R, d_L, d_R);

			printf("%d %d %d %d", count1, count2, count3, count4);
			bad_hllc_flag = true;
		}
#endif

		//const Type p_L = pressure[num_cell]; 
		//const Type p_R = pressure[neig];

		const Type v_L = U_L[1] / d_L; // sqrt(U_L[1] * U_L[1] + U_L[2] * U_L[2] + U_L[3] * U_L[3]);  //velocity[num_cell].norm();
		const Type v_R = U_R[1] / d_R; //sqrt(U_R[1] * U_R[1] + U_R[2] * U_R[2] + U_R[3] * U_R[3]); //velocity[neig].norm();

		const Type a_L = sqrt(gamma1 * p_L / d_L);
		const Type a_R = sqrt(gamma1 * p_R / d_R);

		const Type P = (p_L + p_R) / 2;
		const Type Den = (d_L + d_R) / 2;
		const Type A = (a_L + a_R) / 2;
		const Type _p = max(0.0, P - (v_R - v_L) * Den * A / 2);

		//const Type _p = max(0.0, ((p_L + p_R) - ((v_L - v_R) * (a_L + a_R) * (d_L + d_R) / 4.0)) / 2.0);

		// pressure-based wave speed estimates
		Type q_L = 1;
		const Type G = (gamma1 + 1) / (2 * gamma1);
		if (_p > p_L) q_L = sqrt(1 + G * (_p / p_L - 1));

		Type q_R = 1;
		if (_p > p_R) q_R = sqrt(1 + G * (_p / p_R - 1));

		const Type S_L = v_L - a_L * q_L;
		const Type S_R = v_R + a_R * q_R;

#ifdef DBG_OUTPUT
		printf("(%d), S_L=%lf,  S_R=%lf\n", num_cell, S_L, S_R);
#endif
		flag_bound = false;
		if (S_R <= 0) // если верно выполнить всегда
		{
			F(0) = U_R[1];
			F(1) = U_R[1] * U_R[1] / d_R + p_R;
			F(2) = U_R[1] * U_R[2] / d_R;
			F(3) = U_R[1] * U_R[3] / d_R;
			F(4) = (U_R[4] + p_R) * U_R[1] / d_R;
			//continue;
			c0++;
		}
		else if (S_L >= 0 || flag_bound) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[1];
			F(1) = U_L[1] * U_L[1] / d_L + p_L;
			F(2) = U_L[1] * U_L[2] / d_L;
			F(3) = U_L[1] * U_L[3] / d_L;
			F(4) = (U_L[4] + p_L) * U_L[1] / d_L;  // TU[4]*d_L
			//continue;
			c1++;
		}
		else
		{
			const Type roS_L = d_L * (S_L - v_L);
			const Type roS_R = d_R * (S_R - v_R);
			const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

			const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;

			Eigen::VectorXd D(5); D << 0, 1, 0, 0, _S;
		
			if (_S >= 0)
			{
				F(0) = U_L[1];
				F(1) = U_L[1] * U_L[1] / d_L + p_L;
				F(2) = U_L[1] * U_L[2] / d_L;
				F(3) = U_L[1] * U_L[3] / d_L;
				F(4) = (U_L[4] + p_L) * U_L[1] / d_L;

				VectorX buf = F;
				F = (_S * (S_L * U_L - buf) + S_L * P_LR * D) / (S_L - _S);
				c2++;
			}
			else //(_S <= 0)
			{
				F(0) = U_R[1];
				F(1) = U_R[1] * U_R[1] / d_R + p_R;
				F(2) = U_R[1] * U_R[2] / d_R;
				F(3) = U_R[1] * U_R[3] / d_R;
				F(4) = (U_R[4] + p_R) * U_R[1] / d_R;

				VectorX buf = F;
				F = (_S * (S_R * U_R - buf) + S_R * P_LR * D) / (S_R - _S);
				c3++;
			}
		}

#ifdef DBG_OUTPUT

		std::cout << "F= " << F << '\n';

		ofstream ofile;
		const std::string name_file = "D:\\Desktop\\FilesCourse\\dbg_hllc.txt";
		static bool start = false;
		if (!start)
		{
			start = true;
			ofile.open(name_file);
			ofile.close();
		}

		//if (fabs(F(0)) > 1e-1)
		{
			ofile.open(name_file, std::ios::app);
			if (!ofile.is_open())
			{
				printf("No open file\n");
				exit(1);
			}
			int i = num_cell;

			ofile << "S_L[" << num_cell << "]= " <<S_L  << " S_R= " << S_R << '\n';
			ofile << "U[" << num_cell << "]= " << U(0) << ", " << U(1) << ", " << U(4) << '\n';
			ofile << "F[" << num_cell << "]= " << F(0) << ", " << F(1) << ", " << F(4) << '\n';

			ofile.close();
		}
#endif

		VectorX buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[4 * num_cell + i];

	}// for


	/*U_full[num_cell] =*/ return (U - SumF * tau / volume[num_cell]);
}


int HLLC_3d(const int N, const Type tau, const std::vector<int>& bound_cells, const std::vector<int>& neighbours_id_face_ts,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	const std::vector<VectorX>& U_full_prev, const std::vector<VectorX>& U_full) {

#pragma omp parallel default(none) shared(N, tau, neighbours_id_face_ts, normals, squares_cell, volume, U_full_prev, U_full)
		{
			VectorX buf(5);
#pragma omp for
			for (int i = 0; i < N; i++)
			{
				buf = HLLC_stepToOMP(i, tau, neighbours_id_face_ts, normals, squares_cell, volume, U_full_prev);

				//buf[1] = Vector3(buf[1], buf[2], buf[3]).norm();
				//buf[2] = 0;
				//buf[3] = 0;
#pragma omp critical
				{
					U_full[i] = buf;
					// Для задачи SODA:
						//U_full[i][2] = 0;
						//U_full[i][3] = 0;
				}
			}
		}		
		return 0;
}

#else //NewCLASS

static int flux_t_calc(const flux_t& val_l, const flux_t& val_r, face_t& f)
{
	Matrix3 T;
	MakeRotationMatrix(f.geo.n, T);
	
	flux_t U_L = val_l; 
	U_L.v = T * val_l.v;	
	
	flux_t U_R= val_r;  
	U_R.v = T * val_r.v;

	Type d_L = U_L[0]; //density[num_cell];
	Type d_R = U_R[0]; // density[neig];

	Vector3 vel(U_L(1) / d_L, U_L(2) / d_L, U_L(3) / d_L);
	Type v = vel.dot(vel);
	Type p_L = (U_L(4) - v * d_L / 2.) * (gamma1 - 1);

	vel << U_R(1) / d_R, U_R(2) / d_R, U_R(3) / d_R;
	v = vel.dot(vel);
	Type p_R = (U_R(4) - v * d_R / 2.) * (gamma1 - 1);


	const Type v_L = U_L[1] / d_L; // sqrt(U_L[1] * U_L[1] + U_L[2] * U_L[2] + U_L[3] * U_L[3]);  //velocity[num_cell].norm();
	const Type v_R = U_R[1] / d_R; //sqrt(U_R[1] * U_R[1] + U_R[2] * U_R[2] + U_R[3] * U_R[3]); //velocity[neig].norm();

	const Type a_L = sqrt(gamma1 * p_L / d_L);
	const Type a_R = sqrt(gamma1 * p_R / d_R);

	const Type P = (p_L + p_R) / 2;
	const Type Den = (d_L + d_R) / 2;
	const Type A = (a_L + a_R) / 2;
	const Type _p = max(0.0, P - (v_R - v_L) * Den * A / 2);
	
	// pressure-based wave speed estimates
	Type q_L = 1;
	const Type G = (gamma1 + 1) / (2 * gamma1);
	if (_p > p_L) q_L = sqrt(1 + G * (_p / p_L - 1));

	Type q_R = 1;
	if (_p > p_R) q_R = sqrt(1 + G * (_p / p_R - 1));

	const Type S_L = v_L - a_L * q_L;
	const Type S_R = v_R + a_R * q_R;	


	flux_t F;
	if (S_R <= 0) // если верно выполнить всегда
	{
		F.d = U_R[1];
		F.v(0) = U_R[1] * U_R[1] / d_R + p_R;
		F.v(1) = U_R[1] * U_R[2] / d_R;
		F.v(2) = U_R[1] * U_R[3] / d_R;
		F.p = (U_R[4] + p_R) * U_R[1] / d_R;
		//continue;		
	}
	else if (S_L >= 0) // выполнить либо по условию либо для всех границ
	{
		F.d = U_L[1];
		F.v(0) = U_L[1] * U_L[1] / d_L + p_L;
		F.v(1) = U_L[1] * U_L[2] / d_L;
		F.v(2) = U_L[1] * U_L[3] / d_L;
		F.p = (U_L[4] + p_L) * U_L[1] / d_L;  // TU[4]*d_L
		//continue;		
	}
	else
	{
		const Type roS_L = d_L * (S_L - v_L);
		const Type roS_R = d_R * (S_R - v_R);
		const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

		const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;
		
		flux_t D(0, 1, 0, 0, _S);

		if (_S >= 0)
		{
			F.d = U_L[1];
			F.v(0) = U_L[1] * U_L[1] / d_L + p_L;
			F.v(1) = U_L[1] * U_L[2] / d_L;
			F.v(2) = U_L[1] * U_L[3] / d_L;
			F.p = (U_L[4] + p_L) * U_L[1] / d_L;

			flux_t buf = F;
			F = ((U_L * S_L - buf) *_S + D * P_LR * S_L) / (S_L - _S);			
		}
		else //(_S <= 0)
		{
			F.d = U_R[1];
			F.v(0) = U_R[1] * U_R[1] / d_R + p_R;
			F.v(1) = U_R[1] * U_R[2] / d_R;
			F.v(2) = U_R[1] * U_R[3] / d_R;
			F.p = (U_R[4] + p_R) * U_R[1] / d_R;

			flux_t buf = F;
			F = ( ( U_R* S_R - buf) *_S +  D* S_R * P_LR) / (S_R - _S);			
		}
	}

	f.f = F;
	f.f.v = (T.transpose()) * F.v;
	f.f = f.f * f.geo.S;

	return 0;
}

static void hllc_get_phys_value(const flux_t& U, flux_t& W)
{
	const Type d = U.d;
	W.d = d;
	W.v = U.v / d;
	const Type vv = W.v.dot(W.v);
	W.p = (U.p - vv * d / 2.) * (gamma1 - 1);
}

int HLLC_3d(const Type tau, grid_t& grid)
{
	flux_t bound_val;
	MatrixX T(5, 5);
	VectorX U(5);
	elem_t* cell;
	// потоки
	for(auto &f : grid.faces)
	{			
		cell = &grid.cells[f.geo.id_l];
		switch (f.geo.id_r)// id соседа она же признак ГУ
		{
		case eBound_FreeBound:
			bound_val = cell->conv_val;
			break;
		case eBound_InnerSource:
			bound_val = cell->conv_val;
			break;
		case eBound_OutSource:			
#if 0
			Matrix3 TT;
			MakeRotationMatrix(f.geo.n, TT);
			Vector3 UU = cells[f.geo.id_l].conv_val.v;
			UU = T * UU;
			UU[1] = -UU[1];
			UU = (TT.transpose()) * UU;
			bound_val.d = cells[f.geo.id_l].conv_val.d;
			bound_val.v = UU;
			bound_val.p = cells[f.geo.id_l].conv_val.p;
#endif
			MakeRotationMatrix(f.geo.n, T);
			U << grid.cells[f.geo.id_l].conv_val.d, grid.cells[f.geo.id_l].conv_val.v(0),
				grid.cells[f.geo.id_l].conv_val.v(1), grid.cells[f.geo.id_l].conv_val.v(2),
				grid.cells[f.geo.id_l].conv_val.p;
			U = T * U;
			U[1] = -U[1];
			U = (T.transpose()) * U;

			bound_val.d = cell->conv_val.d;
			bound_val.v(0) = U(1);
			bound_val.v(1) = U(2);
			bound_val.v(2) = U(3);			
			bound_val.p = cell->conv_val.p;

			break;
		case eBound_LockBound:

			MakeRotationMatrix(f.geo.n, T);
			U << cell->conv_val.d,
				cell->conv_val.v(0),
				cell->conv_val.v(1),
				cell->conv_val.v(2),
				cell->conv_val.p;
			U = T * U;
			U[1] = -U[1];
			U = (T.transpose()) * U;

			bound_val.d = cell->conv_val.d;
			bound_val.v(0) = U(1);
			bound_val.v(1) = U(2);
			bound_val.v(2) = U(3);
			bound_val.p = cell->conv_val.p;

			break;

		default:
			if (f.geo.id_r < 0)
			{
				printf("Err bound in HLLC_3d\n");
				EXIT(1);				
			}

			bound_val = grid.cells[f.geo.id_r].conv_val;
			break;
		}

		flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val, f);
	}

	// ячейки
	for (auto& el : grid.cells)
	{		
		flux_t sumF;
		for (int j = 0; j < base; j++)
		{
			if (el.geo.sign_n[j])
			{
				sumF += grid.faces[el.geo.id_faces[j]].f;
			}
			else
			{
				sumF -= grid.faces[el.geo.id_faces[j]].f;
			}
		}
		el.conv_val -= sumF * (tau / el.geo.V);
	}

	// востановление физических переменных
	for (auto& el : grid.cells)
	{
		hllc_get_phys_value(el.conv_val, el.phys_val);
	}

	return 0;

}
#endif //NEW_CLASS


#endif // HLLC3d

