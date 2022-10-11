#include "solve_short_characteristics_hllc.h"

#ifdef USE_VTK

#else
int count1=0, count2=0, count3=0, count4=0;

std::vector<int> b_right;
std::vector<int> b_left;

int ReBuildDataForHLLC(const int N, std::vector<VectorX>& data) {

	data.resize(N);
	VectorX cell(5);
	
	for (size_t i = 0; i < N; i++)
	{
		const Type d = density[i];
		cell[0] = d;
		cell[1] = d * velocity[i][0];
		cell[2] = d * velocity[i][1];
		cell[3] = d * velocity[i][2];

		const Type v = velocity[i].dot(velocity[i]);
		cell[4] = pressure[i] / (gamma1 - 1) + d * v / 2;
		//cell[4] = e_substance[i];

		data[i] = cell;
	}

	velocity.clear();
	density.clear();
	e_substance.clear();

	printf("Array velocity, density was cleared\n");
	
	return 0;
}

inline void MakeRotationMatrix(const Vector3& n, Eigen::MatrixXd& T) {
	
	//T = Eigen::MatrixXd::Zero(5,5);
	T(0, 0) = T(4, 4) = 1;

	if (fabs(n[2] * n[2] - 1) > eps)
	{

		T(1, 1) = n[0];
		T(1, 2) = n[1];
		T(1, 3) = n[2];

		Type sqr = sqrt(1 - n[2] * n[2]);

		if (sqr < eps * eps - eps / 10)
			printf("Err T\n");

		T(2, 1) = -n[1] / sqr;
		T(2, 2) = n[0] / sqr;

		T(3, 1) = -n[0] * n[2] / sqr;
		T(3, 2) = -n[1] * n[2] / sqr;
		T(3, 3) = sqr;
	}
	else if (n[2] > 0)  // n_z == 1
	{
		T(1, 3) = 1;
		T(2, 2) = 1;
		T(3, 1) = -1;
	}
	else  // n_z == -1
	{
		T(1, 3) = -1;
		T(2, 2) = -1;
		T(3, 1) = 1;
	}
	
}


VectorX HLLC_stepToOMP(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev) {

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);

	bool flag_bound = 0; // граничная ячейка или нет	

	//if (num_cell == 10766)
	//{
	//	//for (size_t k = 0; k < 5; k++)
	//	//{
	//	//	printf("%.10lf\n", U[k]);
	//	//}
	//	//printf("\n\n");
	//}

	for (size_t i = 0; i < 4; i++) // по граням
	{
		const int neig = neighbours_id_faces[4 * num_cell + i];

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

		//U_R << d, d* velocity[neig][0], d* velocity[neig][1], d* velocity[neig][2], e_substance[neig];		

		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????
		U_L = T * U;
		U_R = T * U_R;

		// if(nan -> save_dump_prev_state and restart to debug)
		if (flag_bound) 
		{

#if 0// defined Sphere
			if (neighbours_id_faces[4 * num_cell + i] == -1) // внешняя  граница
			{
				if (normals[num_cell].n[i].dot(Vector3(1, 0, 0)) > 0)
				{
					U_R[1] *= -1; // отражение на полу сфере
				}
			}			
#endif

#if defined Cone
			if (fabs((normals[num_cell].n[i] - Vector3(1, 0, 0)).norm()) > eps)  // если не открытый конец, то отражение				
			{
				//U_R[1] *= -1;
			}
			else // открытый конец -> условие разрежениия
			{
#if 1
				Type d = U_R[0];
				Vector3 vel(U_R[1] / d, U_R[2] / d, U_R[3] / d);
				Type v = vel.norm();
				double pressure = (U_R[4] - v * v * d / 2.) * (gamma1 - 1);

				d *= 0.5;
				vel[0] *= 1.1;
				
				v = vel.dot(vel);
				U_R[0] = d;
				U_R[1] = d * vel[0];
				U_R[2] = d * vel[1];
				U_R[3] = d * vel[2];				
				U_R[4] = pressure / (gamma1 - 1) + d * v / 2;
#endif
			}
#endif
		}

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

		//if (_Is_nan(p_L)) 
		//{
		//	printf("NaN\n"); exit(1);
		//}
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
		flag_bound = false;
		if (S_R <= 0) // если верно выполнить всегда
		{
			F(0) = U_R[1];
			F(1) = U_R[1] * U_R[1] / d_R + p_R;
			F(2) = U_R[1] * U_R[2] / d_R;
			F(3) = U_R[1] * U_R[3] / d_R;
			F(4) = (U_R[4] + p_R) * U_R[1] / d_R;
			//continue;
			if(num_cell==2000)
				count1++;
		}
		else if (S_L >= 0 || flag_bound) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[1];
			F(1) = U_L[1] * U_L[1] / d_L + p_L;
			F(2) = U_L[1] * U_L[2] / d_L;
			F(3) = U_L[1] * U_L[3] / d_L;
			F(4) = (U_L[4] + p_L) * U_L[1] / d_L;  // TU[4]*d_L
			//continue;
			if (num_cell == 2000)
				count2++;
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
				if (num_cell == 2000)
					count3++;

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
				if (num_cell == 2000)
					count4++;

			}
		}

		VectorX buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[4 * num_cell + i];

	}// for

	
	/*U_full[num_cell] =*/ return (U - SumF * tau / volume[num_cell]);
}


int HLLC_step(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	std::vector<Eigen::VectorXd>& U_full, const std::vector<VectorX>& U_full_prev) {

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);

	bool flag_bound = 0; // граничная ячейка или нет

	for (size_t i = 0; i < 4; i++) // по граням
	{
		const int neig = neighbours_id_faces[4 * num_cell + i];

		if (neig >= 0) 
		{
			flag_bound = 0;
			U_R = U_full_prev[neig / 4];
		}
		else 
		{
			flag_bound = 1;
			U_R = U;			
		}

		//U_R << d, d* velocity[neig][0], d* velocity[neig][1], d* velocity[neig][2], e_substance[neig];		
		
		MakeRotationMatrix(normals[num_cell].n[i], T);
		U_L = T * U;
		U_R = T * U_R;
		

		 Type d_L = U_L[0]; //density[num_cell];
		 Type d_R = U_R[0]; // density[neig];

		Vector3 vel(U_L(1) / d_L, U_L(2) / d_L, U_L(3) / d_L);
		Type v = vel.dot(vel);
		 Type p_L = (U_L(4) - v * d_L / 2.) * (gamma1 - 1);

		vel << U_R(1) / d_R, U_R(2) / d_R, U_R(3) / d_R;
		v = vel.dot(vel);
		Type p_R = (U_R(4) - v * d_R / 2.) * (gamma1 - 1);

		//if (_Is_nan(p_L)) 
		//{
		//	printf("NaN\n"); exit(1);
		//}
		
		//if (p_L < 0) //printf("(%d) p_L= %f\n", num_cell, p_L);

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

		if (S_R <= 0) // если верно выполнить всегда
		{
			F(0) = U_R[1];
			F(1) = U_R[1] * U_R[1] / d_R + p_R;
			F(2) = U_R[1] * U_R[2] / d_R;
			F(3) = U_R[1] * U_R[3] / d_R;
			F(4) = (U_R[4] + p_R) * U_R[1] / d_R;
			//continue;
		}
		else if (S_L >= 0 || flag_bound) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[1];
			F(1) = U_L[1] * U_L[1] / d_L + p_L;
			F(2) = U_L[1] * U_L[2] / d_L;
			F(3) = U_L[1] * U_L[3] / d_L;
			F(4) = (U_L[4] + p_L) * U_L[1] / d_L;  // TU[4]*d_L
			//continue;
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
			}
		}

		VectorX buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[4 * num_cell + i];

	}// for

	U_full[num_cell] = U - SumF * tau / volume[num_cell];
// скорости:	
 
	// Для задачи SODA
 U_full[num_cell][2] = 0; 
 U_full[num_cell][3] = 0;

	return 0;
}



int HLLC(const int N, const Type tau, const std::vector<int>& bound_cells, const std::vector<int>& neighbours_id_faces,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	const std::vector<VectorX>& U_full_prev) {

#pragma omp parallel default(none) shared(N, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, U_full)
		{
			VectorX buf(5);
#pragma omp for
			for (size_t i = 0; i < N; i++)
			{
				buf = HLLC_stepToOMP(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);
#pragma omp critical
				{
					U_full[i] = buf;
					// Для задачи SODA:
						//U_full[i][2] = 0;
						//U_full[i][3] = 0;
				}
			}
		}
#if 0 //defined Sphere // условие отражение на полу-сфере
		for (size_t num_cell = 0; num_cell < N; num_cell++)
		{
			for (size_t i = 0; i < 4; i++)			
			if (neighbours_id_faces[4 * num_cell + i] == -1) // внешняя  граница
			{
				if (normals[num_cell].n[i].dot(Vector3(1, 0, 0)) > 0)
				{
					Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
					MakeRotationMatrix(normals[num_cell].n[i], T);
					U_full[num_cell] = T * U_full[num_cell];
					U_full[num_cell][1] *= -1;
					U_full[num_cell] = (T.transpose()) * U_full[num_cell];
				}
			}
		}
#endif

#if defined Cone // условие отражение на боковой поверхности конуса
		for (size_t num_cell = 0; num_cell < N; num_cell++)
		{
			for (size_t i = 0; i < 4; i++)
				if (neighbours_id_faces[4 * num_cell + i] == -1) // внешняя  граница
				{
					if ((normals[num_cell].n[i] - (Vector3(1, 0, 0))).norm() > eps)
					{
						Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
						MakeRotationMatrix(normals[num_cell].n[i], T);
						U_full[num_cell] = T * U_full[num_cell];
						U_full[num_cell][1] *= -1;
						U_full[num_cell] = (T.transpose()) * U_full[num_cell];
					}
				}
		}
#endif
		// расчет ГУ путем усреднения по соседним к граничным ячейкам 
#if 0	

		int neighb = 0; // номер внешней грани
		register int count = 0;  // число граничных граней у ячейки
		register int cc = 0;		// число соседних граней  {4-count} 

		for (auto num_cell : bound_cells)
		{
			VectorX U = VectorX::Zero(5);   
			cc = 0;
			count = 0;
			
			for (size_t j = 0; j < 4; j++)
			{
				int n = neighbours_id_faces[4 * num_cell + j];
				if (n >= 0)
				{
					U += U_full[n / 4];
					cc++;
				}
				else 
				{
					neighb = j; 
					count++; 
				}
			}
			U_full[num_cell] = U / cc;

			if (count == 1) // если это боковая поверхность(не угол), то условие отражения
			{
				Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
				MakeRotationMatrix(normals[num_cell].n[neighb], T);
				U_full[num_cell] = T * U_full[num_cell];
				U_full[num_cell][1] *= -1;
				U_full[num_cell] = (T.transpose()) * U_full[num_cell];
				//U_full[num_cell][2] = 0;
				//U_full[num_cell][3] = 0;
			}			
		}

#endif
		
		return 0;
}
#endif // USE_VTK