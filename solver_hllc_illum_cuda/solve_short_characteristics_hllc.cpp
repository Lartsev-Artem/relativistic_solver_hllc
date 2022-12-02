#include "solve_short_characteristics_hllc.h"

#ifdef USE_VTK

#else
#define SIGN(a) (a < 0.0 ? -1.0 : 1.0) 

std::vector<int> b_right;
std::vector<int> b_left;

static int c0 = 0;
static int c1 = 0;
static int c2 = 0;
static int c3 = 0;
VectorX GetPhysVal( Type p0, const VectorX& U);

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
	
	T = Eigen::MatrixXd::Zero(5,5);
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
	Eigen::VectorXd tU;

	bool flag_bound = 1; // граничная ячейка или нет	

	for (size_t i = 0; i < 4; i++) // по граням
	{
		const int neig = neighbours_id_faces[4 * num_cell + i];
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
		ofstream ofile;
		const std::string name_file = "D:\\Desktop\\FilesCourse\\dbg_hllc.txt";
		static bool start = false;
		if (!start)
		{
			start = true;
			ofile.open(name_file);
			ofile.close();
		}

		if (fabs(F(0)) > 1e-1)
		{
			ofile.open(name_file, std::ios::app);
			if (!ofile.is_open())
			{
				printf("No open file\n");
				exit(1);
			}
			int i = num_cell;
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

#ifdef DBG_OUTPUT
	omp_set_num_threads(1);
#endif
#pragma omp parallel default(none) shared(N, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, U_full)
		{
			VectorX buf(5);
#pragma omp for
			for (int i = 0; i < N; i++)
			{
				buf = HLLC_stepToOMP(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev);

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

#ifdef DBG_OUTPUT
		printf("%d %d %d %d\n", c0, c1, c2, c3);
		printf("End output\n");
		static int count = 0;
		if (count++ > 2)
		{
			exit(1);
		}
#endif
		//printf("counts %d, %d, %d", c0, c1, c2);
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

#if 0 //defined Cone // условие отражение на боковой поверхности конуса
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






#ifdef RHLLC

static std::vector<VectorX> W_full;

static inline void MakeRotationMatrix(const Vector3& n, Matrix3& T)
{
	const Type x = n[0];
	const Type y = n[1];
	const Type theta = atan2(y, x);
	const Type phi = atan2(sqrt(x * x + y * y), n[2]);

	T(0, 0) = cos(theta) * sin(phi);
	T(0, 1) = sin(theta) * sin(phi);
	T(0, 2) = cos(phi);

	T(1, 0) = -sin(theta);
	T(1, 1) = cos(theta);
	T(1, 2) = 0;

	T(2, 0) = -cos(theta) * cos(phi);
	T(2, 1) = -sin(theta) * cos(phi);
	T(2, 2) = sin(phi);

	return;
}

int ReBuildDataForHLLCRel(const int N, std::vector<VectorX>& data) {

	data.resize(N);
	VectorX cell(5);

	W_full.resize(N);
	VectorX PhysCell(5);

	for (size_t i = 0; i < N; i++)
	{
		const Type v = velocity[i].dot(velocity[i]);
		const Type d = density[i];
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * pressure[i] / d;
		const Type dhGG = d * h * Gamma * Gamma;

		cell[0] = Gamma * d;
		cell[1] = dhGG * velocity[i][0];
		cell[2] = dhGG * velocity[i][1];
		cell[3] = dhGG * velocity[i][2];

		cell[4] = dhGG - pressure[i];
		//cell[4] = pressure[i] / (gamma1 - 1) + d * v / 2;
		//cell[4] = e_substance[i];

		data[i] = cell;

		PhysCell << density[i], velocity[i][0], velocity[i][1], velocity[i][2], pressure[i];
		W_full[i] = PhysCell;
	}

//	velocity.clear();
//	density.clear();
//	e_substance.clear();

//	printf("Array velocity, density was cleared\n");

	return 0;
}

VectorX NewStyleMatrix_HLLC_stepToOMPRel(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev) {

	//const Type h = 1 + (gamma1) / (gamma1 - 1) * p / rho;
	//const Type Gamma = 1. / sqrt(1 - v * v);

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	Eigen::Matrix3d T;
	Eigen::Vector3d tU;


	for (size_t i = 0; i < 4; i++) // по граням
	{
		const int neig = neighbours_id_faces[4 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;

#if defined Cone && defined Jet && !defined RHLLC // условие разрежения. См как заданы границы в make_struct.prj

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

				vel << U[1], U[2], U[3];
				tU = T * vel;
				tU[0] = -tU[0];
				vel = (T.transpose()) * tU;

				U_R[0] = U[0];
				U_R[1] = vel[0];
				U_R[2] = vel[1];
				U_R[3] = vel[2];
				U_R[4] = U[4];
				break;
			case eBound_LockBound:

				MakeRotationMatrix(normals[num_cell].n[i], T);

				vel << U[1], U[2], U[3];
				tU = T * vel;
				tU[0] = -tU[0];
				vel = (T.transpose()) * tU;

				U_R[0] = U[0];
				U_R[1] = vel[0];
				U_R[2] = vel[1];
				U_R[3] = vel[2];
				U_R[4] = U[4];

				break;
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				U_R = U_full_prev[neig / 4];
				break;
			}
		}

		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????
		Vector3 RotU(U[1], U[2], U[3]);
		RotU = T * RotU;
		U_L[0] = U[0];
		U_L[1] = RotU[0];
		U_L[2] = RotU[1];
		U_L[3] = RotU[2];
		U_L[4] = U[4];
		//U_L = T * U;

		RotU << U_R[1], U_R[2], U_R[3];
		RotU = T * RotU;
		U_R[0] = U_R[0];
		U_R[1] = RotU[0];
		U_R[2] = RotU[1];
		U_R[3] = RotU[2];
		U_R[4] = U_R[4];
		//U_R = T * U_R;

		// нормальная сокорость
		const Vector3 Vel_L = T * velocity[num_cell];
		const Vector3 Vel_R = T * velocity[neig / 4];


		Type d_L = density[num_cell];//U_L[0]; 
		Type d_R = density[neig / 4]; //U_R[0]; 		

		//Vector3 vel(U_L(1) / d_L, U_L(2) / d_L, U_L(3) / d_L);
		Type v = Vel_L.dot(Vel_L); // vel.dot(vel);

		const Type h_L = 1 + (gamma1 / (gamma1 - 1)) * pressure[num_cell] / d_L;
		Type Gamma_L = 1. / sqrt(1 - v);
		Type W_L = d_L * h_L * Gamma_L * Gamma_L;
		/////////Type p_L = (W_L - U_L[0] * Gamma_L) / ((gamma1 / (gamma1 - 1)) * Gamma_L * Gamma_L);
		Type p_L = pressure[num_cell];

		//Type p_L = (U_L(4) - v * d_L / 2.) * (gamma1 - 1);

		//vel << U_R(1) / d_R, U_R(2) / d_R, U_R(3) / d_R;
		v = Vel_R.dot(Vel_R);//vel.dot(vel);
		const Type h_R = 1 + (gamma1 / (gamma1 - 1)) * pressure[neig / 4] / d_R;
		Type Gamma_R = 1. / sqrt(1 - v);
		Type W_R = d_R * h_R * Gamma_R * Gamma_R;
		////////Type p_R = (W_R - U_R[0] * Gamma_R) / ((gamma1 / (gamma1 - 1)) * Gamma_R * Gamma_R);
		Type p_R = pressure[neig / 4];

		//Type p_R = (U_R(4) - v * d_R / 2.) * (gamma1 - 1);


		//const Type p_L = pressure[num_cell]; 
		//const Type p_R = pressure[neig];

		const Type v_L = Vel_L[0];//U_L[1] / d_L; // sqrt(U_L[1] * U_L[1] + U_L[2] * U_L[2] + U_L[3] * U_L[3]);  //velocity[num_cell].norm();
		const Type v_R = Vel_R[0];//U_R[1] / d_R; //sqrt(U_R[1] * U_R[1] + U_R[2] * U_R[2] + U_R[3] * U_R[3]); //velocity[neig].norm();

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
		if (S_R <= 0) // если верно выполнить всегда
		{
			F(0) = U_R[0] * Vel_R[0]; //D*v_x
			F(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
			F(2) = U_R[2] * Vel_R[0];
			F(3) = U_R[3] * Vel_R[0];
			F(4) = U_R[1];
			//continue;
		}
		else if (S_L >= 0) // выполнить либо по условию либо для всех границ
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
			const Type roS_L = d_L * (S_L - v_L);
			const Type roS_R = d_R * (S_R - v_R);
			const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

			const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;

			Eigen::VectorXd D(5); D << 0, 1, 0, 0, _S;

			if (_S >= 0)
			{
				F(0) = U_L[0] * Vel_L[0]; //D*v_x
				F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
				F(2) = U_L[2] * Vel_L[0];
				F(3) = U_L[3] * Vel_L[0];
				F(4) = U_L[1];

				VectorX buf = F;
				F = (_S * (S_L * U_L - buf) + S_L * P_LR * D) / (S_L - _S);
			}
			else //(_S <= 0)
			{
				F(0) = U_R[0] * Vel_R[0]; //D*v_x
				F(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
				F(2) = U_R[2] * Vel_R[0];
				F(3) = U_R[3] * Vel_R[0];
				F(4) = U_R[1];

				VectorX buf = F;
				F = (_S * (S_R * U_R - buf) + S_R * P_LR * D) / (S_R - _S);
			}
		}

		Vector3 RotF(F[1], F[2], F[3]);
		RotF = (T.transpose()) * RotF;
		//F[0] = F[0];
		F[1] = RotF[0];
		F[2] = RotF[1];
		F[3] = RotF[2];
		//F[4] = F[4];
		//F = (T.transpose()) * buf;

		SumF += F * squares_cell[4 * num_cell + i];

	}// for


	return (U - SumF * tau / volume[num_cell]);
}


VectorX HLLC_stepToOMPRel(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev) {

	/*
	An HLLC Riemann solver for relativistic flows – I. Hydrodynamics 
	A. Mignone and G. Bodo
	INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
	Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
	Mon. Not. R. Astron. Soc. 364, 126–136 (2005)
	*/

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5,5);
	Eigen::VectorXd tU(5);
	
	for (size_t i = 0; i < 4; i++) // по граням
	{
		int neig = neighbours_id_faces[4 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;
				neig = num_cell * 4 + i;
				break;
			case eBound_InnerSource:
				U_R = U;
				neig = num_cell * 4 + i;
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

				U_R = U_full_prev[neig / 4];
				break;
			}
		}
			
		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????

		U_L = T * U;
		U_R = T * U_R;
//==================== Кэшируем физические переменные слева и справа============================//
		// нормальная сокорость
		Eigen::VectorXd bufVL(5), bufVR(5);
		bufVL << 0, velocity[num_cell][0], velocity[num_cell][1], velocity[num_cell][2], 0;
		bufVL = T * bufVL;
		const Vector3 Vel_L(bufVL[1], bufVL[2], bufVL[3]);  //T * velocity[num_cell];

		bufVR << 0, velocity[neig / 4][0], velocity[neig / 4][1], velocity[neig / 4][2], 0;
		bufVR = T * bufVR;
		const Vector3 Vel_R(bufVR[1], bufVR[2], bufVR[3]);  //T * velocity[neig / 4];


		const Type d_L = density[num_cell];
		const Type d_R = density[neig/4]; 

		const Type p_L = pressure[num_cell];  
		const Type p_R = pressure[neig / 4]; 

//========================================================================================//

//=========================Вычисляем релятивистикие параметры============================//
		const Type g_L = 1. / sqrt(1 - (Vel_L.dot(Vel_L)));	// фактор Лоренца
		const Type g_R = 1. / sqrt(1 - (Vel_R.dot(Vel_R)));

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
			c0++; // для hllc этот код срабатывает на 8 принте для задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[0] * Vel_L[0]; //D*v_x
			F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F(2) = U_L[2] * Vel_L[0];
			F(3) = U_L[3] * Vel_L[0];
			F(4) = U_L[1];
			//continue;
			c1++;
		}
		else
		{
//====================Расчёт потоков и приближений hll=========================================//
			Eigen::VectorXd F_L(5);
			Eigen::VectorXd F_R(5);
			Eigen::VectorXd U_hll(5);
			Eigen::VectorXd F_hll(5);

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

			F_hll = (lambda_R * F_L - lambda_L * F_R + lambda_R * lambda_L * (U_R - U_L)) / (lambda_R - lambda_L);
			U_hll = (lambda_R * U_R - lambda_L * U_L + (F_L - F_R)) / (lambda_R - lambda_L);

		//	cout << "F_hll\n" << F_hll << "\nU_hll\n" << U_hll << '\n';

		//F = F_hll;
			
//============================================================================================//
#if 1		
//=========================Поиск скорости промежуточной волны===============================//
			const Type a = F_hll[4];			//F_E^hll
			const Type b = U_hll[4] + F_hll[1]; // (E_hll + F_mx^hll)
			const Type c = U_hll[1];			//mx_hll

#if 1 // как описано в Mignone...
			bool flag_hll = false;

			Type _lambda = 0;
			if (fabs(a) > 1e-10 /*&& fabs(c) > eps*/)
			{
				Type D = b * b - 4 * a * c;

				if (D < 0)
				{
					printf("Err D %d (a=%lf, b=%lf, c=%lf)\n", num_cell, a, b, c);
					exit(1);
				}

				_lambda = (b - sqrt(D)) / (2 * a); // здесь b без минуса в силу уравнения	
			}
			else
			{
				//VectorX physVal = GetPhysVal(p_L, U_hll);
				//_lambda = c / physVal[0];
				//flag_hll = true;

				_lambda = c / b;				
			}

/*			if (flag_hll) 
			{
				F = F_hll;
			}
			else*/			
#endif		
			{
				if (_lambda >= 0)
				{
					c2++;
					//============================Поиск промежуточного давления ===================================//
					const Type A = lambda_L * U_L[4] - U_L[1];
					const Type B = U_L[1] * (lambda_L - Vel_L[0]) - p_L;
					const Type _p = (-B + A * _lambda) / (1 + _lambda * lambda_L);
					//============================================================================================//

					//==========================Финальный поток HLLC=============================================//
					Eigen::VectorXd _U_L(5);
					const Type dif_L = lambda_L - _lambda;

					_U_L = U_L * ((lambda_L - Vel_L[0]));

					_U_L[1] += (_p - p_L);
					_U_L[4] += (_p * _lambda - p_L * Vel_L[0]);

					_U_L /= (lambda_L - _lambda);


					//_U_L[0] = (U_L[0] * (lambda_L - Vel_L[0])) / dif_L;
					//_U_L[1] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) / dif_L;
					//_U_L[2] = (U_L[2] * (lambda_L - Vel_L[0]) ) / dif_L;
					//_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0]) ) / dif_L;
					//_U_L[4] = (U_L[4] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) / dif_L;

					F = F_L + lambda_L * (_U_L - U_L);

					//============================================================================================//
				}
				else //(_S <= 0)
				{
					c3++;
					//============================Поиск промежуточного давления ===================================//
					const Type A = lambda_R * U_R[4] - U_R[1];
					const Type B = U_R[1] * (lambda_R - Vel_R[0]) - p_R;
					const Type _p = (-B + A * _lambda) / (1 + _lambda * lambda_R);
					//============================================================================================//
					Eigen::VectorXd _U_R(5);
					const Type dif_R = lambda_R - _lambda;

					_U_R = U_R * ((lambda_R - Vel_R[0]));

					_U_R[1] += (_p - p_R);
					_U_R[4] += (_p * _lambda - p_R * Vel_R[0]);

					_U_R /= (lambda_R - _lambda);

					//_U_R[0] = (U_R[0] * (lambda_R - Vel_R[0])) / dif_R;
					//_U_R[1] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) / dif_R;
					//_U_R[2] = (U_R[2] * (lambda_R - Vel_R[0])) / dif_R;
					//_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0])) / dif_R;
					//_U_R[4] = (U_R[4] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) / dif_R;

					F = F_R + lambda_R * (_U_R - U_R);
				}
			}
#endif
		}


		VectorX buf = F;
		F = (T.transpose()) * buf;
		SumF += F * squares_cell[4 * num_cell + i];	
	//	cout << "flew "<<i <<"\n" << F << '\n';

	}// for

	//cout << "ready flew\n" << SumF << '\n';
	return (U - SumF * tau / volume[num_cell]);
}


VectorX RHLLC_stepToOMPGit(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev) {

	/*
	An HLLC Riemann solver for relativistic flows – I. Hydrodynamics
	A. Mignone and G. Bodo
	INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
	Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
	Mon. Not. R. Astron. Soc. 364, 126–136 (2005)

	https://github.com/PrincetonUniversity/Athena-Cversion/blob/master/src/rsolvers/hllc_sr.c

	\note: комбинация кода из mignone 2005 и 2006. в части hllc
	*/

	Eigen::VectorXd SumF = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
#ifdef DBG_OUTPUT
	Eigen::VectorXd SumFdbg = Eigen::VectorXd::Zero(5);  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
#endif

	const Eigen::VectorXd U = U_full_prev[num_cell];
	Eigen::VectorXd F = Eigen::VectorXd::Zero(5);

	Eigen::VectorXd U_R(5);
	Eigen::VectorXd U_L(5);

	const Eigen::VectorXd W = W_full[num_cell];

	Eigen::VectorXd W_R(5);
	Eigen::VectorXd W_L(5);

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(5, 5);
	Eigen::VectorXd tU(5);

	for (size_t i = 0; i < 4; i++) // по граням
	{
		int neig = neighbours_id_faces[4 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;
				W_R = W;
				neig = num_cell * 4 + i;
				break;
			case eBound_OutSource://eBound_InnerSource:
				U_R = U;
				W_R = W;
				neig = num_cell * 4 + i;
				break;
			case eBound_InnerSource://eBound_OutSource:
				MakeRotationMatrix(normals[num_cell].n[i], T);
				tU = T * U;
				tU[1] = -tU[1];
				U_R = (T.transpose()) * tU;

				tU = T * W;
				tU[1] = -tU[1];
				W_R = (T.transpose()) * tU;				
				break;
			case eBound_LockBound:

				MakeRotationMatrix(normals[num_cell].n[i], T);
				tU = T * U;
				tU[1] = -tU[1];
				U_R = (T.transpose()) * tU;


				tU = T * W;
				tU[1] = -tU[1];
				W_R = (T.transpose()) * tU;
				break;
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				U_R = U_full_prev[neig / 4];
				W_R = W_full[neig / 4];
				break;
			}
		}

		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????

		U_L = T * U;
		U_R = T * U_R;

		W_L = T * W;
		W_R = T * W_R;

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
			c0++; // для hllc этот код срабатывает на 8 принте для задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[0] * Vel_L[0]; //D*v_x
			F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F(2) = U_L[2] * Vel_L[0];
			F(3) = U_L[3] * Vel_L[0];
			F(4) = U_L[1];
			//continue;
			c1++;
		}
		else
		{
			//====================Расчёт потоков и приближений hll=========================================//
			Eigen::VectorXd F_L(5);
			Eigen::VectorXd F_R(5);
			Eigen::VectorXd U_hll(5);
			Eigen::VectorXd F_hll(5);

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
			{
				if (_lambda >= 0.0)
				{
					c2++;
					//============================Поиск промежуточного давления ===================================//
					const Type _p = -F_hll[4] * _lambda + F_hll[1];
					//============================================================================================//

					//==========================Финальный поток HLLC=============================================//
					Eigen::VectorXd _U_L(5);
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
					c3++;
					//============================Поиск промежуточного давления ===================================//
					const Type _p = -F_hll[4] * _lambda + F_hll[1];
					//============================================================================================//
					Eigen::VectorXd _U_R(5);
					const Type dif_R = 1.0/(lambda_R - _lambda);

					_U_R[0] = (U_R[0] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[1] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
					_U_R[2] = (U_R[2] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[4] = (U_R[4] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

					F = F_R + lambda_R * (_U_R - U_R);
				}
			}
#endif
		}


		VectorX buf = F;
		F = (T.transpose()) * buf;
		SumF += F * squares_cell[4 * num_cell + i];
	
#ifdef DBG_OUTPUT
		SumFdbg += F;
#endif

	}// for

#ifdef DBG_OUTPUT
	ofstream ofile;
	const std::string name_file = "D:\\Desktop\\FilesCourse\\dbg_rhllc_flux.txt";
	static bool flag_init = true;
	if (flag_init)
	{
		ofile.open(name_file);
		ofile.close();
		flag_init = false;
	}
	ofile.open(name_file, std::ios::app);
	if (!ofile.is_open())
	{
		printf("No open file\n");
		exit(1);
	}
	int i = num_cell;
	ofile << "F[" << i << "]= " << SumFdbg(0) << ", "
		<< SumFdbg(1) << ", " << SumFdbg(2) << ", " << SumFdbg(3) << ", " << SumFdbg(4) << '\n';
	ofile.close();
#endif
	//cout << "ready flew\n" << SumF << '\n';
	return (U - SumF * tau / volume[num_cell]);
}


int ReBuildPhysicValue_ost1098(const int num_cell, const VectorX& U) {

	const Type vv = velocity[num_cell].dot(velocity[num_cell]);
	const Type d = density[num_cell];
	Type Gamma0 = 1. / sqrt(1 - vv);
	const Type h = 1 + gamma_g * pressure[num_cell] / d;

	Type W0 = d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m(U[1], U[2], U[3]);
	Type mm = m.dot(m);

	Type p = pressure[num_cell];
	Vector3 v = velocity[num_cell];

	Type D = U[0];
	Type E = U[4];

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0* Gamma0* Gamma0)*mm/(2*W0*W0*W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
		W0 -= (fW / dFdW);

		Gamma0 = 1./sqrt(1 - mm / (W0 * W0));
		
		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v[0] = m[0] / W0;
		v[1] = m[1] / W0;
		v[2] = m[2] / W0;

		err -= W0;
		cc++;
	} while (fabs(err/W0) > 1e-14);

	if (p < 0 || U[0] < 0 || std::isnan(p) || std::isnan(U[0]))
	{
		printf("Error cell %d (p = %lf, d= %lf)", num_cell, p, D / Gamma0);
		exit(1);
	}

#ifdef DBG_OUTPUT
	ofstream ofile;
	const std::string name_file = "D:\\Desktop\\FilesCourse\\dbg_rhllc2.txt";
	if (num_cell == 0)
	{
		ofile.open(name_file);
		ofile.close();
	}
	
	if (cc > 2)
	{		
		ofile.open(name_file, std::ios::app);
		if (!ofile.is_open())
		{
			printf("No open file\n");
			exit(1);
		}
		int i = num_cell;
		ofile << "U[" << i << "]= " << U(0) << ", " << U(1) << ", " << U(4) << '\n';
		ofile << "W[" << i << "]= " << D / Gamma0 << ", " << v(0) << ", " << p << '\n';
		
		ofile.close();		
	}
#endif
	pressure[num_cell] = p;
	velocity[num_cell] = v;
	//Gamma0;
	density[num_cell] = D / Gamma0;

	return 0;
}

int ReBuildPhysicValue(const int num_cell, const VectorX& U)
{
	const Type E = U[4];
	const Type D = U[0];
	const Vector3 M(U[1], U[2], U[3]);

	const Type mm = M.dot(M);
	Type p1 = pressure[num_cell];
	Type p0 = p1;

	do
	{
		p0 = p1;

		Type gp = 1. / sqrt(1 - (mm) / ((E + p0) * (E + p0)));
		Type Fp = E + p0 - D * gp - gamma_g * p0 * gp * gp;

		Type div = mm - (E + p0) * (E + p0);
		Type dFdP = 1 - gamma_g * gp * gp + (mm * (E + p0) * (2 * gamma_g * p0 + D * (1. / gp))) / (div * div);

		p1 = p0 - Fp / dFdP;
	} while (fabs((p0 - p1) / p1) > 1e-13);

	Type Dhg = E + p1;
	Vector3 V = M / Dhg;

	Type g = 1. / sqrt(1 - V.dot(V));

	Type rho = D / g;

	if (p1 < 0 || rho < 0 || std::isnan(p1) || std::isnan(rho))
	{
		printf("Error cell  (p = %lf, d= %lf)", p1, rho);
		exit(1);
	}

	density[num_cell] = rho;
	velocity[num_cell] = V;
	pressure[num_cell] = p1;
	
	return 0;
}
VectorX GetPhysVal(Type p0, const VectorX& U)
{
	const Type E = U[4];
	const Type D = U[0];
	const Vector3 M(U[1], U[2], U[3]);

	const Type mm = M.dot(M);	
	Type p1 = p0;

	do
	{
		p0 = p1;

		Type gp = 1. / sqrt(1 - (mm) / ((E + p0) * (E + p0)));
		Type Fp = E + p0 - D * gp - gamma_g * p0 * gp * gp;

		Type div = mm - (E + p0) * (E + p0);
		Type dFdP = 1 - gamma_g * gp * gp + (mm * (E + p0) * (2 * gamma_g * p0 + D * (1. / gp))) / (div * div);

		p1 = p0 - Fp / dFdP;
	} while (fabs((p0 - p1) / p1) > 1e-5);
	
	Type Dhg = E + p1;
	Vector3 V = M / Dhg;

	Type g = 1. / sqrt(1 - V.dot(V));

	Type rho = D / g;

	if (p1 < 0 || rho < 0 || std::isnan(p1) || std::isnan(rho))
	{
		printf("Error cell  (p = %lf, d= %lf)", p1, rho);
		exit(1);
	}

	VectorX res(5);
	res << rho, V[0], V[1], V[2], p1;
		
	return res;
}

void HLLC_Rel(const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, std::vector<VectorX>& U_full_prev)
{	
#ifdef DBG_OUTPUT
	omp_set_num_threads(1);
#endif
	ReBuildDataForHLLCRel(size_grid, U_full_prev); // p0, v0, rho0 -> U_full_prev

#pragma omp parallel default(none) shared(pressure, velocity, density, size_grid, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, U_full)
	{
		VectorX buf(5);
#pragma omp for
		for (int i = 0; i < size_grid; i++)
		{
			//buf = HLLC_stepToOMPRel(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev); // -> U_full		
			buf = RHLLC_stepToOMPGit(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev); // -> U_full		
#ifdef Cube
			// Для 1d задачи Сода. 
			// СЛИШКОМ СИЛЬНЫЕ КОЛЕБАНИЯ ПО Y,Z. что-то не так!!!!!!!!!!!!!!!!!!!!
			//buf[1] = Vector3(buf[1], buf[2], buf[3]).norm();
			buf[2] = 0;
			buf[3] = 0; 
#endif

#pragma omp critical
			{
 				U_full[i] = buf;		
			}
		}
	}

	for (size_t i = 0; i < size_grid; i++)
	{		
		ReBuildPhysicValue_ost1098(i, U_full[i]); // U_full -> p1, v1, rho1
		//ReBuildPhysicValue(i, U_full[i]); // U_full -> p1, v1, rho1
	}

#ifdef DBG_OUTPUT
	ofstream ofile;
	ofile.open("D:\\Desktop\\FilesCourse\\dbg_rhllc.txt");
	if (!ofile.is_open())
	{
		printf("No open file\n");
		exit(1);
	}
	for (size_t i = 0; i < size_grid; i++)
	{
		if (fabs(velocity[i](1)) > 1e-8 || fabs(velocity[i](2)) > 1e-8)
		{
			printf("Error velocity for 1d soda\n");
			exit(2);
		}
		ofile << "U= " << U_full[i](0) << ", " << U_full[i](1) << ", " << U_full[i](4) << '\n';
		ofile << "W= " << density[i] << ", " << velocity[i](0) << ", " << pressure[i] << '\n';
	}
	ofile.close();
	printf("dbg File writed\n");
	exit(1);
#endif
	
	//printf("couts: %d, %d, %d, %d\n", c0, c1, c2,c3);
	U_full.swap(U_full_prev);
}


Type FormTimeStepToRHLLC(const int n, const Type h, const Type k) {

	// формирование шага по времени		
	Type c_max = -1;
	for (size_t i = 0; i < n; i++)
	{		
		const Type d = density[i]; 
		const Vector3 vel = velocity[i];
		const Type v = vel.dot(vel);
		const Type p = pressure[i]; 

		const Type gg = 1. / (sqrt(1 - v));
		const Type ent = 1 + gamma_g * p / d;
		const Type cs = sqrt(gamma1 * p / (d * ent));
		const Type sigma = cs * cs / (gg * gg * (1 - cs * cs));
		
		const Type c = (sqrt(v) + sqrt(sigma * (1 - v + sigma))) / (1 + sigma);
		if (c > c_max) c_max = c;
	}
	//printf("c=%lf\n", c_max);
	return k * h * c_max;
}
#endif //RHLLC

#endif // USE_VTK