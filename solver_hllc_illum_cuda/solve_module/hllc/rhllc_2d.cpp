#include "solve_short_characteristics_hllc.h"




#define SIGN(a) (a < 0.0 ? -1.0 : 1.0) 

#ifdef USE_VTK

#else

//#define SODA_2d

#define Jet_2d

static int c0 = 0;
static int c1 = 0;
static int c2 = 0;
static int c3 = 0;

typedef Eigen::Vector4d Vector4;
typedef Eigen::Vector2d Vector2;
typedef Eigen::Matrix4d Matrix4;

//#ifndef USE_MPI
static std::vector<Vector4> W_full_2d;
static std::vector<Vector4> U_full_2d;
static std::vector<Vector4> U_full_2d_prev;
//#endif


static int WriteFileSolution(const std::string& main_dir, const std::vector<Type>& density,
	const std::vector<Type>& pressure, const std::vector<Vector3>& velocity)
{
	FILE* f;

	f = fopen((main_dir + "density.bin").c_str(), "wb");
	int n = density.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(density.data(), sizeof(Type), n, f);
	fclose(f);

	f = fopen((main_dir + "pressure.bin").c_str(), "wb");
	n = pressure.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(pressure.data(), sizeof(Type), n, f);
	fclose(f);

	f = fopen((main_dir + "velocity.bin").c_str(), "wb");
	n = velocity.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(velocity.data(), sizeof(Vector3), n, f);
	fclose(f);

	return 0;
}

static  int WriteSolution(const int n, const std::string main_dir, const std::vector<Vector4>& W)
{
	std::vector<Type> den(size_grid);
	std::vector<Type> press(size_grid);
	std::vector<Vector3> vel(size_grid);

	for (size_t i = 0; i < size_grid; i++)
	{		
		den[i] = W[i](0);

		vel[i](0) = W[i](1);
		vel[i](1) = W[i](2);
		vel[i](2) = 0;
		
		press[i] = W[i](3);
	}

	WriteFileSolution(main_dir + "Solve" + std::to_string(n), den, press, vel);

	return 0;
}

static int InitNeigh(const int N, const std::vector<Vector3>& centerts, const std::vector<Normals>& normals,
	std::vector<int>& neighbours_id_faces)
{
#ifdef Jet_2d
	for (size_t i = 0; i < N; i++)
	{
		Vector2 x(centerts[i][0], centerts[i][1]);

		for (size_t j = 0; j < 3; j++)
		{
			if (neighbours_id_faces[3 * i + j] < 0)
				if (fabs(normals[i].n[j][0] - 1) < 0.0001)
					if (x[1] < 1 && x[1]>-1 && x[0] < 1)
				{
					neighbours_id_faces[3 * i + j] = eBound_OutSource;
				}
		}
	}
#endif

	return 0;
}
static int RHLLC_Init(const int N, const std::vector<Vector3>& centerts, std::vector<Vector4>& W) {

	W.resize(N);
	Vector4 cell;

	for (size_t i = 0; i < N; i++)
	{
		Vector2 x(centerts[i][0], centerts[i][1]);

#ifdef Jet_2d
		if (x[0] < 1 && x[1] < 1 && x[1]>-1)
		{
			cell(0) = 0.1;			
			cell(1) = 0.99;
			cell(2) = 0;
			cell(3) = 0.01;			
		}
		else
		{
			cell(0) = 10;			
			cell(1) = 0;
			cell(2) = 0;
			cell(3) = 0.01;
		}
#endif

#ifdef SODA_2d
		if (x[0] < 0.499)
		{
			cell(0) = 1;
			cell(1) = 0.9;
			cell(2) = 0;
			cell(3) = 1;
		}
		else
		{
			cell(0) = 1;
			cell(1) = 0;
			cell(2) = 0;
			cell(3) = 10;
		}
#endif

		W[i] = cell;
	}

	return 0;
}

static inline void MakeRotationMatrix(const Vector3& n, Matrix4& T) {

	// n=(x,y,0)!!!!

	T = Matrix4::Zero();
	T(0, 0) = T(3, 3) = 1;

	//T(1, 1) = n[0];
	//T(1, 2) = n[1];

	//T(2, 1) = -n[1];
	//T(2, 2) = n[0];

	T(1, 1) = -n[0];
	T(1, 2) = -n[1];

	T(2, 1) = n[1];
	T(2, 2) = -n[0];

}

Vector4 RHLLC_stepToOMP2d(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume) {

	Vector4 SumF = Vector4::Zero();  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Vector4 U = U_full_2d_prev[num_cell];
	Vector4 F = Vector4::Zero();

	Vector4 U_R;
	Vector4 U_L;

	Vector4 W = W_full_2d[num_cell];
	Vector4 W_R;
	Vector4 W_L;

	Matrix4 T = Matrix4::Zero();

	for (size_t i = 0; i < 3; i++) // по граням
	{
		const int neig = neighbours_id_faces[3 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;
				W_R = W;
				break;
			case eBound_OutSource://eBound_InnerSource:
#ifdef Jet_2d				
				W_R << 0.1, 0.99, 0, 0.01;

				U_R << 0.7088812050083355,
					6.218592964824112,					
					0,
					6.271407035175871;
#else
				U_R = U;
				W_R = W;
#endif
				break;
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				U_R = U_full_2d_prev[neig / 3];
				W_R = W_full_2d[neig / 3];

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
		const Vector2 Vel_L(W_L[1], W_L[2]);// , 0);  //T * velocity[num_cell];
		const Vector2 Vel_R(W_R[1], W_R[2]);// , 0);  //T * velocity[neig / 4];

		const Type d_L = W_L(0);
		const Type d_R = W_R(0);

		const Type p_L = W_L(3);
		const Type p_R = W_R(3);

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
			//F(3) = U_R[3] * Vel_R[0];
			F(3) = U_R[1];
			//continue;
			c0++; // для hllc этот код срабатывает на 8 принте для задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[0] * Vel_L[0]; //D*v_x
			F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F(2) = U_L[2] * Vel_L[0];
			//F(3) = U_L[3] * Vel_L[0];
			F(3) = U_L[1];
			//continue;
			c1++;
		}
		else
		{
			//====================Расчёт потоков и приближений hll=========================================//
			Vector4 F_L;
			Vector4 F_R;
			Vector4 U_hll;
			Vector4 F_hll;

			F_R(0) = U_R[0] * Vel_R[0]; //D*v_x
			F_R(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
			F_R(2) = U_R[2] * Vel_R[0];
			//F_R(3) = U_R[3] * Vel_R[0];
			F_R(3) = U_R[1];

			F_L(0) = U_L[0] * Vel_L[0]; //D*v_x
			F_L(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F_L(2) = U_L[2] * Vel_L[0];
			//F_L(3) = U_L[3] * Vel_L[0];
			F_L(3) = U_L[1];

			F_hll = (lambda_R * F_L - lambda_L * F_R + (lambda_R * lambda_L * (U_R - U_L))) / (lambda_R - lambda_L);
			U_hll = (lambda_R * U_R - lambda_L * U_L + (F_L - F_R)) / (lambda_R - lambda_L);
			
#ifdef ONLY_RHLL
			F = F_hll;
#endif

			//============================================================================================//
#ifndef ONLY_RHLL		
//=========================Поиск скорости промежуточной волны===============================//
			const Type a = F_hll[3];			//F_E^hll
			const Type b = -U_hll[3] - F_hll[1]; // (E_hll + F_mx^hll)
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
					const Type _p = -F_hll[3] * _lambda + F_hll[1];
					//============================================================================================//

					//==========================Финальный поток HLLC=============================================//
					Vector4 _U_L;
					const Type dif_L = 1.0 / (lambda_L - _lambda);

					_U_L[0] = (U_L[0] * (lambda_L - Vel_L[0])) * dif_L;
					_U_L[1] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
					_U_L[2] = (U_L[2] * (lambda_L - Vel_L[0])) * dif_L;
					//_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0])) * dif_L;
					_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

					F = F_L + lambda_L * (_U_L - U_L);

					//============================================================================================//
				}
				else //(_S <= 0)
				{
					c3++;
					//============================Поиск промежуточного давления ===================================//
					const Type _p = -F_hll[3] * _lambda + F_hll[1];
					//============================================================================================//
					Vector4 _U_R;
					const Type dif_R = 1.0 / (lambda_R - _lambda);

					_U_R[0] = (U_R[0] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[1] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
					_U_R[2] = (U_R[2] * (lambda_R - Vel_R[0])) * dif_R;
					//_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

					F = F_R + lambda_R * (_U_R - U_R);
				}
			}
#endif
		}

		Vector4 buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[3 * num_cell + i];

	}// for


	/*U_full[num_cell] =*/ return (U - SumF * tau / volume[num_cell]);
}


Vector4 RHLLC_stepToMpi2d(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume,
	const std::vector<Vector4>&U_full_2d_prev, const std::vector<Vector4>& W_full_2d) {

	Vector4 SumF = Vector4::Zero();  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Vector4 U = U_full_2d_prev[num_cell];
	Vector4 F = Vector4::Zero();

	Vector4 U_R;
	Vector4 U_L;

	Vector4 W = W_full_2d[num_cell];
	Vector4 W_R;
	Vector4 W_L;

	Matrix4 T = Matrix4::Zero();

	for (size_t i = 0; i < 3; i++) // по граням
	{
		const int neig = neighbours_id_faces[3 * num_cell + i];
		{ // эта скобочка нужна. Т.к. далее могут быть наложения имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизация
			Type d, v, pressure;
			Vector3 vel;
			switch (neig)
			{
			case eBound_FreeBound:
				U_R = U;
				W_R = W;
				break;
			//case -10:
				//return Vector4(0, 0, 0, 0);  // плохо. Вообще не надо допускать такую ячейку
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC %d, %d\n", num_cell, neig);
					exit(1);
				}

				U_R = U_full_2d_prev[neig / 3];
				W_R = W_full_2d[neig / 3];

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
		const Vector2 Vel_L(W_L[1], W_L[2]);// , 0);  //T * velocity[num_cell];
		const Vector2 Vel_R(W_R[1], W_R[2]);// , 0);  //T * velocity[neig / 4];

		const Type d_L = W_L(0);
		const Type d_R = W_R(0);

		const Type p_L = W_L(3);
		const Type p_R = W_R(3);

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
			//F(3) = U_R[3] * Vel_R[0];
			F(3) = U_R[1];
			//continue;
			c0++; // для hllc этот код срабатывает на 8 принте для задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[0] * Vel_L[0]; //D*v_x
			F(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F(2) = U_L[2] * Vel_L[0];
			//F(3) = U_L[3] * Vel_L[0];
			F(3) = U_L[1];
			//continue;
			c1++;
		}
		else
		{
			//====================Расчёт потоков и приближений hll=========================================//
			Vector4 F_L;
			Vector4 F_R;
			Vector4 U_hll;
			Vector4 F_hll;

			F_R(0) = U_R[0] * Vel_R[0]; //D*v_x
			F_R(1) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
			F_R(2) = U_R[2] * Vel_R[0];
			//F_R(3) = U_R[3] * Vel_R[0];
			F_R(3) = U_R[1];

			F_L(0) = U_L[0] * Vel_L[0]; //D*v_x
			F_L(1) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
			F_L(2) = U_L[2] * Vel_L[0];
			//F_L(3) = U_L[3] * Vel_L[0];
			F_L(3) = U_L[1];

			F_hll = (lambda_R * F_L - lambda_L * F_R + (lambda_R * lambda_L * (U_R - U_L))) / (lambda_R - lambda_L);
			U_hll = (lambda_R * U_R - lambda_L * U_L + (F_L - F_R)) / (lambda_R - lambda_L);

#ifdef ONLY_RHLL
			F = F_hll;
#endif

			//============================================================================================//
#ifndef ONLY_RHLL		
//=========================Поиск скорости промежуточной волны===============================//
			const Type a = F_hll[3];			//F_E^hll
			const Type b = -U_hll[3] - F_hll[1]; // (E_hll + F_mx^hll)
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
					const Type _p = -F_hll[3] * _lambda + F_hll[1];
					//============================================================================================//

					//==========================Финальный поток HLLC=============================================//
					Vector4 _U_L;
					const Type dif_L = 1.0 / (lambda_L - _lambda);

					_U_L[0] = (U_L[0] * (lambda_L - Vel_L[0])) * dif_L;
					_U_L[1] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
					_U_L[2] = (U_L[2] * (lambda_L - Vel_L[0])) * dif_L;
					//_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0])) * dif_L;
					_U_L[3] = (U_L[3] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

					F = F_L + lambda_L * (_U_L - U_L);

					//============================================================================================//
				}
				else //(_S <= 0)
				{
					c3++;
					//============================Поиск промежуточного давления ===================================//
					const Type _p = -F_hll[3] * _lambda + F_hll[1];
					//============================================================================================//
					Vector4 _U_R;
					const Type dif_R = 1.0 / (lambda_R - _lambda);

					_U_R[0] = (U_R[0] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[1] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
					_U_R[2] = (U_R[2] * (lambda_R - Vel_R[0])) * dif_R;
					//_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0])) * dif_R;
					_U_R[3] = (U_R[3] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

					F = F_R + lambda_R * (_U_R - U_R);
				}
			}
#endif
		}

		Vector4 buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[3 * num_cell + i];

	}// for


	/*U_full[num_cell] =*/ return (U - SumF * tau / volume[num_cell]);
}


static int CheckState(const std::vector<Vector4>& W)
{
	for (size_t i = 0; i < size_grid; i++)
	{
		const Type d = W[i](0);
		const Vector3 v(W[i](1) , W[i](2), 0);		
		const Type p = W[i](3);

		if (d < 0 || isnan(d) || p < 0 || isnan(p) || v[0] < -0.001)
		{
			printf("cell = %d, p= %lf, d=%lf, v=%lf\n", i, p, d, v[0]);
			exit(1);
		}

	}

	static int cnt = 0;
	printf("Good check: %d\n", cnt++);

}

static int ReBuildPhysicValue(const Vector4& U, Vector4& W)
{
	Vector2 v(W(1), W(2));

	const Type vv = v.dot(v);
	const Type d = W(0);
	Type Gamma0 = 1. / sqrt(1 - vv);
	const Type h = 1 + gamma_g * W(3) / d;

	Type W0 = d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector2 m(U[1], U[2]);
	Type mm = m.dot(m);

	Type p = W(3);


	Type D = U[0];
	Type E = U[3];

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
		W0 -= (fW / dFdW);

		Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v[0] = m[0] / W0;
		v[1] = m[1] / W0;
		//v[2] = m[2] / W0;

		err -= W0;
		cc++;
	} while (fabs(err / W0) > 1e-14);

	if (p < 0 || D < 0 || std::isnan(p) || std::isnan(D))
	{
		printf("Error (p = %lf, d= %lf)",  p, D / Gamma0);
		return 1;
	}

	W(0) = D / Gamma0;
	W(1) = v(0);
	W(2) = v(1);
	W(3) = p;

	return 0;
}

static int ReBuildPhysicValue(const  std::vector<Vector4>& U, std::vector<Vector4>& W) {
	
	bool flag = false;
#pragma omp parallel default(none) shared(U, W, flag)
	{	
		const int size = U.size();

#pragma omp for
		for (int num_cell = 0; num_cell < size; num_cell++)
		{
			if (!flag && ReBuildPhysicValue(U[num_cell], W[num_cell]))
			{
#pragma omp critical
				{
					flag = true;
					printf("Error cell= %d\n", num_cell);
					//MPI_RETURN(1);
				}
			}
		}
	}

	return flag;
}

static int ReBuildPhysicValueold(const  std::vector<Vector4>& U, std::vector<Vector4>& W) {

	for (size_t num_cell = 0; num_cell < size_grid; num_cell++)
	{
		Vector2 v(W[num_cell](1), W[num_cell](2));

		const Type vv = v.dot(v);
		const Type d = W[num_cell](0);
		Type Gamma0 = 1. / sqrt(1 - vv);
		const Type h = 1 + gamma_g * W[num_cell](3) / d;

		Type W0 = d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

		Vector2 m(U[num_cell][1], U[num_cell][2]);
		Type mm = m.dot(m);

		Type p = W[num_cell](3);


		Type D = U[num_cell][0];
		Type E = U[num_cell][3];

		int  cc = 0;

		Type err = 1;
		do
		{
			err = W0;

			Type fW = W0 - p - E;

			Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
			Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
			W0 -= (fW / dFdW);

			Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

			p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

			v[0] = m[0] / W0;
			v[1] = m[1] / W0;
			//v[2] = m[2] / W0;

			err -= W0;
			cc++;
		} while (fabs(err / W0) > 1e-14);

		if (p < 0 || D < 0 || std::isnan(p) || std::isnan(D))
		{
			printf("Error cell %d (p = %lf, d= %lf)", num_cell, p, D / Gamma0);
			exit(1);
		}



		W[num_cell](0) = D / Gamma0;
		W[num_cell](1) = v(0);
		W[num_cell](2) = v(1);
		W[num_cell](3) = p;

	}

	return 0;
}
static int ReBuildConvValue(const std::vector<Vector4>& W, std::vector<Vector4>& U) {

	U.resize(size_grid);
	Vector4 cell;
	
	for (size_t i = 0; i < size_grid; i++)
	{
		const Type v = (W[i](1) * W[i](1) + W[i](2) * W[i](2));
		const Type d = W[i](0);
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * W[i](3) / d;
		const Type dhGG = d * h * Gamma * Gamma;

		cell[0] = Gamma * d;
		cell[1] = dhGG * W[i][1];
		cell[2] = dhGG * W[i][2];		

		cell[3] = dhGG - W[i](3);
		
		U[i] = cell;
	}

	return 0;
}


int RHLLC2d(std::string& main_dir,
	const std::vector<Vector3>& centerts, std::vector<int>& neighbours_id_faces,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume) {

	printf("start 2d rhllc task\n");

#ifdef Jet_2d
	Type t = 0;
	const Type T = 100;
	const Type h = 0.06610106322896746; // jet
	Type tau = h *h* 0.1;
	Type print_time = 1;
	Type cur_time = 10;
#else
	Type t = 0;
	const Type T = 0.5;	
	Type tau = 1e-4;
	Type print_time = 0.01;
	Type cur_time = 10;
#endif

	RHLLC_Init(size_grid, centerts, W_full_2d);

	InitNeigh(size_grid, centerts, normals, neighbours_id_faces);

	U_full_2d.resize(size_grid);
	U_full_2d_prev.resize(size_grid);

	ReBuildConvValue(W_full_2d, U_full_2d_prev);

	int count = 0;
	int sol_cnt = 0;	

	while (t < T)
	{
		if (cur_time >= print_time)
		{
			if (WriteSolution(sol_cnt, main_dir + "Solve\\", W_full_2d))
			{
				return 1;
			}
			cur_time = 0;
			printf("File sol: %d. time= %lf\n", sol_cnt++, t);
		}

		//ReBuildConvValue(W_full_2d, U_full_2d_prev);

#pragma omp parallel default(none) shared(size_grid, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_2d, U_full_2d_prev,W_full_2d)
		{
			Vector4 buf;
#pragma omp for
			for (int i = 0; i < size_grid; i++)
			{
				U_full_2d[i] = RHLLC_stepToOMP2d(i, tau, neighbours_id_faces, normals, squares_cell, volume);

				//buf[1] = Vector3(buf[1], buf[2], buf[3]).norm();
				//buf[2] = 0;
				//buf[3] = 0;
			}
		}


		ReBuildPhysicValue(U_full_2d, W_full_2d);
		//CheckState(U_full_2d);

#ifdef WRITE_LOG
		ofstream ofile;
		ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
		ofile << "\nt= " << t << "; tau= " << tau << "; step= " << count << '\n';
		ofile.close();
#endif

#ifdef WRITE_LOG_NO
		{
			ofstream ofile;
			ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);

			ofile << "U_full_2d\n\n";
			for (auto el : U_full_2d)
				ofile << el << '\n';

			ofile << "\n\n";
			ofile.close();
		}
#endif

		U_full_2d.swap(U_full_2d_prev);
		t += tau;
		cur_time += tau;
		count++;

	}//while

	return 0;
}


#define base 3// число граней
#ifdef USE_MPI

void InitMPI_Struct()
{

}

int MPI_RHLLC(std::string& main_dir,
	 std::vector<Vector3>& centerts,  std::vector<int>& neighbours_id_faces,
	 std::vector<Normals>& normals,  std::vector<Type>& squares_cell,  std::vector<Type>& volume)
{

	int np, myid;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef WRITE_LOG
	{		
		remove((main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt").c_str());
	}
#endif

	//np = 3; myid = 2;

	Vector4 a;	
	MPI_Datatype MPI_EigenVector4;
	int len[2] = { base + 1, 1 };
	const int offset = ((uint32_t*)a.data()) - ((uint32_t*)&a);
	MPI_Aint pos[2] = { offset ,sizeof(Vector4) };
	MPI_Datatype typ[2] = { MPI_DOUBLE, MPI_UB };
	MPI_Type_struct(2, len, pos, typ, &MPI_EigenVector4);
	MPI_Type_commit(&MPI_EigenVector4);
	//
	printf("Run\n");

	std::vector<Vector4> W_full_2d;
	std::vector<Vector4> U_full_2d;
	std::vector<Vector4> U_full_2d_prev;
	
	std::vector<Vector4> W_full_2d_local;
	std::vector<Vector4> U_full_2d_local;
	std::vector<Vector4> U_full_2d_prev_local;


	// вычисление сдвигов и числа тел на узел
	std::vector<int> send_count;
	std::vector<int> disp;
	//if(np > 1)
	{
		const int npp = np > 1 ? np - 1 : 1;
		send_count.resize(np, size_grid / npp); send_count[0] = 0;
		disp.resize(np, 0);
		for (int i = 2; i < np; i++)
			disp[i] = (i-1) * (size_grid / npp);

		//int wtf = size_grid % npp;
		if (size_grid % npp != 0) // если число процессов не кратно размерности задачи 
		{ 
			for (int i = 0; i < size_grid % npp; i++) // первые процессы берут на единицу больше тел
				++send_count[i+1];

			for (int i = 2; i < np; i++) // смещения для процессов за ними увеличивается начиная со второго
				++disp[i];
		}
	}

#ifdef WRITE_LOG
	std::vector<int> orig_id_cell_nodes;
	{
		

		//if (myid == 0)
		{
			printf("disp= ");
			for (auto el : disp)
				printf("%d ", el);

			printf("\nsend_count= ");
			for (auto el : send_count)
				printf("%d ", el);
			printf("\n");

		}
	}

	printf("np= %d, myid= %d\n", np, myid);
	//np = 4, disp= 0 0 5 8 send_count = 0 5 4 4
#endif
	
	std::vector<int> neighbours_id_faces_local;
	std::vector<Normals> normals_local;
	std::vector<Type> squares_cell_local;
	std::vector<Type> volume_local;


	if (myid == 0)
	{
		// формируем массивы на управляющем узле
		RHLLC_Init(size_grid, centerts, W_full_2d);

		U_full_2d.resize(size_grid);
		U_full_2d_prev.resize(size_grid);

		ReBuildConvValue(W_full_2d, U_full_2d_prev);
	}
	else
	{
		neighbours_id_faces_local.resize(base * send_count[myid]);
		normals_local.resize(send_count[myid]);
		squares_cell_local.resize(base * send_count[myid]);
		volume_local.resize(send_count[myid]);

		U_full_2d_local.resize(send_count[myid]);
		U_full_2d_prev_local.resize(send_count[myid]);
		W_full_2d_local.resize(send_count[myid]);
	}	
		
		/// На других узлах переменные так normals... так же будут заполнены. На каждом узле в отдельности 
		// мы проводим перенумерацию и затем на всех, кроме управляющего удаляем исходные массивы
	if (myid != 0 && np > 1)
	{
		const int size_node = send_count[myid];
		const int shift_node = disp[myid];
		int id = 0;

		for (int i = shift_node; i < shift_node+size_node; i++)
		{
			Normals nn(base);

			for (int k = 0; k < base; k++)
			{
				const int neigh = neighbours_id_faces[base * i + k];
				const int neigh_cell = neigh / base;
				if (neigh >= 0)
				{
					const int new_idx_cell = neigh_cell - shift_node;
					const int new_idx_face = new_idx_cell * base + neigh % base;

					if ((neigh_cell > shift_node) && (neigh_cell < shift_node + size_node)) // данные на этом узле
					{
						neighbours_id_faces_local[base * id + k] = new_idx_face;
						orig_id_cell_nodes.push_back(i);
					}
					else
					{
						neighbours_id_faces_local[base * id + k] = -10;  // флаг на то, что данные на другом узле
						//orig_id_cell_nodes.erase(orig_id_cell_nodes.end() - k, orig_id_cell_nodes.end());
							//break;
					}
				}
				else
				{
					neighbours_id_faces_local[base * id + k] = neigh;
					orig_id_cell_nodes.push_back(i);
				}

				nn.n[k] = normals[i].n[k];
				squares_cell_local[base * id + k] = squares_cell[base * i + k];
			}
			normals_local[id] = nn;
			volume_local[id] = volume[i];

			id++;
		}
		auto last = std::unique(orig_id_cell_nodes.begin(), orig_id_cell_nodes.end());
		orig_id_cell_nodes.erase(last, orig_id_cell_nodes.end());
		// очищаем память на вычислительном узле
		neighbours_id_faces.clear();
		normals.clear();
		squares_cell.clear();
		volume.clear();
	
	}

	 // поиск ячеек, которые должны быть расчитаны на главном узле (один раз)
	 //--------------------------------------------до while(T)_----
	 std::vector<int> id_cells_node0;	 	 
	 if (myid == 0 && np > 1) // управляющий узел
	 {
		 id_cells_node0.reserve(size_grid);	

		 for (size_t id = 1; id < np; id++) // по всем узлам
		 {
			 const int shift_node = disp[id];
			 const int size_node = send_count[id];

			 for (int i = shift_node; i < shift_node+ size_node; i++)// все ячеки узла
			 {
				 for (size_t j = 0; j < base; j++) // все грани
				 {
					 const int neigh_face = neighbours_id_faces[i * base + j]; // -> грань
					 if (neigh_face >= 0)
					 {
						 const int neigh_cell = neigh_face / base; // -> ячейка
						 if ((neigh_cell <= shift_node) || (neigh_cell >= shift_node + size_node)) // ячейка не на данном узле					 
						 {
							 id_cells_node0.push_back(i);
							 orig_id_cell_nodes.push_back(i);
							 break;
						 }
					 }
				 }				 
			 }
		 }
		 
		 std::sort(id_cells_node0.begin(), id_cells_node0.end());
		 auto last = std::unique(id_cells_node0.begin(), id_cells_node0.end());  // могуть быт повторы	
		 id_cells_node0.erase(last, id_cells_node0.end());

	 }
	 else if (np == 1)
	 {
		 id_cells_node0.resize(size_grid);
		 for (int i = 0; i < size_grid; i++)
		 {
			 id_cells_node0[i] = i;  // вся сетка
		 }		 
	 }

	 if (myid == 0)
	 {
		 U_full_2d_local.resize(id_cells_node0.size());
		 W_full_2d_local.resize(id_cells_node0.size());
	 }
	 //--------------------------------------------до while(T)_--------------------

	
	
#ifdef WRITE_LOG
	 {
		 ofstream ofile;
		 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);

		 if (myid == 0)
		 {
			 ofile << "neighbours_id_faces= " << id_cells_node0.size() << "\n\n";
			 for (auto& el : id_cells_node0)
				 ofile << el << ' ';
			 ofile << "\n\n";
		 }
		 else
		 {
			 ofile << "neighbours_id_faces= "<< neighbours_id_faces_local.size() << "\n\n";
			 for (auto& el : neighbours_id_faces_local)
				 ofile << el/base << ' ';
			 ofile << "\n\n";
		 }
		 		 
		 ofile.close();
	 }
#endif
	 // Расчет

	 Type t = 0;
	 const Type T = 0.5;	 
	 Type print_time = 0.01;
	 Type cur_time = 10;	 
	 Type tau = 0;

	 // ------------------------------------------ Основной расчёт------------------------------------------
	 int sol_cnt = 0;
	 int count = 0;

	 while (t < T)
	 {

		 {
			 // и себе тоже-> блок себя удалить, если np!= 1? 
			 MPI_Scatterv(W_full_2d.data(), send_count.data(), disp.data(), MPI_EigenVector4, W_full_2d_local.data(), send_count[myid], MPI_EigenVector4, 0, MPI_COMM_WORLD);
			 MPI_Scatterv(U_full_2d_prev.data(), send_count.data(), disp.data(), MPI_EigenVector4, U_full_2d_prev_local.data(), send_count[myid], MPI_EigenVector4, 0, MPI_COMM_WORLD);
			 //MPI_Scatterv(U_full_2d.data(), send_count.data(), disp.data(), MPI_EigenVector4, U_full_2d_local.data(), send_count[myid], MPI_EigenVector4, 0, MPI_COMM_WORLD);
		 }

#ifndef WRITE_LOG
		 {
			 ofstream ofile;
			 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);

			 //if (myid == 0)
			 {
				 ofile << "U_full_2d_prev_local= " << U_full_2d_prev_local.size() << "\n\n";
				 for (auto& el : U_full_2d_prev_local)
				 {
					 ofile << el[0] << '\n';
				 }
				 ofile << "\n\n";
			 }
			 ofile.close();

			 int ii = 0;
			 if (myid != 0)
			 {
				 for (size_t i = 0; i < U_full_2d_local.size(); i++)
				 {
					 if (neighbours_id_faces_local[i * base + 0] == -10 ||
						 neighbours_id_faces_local[i * base + 1] == -10 ||
						 neighbours_id_faces_local[i * base + 2] == -10) continue;
					 
						 int id = neighbours_id_faces_local[i] / base;
						 U_full_2d_local[i] = -U_full_2d_prev_local[i];
					 
				 }				 
			 }
			 else
			 {
				 for (size_t i = 0; i < id_cells_node0.size(); i++)
				 {
					 U_full_2d_local[i][0] = -id_cells_node0[i];
				 }
			 }
			 MPI_Gatherv(U_full_2d_local.data(), send_count[myid], MPI_EigenVector4, U_full_2d.data(), send_count.data(), disp.data(), MPI_EigenVector4, 0, MPI_COMM_WORLD);

			 if (myid == 0)
			 {
				 for (size_t i = 0; i < id_cells_node0.size(); i++)
				 {
					 U_full_2d[id_cells_node0[i]][0] = U_full_2d_local[i][0];
				 }
			 }
			 //ofstream ofile;
			 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);

			 if (myid == 0)
			 {
				 ofile << "U_full_2d= " << U_full_2d.size() << "\n\n";
				 for (auto& el : U_full_2d)
				 {
					 ofile << el[0] << '\n';
				 }
				 ofile << "\n\n";
			 }
			 ofile.close();
		 }
#endif


		 if (myid == 0)
		 {
			 tau = 1e-4;
		 }

		 MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  // заготовка под переменный шаг
		 
		 int cc = 0;

		 if (myid == 0) // на главном узле хранится полный объем данных
		 {

			 // будет время от других узлов
			 if (cur_time >= print_time)
			 {
				 if (WriteSolution(sol_cnt, main_dir + "Solve\\", W_full_2d)) //W- не перезаписана!
				 {
					 MPI_RETURN(1);
				 }
				 cur_time = 0;
				 printf("File sol: %d. time= %lf\n", sol_cnt++, t);
#ifdef WRITE_LOG
				 {
					 ofstream ofile;
					 ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);

					 ofile << "t= " << t << " sol_cnt= " << sol_cnt << '\n';
					 ofile.close();
				 }
#endif
			 }

			 const int N = id_cells_node0.size();
			 for (int id = 0; id < N; id++)// можно omp
			 {
				 const int i = id_cells_node0[id];

				 U_full_2d_local[id] = RHLLC_stepToMpi2d(i, tau, neighbours_id_faces, normals, squares_cell, volume,
					 U_full_2d_prev, W_full_2d);
				 cc++;

				 if (ReBuildPhysicValue(U_full_2d_local[id], W_full_2d[i]))
				 {
#ifdef WRITE_LOG
					 {
						 ofstream ofile;
						 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);
						 ofile << "Err cell " << i << '\n';
						 ofile.close();
					 }
#endif
					 MPI_RETURN(1);
				 }
				 W_full_2d_local[id] = W_full_2d[i];
			 }
			 
		 }
		 else // вычислительные узлы
		 {
			 const int size = U_full_2d_local.size();
			 for (int i = 0; i < size; i++)
			 {
				 
				 if (neighbours_id_faces_local[i*base] == -10 || neighbours_id_faces_local[i * base +1] == -10
					 || neighbours_id_faces_local[i * base +2] == -10) // другой узел
				 {
					 continue; // этот узел будет выпоолнен на основном узле
				 }
				 // U_full_2d_prev,  W_full_2d -> local
				 U_full_2d_local[i] = RHLLC_stepToMpi2d(i, tau, neighbours_id_faces_local, normals_local, squares_cell_local, volume_local,
					 U_full_2d_prev_local, W_full_2d_local);

				 cc++;
				 if (ReBuildPhysicValue(U_full_2d_local[i], W_full_2d_local[i]))
				 {
#ifdef WRITE_LOG
					 {
						 ofstream ofile;
						 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);
						 ofile << "Err cell " << i << '\n';
						 ofile.close();
					 }
#endif
					 MPI_RETURN(1);
				 }

			 }
		 }

		 static int one_print = 0;
		 if (!one_print)
		 {
			 printf("[%d] CC=%d\n", myid, cc);
			 one_print = 1;
		 }
		 //if (myid != 0 && np > 1)
		 {
			 MPI_Gatherv(W_full_2d_local.data(), send_count[myid], MPI_EigenVector4, W_full_2d.data(), send_count.data(), disp.data(), MPI_EigenVector4, 0, MPI_COMM_WORLD);
			 //MPI_Gatherv(U_full_2d_prev_local.data(), send_count[myid], MPI_EigenVector4, U_full_2d_prev.data(), send_count.data(), disp.data(), MPI_EigenVector4, 0, MPI_COMM_WORLD);
			 MPI_Gatherv(U_full_2d_local.data(), send_count[myid], MPI_EigenVector4, U_full_2d.data(), send_count.data(), disp.data(), MPI_EigenVector4, 0, MPI_COMM_WORLD);

		 }

		 if (myid == 0)
		 {
			 int ii = 0;
			 for (auto id : id_cells_node0)
			 {
				 U_full_2d[id] = U_full_2d_local[ii];  // сохраняем ячейки расчитаные отдельно
				 W_full_2d[id] = W_full_2d_local[ii++];  // сохраняем ячейки расчитаные отдельно
			 }			 

			 std::swap(U_full_2d, U_full_2d_prev);
			 //U_full_2d.swap(U_full_2d_prev); // здесь он переставит указатели, mpi может сойти с ума
		 }

#ifndef WRITE_LOG
		 {
			 ofstream ofile;
			 ofile.open(main_dir + "File_with_Logs_solve" + to_string(myid) + ".txt", std::ios::app);
			

			 ofile<<"t= "<<t << " U_full_2d_local= " << U_full_2d_local.size() << "\n\n";
			 for (auto& el : U_full_2d_local)
				 ofile << el << '\n';
			 ofile << "\n\n";

			 ofile << "\n\n";
			 ofile.close();
	 }
#endif

		 t += tau;
		 cur_time += tau;
		 count++;
	 }// while

	MPI_RETURN(0);
}

#endif
#endif // USE_VTK