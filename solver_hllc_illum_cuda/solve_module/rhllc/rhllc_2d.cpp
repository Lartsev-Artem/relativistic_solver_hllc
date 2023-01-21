#include "../solve_config.h"
#if defined RHLLC && NUMBER_OF_MEASUREMENTS == 2 && !defined USE_MPI
#include "../solve_global_struct.h"
#include "../../file_module/writer_bin.h"
#include "../../utils/grid_geometry/geometry_solve.h"
#include "rhllc_utils.h"

//#define SODA_2d
#define Jet_2d

static int c0 = 0;
static int c1 = 0;
static int c2 = 0;
static int c3 = 0;

static std::vector<Vector4> W_full_2d;
static std::vector<Vector4> U_full_2d;
static std::vector<Vector4> U_full_2d_prev;

static  int WriteSolution(const int n, const int size_grid, const std::string main_dir, const std::vector<Vector4>& W)
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

	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "density.bin", den);
	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "pressure.bin", press);
	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "velocity.bin", vel);	

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

static Vector4 RHLLC_stepToOMP2d(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
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

	for (size_t i = 0; i < 3; i++) // по гран€м
	{
		const int neig = neighbours_id_faces[3 * num_cell + i];
		{ // эта скобочка нужна. “.к. далее могут быть наложени€ имен переменных. 
			// todo: переписать lock bound и free_bound на функции + оптимизаци€
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

		//====================  эшируем физические переменные слева и справа============================//
		// нормальна€ сокорость
		const Vector2 Vel_L(W_L[1], W_L[2]);// , 0);  //T * velocity[num_cell];
		const Vector2 Vel_R(W_R[1], W_R[2]);// , 0);  //T * velocity[neig / 4];

		const Type d_L = W_L(0);
		const Type d_R = W_R(0);

		const Type p_L = W_L(3);
		const Type p_R = W_R(3);

		const Type VV_L = Vel_L.dot(Vel_L);
		const Type VV_R = Vel_R.dot(Vel_R);

		//========================================================================================//

		//=========================¬ычисл€ем рел€тивистикие параметры============================//				
		const Type g_L = 1. / sqrt(1 - VV_L);	// фактор Ћоренца
		const Type g_R = 1. / sqrt(1 - VV_R);

		const Type h_L = 1 + gamma_g * p_L / d_L; // энтальпи€
		const Type h_R = 1 + gamma_g * p_R / d_R;

		const Type cs_L = sqrt((gamma1 * p_L) / (d_L * h_L)); // скорость звука
		const Type cs_R = sqrt((gamma1 * p_R) / (d_R * h_R));

		const Type sigmaS_L = (cs_L * cs_L) / (g_L * g_L * (1 - cs_L * cs_L)); // что-то дл€ расчета собственных чисел HHL
		const Type sigmaS_R = (cs_R * cs_R) / (g_R * g_R * (1 - cs_R * cs_R));

		//========================================================================================//

		const Type sqr_L = sqrt(sigmaS_L * (1 - Vel_L[0] * Vel_L[0] + sigmaS_L));
		const Type sqr_R = sqrt(sigmaS_R * (1 - Vel_R[0] * Vel_R[0] + sigmaS_R));

		// здесь встречалась альтернатива сравнени€ с нулем min(0,L), max(0,R)
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
			c0++; // дл€ hllc этот код срабатывает на 8 принте дл€ задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо дл€ всех границ
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
			//====================–асчЄт потоков и приближений hll=========================================//
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
//=========================ѕоиск скорости промежуточной волны===============================//
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
					//============================ѕоиск промежуточного давлени€ ===================================//
					const Type _p = -F_hll[3] * _lambda + F_hll[1];
					//============================================================================================//

					//==========================‘инальный поток HLLC=============================================//
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
					//============================ѕоиск промежуточного давлени€ ===================================//
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

int RHLLC2d(std::string& solve_dir,
	const std::vector<Vector3>& centerts, std::vector<int>& neighbours_id_faces,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume) 
{

	printf("start 2d rhllc task\n");

	const int size_grid = centerts.size();

#ifdef Jet_2d
	Type t = 0;
	const Type T = 100;
	const Type h = 0.06610106322896746; // jet
	Type tau = h * h * 0.1;
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
			if (WriteSolution(sol_cnt, size_grid, solve_dir, W_full_2d))
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

		WRITE_LOG("\nt= " << t << "; tau= " << tau << "; step= " << count << '\n');

		U_full_2d.swap(U_full_2d_prev);
		t += tau;
		cur_time += tau;
		count++;

	}//while

	return 0;
}



#endif // RHLLC2d