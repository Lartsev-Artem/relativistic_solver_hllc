#include "../solve_config.h"
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 2 && !defined USE_MPI
#include "../solve_global_struct.h"
#include "../../file_module/writer_bin.h"
#include "../../utils/grid_geometry/geometry_solve.h"
#include "hllc_utils.h"


#define JET_2D

static int c0 = 0;
static int c1 = 0;
static int c2 = 0;
static int c3 = 0;

//static std::vector<Vector4> W_2d;
static std::vector<Vector4> U_full_2d;
static std::vector<Vector4> U_full_2d_prev;

static int WriteSolution(const int n, const std::string main_dir, const std::vector<Vector4>& U_full)
{
	std::vector<Type> den(size_grid);
	std::vector<Type> press(size_grid);
	std::vector<Vector3> vel(size_grid);

	for (size_t i = 0; i < size_grid; i++)
	{
		const Type d = U_full[i](0);
		den[i] = d;

		vel[i](0) = U_full[i](1) / d;
		vel[i](1) = U_full[i](2) / d;		
		vel[i](2) = 0;

		const Type v = vel[i].norm();
		press[i] = (U_full[i](3) - v * v * d / 2.) * (gamma1 - 1);
	}

	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "density.bin", den);
	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "pressure.bin", press);
	WriteSimpleFileBin(main_dir + "Solve" + std::to_string(n) + "velocity.bin", vel);


	return 0;
}

static int ReBuildDataForHLLC2d(const int N, const std::vector<Vector3>& centerts, std::vector<Vector4>& data) {
	
	data.resize(N);
	Vector4 cell;

	for (size_t i = 0; i < N; i++)
	{		
		Vector2 x(centerts[i][0], centerts[i][1]);
#ifdef JET_2D
		if (x[0] < 1 && x[1] < 1 && x[1]>-1)
		{
			const Type d = 0.1;
			const Vector3 v(0.99, 0, 0);
			const Type p = 0.01;
			
			cell[0] = d;
			cell[1] = d * v[0];
			cell[2] = d * v[1];

			const Type vv = v.dot(v);
			cell[3] = p / (gamma1 - 1) + d * vv / 2;
		}
		else
		{
			const Type d = 10;
			const Vector3 v(0, 0, 0);
			const Type p = 0.01;

			cell[0] = d;
			cell[1] = d * v[0];
			cell[2] = d * v[1];

			const Type vv = v.dot(v);
			cell[3] = p / (gamma1 - 1) + d * vv / 2;
		}
#else

		if (x[0] < 0.5)
		{
			const Type d = 1;
			const Vector3 v(0, 0, 0);
			const Type p = 1;

			cell[0] = d;
			cell[1] = d * v[0];
			cell[2] = d * v[1];

			const Type vv = v.dot(v);
			cell[3] = p / (gamma1 - 1) + d * vv / 2;
		}
		else
		{
			const Type d = 0.125;
			const Vector3 v(0, 0, 0);
			const Type p = 0.1;

			cell[0] = d;
			cell[1] = d * v[0];
			cell[2] = d * v[1];

			const Type vv = v.dot(v);
			cell[3] = p / (gamma1 - 1) + d * vv / 2;
		}
#endif // JET_2D
		data[i] = cell;
	}

	return 0;
}

static inline void MakeRotationMatrix(const Vector3& n, Matrix4& T) {

	// n=(x,y,0)!!!!

	T = Matrix4::Zero();
	T(0, 0) = T(3, 3) = 1; 

	T(1, 1) = n[0];
	T(1, 2) = n[1];

	T(2, 1) = -n[1];
	T(2, 2) = n[0];

}

static Vector4 HLLC_stepToOMP2d(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume) {

	Vector4 SumF = Vector4::Zero();  // интеграл по поверхности от F (т.е. Sum{F*n*dS}
	const Vector4 U = U_full_2d_prev[num_cell];
	Vector4 F = Vector4::Zero();

	Vector4 U_R;
	Vector4 U_L;

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
				break;		
			case eBound_OutSource://eBound_InnerSource:
#ifdef JET_2D	
				 d = 0.1;
				 vel << 0.99, 0, 0;
				 pressure = 0.01;
				v = vel.dot(vel);				
				U_R << d, d* vel[0], d* vel[1], pressure / (gamma1 - 1) + d * v / 2;
					
#else
				U_R = U;				
#endif
				break;
			default:
				if (neig < 0)
				{
					printf("Err bound in HLLC\n");
					exit(1);
				}

				U_R = U_full_2d_prev[neig / 3];
				break;
			}
		}
			
		MakeRotationMatrix(normals[num_cell].n[i], T);
		// здесь одна матрица поворота или с разным знаком?????
		U_L = T * U;
		U_R = T * U_R;

		Type d_L = U_L[0]; //density[num_cell];
		Type d_R = U_R[0]; // density[neig];

		Vector2 vel(U_L(1) / d_L, U_L(2) / d_L);
		Type v = vel.dot(vel);
		Type p_L = (U_L(3) - v * d_L / 2.) * (gamma1 - 1);

		vel << U_R(1) / d_R, U_R(2) / d_R;
		v = vel.dot(vel);
		Type p_R = (U_R(3) - v * d_R / 2.) * (gamma1 - 1);


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
		
		if (S_R <= 0) // если верно выполнить всегда
		{
			F(0) = U_R[1];
			F(1) = U_R[1] * U_R[1] / d_R + p_R;
			F(2) = U_R[1] * U_R[2] / d_R;			
			F(3) = (U_R[3] + p_R) * U_R[1] / d_R;
			//continue;
			c0++;
		}
		else if (S_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F(0) = U_L[1];
			F(1) = U_L[1] * U_L[1] / d_L + p_L;
			F(2) = U_L[1] * U_L[2] / d_L;			
			F(3) = (U_L[3] + p_L) * U_L[1] / d_L;  // TU[4]*d_L
			//continue;
			c1++;
		}
		else
		{
			const Type roS_L = d_L * (S_L - v_L);
			const Type roS_R = d_R * (S_R - v_R);
			const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

			const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;

			Vector4 D; D << 0, 1, 0, _S;

			if (_S >= 0)
			{
				F(0) = U_L[1];
				F(1) = U_L[1] * U_L[1] / d_L + p_L;
				F(2) = U_L[1] * U_L[2] / d_L;				
				F(3) = (U_L[3] + p_L) * U_L[1] / d_L;

				Vector4 buf = F;
				F = (_S * (S_L * U_L - buf) + S_L * P_LR * D) / (S_L - _S);
				c2++;
			}
			else //(_S <= 0)
			{
				F(0) = U_R[1];
				F(1) = U_R[1] * U_R[1] / d_R + p_R;
				F(2) = U_R[1] * U_R[2] / d_R;				
				F(3) = (U_R[3] + p_R) * U_R[1] / d_R;

				Vector4 buf = F;
				F = (_S * (S_R * U_R - buf) + S_R * P_LR * D) / (S_R - _S);
				c3++;
			}
		}

		Vector4 buf = F;
		F = (T.transpose()) * buf;

		SumF += F * squares_cell[3 * num_cell + i];

	}// for


	/*U_full[num_cell] =*/ return (U - SumF * tau / volume[num_cell]);
}

static int InitNeigh(const int N, const std::vector<Vector3>& centerts, const std::vector<Normals>& normals,
	std::vector<int>& neighbours_id_faces)
{
#ifdef JET_2D
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

int HLLC_2d(std::string& solve_dir,
	const std::vector<Vector3>& centerts, std::vector<int>& neighbours_id_faces,
	const std::vector<Normals>& normals, const std::vector<Type>& squares_cell, const std::vector<Type>& volume) {

	printf("start 2d task\n");
	const int size_grid = centerts.size();

	Type t = 0;
#ifdef JET_2D
	const Type T = 100;
	const Type h = 0.06610106322896746; // jet
	Type tau = h * h * 0.1;
	Type print_time = 1;
	Type cur_time = 10;
#else	
	const Type T = 0.4;
	Type tau = 1e-4;
	Type print_time = 0.01;
	Type cur_time = 10;
#endif


	ReBuildDataForHLLC2d(size_grid, centerts, U_full_2d_prev);
	U_full_2d.resize(size_grid);

	InitNeigh(size_grid, centerts, normals, neighbours_id_faces);

	int count = 0;
	int sol_cnt = 0;

	while (t < T)
	{
		if (cur_time >= print_time)
		{
			if (WriteSolution(sol_cnt, solve_dir, U_full_2d_prev))
			{
				return 1;
			}
			cur_time = 0;
			printf("File sol: %d\n", sol_cnt++);

			WRITE_LOG("\nt= " << t << "; tau= " << tau << "; step= " << count << '\n');			
		}

#pragma omp parallel default(none) shared(size_grid, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_2d, U_full_2d_prev)
		{
#pragma omp for
			for (int i = 0; i < size_grid; i++)
			{
				U_full_2d[i] = HLLC_stepToOMP2d(i,tau, neighbours_id_faces, normals, squares_cell, volume);

				//buf[1] = Vector3(buf[1], buf[2], buf[3]).norm();
				//buf[2] = 0;
				//buf[3] = 0;

			}
		}

		U_full_2d.swap(U_full_2d_prev);
		t += tau;
		cur_time += tau;
		count++;

	}//while

	return 0;
}


#endif // USE_VTK