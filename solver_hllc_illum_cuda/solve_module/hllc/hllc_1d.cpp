#include "../solve_config.h"
#if defined HLLC && NUMBER_OF_MEASUREMENTS == 1

#include "../solve_global_struct.h"

struct PhysVal 
{
	Type rho;
	Type v;
	Type p;
};

struct ConVal 
{
	Type d;
	Type v;
	Type p;
};

struct Flux 
{
	Type d;
	Type v;
	Type p;

	Flux() 
	{
		d = 0; v = 0; p = 0;
	}
};

static int size_g= 101;
static const Type h = 1./(size_g-1);
static Type tau = 1e-5; // h * 0.1;


static std::vector<PhysVal> W(size_g);

static std::vector<ConVal> U_full_1d(size_g);
static std::vector<ConVal> U_full_1d_prev(size_g);

static std::vector<Flux> F_full(size_g);

static int c0,c1, c2, c3, c4;

static int WriteSolution(const std::string& dir, const std::vector<ConVal>& vector_U, bool flag_init)
{
	int n = vector_U.size();
	std::vector<Type> density(n);
	std::vector<Type> pressure(n);
	std::vector<Type> velocity(n);

	for (size_t i = 0; i < n; i++)
	{
		const Type d = vector_U[i].d;
		density[i] = d;
		velocity[i] = vector_U[i].v / d;
		const Type v = velocity[i];
		pressure[i] = (vector_U[i].p - v * v * d / 2.) * (gamma1 - 1);
	}

	std::ofstream ofile_1;

	if (flag_init)
	{
		ofile_1.open(dir + "density.txt");
		ofile_1 << n << '\n';
	}
	else
		ofile_1.open((dir + "density.txt"), std::ios::app);


	if (!ofile_1.is_open())
	{
		printf("density not opened\n");
		return 1;
	}

	for (size_t i = 0; i < n; i++)
	{
		ofile_1 << density[i] << '\n';
	}
	ofile_1.close();

	if (flag_init)
	{
		ofile_1.open((dir + "velocity.txt"));
		ofile_1 << n << '\n';
	}
	else
		ofile_1.open((dir + "velocity.txt"), std::ios::app);

	if (!ofile_1.is_open())
	{
		printf("densit not opened\n");
		return 1;
	}

	for (size_t i = 0; i < n; i++)
	{
		ofile_1 << velocity[i] << '\n';
	}
	ofile_1.close();

	if (flag_init)
	{
		ofile_1.open((dir + "pressure.txt"));
		ofile_1 << n << '\n';
	}
	else
		ofile_1.open((dir + "pressure.txt"), std::ios::app);

	if (!ofile_1.is_open())
	{
		printf("densit not opened\n");
		return 1;
	}

	for (size_t i = 0; i < n; i++)
	{
		ofile_1 << pressure[i] << '\n';
	}
	ofile_1.close();

	return 0;
}

static inline double MakeRotateCoef(const int num_face)
{
	if (num_face == 0)
	{
		return -1.0; //лево
	}
	else
	{
		return 1.0; // право
	}
}

static ConVal HLLC_flux(int num_cell)
{

	Flux SumF;
	ConVal U = U_full_1d_prev[num_cell];
	Flux F;

	ConVal U_R;
	ConVal U_L;
	
	for (size_t i = 0; i < 2; i++) // по граням
	{
		const int neig = i == 0 ? num_cell - 1 : num_cell + 1;
		{ 
			//todo: поворот
			if (neig < 0) // левая граница
			{
				U_R = U;
			}
			else if (neig == size_g) // правая граница
			{
				U_R = U;
			}
			else
			{
				U_R = U_full_1d_prev[neig];
			}
		}

		const Type T = MakeRotateCoef(i);
		U_L = U;  U_L.v *= T;
		U_R = U_R; U_R.v *= T;

		Type d_L = U_L.d; //density[num_cell];
		Type d_R = U_R.d; // density[neig];

		Type vel = U_L.v / d_L; ////(U_L(1) / d_L, U_L(2) / d_L, U_L(3) / d_L);
		Type v = vel * vel;
		Type p_L = (U_L.p - v * d_L / 2.) * (gamma1 - 1);

		vel = U_R.v / d_R;
		v = vel*vel;
		Type p_R = (U_R.p - v * d_R / 2.) * (gamma1 - 1);

		const Type v_L = U_L.v / d_L; // sqrt(U_L[1] * U_L[1] + U_L[2] * U_L[2] + U_L[3] * U_L[3]);  //velocity[num_cell].norm();
		const Type v_R = U_R.v / d_R; //sqrt(U_R[1] * U_R[1] + U_R[2] * U_R[2] + U_R[3] * U_R[3]); //velocity[neig].norm();

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
			F.d = U_R.v;
			F.v = U_R.v * U_R.v / d_R + p_R;			
			F.p = (U_R.p + p_R) * U_R.v / d_R;
			//continue;
			c0++;
		}
		else if (S_L >= 0) // выполнить либо по условию либо для всех границ
		{
			F.d = U_L.v;
			F.v = U_L.v * U_L.v / d_L + p_L;		
			F.p = (U_L.p + p_L) * U_L.v / d_L;  // TU[4]*d_L
			//continue;
			c1++;
		}
		else
		{
			const Type roS_L = d_L * (S_L - v_L);
			const Type roS_R = d_R * (S_R - v_R);
			const Type _S = (p_R - p_L + roS_L * v_L - roS_R * v_R) / (roS_L - roS_R);

			const Type P_LR = (p_L + p_R + roS_L * (_S - v_L) + roS_R * (_S - v_R)) / 2.0;

			Eigen::Vector3d D(0, 1, _S);

			if (_S >= 0)
			{
				F.d = U_L.v;
				F.v = U_L.v * U_L.v / d_L + p_L;				
				F.p = (U_L.p + p_L) * U_L.v / d_L;

				Flux buf = F;
				F.d = (_S * (S_L * U_L.d - buf.d) + S_L * P_LR * D(0)) / (S_L - _S);
				F.v = (_S * (S_L * U_L.v - buf.v) + S_L * P_LR * D(1)) / (S_L - _S);
				F.p = (_S * (S_L * U_L.p - buf.p) + S_L * P_LR * D(2)) / (S_L - _S);
				c2++;
			}
			else //(_S <= 0)
			{
				F.d = U_R.v;
				F.v = U_R.v * U_R.v / d_R + p_R;				
				F.p = (U_R.p + p_R) * U_R.v / d_R;

				Flux buf = F;
				F.d = (_S * (S_R * U_R.d - buf.d) + S_R * P_LR * D(0)) / (S_R - _S);
				F.v = (_S * (S_R * U_R.v - buf.v) + S_R * P_LR * D(1)) / (S_R - _S);
				F.p = (_S * (S_R * U_R.p - buf.p) + S_R * P_LR * D(2)) / (S_R - _S);
				c3++;
			}
		}
		
		F.v *= T;  //???

		SumF.d += F.d;
		SumF.v += F.v;
		SumF.p += F.p;

	}// for


	U.d -= (SumF.d * tau / h);
	U.v -= (SumF.v * tau / h);
	U.p -= (SumF.p * tau / h);
	
	/*U_full_1d[num_cell] =*/ return U;
}

static int HLLC_Init(std::vector<ConVal>& U) {
	
	int n = U.size();	
	for (size_t i = 0; i < n; i++)
	{
		if (i * h < 0.5)
		{
			U[i].d = 1;
			U[i].v = 0;
			U[i].p = 1;
		}
		else
		{
			U[i].d = 0.125;
			U[i].v = 0;
			U[i].p = 0.1;
		}
	}
	return 0;
}

int RunHllc_1d(std::string& main_dir) 
{	

	Type t = 0;
	const Type T = 0.4;
	const Type print_time = 0.1;
	Type cur_time = 0;
	
	HLLC_Init(U_full_1d_prev);
	if (WriteSolution(main_dir, U_full_1d_prev, 1))
	{
		return 1;
	}

	int count = 0;
	while (t < T)	
	{	
		for (size_t i = 0; i < size_g; i++)
		{
			U_full_1d[i] = HLLC_flux(i);
		}		

		if (cur_time >= print_time)
		{
			if (WriteSolution(main_dir, U_full_1d, 0))
			{
				return 1;
			}
			cur_time = 0;
		}


		U_full_1d.swap(U_full_1d_prev);
		t += tau;
		cur_time += tau;
		count++;
	}

	printf("timer %d\n", count);
	return 0;
}

#endif // HLLC 1d