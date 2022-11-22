#include "solve_short_characteristics_headers.h"
#include "solve_short_characteristics_global_structure.h"

#define SIGN(a) (a < 0.0 ? -1.0 : 1.0) 
struct PhysVal {
	Type rho;
	Type v;
	Type p;
};

struct ConVal {
	Type d;
	Type v;
	Type p;
};

struct Flux {
	Type d;
	Type v;
	Type p;

	Flux()
	{
		d = 0; v = 0; p = 0;
	}
};

static int size_g = 101;
static const Type h = 1. / (size_g - 1);
static Type tau = h * 0.1;


static std::vector<PhysVal> W_full_1d(size_g);
static std::vector<PhysVal> W_full_1d_prev(size_g);

static std::vector<ConVal> U_full_1d(size_g);
static std::vector<ConVal> U_full_1d_prev(size_g);

static std::vector<Flux> F_full(size_g*2);

static int c0, c1, c2, c3, c4;

int WriteSolution(const std::string& dir, const std::vector<PhysVal>& vector_U, bool flag_init);

static double MakeRotateCoef(int num_face)
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

ConVal RHLLC_flux(const int num_cell) {
	
	Vector3 SumF;
	ConVal U = U_full_1d_prev[num_cell];
	PhysVal W = W_full_1d_prev[num_cell];

	Vector3 F;

	ConVal U_R;
	ConVal U_L;

	PhysVal W_R;
	PhysVal W_L;

	for (size_t i = 0; i < 2; i++) // по гран€м
	{
		const int neig = i == 0 ? num_cell - 1 : num_cell + 1;
		{
			//todo: поворот
			if (neig < 0) // лева€ граница
			{
				U_R = U;
				W_R = W;
			}
			else if (neig == size_g) // права€ граница
			{
				U_R = U;
				W_R = W;
			}
			else
			{
				U_R = U_full_1d_prev[neig];
				W_R = W_full_1d_prev[neig];
			}
		}

		const Type T = MakeRotateCoef(i);
		U_L = U;  U_L.v *= T;
		U_R = U_R; U_R.v *= T;

		W_L = W;  W_L.v *= T;
		W_R = W_R; W_R.v *= T;


		//====================  эшируем физические переменные слева и справа============================//
				// нормальна€ сокорость			
		const Vector3 Vel_L (W_L.v,0,0); 
		const Vector3 Vel_R (W_R.v,0,0);
		
		const Type d_L = W_L.rho;
		const Type d_R = W_R.rho;

		const Type p_L = W_L.p;
		const Type p_R = W_R.p;

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
/// \warning: в 3d случае здесь нормальна€ компонента. “ак ли это?


		// здесь встречалась альтернатива сравнени€ с нулем min(0,L), max(0,R)
		const Type lambda_L = min((Vel_L[0] - sqr_L) / (1 + sigmaS_L), (Vel_R[0] - sqr_R) / (1 + sigmaS_R));
		const Type lambda_R = max((Vel_L[0] + sqr_L) / (1 + sigmaS_L), (Vel_R[0] + sqr_R) / (1 + sigmaS_R));

		if (lambda_R <= 0) // если верно выполнить всегда
		{
			F.d = U_R.d * Vel_R; //D*v_x
			F.v = U_R.v * Vel_R + p_R; //mx*vx+p			
			F.p = U_R.v;
			//continue;
			c0++; // дл€ hllc этот код срабатывает на 8 принте дл€ задачи soda
		}
		else if (lambda_L >= 0) // выполнить либо по условию либо дл€ всех границ
		{
			F.d = U_L.d * Vel_L; //D*v_x
			F.v = U_L.v * Vel_L + p_L; //mx*vx+p
			F.p = U_L.v;
			//continue;
			c1++;
		}
		else
		{
			//====================–асчЄт потоков и приближений hll=========================================//
			Flux F_L;
			Flux F_R;
			ConVal U_hll;
			Flux F_hll;

			F_R.d = U_R.d * Vel_R; //D*v_x
			F_R.v = U_R.v * Vel_R + p_R; //mx*vx+p			
			F_R.p = U_R.v;

			F_L.d = U_L.d * Vel_L; //D*v_x
			F_L.v = U_L.v * Vel_L + p_L; //mx*vx+p			
			F_L.p = U_L.v;

			F_hll.d = (lambda_R * F_L.d - lambda_L * F_R.d + (lambda_R * lambda_L * (U_R.d - U_L.d))) / (lambda_R - lambda_L);
			F_hll.v = (lambda_R * F_L.v - lambda_L * F_R.v + (lambda_R * lambda_L * (U_R.v - U_L.v))) / (lambda_R - lambda_L);
			F_hll.p = (lambda_R * F_L.p - lambda_L * F_R.p + (lambda_R * lambda_L * (U_R.p - U_L.p))) / (lambda_R - lambda_L);

			U_hll.d = (lambda_R * U_R.d - lambda_L * U_L.d + (F_L.d - F_R.d)) / (lambda_R - lambda_L);
			U_hll.v = (lambda_R * U_R.v - lambda_L * U_L.v + (F_L.v - F_R.v)) / (lambda_R - lambda_L);
			U_hll.p = (lambda_R * U_R.p- lambda_L * U_L.p + (F_L.p - F_R.p)) / (lambda_R - lambda_L);

			
#ifdef ONLY_RHLL
			F = F_hll;
#endif

			//============================================================================================//
#ifndef ONLY_RHLL		
//=========================ѕоиск скорости промежуточной волны===============================//
			const Type a = F_hll.p;			//F_E^hll
			const Type b = -U_hll.p - F_hll.v; // (E_hll + F_mx^hll)
			const Type c = U_hll.v;			//mx_hll

#if 1 // как описано в Mignone...
			Type quad = -0.5 * (b + SIGN(b) * sqrt(b * b - 4 * a * c));
			Type _lambda = c / quad;

#endif		
			{
				if (_lambda >= 0.0)
				{
					c2++;
					//============================ѕоиск промежуточного давлени€ ===================================//
					const Type _p = -F_hll.p * _lambda + F_hll.v;
					//============================================================================================//

					//==========================‘инальный поток HLLC=============================================//
					ConVal _U_L;
					const Type dif_L = 1.0 / (lambda_L - _lambda);

					_U_L.d = (U_L.d * (lambda_L - Vel_L)) * dif_L;
					_U_L.v = (U_L.v * (lambda_L - Vel_L) + _p - p_L) * dif_L;					
					_U_L.p = (U_L.p * (lambda_L - Vel_L) + _p * _lambda - p_L * Vel_L) * dif_L;

					F.d = F_L.d + lambda_L * (_U_L.d - U_L.d);
					F.v = F_L.v + lambda_L * (_U_L.v - U_L.v);
					F.p = F_L.p + lambda_L * (_U_L.p - U_L.p);

					//============================================================================================//
				}
				else //(_S <= 0)
				{
					c3++;
					//============================ѕоиск промежуточного давлени€ ===================================//
					const Type _p = -F_hll.p * _lambda + F_hll.v;
					//============================================================================================//
					ConVal _U_R;
					const Type dif_R = 1.0 / (lambda_R - _lambda);

					_U_R.d = (U_R.d * (lambda_R - Vel_R)) * dif_R;
					_U_R.v = (U_R.v * (lambda_R - Vel_R) + _p - p_R) * dif_R;					
					_U_R.p = (U_R.p * (lambda_R - Vel_R) + _p * _lambda - p_R * Vel_R) * dif_R;

					F.d = F_R.d + lambda_R * (_U_R.d - U_R.d);
					F.v = F_R.v + lambda_R * (_U_R.v - U_R.v);
					F.p = F_R.p + lambda_R * (_U_R.p - U_R.p);
				}
			}
#endif
		}


		
		F.v *= T;
		SumF.d += F.d;
		SumF.v += F.v;
		SumF.p += F.p;

		
		F_full[num_cell * 2 + i] = F;

	}// for

	//cout << "ready flew\n" << SumF << '\n';
	U.d -= (SumF.d * tau / h);
	U.v -= (SumF.v * tau / h);
	U.p -= (SumF.p * tau / h);

	return U;
}

int RHLLC_Init(std::vector<PhysVal>& W) {

	int n = W.size();
	//todo: check left 
	for (size_t i = 0; i < n; i++)
	{
		if (i * h < 0.5)
		{			
			W[i].rho = 1;
			W[i].v = 0.9;
			W[i].p = 1;
		}
		else
		{
			W[i].rho = 1;
			W[i].v = 0;
			W[i].p = 10;
		}
	}
	return 0;
}


int ReBuildCon2Phys(const std::vector<ConVal>& U, std::vector<PhysVal>& W)
{
	int n = U.size();
	for (size_t num_cell = 0; num_cell < n; num_cell++)
	{
		const Type vv = W[num_cell].v * W[num_cell].v;
		const Type d = W[num_cell].rho;
		Type Gamma0 = 1. / sqrt(1 - vv);
		const Type h = 1 + gamma_g * W[num_cell].p / d;

		Type W0 = d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;
		
		Type m = U[num_cell].v;
		Type mm = U[num_cell].v * U[num_cell].v;

		Type p = W[num_cell].p;
		Type v = W[num_cell].v;

		Type D = U[num_cell].d;
		Type E = U[num_cell].p;

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

			v = m / W0;
			
			err -= W0;
			cc++;
		} while (fabs(err / W0) > 1e-14);

		if (p < 0 || U[num_cell].d < 0 || std::isnan(p) || std::isnan(U[num_cell].d))
		{
			printf("Error cell %d (p = %lf, d= %lf)", num_cell, p, D / Gamma0);
			exit(1);
		}

		W[num_cell].p = p;
		W[num_cell].v = v;		
		W[num_cell].rho = D / Gamma0;
	}
	return 0;
}
int ReBuildPhys2Con(const std::vector<PhysVal>& W, std::vector<ConVal>& U)
{	
	ConVal cell;
	int n = W.size();

	for (size_t i = 0; i < n; i++)
	{
		const Type v = W[i].v * W[i].v;
		const Type d = W[i].rho;
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * W[i].p / d;
		const Type dhGG = d * h * Gamma * Gamma;

		cell.d = Gamma * d;
		cell.v = dhGG * W[i].v;		
		cell.p = dhGG - W[i].p;
		
		U[i] = cell;
	}

	return 0;
}

int RHLLC_1d(std::string& main_dir) {

	Type t = 0;
	const Type T = 0.5;

	RHLLC_Init(W_full_1d_prev);
	if (WriteSolution(main_dir, W_full_1d_prev, 1))
	{
		return 1;
	}

	int count = 0;
	while (t < T)
		//for (size_t k = 0; k < ; k++)		
	{
		ReBuildPhys2Con(W_full_1d_prev, U_full_1d_prev);

		for (size_t i = 0; i < size_g; i++)
		{
			U_full_1d[i] = RHLLC_flux(i);
		}
	
		ReBuildCon2Phys(U_full_1d, W_full_1d_prev);

		if (WriteSolution(main_dir, W_full_1d_prev, 0))
		{
			return 1;
		}

#ifdef DBG_OUTPUT
		ofstream ofile;
		static bool flag_clear = true;
		const std::string name_file = "D:\\Desktop\\FilesCourse\\dbg_rhllc_1d.txt";
		if (flag_clear)
		{
			ofile.open(name_file);
			ofile.close();
			flag_clear = false;
		}

		ofile.open(name_file, std::ios::app);
		if (!ofile.is_open())
		{
			printf("No open file\n");
			exit(1);
		}
		for (size_t i = 0; i < size_g; i++)
		{		
			ofile << "U[" << i << "]= " << U_full_1d[i].d << ", " << U_full_1d[i].v<< ", " << U_full_1d[i].p << '\n';
			ofile << "W[" << i << "]= " << W_full_1d_prev[i].rho << ", " << W_full_1d_prev[i].v << ", " << W_full_1d_prev[i].p << '\n';
			ofile << "F[" << i << "]= " << F_full[2*i].d + F_full[2*i+1].d << ", " 
				<< F_full[2 * i].v+ F_full[2 * i + 1].v << ", " << F_full[2 * i].p+ F_full[2 * i + 1].p << '\n';
		}
		ofile << "\n\n======================================\n\n";
		ofile.close();

#endif
		U_full_1d.swap(U_full_1d_prev);
		//W_full_1d.swap(W_full_1d_prev);

		t += tau;
		count++;
	}

	printf("timer %d\n", count);
	return 0;
}

int WriteSolution(const std::string& dir, const std::vector<PhysVal>& vector_U, bool flag_init)
{
	int n = vector_U.size();
	std::vector<Type> density(n);
	std::vector<Type> pressure(n);
	std::vector<Type> velocity(n);

	for (size_t i = 0; i < n; i++)
	{		
		density[i] = vector_U[i].rho;
		velocity[i] = vector_U[i].v;		
		pressure[i] = vector_U[i].p;
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