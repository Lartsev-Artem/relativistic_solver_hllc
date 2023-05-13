#if defined SOLVE
#include "../solve_config.h"
#if defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 && defined RHLLC_MPI
#include "../solve_global_struct.h"
#include "../../global_value.h"
#include "../../utils/grid_geometry/geometry_solve.h"
#include "../solve_utils.h"
#include "rhllc_utils.h"

extern std::vector<flux_t> phys_local;
extern std::vector<int> send_hllc;
extern std::vector<int> disp_hllc;


static int rhllc_get_conv_value_ost1098(const flux_t& W, flux_t& U)
{
	const Type d = W.d;
	const Type Gamma = 1. / sqrt(1 - W.v.dot(W.v));
	const Type h = 1 + gamma_g * W.p / d;
	const Type dhGG = d * h * Gamma * Gamma;

	U.d = Gamma * W.d;
	U.v = dhGG * W.v;
	U.p = dhGG - W.p;

	return 0;
}

static int flux_t_calc(const flux_t& conv_val_l, const flux_t& conv_val_r,
	const flux_t& phys_val_l, const flux_t& phys_val_r,
	face_t& f)
{
	/*
An HLLC Riemann solver for relativistic flows Ц I. Hydrodynamics
A. Mignone and G. Bodo
INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
Mon. Not. R. Astron. Soc. 364, 126Ц136 (2005)

https://github.com/PrincetonUniversity/Athena-Cversion/blob/master/src/rsolvers/hllc_sr.c

\note: комбинаци€ кода из mignone 2005 и 2006. в части hllc
*/

	Matrix3 T;
	MakeRotationMatrix(f.geo.n, T);

	flux_t U_L = conv_val_l;
	U_L.v = T * conv_val_l.v;

	flux_t U_R = conv_val_r;
	U_R.v = T * conv_val_r.v;

	flux_t W_L = phys_val_l;
	W_L.v = T * phys_val_l.v;

	flux_t W_R = phys_val_r;
	W_R.v = T * phys_val_r.v;

	//====================  эшируем физические переменные слева и справа============================//
	// нормальна€ сокорость
	const Vector3 Vel_L(W_L[1], W_L[2], W_L[3]);  //T * velocity[num_cell];
	const Vector3 Vel_R(W_R[1], W_R[2], W_R[3]);  //T * velocity[neig / 4];

	const Type d_L = W_L(0);
	const Type d_R = W_R(0);

	const Type p_L = W_L(4);
	const Type p_R = W_R(4);

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

	flux_t F;
	if (lambda_R <= 0) // если верно выполнить всегда
	{
		F.d = U_R[0] * Vel_R[0]; //D*v_x
		F.v(0) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
		F.v(1) = U_R[2] * Vel_R[0];
		F.v(2) = U_R[3] * Vel_R[0];
		F.p = U_R[1];
		//continue;			
	}
	else if (lambda_L >= 0) // выполнить либо по условию либо дл€ всех границ
	{
		F.d = U_L[0] * Vel_L[0]; //D*v_x
		F.v(0) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
		F.v(1) = U_L[2] * Vel_L[0];
		F.v(2) = U_L[3] * Vel_L[0];
		F.p = U_L[1];
		//continue;			
	}
	else
	{
		//====================–асчЄт потоков и приближений hll=========================================//
		flux_t F_L;
		flux_t F_R;
		flux_t U_hll;
		flux_t F_hll;

		F_R.d = U_R[0] * Vel_R[0]; //D*v_x
		F_R.v(0) = U_R[1] * Vel_R[0] + p_R; //mx*vx+p
		F_R.v(1) = U_R[2] * Vel_R[0];
		F_R.v(2) = U_R[3] * Vel_R[0];
		F_R.p = U_R[1];

		F_L.d = U_L[0] * Vel_L[0]; //D*v_x
		F_L.v(0) = U_L[1] * Vel_L[0] + p_L; //mx*vx+p
		F_L.v(1) = U_L[2] * Vel_L[0];
		F_L.v(2) = U_L[3] * Vel_L[0];
		F_L.p = U_L[1];

		for (int i = 0; i < base + 1; i++)
		{
			F_hll[i] = (F_L[i] * lambda_R - F_R[i] * lambda_L + ((U_R[i] - U_L[i]) * lambda_R * lambda_L)) / (lambda_R - lambda_L);
			U_hll[i] = ((U_R[i] * lambda_R) - (U_L[i] * lambda_L) + (F_L[i] - F_R[i])) / (lambda_R - lambda_L);
		}
		/*F_hll = (F_L * lambda_R - F_R * lambda_L + ((U_R - U_L) * lambda_R * lambda_L)) / (lambda_R - lambda_L);
		U_hll = ((U_R * lambda_R) - (U_L * lambda_L) + (F_L - F_R)) / (lambda_R - lambda_L);*/

#ifdef ONLY_RHLL
		F = F_hll;
#endif

		//============================================================================================//
#ifndef ONLY_RHLL		
//=========================ѕоиск скорости промежуточной волны===============================//
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
				//============================ѕоиск промежуточного давлени€ ===================================//
				const Type _p = -F_hll[4] * _lambda + F_hll[1];
				//============================================================================================//

				//==========================‘инальный поток HLLC=============================================//
				flux_t _U_L;
				const Type dif_L = 1.0 / (lambda_L - _lambda);

				_U_L.d = (U_L[0] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.v[0] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
				_U_L.v[1] = (U_L[2] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.v[2] = (U_L[3] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.p = (U_L[4] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

				F = F_L + (_U_L - U_L) * lambda_L;

				//============================================================================================//
			}
			else //(_S <= 0)
			{
				//============================ѕоиск промежуточного давлени€ ===================================//
				const Type _p = -F_hll[4] * _lambda + F_hll[1];
				//============================================================================================//
				flux_t _U_R;
				const Type dif_R = 1.0 / (lambda_R - _lambda);

				_U_R.d = (U_R[0] * (lambda_R - Vel_R[0])) * dif_R;
				_U_R.v[0] = (U_R[1] * (lambda_R - Vel_R[0]) + _p - p_R) * dif_R;
				_U_R.v[1] = (U_R[2] * (lambda_R - Vel_R[0])) * dif_R;
				_U_R.v[2] = (U_R[3] * (lambda_R - Vel_R[0])) * dif_R;
				_U_R.p = (U_R[4] * (lambda_R - Vel_R[0]) + _p * _lambda - p_R * Vel_R[0]) * dif_R;

				F = F_R + (_U_R - U_R) * lambda_R;
			}
		}
#endif
	}
	//WRITE_LOG(" F0= " << F[0] << ' ' << F[1] << ' ' << F[2] << ' ' << F[3] << ' ' << F[4] << '\n');
	f.f = F;
	f.f.v = (T.transpose()) * F.v;
	//WRITE_LOG(" F1= " << f.f[0] << ' ' << f.f[1] << ' ' << f.f[2] << ' ' << f.f[3] << ' ' << f.f[4] << '\n');
	f.f = f.f * f.geo.S;

	return 0;
}


static int rhllc_get_phys_value_ost1098(const flux_t& U, flux_t& W)
{
	int call_back_flag = 0;
	Type Gamma0 = 1. / sqrt(1 - (W.v.dot(W.v)));
	Type p = W.p;

	const Type h = 1 + gamma_g * p / W.d;

	Type W0 = W.d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m = U.v;
	Type mm = m.dot(m);


	Vector3 v = W.v;

	Type D = U.d;
	Type E = U.p;

	int  cc = 0;

	Type err = 1;
	do
	{
		err = W0;

		//W.d * h * Gamma0 * Gamma0 - p - E;
		Type fW = W0 - p - E;

		Type dGdW = -(Gamma0 * Gamma0 * Gamma0) * mm / (2 * W0 * W0 * W0);
		Type dFdW = 1 - ((Gamma0 * (1 + D * dGdW) - 2 * W0 * dGdW) / (Gamma0 * Gamma0 * Gamma0 * gamma_g));
		W0 -= (fW / dFdW);

		Gamma0 = 1. / sqrt(1 - mm / (W0 * W0));

		p = (W0 - D * Gamma0) / (Gamma0 * Gamma0 * gamma_g);

		v = m / W0;

		err -= W0;
		cc++;

		if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d))
		{
			//	EXIT(1);
			//	p = max(sqrt(mm) - E, 1e-20);
			//	if (std::isnan(p)) p = 1e-20;
			//	v = m / (E + p);
			//	if (v.norm() > 0.9999999995)
			//	{
			//		v /= 1.0005;
			//	}
			//	Gamma0 = 1. / sqrt(1 - v.dot(v));
			////	printf("try resolve %.16lf\n", v.norm());
			//	//printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
			//	//D_LD;
			call_back_flag = 1;
			break;
		}

	} while (fabs(err / W0) > 1e-14);

	/*if (p < 0 || U.d < 0 || std::isnan(p) || std::isnan(U.d) || v.norm()>1)
	{
		printf("W= %lf,(%lf), %lf\n", W.d, W.v.norm(), W.p);
		printf("Error cell (p = %lf, d= %lf, %lf)\n", p, D , Gamma0);
		EXIT(1);
		call_back_flag = 2;
	}*/

	W.d = D / Gamma0;
	if (std::isnan(Gamma0) || std::isinf(Gamma0) || (D / Gamma0) < 1e-10)
	{
		W.d = 1e-10;
		call_back_flag = 1;
	}
	W.v = v;

	W.p = p;
	if (std::isnan(p) || std::isinf(p) || p < 1e-20)
	{
		W.p = 1e-20;
		call_back_flag = 1;
	}

	return call_back_flag;
}


static std::vector<flux_all_t> left_neig;
static std::vector<flux_all_t> right_neig;
std::vector<mpi_conf_t> mpi_conf;

static std::vector<MPI_Request> send_rq_rhllc;
static std::vector<MPI_Request> recv_rq_rhllc;

int InitMPI_RHllc(const std::vector<elem_t>& cells)
{
	int np, myid;
	MPI_GET_INF(np, myid);
	
	if (myid == 0 || myid == np - 1)
	{
		send_rq_rhllc.resize(1, MPI_REQUEST_NULL);
		recv_rq_rhllc.resize(1, MPI_REQUEST_NULL);
	}
	else
	{
		send_rq_rhllc.resize(2, MPI_REQUEST_NULL);
		recv_rq_rhllc.resize(2, MPI_REQUEST_NULL);
	}

	int cc = 0;
	if (myid != np - 1)
	{
		int size_l = mpi_conf[myid + 1].left_cell_max - mpi_conf[myid + 1].left_cell_min + 1;
		MPI_Send_init(cells.data() + mpi_conf[myid + 1].left_cell_min, size_l, MPI_flux_elem_t, myid + 1, 0, MPI_COMM_WORLD, &send_rq_rhllc[cc]);
		MPI_Recv_init(right_neig.data(), right_neig.size(), MPI_flux_all_t, myid + 1, 0, MPI_COMM_WORLD, &recv_rq_rhllc[cc++]);
	}

	if (myid != 0)
	{
		int size_r = mpi_conf[myid - 1].right_cell_max - mpi_conf[myid - 1].right_cell_min + 1;
		MPI_Send_init(cells.data() + mpi_conf[myid - 1].right_cell_min, size_r, MPI_flux_elem_t, myid - 1, 0, MPI_COMM_WORLD, &send_rq_rhllc[cc]);
		MPI_Recv_init(left_neig.data(), left_neig.size(), MPI_flux_all_t, myid - 1, 0, MPI_COMM_WORLD, &recv_rq_rhllc[cc++]);
	}

	return 0;
}

int ReadMpiConf(const std::string& file, int myid, int np)
{
	std::ifstream ifile;
	OPEN_FSTREAM(ifile, file.c_str());

	int size;
	ifile >> size;

	if (np != size)
	{
		printf("Bad size conf file. Need %d nodes\n", size);
		D_LD;
	}

	mpi_conf.resize(np);
	int id = 0;

	for (int i = 0; i <= np; i++)
	{
		ifile >> id;

		ifile >> mpi_conf[id].left;
		ifile >> mpi_conf[id].right;

		ifile >> mpi_conf[id].left_cell_min;
		ifile >> mpi_conf[id].left_cell_max;

		ifile >> mpi_conf[id].right_cell_min;
		ifile >> mpi_conf[id].right_cell_max;
	}
	ifile.close();

	if(mpi_conf[myid].left_cell_max != -1)
	left_neig.resize(mpi_conf[myid].left_cell_max - mpi_conf[myid].left_cell_min + 1);

	if (mpi_conf[myid].right_cell_max != -1)
	right_neig.resize(mpi_conf[myid].right_cell_max - mpi_conf[myid].right_cell_min + 1);
}

int RHLLC_3d_MPI(const Type tau, grid_t& grid)
{
	const int size_grid = grid.size;

	int np, myid;
	MPI_GET_INF(np, myid);

#ifdef ILLUM

#pragma omp parallel default(none) shared(tau,myid, grid, glb_files, mpi_conf)
	{
		const int left = mpi_conf[myid].left;
		const int right = mpi_conf[myid].right;
		// востановление физических переменных
#pragma omp for
		for (int i = left; i < right; i++)
		{
			int back = rhllc_get_phys_value_ost1098(grid.cells[i].conv_val, grid.cells[i].phys_val);
			if (back)
			{
				if (back == 2)
				{
					D_LD;
				}
				else
				{
					//printf("try id= %d\n", i);
				}
				// если был пернсчЄт
				rhllc_get_conv_value_ost1098(grid.cells[i].phys_val, grid.cells[i].conv_val);
			}
		}
	}
#endif

	MPI_Startall(recv_rq_rhllc.size(), recv_rq_rhllc.data());

	MPI_Waitall(send_rq_rhllc.size(), send_rq_rhllc.data(), MPI_STATUSES_IGNORE);
	MPI_Startall(send_rq_rhllc.size(), send_rq_rhllc.data());
	

#pragma omp parallel default(none) shared(tau, grid, glb_files, mpi_conf, myid)
	{
		flux_t bound_val;
		flux_t phys_bound_val;
		Matrix3 T;
		Matrix3 TT;
		elem_t* cell;
		const int size_face = grid.faces.size();

#pragma omp for
		for (int i = 0; i < size_face; i++)
		{
			face_t& f = grid.faces[i];

			if (f.geo.id_l < mpi_conf[myid].left || f.geo.id_l >= mpi_conf[myid].right)
			{
				continue;
			}

			cell = &grid.cells[f.geo.id_l];
			switch (f.geo.id_r)// id соседа она же признак √”
			{
			case eBound_FreeBound:
				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
				break;
			case eBound_InnerSource:
#ifdef Sphere
				phys_bound_val.d = 0.1; phys_bound_val.v << 0, 0, 0; phys_bound_val.p = 0.1;
				rhllc_get_conv_value_ost1098(phys_bound_val, bound_val);
#elif defined Cone_JET
				//phys_bound_val.d = 0.1; // (3 * 1e-8 + 1e-12) / DENSITY;
				//phys_bound_val.p = 1; // (100 + (1e-2)) / PRESSURE;
				//phys_bound_val.v = Vector3(1e-4, 0, 0);// (Vector3(1e4, 0, 0)) / VELOCITY;

				//phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				//phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
				//phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;

				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
				rhllc_get_conv_value_ost1098(phys_bound_val, bound_val);
#else
				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
#endif
				break;
			case eBound_OutSource:
#if defined Cylinder 
				phys_bound_val.d = 0.1;
				phys_bound_val.v << 0.99, 0, 0;
				phys_bound_val.p = 0.01;
				rhllc_get_conv_value_ost1098(phys_bound_val, bound_val);
#elif defined Cone
				//phys_bound_val.d = 0.1; // (3 * 1e-8 + 1e-12) / DENSITY;
				//phys_bound_val.p = 1; // (100 + (1e-2)) / PRESSURE;
				//phys_bound_val.v = Vector3(1e-4, 0, 0);// (Vector3(1e4, 0, 0)) / VELOCITY;

				//phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				//phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
				//phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;

				phys_bound_val.d = 0.1;
				phys_bound_val.v << 0.99, 0, 0;
				phys_bound_val.p = 0.01;
				rhllc_get_conv_value_ost1098(phys_bound_val, bound_val);
#else
				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
#endif	
				break;
			case eBound_LockBound:
#if 1
				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
				MakeRotationMatrix(f.geo.n, T);

				bound_val.v = T * bound_val.v;
				phys_bound_val.v = T * phys_bound_val.v;

				bound_val.v[0] = -bound_val.v[0];
				phys_bound_val.v[0] = -phys_bound_val.v[0];

				TT = T.transpose();

				bound_val.v = TT * bound_val.v;
				phys_bound_val.v = TT * phys_bound_val.v;
#else
				bound_val = cell->conv_val;
				phys_bound_val = cell->phys_val;
#endif
				break;

			default:
				DIE_IF(f.geo.id_r < 0); //Err bound in RHLLC_3d

				if (f.geo.id_r < mpi_conf[myid].left || f.geo.id_r >= mpi_conf[myid].right)
				{
					continue;
				}

				bound_val = grid.cells[f.geo.id_r].conv_val;
				phys_bound_val = grid.cells[f.geo.id_r].phys_val;
				break;
			}

			flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val, grid.cells[f.geo.id_l].phys_val, phys_bound_val, f);
		}

		//#pragma omp barrier
	}

	MPI_Waitall(recv_rq_rhllc.size(), recv_rq_rhllc.data(), MPI_STATUSES_IGNORE);
	int flags[2];
	MPI_Testall(send_rq_rhllc.size(), send_rq_rhllc.data(), flags, MPI_STATUSES_IGNORE);

#pragma omp parallel default(none) shared(tau, grid, glb_files, mpi_conf, myid, right_neig, left_neig)
	{
		// after mpi
		flux_t bound_val;
		flux_t phys_bound_val;
		const int size_face = grid.faces.size();
		int r_id = 0;
		int l_id = 0;

		if (right_neig.size() != 0)//(myid == 0)
		{
#pragma omp for
			for (int i = 0; i < size_face; i++)
			{
				face_t& f = grid.faces[i];

				if (f.geo.id_r < 0) continue;

				if (f.geo.id_l < mpi_conf[myid].left || f.geo.id_l >= mpi_conf[myid].right)//текущий узел
				{
					continue;
				}

				if (f.geo.id_r >= mpi_conf[myid].right_cell_min && f.geo.id_r <= mpi_conf[myid].right_cell_max)
				{
					r_id = f.geo.id_r - mpi_conf[myid].right_cell_min;
					bound_val = right_neig[r_id].conv_val;
					phys_bound_val = right_neig[r_id].phys_val;
				}
				else
				{
					continue;
				}

				flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val, grid.cells[f.geo.id_l].phys_val, phys_bound_val, f);
			}
		}

		if (left_neig.size() != 0)
		{

#pragma omp for
			for (int i = 0; i < size_face; i++)
			{
				face_t& f = grid.faces[i];

				if (f.geo.id_r < 0) continue;

				if (f.geo.id_r < mpi_conf[myid].left || f.geo.id_r >= mpi_conf[myid].right)//текущий узел
				{
					continue;
				}

				if (f.geo.id_l >= mpi_conf[myid].left_cell_min && f.geo.id_l <= mpi_conf[myid].left_cell_max)
				{
					l_id = f.geo.id_l - mpi_conf[myid].left_cell_min;
					bound_val = left_neig[l_id].conv_val;
					phys_bound_val = left_neig[l_id].phys_val;
					//l_id++;
				}
				else
				{
					continue;
				}

				flux_t_calc(bound_val, grid.cells[f.geo.id_r].conv_val, phys_bound_val, grid.cells[f.geo.id_r].phys_val, f);
			}
		}
	}


#ifdef ILLUM
	auto calc{ [&grid, myid, tau](const int left, const int right)
	{
#pragma omp for
		for (int i = left; i < right; i++)
		{
			elem_t& el = grid.cells[i];
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

			rhllc_get_phys_value_ost1098(el.conv_val, el.phys_val); // востановление физических переменных
			phys_local[i] = el.phys_val;
		}
	}};
	
	const int left = mpi_conf[myid].left;

	for (int i = 0; i < disp_hllc.size(); i++)
	{
#pragma omp parallel default(none)firstprivate(i,tau, left) shared(grid, send_hllc, disp_hllc, phys_local, calc)
		{
			calc(left + disp_hllc[i], left + disp_hllc[i] + send_hllc[i]);
		}

//#pragma omp single
		{
			SendPhysValue(phys_local.data() + left + disp_hllc[i], send_hllc[i], i);
		}
	}

#else

#pragma omp parallel default(none) shared(tau, grid, glb_files, mpi_conf, myid)
	{
		const int size_grid = grid.size;

#pragma omp for
		for (int i = mpi_conf[myid].left; i < mpi_conf[myid].right; i++)//for (int i = 0; i < size_grid; i++)
		{
			elem_t& el = grid.cells[i];
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

			rhllc_get_phys_value_ost1098(el.conv_val, el.phys_val); // востановление физических переменных
		}
	} //omp

#endif

	return 0;
}

#endif 
#endif // SOLVE && USE_MPI