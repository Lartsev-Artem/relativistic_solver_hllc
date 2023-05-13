#include "../solve_config.h"
#if defined RHLLC && NUMBER_OF_MEASUREMENTS == 3 && !defined RHLLC_MPI && defined SOLVE
#include "../solve_global_struct.h"
#include "../../file_module/reader_bin.h"
#include "../../utils/grid_geometry/geometry_solve.h"
#include "../solve_utils.h"
#include "rhllc_utils.h"

#ifdef USE_MPI
extern std::vector<flux_t> phys_local;
extern std::vector<int> send_hllc;
extern std::vector<int> disp_hllc;
#endif

#ifdef OLD_CLASS
static int ReBuildDataForRHLLC(const std::vector<VectorX>& W_full, std::vector<VectorX>& data)
{
	const int N = W_full.size();
	data.resize(N);
	VectorX cell(5);

	for (size_t i = 0; i < N; i++)
	{
		const Vector3 velocity = Vector3(W_full[i](1), W_full[i](2), W_full[i](3));
		const Type v = velocity.dot(velocity);
		const Type d = W_full[i](0);
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * W_full[i](4) / d;
		const Type dhGG = d * h * Gamma * Gamma;

		cell[0] = Gamma * d;
		cell[1] = dhGG * velocity[0];
		cell[2] = dhGG * velocity[1];
		cell[3] = dhGG * velocity[2];

		cell[4] = dhGG - W_full[i](4);
		//cell[4] = pressure[i] / (gamma1 - 1) + d * v / 2;
		//cell[4] = e_substance[i];

		data[i] = cell;		
	}

	return 0;
}

static VectorX RHLLC_stepToOMPGit(const int num_cell, const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, const std::vector<VectorX>& U_full_prev,
	const std::vector<VectorX>& W_full) {

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
#ifdef Cylinder				
				W_R << 0.1, 0.99, 0, 0, 0.01;

				U_R << 0.7088812050083355,
					6.218592964824112,
					0,
					0,
					6.271407035175871;
#else
				U_R = U;
				W_R = W;
#endif
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

		WRITE_LOG(i<<" F0= "<< F[0] << ' ' << F[1] << ' ' << F[2] << ' ' << F[3] << ' ' << F[4] << '\n');
		VectorX buf = F;
		F = (T.transpose()) * buf;

		WRITE_LOG(i << " F1= " << F[0] << ' ' << F[1] << ' ' << F[2] << ' ' << F[3] << ' ' << F[4] << '\n');
		SumF += F * squares_cell[4 * num_cell + i];
	

	}// for
	WRITE_LOG(num_cell<<" SumF= " << SumF[0] << ' ' << SumF[1] << ' ' << SumF[2] << ' ' << SumF[3] << ' ' << SumF[4] << '\n');
	return (U - SumF * tau / volume[num_cell]);
}


static int ReBuildPhysicValue_ost1098(const VectorX& U, VectorX& W) {

	const Vector3 velocity_ = Vector3(W(1), W(2), W(3));
	const Type vv = velocity_.dot(velocity_);
	const Type d = W(0);
	Type Gamma0 = 1. / sqrt(1 - vv);
	const Type h = 1 + gamma_g * W(4) / d;

	Type W0 = d * h * Gamma0 * Gamma0; //U[0] * Gamma0 * h;

	Vector3 m(U[1], U[2], U[3]);
	Type mm = m.dot(m);

	Type p =W(4);
	Vector3 v = velocity_;

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
		printf("Error cell (p = %lf, d= %lf)", p, D / Gamma0);
		exit(1);
	}

	W(4) = p;

	W(1) = v(0);
	W(2) = v(1);
	W(3) = v(2);
	//Gamma0;
	W(0) = D / Gamma0;

	return 0;
}

int RHLLC_Init_3d(const int N, const std::vector<Vector3>& centerts, std::vector<VectorX>& W) {

	W.resize(N);
	VectorX cell(5);

	for (size_t i = 0; i < N; i++)
	{
		Vector3 x = centerts[i];

#ifdef Jet_3d		
		if (Vector2(x[1], x[2]).norm() < 0.2 && x[0] < 0.5)
		{
			cell(0) = 0.1;
			cell(4) = 0.01;

			cell(1) = 0.99;
			cell(2) = 0;
			cell(3) = 0;			
		}
		else
		{
			cell(0) = 10;
			cell(4) = 0.01;

			cell(1) = 0;
			cell(2) = 0;
			cell(3) = 0;			
		}
#endif

#if 1 //def SODA
		if (x[0] < 0.5)
		{
			cell(0) = 1;
			cell(1) = 0.9;
			cell(2) = 0;
			cell(3) = 0;
			cell(4) = 1;
		}
		else
		{
			cell(0) = 1;
			cell(1) = 0;
			cell(2) = 0;
			cell(3) = 0;
			cell(4) = 10;
		}
#endif

		W[i] = cell;
	}

	return 0;
}

void RHLLC_3d_old(const Type tau, const std::vector<int>& neighbours_id_faces, const std::vector<Normals>& normals,
	const std::vector<Type>& squares_cell, const std::vector<Type>& volume, std::vector<VectorX>& U_full_prev, std::vector<VectorX>& U_full,
	std::vector<VectorX> W_full)
{
	const int size_grid = U_full.size();

	ReBuildDataForRHLLC(W_full, U_full_prev); // p0, v0, rho0 -> U_full_prev
	omp_set_num_threads(1);
#pragma omp parallel default(none) shared(size_grid, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, U_full)
	{
		VectorX buf(5);
#pragma omp for
		for (int i = 0; i < size_grid; i++)
		{
			//buf = HLLC_stepToOMPRel(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev); // -> U_full		
			U_full[i] = RHLLC_stepToOMPGit(i, tau, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, W_full); // -> U_full

			WRITE_LOG(i << " U_full= " << U_full[i][0] << ' ' << U_full[i][1] << ' ' << U_full[i][2] << ' ' << U_full[i][3] << ' ' << U_full[i][4] << "\n\n");
#if 0//def Cube
			// Для 1d задачи Сода. 
			// СЛИШКОМ СИЛЬНЫЕ КОЛЕБАНИЯ ПО Y,Z. что-то не так!!!!!!!!!!!!!!!!!!!!
			//buf[1] = Vector3(buf[1], buf[2], buf[3]).norm();
			buf[2] = 0;
			buf[3] = 0;
#endif
		}
	}

#pragma omp parallel default(none) shared(size_grid,  gamma_g, U_full)
	{
#pragma omp for
		for (int i = 0; i < size_grid; i++)
		{
			ReBuildPhysicValue_ost1098(U_full[i], W_full[i]); // U_full -> p1, v1, rho1			
		}
	}

	U_full.swap(U_full_prev);
}

#include "../../file_module/writer_bin.h"
int RHLLC_3d(const Type tau, grid_t& grid)
{
	std::vector<int> neighbours_id_faces;
	std::vector<Normals> normals;
	std::vector<Type> squares_cell;
	std::vector<Type>volume; 
	std::vector<Vector3> centers;

	std::vector<VectorX> U_full_prev;
	std::vector<VectorX> U_full;
	std::vector<VectorX> W_full;

	const std::string name_file_id_neighbors = glb_files.base_adress + "pairs.bin";
	const std::string name_file_normals = glb_files.base_adress + "normals.bin";
	const std::string name_file_centers = glb_files.base_adress + "centers.bin";
	const std::string name_file_squares = glb_files.base_adress + "squares.bin";
	const std::string name_file_volume = glb_files.base_adress + "volume.bin";

	if (ReadSimpleFileBin(name_file_id_neighbors, neighbours_id_faces)) RETURN_ERR("Error reading file neighbours\n");
	if (ReadNormalFile(name_file_normals, normals)) RETURN_ERR("Error reading file normals\n");
	if (ReadSimpleFileBin(name_file_squares, squares_cell)) RETURN_ERR("Error reading file squares_faces\n");
	if (ReadSimpleFileBin(name_file_volume, volume)) RETURN_ERR("Error reading file volume\n");
	if (ReadSimpleFileBin(name_file_centers, centers)) RETURN_ERR("Error reading file centers\n");

	RHLLC_Init_3d(centers.size(), centers, W_full);
	ReBuildDataForRHLLC(W_full, U_full_prev);
	U_full.assign(U_full_prev.begin(), U_full_prev.end());

	RHLLC_3d_old(1e-5, neighbours_id_faces, normals, squares_cell, volume, U_full_prev, U_full, W_full);

	//for (size_t i = 0; i < grid.size; i++)	
	//{
	//	for (size_t j = 0; j < base+1; j++)
	//	{	
	//		grid.cells[i].conv_val(j) = U_full[i][j];
	//		grid.cells[i].phys_val[j] = W_full[i][j];
	//	}
	//	
	//}
	
	return 0;
}

#else //OLD_CLASS

int rhllc_get_conv_value_ost1098(const flux_t& W, flux_t& U)
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
An HLLC Riemann solver for relativistic flows – I. Hydrodynamics
A. Mignone and G. Bodo
INAF Osservatorio Astronomico di Torino, 10025 Pino Torinese, Italy
Accepted 2005 August 22. Received 2005 August 16; in original form 2005 May 16
Mon. Not. R. Astron. Soc. 364, 126–136 (2005)

https://github.com/PrincetonUniversity/Athena-Cversion/blob/master/src/rsolvers/hllc_sr.c

\note: комбинация кода из mignone 2005 и 2006. в части hllc
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
	else if (lambda_L >= 0) // выполнить либо по условию либо для всех границ
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
		//====================Расчёт потоков и приближений hll=========================================//
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

		for (int i = 0; i < base+1; i++)
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
				//============================Поиск промежуточного давления ===================================//
				const Type _p = -F_hll[4] * _lambda + F_hll[1];
				//============================================================================================//

				//==========================Финальный поток HLLC=============================================//
				flux_t _U_L;
				const Type dif_L = 1.0 / (lambda_L - _lambda);

				_U_L.d = (U_L[0] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.v[0] = (U_L[1] * (lambda_L - Vel_L[0]) + _p - p_L) * dif_L;
				_U_L.v[1] = (U_L[2] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.v[2] = (U_L[3] * (lambda_L - Vel_L[0])) * dif_L;
				_U_L.p = (U_L[4] * (lambda_L - Vel_L[0]) + _p * _lambda - p_L * Vel_L[0]) * dif_L;

				F = F_L + (_U_L - U_L)* lambda_L;

				//============================================================================================//
			}
			else //(_S <= 0)
			{
				//============================Поиск промежуточного давления ===================================//
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

int RHLLC_3d(const Type tau, grid_t& grid)
{
#pragma omp parallel default(none) shared(tau, grid, glb_files)
	{
		const int size_grid = grid.size;	

#ifdef ILLUM
		// востановление физических переменных
#pragma omp for
		for (int i = 0; i < size_grid; i++)
		{
			int back = rhllc_get_phys_value_ost1098(grid.cells[i].conv_val, grid.cells[i].phys_val);
			if (back)
			{		
				D_LD;
				if (back == 2)
				{
					
				}
				else
				{
					//printf("try id= %d\n", i);
				}
				// если был пернсчёт
				rhllc_get_conv_value_ost1098(grid.cells[i].phys_val, grid.cells[i].conv_val);
			}			
		}		
#endif
	}

#pragma omp parallel default(none) shared(tau, grid, glb_files)
	{
		flux_t bound_val;
		flux_t phys_bound_val;
		Matrix3 T;
		Matrix3 TT;
		elem_t* cell;		
		const int size_face = grid.faces.size();

		// потоки
#pragma omp for
		for (int i = 0; i < size_face; i++)
		{			
			face_t& f = grid.faces[i];
			cell = &grid.cells[f.geo.id_l];
			switch (f.geo.id_r)// id соседа она же признак ГУ
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
				//phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				//phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
				//phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;
				phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
				phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;
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
			/*	phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
				phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;*/

				//phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
				//phys_bound_val.p = /*Pressure(Vector3::Zero())*/5000. / PRESSURE;
				//phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;

				phys_bound_val.d = 1;
				phys_bound_val.p = 1;
				phys_bound_val.v = Vector3(1e-2, 0, 0);

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
				
				bound_val = grid.cells[f.geo.id_r].conv_val;
				phys_bound_val = grid.cells[f.geo.id_r].phys_val;
				break;
			}

			flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val,
				grid.cells[f.geo.id_l].phys_val, phys_bound_val, f);
		}

		//#pragma omp barrier 
	}

#pragma omp parallel default(none) shared(tau, grid)
	{
		const int size_grid = grid.size;		
#pragma omp for
		for (int i = 0; i < size_grid; i++)
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

	return 0;

}


int MPI_RHLLC_3d(const int myid, const Type tau, grid_t& grid)
{
	//todo: shadule (static 1) (dynamic 1)
#ifdef ILLUM

	if (myid == 0)
	{
#pragma omp parallel default(none) shared(tau, grid, glb_files)
		{
			const int size_grid = grid.size;
#ifdef ILLUM
			// востановление физических переменных
#pragma omp for
			for (int i = 0; i < size_grid; i++)
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
					// если был пернсчёт
					rhllc_get_conv_value_ost1098(grid.cells[i].phys_val, grid.cells[i].conv_val);
				}
			}
#endif
		}

#pragma omp parallel default(none) shared(tau, grid, glb_files)
		{
			flux_t bound_val;
			flux_t phys_bound_val;
			Matrix3 T;
			Matrix3 TT;
			elem_t* cell;
			const int size_face = grid.faces.size();

			// потоки
#pragma omp for
			for (int i = 0; i < size_face; i++)
			{
				face_t& f = grid.faces[i];
				cell = &grid.cells[f.geo.id_l];
				switch (f.geo.id_r)// id соседа она же признак ГУ
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
					
					phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
					phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
					phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;
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
					
					phys_bound_val.d = Density(Vector3::Zero()) / DENSITY;
					phys_bound_val.p = Pressure(Vector3::Zero()) / PRESSURE;
					phys_bound_val.v = Velocity(Vector3::Zero()) / VELOCITY;
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

					bound_val = grid.cells[f.geo.id_r].conv_val;
					phys_bound_val = grid.cells[f.geo.id_r].phys_val;
					break;
				}

				flux_t_calc(grid.cells[f.geo.id_l].conv_val, bound_val,
					grid.cells[f.geo.id_l].phys_val, phys_bound_val, f);
			}

			//#pragma omp barrier 
		}


	} //myid ==0
	

		auto calc{ [&grid, myid, tau](const int left, const int right)
		{
			if (myid == 0)
			{
//#pragma omp parallel default(none) shared(tau, grid, send_hllc, disp_hllc, phys_local,left,right)
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
			}
		}
		} };

		const int size_grid = grid.size;

		for (int i = 0; i < disp_hllc.size(); i++)
		{
#pragma omp parallel default(none)firstprivate(i,tau) shared(grid, send_hllc, disp_hllc, phys_local, calc)
			{
				calc(disp_hllc[i], disp_hllc[i] + send_hllc[i]);
			}

//#pragma omp single
			{
				SendPhysValue(phys_local.data()+ disp_hllc[i], send_hllc[i], i);
			}
		}
	

#else
#pragma error "bad prj config"
#endif //USE_MPI


	return 0;

}


#endif //NEW_CLASS
#endif // RHLLC_3d
