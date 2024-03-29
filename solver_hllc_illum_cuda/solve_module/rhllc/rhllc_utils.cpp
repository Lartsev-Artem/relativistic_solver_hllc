#include "rhllc_utils.h"

#ifdef RHLLC
#include "../../file_module/reader_bin.h"
#include "../../file_module/reader_txt.h"

int RHllcGetTimeStep(hllc_value_t& hllc_set, const std::vector<elem_t>& cells) 
{

	hllc_set.tau = hllc_set.CFL * hllc_set.h;

	hllc_set.print_timer = hllc_set.tau*0.1;
	hllc_set.T = hllc_set.tau * 3;

	if (hllc_set.tau < 0)
	{
		EXIT_ERRS("Error tau = %lf\n", hllc_set.tau);
	}

	return  0;
}

int RHllcPhysToConv(std::vector<elem_t>& cells)
{
	for (auto& el : cells)
	{
		const Type v = el.phys_val.v.dot(el.phys_val.v);
		const Type d = el.phys_val.d;
		const Type Gamma = 1. / sqrt(1 - v);
		const Type h = 1 + gamma_g * el.phys_val.p / d;
		const Type dhGG = d * h * Gamma * Gamma;

		el.conv_val.d = Gamma * d;
		el.conv_val.v = el.phys_val.v*dhGG;
		el.conv_val.p = dhGG - el.phys_val.p;
	}
	return 0;
}

Type Density(const Vector3& p)
{
	return 1e-8;
	const Type x = p[0];	
	const Type x1 = 0.04;
	const Type x2 = 0.06;
	const Type x3 = 0.1;

	const Type a = 1e-9;
	const Type b = 1e-14;
	const Type betta0 = 0.012;
	const Type betta1 = 0.1;
	if (x < x1)
	{
		const Type val = (x1 - x2) / betta0;
		return a * exp(-val * val);
	}

	if (x < x2)
	{
		const Type val = (x - x2) / betta0;
		return a * exp(-val * val);
	}

	if (x < x3) 
	{
		return a;
	}

	const Type val = (x - x3) / betta1;
	return a * exp(-val * val) + b;	
}
Type Pressure(const Vector3& p)
{
	return 100;
	const Type x = p[0];
	const Type x1 = 0.04;
	const Type x2 = 0.06;
	const Type x3 = 0.1;

	const Type A = 200;
	const Type betta2 = 0.03;
	const Type betta3 = 0.02;
	const Type b = 1e-2;

	if (x < x1)
	{
		return A;
	}

	if (x < x2)
	{
		const Type val = (x - x1) / betta3;
		return A * exp(-val * val);
	}

	if (x < x3)
	{
		const Type val = (x2 - x1) / betta3;
		return A * exp(-val * val);
	}

	const Type val0 = (x2 - x1) / betta3;
	const Type coef = A * exp(-val0 * val0);

	const Type val = (x - x3) / betta2;	
	return coef * exp(-val * val)+b;
}
Vector3 Velocity(const Vector3& p)
{
	return Vector3(C_LIGHT * 1e-5, 0, 0);
}

static int SetRHllcValueDefault(std::vector<elem_t>& cells)
{
	std::vector<Vector3> centers;
	if (ReadSimpleFileBin(glb_files.base_adress + "centers.bin", centers)) RETURN_ERR("Default rhllc value not set\n");

	int i = 0;
	for (auto& el : cells)
	{
#if defined Cylinder
		Vector3 x = centers[i];
		if (Vector2(x[1], x[2]).norm() < 0.03 && x[0] < 0.1)
		{
			el.phys_val.d = 0.1;
			el.phys_val.p = 0.01;
			el.phys_val.v = Vector3(0.99, 0, 0);
		}
		else
		{
			el.phys_val.d = 10;
			el.phys_val.p = 0.01;
			el.phys_val.v = Vector3(0, 0, 0);
		}
#elif defined Cube  //Cylinder
		Type x = centers[i][0];
		if (x < 0.5)
		{
			el.phys_val.d = 1;
			el.phys_val.p = 1;
			el.phys_val.v = Vector3(0.9, 0, 0);
		}
		else
		{
			el.phys_val.d = 1;
			el.phys_val.p = 10;
			el.phys_val.v = Vector3(0, 0, 0);
		}
#elif defined Cone_JET
		Vector3 x = centers[i];
		if (Vector2(x[1], x[2]).norm() < 0.01 && x[0] < 0.05)
		{
			el.phys_val.d = 0.1;
			el.phys_val.p = 0.01;
			el.phys_val.v = Vector3(0.99, 0, 0);
		}
		else
		{
			el.phys_val.d = 10;
			el.phys_val.p = 0.01;
			el.phys_val.v = Vector3(0, 0, 0);
		}

		//const Vector3 x = centers[i];
		//el.phys_val.d = Density(x) / DENSITY;
		//el.phys_val.p = Pressure(x) / PRESSURE;
		//el.phys_val.v = Velocity(x) / VELOCITY;
#elif defined Cone

		{
			el.phys_val.d = 1;
			el.phys_val.p = 0.1;
			el.phys_val.v = Vector3(0, 0, 0);
		}

#elif defined Sphere
		el.phys_val.d = 0.1;
		el.phys_val.p = 0.1;
		el.phys_val.v = Vector3(0, 0, 0);
#else		
		el.phys_val.d = 0.1;
		el.phys_val.p = 1;
		el.phys_val.v = Vector3(1e-4, 0, 0);
#endif // Cube

		i++;
	} //for

	return 0;
}

// ------------------
static int SetHllcSettingDefault(hllc_value_t& hllc_set)
{
	hllc_set.h = 0.004; // default
#ifdef ILLUM
	///hllc_set.h = 0.0007166575761593; //jet
	hllc_set.tau = 1e-7;
	hllc_set.CFL = 0.01;
	hllc_set.print_timer = 0.0001;
	hllc_set.T = 0.0002;// 1;
#else //ILUM

	hllc_set.h = 0.0005369329546790;
	hllc_set.tau = 1e-5;
	hllc_set.CFL = 0.01;
	hllc_set.print_timer = 1e-5;
	hllc_set.T = 1e-5*2-1e-6;

#if defined Cube
	//const Type h = 0.0007123669658939; // Soda1d_2
	//const Type h  = 0.0010828369115320; // Soda1d
	hllc_set.h = 0.0010307259619874; // Soda1d_3
#elif defined Step
	hllc_set.h = 0.0018751819368151;
	hllc_set.T = 5;
	hllc_set.CFL = 0.5;
	hllc_set.print_timer = 0.1;
#endif

#endif //ILUM

	return 0;
}

int InitRHLLC(file_name file_settings_hllc, hllc_value_t& hllc_set,
	file_name file_init_value, std::vector<elem_t>& cells)
{
	if (ReadHllcSetting(file_settings_hllc, hllc_set))
	{
		SetHllcSettingDefault(hllc_set);
		WRITE_LOG("SetDefault hllc settings\n");
	}

	if (ReadHllcInit(file_init_value, cells))
	{
		if (SetRHllcValueDefault(cells))
		{
			WRITE_LOG("Default rhllc value not set\n");
			return 1;
		}
		WRITE_LOG("SetDefault hllc value\n");
	}

	RHllcPhysToConv(cells);

	RHllcGetTimeStep(hllc_set, cells);

	return 0;
}



#if NUMBER_OF_MEASUREMENTS == 2
int ReBuildPhysicValue(const Vector4& U, Vector4& W)
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
		printf("Error (p = %lf, d= %lf)", p, D / Gamma0);
		return 1;
	}

	W(0) = D / Gamma0;
	W(1) = v(0);
	W(2) = v(1);
	W(3) = p;

	return 0;
}

int ReBuildPhysicValue(const  std::vector<Vector4>& U, std::vector<Vector4>& W) {

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

int ReBuildConvValue(const std::vector<Vector4>& W, std::vector<Vector4>& U) {

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
#endif

#endif //RHLLC