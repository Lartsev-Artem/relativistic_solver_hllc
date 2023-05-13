#include "hllc_utils.h"
#include "../../global_def.h"
#if defined HLLC && defined SOLVE
#include "../../file_module/reader_txt.h"
#include "../../file_module/reader_bin.h"


int HllcPhysToConv(std::vector<elem_t>& cells)
{
	for (auto& el : cells)
	{
		const Type v = el.phys_val.v.dot(el.phys_val.v);
		const Type d = el.phys_val.d;

		el.conv_val.d = d;
		el.conv_val.v *= d;
		el.conv_val.p = el.phys_val.p / (gamma1 - 1) + d * v / 2;
	}

	return 0;
}

int HllcConvToPhys(std::vector<elem_t>& cells)
{
	for (auto& el : cells)
	{
		const Type d = el.conv_val.d;
		el.phys_val.d = d;
		el.phys_val.v = el.conv_val.v / d;
		const Type vv = el.phys_val.v.dot(el.phys_val.v);
		el.phys_val.p = (el.conv_val.p - vv * d / 2.) * (gamma1 - 1);
	}
	return 0;
}

int HllcGetTimeStep(hllc_value_t& hllc_set, const std::vector<elem_t>& cells) {

	Type c_max = -1;
	for (auto& el : cells)
	{				
		const Type a = sqrt(gamma1 * el.phys_val.p / el.phys_val.d);
		const Type c = el.phys_val.v.norm() + a;
		if (c > c_max) c_max = c;
	}
	hllc_set.tau = hllc_set.CFL * hllc_set.h / c_max;

	if (hllc_set.tau < 0)
	{
		EXIT_ERRS("Error tau = %lf\n", hllc_set.tau);
	}

	return  0;
}


static int SetHllcValueDefault(std::vector<elem_t>& cells)
{
	std::vector<Vector3> centers;
	if (ReadSimpleFileBin(glb_files.base_adress + "centers.bin", centers)) RETURN_ERR("Default hllc value not set\n");

	int i = 0;
	for (auto& el : cells)
	{
#ifdef Jet
		const Type betta = 0.1;
		const Type a = 1;
		const Type b = 0.001;
		Type x = centers[i][0];
		el.phys_val.d = a * exp(-x * x / betta) + b;
		el.phys_val.p = a * exp(-x * x / betta) + (1e-5);
		el.phys_val.v = Vector3(1e-4, 0, 0);
#endif

#ifdef Cube 
		Type x = centers[i][0];
		if (x < 0.5)
		{
			el.phys_val.d = 1;
			el.phys_val.p = 1;
			el.phys_val.v = Vector3(0, 0, 0);
		}
		else
		{
			el.phys_val.d = 0.125;
			el.phys_val.p = 0.1;
			el.phys_val.v = Vector3(0, 0, 0);
		}
#endif // Cube

#if defined Sphere || defined Step
		el.phys_val.d = 0.1;
		el.phys_val.p = 0.01;
		el.phys_val.v = Vector3(0, 0, 0);
#endif //Sphere

		i++;
	}

	return 0;
}
static int SetHllcSettingDefault(hllc_value_t& hllc_set)
{
	hllc_set.h = 0.005; // default
#ifdef ILLUM
	hllc_set.tau = 1e-8;
	hllc_set.CFL = 0.001;
	hllc_set.print_timer = 1e-5;
	hllc_set.T = 0.5;
#else //ILUM

	hllc_set.tau = 1e-5;
	hllc_set.CFL = 0.5;
	hllc_set.print_timer = 0.01;
	hllc_set.T = 0.4;

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

int InitHLLC(file_name file_settings_hllc, hllc_value_t& hllc_set, 
	file_name file_init_value, std::vector<elem_t>& cells)
{
	if (ReadHllcSetting(file_settings_hllc, hllc_set))
	{
		SetHllcSettingDefault(hllc_set);
		WRITE_LOG("SetDefault hllc settings\n");
	}

	if (ReadHllcInit(file_init_value, cells))
	{
		if (SetHllcValueDefault(cells))
		{
			WRITE_LOG("Default rhllc value not set\n");
			return 1;
		}
		WRITE_LOG("SetDefault hllc value\n");
	}

	HllcPhysToConv(cells);

	HllcGetTimeStep(hllc_set, cells);

	return 0;
}

#endif //HLLC