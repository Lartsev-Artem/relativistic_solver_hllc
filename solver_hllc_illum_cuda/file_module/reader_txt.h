#include "../prj_config.h"

#if (!defined READER_TXT)

#define READER_TXT

#include "../global_def.h"

template<typename Str_Type, typename T>
size_t ReadSimpleFileTxt(const Str_Type name_file, std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	std::ifstream ifile;
	OPEN_FSTREAM(ifile, name_file.c_str());

	int size;
	ifile >> size;
	data_array.resize(size);

	for (int i = 0; i < size; i++)
	{
		ifile >> data_array[i];
	}

	ifile.close();
	return 0;
}

template<typename Str_Type, typename T>
size_t WriteSimpleFileTxt(const Str_Type name_file, std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, name_file.c_str());
	
	ofile << data_array.size()<<'\n';

	for (int i = 0; i < data_array.size(); i++)
	{
		ofile << data_array[i] << '\n';
	}

	ofile.close();
	return 0;
}

int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all);
size_t ReadSphereDirectionDecartToSpherical(const std::string& name_file_sphere_direction, grid_directions_t& grid_direction);


#include "../solve_module/solve_config.h"
#if defined HLLC || defined RHLLC
#include "../solve_module/solve_config.h"
#include "../solve_module/solve_global_struct.h"
int ReadHllcSetting(file_name file_settings_hllc, hllc_value_t& hllc_set);

#endif //SOLVE
#endif // !READER_TXT
