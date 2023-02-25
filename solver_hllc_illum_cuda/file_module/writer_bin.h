#ifndef WRITER_BIN
#define WRITER_BIN

#include "../prj_config.h"

#include "../global_def.h"
#include "../solve_module/solve_global_struct.h"

template<typename Str_Type, typename T>
size_t WriteSimpleFileBin(const Str_Type name_file, const std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	FILE* f;
	f = fopen(std::string(name_file).c_str(), "wb");

	if (!f)
	{
		RETURN_ERRS("file %s wasn't writed\n", std::string(name_file).c_str());
	}

	int n = data_array.size();
	fwrite(&n, sizeof(int), 1, f);	
	fwrite(data_array.data(), sizeof(T), n, f);
	
	fclose(f);
	return 0;
}

size_t WriteFileSolution(const std::string& main_dir, const grid_t& grid);
int WriteGeometryGrid(const std::string& file_cells, const std::string& file_faces, grid_t& grid);

#endif // !WRITER_BIN
