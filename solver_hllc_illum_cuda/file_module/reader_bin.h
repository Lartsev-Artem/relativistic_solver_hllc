#include "../prj_config.h"

#if (!defined READER_BIN)

#define READER_BIN

#include "../global_def.h"
template<typename Str_Type, typename T>
size_t ReadSimpleFileBin(const Str_Type name_file, std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	FILE* f;
	f = fopen(std::string(name_file).c_str(), "rb");

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	data_array.resize(n);
	fread_unlocked(data_array.data(), sizeof(T), n, f);

	fclose(f);

	printf("read simple data: %s \n", name_file.c_str());
	return 0;
}


int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals);
#endif // !READER_BIN
