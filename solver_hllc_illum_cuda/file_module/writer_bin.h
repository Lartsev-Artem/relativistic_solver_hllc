#include "../prj_config.h"

#if (!defined WRITER_BIN)

#define WRITER_BIN

#include "../global_def.h"
template<typename Str_Type, typename T>
size_t WriteSimpleFileBin(const Str_Type name_file, const std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	FILE* f;
	f = fopen(std::string(name_file).c_str(), "rb");

	if (!f)
	{
		RETURN_ERR(("file %s wasn't writed\n", std::string(name_file).c_str()));
	}

	int n = data_array.size();
	fwrite(&n, sizeof(int), 1, f);	
	fwrite(data_array.data(), sizeof(T), n, f);
	
	fclose(f);
	return 0;
}

#define WRITE_FILE(name_file, data, value) \
{ \
FILE* f;\
int n = cells.size(); \
f = fopen(name_file, "wb"); \
fwrite(&n, sizeof(int), 1, f); \
for (auto& el : data) \
{	\
	fwrite(&el.value, sizeof(el.value), 1, f);	\
}	\
fclose(f);\
}


#endif // !WRITER_BIN
