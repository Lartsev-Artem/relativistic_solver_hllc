#include "../prj_config.h"

#if (!defined READER_BIN)

#define READER_BIN

#include "../global_def.h"
template<typename T>
size_t ReadSimpleFileBin(file_name name_file, std::vector<T>& data_array) {
	// Файл должен содержать в первой строке число элементов. Далее последовательные данные

	FILE* f;
	f = fopen(name_file.c_str(), "rb");
	if (!f) RETURN_ERRS("file %s not open\n", name_file.c_str());

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	data_array.resize(n);
	fread_unlocked(data_array.data(), sizeof(T), n, f);

	fclose(f);

	printf("read simple data: %s \n", name_file.c_str());
	return 0;
}

#define READ_FILE(name_file, data, value) \
{ \
FILE* f;\
f = fopen(name_file, "rb"); \
if(!f) RETURN_ERRS("file %s not open\n",name_file); \
int n; \
fread(&n, sizeof(int), 1, f); \
data.resize(n); \
for (auto& el : data) \
{	\
	fread(&el.value, sizeof(el.value), 1, f);	\
}	\
fclose(f);\
}


int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals);


int ReadDataArray(const size_t class_file_vtk, const std::string& main_dir,
	std::vector<Type>& density, std::vector<Type>& absorp_coef, std::vector<Type>& rad_en_loose_rate,
	std::vector<Vector3>& velocity, std::vector<Type>& pressure, const bool is_print=false);

#ifdef SOLVE
#include "../solve_module/solve_global_struct.h"
int ReadDataArray(const solve_mode_t& mode, file_name main_dir, grid_t& grid);

int ReadGeometryGrid(const std::string& file_cells, const std::string& file_faces, grid_t& grid);
int ReadValueGrid(const std::string& main_dir, grid_t& grid);
int ReadIllumGeometry(const int count_dir,
	file_name file_x, file_name file_state, file_name file_x0, file_name file_graph, file_name file_res,
	std::vector<BasePointTetra>& vec_x,
	std::vector <std::vector<int>>& face_states,
	std::vector <std::vector<cell_local>>& vec_x0,
	std::vector < std::vector<int>>& sorted_id_cell,
	std::vector<Type>& vec_res_bound);

#if defined HLLC || defined RHLLC
inline int ReadHllcInit(file_name file_init_value, std::vector<elem_t>& cells)
{
	READ_FILE(file_init_value.c_str(), cells, phys_val);
	return 0;
}
#endif

#endif // Solve
#endif // !READER_BIN
