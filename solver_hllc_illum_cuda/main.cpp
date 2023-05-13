#include "prj_config.h"
#include "global_def.h"
#include "global_value.h"
#include "global_headers.h"

global_files_t glb_files;
int main(int argc, char* argv[])
{	
	MPI_START(argc, argv);	

	std::string name_file_settings;
	if (argc <= 1)
	{
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file.txt";
		std::string foo; int a;
		if (ReadStartSettings(name_file_settings, a, foo, foo, foo, foo, glb_files.base_adress, foo, a))
		{
			RETURN_ERR("Err reading default settings file\n");
		}		
	}
	else
		name_file_settings = argv[1];
	
	int id, np;
	MPI_GET_INF(np, id);
	if(id == 0) std::remove(Files_log);
		
	std::remove((Files_log + std::to_string(id) + ".txt").c_str());

#if defined BUILD
	RunBuildModule(name_file_settings);
#endif //BUILD

#if defined MAKE
	RunMakeModule(name_file_settings, 0, 0);
#endif //MAKE

#if defined SOLVE	
	RunSolveModule(argc, argv, name_file_settings);
#endif //SOLVE

#if defined UTILS	
	RunUtilsModule(argc, argv, name_file_settings);	
#endif //UTILS

	MPI_END;
	return EXIT_SUCCESS;
}
