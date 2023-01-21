#include "prj_config.h"
#include "global_def.h"
#include "global_value.h"
#include "global_headers.h"

std::string BASE_ADRESS;

int main(int argc, char* argv[])
{
	MPI_START(argc, argv);

	std::string name_file_settings;
	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file.txt";
	else
		name_file_settings = argv[1];

#if defined BUILD
	RunBuildModule(name_file_settings);
#endif //BUILD

#if defined MAKE
	RunMakeModule(name_file_settings, 0, 0);
#endif //MAKE

#if defined SOLVE
	RunSolveModule(name_file_settings);
#endif //SOLVE

#if defined UTILS
	RunUtilsModule(argc, argv);
#endif

	MPI_END;
	return EXIT_SUCCESS;
}
