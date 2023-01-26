#include "prj_config.h"
#include "global_def.h"
#include "global_value.h"
#include "global_headers.h"

std::string BASE_ADRESS;
#include "utils/rebuild_solve.h"
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
	Type L = 1 * 1e12;
	Type V = 3 * 1e8;
	Type M = 1 * 1e21;

	Type t = L / V;

	Type d = M / (L * L * L);
	Type p = M / (L * t * t);

	Type e = M * L * L / (t * t);
	Type I = e / (t * L * L);
	Type data[10] = {L,V,M,t,0,d,p,e,I};

	for (size_t i = 0; i < 10; i++)
	{
		std::cout << std::setprecision(16) << data[i] << '\n';
	}
	return 0;


	rebuild_solve(name_file_settings);//RunUtilsModule(argc, argv, name_file_settings);
#endif

	MPI_END;
	return EXIT_SUCCESS;
}
