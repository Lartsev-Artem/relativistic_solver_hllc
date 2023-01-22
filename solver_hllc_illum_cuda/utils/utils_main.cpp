#include "../prj_config.h"
#ifdef UTILS
#include "rebuild_solve.h"

int RunUtilsModule(int argc, char* argv[], const std::string& settings)
{
	rebuild_solve(settings);
	return 0;
}
#endif //UTILS
