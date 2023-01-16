#include "../prj_config.h"

#if (!defined READER_TXT)

#define READER_TXT

#include "../global_def.h"

int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all);
size_t ReadSphereDirectionDecartToSpherical(const std::string& name_file_sphere_direction, vector<direction_s>& directions_all, Type& square_surface);
#endif // !READER_TXT
