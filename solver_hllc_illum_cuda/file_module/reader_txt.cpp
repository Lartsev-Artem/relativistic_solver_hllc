#include "../prj_config.h"

#include "../global_headers.h"
#include "../global_def.h"


int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all) {

	std::ifstream ifile;

	OPEN_FSTREAM(ifile, name_file_sphere_direction.c_str());

	int N = 0;
	ifile >> N;
	directions_all.resize(N);

	Type buf_s;
	int i = 0;

	for (int i = 0; i < N; i++) {
		ifile >> buf_s;
		ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
	}
	ifile >> buf_s;
	ifile.close();
	return 0;
}

size_t ReadSphereDirectionDecartToSpherical(const std::string& name_file_sphere_direction, vector<Vector3>& directions_all, vector<Type>& squares, Type& square_surface) {

	ifstream ifile;

	OPEN_FSTREAM(ifile, name_file_sphere_direction.c_str());

	int N = 0;
	ifile >> N;
	directions_all.resize(N);
	squares.resize(N);

	for (int i = 0; i < N; i++) {
		ifile >> squares[i];
		ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
	}
	ifile >> square_surface;
	ifile.close();

	return 0;
}
