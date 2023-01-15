#include "../prj_config.h"

#if defined BUILD

#include "../global_headers.h"
#include "../global_def.h"

#include "reader_vtk.h"

#ifdef USE_VTK

#endif // USE_VTK


int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals) {

	FILE* f;
	f = fopen(name_file_normals.c_str(), "rb");

	int n;
	fread_unlocked(&n, sizeof(int), 1, f);
	normals.resize(n);

	Normals norm(4);
	for (size_t i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++)
			fread_unlocked(norm.n[j].data(), sizeof(Type), 3, f);
		normals[i] = norm;
	}
	fclose(f);

	printf("Normals read\n");


	//std::ifstream ifile;

	//ifile.open(name_file_normals);
	//if (!ifile.is_open()) {
	//	std::cout << "Error read file normals\n";
	//	return 1;
	//}

	//int N;
	//ifile >> N;
	//normals.resize(N);

	//Normals norm(4);
	//for (size_t i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//		ifile >> norm.n[j][0] >> norm.n[j][1] >> norm.n[j][2];
	//	normals[i] = norm;
	//}

	//ifile.close();
	return 0;
}

#endif //BUILD