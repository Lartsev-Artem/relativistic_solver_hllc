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
	fread(&n, sizeof(int), 1, f);
	normals.resize(n);

	Normals norm(base);
	for (size_t i = 0; i < n; i++)
	{
		for (int j = 0; j < base; j++)
			fread(&norm.n[j], sizeof(Vector3), 1, f);

		normals[i] = norm;
	}
	fclose(f);

	return 0;
}

#endif //BUILD