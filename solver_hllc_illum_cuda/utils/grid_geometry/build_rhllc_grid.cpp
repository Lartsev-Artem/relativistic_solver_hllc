#if defined UTILS && defined BUILD
#include "build_rhllc_grid.h"

#if defined USE_VTK
#include "../../file_module/reader_vtk.h"

#include "../../file_module/reader_txt.h"
#include "../../file_module/reader_bin.h"
int aaaaa(int np, const std::string& base_adress)
{
	const std::string name_file_id_neighbors = base_adress + F_NEIB;
	const std::string name_file_normals = base_adress + F_NORMALS;
	const std::string name_file_centers = base_adress + F_CENTERS;
	const std::string name_file_squares = base_adress + F_SQUARES;
	const std::string name_file_volume = base_adress + F_VOLUME;
	const std::string name_file_metis = "D:\\Desktop\\FilesCourse\\Metis\\plane.txt.epart.txt";

	std::vector<int> neighbours_id_faces;
	std::vector<Normals> normals;
	std::vector<Type> squares_faces;
	std::vector<Type> volume;
	std::vector<Vector3> centers;
	std::vector<uint8_t> cells_nodes;

	if (ReadSimpleFileBin(name_file_id_neighbors, neighbours_id_faces)) RETURN_ERR("Error reading file neighbours\n");
	if (ReadNormalFile(name_file_normals, normals)) RETURN_ERR("Error reading file normals\n");
	if (ReadSimpleFileBin(name_file_squares, squares_faces)) RETURN_ERR("Error reading file squares_faces\n");
	if (ReadSimpleFileBin(name_file_volume, volume)) RETURN_ERR("Error reading file volume\n");
	if (ReadSimpleFileBin(name_file_centers, centers)) RETURN_ERR("Error reading file centers\n");	
	if (ReadSimpleFileTxt(name_file_metis, cells_nodes))RETURN_ERR("Error reading file metis\n");

	std::vector<int> loc_size(np, 0);
	for (auto id: cells_nodes)
	{		
		if (id >= np)
		{
			WRITE_LOG_ERR("Error config claster! Need " << np << "nodes\n");
			D_LD; //если файл metis собран под другую конфинурацию будет ошибка!		
		}
		loc_size[id]++;
	}
#ifdef DEBUG
	for (auto N : loc_size) if (N == 0)
	{
		WRITE_LOG_ERR("empty node in mpi_rhllc\n");
		D_LD;
	}
#endif

	std::vector<std::vector<face_t>> faces(np);
	std::vector<std::vector<elem_t>> cells(np);

	for (int id = 0; id < np; id++)
	{
		const int N = loc_size[id];
		cells[id].resize(N);
	}

	int cc = 0;
	const int N = centers.size();

	std::vector<int> idx_cells(np, 0);

	for (int i = 0; i < N * base; i++)
	{
		int idx = neighbours_id_faces[i];

		if (idx != -10)
		{
			face_t f;
			f.geo.id_l = idx_cells[i / base]++; //€чейка

			f.geo.n = normals[i / base].n[i % base];
			f.geo.S = squares_faces[i];

			cells[i / base].geo.sign_n[i % base] = true;
			cells[i / base].geo.id_faces[i % base] = cc;

			neighbours_id_faces[i] = -10;
			if (idx >= 0)
			{
				f.geo.id_r = idx / base; //сосед

				cells[idx / base].geo.sign_n[idx % base] = false;
				cells[idx / base].geo.id_faces[idx % base] = cc;

				neighbours_id_faces[idx] = -10;
			}
			else
			{
				f.geo.id_r = idx; // код границы
			}

			faces[cells_nodes[i/base]].push_back(f); // как потом искать с €чейками?
			cc++;
		}
	

	return EXIT_SUCCESS;
}

#endif //USE_VTK
#endif //UTILS
