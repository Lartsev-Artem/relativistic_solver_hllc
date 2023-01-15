#include "Header.h"
#define base 4
struct flux
{
	Type d;
	Vector3 v;
	Type p;

	flux()
	{
		d = 0;
		v = Vector3::Zero();
		p = 0;
	}

	flux(const Type a, const Type b, const Type c, const Type dd, const Type e)
	{
		d = a;
		v = Vector3(b, c, dd);
		p = e;
	}

	flux operator+ (const flux& x) { d += x.d; v += x.v; p += x.p; return *this; }
	flux operator+= (const flux& x) { d += x.d; v += x.v; p += x.p; return *this; }
	flux operator-= (const flux& x) { d -= x.d; v -= x.v; p -= x.p; return *this; }
	flux operator* (const Type x) { d *= x; v *= x; p *= x; return *this; }
	flux operator- (const flux& x) { d -= x.d; v -= x.v; p -= x.p; return *this; }
	flux operator/ (const Type x) { d /= x; v /= x; p /= x; return *this; }

	// это временно дл€ свз€и со старым кодом
	Type operator[](const int i)
	{
		switch (i)
		{
		case 0: return d;
		case 1: return v[0];
		case 2: return v[1];
		case 3: return v[2];
		case 4: return p;

		default:
			printf("err idx in flux opreator []\n");
			return 0;
			break;
		}
	}
	Type operator()(const int i)
	{
		switch (i)
		{
		case 0: return d;
		case 1: return v[0];
		case 2: return v[1];
		case 3: return v[2];
		case 4: return p;

		default:
			printf("err idx in flux opreator ()\n");
			return 0;
			break;
		}
	}

	//private:
	//	flux(const flux& f) {};
};
struct face
{
	flux f;
	int id_l;
	int id_r;

	Vector3 n;
	Type S;

	face()
	{
		id_l = 0;
		id_r = 0;
		n = Vector3::Zero();
		S = 0;
	}
};
struct elem// пока спорно
{
	flux val;
	int id_faces[base];
	Type V;
	bool sign_n[base];

	Vector3 center; //if define ??

	elem()
	{
		for (int i = 0; i < base; i++)
		{
			id_faces[i] = -1;
			sign_n[i] = 0;
		}
		V = 0;
	}

//private:
	elem(const elem& el) {};
};

int WriteNewStructFile(const std::string file_dir, const std::vector<face>& faces, const std::vector<elem>& cells);
int ReBuildFromOldToNewStruct(int argc, char* argv[])
{
	/*std::string name_file_normals = "";
	std::string name_file_squares = "";
	std::string name_file_volume = "";
	std::string name_file_id_neighbors = "";*/

	std::string main_dir = "D:\\Desktop\\FilesCourse\\Test\\";
	std::string name_file_id_neighbors = main_dir + "pairs";
	std::string name_file_normals = main_dir + "normals.bin";
	std::string name_file_centers = main_dir + "centers.bin";
	std::string name_file_squares = main_dir + "squares.bin";
	std::string name_file_volume = main_dir + "volume.bin";
	


	std::vector<Normals> normals;
	std::vector<Type> squares_faces;
	std::vector<Type> volume;
	std::vector<int> neighbours_id_faces;
	std::vector<Vector3> centers;

	if (ReadNormalFile(name_file_normals, normals)) ERR_RETURN("Error reading file normals\n")
		if (ReadSimpleFileBin(name_file_squares, squares_faces)) ERR_RETURN("Error reading file squares_cell\n")
			if (ReadSimpleFileBin(name_file_volume, volume)) ERR_RETURN("Error reading file volume\n")
				if (ReadSimpleFileBin(name_file_id_neighbors + ".bin", neighbours_id_faces)) ERR_RETURN("Error reading file neighbours\n")

	if (ReadSimpleFileBin(name_file_centers, centers)) ERR_RETURN("Error reading file centers\n")
	
	const int N = neighbours_id_faces.size() / base;
	
	std::vector<face> faces; faces.reserve(N);
	std::vector<elem> cells(N);
	

	cells.resize(N);

	int cc = 0;
	for (int i = 0; i < N * base; i++)
	{
		int idx = neighbours_id_faces[i];

		if (idx != -10)
		{
			face f;
			f.id_l = i / base; //€чейка

			f.n = normals[i / base].n[i % base];
			f.S = squares_faces[i];

			cells[i / base].sign_n[i % base] = true;
			cells[i / base].id_faces[i % base] = cc;

			neighbours_id_faces[i] = -10;
			if (idx >= 0)
			{
				f.id_r = idx / base; //сосед

				cells[idx / base].sign_n[idx % base] = false;
				cells[idx / base].id_faces[idx % base] = cc;

				neighbours_id_faces[idx] = -10;
			}
			else
			{
				f.id_r = idx; // код границы
			}

			faces.push_back(f); // как потом искать с €чейками?
			cc++;
		}
	}

	for (size_t i = 0; i < N; i++)
	{
		cells[i].V = volume[i];
		cells[i].center = centers[i];
	}
	volume.clear();

	neighbours_id_faces.clear();
	normals.clear();
	squares_faces.clear();

	//int id_face = 0;
	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < base; j++)
	//	{
	//		const int idx = neighbours_id_faces[i * N + j];

	//		if (idx != -10)
	//		{
	//			face f;
	//			f.id_l = i; //€чейка

	//			f.n = normals[i].n[j];
	//			f.S = squares_faces[i * N + j];

	//			cells[i].sign_n[j] = true;
	//			cells[i].id_faces[j] = id_face;

	//			neighbours_id_faces[i * N + j] = -10;
	//			if (idx >= 0)
	//			{
	//				f.id_r = idx / base; //сосед

	//				cells[idx / base].sign_n[idx % base] = false;
	//				cells[idx / base].id_faces[idx % base] = id_face;

	//				neighbours_id_faces[idx] = -10;
	//			}
	//			else
	//			{
	//				f.id_r = idx; // код границы
	//			}

	//			faces.push_back(f); // как потом искать с €чейками?
	//			id_face++;
	//		}
	//	}

	//	cells[i].V = volume[i];
	//}


	neighbours_id_faces.clear();
	normals.clear();
	squares_faces.clear();
	volume.clear();
	centers.clear();


	WriteNewStructFile(main_dir, faces, cells);

	std::vector<face> faces2;
	std::vector<elem> cells2;

	
	ReadSimpleFileBin(main_dir + "cells.bin", cells2);
	ReadSimpleFileBin(main_dir + "faces.bin", faces2);

	//compare
	{
		for (int i = 0; i < N; i++)
		{
			if(fabs(cells[i].V - cells2[i].V) > 1e-10) 
			{
				printf("err\n");
				exit(1);
			}

			if ((faces[i].n - faces2[i].n).norm() > 1e-10)
			{
				printf("err2\n");
				exit(1);
			}

			if (faces[i].id_l != faces2[i].id_l)
			{
				printf("err3\n");
				exit(1);
			}
		}
	}

	return 0;
}

int WriteNewStructFile( const std::string file_dir, const std::vector<face>& faces, const std::vector<elem>& cells)
{
	FILE* f;
	f = fopen((file_dir+"faces.bin").c_str(), "wb");

	int n = faces.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(faces.data(), sizeof(faces[0]), n, f);
	fclose(f);

	f = fopen((file_dir + "cells.bin").c_str(), "wb");

	n = cells.size();
	fwrite(&n, sizeof(int), 1, f);
	fwrite(cells.data(), sizeof(cells[0]), n, f);
	fclose(f);

	return 0;
}