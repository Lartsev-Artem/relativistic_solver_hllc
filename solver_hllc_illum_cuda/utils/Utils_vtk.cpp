#include "Header.h"



size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	/*Чтение исходного файла и запись в vtkUnstructuredGrid*/

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();


	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	switch (class_file_vtk) {
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	case 5:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("pressure");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("velocity");  // GetScalars("radEnLooseRate");
		break;
	}

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}



inline void MakeRotationMatrix(const Vector3& n, Eigen::MatrixXd& T) {

	T = Eigen::MatrixXd::Zero(5, 5);
	T(0, 0) = T(4, 4) = 1;

	if (fabs(n[2] * n[2] - 1) > 1e-10)
	{

		T(1, 1) = n[0];
		T(1, 2) = n[1];
		T(1, 3) = n[2];

		Type sqr = sqrt(1 - n[2] * n[2]);

		if (sqr < eps * eps - eps / 10)
			printf("Err T\n");

		T(2, 1) = -n[1] / sqr;
		T(2, 2) = n[0] / sqr;

		T(3, 1) = -n[0] * n[2] / sqr;
		T(3, 2) = -n[1] * n[2] / sqr;
		T(3, 3) = sqr;
	}
	else if (n[2] > 0)  // n_z == 1
	{
		T(1, 3) = 1;
		T(2, 2) = 1;
		T(3, 1) = -1;
	}
	else  // n_z == -1
	{
		T(1, 3) = -1;
		T(2, 2) = -1;
		T(3, 1) = 1;
	}

}

void NewMakeRotationMatrix(const Vector3& n, Matrix3& T)
{	
	const Type x = n[0];
	const Type y = n[1];
	const Type theta = atan2(y, x) ;
	const Type phi = atan2(sqrt(x * x + y * y), n[2]);

	T(0, 0) = cos(theta) * sin(phi);
	T(0, 1) = sin(theta) * sin(phi);
	T(0, 2) = cos(phi);

	T(1, 0) = -sin(theta);
	T(1, 1) = cos(theta);
	T(1, 2) = 0;

	T(2, 0) = -cos(theta) * cos(phi);
	T(2, 1) = -sin(theta) * cos(phi);
	T(2, 2) = sin(phi);

	return;
}

int ReadNormalFile(const std::string& name_file_normals, std::vector<Normals>& normals) {

	FILE* f;
	f = fopen(name_file_normals.c_str(), "rb");

	int n;
	fread(&n, sizeof(int), 1, f);
	normals.resize(n);

	Normals norm(4);
	for (size_t i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++)
			fread(norm.n[j].data(), sizeof(Type), 3, f);
		normals[i] = norm;
	}
	fclose(f);
	return 0;
}

//#define CHECK_NEW_ROTATE


#define SIGN(a) (a < 0.0 ? -1.0 : 1.0) 

int ReBuildNetgenToMetis(int argc, char* argv[]) 
{
	std::string name_file_in = "";
	std::string name_file_out = "";
	int size = 3; // по умолчанию 3d сетка

	if (argc <= 2)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input: path\\input_file, path\\output_file.vtk\n");
		return 1;
#else
		name_file_in = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.txt";
		name_file_out = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.vtk";
#endif
	}
	else
	{
		name_file_in = argv[1];
		name_file_out = argv[2];

		if (argc == 4)
		{
			size = std::stoi(argv[3]);
		}
	}


	std::ifstream ifile(name_file_in);

	if (!ifile.is_open())
	{
		printf("Error. Input file not opened\n");
		return 1;
	}

	int n;
	ifile >> n;
	std::vector<Eigen::Vector3d> point(n);
	
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			ifile >> point[i][j];
		}
	}

	ifile >> n;
	std::vector<Eigen::Vector4i> cell(n);
	for (size_t i = 0; i < n; i++)
	{
		ifile >> cell[i][0];
		for (size_t j = 0; j < size+1; j++)
		{
			ifile >> cell[i][j];
		}
	}


	std::ofstream ofile(name_file_out);

	if (!ofile.is_open())
	{
		printf("Error. Output file not opened\n\n");
		return 1;
	}

	ofile  << cell.size() << "\n";
	for (size_t i = 0; i < cell.size(); i++)
	{
		for (size_t j = 0; j < size + 1; j++)
		{
			ofile << cell[i][j] << ' ';
		}
		ofile << '\n';
	}	

	ofile.close();
	return 0;
}

#include <filesystem>
#include <fstream>
#include <string>

int RenameFile(int argc, char* argv[])
{
	std::string name_file = "";
	std::string name_file_in = "";
	std::string name_file_out = "";
	if (argc <= 1)
	{
#ifdef RELEASE_VER
		printf("Error input data!\n");
		printf("Input: path\\file, with \n path\\output_file_out \n path\\output_file\n");
		return 1;
#else
		name_file_in = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.txt";
		name_file_out = "D:\\Desktop\\FilesCourse\\Grids\\plane2d.vtk";
#endif
	}
	else
	{
		name_file = argv[1];
	}

	std::ifstream ifile(name_file);

	if (!ifile.is_open())
	{
		printf("Error. Input file not opened\n");
		return 1;
	}
	std::getline(ifile, name_file_in);
	std::getline(ifile, name_file_out);
	ifile.close();

	std::remove(name_file_out.c_str());
	//std::filesystem::copy_file(name_file_in.c_str(), name_file_out.c_str());

	return 0;
}

int RenameFile(std::string name_file_in, std::string name_file_out)
{
	std::remove(name_file_out.c_str());
	std::rename(name_file_in.c_str(), name_file_out.c_str());
	return 0;
}
int SetFile( const std::string& file_set, const std::string& file_in, const std::string& file_out)
{
	std::ofstream ofile(file_set);

	if (!ofile.is_open())
	{
		printf("Error. Input file not opened\n");
		return 1;
	}

	ofile << file_in << '\n';
	ofile << file_out;

	ofile.close();
	return 0;
}




int main(int argc, char* argv[])
{	
	ReBuildFromOldToNewStruct(argc, argv);
	return 0;
#ifdef CHECK_NEW_ROTATE
	Eigen::MatrixXd T;
	Matrix3 T1;
	Vector3 n(0, 1, 0);
	MakeRotationMatrix(n, T);
	NewMakeRotationMatrix(n, T1);

	Vector3 v(0, 0, 5);
	Eigen::VectorXd V(5); V << 1, 0, 0, 5, 1;

	std::cout << T * V << '\n';
	std::cout << T1 * v << '\n';


	return 0;
#endif

#ifdef CHECK_NEW_ROTATE
	std::vector<Normals> normals;
	ReadNormalFile("D:\\Desktop\\FilesCourse\\Soda\\normals.bin", normals);
	Eigen::MatrixXd T;
	Matrix3 T1;

	for (size_t i = 0; i < normals.size(); i++)
	{
		for (size_t j = 0; j < 4; j++) 
		{
			MakeRotationMatrix(normals[i].n[j], T);
			NewMakeRotationMatrix(normals[i].n[j], T1);

			Eigen::VectorXd U = Eigen::VectorXd::Random(5);
			Eigen::VectorXd U1 = U;
			Vector3 UU(U[1], U[2], U[3]);

			U=T*U;
			UU = T1 * UU; 
			U1[1] = UU[0]; U1[2] = UU[1]; U1[3] = UU[2];

			if ((U - U1).norm() > 1e-9)
			{
				printf("Error new mat %d - %d\n", i, j);
				std::cout << "T= " << U << "\n\n T1= " << U1 << '\n';
				std::cout << "n= " << normals[i].n[j] << '\n';
			}
		}
	}
	

	return 0;

#endif

	//GetAverageArea(argc, argv);
	//ReNumberingGrid(argc, argv);

	//return 0;

	//Trace2D(argc, argv);

	//RunMake1d(argc, argv);
	//return 0;

	/*std::string file_trace = "D:\\Desktop\\FilesCourse\\Plane2D\\trace0.txt";
	std::string file_centers = "D:\\Desktop\\FilesCourse\\Plane2D\\centers.bin"; 
	std::string file_solve = "D:\\Desktop\\FilesCourse\\Plane2D\\Solve\\Solve";
	int max_iter = 14;

	Make1dFromTrace(file_trace, file_centers, file_solve, max_iter);

	return 0;*/

	//ReBuildNetgen2dToVTK(argc, argv);
	//GetAverageArea(argc, argv);
	
	//SetScalarDataVtkFromFile(argc, argv);
	//return 0;
	
	//SetFile("SetRename.txt", "Utils_vtk.exe", "ReBuildNetgenToMetis.exe");
	//RenameFile("D:\\BMSTU\\Utils_vtk\\x64\Release\\Utils_vtk.exe", "D:\\BMSTU\\Utils_vtk\\x64\\Release\\ReBuildNetgenToMetis.exe");
	ReBuildNetgenToMetis(argc, argv);
	return 0;
	
	BuildHLLC_1dTime(argc, argv);

	//SetBoundCells(argc, argv);

	return 0;
}



#ifdef HYDRO
/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *      const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   - Wl,Wr = L/R-states of PRIMITIVE variables at cell interface
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
	const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS* pFlux)
{
	Cons1DS Fl, Fr, Fhll, Uhll, Usl, Usr;
	Real rhl, rhr, csl, csr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
	Real ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
	Real lmdal, lmdar; /* Left and Right wave speeds */
	Real lmdas; /* Contact wave speed */
	Real ovlrmll;
	Real a, b, c, quad, rad;
	Real den, ps; /* PressUre in inner region */

	/*--- Step 1. ------------------------------------------------------------------
	 * Compute the max and min wave speeds used in Mignone
	 */

	rhl = Wl.d + Wl.P * Gamma / Gamma_1; /* Mignone Eq 3.5 */
	rhr = Wr.d + Wr.P * Gamma / Gamma_1;

	csl = sqrt(Gamma * Wl.P / rhl); /* Mignone Eq 4 */
	csr = sqrt(Gamma * Wr.P / rhr);

	cslsq = SQR(csl);
	csrsq = SQR(csr);

	vsql = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
	vsqr = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);

	gammasql = 1.0 / (1.0 - vsql);
	gammasqr = 1.0 / (1.0 - vsqr);

	ssl = cslsq / (gammasql * (1.0 - cslsq)); /* Mignone Eq 22.5 */
	ssr = csrsq / (gammasqr * (1.0 - csrsq));

	radl = sqrt(ssl * (1.0 - SQR(Wl.Vx) + ssl)); /* Mignone Eq 23 (radical part) */
	radr = sqrt(ssr * (1.0 - SQR(Wr.Vx) + ssr));

	lmdapl = (Wl.Vx + radl) / (1.0 + ssl); /* Mignone Eq 23 */
	lmdapr = (Wr.Vx + radr) / (1.0 + ssr);
	lmdaml = (Wl.Vx - radl) / (1.0 + ssl);
	lmdamr = (Wr.Vx - radr) / (1.0 + ssr);

	lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
	lmdar = MAX(lmdapl, lmdapr);


	/*--- Step 2. ------------------------------------------------------------------
	 * Compute L/R fluxes according to Mignone 2
	 */

	Fl.d = Ul.d * Wl.Vx;
	Fl.Mx = Ul.Mx * Wl.Vx + Wl.P;
	Fl.My = Ul.My * Wl.Vx;
	Fl.Mz = Ul.Mz * Wl.Vx;
	Fl.E = Ul.Mx;

	Fr.d = Ur.d * Wr.Vx;
	Fr.Mx = Ur.Mx * Wr.Vx + Wr.P;
	Fr.My = Ur.My * Wr.Vx;
	Fr.Mz = Ur.Mz * Wr.Vx;
	Fr.E = Ur.Mx;

	/*--- Step 3. ------------------------------------------------------------------
	 * Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
	 * Compute HLL conserved quantities using Mignone eq 9
	 */

	ovlrmll = 1.0 / (lmdar - lmdal);
	lmdatlmda = lmdal * lmdar;

	Fhll.d = (lmdar * Fl.d - lmdal * Fr.d + lmdatlmda * (Ur.d - Ul.d)) * ovlrmll;
	Fhll.Mx = (lmdar * Fl.Mx - lmdal * Fr.Mx + lmdatlmda * (Ur.Mx - Ul.Mx)) * ovlrmll;
	Fhll.My = (lmdar * Fl.My - lmdal * Fr.My + lmdatlmda * (Ur.My - Ul.My)) * ovlrmll;
	Fhll.Mz = (lmdar * Fl.Mz - lmdal * Fr.Mz + lmdatlmda * (Ur.Mz - Ul.Mz)) * ovlrmll;
	Fhll.E = (lmdar * Fl.E - lmdal * Fr.E + lmdatlmda * (Ur.E - Ul.E)) * ovlrmll;

	Uhll.d = (lmdar * Ur.d - lmdal * Ul.d + Fl.d - Fr.d) * ovlrmll;
	Uhll.Mx = (lmdar * Ur.Mx - lmdal * Ul.Mx + Fl.Mx - Fr.Mx) * ovlrmll;
	Uhll.My = (lmdar * Ur.My - lmdal * Ul.My + Fl.My - Fr.My) * ovlrmll;
	Uhll.Mz = (lmdar * Ur.Mz - lmdal * Ul.Mz + Fl.Mz - Fr.Mz) * ovlrmll;
	Uhll.E = (lmdar * Ur.E - lmdal * Ul.E + Fl.E - Fr.E) * ovlrmll;

	/*--- Step 4. ------------------------------------------------------------------
	 * Compute contact wave speed using larger root from Mignone Eq 18
	 * Physical root is the root with the minus sign
	 */

	 /* quadratic formUla calcUlation */

	a = Fhll.E;
	b = -(Uhll.E + Fhll.Mx);
	c = Uhll.Mx;


	quad = -0.5 * (b + SIGN(b) * sqrt(b * b - 4.0 * a * c));
	lmdas = c / quad;

	/*--- Step 5. ------------------------------------------------------------------
	 * Determine intercell flux according to Mignone 13
	 */

	if (lmdal >= 0.0) { /* Fl */
		/* intercell flux is left flux */
		pFlux->d = Fl.d;
		pFlux->Mx = Fl.Mx;
		pFlux->My = Fl.My;
		pFlux->Mz = Fl.Mz;
		pFlux->E = Fl.E;

		return;
	}
	else if (lmdas >= 0.0) { /* Fls */

		/* Mignone 2006 Eq 48 */
		ps = -Fhll.E * lmdas + Fhll.Mx;

		/* now calcUlate Usl with Mignone Eq 16 */
		den = 1.0 / (lmdal - lmdas);

		Usl.d = Ul.d * (lmdal - Wl.Vx) * den;
		Usl.Mx = (Ul.Mx * (lmdal - Wl.Vx) + ps - Wl.P) * den;
		Usl.My = Ul.My * (lmdal - Wl.Vx) * den;
		Usl.Mz = Ul.Mz * (lmdal - Wl.Vx) * den;
		Usl.E = (Ul.E * (lmdal - Wl.Vx) + ps * lmdas - Wl.P * Wl.Vx) * den;

		/* now calcUlate Fsr using Mignone Eq 14 */

		pFlux->d = lmdal * (Usl.d - Ul.d) + Fl.d;
		pFlux->Mx = lmdal * (Usl.Mx - Ul.Mx) + Fl.Mx;
		pFlux->My = lmdal * (Usl.My - Ul.My) + Fl.My;
		pFlux->Mz = lmdal * (Usl.Mz - Ul.Mz) + Fl.Mz;
		pFlux->E = lmdal * (Usl.E - Ul.E) + Fl.E;

		return;
	}
	else if (lmdar >= 0.0) { /* Frs */

		/* Mignone 2006 Eq 48 */
		ps = -Fhll.E * lmdas + Fhll.Mx;

		/* now calcUlate Usr with Mignone Eq 16 */
		den = 1.0 / (lmdar - lmdas);

		Usr.d = Ur.d * (lmdar - Wr.Vx) * den;
		Usr.Mx = (Ur.Mx * (lmdar - Wr.Vx) + ps - Wr.P) * den;
		Usr.My = Ur.My * (lmdar - Wr.Vx) * den;
		Usr.Mz = Ur.Mz * (lmdar - Wr.Vx) * den;
		Usr.E = (Ur.E * (lmdar - Wr.Vx) + ps * lmdas - Wr.P * Wr.Vx) * den;

		/* now calcUlate Fsr using Mignone Eq 14 */

		pFlux->d = lmdar * (Usr.d - Ur.d) + Fr.d;
		pFlux->Mx = lmdar * (Usr.Mx - Ur.Mx) + Fr.Mx;
		pFlux->My = lmdar * (Usr.My - Ur.My) + Fr.My;
		pFlux->Mz = lmdar * (Usr.Mz - Ur.Mz) + Fr.Mz;
		pFlux->E = lmdar * (Usr.E - Ur.E) + Fr.E;

		return;
	}
	else { /* Fr */
		/* intercell flux is right flux */
		pFlux->d = Fr.d;
		pFlux->Mx = Fr.Mx;
		pFlux->My = Fr.My;
		pFlux->Mz = Fr.Mz;
		pFlux->E = Fr.E;

		return;
	}

	/* need to deal with scalar fluxes */
}

#endif /* HYDRO */