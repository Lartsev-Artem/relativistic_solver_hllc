#ifndef GLOBAL_HEADERS
#define GLOBAL_HEADERS

#include "prj_config.h"


#include <iomanip>
#include <algorithm>

//#include <execution>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <stdlib.h>

#include <chrono>
#include <cstdio>
#include <inttypes.h>
#include <memory>  

#include <stdio.h>

#include <omp.h>

#if !__NVCC__
#if defined USE_MPI
#include "mpi.h"
#endif
#endif

#if !defined CLASTER

#include <sys/stat.h>
#include <direct.h>

#ifdef USE_VTK
#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkDataSet.h>
#include <vtk-9.0\vtkDataSetAttributes.h>
#include <vtk-9.0\vtkDataObject.h>
#include <vtk-9.0\vtkDoubleArray.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkLine.h>
#include <vtk-9.0\vtkLineSource.h>
#include <vtk-9.0\vtkMath.h>
#include <vtk-9.0\vtkNamedColors.h>
#include <vtk-9.0\vtkPointData.h>
#include <vtk-9.0\vtkPoints.h>
#include <vtk-9.0\vtkQuad.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkTetra.h>
#include <vtk-9.0\vtkTriangle.h>
#include <vtk-9.0\vtkUnsignedCharArray.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>
#endif // USE_VTK

//#include <eigen3/Eigen/Dense>
#include "C:/DEV/vcpkg-master/installed/x64-windows/include/eigen3/Eigen/Dense"
#else

#include </nethome/student/FS18/FS2-x1/Lartsev/Eigen/Dense>

#endif // CLASTER

#if defined BUILD
	int RunBuildModule(const std::string& name_file_settings);
#endif //BUILD

#if defined MAKE
	int RunMakeModule(std::string name_file_settings, int a, int b);
#endif //MAKE

#if defined SOLVE
	int RunSolveModule(const std::string& name_file_settings);
#endif //SOLVE

#if defined UTILS
	int RunUtilsModule(int argc, char* argv[], const std::string& settings);
#endif //SOLVE
	

#include "global_def.h"

template<typename Str>
int ReadStartSettings(const std::string& name_file_settings, int& class_file_vtk, Str& name_file_vtk, Str& name_file_sphere_direction, 
	Str& graph_adress, Str& illum_geo_adress, Str& base_adress, Str& solve_adress, int& number_of_iter)
{
	std::ifstream ifile(name_file_settings);
	if (!ifile.is_open()) 
	{
		printf("Error : file %s is not open\n", name_file_settings.c_str());
		return 1;
	}	
	
	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	std::getline(ifile, str);

	std::getline(ifile, name_file_vtk);
	std::getline(ifile, name_file_sphere_direction);
	std::getline(ifile, graph_adress);
	std::getline(ifile, illum_geo_adress);
	std::getline(ifile, base_adress);
	std::getline(ifile, solve_adress);
	ifile >> number_of_iter;
	
	ifile.close();

	return 0;
}


template<typename Str>
int ReadStartSettings(const std::string& name_file_settings, int& class_file_vtk, Str& name_file_vtk, Str& name_file_sphere_direction,
	Str& graph_adress, Str& illum_geo_adress, Str& base_adress, Str& solve_adress, int& number_of_iter, Str& hllc_init_value)
{
	std::ifstream ifile(name_file_settings);
	if (!ifile.is_open())
	{
		printf("Error : file %s is not open\n", name_file_settings.c_str());
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	std::getline(ifile, str);

	std::getline(ifile, name_file_vtk);
	std::getline(ifile, name_file_sphere_direction);
	std::getline(ifile, graph_adress);
	std::getline(ifile, illum_geo_adress);
	std::getline(ifile, base_adress);
	std::getline(ifile, solve_adress);	
	ifile >> number_of_iter; std::getline(ifile, str);
	std::getline(ifile, hllc_init_value);

	ifile.close();
	return 0;
}

#endif //GLOBAL_HEADERS
