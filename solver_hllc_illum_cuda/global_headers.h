#ifndef GLOBAL_HEADERS
#define GLOBAL_HEADERS

#include "prj_config.h"


#include <iomanip>
#include <algorithm>

#include <execution>
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

#if defined USE_MPI
#include "mpi.h"
#endif

#if !defined CLASTER

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

#include <eigen3/Eigen/Dense>
#else

#include </nethome/student/FS18/FS2-x1/Lartsev/Eigen/Dense>

#endif // CLASTER

#if defined BUILD
	int RunBuildModule(const std::string& name_file_settings);
#endif //BUILD

#if defined MAKE

#endif //MAKE

#if defined SOLVE

#endif //SOLVE

#include "global_def.h"

template<typename Str>
int ReadStartSettings(const std::string& name_file_settings, int& class_file_vtk, Str& name_file_vtk, Str& name_file_sphere_direction, 
	Str& graph_adress, Str& base_adress, Str& solve_adress)
{
	std::ifstream ifile;

	OPEN_FSTREAM(ifile, name_file_settings.c_str())
	
	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	getline(ifile, str);

	std::getline(ifile, name_file_vtk);
	std::getline(ifile, name_file_sphere_direction);
	std::getline(ifile, graph_adress);
	std::getline(ifile, base_adress);
	std::getline(ifile, solve_adress);
	
	ifile.close();
	return 0;
}

#endif //GLOBAL_HEADERS
