#ifndef SHORT_CHARACTERISTICS_HEADER_H
#define SHORT_CHARACTERISTICS_HEADER_H

#define ONLY_HLLC
//#define ONLY_ILLUM

#define WRITE_LOG

//#define USE_VTK

//#define SAVE_DUMP_HLLC

// виды геометрии
//#define Cone
//#define Const_1d
#define Sphere
//#define Ellips


#ifdef USE_VTK
#define ReBuildSolve
//#define ReWriteSolve // записать все массивы из grid в bin
#endif

#include <iomanip>
#include <algorithm>
//#include <ctime>
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
//#include </nethome/student/FS18/FS2-x1/Lartsev/Eigen/Dense>

#include <omp.h>

#endif
