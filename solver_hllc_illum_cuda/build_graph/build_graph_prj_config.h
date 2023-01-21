#ifndef BUILD_GRAPH_PRJ_CONFIG
#define BUILD_GRAPH_PRJ_CONFIG

#include "../prj_config.h"

#if !defined CLASTER && defined USE_VTK
#define WriteFiles  // ��������� ������ vtk � �������� ������ ����������� ��� ��������� ������
#endif

#ifndef ONLY_GEO_DATA
#define ReadFiles     // ������ �������������� ������ ��� ���������� ������ � �� ���������� ������
#endif

#define FastWriteFile

//#define USE_OMP        // ����������� ���������� omp (�������������� / ������ mpi) 

#define GRID_WITH_INNER_BOUNDARY  // ���� ��� ����� � ���������� �������� (��� ����������� ������� ���������, ��� ������������� �������� �����)

#if NUMBER_OF_MEASUREMENTS == 3
#define TASK_3D  // ���������� �����
#include <set>
#elif NUMBER_OF_MEASUREMENTS == 2
#define TASK_2D  // ����������� ������. ������ �������, �������, �����, ����� � �������� (��� hllc) + ������
#endif

#define USE_STRANGE_FUNCTION // ����� �� ������ �������, ������� � ������� ������� �� ������������

// DEBUG
//#define ONLY_ONE_DIRECTION  // ������  � ���������������� ������ 

#endif