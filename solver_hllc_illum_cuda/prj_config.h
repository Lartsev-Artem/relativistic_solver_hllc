#ifndef PRJ_CONFIG
#define PRJ_CONFIG

//#define CLASTER	// ������ ��� �������

#ifndef CLASTER

#define USE_VTK	  // ������������� vtk ��� ������ ����������� � ���� �����

//#define UTILS

#endif

#define NUMBER_OF_MEASUREMENTS 3

#define BUILD  // ������ ������ ����������� ������

//#define MAKE   // ������ ������ ���������������� �������� ��� ���������

#define SOLVE  // ������ ������ �������

//#define USE_MPI  // ����������� ���������� mpi

#define USE_CUDA  // ����������� ���������� cuda

#define WRITE_GLOBAL_LOG	// ������ ��� ����

#if defined BUILD || defined MAKE
#define ONLY_GEO_DATA
#endif

#ifdef UTILS

#define ILLUM
#define HLLC

#endif //UTILS


// ���������

#define Cube

//#define Step

//#define Cone

//#define Sphere

//#define Cylinder

#endif //PRJ_CONFIG
