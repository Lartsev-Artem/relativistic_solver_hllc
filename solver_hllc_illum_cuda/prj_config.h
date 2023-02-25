#ifndef PRJ_CONFIG
#define PRJ_CONFIG

//#define CLASTER	// ������ ��� �������

#ifndef CLASTER

#if !__NVCC__
#define USE_VTK	  // ������������� vtk ��� ������ ����������� � ���� �����
#endif

//#define UTILS

//#define WRITE_LOG_ON_SCREAN //������ ��� �� � ����, � �� �����

#endif

#define NUMBER_OF_MEASUREMENTS 3

//#define BUILD  // ������ ������ ����������� ������

//#define MAKE   // ������ ������ ���������������� �������� ��� ���������

#define SOLVE  // ������ ������ �������

#if !__NVCC__
#define USE_MPI  // ����������� ���������� mpi
#endif

//#define USE_CUDA  // ����������� ���������� cuda

#define WRITE_GLOBAL_LOG	// ������ ��� ����

#if defined BUILD || defined MAKE  //todo rename on BUILD_DATA_TO_ILLUM
//#define ONLY_GEO_DATA
#endif

#if defined UTILS && !defined SOLVE

#define ILLUM
#define HLLC

#endif //UTILS


// ���������

//#define Cube

//#define Step

#define Cone

#define Cone_JET

//#define Sphere

//#define Cylinder

#endif //PRJ_CONFIG
