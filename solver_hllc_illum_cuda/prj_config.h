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

//#define SOLVE  // ������ ������ �������

//#define USE_CUDA  // ����������� ���������� cuda


#if !__NVCC__
#define USE_MPI  // ����������� ���������� mpi
#endif


#define WRITE_GLOBAL_LOG	// ������ ��� ����

#ifdef WRITE_GLOBAL_LOG
#define WRITE_MPI_LOG	// ������ mpi ��� ����
#endif

#if defined BUILD || defined MAKE  //todo rename on BUILD_DATA_TO_ILLUM
#define ONLY_GEO_DATA
#endif

#if defined UTILS && !defined SOLVE

#define ILLUM
#define HLLC

#endif //UTILS

// ���������

//#define Cube

//#define Step

#define Cone

//#define Cone_JET

//#define Sphere

//#define Cylinder

#ifdef DEBUG
//#define DEBUG_MPI_RHLLC

#endif

#endif //PRJ_CONFIG
