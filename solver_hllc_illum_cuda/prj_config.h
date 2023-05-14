#ifndef PRJ_CONFIG
#define PRJ_CONFIG

//#define CLASTER	// сборка под кластер

#ifndef CLASTER

#if !__NVCC__
#define USE_VTK	  // использование vtk для вывода результатов в виде сетки
#endif

//#define UTILS

//#define WRITE_LOG_ON_SCREAN //писать лог не в файл, а на экран

#endif

#define NUMBER_OF_MEASUREMENTS 3

//#define BUILD  // запуск модуля построителя графов

//#define MAKE   // запуск модуля подготовительных расчётов для излучения

//#define SOLVE  // запуск модуля решения

//#define USE_CUDA  // подключение технологии cuda


#if !__NVCC__
#define USE_MPI  // подключение технологии mpi
#endif


#define WRITE_GLOBAL_LOG	// писать лог файл

#ifdef WRITE_GLOBAL_LOG
#define WRITE_MPI_LOG	// писать mpi лог файл
#endif

#if defined BUILD || defined MAKE  //todo rename on BUILD_DATA_TO_ILLUM
#define ONLY_GEO_DATA
#endif

#if defined UTILS && !defined SOLVE

#define ILLUM
#define HLLC

#endif //UTILS

// Геометрия

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
