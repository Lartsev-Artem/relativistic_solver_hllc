#ifndef PRJ_CONFIG
#define PRJ_CONFIG

//#define CLASTER	// сборка под кластер

#ifndef CLASTER

#define USE_VTK	  // использование vtk для вывода результатов в виде сетки

//#define UTILS

#endif

#define NUMBER_OF_MEASUREMENTS 3

//#define BUILD  // запуск модуля построителя графов

//#define MAKE   // запуск модуля подготовительных расчётов для излучения

#define SOLVE  // запуск модуля решения

//#define USE_MPI  // подключение технологии mpi

#define USE_CUDA  // подключение технологии cuda

#define WRITE_GLOBAL_LOG	// писать лог файл

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

#endif //PRJ_CONFIG
