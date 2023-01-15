#ifndef PRJ_CONFIG
#define PRJ_CONFIG

//#define CLASTER	// сборка под кластер

#ifndef CLASTER

#define USE_VTK	  // использование vtk для вывода результатов в виде сетки

//#define ONLY_BUILD_DATA

#endif

#define NUMBER_OF_MEASUREMENTS 3

#define BUILD  // запуск модуля построителя графов

#define MAKE   // запуск модуля подготовительных расчётов для излучения

#define SOLVE  // запуск модуля решения

#define USE_MPI  // подключение технологии mpi

#define WRITE_GLOBAL_LOG	// писать лог файл


// Геометрия
// 
//#define Cube

//#define Step

//#define Cone

//#define Sphere

#define Cylinder


#endif //PRJ_CONFIG
