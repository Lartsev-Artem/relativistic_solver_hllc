#ifndef BUILD_GRAPH_PRJ_CONFIG
#define BUILD_GRAPH_PRJ_CONFIG

#include "../prj_config.h"

#if !defined CLASTER && defined USE_VTK
#define WriteFiles  // включение файлов vtk и создание файлов необходимых для постоения графов
#endif

#ifndef ONLY_GEO_DATA
#define ReadFiles     // чтение сформированных файлов для построения графов и из дальнейший расчёт
#endif

#define FastWriteFile

//#define USE_OMP        // подключение технологии omp (самостоятельно / вместе mpi) 

#define GRID_WITH_INNER_BOUNDARY  // граф для сетки с внутренней границей (для оптимизации следует отключить, при использование сплошной сетки)

#if NUMBER_OF_MEASUREMENTS == 3
#define TASK_3D  // трехмерная сетка
#include <set>
#elif NUMBER_OF_MEASUREMENTS == 2
#define TASK_2D  // сокращенная версия. Толкьо нормали, площадь, объем, связи с соседями (для hllc) + центры
#endif

#define USE_STRANGE_FUNCTION // какие то старые функции, которые в текущем конфиге не используются

// DEBUG
//#define ONLY_ONE_DIRECTION  // только  в последовательном режиме 

#endif