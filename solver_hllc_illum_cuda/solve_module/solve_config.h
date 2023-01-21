#ifndef SOLVE_CONFIG_H
#define SOLVE_CONFIG_H

#include "../prj_config.h"
#include "../global_def.h"
#ifdef SOLVE

#include "../global_def.h"
#include "../global_value.h"

//#define HLLC

#define RHLLC

//#define ILLUM

#ifdef ILLUM
//#define SORT_ILLUM  //сорьтровать ли суммирование интегралов по возростанию
#endif

#ifdef USE_CUDA
#include "../solve_short_characteristic_cuda.cuh"
#endif


#endif //SOLVE
#endif //SOLVE_CONFIG_H