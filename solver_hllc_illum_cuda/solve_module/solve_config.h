#ifndef SOLVE_CONFIG_H
#define SOLVE_CONFIG_H

#include "../prj_config.h"
#include "../global_def.h"
#include "../global_value.h"

#ifdef SOLVE

//#define HLLC

//#define RHLLC

#if NUMBER_OF_MEASUREMENTS == 3 // излучение доступно только для 3d

#define ILLUM

#ifdef ILLUM
//#define SORT_ILLUM  //сорьтровать ли суммирование интегралов по возростанию
#endif

#endif //3d

#ifdef USE_CUDA
#include "../solve_short_characteristic_cuda.cuh"
#endif

#ifdef USE_MPI

#ifdef ILLUM
#define ILLUM_MPI
#endif

#ifdef HLLC
#define HLLC_MPI
#endif

#ifdef RHLLC
//#define RHLLC_MPI
#endif

#endif

#ifndef CLASTER
//#define RUN_TEST
#endif

#endif //SOLVE
#endif //SOLVE_CONFIG_H