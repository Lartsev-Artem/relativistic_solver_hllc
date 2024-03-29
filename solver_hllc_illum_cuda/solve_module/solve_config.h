#if !defined SOLVE_CONFIG_H && defined SOLVE
#define SOLVE_CONFIG_H

#include "../prj_config.h"
#include "../global_def.h"
#include "../global_value.h"

#ifdef SOLVE

//#define HLLC

#define RHLLC

#if NUMBER_OF_MEASUREMENTS == 3 // ��������� �������� ������ ��� 3d

#define ILLUM

#ifdef ILLUM
//#define SORT_ILLUM  //����������� �� ������������ ���������� �� �����������
#endif

#endif //3d

#ifdef USE_MPI

#include "mpi.h"

#ifdef ILLUM
#define ILLUM_MPI
#endif

#ifdef HLLC
#define HLLC_MPI
#endif

#ifdef RHLLC
#define RHLLC_MPI
#endif

#endif

#ifndef CLASTER
//#define RUN_TEST
#endif

#if defined USE_CUDA && !defined ILLUM
#error "Bad config. CUDA only with ILLUM!!!"
#endif


#define ON_FULL_ILLUM_ARRAYS

#endif //SOLVE
#endif //SOLVE_CONFIG_H