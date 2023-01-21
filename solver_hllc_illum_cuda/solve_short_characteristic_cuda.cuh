#ifndef SHORT_CHARACTERISTICS_CUDA_H
#define SHORT_CHARACTERISTICS_CUDA_H

#include "prj_config.h"
#ifdef USE_CUDA
#include "global_def.h"

#include "solve_module/solve_config.h"


//************Global Value**********************
extern Type* dev_directions;
extern Type* dev_squares;
extern Type* dev_illum;
extern Type* dev_int_scattering;

extern Type* dev_energy;
extern Type* dev_stream;
extern Type* dev_impuls;

//**********************************************

#define BS 32
#define CUDA_RETURN(s) {printf(s); return 1;}


int CheckDevice();

int InitDevice(const int num_dir, const int num_cell, const int mod = 0);
int HostToDevice(const grid_directions_t& host_directions, std::vector<Type>& host_illum, const int mod = 0);

int CalculateIntScattering(const int b_size, const int N, const int M, std::vector<Type>& host_illum, std::vector<Type>& host_int_scattering);
int CalculateEnergy(const int b_size, const int N, const int M, std::vector<Type>& host_energy);
int CalculateStream(const int b_size, const int N, const int M, std::vector<Vector3>& host_stream);
int CalculateImpuls(const int b_size, const int N, const int M, std::vector<Matrix3>& host_impuls);

int ClearDevice(const int mod = 0);

#endif //CUDA
#endif