#if !defined CUDA_UTILS_H && defined USE_CUDA
#define CUDA_UTILS_H

#include "cuda_struct.cuh"

__device__ Type Gamma(const cuda_vector_t<Type, 3>& direction, const cuda_vector_t<Type, 3>& direction2);

__device__ Type IntegarteDirection(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid);

__device__ cuda_vector_t<Type, 3> IntegarteDirection3(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid);

__device__ cuda_vector_t<Type, 9>  IntegarteDirection9(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid);

__device__ void IntegarteDirection3Faces(const int num_cell, const grid_directions_device_t* dir_grid, grid_device_t* grid, cuda_vector_t<Type, 3>* Stream);

__device__ void IntegarteDirection9Faces(const int num_cell, const grid_directions_device_t* dir_grid, grid_device_t* grid, cuda_vector_t<Type, 9>* Impuls);

#endif //CUDA_UTILS_H