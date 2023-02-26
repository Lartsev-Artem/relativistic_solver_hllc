#if !defined CUDA_KERNEL_H && USE_CUDA
#define CUDA_KERNEL_H
#include "cuda_struct.cuh"

__global__ void d_GetS(const grid_directions_device_t* dir, grid_device_t* grid);

__global__ void d_MakeEnergy(const grid_directions_device_t* dir, grid_device_t* grid);

__global__ void d_MakeStream(const grid_directions_device_t* dir, grid_device_t* grid);

__global__ void d_MakeImpuls(const grid_directions_device_t* dir, grid_device_t* grid);

__global__ void d_MakeDivStream(const grid_directions_device_t* dir, grid_device_t* grid);

__global__ void d_MakeDivImpuls(const grid_directions_device_t* dir, grid_device_t* grid);

#endif //CUDA_KERNEL_H