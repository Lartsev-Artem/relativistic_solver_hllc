#if !defined CUDA_STRUCT_H && defined USE_CUDA
#define CUDA_STRUCT_H

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "cuda_def.cuh"

#ifdef CUDA_DEBUG_MODE
#define CHECK_SCALE(i, size) if (i >= size)  CUDA_EXIT_ERR(cuda: Error CHECK_SCALE  i >= size)
#else
#define CHECK_SCALE(i, size) {}
#endif

template<class type_t, int size>
struct cuda_vector_t
{
    type_t data[size];

    __host__ __device__  cuda_vector_t() : data{ 0 } {}

    __host__ __device__ Type operator [](const int i)  const
    {
        CHECK_SCALE(i, size);
        return data[i];
    }

    __host__ __device__ Type& operator [](const int i)
    {
        CHECK_SCALE(i, size);
        return *(data + i);
    }
};

struct direction_device_t
{
    cuda_vector_t<double, 3> dir;
    double area;
};
struct grid_directions_device_t
{
    int size;
    direction_device_t* directions;
    double full_area;    

    __host__ __device__  grid_directions_device_t() : size(0), directions(nullptr), full_area(0) {}
    __host__ __device__ ~grid_directions_device_t() { }
};

struct grid_device_t
{
    int size;
    int local_scattering_size;
    int local_scattering_disp;
    Type* illum;
    Type* int_scattering;

    Type* divstream;
    cuda_vector_t<Type, 3>* divimpuls;

    cuda_vector_t<Type, 3>* normals;
    Type* areas;
    Type* volume;

#ifdef CUDA_FULL_ARRAYS
    Type* energy;
    cuda_vector_t<Type, 3>* stream;
    cuda_vector_t<Type, 9>* impuls;
#endif

    grid_device_t() : size(0), illum(nullptr), int_scattering(nullptr), divstream(nullptr), divimpuls(nullptr),
        normals(nullptr), areas(nullptr), volume(nullptr)
#ifdef CUDA_FULL_ARRAYS       
        , energy(nullptr), stream(nullptr), impuls(nullptr) {}
#endif

};


struct device_host_ptr_t
{
    direction_device_t* directions;

    Type* illum;
    Type* int_scattering;

    Type* divstream;
    cuda_vector_t<Type, 3>* divimpuls;

    cuda_vector_t<Type, 3>* normals;
    Type* areas;
    Type* volume;

#ifdef CUDA_FULL_ARRAYS
    Type* energy;
    cuda_vector_t<Type, 3>* stream;
    cuda_vector_t<Type, 9>* impuls;
#endif
};
extern device_host_ptr_t device_host_ptr;
extern grid_directions_device_t* grid_dir_device_ptr;
extern grid_device_t* grid_cell_device_ptr;

#endif //CUDA_STRUCT_H