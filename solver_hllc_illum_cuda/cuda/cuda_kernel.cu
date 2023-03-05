#ifdef USE_CUDA
//***********************************************************************//
//********************Functions calls kernel*****************************//
//***********************************************************************//

#include "cuda_utils.cuh"
#include "cuda_struct.cuh"

#ifdef CUDA_FULL_ARRAYS
#define CUDA_CONVERT_FATE_TO_CELL(val, size, src) \
for (int k = 0; k < size; k++) \
{ \
    val[k] = 0; \
    for (int f = 0; f < base; f++) \
        val[k] += src[f][k]; \
    val[k] /= base; \
}
#else
    CUDA_CONVERT_FATE_TO_CELL(val, size, src) {}
#endif


__global__ void d_GetS(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N) return;
    if (k >= M) return;

    grid->int_scattering[k * N + i] = 0;


    for (int num_direction = 0; num_direction < M; num_direction++)
    {
        Type I = (
            grid->illum[4 * num_direction * N + 4 * i] +
            grid->illum[4 * num_direction * N + 4 * i + 1] +
            grid->illum[4 * num_direction * N + 4 * i + 2] +
            grid->illum[4 * num_direction * N + 4 * i + 3]) / 4;

        grid->int_scattering[k * N + i] += Gamma(dir->directions[num_direction].dir, dir->directions[k].dir) *
            I * dir->directions[num_direction].area;
    }

    grid->int_scattering[k * N + i] /= dir->full_area;
}

__global__ void d_GetS_MPI(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N) return;
    //if (k >= M) return;    
    if (k >= grid->local_scattering_size) return;

    grid->int_scattering[k * N + i] = 0;


    for (int num_direction = 0; num_direction < M; num_direction++)
    {
        Type I = (
            grid->illum[4 * num_direction * N + 4 * i] +
            grid->illum[4 * num_direction * N + 4 * i + 1] +
            grid->illum[4 * num_direction * N + 4 * i + 2] +
            grid->illum[4 * num_direction * N + 4 * i + 3]) / 4;

        grid->int_scattering[k * N + i] += Gamma(dir->directions[num_direction].dir, dir->directions[grid->local_scattering_disp + k].dir) *
            I * dir->directions[num_direction].area;
    }

    grid->int_scattering[k * N + i] /= dir->full_area;
}

__global__ void d_GetS_MPI_Stream(const grid_directions_device_t* dir, grid_device_t* grid, const int start, const int end)
{
    const int M = dir->size;
    const int N = grid->size;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N) return;
    //if (k >= M) return;    
    if (k >= end/*grid->local_scattering_size*/) return;
    if (k < start) return;

    grid->int_scattering[k * N + i] = 0;


    for (int num_direction = 0; num_direction < M; num_direction++)
    {
        Type I = (
            grid->illum[4 * num_direction * N + 4 * i] +
            grid->illum[4 * num_direction * N + 4 * i + 1] +
            grid->illum[4 * num_direction * N + 4 * i + 2] +
            grid->illum[4 * num_direction * N + 4 * i + 3]) / 4;

        grid->int_scattering[k * N + i] += Gamma(dir->directions[num_direction].dir, dir->directions[grid->local_scattering_disp + k].dir) *
            I * dir->directions[num_direction].area;
    }

    grid->int_scattering[k * N + i] /= dir->full_area;
}

__device__ void d_MakeEnergy(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int N = grid->size;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    grid->energy[i] = IntegarteDirection(i, dir, grid);
}

__global__ void d_MakeStream(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    grid->stream[i] = IntegarteDirection3(i, dir, grid);

    return;
}

__global__ void d_MakeImpuls(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    grid->impuls[i] = IntegarteDirection9(i, dir, grid);

    return;
}

__device__ void d_MakeDivStream(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    cuda_vector_t<Type, 3> Stream[base];
    IntegarteDirection3Faces(i, dir, grid, Stream);

    CUDA_CONVERT_FATE_TO_CELL(grid->stream[i], 3, Stream);

    //Type& d_divstream = grid->divstream[i];
    
    grid->divstream[i] = 0;
    for (size_t f = 0; f < base; f++)
    {
        Type sum = 0;
        for (size_t k = 0; k < 3; k++)
        {
            sum += Stream[f][k] * grid->normals[i * base + f][k];
        } 
        grid->divstream[i] += sum * grid->areas[i * base + f];
    }

    grid->divstream[i] /= grid->volume[i];
    return;
}

__device__ void d_MakeDivImpuls(const grid_directions_device_t* dir, grid_device_t* grid)
{
    const int M = dir->size;
    const int N = grid->size;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    cuda_vector_t<Type, 9> impuls[base];

    IntegarteDirection9Faces(i, dir, grid, impuls);

    CUDA_CONVERT_FATE_TO_CELL(grid->impuls[i], 9, impuls);

    //cuda_vector_t<Type, 3>& divimpuls = grid->divimpuls[i];

    for (size_t j = 0; j < 3; j++)
        grid->divimpuls[i][j] = 0;

    for (size_t j = 0; j < base; j++)
    {
        for (size_t h = 0; h < 3; h++)
        {
            Type sum = 0;
            for (size_t k = 0; k < 3; k++)
            {
                sum += impuls[j][h * 3 + k] * grid->normals[i * base + j][k];
            }

            grid->divimpuls[i][h] += sum * grid->areas[i * base + j];
        }
    }

    for (size_t j = 0; j < 3; j++)
        grid->divimpuls[i][j] /= grid->volume[i];

    return;
}


__global__ void d_MakeIllumParam(const grid_directions_device_t* dir, grid_device_t* grid)
{
    // эти функции можно объденить в одну. Тогда будет одно общее обращение в память к illum
    d_MakeEnergy(dir, grid);
    d_MakeDivStream(dir, grid);
    d_MakeDivImpuls(dir, grid);
}

#endif //USE_CUDA