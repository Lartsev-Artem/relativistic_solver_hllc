#include "../prj_config.h"
#ifdef USE_CUDA
#include "../solve_module/solve_global_struct.h"

#include "cuda_kernel.cuh"
#include "cuda_memory.cuh"

#define CUDA_BLOCKS_2D(val, cells, dirs) dim3 val((cells + BS - 1) / BS, (dirs + BS - 1) / BS);
#define CUDA_TREADS_2D(val)              dim3 threads(BS, BS);

#define CUDA_BLOCKS_1D(val, cells) dim3 val((cells + BS - 1) / BS);
#define CUDA_TREADS_1D(val)              dim3 threads(BS);


void SetDevice(const int num_dev)
{
    CUDA_CALL_FUNC(cudaSetDevice, num_dev);        
}

void CopyIllumOnDevice(const int size, const Type* Illum_host)
{
    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum, Illum_host, size * sizeof(Illum_host[0]));
}

int CalculateIntScattering(const grid_directions_t& grid_dir, grid_t& grid)
{
    const int M = grid_dir.size;
    const int N = grid.size;
    
    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum, grid.Illum, N * M * base * sizeof(grid.Illum[0]));

    CUDA_TREADS_2D(threads);
    CUDA_BLOCKS_2D(blocks, N, M);

    CUDA_CALL_KERNEL(d_GetS, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);
       
    CUDA_MEMCPY_TO_HOST(grid.scattering, device_host_ptr.int_scattering, N * M * sizeof(grid.scattering[0]));
       
    return 0;
}

int CalculateEnergy(const grid_directions_t& grid_dir, grid_t& grid)
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.energy, grid.energy, N * sizeof(grid.energy[0]));

    CUDA_TREADS_1D(threads);
    CUDA_BLOCKS_1D(blocks, N);

    CUDA_CALL_KERNEL(d_MakeEnergy, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.energy, device_host_ptr.energy, N * sizeof(grid.energy[0]));

    return 0;
}

int CalculateStream(const grid_directions_t& grid_dir, grid_t& grid) 
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.stream, grid.stream, N * sizeof(grid.stream[0]));

    CUDA_TREADS_1D(threads);
    CUDA_BLOCKS_1D(blocks, N);
    CUDA_CALL_KERNEL(d_MakeStream, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.stream, device_host_ptr.stream, N * sizeof(grid.stream[0]));

    return 0;
}

int CalculateImpuls(const grid_directions_t& grid_dir, grid_t& grid) 
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.impuls, grid.impuls, N * sizeof(grid.impuls[0]));

    CUDA_TREADS_1D(threads);
    CUDA_BLOCKS_1D(blocks, N);
    CUDA_CALL_KERNEL(d_MakeImpuls, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.impuls, device_host_ptr.impuls, N * sizeof(grid.impuls[0]));

    return 0;    
}

int CalculateDivImpuls(const grid_directions_t& grid_dir, grid_t& grid) 
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.divimpuls, grid.divimpuls, N * sizeof(grid.divimpuls[0]));

    CUDA_TREADS_1D(threads);
    CUDA_BLOCKS_1D(blocks, N);
    CUDA_CALL_KERNEL(d_MakeDivImpuls, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.divimpuls, device_host_ptr.divimpuls, N * sizeof(grid.divimpuls[0]));

#ifdef CUDA_FULL_ARRAYS
    CUDA_MEMCPY_TO_HOST(grid.impuls, device_host_ptr.impuls, N * sizeof(grid.impuls[0]));
#endif

    return 0;    
}

int CalculateDivStream(const grid_directions_t& grid_dir, grid_t& grid) 
{
    const int M = grid_dir.size;
    const int N = grid.size;    
    
    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.divstream, grid.divstream, N * sizeof(grid.divstream[0]));

    CUDA_TREADS_1D(threads);
    CUDA_BLOCKS_1D(blocks, N);
    
    CUDA_CALL_KERNEL(d_MakeDivStream, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.divstream, device_host_ptr.divstream, N * sizeof(grid.divstream[0]));

#ifdef CUDA_FULL_ARRAYS
    CUDA_MEMCPY_TO_HOST(grid.stream, device_host_ptr.stream, N * sizeof(grid.stream[0]));
#endif

    return 0;   
}

#endif //USE_CUDA