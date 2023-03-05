#ifdef USE_CUDA
#include "../solve_module/solve_global_struct.h"

#include "cuda_kernel.cuh"
#include "cuda_memory.cuh"

#define CUDA_BLOCKS_2D(val, cells, dirs) dim3 val((cells + BS - 1) / BS, (dirs + BS - 1) / BS);
#define CUDA_TREADS_2D(val)              dim3 threads(BS, BS);

#define CUDA_BLOCKS_1D(val, cells, bs) dim3 val((cells + bs - 1) / bs);
#define CUDA_TREADS_1D(val, bs)              dim3 threads(bs);

static cudaStream_t cuda_streams[eCuda_count];

static cudaEvent_t cuda_events[eCuda_count*3];
void CudaEventPrint()
{
    for (int i = 0; i < eCuda_count * 3; i++)
    {
        //synchronize
        cudaEventSynchronize(cuda_events[i]); //optional        
    }
   
    float dt_ms;
    cudaEventElapsedTime(&dt_ms, cuda_events[0], cuda_events[1]);
    printf("time kernel1 %f\n", dt_ms);

    cudaEventElapsedTime(&dt_ms, cuda_events[1], cuda_events[2]);
    printf("time send1 %f\n", dt_ms);

    cudaEventElapsedTime(&dt_ms, cuda_events[3], cuda_events[4]);
    printf("time kernel2 %f\n", dt_ms);

    cudaEventElapsedTime(&dt_ms, cuda_events[4], cuda_events[5]);
    printf("time send2 %f\n\n", dt_ms);
}
void CudaSyncStream(const e_cuda_stream_id_t stream_id)
{   
    cudaStreamSynchronize(cuda_streams[stream_id]);
}

void CudaWait()
{
    CUDA_CALL_FUNC(cudaDeviceSynchronize);
}
void SetDevice(const int num_dev)
{
    CUDA_CALL_FUNC(cudaSetDevice, num_dev);  

    //for (int i = 0; i < eCuda_count; i++)
    //{
    //    CUDA_CALL_FUNC(cudaStreamCreate, &cuda_streams[i]);
    //}

    //for (int i = 0; i < eCuda_count*3; i++)
    //{
    //    CUDA_CALL_FUNC(cudaEventCreate, &cuda_events[i]);
    //}

    // get the range of stream priorities for this device
    int priority_high, priority_low;
    cudaDeviceGetStreamPriorityRange(&priority_low, &priority_high);
    
    CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[eCuda_params], cudaStreamNonBlocking, priority_high);
    CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[eCuda_scattering_1], cudaStreamNonBlocking, (priority_high+ priority_low)/2);
    CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[eCuda_scattering_2], cudaStreamNonBlocking, priority_low);
}

void CudaSendIllumAsync(const int size,const int shift ,const Type* Illum_host)
{
    //CUDA_MEMCPY_TO_DEVICE_ASYNC(device_host_ptr.illum + shift, Illum_host + shift, size*sizeof(Illum_host[0]));
    CUDA_CALL_FUNC(cudaMemcpyAsync, device_host_ptr.illum + shift, Illum_host + shift, size * sizeof(Illum_host[0]),
        cudaMemcpyHostToDevice, cuda_streams[0]);

    //CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum + shift, Illum_host + shift, size * sizeof(Illum_host[0]));
}
void CudaSendIllumAsync(const int size, const int shift_dev, const int shift_host, const Type* Illum_host)
{
    //CUDA_MEMCPY_TO_DEVICE_ASYNC(device_host_ptr.illum + shift_dev, Illum_host + shift_host, size * sizeof(Illum_host[0]));
    CUDA_CALL_FUNC(cudaMemcpyAsync, device_host_ptr.illum + shift_dev, Illum_host + shift_host, size * sizeof(Illum_host[0]),
        cudaMemcpyHostToDevice, cuda_streams[0]);
    //CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum + shift, Illum_host + shift, size * sizeof(Illum_host[0]));
}

void CopyIllumOnDevice(const int size, const Type* Illum_host)
{
    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum, Illum_host, size * sizeof(Illum_host[0]));
}
int CalculateIntScatteringAsyncMPIStream(const grid_directions_t& grid_dir, grid_t& grid,
    const int start, const int end, const e_cuda_stream_id_t stream)
{
    const int M = grid_dir.size;
    const int N = grid.size;
  
    //  CUDA_TREADS_2D(threads);
    //  CUDA_BLOCKS_2D(blocks, N, M);

    dim3 threads(32, 16);
    dim3 blocks((N + 32 - 1) / 32, (M + 16 - 1) / 16);

    // надо как то дать задержку второму потоку, относительно первого

    CUDA_CALL_KERNEL_STREAM(d_GetS_MPI_Stream, blocks, threads, cuda_streams[stream], grid_dir_device_ptr, grid_cell_device_ptr, start, end);

    // CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.scattering+start, device_host_ptr.int_scattering + start, N * (end-start) * sizeof(grid.scattering[0]),
        cuda_streams[stream]);    
    
    return 0;
}
int CalculateIntScatteringAsyncMPI(const grid_directions_t& grid_dir, grid_t& grid, int local_size)
{
    const int M = grid_dir.size;
    const int N = grid.size;

    //  CUDA_TREADS_2D(threads);
    //  CUDA_BLOCKS_2D(blocks, N, M);

    dim3 threads(32, 16);
    dim3 blocks((N + 32 - 1) / 32, (M + 16 - 1) / 16);

    CUDA_CALL_KERNEL(d_GetS_MPI, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    // CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST_ASYNC(grid.scattering, device_host_ptr.int_scattering, N * local_size * sizeof(grid.scattering[0]));

    return 0;
}

int CalculateIntScatteringAsync(const grid_directions_t& grid_dir, grid_t& grid)
{
    const int M = grid_dir.size;
    const int N = grid.size;
    
  //  CUDA_TREADS_2D(threads);
  //  CUDA_BLOCKS_2D(blocks, N, M);

    dim3 threads(32, 16);
    dim3 blocks((N + 32 - 1) / 32, (M + 16 - 1) / 16);
    

    CUDA_CALL_KERNEL(d_GetS, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

   // CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST_ASYNC(grid.scattering, device_host_ptr.int_scattering, N * M * sizeof(grid.scattering[0]));

    return 0;
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

    CUDA_TREADS_1D(threads, BS);
    CUDA_BLOCKS_1D(blocks, N, BS);

  //  CUDA_CALL_KERNEL(d_MakeEnergy, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.energy, device_host_ptr.energy, N * sizeof(grid.energy[0]));

    return 0;
}

int CalculateStream(const grid_directions_t& grid_dir, grid_t& grid) 
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.stream, grid.stream, N * sizeof(grid.stream[0]));

    CUDA_TREADS_1D(threads, 512);
    CUDA_BLOCKS_1D(blocks, N, 512);
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

    CUDA_TREADS_1D(threads, 512);
    CUDA_BLOCKS_1D(blocks, N, 512);
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

    CUDA_TREADS_1D(threads, BS);
    CUDA_BLOCKS_1D(blocks, N, BS);
 //   CUDA_CALL_KERNEL(d_MakeDivImpuls, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

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

    CUDA_TREADS_1D(threads, BS);
    CUDA_BLOCKS_1D(blocks, N, BS);
    
   // CUDA_CALL_KERNEL(d_MakeDivStream, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    CUDA_CALL_FUNC(cudaGetLastError);

    CUDA_MEMCPY_TO_HOST(grid.divstream, device_host_ptr.divstream, N * sizeof(grid.divstream[0]));

#ifdef CUDA_FULL_ARRAYS
    CUDA_MEMCPY_TO_HOST(grid.stream, device_host_ptr.stream, N * sizeof(grid.stream[0]));
#endif

    return 0;   
}

int CalculateAllParam(const grid_directions_t& grid_dir, grid_t& grid)
{
    const int M = grid_dir.size;
    const int N = grid.size;
    
    CUDA_TREADS_1D(threads, BS);
    CUDA_BLOCKS_1D(blocks, N, BS);

    CUDA_CALL_KERNEL(d_MakeIllumParam, blocks, threads, grid_dir_device_ptr, grid_cell_device_ptr);

    //CUDA_CALL_FUNC(cudaGetLastError);

    //не асинхронно т.к. сразу затем идет расчет газовый
    CUDA_MEMCPY_TO_HOST(grid.divstream, device_host_ptr.divstream, N * sizeof(grid.divstream[0]));
    CUDA_MEMCPY_TO_HOST(grid.divimpuls, device_host_ptr.divimpuls, N * sizeof(grid.divimpuls[0]));

#ifdef CUDA_FULL_ARRAYS
    CUDA_MEMCPY_TO_HOST_ASYNC(grid.energy, device_host_ptr.energy, N * sizeof(grid.energy[0]));
    CUDA_MEMCPY_TO_HOST_ASYNC(grid.stream, device_host_ptr.stream, N * sizeof(grid.stream[0]));
    CUDA_MEMCPY_TO_HOST_ASYNC(grid.impuls, device_host_ptr.impuls, N * sizeof(grid.impuls[0]));
#endif

    return 0;
}

int CalculateAllParamStream(const grid_directions_t& grid_dir, grid_t& grid, e_cuda_stream_id_t st)
{
    const int M = grid_dir.size;
    const int N = grid.size;

    CUDA_TREADS_1D(threads, BS);
    CUDA_BLOCKS_1D(blocks, N, BS);

    CUDA_CALL_KERNEL_STREAM(d_MakeIllumParam, blocks, threads,cuda_streams[st], grid_dir_device_ptr, grid_cell_device_ptr);

    //CUDA_CALL_FUNC(cudaGetLastError);

    //не асинхронно т.к. сразу затем идет расчет газовый
    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.divstream, device_host_ptr.divstream, N * sizeof(grid.divstream[0]), cuda_streams[st]);
    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.divimpuls, device_host_ptr.divimpuls, N * sizeof(grid.divimpuls[0]), cuda_streams[st]);

#ifdef CUDA_FULL_ARRAYS
    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.energy, device_host_ptr.energy, N * sizeof(grid.energy[0]), cuda_streams[st]);
    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.stream, device_host_ptr.stream, N * sizeof(grid.stream[0]), cuda_streams[st]);
    CUDA_MEMCPY_TO_HOST_ASYNC_STREAM(grid.impuls, device_host_ptr.impuls, N * sizeof(grid.impuls[0]), cuda_streams[st]);
#endif

    return 0;
}
#endif //USE_CUDA