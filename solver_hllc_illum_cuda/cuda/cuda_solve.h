#if !defined CUDA_SOLVE_H && defined USE_CUDA
#define CUDA_SOLVE_H
#include "../solve_module/solve_global_struct.h"

void CudaEventPrint();
void CudaSyncStream(const e_cuda_stream_id_t stream_id);
int CalculateIntScatteringAsyncMPIStream(const grid_directions_t& grid_dir, grid_t& grid,
    const int start, const int end, const e_cuda_stream_id_t stream);

void CudaSendIllumAsync(const int size, const int shift, const Type* Illum_host);
void CudaSendIllumAsync(const int size, const int shift_dev, const int shift_host, const Type* Illum_host);
int CalculateIntScatteringAsyncMPI(const grid_directions_t& grid_dir, grid_t& grid, int local_size);
void CopyIllumOnDevice(const int size, const Type* Illum_host);

int CalculateIntScattering(const grid_directions_t& grid_dir, grid_t& grid);
int CalculateIntScatteringAsync(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateEnergy(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateStream(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateImpuls(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateDivImpuls(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateDivStream(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateAllParam(const grid_directions_t& grid_dir, grid_t& grid);
int CalculateAllParamStream(const grid_directions_t& grid_dir, grid_t& grid, e_cuda_stream_id_t st);


void SetDevice(const int num_dev);
void InitDevice(const grid_directions_t& grid_dir_host, grid_t& grid_host, const int start, const int end);
void ClearDevice();

void ClearHost(grid_t& grid_host);
void CudaWait();

#endif //CUDA_SOLVE_H