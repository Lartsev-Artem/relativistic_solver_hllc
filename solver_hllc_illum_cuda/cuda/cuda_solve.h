#if !defined CUDA_SOLVE_H && defined USE_CUDA
#define CUDA_SOLVE_H
#include "../solve_module/solve_global_struct.h"

void CudaSendIllumAsync(const int size, const int shift, const Type* Illum_host);
void CopyIllumOnDevice(const int size, const Type* Illum_host);

int CalculateIntScattering(const grid_directions_t& grid_dir, grid_t& grid);
int CalculateIntScatteringAsync(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateEnergy(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateStream(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateImpuls(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateDivImpuls(const grid_directions_t& grid_dir, grid_t& grid);

int CalculateDivStream(const grid_directions_t& grid_dir, grid_t& grid);

void SetDevice(const int num_dev);
void InitDevice(const grid_directions_t& grid_dir_host, grid_t& grid_host);
void ClearDevice();

void ClearHost(grid_t& grid_host);
void CudaWait();

#endif //CUDA_SOLVE_H