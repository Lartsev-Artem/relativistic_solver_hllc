#if !defined CUDA_MEMORY_H && defined USE_CUDA
#define CUDA_MEMORY_H
#include "cuda_def.cuh"

#define CUDA_MEMCPY_TO_DEVICE(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST(dist, src, size) CUDA_CALL_FUNC(cudaMemcpy, dist, src, size, cudaMemcpyDeviceToHost);


#define CUDA_MEMCPY_TO_DEVICE_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyHostToDevice);
#define CUDA_MEMCPY_TO_HOST_ASYNC(dist, src, size) CUDA_CALL_FUNC(cudaMemcpyAsync, dist, src, size, cudaMemcpyDeviceToHost);

#define CUDA_MALLOC(src, size) CUDA_CALL_FUNC(cudaMalloc, (void**)src, size);
#define CUDA_MALLOC_HOST(src, size) CUDA_CALL_FUNC(cudaMallocHost, (void**)src, size);


#define CUDA_FREE_MEMORY(val)  CUDA_CALL_FUNC(cudaFree, val)
#define CUDA_FREE_HOST_MEMORY(val)  CUDA_CALL_FUNC(cudaFreeHost, val)

#endif //CUDA_MEMORY_H