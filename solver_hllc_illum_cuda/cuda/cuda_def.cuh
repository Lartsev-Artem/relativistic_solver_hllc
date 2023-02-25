#ifndef CUDA_DEF_H
#define CUDA_DEF_H

#include "../prj_config.h"
#ifdef USE_CUDA

#include <string>
#include "../global_def.h"
#include "../solve_module/solve_config.h"

#define BS 32
#define CUDA_DEBUG_MODE //отладка в структурах

//#define CUDA_EXIT(a){printf("\nCalls from: ") ;printf(__FUNCTION__); EXIT(a)} //если не сработает mpi, то заменить на exit(a)
//#define CUDA_EXIT_ERR(...) {printf(CONVERT_TO_STRING(__VA_ARGS__)); CUDA_EXIT(1);}

#define CUDA_EXIT(a){printf("\nCalls from: %s", __FUNCTION__); EXIT(a)} //если не сработает mpi, то заменить на exit(a)
#define CUDA_EXIT_ERR(...) {printf(CONVERT_TO_STRING(__VA_ARGS__)); CUDA_EXIT(1);}


#define CUDA_CALL_FUNC(_func, ...) \
if (CheckError(_func(__VA_ARGS__))){ \
CUDA_EXIT_ERR(Error in function: _func \nargs:  __VA_ARGS__ ); }

#define CUDA_CALL_KERNEL(_func, _blocks, _threads, ...) _func << <_blocks, _threads >> > (__VA_ARGS__) //\
CUDA_CALL_FUNC(cudaGetLastError)


#ifdef ON_FULL_ILLUM_ARRAYS
    #define CUDA_FULL_ARRAYS 
#endif

template<typename Err_t>
int CheckError(Err_t cudaStatus, const std::string my_text = "")
{
    if (cudaStatus != cudaSuccess)
    {
        printf("Error cuda code: %d \n%s \n%s\n", cudaStatus, cudaGetErrorString(cudaStatus), my_text.c_str());
        return 1;
    }
    return 0;
}

#endif //USE_CUDA
#endif //CUDA_DEF_H