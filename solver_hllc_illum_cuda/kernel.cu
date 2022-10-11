#include "solve_short_characteristic_cuda.cuh"

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

struct dev_Vector3 {
    Type data[3];

    __host__ __device__  dev_Vector3() {
        data[0] = 0;
        data[1] = 0;
        data[2] = 0;
    }

    __host__ __device__ Type operator [](const int i)  const {
        if (i >= 3) {
            printf("Error dev_Vector3\n");
            return 0;
        }

        return data[i];
    }

    __host__ __device__ Type& operator [](const int i) {
        if (i >= 3) {
            printf("Error dev_Vector3\n");
            return *data;
        }
        return *(data + i);
    }
};

struct dev_Matrix3 {
    Type data[9];

    __host__ __device__ dev_Matrix3() {
        data[0] = 0; data[1] = 0; data[2] = 0;
        data[3] = 0; data[4] = 0; data[5] = 0;
        data[6] = 0; data[7] = 0; data[8] = 0;
    }

    __host__ __device__ Type operator [](const int i)  const {
        if (i >= 9) {
            printf("Error dev_Matrix3\n");
            return 0;
        }

        return data[i];
    }

    __host__ __device__  Type& operator [](const int i) {
        if (i >= 9) {
            printf("Error dev_Matrix3\n");
            return *data;
        }
        return *(data + i);
    }
};




//************Global Value**********************
Type* dev_directions;
Type* dev_squares;
Type* dev_illum;
Type* dev_int_scattering;

Type* dev_energy;
Type* dev_stream;
Type* dev_impuls;

//**********************************************


int CheckError(cudaError_t cudaStatus, const char* my_text = "") {
    if (cudaStatus != cudaSuccess) {
        printf("Error code: %d \n%s \n%s\n", cudaStatus, cudaGetErrorString(cudaStatus), my_text);
        return 1;
    }
    return 0;
}


__device__ Type Gamma(const dev_Vector3& direction, const dev_Vector3& direction2) {

    Type sum = 0;
    for (size_t i = 0; i < 3; i++)
    {
        sum += direction[i] * direction2[i];
    }

    return (3. * (1 + sum * sum)) / 4;
}
__device__ Type IntegarteDirection(const int N, const int M, int num_cell,
    Type* illum, Type* squares, const Type square_surface) {

    Type res = 0;

    for (size_t i = 0; i < M; i++) {
        res += illum[N * i + num_cell] * squares[i];
    }

    return res / square_surface;
}
__device__ dev_Vector3 IntegarteDirection3(const int N, const int M, int num_cell,
    Type* illum, const dev_Vector3* directions, Type* squares, const Type square_surface) {

    dev_Vector3 res;

    for (size_t i = 0; i < M; i++) {
        res[0] += directions[i][0] * illum[N * i + num_cell] * squares[i];
        res[1] += directions[i][1] * illum[N * i + num_cell] * squares[i];
        res[2] += directions[i][2] * illum[N * i + num_cell] * squares[i];
    }

    res[0] /= square_surface;
    res[1] /= square_surface;
    res[2] /= square_surface;

    return res;
}
__device__ dev_Matrix3 IntegarteDirection9(const int N, const int M, int num_cell,
    Type* illum, const dev_Vector3* directions, Type* squares, const Type square_surface) {

    dev_Matrix3 res;

    for (size_t i = 0; i < 3; i++)
        for (size_t k = 0; k < 3; k++)

            for (size_t j = 0; j < M; j++) {
                res[i*3 +k] += directions[j][i] * directions[j][k] * (illum[N * j + num_cell] * squares[j]);
            }


    for (size_t i = 0; i < 9; i++)
        res[i] /= square_surface;
    
    return res;
}


__global__ void d_GetS(const int N, const int M, Type* illum_old, Type* Integ,
    const Type* d_directions, Type* squares, const Type square_surface) {

    const dev_Vector3* directions = (dev_Vector3*)d_directions;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= N) return;
    if (k >= M) return;
    Integ[k * N + i] = 0;


    for (int num_direction = 0; num_direction < M; num_direction++)
    {
        //illum_old[i * M + num_direction]
        Type I = (illum_old[4*num_direction * N + 4*i] + illum_old[4*num_direction * N + 4 * i + 1] +
            illum_old[4 * num_direction * N + 4 * i + 2] + illum_old[4 * num_direction * N + 4 * i + 3]) / 4;

        Integ[k * N + i] += Gamma(directions[num_direction], directions[k]) *
            I * squares[num_direction];
    }

    Integ[k * N + i] /= square_surface;
}

__global__ void d_MakeEnergy(const int N, const int M, Type* illum, Type* energy, Type* squares, const Type square_surface) {


    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    //for (size_t i = 0; i < n; ++i)

    energy[i] = IntegarteDirection(N, M, i, illum, squares, square_surface);

    return;
}

__global__ void d_MakeStream(const int N, const int M, Type* illum, Type* d_stream, Type* d_directions, Type* squares, const Type square_surface) {


    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    const dev_Vector3* directions = (dev_Vector3*)d_directions;
    dev_Vector3* stream = (dev_Vector3*)d_stream;

    //for (size_t i = 0; i < n; ++i)

    stream[i] = IntegarteDirection3(N, M, i, illum, directions, squares, square_surface);

    return;
}

__global__ void d_MakeImpuls(const int N, const int M, Type* illum, Type* d_impuls, Type* d_directions, Type* squares, const Type square_surface) {


    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    const dev_Vector3* directions = (dev_Vector3*)d_directions;
    dev_Matrix3* impuls = (dev_Matrix3*)d_impuls;

    //for (size_t i = 0; i < n; ++i)

    impuls[i] = IntegarteDirection9(N, M, i, illum, directions, squares, square_surface);

    return;
}



//***********************************************************************//
//*********************Functions from host*******************************//
//***********************************************************************//

int CheckDevice() {
    //проверка наличия карты 
    if (CheckError(cudaSetDevice(0), "cudaSetDevice failed!Do you have a CUDA - capable GPU installed ?")) return 1;
    return 0;
}

int InitDevice(const int num_dir, const int num_cell, const int mod) {

    if (CheckError(cudaMalloc(&dev_directions, num_dir * 3 * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_squares, num_dir * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_illum, 4 * num_dir * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_int_scattering, num_dir * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;

    if (mod) //extended
    {
        if (CheckError(cudaMalloc(&dev_energy, num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_stream, 3 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_impuls, 9 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
    }

    return 0;
}

int HostToDevice(std::vector<Vector3>& host_directions, std::vector<Type>& host_squares, std::vector<Type>& host_illum, const int mod) {
    
    if (CheckError(cudaMemcpy(dev_illum, host_illum.data(), host_illum.size() * sizeof(host_illum[0]), 
        cudaMemcpyHostToDevice), "cudaMemcpy failed illum!")) return 1;

    if (mod)  // first init
    {
        std::vector<dev_Vector3> dev_dir(host_directions.size());  // может можно сразу?
        for (size_t i = 0; i < host_directions.size(); i++)
        {
            dev_dir[i].data[0] = host_directions[i][0];
            dev_dir[i].data[1] = host_directions[i][1];
            dev_dir[i].data[2] = host_directions[i][2];
        }


        if (CheckError(cudaMemcpy( dev_directions, dev_dir.data(), dev_dir.size() * sizeof(dev_Vector3), cudaMemcpyHostToDevice),
            "cudaMemcpy failed dir!")) return 1;
  
        if (CheckError(cudaMemcpy(dev_squares, host_squares.data(), host_squares.size() * sizeof(host_squares[0]), cudaMemcpyHostToDevice), 
            "cudaMemcpy failed! squares")) return 1;
    }

    return 0;
}

int CalculateIntScattering(const int b_size, const int N, const int M, std::vector<Type>& host_illum, std::vector<Type>& host_int_scattering) {
    
    if (CheckError(cudaMemcpy(dev_illum, host_illum.data(), host_illum.size() *sizeof(host_illum[0]), cudaMemcpyHostToDevice), 
        "cudaMemcpy failed!, Illum in scattering")) return 1;

    dim3 threads(b_size, b_size);
    dim3 blocks((N + b_size - 1) / b_size, (M + b_size - 1) / b_size);

    d_GetS << <blocks, threads >> > (N, M, dev_illum, dev_int_scattering, dev_directions, dev_squares, square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_GetS failed!")) return 1;

    if (CheckError(cudaMemcpy(host_int_scattering.data(), dev_int_scattering, host_int_scattering.size() * sizeof(host_int_scattering[0]), 
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, int_scattering")) return 1;
    
    return 0;
}
int CalculateEnergy(const int b_size, const int N, const int M, std::vector<Type>& host_energy) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeEnergy << <blocks, threads >> > (N, M, dev_illum, dev_energy, dev_squares, square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeEnergy failed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_energy.data(), dev_energy, host_energy.size() * sizeof(host_energy[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, energy")) return 1;

    return 0;
}

int CalculateStream(const int b_size, const int N, const int M, std::vector<Vector3>& host_stream) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeStream << <blocks, threads >> > (N, M, dev_illum, dev_stream, dev_directions, dev_squares, square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeStreamfailed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_stream.data(), dev_stream, host_stream.size() * sizeof(host_stream[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, energy")) return 1;

    return 0;
}

int CalculateImpuls(const int b_size, const int N, const int M, std::vector<Matrix3>& host_impuls) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeImpuls << <blocks, threads >> > (N, M, dev_illum, dev_impuls, dev_directions, dev_squares, square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeImpuls failed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_impuls.data(), dev_impuls, host_impuls.size() * sizeof(host_impuls[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, impuls")) return 1;

    return 0;
}

int ClearDevice(const int mod) {

    cudaFree(dev_directions);
    cudaFree(dev_squares);
    cudaFree(dev_illum);
    cudaFree(dev_int_scattering);

    if (mod)
    {
        cudaFree(dev_energy);
        cudaFree(dev_stream);
        cudaFree(dev_impuls);
    }

    // cudaDeviceReset must be called before exiting in order for profiling and tracing tools such as Nsight and Visual Profiler to show complete traces.
    if (CheckError(cudaDeviceReset(), "cudaDeviceReset failed!")) return 1;

    return 0;
}
