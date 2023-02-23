#include "solve_short_characteristic_cuda.cuh"
#ifdef USE_CUDA

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "global_value.h"

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
Type dev_square_surface;
Type* dev_directions;
Type* dev_squares;
Type* dev_illum;
Type* dev_int_scattering;

Type* dev_energy;
Type* dev_stream;
Type* dev_impuls;

Type* dev_divimpuls;
Type* dev_divstream;

Type* dev_normals;
Type* dev_areas;
Type* dev_volume;

//**********************************************


static inline int CheckError(cudaError_t cudaStatus, const char* my_text = "") {
    if (cudaStatus != cudaSuccess) 
    {        
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

    for (size_t i = 0; i < M; i++) 
    {
        Type I = 0;
        for (size_t k = 0; k < base; k++)
        {
            I += illum[base*(N * i + num_cell) + k];
        }
        I /= base;        

        res += I * squares[i];
    }

    return res / square_surface;
}
__device__ dev_Vector3 IntegarteDirection3(const int N, const int M, int num_cell,
    Type* illum, const dev_Vector3* directions, Type* squares, const Type square_surface) {

    dev_Vector3 res;

    for (size_t i = 0; i < M; i++) 
    {
        Type I = 0;
        for (size_t k = 0; k < base; k++)
        {
            I += illum[base * (N * i + num_cell) + k];
        }
        I /= base;

        res[0] += directions[i][0] * I * squares[i];
        res[1] += directions[i][1] * I * squares[i];
        res[2] += directions[i][2] * I * squares[i];
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

            for (size_t j = 0; j < M; j++) 
            {
                Type I = 0;
                for (size_t h = 0; h < base; h++)
                {
                    I += illum[base * (N * j + num_cell) + h];
                }
                I /= base;

                res[i*3 +k] += directions[j][i] * directions[j][k] * (I * squares[j]);
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


__device__ void IntegarteDirection9Faces(const int N, const int M, int num_cell,
    Type* illum, const dev_Vector3* directions, Type* squares, const Type square_surface, dev_Matrix3* Impuls)
{

    for (size_t h = 0; h < base; h++)
    {
        for (size_t i = 0; i < 9; i++)
        {
            Impuls[h].data[i] = 0;
        }
    }

    for (size_t dir = 0; dir < M; dir++)
    {
        for (size_t f = 0; f < base; f++)
        {
            for (size_t i = 0; i < 3; i++)
                for (size_t k = 0; k < 3; k++)
                {
                    Type I = illum[base * (N * dir + num_cell) + f];
                    Impuls[f][i * 3 + k] += directions[dir][i] * directions[dir][k] * (I * squares[dir]);
                }
        }
    }

    for (size_t h = 0; h < base; h++)
        for (size_t i = 0; i < 9; i++)
            Impuls[h][i] /= square_surface;

    return;
}

__global__ void d_MakeDivImpuls(const int N, const int M, Type* illum, Type* d_divimpuls, Type* d_directions, Type* squares, const Type square_surface,
    const Type* d_normals, const Type* areas, const Type* volume) {


    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    const dev_Vector3* directions = (dev_Vector3*)d_directions;
    dev_Matrix3 impuls[base];

    //for (size_t i = 0; i < n; ++i)
    
    IntegarteDirection9Faces(N, M, i, illum, directions, squares, square_surface, impuls);

    dev_Vector3& divimpuls = (((dev_Vector3*)d_divimpuls)[i]);

    const dev_Vector3* normals = (dev_Vector3*)d_normals;

    for (size_t j = 0; j < 3; j++)
        divimpuls[j] = 0;

    for (size_t j = 0; j < base; j++)
    {
        for (size_t h = 0; h < 3; h++)
        {
            Type sum = 0;
            for (size_t k = 0; k < 3; k++)
            {
                sum += impuls[j][h * 3 + k] * normals[i * base + j][k];
            }

            divimpuls[h] += sum * areas[i * base + j];
        }        
    }

    for (size_t j = 0; j < 3; j++)
        divimpuls[j] /= volume[i];

    return;
}


__device__ void IntegarteDirection3Faces(const int N, const int M, int num_cell,
    Type* illum, const dev_Vector3* directions, Type* squares, const Type square_surface, dev_Vector3* Stream) {

    for (size_t h = 0; h < base; h++)
    {
        for (size_t i = 0; i < 3; i++)
        {
            Stream[h].data[i] = 0;
        }
    }

    for (size_t i = 0; i < M; i++)
    {
        for (size_t f = 0; f < base; f++)
        {
            Type I = illum[base * (N * i + num_cell) + f];         

            Stream[f][0] += directions[i][0] * I * squares[i];
            Stream[f][1] += directions[i][1] * I * squares[i];
            Stream[f][2] += directions[i][2] * I * squares[i];
        }
    }


    for (size_t h = 0; h < base; h++)
        for (size_t i = 0; i < 3; i++)
            Stream[h][i] /= square_surface;

    return;
}

__global__ void d_MakeDivStream(const int N, const int M, Type* illum, Type* d_divstream, Type* d_directions, Type* squares, const Type square_surface,
    const Type* d_normals, const Type* areas, const Type* volume) 
{


    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= N) return;

    const dev_Vector3* directions = (dev_Vector3*)d_directions;
    dev_Vector3 Stream[base];

    //for (size_t i = 0; i < n; ++i)

    IntegarteDirection3Faces(N, M, i, illum, directions, squares, square_surface, Stream);
    
    const dev_Vector3* normals = (dev_Vector3*)d_normals;

    d_divstream[i] = 0;

    for (size_t f = 0; f < base; f++)
    {
        Type sum = 0;
        for (size_t k = 0; k < 3; k++)
        {
            sum += Stream[f][k] * normals[i * base + f][k];
        } 
        d_divstream[i] += sum * areas[i * base + f];
    }

    d_divstream[i] /= volume[i];
    return;
}

//***********************************************************************//
//*********************Functions from host*******************************//
//***********************************************************************//

int CheckDevice()     //проверка наличия карты 
{
    return CheckError(cudaSetDevice(0), "cudaSetDevice failed!Do you have a CUDA - capable GPU installed ?");
}

int InitDevice(const int num_dir, const int num_cell, const int mod) {

    if (CheckError(cudaMalloc(&dev_directions, num_dir * NUMBER_OF_MEASUREMENTS * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_squares, num_dir * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_illum, 4 * num_dir * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
    if (CheckError(cudaMalloc(&dev_int_scattering, num_dir * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;

    if (mod) //extended
    {
        if (CheckError(cudaMalloc(&dev_energy, num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_stream, 3 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_impuls, 9 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;

        if (CheckError(cudaMalloc(&dev_divstream, num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_divimpuls, 3 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;


        if (CheckError(cudaMalloc(&dev_volume, num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_areas, base*num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
        if (CheckError(cudaMalloc(&dev_normals, base*3 * num_cell * sizeof(double)), "cudaMalloc failed!")) return 1;
    }

    return 0;
}
int HostToDeviceInit(const grid_directions_t& host_directions, std::vector<Type>& host_illum,
    std::vector<Type>& host_volume, std::vector<Type>& host_areas, std::vector<Normals>& host_normals)
{

    if (CheckError(cudaMemcpy(dev_illum, host_illum.data(), host_illum.size() * sizeof(host_illum[0]),
        cudaMemcpyHostToDevice), "cudaMemcpy failed illum!")) return 1;


    std::vector<dev_Vector3> dev_dir(host_directions.size);  // может можно сразу?
    std::vector<Type>host_squares(host_directions.size);

    int i = 0;
    for (auto& el : host_directions.directions)
    {
        dev_dir[i].data[0] = el.dir[0];
        dev_dir[i].data[1] = el.dir[1];
        dev_dir[i].data[2] = el.dir[2];

        host_squares[i] = el.area;
        i++;
    }

    dev_square_surface = host_directions.full_area;
    /*if (CheckError(cudaMemcpy(&dev_square_surface, &host_directions.full_area, 1 * sizeof(double), cudaMemcpyHostToDevice),
        "cudaMemcpy failed square_surface!")) return 1;*/

    if (CheckError(cudaMemcpy(dev_directions, dev_dir.data(), dev_dir.size() * sizeof(dev_Vector3), cudaMemcpyHostToDevice),
        "cudaMemcpy failed dir!")) return 1;

    if (CheckError(cudaMemcpy(dev_squares, host_squares.data(), host_squares.size() * sizeof(host_squares[0]), cudaMemcpyHostToDevice),
        "cudaMemcpy failed! squares")) return 1;
    

    WRITE_LOG("Start Init dev_norm\n");

    std::vector<dev_Vector3> dev_norm(host_normals.size() * base);
    for (size_t i = 0; i < host_normals.size(); i++)
    {
        for (size_t j = 0; j < base; j++)
        {
            dev_norm[i * base + j].data[0] = host_normals[i].n[j][0];
            dev_norm[i * base + j].data[1] = host_normals[i].n[j][1];
            dev_norm[i * base + j].data[2] = host_normals[i].n[j][2];
        }
    }
    WRITE_LOG("Start Init normals\n");

    if (CheckError(cudaMemcpy(dev_normals, dev_norm.data(), dev_norm.size() * sizeof(dev_norm[0]), cudaMemcpyHostToDevice),
        "cudaMemcpy failed! normals")) return 1;


    WRITE_LOG("Start Init volume\n");

    if (CheckError(cudaMemcpy(dev_volume, host_volume.data(), host_volume.size() * sizeof(host_volume[0]), cudaMemcpyHostToDevice),
        "cudaMemcpy failed! volume")) return 1;

    WRITE_LOG("Start Init areas\n");

    if (CheckError(cudaMemcpy(dev_areas, host_areas.data(), host_areas.size() * sizeof(host_areas[0]), cudaMemcpyHostToDevice),
        "cudaMemcpy failed! areas")) return 1;



    WRITE_LOG("final Init\n");

    return 0;
}

int HostToDevice(const grid_directions_t& host_directions, std::vector<Type>& host_illum, const int mod) {
    
    if (CheckError(cudaMemcpy(dev_illum, host_illum.data(), host_illum.size() * sizeof(host_illum[0]), 
        cudaMemcpyHostToDevice), "cudaMemcpy failed illum!")) return 1;

    if (mod)  // first init
    {
        std::vector<dev_Vector3> dev_dir(host_directions.size);  // может можно сразу?
        std::vector<Type>host_squares(host_directions.size);

        int i = 0;
        for (auto &el:host_directions.directions)
        {
            dev_dir[i].data[0] = el.dir[0];
            dev_dir[i].data[1] = el.dir[1];
            dev_dir[i].data[2] = el.dir[2];

            host_squares[i] = el.area;
            i++;
        }
                
        dev_square_surface = host_directions.full_area;
        /*if (CheckError(cudaMemcpy(&dev_square_surface, &host_directions.full_area, 1 * sizeof(double), cudaMemcpyHostToDevice),
            "cudaMemcpy failed square_surface!")) return 1;*/

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

    d_GetS << <blocks, threads >> > (N, M, dev_illum, dev_int_scattering, dev_directions, dev_squares, dev_square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_GetS failed!")) return 1;

    if (CheckError(cudaMemcpy(host_int_scattering.data(), dev_int_scattering, host_int_scattering.size() * sizeof(host_int_scattering[0]), 
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, int_scattering")) return 1;
    
    return 0;
}

int CalculateEnergy(const int b_size, const int N, const int M, std::vector<Type>& host_energy) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeEnergy << <blocks, threads >> > (N, M, dev_illum, dev_energy, dev_squares, dev_square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeEnergy failed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_energy.data(), dev_energy, host_energy.size() * sizeof(host_energy[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, energy")) return 1;

    return 0;
}

int CalculateStream(const int b_size, const int N, const int M, std::vector<Vector3>& host_stream) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeStream << <blocks, threads >> > (N, M, dev_illum, dev_stream, dev_directions, dev_squares, dev_square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeStreamfailed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_stream.data(), dev_stream, host_stream.size() * sizeof(host_stream[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, energy")) return 1;

    return 0;
}

int CalculateImpuls(const int b_size, const int N, const int M, std::vector<Matrix3>& host_impuls) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeImpuls << <blocks, threads >> > (N, M, dev_illum, dev_impuls, dev_directions, dev_squares, dev_square_surface);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeImpuls failed!")) return 1;
    
    if (CheckError(cudaMemcpy(host_impuls.data(), dev_impuls, host_impuls.size() * sizeof(host_impuls[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, impuls")) return 1;

    return 0;
}

int CalculateDivImpuls(const int b_size, const int N, const int M, std::vector<Vector3>& host_div_impuls) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeDivImpuls << <blocks, threads >> > (N, M, dev_illum, dev_divimpuls, dev_directions, dev_squares, dev_square_surface,
        dev_normals, dev_areas,dev_volume);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeDivImpuls failed!")) return 1;

    if (CheckError(cudaMemcpy(host_div_impuls.data(), dev_divimpuls, host_div_impuls.size() * sizeof(host_div_impuls[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, div_impuls")) return 1;

    return 0;
}
int CalculateDivStream(const int b_size, const int N, const int M, std::vector<Type>& host_divstream) {

    dim3 threads(b_size);
    dim3 blocks((N + b_size - 1) / b_size);

    d_MakeDivStream << <blocks, threads >> > (N, M, dev_illum, dev_divstream, dev_directions, dev_squares, dev_square_surface,
        dev_normals, dev_areas, dev_volume);

    // Check for any errors launching the kernel
    if (CheckError(cudaGetLastError(), "d_MakeDivStream failed!")) return 1;

    if (CheckError(cudaMemcpy(host_divstream.data(), dev_divstream, host_divstream.size() * sizeof(host_divstream[0]),
        cudaMemcpyDeviceToHost), "cudaMemcpy failed!, d_MakeDivStream")) return 1;

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
        cudaFree(dev_divstream);
        cudaFree(dev_divimpuls);

        cudaFree(dev_areas);
        cudaFree(dev_volume);
        cudaFree(dev_normals);
    }

    // cudaDeviceReset must be called before exiting in order for profiling and tracing tools such as Nsight and Visual Profiler to show complete traces.
    if (CheckError(cudaDeviceReset(), "cudaDeviceReset failed!")) return 1;

    return 0;
}
#endif //USE_CUDA