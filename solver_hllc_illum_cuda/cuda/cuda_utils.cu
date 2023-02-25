//***********************************************************************//
//*********************Functions from device*****************************//
//***********************************************************************//

#include "../prj_config.h"
#ifdef USE_CUDA
#include "cuda_utils.cuh"

__device__ Type Gamma(const cuda_vector_t<Type, 3>& direction, const cuda_vector_t<Type, 3>& direction2)
{    
    Type sum = 0;
    for (int i = 0; i < 3; i++)
    {
        sum += direction[i] * direction2[i];
    }

    return (3. * (1 + sum * sum)) / 4;
}

__device__ Type IntegarteDirection(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid)    
{
    const int M = dir->size;
    const int N = grid->size;

    Type res = 0;
    for (int i = 0; i < M; i++)
    {
        Type I = 0;
        for (int k = 0; k < base; k++)
        {
            I += grid->illum[base*(N * i + num_cell) + k];
        }
        I /= base;        

        res += I * dir->directions[i].area;
    }

    return res / dir->full_area;
}

__device__ cuda_vector_t<Type,3> IntegarteDirection3(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid) 
{
    const int M = dir->size;
    const int N = grid->size;
    cuda_vector_t<Type, 3> res;

    for (int i = 0; i < M; i++) 
    {
        Type I = 0;
        for (int k = 0; k < base; k++)
        {
            I += grid->illum[base * (N * i + num_cell) + k];
        }
        I /= base;

        direction_device_t& buf = dir->directions[i];
        Type i_s = I * buf.area;

        res[0] += buf.dir[0] * i_s;
        res[1] += buf.dir[1] * i_s;
        res[2] += buf.dir[2] * i_s;
    }

    res[0] /= dir->full_area;
    res[1] /= dir->full_area;
    res[2] /= dir->full_area;

    return res;
}

__device__ cuda_vector_t<Type, 9>  IntegarteDirection9(const int num_cell, const grid_directions_device_t* dir, grid_device_t* grid) {

    const int M = dir->size;
    const int N = grid->size;
    cuda_vector_t<Type, 9> res;

    for (int i = 0; i < 3; i++)
        for (int k = 0; k < 3; k++)

            for (int j = 0; j < M; j++)
            {
                Type I = 0;
                for (int h = 0; h < base; h++)
                {
                    I += grid->illum[base * (N * j + num_cell) + h];
                }
                I /= base;

                direction_device_t& buf = dir->directions[j];

                res[i * 3 + k] += buf.dir[i] * buf.dir[k] * (I * buf.area);
            }


    for (int i = 0; i < 9; i++)
        res[i] /= dir->full_area;

    return res;
}

__device__ void IntegarteDirection9Faces(const int num_cell, const grid_directions_device_t* dir_grid, grid_device_t* grid, cuda_vector_t<Type, 9>* Impuls)
{
    const int M = dir_grid->size;
    const int N = grid->size;

    for (size_t h = 0; h < base; h++)
    {
        for (size_t i = 0; i < 9; i++)
        {
            Impuls[h][i] = 0;
        }
    }

    for (size_t dir = 0; dir < M; dir++)
    {
        for (size_t f = 0; f < base; f++)
        {
            for (size_t i = 0; i < 3; i++)
                for (size_t k = 0; k < 3; k++)
                {
                    Type I = grid->illum[base * (N * dir + num_cell) + f];
                    Impuls[f][i * 3 + k] += dir_grid->directions[dir].dir[i] * dir_grid->directions[dir].dir[k] * (I * dir_grid->directions[dir].area);
                }
        }
    }

    for (size_t h = 0; h < base; h++)
        for (size_t i = 0; i < 9; i++)
            Impuls[h][i] /= dir_grid->full_area;

    return;
}

__device__ void IntegarteDirection3Faces(const int num_cell, const grid_directions_device_t* dir_grid, grid_device_t* grid, cuda_vector_t<Type, 3>* Stream) {

    const int M = dir_grid->size;
    const int N = grid->size;


    for (size_t h = 0; h < base; h++)
    {
        for (size_t i = 0; i < 3; i++)
        {
            Stream[h][i] = 0;
        }
    }

    for (size_t i = 0; i < M; i++)
    {
        direction_device_t& dir = dir_grid->directions[i];

        for (size_t f = 0; f < base; f++)
        {
            Type I = grid->illum[base * (N * i + num_cell) + f];                     
            Type i_a = I * dir.area;

            Stream[f][0] += dir.dir[0] * i_a;
            Stream[f][1] += dir.dir[1] * i_a;
            Stream[f][2] += dir.dir[2] * i_a;
        }
    }


    for (size_t h = 0; h < base; h++)
        for (size_t i = 0; i < 3; i++)
            Stream[h][i] /= dir_grid->full_area;

    return;
}

#endif //USE_CUDA