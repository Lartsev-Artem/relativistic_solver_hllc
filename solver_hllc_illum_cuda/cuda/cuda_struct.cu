#include "../prj_config.h"
#ifdef USE_CUDA

#include "cuda_struct.cuh"

device_host_ptr_t device_host_ptr;
/*static*/ grid_directions_device_t* grid_dir_device_ptr;
/*static*/ grid_device_t* grid_cell_device_ptr;

#include "cuda_memory.cuh"

#define CUDA_INIT_GRID_PTR_DEVICE(dist, src, N) \
    CUDA_MALLOC(&src, N)\
    CUDA_MEMCPY_TO_DEVICE(&dist, &src, sizeof(src))


//**********************************************
static int grid_directions_init_copy_host_to_device(grid_directions_device_t*& grid_device, const grid_directions_t& grid_host)
{
    const int n = grid_host.size;

    CUDA_MALLOC((void**)&grid_device, sizeof(grid_directions_device_t));
    
    CUDA_MEMCPY_TO_DEVICE(&grid_device->size, &n, sizeof(grid_host.size));
    CUDA_MEMCPY_TO_DEVICE(&grid_device->full_area, &grid_host.full_area, sizeof(grid_host.full_area));
    
    std::vector<direction_device_t> directions(n);
    for (size_t i = 0; i < n; i++)
    {
        directions[i].area = grid_host.directions[i].area;
        directions[i].dir[0] = grid_host.directions[i].dir[0];
        directions[i].dir[1] = grid_host.directions[i].dir[1];
        directions[i].dir[2] = grid_host.directions[i].dir[2];
    }
        
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->directions, device_host_ptr.directions, n * sizeof(direction_device_t));    
    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.directions, directions.data(), n * sizeof(directions[0]));    
    return 0;
}

static int grid_directions_clear_device(grid_directions_device_t*& grid_device)
{       
    CUDA_FREE_MEMORY(device_host_ptr.directions);
    CUDA_FREE_MEMORY(grid_device);
    return 0;
}

//**********************************************

static int grid_init_device(grid_device_t*& grid_device, const int num_cell, const int num_dir)
{    
    CUDA_MALLOC((void**)&grid_device, sizeof(grid_device_t));
    
    CUDA_MEMCPY_TO_DEVICE(&grid_device->size, &num_cell, sizeof(num_cell));
        
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->illum, device_host_ptr.illum, (base * num_dir * num_cell * sizeof(Type)) );

    CUDA_INIT_GRID_PTR_DEVICE(grid_device->int_scattering, device_host_ptr.int_scattering, num_dir * num_cell * sizeof(Type));

    CUDA_INIT_GRID_PTR_DEVICE(grid_device->energy, device_host_ptr.energy, num_cell * sizeof(Type));

    CUDA_INIT_GRID_PTR_DEVICE(grid_device->divstream, device_host_ptr.divstream, num_cell * sizeof(Type));
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->divimpuls, device_host_ptr.divimpuls, num_cell * sizeof(cuda_vector_t<Type, 3>));

    CUDA_INIT_GRID_PTR_DEVICE(grid_device->volume, device_host_ptr.volume, num_cell * sizeof(Type));
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->areas, device_host_ptr.areas, base * num_cell * sizeof(Type));
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->normals, device_host_ptr.normals, base * num_cell * sizeof(cuda_vector_t<Type, 3>));

#ifdef CUDA_FULL_ARRAYS
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->energy, device_host_ptr.energy, num_cell * sizeof(Type));
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->stream, device_host_ptr.stream, num_cell * sizeof(cuda_vector_t<Type, 3>));
    CUDA_INIT_GRID_PTR_DEVICE(grid_device->impuls, device_host_ptr.impuls, num_cell * sizeof(cuda_vector_t<Type, 9>));
#endif

    return 0;
}


#include "../file_module/reader_bin.h"

static int grid_init_copy_host_to_device(const int num_dir, const std::vector<Type>& host_volume, const std::vector<Type>& host_areas, const std::vector<Normals>& host_normals,
    const grid_t& grid)
{       
    const int N = grid.size;
    
   // CUDA_MEMCPY_TO_DEVICE(device_host_ptr.illum, grid.Illum, (base * num_dir * N * sizeof(Type)));
   
    std::vector<cuda_vector_t<Type,3>> dev_norm(host_normals.size() * base);
    for (size_t i = 0; i < host_normals.size(); i++)
    {
        for (size_t j = 0; j < base; j++)
        {
            dev_norm[i * base + j].data[0] = host_normals[i].n[j][0];
            dev_norm[i * base + j].data[1] = host_normals[i].n[j][1];
            dev_norm[i * base + j].data[2] = host_normals[i].n[j][2];
        }
    }    

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.normals, dev_norm.data(), dev_norm.size() * sizeof(dev_norm[0]));

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.volume, host_volume.data(), N * sizeof(host_volume[0]));

    CUDA_MEMCPY_TO_DEVICE(device_host_ptr.areas, host_areas.data(), base * N * sizeof(host_areas[0]));

    return 0;
}

static int grid_clear_device(grid_device_t*& grid_device)
{       
    CUDA_FREE_MEMORY(device_host_ptr.illum);
    CUDA_FREE_MEMORY(device_host_ptr.int_scattering);
    CUDA_FREE_MEMORY(device_host_ptr.divimpuls);
    CUDA_FREE_MEMORY(device_host_ptr.divstream);

    CUDA_FREE_MEMORY(device_host_ptr.normals);
    CUDA_FREE_MEMORY(device_host_ptr.volume);
    CUDA_FREE_MEMORY(device_host_ptr.areas);

#ifdef CUDA_FULL_ARRAYS
    CUDA_FREE_MEMORY(device_host_ptr.energy);
    CUDA_FREE_MEMORY(device_host_ptr.stream);
    CUDA_FREE_MEMORY(device_host_ptr.impuls);
#endif
    
    CUDA_FREE_MEMORY(grid_device);

    return 0;
}


static int InitDevice(grid_directions_device_t*& grid_dir_device, grid_device_t*& grid_device, 
    const grid_directions_t& grid_dir_host, grid_t& grid_host)
{
    const std::string name_file_normals = BASE_ADRESS + "normals.bin";
    const std::string name_file_squares = BASE_ADRESS + "squares.bin";
    const std::string name_file_volume = BASE_ADRESS + "volume.bin";

    std::vector<Normals> normals;
    std::vector<Type> squares_faces;
    std::vector<Type> volume;

    if (ReadNormalFile(name_file_normals, normals)) RETURN_ERR("Error reading file normals\n");
    if (ReadSimpleFileBin(name_file_squares, squares_faces)) RETURN_ERR("Error reading file squares_faces\n");
    if (ReadSimpleFileBin(name_file_volume, volume)) RETURN_ERR("Error reading file volume\n");


    grid_directions_init_copy_host_to_device(grid_dir_device, grid_dir_host);

    grid_init_device(grid_device, grid_host.size, grid_dir_host.size);

    CUDA_MALLOC_HOST(&grid_host.Illum, (base * grid_dir_host.size * grid_host.size * sizeof(Type)));
    CUDA_MALLOC_HOST(&grid_host.scattering, (grid_dir_host.size * grid_host.size * sizeof(Type)));

    CUDA_MALLOC_HOST(&grid_host.divstream, (grid_host.size * sizeof(Type)));
    CUDA_MALLOC_HOST(&grid_host.divimpuls, (grid_host.size * sizeof(Vector3)));  

#ifdef ON_FULL_ILLUM_ARRAYS
    grid_host.energy = new Type[grid_host.size];
    grid_host.stream = new Vector3[grid_host.size];
    grid_host.impuls = new Matrix3[grid_host.size];
#endif

    grid_init_copy_host_to_device(grid_dir_host.size, volume, squares_faces, normals, grid_host);

    return 0;
}

static int ClearDevice(grid_directions_device_t*& grid_dir_device, grid_device_t*& grid_device)
{
    grid_directions_clear_device(grid_dir_device);
    grid_clear_device(grid_device);
    return 0;
}

void InitDevice(const grid_directions_t& grid_dir_host, grid_t& grid_host)
{
    InitDevice(grid_dir_device_ptr, grid_cell_device_ptr, grid_dir_host, grid_host);
}
void ClearDevice()
{
    grid_directions_clear_device(grid_dir_device_ptr);
    grid_clear_device(grid_cell_device_ptr);
}

void ClearHost(grid_t& grid_host)
{    
    CUDA_FREE_HOST_MEMORY(grid_host.Illum);
    CUDA_FREE_HOST_MEMORY(grid_host.scattering);

    CUDA_FREE_HOST_MEMORY(grid_host.divstream);
    CUDA_FREE_HOST_MEMORY(grid_host.divimpuls);

    delete[] grid_host.energy;
    delete[] grid_host.stream;
    delete[] grid_host.impuls;

    WRITE_LOG("Free host arrays\n");
}
#endif //USE_CUDA