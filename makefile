CPC = mpiicpc
NVCC = nvcc
OUT = run

DEFINES = CLASTER SOLVE USE_CUDA #BUILD MAKE

ONDEF =  $(addprefix -D , $(DEFINES))

CPP_FLAGS  = $(ONDEF) -fPIE -Ofast -fopenmp

NVCC_OTHER_FLAGS = -Xcompiler "-fopenmp"
NVCC_FLAGS = $(ONDEF) -O2 -gencode arch=compute_70,code=sm_70 -dc $(NVCC_OTHER_FLAGS)

SRC_BASE = ./..

SRC_FILE = $(SRC_BASE)/file_module

SRC_UTILS = $(SRC_BASE)/utils

SRC_GRAPH = $(SRC_BASE)/build_graph

SRC_MAKE = $(SRC_BASE)/make_illum_struct

SRC_SOLVE = $(SRC_BASE)/solve_module

SRC_HLCC = $(SRC_SOLVE)/hllc

SRC_RHLCC = $(SRC_SOLVE)/rhllc

SRC_ILLUM = $(SRC_SOLVE)/illum

SRC_CUDA = $(SRC_BASE)/cuda


CUDA_OBJ = cuda_solve.o cuda_kernel.o cuda_struct.o cuda_utils.o

FILE_OBJ = writer_bin.o reader_bin.o reader_txt.o

UTILS_OBJ = geometry_data.o geometry_solve.o get_grid_data.o rebuild_grid.o rebuild_solve.o utils_main.o

STRUCT_OBJ = struct_short_characteristics_calculations.o struct_short_characteristics_logic_function.o struct_short_characteristics_main.o

BUILD_OBJ = build_graph_read_write.o build_graph_main.o build_graph_calculation.o

HLLC_OBJ = hllc_1d.o hllc_2d.o hllc_3d.o hllc_utils.o

RHLCC_OBJ = rhllc_1d.o rhllc_2d.o rhllc_3d.o rhllc_utils.o rhllc_2d_mpi.o rhllc_3d_mpi.o

ILLUM_OBJ = illum_utils.o Illum_part.o Illum_part_mpi.o Illum_part_mpi_omp.o Illum_async.o

SOLVE_OBJ = solve_utils.o solve_global_struct.o solve_main.o $(HLLC_OBJ) $(ILLUM_OBJ) $(RHLCC_OBJ) $(CUDA_OBJ)

OBJ = main.o  $(FILE_OBJ) $(UTILS_OBJ) $(STRUCT_OBJ) $(BUILD_OBJ) $(SOLVE_OBJ)

all: solution

writer_bin.o: $(SRC_FILE)/writer_bin.cpp illum_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_FILE)/writer_bin.cpp

reader_bin.o: $(SRC_FILE)/reader_bin.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_FILE)/reader_bin.cpp

reader_txt.o: $(SRC_FILE)/reader_txt.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_FILE)/reader_txt.cpp
 
 
geometry_data.o: $(SRC_UTILS)/grid_geometry/geometry_data.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/grid_geometry/geometry_data.cpp

geometry_solve.o: $(SRC_UTILS)/grid_geometry/geometry_solve.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/grid_geometry/geometry_solve.cpp

get_grid_data.o: $(SRC_UTILS)/grid_geometry/get_grid_data.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/grid_geometry/get_grid_data.cpp

rebuild_grid.o: $(SRC_UTILS)/rebuild_grid.cpp writer_bin.o reader_txt.o reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/rebuild_grid.cpp

rebuild_solve.o: $(SRC_UTILS)/rebuild_solve.cpp writer_bin.o reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/rebuild_solve.cpp

utils_main.o: $(SRC_UTILS)/utils_main.cpp rebuild_solve.o rebuild_grid.o get_grid_data.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_UTILS)/utils_main.cpp 


build_graph_read_write.o: $(SRC_GRAPH)/build_graph_read_write.cpp geometry_data.o writer_bin.o reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_GRAPH)/build_graph_read_write.cpp

build_graph_calculation.o: $(SRC_GRAPH)/build_graph_calculation.cpp geometry_data.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_GRAPH)/build_graph_calculation.cpp

build_graph_main.o: $(SRC_GRAPH)/build_graph_main.cpp build_graph_calculation.o build_graph_read_write.o writer_bin.o reader_bin.o reader_txt.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_GRAPH)/build_graph_main.cpp


struct_short_characteristics_calculations.o: $(SRC_MAKE)/struct_short_characteristics_calculations.cpp writer_bin.o geometry_data.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_MAKE)/struct_short_characteristics_calculations.cpp

struct_short_characteristics_logic_function.o: $(SRC_MAKE)/struct_short_characteristics_logic_function.cpp struct_short_characteristics_calculations.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_MAKE)/struct_short_characteristics_logic_function.cpp

struct_short_characteristics_main.o: $(SRC_MAKE)/struct_short_characteristics_main.cpp struct_short_characteristics_calculations.o struct_short_characteristics_logic_function.o writer_bin.o reader_bin.o reader_txt.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_MAKE)/struct_short_characteristics_main.cpp

hllc_utils.o: $(SRC_HLCC)/hllc_utils.cpp reader_txt.o reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_HLCC)/hllc_utils.cpp

hllc_1d.o: $(SRC_HLCC)/hllc_1d.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_HLCC)/hllc_1d.cpp

hllc_2d.o: $(SRC_HLCC)/hllc_2d.cpp writer_bin.o geometry_solve.o hllc_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_HLCC)/hllc_2d.cpp

hllc_3d.o: $(SRC_HLCC)/hllc_3d.cpp reader_bin.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_HLCC)/hllc_3d.cpp


rhllc_utils.o: $(SRC_RHLCC)/rhllc_utils.cpp reader_txt.o reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_utils.cpp

rhllc_1d.o: $(SRC_RHLCC)/rhllc_1d.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_1d.cpp

rhllc_2d.o: $(SRC_RHLCC)/rhllc_2d.cpp writer_bin.o geometry_solve.o rhllc_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_2d.cpp

rhllc_3d.o: $(SRC_RHLCC)/rhllc_3d.cpp reader_bin.o geometry_solve.o writer_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_3d.cpp

rhllc_2d_mpi.o: $(SRC_RHLCC)/rhllc_2d_mpi.cpp writer_bin.o geometry_solve.o rhllc_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_2d_mpi.cpp

rhllc_3d_mpi.o: $(SRC_RHLCC)/rhllc_3d_mpi.cpp reader_bin.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_RHLCC)/rhllc_3d_mpi.cpp

illum_utils.o: $(SRC_ILLUM)/illum_utils.cpp reader_bin.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_ILLUM)/illum_utils.cpp

Illum_part.o: $(SRC_ILLUM)/Illum_part.cpp reader_txt.o illum_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_ILLUM)/Illum_part.cpp


Illum_part_mpi.o: $(SRC_ILLUM)/Illum_part_mpi.cpp reader_bin.o reader_txt.o illum_utils.o writer_bin.o solve_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_ILLUM)/Illum_part_mpi.cpp


Illum_part_mpi_omp.o: $(SRC_ILLUM)/Illum_part_mpi_omp.cpp reader_bin.o reader_txt.o illum_utils.o writer_bin.o solve_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_ILLUM)/Illum_part_mpi_omp.cpp
 
Illum_async.o: $(SRC_ILLUM)/Illum_async.cpp reader_bin.o reader_txt.o illum_utils.o writer_bin.o solve_utils.o cuda_solve.o geometry_solve.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_ILLUM)/Illum_async.cpp


solve_utils.o: $(SRC_SOLVE)/solve_utils.cpp reader_bin.o writer_bin.o hllc_utils.o rhllc_utils.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_SOLVE)/solve_utils.cpp

solve_global_struct.o: $(SRC_SOLVE)/solve_global_struct.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_SOLVE)/solve_global_struct.cpp
 
solve_main.o: $(SRC_SOLVE)/solve_main.cpp reader_bin.o reader_txt.o writer_bin.o geometry_solve.o solve_utils.o hllc_utils.o illum_utils.o solve_global_struct.o
	$(CPC) $(CPP_FLAGS) -c $(SRC_SOLVE)/solve_main.cpp 

cuda_solve.o: $(SRC_CUDA)/cuda_solve.cu cuda_kernel.o
	$(NVCC) $(NVCC_FLAGS) -c $(SRC_CUDA)/cuda_solve.cu
 
cuda_kernel.o: $(SRC_CUDA)/cuda_kernel.cu cuda_struct.o cuda_utils.o
	$(NVCC) $(NVCC_FLAGS) -c $(SRC_CUDA)/cuda_kernel.cu

cuda_struct.o: $(SRC_CUDA)/cuda_struct.cu
	$(NVCC) $(NVCC_FLAGS) -c $(SRC_CUDA)/cuda_struct.cu
  
cuda_utils.o: $(SRC_CUDA)/cuda_utils.cu
	$(NVCC) $(NVCC_FLAGS) -c $(SRC_CUDA)/cuda_utils.cu   
 

main.o: $(SRC_BASE)/main.cpp
	$(CPC) $(CPP_FLAGS) -c $(SRC_BASE)/main.cpp 
 


solution: $(OBJ)
	$(NVCC) -ccbin=$(CPC) $(NVCC_OTHER_FLAGS) -o $(OUT) $(OBJ)
  
graph: main.o $(BUILD_OBJ)
	$(CPC) $(CPP_FLAGS) -o $(OUT) main.o $(BUILD_OBJ) $(FILE_OBJ) geometry_solve.o geometry_data.o
 
illum_struct: main.o $(STRUCT_OBJ)
	$(CPC) $(CPP_FLAGS) -o $(OUT) main.o $(STRUCT_OBJ) $(FILE_OBJ) geometry_solve.o
 
solve: main.o $(SOLVE_OBJ) $(FILE_OBJ) geometry_solve.o
	$(NVCC) -ccbin=$(CPC) $(NVCC_OTHER_FLAGS) -o $(OUT) main.o $(SOLVE_OBJ) $(FILE_OBJ) geometry_solve.o
 
utils: main.o $(UTILS_OBJ)
	$(CPC) $(CPP_FLAGS) -o $(OUT) main.o $(UTILS_OBJ) 
 


clean:
	-rm -f .cppdefs $(OBJ) $(OUT)