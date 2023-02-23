CC = mpiicc
CPC = mpiicpc
NVC = nvcc 
CCFLAGS = -fPIE  -O3 -openmp -c
NVCFLAGS = -O2 -gencode arch=compute_70,code=sm_70 -c  -Xcompiler " -openmp"
NVCFLAGSRUN = -ccbin= $(CPC) -Xcompiler " -fopenmp" -o

#.DEFAULT:
#-touch $@

all: solution

FILE_OBJ = writer_bin.o reader_bin.o reader_txt.o

UTILS_OBJ = geometry_data.o geometry_solve.o get_grid_data.o rebuild_grid.o rebuild_solve.o utils_main.o

GRAPH_OBJ = build_graph_main.o build_graph_calculation.o build_graph_read_write.o

MAKE_OBJ = struct_short_characteristics_main.o struct_short_characteristics_calculations.o struct_short_characteristics_logic_function.o

HLLC_OBJ = hllc_1d.o hllc_2d.o hllc_3d.o hllc_utils.o

RHLCC_OBJ = rhllc_1d.o rhllc_2d.o rhllc_3d.o rhllc_utils.o rhllc_2d_mpi.o rhllc_3d_mpi.o

ILLUM_OBJ = illum_utils.o Illum_part.o Illum_part_mpi.o Illum_part_mpi_omp.o

SOLVE_OBJ = solve_utils.o solve_global_struct.o solve_main.o $(HLLC_OBJ) $(RHLLC_OBJ) $(ILLUM_OBJ)

writer_bin.o: $(SRC_FILE)/writer_bin.cpp illum_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_FILE)/writer_bin.cpp

reader_bin.o: $(SRC_FILE)/reader_bin.cpp \
    $(CPC) $(CCFLAGS) $(SRC_FILE)/reader_bin.cpp

reader_txt.o: $(SRC_FILE)/reader_txt.cpp \
    $(CPC) $(CCFLAGS) $(SRC_FILE)/reader_txt.cpp

geometry_data.o: $(SRC_UTILS)/grid_geometry/geometry_data.cpp \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/grid_geometry/geometry_data.cpp

geometry_solve.o: $(SRC_UTILS)/grid_geometry/geometry_solve.cpp \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/grid_geometry/geometry_solve.cpp

get_grid_data.o: $(SRC_UTILS)/grid_geometry/get_grid_data.cpp $(SRC_FILE)/reader_vtk.o \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/grid_geometry/get_grid_data.cpp

rebuild_grid.o: $(SRC_UTILS)/rebuild_grid.cpp writer_bin.o reader_txt.o reader_bin.o \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/rebuild_grid.cpp

rebuild_solve.o: $(SRC_UTILS)/rebuild_solve.cpp writer_bin.o reader_bin.o \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/rebuild_solve.cpp

utils_main.o: $(SRC_UTILS)/utils_main.cpp rebuild_solve.o rebuild_grid.o get_grid_data.o \
$(CPC) $(CCFLAGS) $(SRC_UTILS)/utils_main.cpp


build_graph_read_write.o: $(SRC_GRAPH)/build_graph_read_write.cpp geometry_data.o writer_bin.o reader_bin.o \
$(CPC) $(CCFLAGS) $(SRC_GRAPH)/build_graph_read_write.cpp

build_graph_calculation.o: $(SRC_GRAPH)/build_graph_calculation.cpp geometry_data.o geometry_solve.o \
$(CPC) $(CCFLAGS) $(SRC_GRAPH)/build_graph_calculation.cpp

build_graph_main.o: $(SRC_GRAPH)/build_graph_main.cpp build_graph_calculation.o build_graph_read_write.o writer_bin.o reader_bin.o reader_txt.o \
$(CPC) $(CCFLAGS) $(SRC_GRAPH)/build_graph_main.cpp

struct_short_characteristics_calculations.o: $(SRC_MAKE)/struct_short_characteristics_calculations.cpp writer_bin.o geometry_data.o \
$(CPC) $(CCFLAGS) $(SRC_MAKE)/struct_short_characteristics_calculations.cpp

struct_short_characteristics_logic_function.o: $(SRC_MAKE)/struct_short_characteristics_logic_function.cpp struct_short_characteristics_calculations.o geometry_solve.o \
$(CPC) $(CCFLAGS) $(SRC_MAKE)/struct_short_characteristics_logic_function.cpp

struct_short_characteristics_main.o: $(SRC_MAKE)/struct_short_characteristics_main.cpp \
struct_short_characteristics_calculations.o  struct_short_characteristics_logic_function.o \
writer_bin.o reader_bin.o reader_txt.o geometry_solve.o \
$(CPC) $(CCFLAGS) $(SRC_MAKE)/struct_short_characteristics_main.cpp


hllc_utils.o: $(SRC_HLCC)/hllc_utils.cpp reader_txt.o reader_bin.o \
    $(CPC) $(CCFLAGS) $(SRC_HLCC)/hllc_utils.cpp

hllc_1d.o: $(SRC_HLCC)/hllc_1d.cpp \
    $(CPC) $(CCFLAGS) $(SRC_HLCC)/hllc_1d.cpp

hllc_2d.o: $(SRC_HLCC)/hllc_2d.cpp writer_bin.o geometry_solve.o hllc_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_HLCC)/hllc_2d.cpp

hllc_3d.o: $(SRC_HLCC)/hllc_3d.cpp reader_bin.o geometry_solve.o \
    $(CPC) $(CCFLAGS) $(SRC_HLCC)/hllc_3d.cpp


rhllc_utils.o: $(SRC_RHLCC)/rhllc_utils.cpp reader_txt.o reader_bin.o \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_utils.cpp

rhllc_1d.o: $(SRC_RHLCC)/rhllc_1d.cpp \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_1d.cpp

rhllc_2d.o: $(SRC_RHLCC)/rhllc_2d.cpp writer_bin.o geometry_solve.o rhllc_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_2d.cpp

rhllc_3d.o: $(SRC_RHLCC)/rhllc_3d.cpp reader_bin.o geometry_solve.o writer_bin.o \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_3d.cpp

rhllc_2d_mpi.o: $(SRC_RHLCC)/rhllc_2d_mpi.cpp writer_bin.o geometry_solve.o rhllc_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_2d_mpi.cpp

rhllc_3d_mpi.o: $(SRC_RHLCC)/rhllc_3d_mpi.cpp reader_bin.o geometry_solve.o \
    $(CPC) $(CCFLAGS) $(SRC_RHLCC)/rhllc_3d_mpi.cpp


illum_utils.o: $(SRC_ILLUM)/illum_utils.cpp reader_bin.o \
    $(CPC) $(CCFLAGS) $(SRC_ILLUM)/illum_utils.cpp

Illum_part.o: $(SRC_ILLUM)/Illum_part.cpp reader_txt.o illum_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_ILLUM)/Illum_part.cpp


Illum_part_mpi.o: $(SRC_ILLUM)/Illum_part_mpi.cpp reader_bin.o reader_txt.o illum_utils.o writer_bin.o solve_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_ILLUM)/Illum_part_mpi.cpp


Illum_part_mpi_omp.o: $(SRC_ILLUM)/Illum_part_mpi_omp.cpp reader_bin.o reader_txt.o illum_utils.o writer_bin.o solve_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_ILLUM)/Illum_part_mpi_omp.cpp


solve_utils.o: $(SRC_SOLVE)/solve_utils.cpp reader_bin.o writer_bin.o hllc_utils.o rhllc_utils.o \
    $(CPC) $(CCFLAGS) $(SRC_SOLVE)/solve_utils.cpp

solve_global_struct.o: $(SRC_SOLVE)/solve_global_struct.cpp \
$(CPC) $(CCFLAGS) $(SRC_SOLVE)/solve_global_struct.cpp

solve_main.o: $(SRC_SOLVE)/solve_main.cpp reader_bin.o reader_txt.o writer_bin.o geometry_solve.o solve_utils.o hllc_utils.o illum_utils.o solve_global_struct.o \
$(CPC) $(CCFLAGS) $(SRC_SOLVE)/solve_main.cpp


main.o: ./../main.cpp \
    $(CPC) $(CCFLAGS) ./../main.cpp

kernel.o: ./../kernel.cu \
    $(NVC) $(NVCFLAGS) ./../kernel.cu


SRC_MAIN = ./..

SRC_FILE = $(SRC_MAIN)/file_module

SRC_UTILS = $(SRC_MAIN)/utils

SRC_GRAPH = $(SRC_MAIN)/build_graph

SRC_MAKE = $(SRC_MAIN)/make_illum_struct

SRC_SOLVE = $(SRC_MAIN)/solve_module

SRC_HLCC = ($SRC_SOLVE)/hllc

SRC_RHLCC = ($SRC_SOLVE)/rhllc

SRC_ILLUM = ($SRC_SOLVE)/illum

OBJ = main.o kernel.o $(FILE_OBJ) $(UTILS_OBJ) $(GRAPH_OBJ) $(MAKE_OBJ) $(SOLVE_OBJ)

solution: $(OBJ) \
    $(NVC) $(OBJ) $(NVCFLAGSRUN) run
