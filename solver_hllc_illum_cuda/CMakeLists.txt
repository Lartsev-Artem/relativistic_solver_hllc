cmake_minimum_required(VERSION 3.15)

project(solution LANGUAGES CXX)


#set(CMAKE_CXX_COMPILER icpc)
#project(test LANGUAGES CXX CUDA)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp")

find_package(OpenMP REQUIRED)


#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -O0 -Wformat -Wstack-protector -Wextra -fstack-protector-all -fno-omit-frame-pointer")


 #Общие настройки сборки
#set(CMAKE_CXX_FLAGS "-O3" )
#set(CMAKE_CXX_FLAGS "-std=c++11" )
#set (CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-fopenmp")

set (MAIN_DIR ${CMAKE_CURRENT_BINARY_DIR})
set (ROOT_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)


set (UTILS_DIR ${MAIN_DIR}/utils)
set (SOLVE_DIR ${MAIN_DIR}/solve_module)
set (GRAPH_DIR ${MAIN_DIR}/build_graph)
set (MAKE_DIR ${MAIN_DIR}/make_illum_struct)
set (FILE_DIR ${MAIN_DIR}/file_module)

message("UTILS_DIR=" ${UTILS_DIR})

include_directories(
        ${MAIN_DIR}
        ${MAIN_DIR}/build_graph/
        ${MAIN_DIR}/file_module/
        ${MAIN_DIR}/make_illum_struct/
        ${MAIN_DIR}/utils/
        ${MAIN_DIR}/utils/grid_geometry/
        ${MAIN_DIR}/solve_module/
        ${MAIN_DIR}/solve_module/hllc/
        ${MAIN_DIR}/solve_module/rhllc/
        ${MAIN_DIR}/solve_module/illum/
)

set (UTILS_SRC
        ${UTILS_DIR}/rebuild_grid.cpp
        ${UTILS_DIR}/rebuild_solve.cpp
        ${UTILS_DIR}/utils_main.cpp
        ${UTILS_DIR}/grid_geometry/geometry_data.cpp
        ${UTILS_DIR}/grid_geometry/geometry_solve.cpp
        ${UTILS_DIR}/grid_geometry/get_grid_data.cpp
)
#add_library(UTILS STATIC ${UTILS_SRC})

set (GRAPH_SRC
        ${GRAPH_DIR}/build_graph_calculation.cpp
        ${GRAPH_DIR}/build_graph_main.cpp
        ${GRAPH_DIR}/build_graph_read_write.cpp
)
#add_library(GRAPH STATIC ${GRAPH_SRC})

set (MAKE_SRC
        ${MAKE_DIR}/struct_short_characteristics_calculations.cpp
        ${MAKE_DIR}/struct_short_characteristics_logic_function.cpp
        ${MAKE_DIR}/struct_short_characteristics_main.cpp
)
#add_library(MAKE STATIC ${MAKE_SRC})

set (FILE_SRC
        ${FILE_DIR}/reader_bin.cpp
        ${FILE_DIR}/reader_txt.cpp
        ${FILE_DIR}/writer_bin.cpp
)
#add_library(FILE STATIC ${FILE_SRC})

set (SOLVE_SRC
        ${SOLVE_DIR}/solve_global_struct.cpp
        ${SOLVE_DIR}/solve_main.cpp
        ${SOLVE_DIR}/solve_utils.cpp

        ${SOLVE_DIR}/hllc/hllc_1d.cpp
        ${SOLVE_DIR}/hllc/hllc_2d.cpp
        ${SOLVE_DIR}/hllc/hllc_3d.cpp
        ${SOLVE_DIR}/hllc/hllc_utils.cpp       

        ${SOLVE_DIR}/rhllc/rhllc_1d.cpp
        ${SOLVE_DIR}/rhllc/rhllc_2d.cpp
        ${SOLVE_DIR}/rhllc/rhllc_3d.cpp
        ${SOLVE_DIR}/rhllc/rhllc_2d_mpi.cpp
        ${SOLVE_DIR}/rhllc/rhllc_3d_mpi.cpp
        ${SOLVE_DIR}/rhllc/rhllc_utils.cpp

        ${SOLVE_DIR}/illum/Illum_part.cpp
        ${SOLVE_DIR}/illum/illum_utils.cpp        
)
#add_library(SOLVE STATIC ${SOLVE_SRC})

set (MAIN_SCR
        ${MAIN_DIR}/main.cpp
)

add_executable(solution ${MAIN_SCR} ${SOLVE_SRC} ${FILE_SRC} ${MAKE_SRC} ${GRAPH_SRC} ${UTILS_SRC})

#target_link_libraries(solution UTILS GRAPH MAKE FILE SOLVE)

#set_target_properties(test PROPERTIES CUDA_ARCHITECTURES "52")