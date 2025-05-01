cmake_minimum_required(VERSION 3.22)
project(MHD2D LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)

# Set CUDA paths explicitly for Arch Linux's non-standard location
set(CMAKE_CUDA_ARCHITECTURES "native")
set(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES "/opt/cuda/targets/x86_64-linux/include")

# Find CUDA (this will find cuSPARSE automatically)
find_package(CUDAToolkit REQUIRED)

# in/out-put data folders include
include_directories(${CMAKE_SOURCE_DIR}/InputData)
include_directories(${CMAKE_SOURCE_DIR}/OutputData)

# Add include directory for local headers (including Eigen)
include_directories(${CMAKE_SOURCE_DIR}/src/include)
include_directories(${CMAKE_SOURCE_DIR}/src/geometry)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Main executable (solver)
add_executable(MHD2D
    src/MHDsolver/MHDSolver1D.cpp
    src/MHDsolver/MHDSolver2D.cpp
    src/MHDsolver/solver.cpp
    src/MHDsolver/MHDgpu.cu
    src/include/json/jsoncpp.cpp
    src/geometry/NetGeometry.cpp
)

# Explicit include paths for Arch Linux
target_include_directories(MHD2D PRIVATE
        ${CMAKE_SOURCE_DIR}/src/json
        ${CMAKE_SOURCE_DIR}/src/FluidSim
        "/opt/cuda/targets/x86_64-linux/include"
)

# Additional linker paths for Arch Linux
target_link_directories(MHD2D PRIVATE
        "/opt/cuda/targets/x86_64-linux/lib"
)

# Set CUDA standard and properties
set_target_properties(MHD2D PROPERTIES
        CUDA_STANDARD 17
        CUDA_SEPARABLE_COMPILATION ON
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(MHD2D PRIVATE OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()

target_compile_options(MHD2D PRIVATE -fopenmp)
target_link_libraries(MHD2D PRIVATE -fopenmp)

#  Geometry executable
add_executable(Geometry
    src/geometry/worldGen.cpp
    src/geometry/NetGeometry.cpp
    src/include/json/jsoncpp.cpp
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Geometry PRIVATE OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()
target_compile_options(Geometry PRIVATE -fopenmp)
target_link_libraries(Geometry PRIVATE -fopenmp)

target_compile_options(MHD2D PRIVATE -fopenmp)
target_link_libraries(MHD2D PRIVATE -fopenmp)

