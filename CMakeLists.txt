cmake_minimum_required(VERSION 3.22)
project(MHD2D LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/InputData)
include_directories(${CMAKE_SOURCE_DIR}/OutputData)
include_directories(${CMAKE_SOURCE_DIR}/src/include)
include_directories(${CMAKE_SOURCE_DIR}/src/geometry)

# Find OpenMP
find_package(OpenMP)

# Main executable
add_executable(MHD2D
    src/MHDsolver/MHDSolver1D.cpp
    src/MHDsolver/MHDSolver2D.cpp
    src/MHDsolver/solver.cpp
    src/include/json/jsoncpp.cpp
    src/geometry/NetGeometry.cpp
    src/include/utility/utility.cpp
)

# Geometry executable
add_executable(Geometry
    src/geometry/worldGen.cpp
    src/geometry/NetGeometry.cpp
    src/include/json/jsoncpp.cpp
    src/include/utility/utility.cpp
        src/include/json/json.h
        src/include/json/jsoncpp.cpp
        src/include/service/FileIO.h
        src/include/service/LinOp.h
)

# Include dirs for specific targets
target_include_directories(MHD2D PRIVATE
    ${CMAKE_SOURCE_DIR}/src/json
)

# OpenMP linking
if(OpenMP_CXX_FOUND)
    target_link_libraries(MHD2D PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(Geometry PRIVATE OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found, continuing without OpenMP support")
endif()

option(USE_CUDA "Enable CUDA support" OFF)

if(USE_CUDA)
    enable_language(CUDA)
    target_sources(MHD2D PRIVATE
        src/MHDsolver/MHDgpu.cu
        src/MHDsolver/GPUFunctions.cu
    )
    target_compile_definitions(MHD2D PRIVATE USE_CUDA)
endif()
