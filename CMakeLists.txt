cmake_minimum_required(VERSION 3.22)
project(MHD2D)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)

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
    src/include/json/jsoncpp.cpp
    src/geometry/NetGeometry.cpp
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

