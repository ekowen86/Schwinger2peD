cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the project name and version
project(Schwinger2peD VERSION 1.0.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(ExternalProject)

# Link OpempMP CMake automatically adds 
# _OPENMP to the CXX flags
option(ENABLE_OPENMP "enable OpenMP" ON)
if(ENABLE_OPENMP)
  find_package(OpenMP)
endif()

add_subdirectory(lib)
add_subdirectory(wilson)
