# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# set pybind
set(pybind11_DIR /Users/franciscojaviernicolas/Library/Python/3.7/lib/python/site-packages/pybind11/)
find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})

# Set the output directory for build artifacts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# project name
project(PyCppTestProgram)


link_directories(${ROOT_LIBRARY_DIR})

set(CMAKE_CXX_STANDARD 11)

# includes and srcs directories
include_directories(src)
file(GLOB SOURCES "src/*.cpp")

# Compile the C++ code into a shared library
add_library(mygeek SHARED ${SOURCES})

pybind11_add_module(pybind_wrapper pybind_wrapper.cpp)
target_link_libraries(pybind_wrapper PRIVATE mygeek)

