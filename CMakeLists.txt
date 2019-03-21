cmake_minimum_required(VERSION 2.8.12)

project(steinhardt)



# Set source directory
set(SOURCE_DIR "src/steinhardt")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})


set(SOURCES "${SOURCE_DIR}/steinhardt.cpp")

add_subdirectory(lib/pybind11)
pybind11_add_module(steinhardt ${SOURCES} "${SOURCE_DIR}/steinhardt_binding.cpp")
