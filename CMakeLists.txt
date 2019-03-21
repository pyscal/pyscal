cmake_minimum_required(VERSION 2.8.12)

project(steinhardt)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#add_library(voro++ lib/voro++/src/voro++.cpp)

# Set source directory
set(SOURCE_DIR "src/steinhardt")
set(SSOURCE_DIR "lib/voro++")
# Tell CMake that headers are also in SOURCE_DIR
include_directories("${PROJECT_SOURCE_DIR}/lib/voro++")
include_directories(${SOURCE_DIR})
set(SOURCES "${SOURCE_DIR}/steinhardt.cpp")

add_subdirectory(lib/pybind11)
pybind11_add_module(steinhardt ${SOURCES} "${SOURCE_DIR}/steinhardt_binding.cpp" "${SSOURCE_DIR}/voro++.cc")
#target_link_libraries(${PROJECT_NAME} PRIVATE voro++)