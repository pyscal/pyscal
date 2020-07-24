cmake_minimum_required(VERSION 2.8.12)

project(pyscal)
#------------------------------------------------------

# Set source directory
#------------------------------------------------------
# main directory for steinhardt cpp modules
set(SOURCE_DIR "src/pyscal")

#sub directory for voro lib
set(SSOURCE_DIR "lib/voro++")


#Headers
#-------------------------------------------------------
# Tell CMake that headers are also in SOURCE_DIR
include_directories("${PROJECT_SOURCE_DIR}/lib/voro++")
include_directories(${SOURCE_DIR})
set(ATOMSOURCES "${SOURCE_DIR}/atom.cpp")
set(SYSSOURCES "${SOURCE_DIR}/system.cpp")

#Now add subdir for the modules
#--------------------------------------------------------
add_subdirectory(lib/pybind11)
#message(STATUS "include dir: ${CMAKE_BUILD_TYPE}")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -DNDEBUG")
#set(CXX_STANDARD 11)

#Finally create the package
#--------------------------------------------------------
pybind11_add_module(catom ${ATOMSOURCES} "${SOURCE_DIR}/atom_binding.cpp")
pybind11_add_module(csystem ${SYSSOURCES} "${SOURCE_DIR}/system_binding.cpp" "${SOURCE_DIR}/atom.cpp" "${SSOURCE_DIR}/voro++.cc")
