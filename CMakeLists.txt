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
set(SOURCES "${SOURCE_DIR}/steinhardt.cpp")

#Now add subdir for the modules
#--------------------------------------------------------
add_subdirectory(lib/pybind11)

#Finally create the package
#--------------------------------------------------------
pybind11_add_module(ccore ${SOURCES} "${SOURCE_DIR}/steinhardt_binding.cpp" "${SSOURCE_DIR}/voro++.cc")