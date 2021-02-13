cmake_minimum_required(VERSION 3.5)
project(RegridIC)
enable_language (Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(ESMF REQUIRED)

set(CMAKE_Fortran_COMPILER  ${ESMF_F90COMPILER})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEPATHS} ${ESMF_F90LINKPATHS} ${ESMF_F90ESMFLINKRPATHS}")

set(TARGET regridic.exe)

set(SOURCE_FILES src/datetime_mod.F90 src/main.F90)

add_executable(${TARGET} ${SOURCE_FILES})

target_link_libraries(${TARGET} ${ESMF_F90ESMFLINKLIBS})

set_target_properties(${TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
