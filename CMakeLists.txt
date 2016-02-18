cmake_minimum_required(VERSION 2.8)

####
# project settings
####
enable_language(Fortran)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)

#set(CMAKE_Fortran_COMPILER "gfortran") 
#set(MY_FFLAGS_DOUBLE $ENV{FFLAGS_DOUBLE})
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

project(haloExtractor C CXX Fortran)
#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
include_directories(include)

####
# sources
####

file(GLOB SOURCES "src/*")

####
# dependencies
####

#include_directories("/usr/local/")
#find_package(CCfits REQUIRED)
#link_directories("/usr/local/lib/")

####
# omp
####
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

#add_compile_options( "-mcmodel=medium" )
####
# build
####
#add_definitions(-std=gnu++0x -L/usr/local/lib/CCfits/.libs -lCCfits -lcfitsio)
set (CMAKE_CXX_FLAGS "-std=gnu++0x -lCCfits -lcfitsio")


add_executable(halo_extractor ${SOURCES})

target_link_libraries(halo_extractor "/usr/local/lib/libCCfits.so")

#target_link_libraries(ulator ${GLAMER_LIBRARIES})
