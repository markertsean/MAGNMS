cmake_minimum_required(VERSION 2.8)

####
# project settings
####

# Compilers
#set (CMAKE_Fortran_COMPILER "/home/local/intel/composerxe-2011.3.174/bin/intel64/ifort")
set (CMAKE_Fortran_COMPILER "/usr/bin/gfortran")


#set (CMAKE_CXX_COMPILER "/home/users/markert/bin/g++")
#set (CMAKE_C_COMPILER "/home/users/markert/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_C_COMPILER "/usr/bin/gcc")


set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


enable_language(Fortran)

project(haloExtractor C CXX Fortran)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (    CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER}     NAME)

####
# sources
####

include_directories(include)
file(GLOB SOURCES "src/*")

####
# set flags
####


set (CMAKE_CXX_FLAGS "-std=gnu++0x -lgfortran -lm")

set (MY_CONVERT_THING big_endian) #Or little_endian? big_endian?


if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fconvert='${MY_CONVERT_THING}'") 
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fconvert='${MY_CONVERT_THING}'")
  message("-- Using gfortran")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-convert ${MY_CONVERT_THING}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-convert ${MY_CONVERT_THING}")
  message("-- Using ifort")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-fconvert='${MY_CONVERT_THING}'") #Or little_endian?
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fconvert='${MY_CONVERT_THING}'")
  message("-- Using gfortran")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


####
# Generate executable
####

add_executable(halo_extractor ${SOURCES})

target_link_libraries(halo_extractor "/usr/local/lib/libCCfits.so")
