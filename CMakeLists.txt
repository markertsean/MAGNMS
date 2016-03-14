cmake_minimum_required(VERSION 2.8)

####
# project settings
####
enable_language(Fortran)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


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

set (CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
set (CMAKE_CXX_FLAGS "-std=gnu++0x -lCCfits -lcfitsio")

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-fconvert='big_endian'") #Or little_endian?
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fconvert='big_endian'")
  message ("Using gfortran...")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-convert big_endian")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-convert big_endian")
  message ("Using ifort...")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-fconvert='big_endian'") #Or little_endian?
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fconvert='big_endian'")
  message ("Using wierd compiler...")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

message ("CMAKE CXX     COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("CMAKE Fortran COMPILER full path: " ${CMAKE_Fortran_COMPILER})

####
# Generate executable
####

add_executable(halo_extractor ${SOURCES})

target_link_libraries(halo_extractor "/usr/local/lib/libCCfits.so")
