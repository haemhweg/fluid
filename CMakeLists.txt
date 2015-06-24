cmake_minimum_required (VERSION 2.6)

project(ComputationalFluidDynamics)

## add debug flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-Wall -Wextra -Wshadow -Wconversion)
endif()

######################################################
## include MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

######################################################
## add C++11 flag
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

set(SRCPATH "${PROJECT_SOURCE_DIR}/src")
set(INCPATH "${PROJECT_SOURCE_DIR}/include")

## add to library
file(GLOB SRC ${SRCPATH}/*.cpp)
file(GLOB INCLUDE ${INCPATH}/*.h)

## add includes
include_directories("${INCPATH}")

add_library(CFD ${INCLUDE} ${SRC})


##########################################################
## executables

## add input

file(GLOB InputFiles ${PROJECT_SOURCE_DIR}/input/*.vtk
                     ${PROJECT_SOURCE_DIR}/input/*.dat
		     ${PROJECT_SOURCE_DIR}/input/config*)

foreach(InputFile ${InputFiles})
  configure_file("${InputFile}" "./" COPYONLY)
endforeach()

add_executable(blatt3 blatt3.cpp)
target_link_libraries(blatt3 CFD ${MPI_LIBRARIES})