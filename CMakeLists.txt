cmake_minimum_required (VERSION 2.6)

project(ComputationalFluidDynamics)

## add debug flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-Wall -Wextra -Wshadow -Wconversion)
endif()

######################################################
## add C++11 flag
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
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
                     ${PROJECT_SOURCE_DIR}/input/*.dat)

foreach(InputFile ${InputFiles})
  configure_file("${InputFile}" "./" COPYONLY)
endforeach()

add_executable(blatt1 blatt1.cpp)
target_link_libraries(blatt1 CFD)