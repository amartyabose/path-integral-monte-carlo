cmake_minimum_required(VERSION 2.6)
project(pimc)

set (pimc_VERSION_MAJOR 0)
set (pimc_VERSION_MINOR 1)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER -Wall -Wextra -Wno-unused-parameter")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

include_directories(${Boost_INCLUDE_DIRS} ${Armadillo_INCLUDE_DIRS})
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/*.cpp)
add_executable(pimc ${SOURCES})
target_link_libraries(pimc ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_C_LIBRARIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(pimc ${MPI_LIBRARIES})
