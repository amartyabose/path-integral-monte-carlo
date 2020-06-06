cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(pimc)

set(pimc_VERSION_MAJOR 0)
set(pimc_VERSION_MINOR 1)

set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/*.cpp)

add_executable(pimc ${SOURCES})
target_include_directories(pimc PUBLIC ${PROJECT_SOURCE_DIR})
target_compile_features(pimc PUBLIC cxx_std_11)

set(default_build_type "Release")

message("Creating build files for ${default_build_type}")

target_compile_definitions(pimc PUBLIC -DSPDLOG_COMPILED_LIB)

set(Boost_NO_BOOST_CMAKE ON)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS system filesystem REQUIRED)

find_package(Armadillo REQUIRED)
target_compile_definitions(pimc PUBLIC -DARMA_NO_DEBUG -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_HDF5)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED)

set(SANITIZERS "-fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak")
set(CMAKE_CXX_FLAGS_DEBUG "$-O0 -g3 ${SANITIZERS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mtune=native -march=native")

include_directories(${Boost_INCLUDE_DIRS} ${Armadillo_INCLUDE_DIRS})

target_link_libraries(pimc
    ${Boost_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    )

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(pimc ${MPI_LIBRARIES})
