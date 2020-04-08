cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(pimc)
set(CMAKE_VERBOSE_MAKEFILE ON)


set (pimc_VERSION_MAJOR 0)
set (pimc_VERSION_MINOR 1)

set (CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set (PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set (default_build_type RELEASE)
# set (CMAKE_BUILD_TYPE Release)

message ("Creating build files for ${default_build_type}")

set(Boost_NO_BOOST_CMAKE ON)

set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package (Boost COMPONENTS system filesystem REQUIRED)
find_package (Armadillo REQUIRED)
find_package (BLAS REQUIRED)
find_package (LAPACK REQUIRED)
find_package (Threads)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PROJECT_SOURCE_DIR} -std=c++11 -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_HDF5 -DBOOST_ALLOW_DEPRECATED_HEADERS -DSPDLOG_COMPILED_LIB")
set (SANITIZERS "-fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g3 ${SANITIZERS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

include_directories (${Boost_INCLUDE_DIRS} ${Armadillo_INCLUDE_DIRS})
file (GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/*.cpp)
add_executable (pimc ${SOURCES})

target_link_libraries (pimc
    ${Boost_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    )

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(pimc ${MPI_LIBRARIES})
# enable_testing()

#MACRO(SUBDIRLIST result curdir)
#  FILE(GLOB children ${curdir}/*)
#  SET(dirlist "")
#  FOREACH(child ${children})
#    IF(IS_DIRECTORY ${child})
#      LIST(APPEND dirlist ${child})
#    ENDIF()
#  ENDFOREACH()
#  SET(${result} ${dirlist})
#ENDMACRO()

#SUBDIRLIST(SUBDIRS "regtests")
#FOREACH(subdir ${SUBDIRS})
#	GET_FILENAME_COMPONENT(dname ${subdir} NAME)
#	message("Adding regtest directory: ${dname}")
#	add_test(NAME "${dname}_compile_test"
#			 COMMAND ${CMAKE_BINARY_DIR}/qcpi in.xml --regtest
#			 WORKING_DIRECTORY ${subdir})
#	add_test(NAME "${dname}_diff_test"
#			 COMMAND ${subdir}/compare.sh
#			 WORKING_DIRECTORY ${subdir})
#ENDFOREACH()
