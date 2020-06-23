cmake_minimum_required (VERSION 2.6)
project(PDC)


find_package(MPI REQUIRED)
if(MPI_FOUND)
message(STATUS "MPI FOUND")
message(STATUS "Using MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
message(STATUS "Using MPI_LIBRARIES: ${MPI_CXX_LIBRARIES}")
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${MPI_EXE_LINKER_FLAGS}")
endif(MPI_FOUND)

add_executable(serial src/serial.cc)
target_link_libraries(serial ${MPI_LIBRARIES})
add_executable(parallel src/parallel.cc)
target_link_libraries(parallel ${MPI_LIBRARIES})


