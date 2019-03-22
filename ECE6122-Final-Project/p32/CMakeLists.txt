cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 11)
project(p32)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(src)
add_executable(${PROJECT_NAME} ${SOURCE_FILES} src/2dfft_mpi.cpp)
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})