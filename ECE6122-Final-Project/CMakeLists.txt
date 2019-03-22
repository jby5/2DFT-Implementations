cmake_minimum_required(VERSION 3.9.1)

#using modern c++
add_definitions(-std=c++11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -std=c++11 -pthread")

#dft: threads
project(p31)
set(SOURCE_FILES_P31
     src/complex.h
     src/complex.cc
     src/input_image.h
     src/input_image.cc
)
add_executable(p31 ${SOURCE_FILES_P31} src/p31.cc)
target_compile_definitions(p31 PRIVATE -DFFT2=threads_dft2)

#fft : mpi
project(p32)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(src)
set(SOURCE_FILES_P32
     src/complex.h
     src/complex.cc
     src/input_image.h
     src/input_image.cc
)
add_executable(${PROJECT_NAME} ${SOURCE_FILES_P32} src/p32.cpp)
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})


# fft cuda
project(p33)
find_package(CUDA REQUIRED)
SET(CUDA_64_BIT_DEVICE_CODE ON)
SET(CUDA_SEPARABLE_COMPILATION ON)
#add_executable(p32 src/main.cpp src/matrix.cpp src/fft_cuda.cu)
SET(CUDA_SRC
        src/p33.cu
        src/complex.cuh
        src/input_image.cuh
        )
CUDA_ADD_EXECUTABLE(
        p33
        src/p33.cu)
target_link_libraries(p33 ${CUDA_LIBRARIES} ${CUDA_CUDART_LIBRARY} ${CUDA_cudadevrt_LIBRARY})
