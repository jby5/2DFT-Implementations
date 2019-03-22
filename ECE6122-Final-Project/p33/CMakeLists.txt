cmake_minimum_required(VERSION 3.9.1)

PROJECT(p33)
FIND_PACKAGE(CUDA QUIET REQUIRED)
SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -rdc=true; -std=c++11;)
SET(CUDA_64_BIT_DEVICE_CODE ON)
SET(CUDA_SEPARABLE_COMPILATION ON)
INCLUDE_DIRECTORIES( src)
SET(CUDA_SRC
        src/input_image.cuh
        src/complexCUDA.cuh
        src/complex.cu
        )
LINK_DIRECTORIES()
CUDA_ADD_EXECUTABLE(
        p33 ${CUDA_SRC} src/main.cu
)
target_link_libraries(p33 ${CUDA_LIBRARIES})
