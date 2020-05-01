cmake_minimum_required(VERSION 3.16)
project(precomputed_FFT)

set(CMAKE_CXX_STANDARD 17)
set(MKL_ROOT /opt/intel/compilers_and_libraries_2020.1.217/linux/mkl)
set(LIBRARY_FILES mkl_core mkl_intel_lp64 mkl_sequential)

add_executable(precomputed_FFT precomputedFFT.cpp)

target_link_directories(precomputed_FFT PRIVATE ${MKL_ROOT}/lib/intel64)
target_link_libraries(precomputed_FFT ${LIBRARY_FILES})