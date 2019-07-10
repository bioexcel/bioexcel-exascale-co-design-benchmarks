# CMake 3.9.6 is probably a bit newer than technically required, so
# try relaxing this if you need to!
cmake_minimum_required(VERSION 3.9.6)

project(GromacsKernelMicrobenchmark)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Create the kernel benchmark target, along with dependencies
add_executable(
    gromacs_kernel_microbenchmark
    main.cpp
    kernel_ref.cpp
    )
target_include_directories(
    gromacs_kernel_microbenchmar
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    )
