cmake_minimum_required(VERSION 3.9)
project(image_diff)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)


set(GCC_COMPILE_OPTIONS "-Wall -Wextra")
set(GCC_COMPILE_DEBUG_OPTIONS "${GCC_COMPILE_OPTIONS};-ggdb;-O0")
set(GCC_COMPILE_RELEASE_OPTIONS "${GCC_COMPILE_OPTIONS};-O3")

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/image_diff.cpp")
add_executable(image_diff ${SRC})