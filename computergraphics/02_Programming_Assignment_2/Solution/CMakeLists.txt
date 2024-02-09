#########################################
#  AdvancedComputerGraphics-PathTracing #
#---------------------------------------#
# @author IGS                           #
# @date    26.4.2022                    #
#---------------------------------------#
#########################################
cmake_minimum_required( VERSION 3.9 )
project( AdvCg-PathTracing)

message("\n * Cpp Coding for Advanced Computer Graphics")

message(STATUS "${PROJECT_NAME} build environment --")
message(STATUS "Build system:     ${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}")

#########################################
#                Options                #
#########################################


#########################################
#              CMake-Stuff              #
#########################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)


#########################################
#     Build/Find External-Libraries     #
#########################################

add_subdirectory(external/tinyobjloader)


#########################################
#        Global Compiler Options        #
#########################################
set(GCC_COMPILE_OPTIONS "-Wall")
set(GCC_COMPILE_DEBUG_OPTIONS "${GCC_COMPILE_OPTIONS};-ggdb;-O0")
set(GCC_COMPILE_RELEASE_OPTIONS "${GCC_COMPILE_OPTIONS};-O3")

add_compile_options("$<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:DEBUG>>:${GCC_COMPILE_DEBUG_OPTIONS}>")
add_compile_options("$<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:RELEASE>>:${GCC_COMPILE_RELEASE_OPTIONS}>")


set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/PathTracing.cpp" )

add_executable( PathTracing ${SRC} )

target_link_libraries(PathTracing tinyobjloader)

target_compile_features( PathTracing PUBLIC cxx_std_17 )
set_target_properties( PathTracing PROPERTIES CXX_EXTENSIONS OFF )

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT PathTracing)