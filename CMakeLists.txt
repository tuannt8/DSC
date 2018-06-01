cmake_minimum_required(VERSION 3.6)
project(DSC)

set(CMAKE_CXX_STANDARD 14)

aux_source_directory(CGLA SRC_CGLA)
aux_source_directory(SOIL SRC_SOIL)
aux_source_directory(src SRC_DSC)
aux_source_directory(is_mesh SRC_MESH)

include_directories(src)
include_directories(is_mesh)
include_directories(./)

add_library(CGLA ${SRC_CGLA})
add_library(SOIL ${SRC_SOIL})
add_library(DSC ${SRC_DSC} ${SRC_MESH})
target_link_libraries(DSC CGLA)
