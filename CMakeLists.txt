cmake_minimum_required(VERSION 3.26)

project(BioRT LANGUAGES C)

if (WIN32)
add_compile_options(/MT /O2)
endif()

add_subdirectory(cvode)
add_subdirectory(src)