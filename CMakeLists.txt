cmake_minimum_required(VERSION 3.26)

project(BioRT LANGUAGES C)

if (WIN32)
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        message("Compiling on MSys2")
    else()
        add_compile_options(/MT /O2)
    endif()
endif()

add_subdirectory(cvode)
add_subdirectory(src)