cmake_minimum_required(VERSION 3.26)

project(biort)


add_executable(
    biort 
    main.c
    custom_io.c
    init.c
    lookup.c 
    optparse.c 
    print.c 
    react.c 
    read_chem.c 
    read_cini.c 
    read_hbv.c 
    read_param.c 
    read_precipchem.c 
    read_soil.c 
    set_numexp.c 
    speciation.c 
    time_func.c 
    transpt.c 
    util_func.c
)

add_dependencies(biort sundials_cvode_static)
add_dependencies(biort sundials_nvecserial_static)

if(WIN32)
    target_include_directories(
        biort PRIVATE 
        ${BioRT_SOURCE_DIR}/build/cvode/include
        ${BioRT_SOURCE_DIR}/cvode/include
        ${BioRT_SOURCE_DIR}/src/include
    )
else()
    target_include_directories(
        biort PRIVATE 
        ${BioRT_SOURCE_DIR}/build/cvode/include
        ${BioRT_SOURCE_DIR}/cvode/include
        ${BioRT_SOURCE_DIR}/src/include
    )
endif()


if(WIN32)
target_link_libraries(
    biort PUBLIC 
    sundials_cvode
    sundials_nvecserial
)
else()
target_link_libraries(
    biort PRIVATE 
    m 
    sundials_cvode
    sundials_nvecserial
)
endif()

if(WIN32)
    target_link_directories(biort PUBLIC
        ${BioRT_SOURCE_DIR}/build/cvode/src/cvode/Release
        ${BioRT_SOURCE_DIR}/build/cvode/src/nvector/serial/Release
    )
else()
    target_link_directories(biort PUBLIC
        ${BioRT_SOURCE_DIR}/build/cvode/src/cvode
        ${BioRT_SOURCE_DIR}/build/cvode/src/nvector/serial
    )
endif()


if(NOT WIN32)
    add_compile_options(
        -O2 -fpermissive
    )
endif()

if(NOT WIN32)
    add_dependencies(biort sundials_cvode_static)
    add_dependencies(biort sundials_nvecserial_static)
endif()

add_custom_command(TARGET biort POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:biort> ${BioRT_SOURCE_DIR})