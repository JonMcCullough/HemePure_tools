cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()

project(HemePure_Tools VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(ZLIB REQUIRED)
find_package(MPI)

add_executable(gmy2gmy+ gmy2gmy+/gmy2gmy+.cc)
target_link_libraries(gmy2gmy+ ZLIB::ZLIB)

add_executable(gmy2lets gmy2inlets/gmy2lets.cc)
target_link_libraries(gmy2lets PRIVATE ZLIB::ZLIB)
target_include_directories(gmy2lets PRIVATE ZLIB::ZLIB)

add_executable(vx2ascii vx2ascii/vx2ascii.cc)

if (MPI_FOUND)
    add_executable(mpivx2gmy vx2gmy/mpivx2gmy.cc)
    target_include_directories(mpivx2gmy PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(mpivx2gmy PRIVATE MPI::MPI_CXX)
    target_link_libraries(mpivx2gmy PRIVATE ZLIB::ZLIB)
endif()

add_executable(vx2paraview vx2paraview/vx2paraview.cc)

add_executable(vx2refine vx2refine/vx2refine.cc)
target_compile_definitions(vx2refine 
    PUBLIC LARGEFILE_SOURCE=1 FILE_OFFSET_BITS=64
)

if (MPI_FOUND)
    add_executable(mpivx2refine vx2refine/parallel/mpivx2refine.cc)
    target_compile_definitions(mpivx2refine 
        PUBLIC LARGEFILE_SOURCE=1 FILE_OFFSET_BITS=64
    )
    target_include_directories(mpivx2refine PRIVATE MPI::MPI_CXX)
    target_link_libraries(mpivx2refine PRIVATE MPI::MPI_CXX)
    target_link_libraries(mpivx2refine PRIVATE ZLIB::ZLIB)
endif()

option(BUILD_VOXELIZER "Build the Palabos-based voxelizer" OFF)

if (BUILD_VOXELIZER)

    include(FetchContent)
    FetchContent_Declare(
        palabos
        URL https://gitlab.com/unigespc/palabos/-/archive/v2.3.0/palabos-v2.3.0.tar.gz
        URL_HASH SHA256=2085de7b06cc9c4a7b3457f7c2d17747d8960c0f861e32fc883acd504dfc1e23
    )
    # Only the palabos headers appear to be needed!
    FetchContent_Populate(palabos)

    #add_library(tinyxml
    #    ${palabos_SOURCE_DIR}/externalLibraries/tinyxml/tinyxml.cpp
    #    ${palabos_SOURCE_DIR}/externalLibraries/tinyxml/tinyxmlerror.cpp
    #    ${palabos_SOURCE_DIR}/externalLibraries/tinyxml/tinyxmlparser.cpp
    #)
    #target_include_directories(tinyxml PUBLIC
    #    ${palabos_SOURCE_DIR}/externalLibraries/tinyxml)

    add_executable(voxelizer voxelizer/source/voxelizer_MultiInput.cpp)
    target_include_directories(voxelizer
        PRIVATE ${palabos_SOURCE_DIR}/src
        PRIVATE ${palabos_SOURCE_DIR}/externalLibraries)

endif()

#
# Installation
#
install(TARGETS gmy2gmy+ gmy2lets vx2ascii vx2paraview vx2refine RUNTIME)

if(MPI_FOUND)
    install(TARGETS mpivx2gmy mpivx2refine RUNTIME)
endif()

if(BUILD_VOXELIZER)
    install(TARGETS voxelizer RUNTIME)
endif()