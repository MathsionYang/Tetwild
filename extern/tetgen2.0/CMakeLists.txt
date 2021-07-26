cmake_minimum_required (VERSION 2.6)

project (TETGEN2)

option(BUILD_STATIC_LIBRARY "Build library." ON)

# The version number.
set (TETGEN_VERSION_MAJOR 2)
set (TETGEN_VERSION_MINOR 0)
set (TETGEN_PATCH_VERSION 0)

set (TETGEN_VERSION "${TETGEN_MAJOR_VERSION}.${TETGEN_MINOR_VERSION}")

add_definitions(-DUSING_GMP)

if(BUILD_STATIC_LIBRARY)
    add_definitions(-DBUILD_LIBRARY)
endif()

set(SOURCES 
tetgen/tetra.cpp
tetgen/io.cpp
tetgen/pred3d.cpp
tetgen/metric.cpp
tetgen/sort.cpp
tetgen/delaunay.cpp
tetgen/flips.cpp
tetgen/constrained.cpp
tetgen/voronoi.cpp
tetgen/shelling.cpp
tetgen/hts.cpp
tetgen/main.cpp)

if(BUILD_STATIC_LIBRARY)
    add_library(tetgen2 STATIC ${SOURCES})
    target_include_directories(tetgen2 PUBLIC .)
else()
    add_executable(tetgen2 ${SOURCES})
endif()

target_link_libraries(tetgen2 gmp gmpxx)

if(BUILD_STATIC_LIBRARY)
    add_library(tetgen::tetgen ALIAS tetgen2)
endif()
