# Top level CMakeLists.txt for AutoGrid.
# Install Boost library before using this project.

cmake_minimum_required (VERSION 3.6)
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project(AutoGrid)
set (CMAKE_CXX_STANDARD 14)

# Basic options
option (TETWILD_WITH_TETGEN         "Boundary recovery with tetgen" ON  )
option (EXPAND_CMAKE_MESSAGE		"Expand output details"			OFF	)
option (AUTOSURF_WITH_ISPC			"Use ISPC"						OFF	)
option (LIBIGL_SKIP_DOWNLOAD		"Skip download external lib"	ON	)
option (LIBIGL_USE_STATIC_LIBRARY	"Use libigl as static library"	OFF	)
option (LIBIGL_WITH_ANTTWEAKBAR		"Use AntTweakBar"				OFF	)
option (LIBIGL_WITH_CGAL			"Use CGAL"						ON	)
option (LIBIGL_WITH_COMISO			"Use CoMiso"					OFF	)
option (LIBIGL_WITH_CORK			"Use Cork"						OFF	)
option (LIBIGL_WITH_LIM				"Use LIM"						OFF	)
option (LIBIGL_WITH_MATLAB			"Use Matlab"					OFF	)
option (LIBIGL_WITH_MOSEK			"Use MOSEK"						OFF	)
option (LIBIGL_WITH_PNG				"Use PNG"						OFF	)
option (LIBIGL_WITH_PYTHON			"Use Python"					OFF	)
option (LIBIGL_WITH_TETGEN			"Use Tetgen"					OFF	)
option (LIBIGL_WITH_TRIANGLE		"Use Triangle"					OFF	)
option (LIBIGL_WITH_XML				"Use XML"						OFF	)
option (LIBIGL_WITH_EMBREE			"Use Embree"					OFF	)
option (LIBIGL_WITH_NANOGUI			"Use Nanogui menu"				OFF	)
option (LIBIGL_WITH_OPENGL			"Use OpenGL"					OFF	)
option (LIBIGL_WITH_OPENGL_GLFW		"Use GLFW"						OFF	)
option (LIBIGL_WITH_VIEWER			"Use OpenGL viewer"				OFF	)
option (GEOGRAM_WITH_TRIANGLE		"Use Triangle"					OFF	)

# Build type
if (NOT CMAKE_BUILD_TYPE)
	message (STATUS "No build type selected, default to Release")
	set (CMAKE_BUILD_TYPE "Release")
endif ()

# Use folder in Visual Studio
if (WIN32)
	set_property (GLOBAL PROPERTY USE_FOLDERS ON)
endif ()

# Extra warnings and color output
include (UseColors)
include (Warnings)

# Unknown property
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

# Collect 3rd party dependencies
set (EXTERNAL_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/extern)
set (EIGEN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen)
include (geogram)
find_package (Boost REQUIRED)
find_package (LIBIGL REQUIRED)
add_subdirectory (${EXTERNAL_DIRS}/pymesh)
add_subdirectory (${EXTERNAL_DIRS}/fmt)
add_library (spdlog INTERFACE)
add_library (spdlog::spdlog ALIAS spdlog)
target_include_directories (spdlog INTERFACE ${EXTERNAL_DIRS}/spdlog/include)
target_compile_definitions (spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
target_link_libraries (spdlog INTERFACE fmt::fmt)

# Project library source directories
set (HEAD_DIR ${PROJECT_SOURCE_DIR}/include/tetwild)
set (CPP_DIR ${PROJECT_SOURCE_DIR}/src/tetwild)
file (GLOB_RECURSE HEAD_FILES LIST_DIRECTORIES false ${HEAD_DIR}/*.*)
file (GLOB_RECURSE CPP_FILES LIST_DIRECTORIES false ${CPP_DIR}/*.* )
if (EXPAND_CMAKE_MESSAGE)
	message (STATUS "Source Files:")
	foreach (FILE IN LISTS HEAD_FILES)
		message (STATUS ${FILE})
	endforeach (FILE)
	foreach (FILE IN LISTS CPP_FILES)
		message (STATUS ${FILE})
	endforeach (FILE)
endif ()

add_subdirectory(${EXTERNAL_DIRS}/tetgen2.0)
if(TETWILD_WITH_TETGEN)
	add_definitions(-DUSE_TETGEN)
endif()
add_library (libAutoGrid STATIC ${HEAD_FILES} ${CPP_FILES})
target_include_directories (libAutoGrid PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories (libAutoGrid PRIVATE ${PROJECT_SOURCE_DIR}/include/tetwild)
target_include_directories (libAutoGrid PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories (libAutoGrid PUBLIC ${EXTERNAL_DIRS})
target_include_directories (libAutoGrid PUBLIC ${EIGEN_ROOT})
target_include_directories (libAutoGrid PUBLIC ${EXTERNAL_DIRS}/tetgen2.0)
if (Boost_FOUND)
	target_include_directories (libAutoGrid PRIVATE ${Boost_INLCUDE_DIRS})
endif ()
target_link_libraries (libAutoGrid PUBLIC geogram)
target_link_libraries (libAutoGrid PUBLIC igl::core)
target_link_libraries (libAutoGrid PUBLIC pymesh::pymesh)
target_link_libraries (libAutoGrid PUBLIC spdlog::spdlog)
target_link_libraries (libAutoGrid PUBLIC igl::cgal)
target_link_libraries (libAutoGrid PRIVATE warnings::all)
if(TETWILD_WITH_TETGEN)
    target_link_libraries (libAutoGrid PUBLIC tetgen::tetgen)
endif()

# Use ispc
if (AUTOSURF_WITH_ISPC)
	message (STATUS "Compiling energy with ISPC")
	add_subdirectory (src/ispc)
	ispc_add_energy (libAutoGrid)
endif ()

# Building executable
add_executable (AutoGrid src/main.cpp)
target_include_directories (AutoGrid PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories (AutoGrid PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries (AutoGrid libAutoGrid)
target_link_libraries (AutoGrid igl::cgal)
target_link_libraries (AutoGrid warnings::all)

# Install
install (TARGETS AutoGrid RUNTIME DESTINATION bin)
