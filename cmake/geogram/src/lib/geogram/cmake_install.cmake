# Install script for directory: /home/gridteam/workplace/TetWild/extern/geogram/src/lib/geogram

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevkitx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1.6.7"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/usr/local/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so.1.6.7"
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so.1"
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1.6.7"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/gridteam/workplace/TetWild/extern/geogram/cmake/lib:"
           NEW_RPATH "/usr/local/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevkit-fullx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1.6.7"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/usr/local/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so.1.6.7"
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so.1"
    "/home/gridteam/workplace/TetWild/cmake/lib/libgeogram.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1.6.7"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/gridteam/workplace/TetWild/extern/geogram/cmake/lib:"
           NEW_RPATH "/usr/local/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevkitx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram" TYPE DIRECTORY FILES "/home/gridteam/workplace/TetWild/extern/geogram/src/lib/geogram/api" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevkit-fullx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram" TYPE DIRECTORY FILES "/home/gridteam/workplace/TetWild/extern/geogram/src/lib/geogram/." FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/license/" EXCLUDE)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/gridteam/workplace/TetWild/cmake/geogram/geogram1.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/gridteam/workplace/TetWild/cmake/geogram/src/lib/geogram/third_party/cmake_install.cmake")

endif()

