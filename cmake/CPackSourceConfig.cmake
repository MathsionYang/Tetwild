# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


SET(CPACK_ARCHIVE_COMPONENT_INSTALL "ON")
SET(CPACK_BINARY_7Z "")
SET(CPACK_BINARY_BUNDLE "")
SET(CPACK_BINARY_CYGWIN "")
SET(CPACK_BINARY_DEB "")
SET(CPACK_BINARY_DRAGNDROP "")
SET(CPACK_BINARY_FREEBSD "")
SET(CPACK_BINARY_IFW "")
SET(CPACK_BINARY_NSIS "")
SET(CPACK_BINARY_OSXX11 "")
SET(CPACK_BINARY_PACKAGEMAKER "")
SET(CPACK_BINARY_PRODUCTBUILD "")
SET(CPACK_BINARY_RPM "")
SET(CPACK_BINARY_STGZ "")
SET(CPACK_BINARY_TBZ2 "")
SET(CPACK_BINARY_TGZ "")
SET(CPACK_BINARY_TXZ "")
SET(CPACK_BINARY_TZ "")
SET(CPACK_BINARY_WIX "")
SET(CPACK_BINARY_ZIP "")
SET(CPACK_BUILD_SOURCE_DIRS "/home/gridteam/workplace/TetWild;/home/gridteam/workplace/TetWild/cmake")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_COMPONENTS_ALL "runtime;devkit;devkit-full;doc-devkit;doc-devkit-full")
SET(CPACK_COMPONENTS_ALL_SET_BY_USER "TRUE")
SET(CPACK_COMPONENTS_GROUPING "IGNORE")
SET(CPACK_COMPONENT_DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit")
SET(CPACK_COMPONENT_DEVKIT-FULL_GROUP "Development")
SET(CPACK_COMPONENT_DEVKIT_DISPLAY_NAME "Vorpaline Developer Kit")
SET(CPACK_COMPONENT_DEVKIT_GROUP "Development")
SET(CPACK_COMPONENT_DOC-DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit Documentation")
SET(CPACK_COMPONENT_DOC-DEVKIT-FULL_GROUP "Documentation")
SET(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_DISPLAY_NAME "Vorpaline Internal Developer Kit Documentation")
SET(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_GROUP "Documentation")
SET(CPACK_COMPONENT_DOC-DEVKIT_DISPLAY_NAME "Vorpaline API Developer Kit Documentation")
SET(CPACK_COMPONENT_DOC-DEVKIT_GROUP "Documentation")
SET(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY "true")
SET(CPACK_COMPONENT_RUNTIME_DISPLAY_NAME "Vorpaline Application")
SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
SET(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
SET(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
SET(CPACK_GENERATOR "TBZ2;TGZ;TXZ;TZ")
SET(CPACK_IGNORE_FILES "/CVS/;/\\.svn/;/\\.bzr/;/\\.hg/;/\\.git/;\\.swp\$;\\.#;/#")
SET(CPACK_INSTALLED_DIRECTORIES "/home/gridteam/workplace/TetWild;/")
SET(CPACK_INSTALL_CMAKE_PROJECTS "")
SET(CPACK_INSTALL_PREFIX "/usr/local")
SET(CPACK_MODULE_PATH "/home/gridteam/workplace/TetWild/cmake")
SET(CPACK_NSIS_DISPLAY_NAME "AutoGrid 1.6.7")
SET(CPACK_NSIS_INSTALLER_ICON_CODE "")
SET(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
SET(CPACK_NSIS_PACKAGE_NAME "AutoGrid 1.6.7")
SET(CPACK_OUTPUT_CONFIG_FILE "/home/gridteam/workplace/TetWild/cmake/CPackConfig.cmake")
SET(CPACK_PACKAGE_DEFAULT_LOCATION "/")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-3.10/Templates/CPack.GenericDescription.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "fast, simple and easy-to-use primitives for geometric programming")
SET(CPACK_PACKAGE_FILE_NAME "AutoGrid-1.6.7-Source")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "AutoGrid 1.6.7")
SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "AutoGrid 1.6.7")
SET(CPACK_PACKAGE_NAME "AutoGrid")
SET(CPACK_PACKAGE_RELOCATABLE "true")
SET(CPACK_PACKAGE_VENDOR "INRIA - ALICE")
SET(CPACK_PACKAGE_VERSION "1.6.7")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "7")
SET(CPACK_RESOURCE_FILE_LICENSE "/usr/share/cmake-3.10/Templates/CPack.GenericLicense.txt")
SET(CPACK_RESOURCE_FILE_README "/usr/share/cmake-3.10/Templates/CPack.GenericDescription.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake-3.10/Templates/CPack.GenericWelcome.txt")
SET(CPACK_RPM_PACKAGE_SOURCES "ON")
SET(CPACK_SET_DESTDIR "OFF")
SET(CPACK_SOURCE_7Z "")
SET(CPACK_SOURCE_CYGWIN "")
SET(CPACK_SOURCE_GENERATOR "TBZ2;TGZ;TXZ;TZ")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;/\\.bzr/;/\\.hg/;/\\.git/;\\.swp\$;\\.#;/#")
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES "/home/gridteam/workplace/TetWild;/")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/gridteam/workplace/TetWild/cmake/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "AutoGrid-1.6.7-Source")
SET(CPACK_SOURCE_RPM "OFF")
SET(CPACK_SOURCE_TBZ2 "ON")
SET(CPACK_SOURCE_TGZ "ON")
SET(CPACK_SOURCE_TOPLEVEL_TAG "Linux64-gcc-Source")
SET(CPACK_SOURCE_TXZ "ON")
SET(CPACK_SOURCE_TZ "ON")
SET(CPACK_SOURCE_ZIP "OFF")
SET(CPACK_STRIP_FILES "")
SET(CPACK_SYSTEM_NAME "Linux64-gcc")
SET(CPACK_TOPLEVEL_TAG "Linux64-gcc-Source")
SET(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/gridteam/workplace/TetWild/cmake/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
