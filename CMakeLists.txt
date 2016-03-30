cmake_minimum_required(VERSION 2.8.9)

# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

project(Fnm C CXX)

set(FNM_MAJOR_VERSION 0)
set(FNM_MINOR_VERSION 7)
set(FNM_PATCH_VERSION 0)
set(FNM_VERSION
  ${FNM_MAJOR_VERSION}.${FNM_MINOR_VERSION}.${FNM_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/Fnm)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# JEM: Allow the developer to select if Dynamic or Static libraries are built
option (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
set(LIB_TYPE STATIC)
if (BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()

option (BUILD_TEST "Build tests" OFF)

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

#include(cmake/gtest.cmake)

#include(cmake/gmock.cmake)

#include(cmake/nanomsg.cmake)

# set up include-directories
include_directories(
  "${PROJECT_SOURCE_DIR}"   # to find fnm/exports_fnm.h
  "${PROJECT_BINARY_DIR}")  # to find fnm/config.h

enable_testing()

# Add sub-directories
add_subdirectory(sps)

add_subdirectory(fnm)

#add_subdirectory(fftw_test)

#add_subdirectory(irq)

#add_subdirectory(msq)

#add_subdirectory(doc)

# add_subdirectory(sps)

# The interesting stuff goes here
# ===============================

# Add all targets to the build-tree export set
export(TARGETS fnm
  FILE "${PROJECT_BINARY_DIR}/FnmTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Fnm)

# Create the FnmConfig.cmake and FnmConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

# JEM: Generation of FnmConfig.cmake needs to know if dynamic linking is used
if (BUILD_SHARED_LIBS)
  set(FNM_DYNAMIC_LINKING TRUE)
else()
  set(FNM_DYNAMIC_LINKING FALSE)
endif()

# JEM: This is used for installations - no hard coding of directories,
file(RELATIVE_PATH REL_BIN_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_BIN_DIR}")

# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "") # For build tree, multiple configurations are used
configure_file(FnmConfig.cmake.in
  "${PROJECT_BINARY_DIR}/FnmConfig.cmake" @ONLY)

# FNM_CMAKE_DIR is just a variable name, used in FnmConfig.cmake.in

# ... for the install tree
set(CONF_INCLUDE_DIRS "\${FNM_CMAKE_DIR}/${REL_INCLUDE_DIR}")
set(CONF_LIBRARY_DIRS "\${FNM_CMAKE_DIR}/${REL_BIN_DIR}")
configure_file(FnmConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FnmConfig.cmake" @ONLY)

# ... for both
configure_file(FnmConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/FnmConfigVersion.cmake" @ONLY)

# Install the FnmConfig.cmake and FnmConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FnmConfig.cmake"
  "${PROJECT_BINARY_DIR}/FnmConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT FnmTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Add uninstall target, in linux simply issue xargs rm < install_manifest.txt
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(UNINSTALL
    COMMAND ${CMAKE_COMMAND} -P
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)
