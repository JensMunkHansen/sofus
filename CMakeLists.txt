set(FNM_CMAKE_MINIMUM_REQUIRED_VERSION 3.7)
set(FNM_CMAKE_MAXIMUM_TESTED_VERSION 3.18)

cmake_minimum_required(VERSION ${FNM_CMAKE_MINIMUM_REQUIRED_VERSION})

if(${CMAKE_VERSION} VERSION_LESS FNM_CMAKE_MAXIMUM_TESTED_VERSION)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION ${FNM_CMAKE_MAXIMUM_TESTED_VERSION})
endif()

# CMP0078: You can switch between LEGACY
# (SWIG_MODULE_<name>_REAL_NAME) or NEW (targets prefix'ed with _) by
# setting UseSWIG_TARGET_NAME_PREFERENCE. If CMP0078=NEW actual
# targets name equal swig target.

# CMP0086: NEW target names can be set using -module flag to the swig
# executable. To use this define SWIG_TARGET_NAME property on .i files

set(project_policies CMP0078 CMP0086)

foreach(policy ${project_policies})
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

project(fnm C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Prefix used for dependencies
set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}/deps-prefix" CACHE STRING "Prefix path for dependencies")

find_package(Doxygen 1.9 QUIET)

include(SetInstallDir)

# Uncommented for Windows Visual Studio 2019
include(DoxygenTarget)

include(SetWarningLevels)

include(DumpVariables)

string(TIMESTAMP timestring %Y%m%d%H%M%S UTC)
set(FNM_MAJOR_VERSION 1)
set(FNM_MINOR_VERSION 0)
set(FNM_PATCH_VERSION 0)
set(FNM_TWEAK_VERSION ${timestring})

set(FNM_VERSION
  ${FNM_MAJOR_VERSION}.${FNM_MINOR_VERSION}.${FNM_PATCH_VERSION}.${FNM_TWEAK_VERSION})

# Setup installation directories and make relative paths absolute 
# ===============================================================
setup_installation_directories()

# JEM: Allow the developer to select if Dynamic or Static libraries are built
option (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
set(LIB_TYPE STATIC)
if (BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()

set(USE_PROGRESS_BAR 0)
option(USE_ProgressBar "Use progress bar" ON)
if (USE_ProgressBar)
  set(USE_PROGRESS_BAR 1)
endif()

option (BUILD_FNM_TEST "Build tests" OFF)

option(USE_MQUEUE "Use message queues" OFF)

option (ENABLE_PYTHON3 "Python 3" ON)

set (_DEBUG_USING_PYTHON_RELEASE_RUNTIME 0)

if (MSVC)
  option(SOFUS_DEBUG_USING_PYTHON_RELEASE_RUNTIME "Debug using Python release runtime" ON)
endif()

if (SOFUS_DEBUG_USING_PYTHON_RELEASE_RUNTIME)
  set (_DEBUG_USING_PYTHON_RELEASE_RUNTIME 1)
endif()

option(BUILD_DOCUMENTATION
    "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if (BUILD_DOCUMENTATION)
  find_program(SED_CMD NAMES sed)
  if (NOT EXISTS "${SED_CMD}")
    message( FATAL_ERROR "sed command is required for building this library")
  endif()
endif()

# set up include-directories
include_directories(
  "${PROJECT_SOURCE_DIR}"   # to find fnm/exports_fnm.h
  "${PROJECT_BINARY_DIR}"   # to find fnm/config.h
  )

if (BUILD_SHARED_LIBS)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif()
endif()

if (BUILD_DOCUMENTATION)
  if(DOXYGEN_FOUND)
    add_doxygen_target(doc)
  else()
    message(FATAL_ERROR "Doxygen is required.")
  endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -msse2 -msse3 -msse4 -mavx -pedantic -Wall")
endif()

# The intention is to include this at a later time
set(FNM_PULSED_WAVE 0)
option(SOFUS_Pulsed_waves "Support pulsed waves" OFF)
if (SOFUS_Pulsed_waves)
  set(FNM_PULSED_WAVE 1)
endif()

set(USE_PROGRESS_BAR 0)
option(USE_ProgressBar "Use progress bar" ON)
if (USE_ProgressBar)
  set(USE_PROGRESS_BAR 1)
endif()

if (1)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(LDFLAGS "${LDFLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

# Add sub-directories
add_subdirectory(sps)

add_subdirectory(gl)

add_subdirectory(sofus)

add_subdirectory(fnm)

add_subdirectory(python)

if (BUILD_FNM_TEST)
  enable_testing()
  add_subdirectory(test)
endif()

# The interesting stuff goes here
# ===============================

# Add all targets to the build-tree export set
export(TARGETS fnm sofus gl sps
  FILE "${PROJECT_BINARY_DIR}/FnmTargets.cmake")

# Export the package for use from the build-tree
export(PACKAGE Fnm)

# Create the FnmConfig.cmake and FnmConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

if (WIN32)
  # Expose location generated DLL's (used on Windows)
  if (BUILD_SHARED_LIBS)
    set(FNM_DYNAMIC_LINKING TRUE)
  else()
    set(FNM_DYNAMIC_LINKING FALSE)
  endif()
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
