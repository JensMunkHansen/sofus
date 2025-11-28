cmake_minimum_required(VERSION 3.20...3.31)

if(NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Add path for custom modules (shared CMake submodule first, then local cmake folder)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Set version from git tag (or default 1.0.0 with commit count as tweak)
include(spsProjectVersion)
sps_set_project_version(
  DEFAULT_VERSION "1.0.0"
  OUT_VERSION SOFUS_VERSION)

project(Sofus
  VERSION ${SOFUS_VERSION}
  LANGUAGES C CXX)

# Ensure PROJECT_VERSION_TWEAK is defined
if(NOT PROJECT_VERSION_TWEAK)
  set(PROJECT_VERSION_TWEAK 0)
endif()

message(STATUS "Sofus version: ${PROJECT_VERSION} (${SOFUS_VERSION_SOURCE})")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# CMP0078: SWIG target names
# CMP0086: SWIG -module flag support
set(project_policies CMP0078 CMP0086)
foreach(policy ${project_policies})
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

# === Installation directories and output layout ===
include(spsInstallDirs)

# === Interface target for build flags ===
add_library(build INTERFACE)
add_library(Sofus::build ALIAS build)

# === Compiler features ===
target_compile_features(build INTERFACE cxx_std_17)

# === Compiler platform flags ===
include(spsCompilerPlatformFlags)

# === Compiler warning flags ===
include(spsCompilerWarningFlags)

# === Options ===
include(CMakeDependentOption)

option(USE_ProgressBar "Use progress bar" ON)
set(USE_PROGRESS_BAR $<BOOL:${USE_ProgressBar}>)
if(USE_ProgressBar)
  set(USE_PROGRESS_BAR 1)
else()
  set(USE_PROGRESS_BAR 0)
endif()

option(BUILD_FNM_TEST "Build tests" OFF)
option(USE_MQUEUE "Use message queues" OFF)
option(ENABLE_PYTHON3 "Python 3" ON)

set(_DEBUG_USING_PYTHON_RELEASE_RUNTIME 0)
if(MSVC)
  option(SOFUS_DEBUG_USING_PYTHON_RELEASE_RUNTIME "Debug using Python release runtime" ON)
  if(SOFUS_DEBUG_USING_PYTHON_RELEASE_RUNTIME)
    set(_DEBUG_USING_PYTHON_RELEASE_RUNTIME 1)
  endif()
endif()

# Pulsed wave support
set(FNM_PULSED_WAVE 0)
option(SOFUS_Pulsed_waves "Support pulsed waves" OFF)
if(SOFUS_Pulsed_waves)
  set(FNM_PULSED_WAVE 1)
endif()

# === Documentation ===
find_package(Doxygen 1.9 QUIET)
option(BUILD_DOCUMENTATION
  "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  find_program(SED_CMD NAMES sed)
  if(NOT EXISTS "${SED_CMD}")
    message(FATAL_ERROR "sed command is required for building documentation")
  endif()
  include(DoxygenTarget)
  if(DOXYGEN_FOUND)
    add_doxygen_target(doc)
  else()
    message(FATAL_ERROR "Doxygen is required for documentation.")
  endif()
endif()

# === SIMD Detection (from local cmake folder) ===
include(simd)

# === Sub-directories ===
add_subdirectory(sps)
add_subdirectory(gl)
add_subdirectory(sofus)
add_subdirectory(fnm)
add_subdirectory(python)

if(BUILD_FNM_TEST)
  enable_testing()
  add_subdirectory(test)
endif()

# === Export targets (both Fnm and Sofus export sets) ===
# Include the build interface target in exports
export(TARGETS fnm sofus gl sps build
  FILE "${PROJECT_BINARY_DIR}/FnmTargets.cmake")

export(TARGETS fnm sofus gl sps build
  FILE "${PROJECT_BINARY_DIR}/SofusTargets.cmake")

export(PACKAGE Fnm)
export(PACKAGE Sofus)

# === Package configuration ===
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# Installation directory variables (compatible with old cmake/SetInstallDir.cmake)
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for headers")

# --- Fnm package config ---
set(INSTALL_CMAKE_DIR_FNM "${CMAKE_INSTALL_LIBDIR}/cmake/Fnm")
file(RELATIVE_PATH REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR_FNM}"
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
file(RELATIVE_PATH REL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR_FNM}"
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")

# For build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "")
configure_file(FnmConfig.cmake.in
  "${PROJECT_BINARY_DIR}/FnmConfig.cmake" @ONLY)

# For install tree
set(CONF_INCLUDE_DIRS "\${FNM_CMAKE_DIR}/${REL_INCLUDE_DIR}")
set(CONF_LIBRARY_DIRS "\${FNM_CMAKE_DIR}/${REL_BIN_DIR}")
configure_file(FnmConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FnmConfig.cmake" @ONLY)

configure_file(FnmConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/FnmConfigVersion.cmake" @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FnmConfig.cmake"
  "${PROJECT_BINARY_DIR}/FnmConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR_FNM}" COMPONENT dev)

# Install the build interface target to both export sets
install(TARGETS build EXPORT FnmTargets)
install(TARGETS build EXPORT SofusTargets)

install(EXPORT FnmTargets
  DESTINATION "${INSTALL_CMAKE_DIR_FNM}" COMPONENT dev)

# --- Sofus package config ---
set(INSTALL_CMAKE_DIR_SOFUS "${CMAKE_INSTALL_LIBDIR}/cmake/Sofus")

# Version config file for Sofus
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/SofusConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(EXPORT SofusTargets
  DESTINATION "${INSTALL_CMAKE_DIR_SOFUS}" COMPONENT dev)

# === Uninstall target ===
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(UNINSTALL
    COMMAND ${CMAKE_COMMAND} -P
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)
endif()
