include(GenerateExportHeader)
include(CheckIncludeFile)

check_include_file(signal.h HAVE_SIGNAL_H)

# === Options ===
option(SOFUS_SHARED_LIBS "Build Shared Libraries" OFF)

set(SOFUS_LIB_TYPE STATIC)
if(SOFUS_SHARED_LIBS)
  set(SOFUS_LIB_TYPE SHARED)
endif()

# Two-dimensional computation (Faster if ON)
set(USE_PROJECTIONS 0)
option(SOFUS_USE_Projections "Use projection for computing radii in planes" ON)
if(SOFUS_USE_Projections)
  set(USE_PROJECTIONS 1)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  @ONLY)

# === Source files ===
set(sofus_HEADERS
  rect_int_limits.hpp
  sofus_version.hpp
  sofus_pulses.hpp
  sofus_calc.hpp
)

set(sofus_SOURCES
  sofus_version.cpp
  sofus_pulses.cpp
  sofus_calc.cpp
)

# === Library target ===
add_library(sofus ${SOFUS_LIB_TYPE} ${sofus_SOURCES} ${sofus_HEADERS} config.h.in)

target_include_directories(sofus PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

# Link to build interface for compiler flags
if(TARGET build)
  target_link_libraries(sofus PRIVATE $<BUILD_INTERFACE:build>)
endif()

generate_export_header(sofus)

# === Installation ===
# Install to both FnmTargets and SofusTargets export sets
install(TARGETS sofus
  EXPORT FnmTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sofus" COMPONENT dev)

install(TARGETS sofus
  EXPORT SofusTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sofus" COMPONENT dev)
