include(GenerateExportHeader)

# === Options ===
option(BUILD_GL_TEST "Build tests" OFF)
option(GL_SHARED_LIBS "Build Shared Libraries" OFF)

set(GL_LIB_TYPE STATIC)
if(GL_SHARED_LIBS)
  set(GL_LIB_TYPE SHARED)
endif()

# === Dependencies ===
find_package(Threads REQUIRED)

# Mathematica is optional - used for generating LUT tables
include(spsMathematica)

# === LUT generation ===
set(_BUILD_LUT OFF)
set(_GL_LUT_TABLE_SIZE 101 CACHE STRING "Max number of abcissa used for Gauss-Legendre integration")
set(_BESSELJZERO_LUT_TABLE_SIZE 20 CACHE STRING "First zeros of BesselJ(0,z)")
set(_BESSELJ_1_SQUARED_LUT_TABLE_SIZE 21 CACHE STRING "Square of BesselJ(1,BesselZero(0,k))")

if(MATHEMATICA_FOUND)
  option(BUILD_LUT "Generate LUT using Mathematica" OFF)
  set(_BUILD_LUT "${BUILD_LUT}")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  @ONLY)

# === Source files ===
set(gl_HEADERS
  gl_lut.hpp
  bessel_lut.hpp
)

set(gl_PUBLIC_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/gl_export.h"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  gl.hpp
)

set(gl_SOURCES
  gl.cpp
)

if(_BUILD_LUT)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp
    COMMENT "Generating Gauss-Legendre LUT" VERBATIM
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.m ${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp ${_GL_LUT_TABLE_SIZE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.m)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp
    COMMENT "Generating Bessel LUT" VERBATIM
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.m ${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp ${_BESSELJZERO_LUT_TABLE_SIZE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.m)

  list(APPEND gl_SOURCES
    "${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp")
else()
  list(APPEND gl_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.cpp")
endif()

# === Library target ===
add_library(gl ${GL_LIB_TYPE} ${gl_SOURCES} ${gl_HEADERS} ${gl_PUBLIC_HEADERS})

target_include_directories(gl PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

# Link to build interface for compiler flags
if(TARGET build)
  target_link_libraries(gl PRIVATE $<BUILD_INTERFACE:build>)
endif()

generate_export_header(gl)

set_target_properties(gl PROPERTIES
  PUBLIC_HEADER "${gl_PUBLIC_HEADERS};${CMAKE_CURRENT_BINARY_DIR}/config.h;${CMAKE_CURRENT_BINARY_DIR}/gl_export.h")

# === Tests ===
if(BUILD_GL_TEST)
  find_package(GTest REQUIRED)

  add_executable(gl_test gl_test.cpp fastgl.cpp)

  target_link_libraries(gl_test PRIVATE
    gl
    GTest::GTest
    GTest::Main
    Threads::Threads)

  if(TARGET build)
    target_link_libraries(gl_test PRIVATE build)
  endif()

  enable_testing()
  add_test(NAME gl_test COMMAND gl_test)
endif()

# === Installation ===
# Install to both FnmTargets and SofusTargets export sets
install(TARGETS gl
  EXPORT FnmTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gl" COMPONENT dev)

install(TARGETS gl
  EXPORT SofusTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gl" COMPONENT dev)
