include(GenerateExportHeader)

option (BUILD_GL_TEST "Build tests" OFF)

option (GL_SHARED_LIBS "Build Shared Libraries" OFF)
set(GL_LIB_TYPE STATIC)
if (GL_SHARED_LIBS)
  set(GL_LIB_TYPE SHARED)
endif()

if (BUILD_GL_TEST)
  find_package( GTest REQUIRED )
endif()

find_package( Threads REQUIRED )
find_package( Mathematica )

include_directories(${GTEST_INCLUDE_DIRS})

set(_BUILD_LUT OFF)
set(_GL_LUT_TABLE_SIZE 101 CACHE STRING "Max number of abcissa used for Gauss-Legendre integration")

if (MATHEMATICA_FOUND)
  option(BUILD_LUT "Generate LUT using Mathematica" OFF)
  set(_BUILD_LUT "${BUILD_LUT}")
else()
endif()

set(_BESSELJZERO_LUT_TABLE_SIZE 20 CACHE STRING "First zeros of BesselJ(0,z)")
set(_BESSELJ_1_SQUARED_LUT_TABLE_SIZE 21 CACHE STRING "Square of BesselJ(1,BesselZero(0,k))")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

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

if (_BUILD_LUT)
  add_custom_command(OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp
    COMMENT "Generating Gauss-Legendre LUT" VERBATIM
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.m ${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp ${_GL_LUT_TABLE_SIZE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.m)
  add_custom_command(OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp
    COMMENT "Generating Bessel LUT" VERBATIM
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.m ${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp ${_BESSELJZERO_LUT_TABLE_SIZE}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.m)
endif()

if (_BUILD_LUT)
  set(gl_SOURCES "${gl_SOURCES}" "${CMAKE_CURRENT_BINARY_DIR}/gl_lut.cpp")
  set(gl_SOURCES "${gl_SOURCES}" "${CMAKE_CURRENT_BINARY_DIR}/bessel_lut.cpp")
else()
  set(gl_SOURCES "${gl_SOURCES}" "${CMAKE_CURRENT_SOURCE_DIR}/gl_lut.cpp")
  set(gl_SOURCES "${gl_SOURCES}" "${CMAKE_CURRENT_SOURCE_DIR}/bessel_lut.cpp")
endif()

add_library(gl ${GL_LIB_TYPE} ${gl_SOURCES} ${gl_HEADERS} ${gl_PUBLIC_HEADERS})

generate_export_header(gl)# DEFINE_NO_DEPRECATED)

set_target_properties(gl PROPERTIES
  PUBLIC_HEADER "${gl_PUBLIC_HEADERS};${CMAKE_CURRENT_BINARY_DIR}/config.h;${CMAKE_CURRENT_BINARY_DIR}/gl_export.h")

if (BUILD_DOCUMENTATION)
#  add_doxygen_target(doc_gl)
endif()

if (BUILD_GL_TEST)
  add_executable(gl_test gl_test.cpp fastgl.cpp)

  # Targets are defined in GoogleTest.cmake
  target_link_libraries(gl_test ${GTEST_LIBRARIES})

  target_link_libraries(gl_test gl)

  # Link to threading library
  set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")

  enable_testing()
  add_test(NAME    gl_test 
           COMMAND gl_test)
endif()

# CMake has screwed up - use export() instead
install(TARGETS gl
  EXPORT SofusTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/gl"
  COMPONENT dev)
