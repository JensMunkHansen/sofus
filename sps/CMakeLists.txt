include(cmake/utils.cmake)
include(cmake/simd.cmake)
include(GenerateExportHeader)

# TODO: Add target SPS

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

include(CheckIncludeFile)
include(CheckSymbolExists)

check_symbol_exists(itoa stdlib.h HAVE_ITOA)
check_symbol_exists(utoa stdlib.h HAVE_UTOA)
check_symbol_exists(ptoa stdlib.h HAVE_PTOA)

#generate_export_header(sps)

if (NOT TARGET signals_test)
  check_include_file(stdint.h HAVE_STDINT_H)
  
  include_directories(
    "${PROJECT_BINARY_DIR}")  # to find sps/config.h
  
  find_package(FFTW)
  if (WIN32)
    # TODO: Support non-default install locations on *NIX
    include_directories(${FFTW_INCLUDES})
  endif()

  # Rotation convention:
  set (RotationConvention "yxy" CACHE STRING "Angle convention used for angular arguments")
  set(RotationConventionValues "zyz;yxy;xyz")
  set_property(CACHE RotationConvention PROPERTY STRINGS ${RotationConventionValues})
  
  if (${RotationConvention} STREQUAL "xyz")
    # Extrinsic (i.e. non-rotating)
    set(ROTATION_CONVENTION "2")
  elseif("${RotationConvention}" STREQUAL "yxy")
    # Intrinsic (i.e. rotating)
    set(ROTATION_CONVENTION "1")
  else() # Convention used is zyz
    # Intrinsic (i.e. rotating)
    set(ROTATION_CONVENTION "0")
  endif()
  
  set (USE_FFTW_THREADS 0)
  option (USE_FFTWThreads "Use threading for FFTW" OFF)
  if (USE_FFTWThreads)
    set(USE_FFTW_THREADS 1)
  endif()
  
  if (MSVC OR CYGWIN)
    add_definitions(-DHAVE_CONFIG_H)
  endif()
  
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)
  
  # Supported from CMake 3.1
  # target_compile_features(signals_test PRIVATE cxx_constexpr)
  use_cxx11()
  
  if (CYGWIN)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -msse2 -msse3 -msse4 -mavx -pedantic -Wall")
  endif()

  set(sps_HEADERS
    smath.hpp
    signals.hpp
    sps_threads.hpp
  )
  
  set(sps_SOURCES
    smath.cpp
    signals.cpp
  )

  add_library(sps ${sps_HEADERS} ${sps_SOURCES})

  generate_export_header(sps)

if (BUILD_TESTS)  
  include_directories(${GTEST_INCLUDE_DIR})
  
  add_executable(signals_test signals_test.cpp signals.cpp)
  
  # Add dependencies
  add_gtest_dependencies(signals_test)
  target_link_libraries(signals_test ${FFTW_LIBRARIES})

  add_executable(trigintrin_test trigintrin_test.cpp)
  
  # Add dependencies
  add_gtest_dependencies(trigintrin_test)

  add_executable(extintrin_test extintrin_test.cpp)
  
  # Add dependencies
  add_gtest_dependencies(extintrin_test)
  
  # Enable the tests
  enable_testing()
  add_test(NAME    signals_test 
           COMMAND signals_test)

  add_test(NAME    trigintrin_test 
           COMMAND trigintrin_test)

  add_test(NAME    extintrin_test 
           COMMAND extintrin_test)
endif()
	 
endif()
