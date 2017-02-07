cmake_minimum_required(VERSION 3.0)

include(cmake/utils.cmake)
include(cmake/simd.cmake)
include(cmake/gtest_macros.cmake)
include(GenerateExportHeader)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  )

include(CheckIncludeFile)
include(CheckSymbolExists)

check_symbol_exists(itoa stdlib.h HAVE_ITOA)
check_symbol_exists(utoa stdlib.h HAVE_UTOA)
check_symbol_exists(ptoa stdlib.h HAVE_PTOA)

option (SPS_Signals "Include signal processing" OFF)

option (BUILD_SPS_TEST "Build SPS tests" OFF)

check_include_file(stdint.h HAVE_STDINT_H)
check_include_file(mqueue.h HAVE_MQUEUE_H)

include_directories(
  "${PROJECT_BINARY_DIR}")  # to find sps/config.h
include_directories(
  "${PROJECT_SOURCE_DIR}")  # to find sps/xxx.h

find_package(FFTW REQUIRED)
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

set(HAVE_CONFIG_H 1)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

use_cxx11()

if (CYGWIN)
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -msse2 -msse3 -msse4 -mavx -pedantic -Wall")
endif()

#add_definitions(-DSPS_DEBUG)

set(sps_HEADERS
  smath.hpp
  smath_types.hpp
  sps_threads.hpp
  multi_malloc.hpp
  memory
)

set(sps_SOURCES
  smath.cpp
  # smath_types.cpp - clean up (it is obsolete as is)
)

if (USE_PROGRESS_BAR)
  # Must be defined for stand-alone sps library
  set(sps_HEADERS "${sps_HEADERS}" "progress.hpp" "progress_if.hpp")
  set(sps_SOURCES "${sps_SOURCES}" "progress.cpp")
endif()

if (SPS_Signals)
  # Must be defined for stand-alone sps library
  set(sps_HEADERS "${sps_HEADERS}" "signals.hpp")
  set(sps_SOURCES "${sps_SOURCES}" "signals.cpp")
endif()

if (HAVE_MQUEUE_H)
  set(sps_HEADERS "${sps_HEADERS}" "sps_mqueue.hpp")
  set(sps_SOURCES "${sps_SOURCES}" "sps_mqueue.cpp")
endif()

if (UNIX AND NOT WIN32)
  set(sps_SOURCES "${sps_SOURCES}" "strace.cpp")
endif()

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
endif()

add_library(sps ${LIB_TYPE} ${sps_HEADERS} ${sps_SOURCES})

if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set_property(TARGET sps PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

if (SPS_Signals)
  if (WIN32)
    link_directories(${FFTW_INCLUDES})
    # Is this necessary???
    if (BUILD_SHARED_LIBS)
      target_link_libraries(sps ${FFTW_LIBRARIES})
    else()
      link_libraries(sps ${FFTW_LIBRARIES})
    endif()
  else()
    link_libraries(sps ${FFTW_LIBRARIES})
  endif()
endif()

generate_export_header(sps)# DEFINE_NO_DEPRECATED)

if (BUILD_SPS_TEST)
  find_package( GoogleTest REQUIRED )
  find_package( Threads REQUIRED )
  
  include_directories(${GTEST_INCLUDE_DIRS})

  if (UNIX)
    add_executable(indexed_types_test indexed_types_test.cpp)
    add_gtest_dependencies(indexed_types_test)
  endif()
  add_executable(extintrin_test extintrin_test.cpp)
  add_gtest_dependencies(extintrin_test)

  add_executable(trigintrin_test trigintrin_test.cpp)
  add_gtest_dependencies(trigintrin_test)

  add_executable(smath_test smath_test.cpp smath.cpp)
  add_gtest_dependencies(smath_test)

  add_executable(memory_test memory_test.cpp)
  add_gtest_dependencies(memory_test)

  if (SPS_Signals)
    add_executable(signals_test signals_test.cpp signals.cpp)
    add_gtest_dependencies(signals_test)
  endif()
  
  # Enable the tests
  enable_testing()
  add_test(NAME    indexed_types_test
           COMMAND indexed_types_test)
    
  add_test(NAME    extintrin_test 
           COMMAND extintrin_test)
  
  add_test(NAME    trigintrin_test 
           COMMAND trigintrin_test)

  add_test(NAME    smath_test 
           COMMAND smath_test)
	 
  add_test(NAME    memory_test 
           COMMAND memory_test)

  if (SPS_Signals)
    add_test(NAME    signals_test 
             COMMAND signals_test)
  endif()
  
  # Link to threading library (general)
  set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
endif()

# SWIG Target

if (SPS_Signals)
  # Need signals.cpp and link to FFTW
  
  # Figure out why cstdarg => FILE__ undefined
  if (MSVC)
    add_definitions(-DSWIG_PYTHON_INTERPRETER_NO_DEBUG)
  endif()

  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  find_package(PythonLibs REQUIRED)

  include_directories(${PYTHON_INCLUDE_PATH})

  if (MSVC)
    set(CMAKE_SWIG_FLAGS "-D_SWIG_WIN32")
  endif()

  set_source_files_properties(swig_signals.i PROPERTIES CPLUSPLUS ON)

  set(swig_signals_HEADERS
    signals.hpp
    )

  set(swig_signals_SOURCES
    )

  swig_add_module(swig_signals python swig_signals.i ${swig_signals_HEADERS} ${swig_signals_SOURCES})

  if (MSVC)
    # Potential uninitialized variable in SWIG_AsVal_
    set_source_files_properties( ${swig_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/wd4701 /wd4127")

    # Allow to debug under windows, if debug versions of Python are missing
    string(REPLACE "_d" "" PYTHON_LIBRARIES "${PYTHON_LIBRARIES}")
  endif()

  swig_link_libraries(swig_signals sps ${FFTW_LIBRARIES} ${PYTHON_LIBRARIES})

  if (WIN32)
    # pyconfig.h is not autogenerated on Windows. To avoid warnings, we
    # add a compiler directive
    get_directory_property(DirDefs COMPILE_DEFINITIONS )
    set_target_properties(_swig_signals PROPERTIES
      COMPILE_DEFINITIONS "${DirDefs};HAVE_ROUND")
  endif()

endif()

# TODO: Create an "export-set" for SPS and import this in FNM and SOFUS using
#       find_package or include(sps.cmake)

if( ${PROJECT_NAME} STREQUAL "Fnm")
  install(TARGETS sps
    # IMPORTANT: Add the fnm library to the "export-set"
    EXPORT FnmTargets
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/sps"
    COMPONENT dev)
endif()

if( ${PROJECT_NAME} STREQUAL "Sofus")
  install(TARGETS sps
    # IMPORTANT: Add the fnm library to the "export-set"
    EXPORT SofusTargets
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/sps"
    COMPONENT dev)
endif()


