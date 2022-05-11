include(CheckIncludeFile)
include(CheckFunctionExists)
include(GenerateExportHeader)
include(GenerateProductVersion)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

if (Threads_FOUND)
  set(HAVE_THREAD 1)
endif()

# DEBUG ISSUES
#ifndef _DEBUG
#define _SECURE_SCL 0
#endif
#add_definitions( -D_SECURE_SCL=0 )
#add_definitions(-DMS_NO_COREDLL)

# On POSIX compliant systems, we would like to support asynchronous execution
check_include_file(pthread.h HAVE_PTHREAD_H)
check_include_file(signal.h HAVE_SIGNAL_H)
check_include_file(mqueue.h HAVE_MQUEUE_H)
check_function_exists(round HAVE_ROUND)

# TODO: Fix why we exit before threads are finished
if (NOT USE_MQUEUE)
  set(HAVE_MQUEUE_H 0)
endif()

#add_definitions(-DSPS_DEBUG)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DSPS_DEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

if (MSVC OR CYGWIN)
  add_definitions(-DHAVE_CONFIG_H)
endif()

option(BUILD_SWIG_INTERFACE "Build SWIG interface" ON)
option(BUILD_SWIG_DOCUMENTATION "Build SWIG documentation" ON)

if (NOT BUILD_DOCUMENTATION)
  set(BUILD_SWIG_DOCUMENTATION OFF)
endif()

set(FNM_CLOSURE_FUNCTIONS 0)
option (FNM_ENABLE_Closures "FNM Enable closure functions" OFF)
if (FNM_ENABLE_Closures)
  set(FNM_CLOSURE_FUNCTIONS 1)
endif()

set(N_MAX_THREADS 8)
set(ACCURATE_TRIGONOMETRICS 0)
option (USE_Accurate_trigonometrics "Use scalar versions" OFF)
if (USE_Accurate_trigonometrics)
  set(ACCURATE_TRIGONOMETRICS 1)
endif()

set (FNM_ENABLE_ATTENUATION 0)
option (FNM_ENABLE_Attenuation "FNM Enable attenuation" ON)
if (FNM_ENABLE_Attenuation)
  set(FNM_ENABLE_ATTENUATION 1)
else()
  set(FNM_ENABLE_ATTENUATION 0)
endif()

set(FNM_DOUBLE_SUPPORT 0)
option(FNM_Double_support "Support double precision" OFF)
if (FNM_Double_support)
  set(FNM_DOUBLE_SUPPORT 1)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

# TODO: Use this everywhere
if (UNIX)
  set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic -Wall ${CMAKE_CXX_FLAGS}")  
  if (ENABLE_SSE2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
  endif()
  if (ENABLE_SSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
  endif()
  if (ENABLE_SSE4)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -msse4.1")
  endif()
  if (ENABLE_AVX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
  endif()
  if (ENABLE_FMA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma")
  endif()
  if (ENABLE_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
  endif()
endif()

# Used by ASA
find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
include_directories(${FFTW_INCLUDES})


set(fnm_HEADERS
  ../sps/progress.hpp
  fnm_types.hpp
  FnmSIMD.hpp
  fnm_calc.hpp
  fnm_data.hpp
  fnm_arrays.hpp
  fnm_basis.hpp
  fnm_common.hpp
  fnm_delays.hpp
  circular_data.hpp
  circular_calc.hpp
  circular.hpp
)

set(fnm_PUBLIC_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/fnm_export.h"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  fnm.hpp
  if_fnm.h
)

set(fnm_SOURCES
  ../sps/progress.cpp
  fnm_types.cpp
  fnm_data.cpp
  fnm_basis.cpp
  fnm_calc.cpp
  fnm_delays.cpp
  fnm.cpp
  fnm_arrays.cpp
  fnm_common.cpp
  fnm_profiling.cpp
  if_fnm.cpp
  circular.cpp
  circular_data.cpp
  circular_calc.cpp
)

set(fnm_HEADERS "${fnm_HEADERS}" "fnm_transient.hpp")
set(fnm_SOURCES "${fnm_SOURCES}" "fnm_transient.cpp")
set(fnm_HEADERS "${fnm_HEADERS}" "fnm_transient_threaded.hpp")
set(fnm_SOURCES "${fnm_SOURCES}" "fnm_transient_threaded.cpp")
set(fnm_HEADERS "${fnm_HEADERS}" "fnm_response.hpp")
set(fnm_SOURCES "${fnm_SOURCES}" "fnm_response.cpp")

if (UNIX AND NOT WIN32)
  set(fnm_HEADERS "${fnm_HEADERS}" "../sps/strace.hpp")
  set(fnm_SOURCES "${fnm_SOURCES}" "../sps/strace.cpp")
endif()

if (0)
generate_product_version(
   VersionFilesOutputVariable
   NAME "${PROJECT_NAME}"
   ICON
   VERSION_MAJOR ${FNM_MAJOR_VERSION}
   VERSION_MINOR ${FNM_MINOR_VERSION}
   VERSION_PATCH ${FNM_PATHC_VERSION} #${BUILD_COUNTER}
   VERSION_REVISION ${FNM_TWEAK_VERSION} # Dev uses data format 20110222
)
endif()

add_library(fnm ${LIB_TYPE} ${fnm_SOURCES} ${fnm_HEADERS} ${fnm_PUBLIC_HEADERS} config.h.in ${VersionFilesOutputVariable})

if (WIN32)
  if (BUILD_SHARED_LIBS)
    # Works with CMake 2.4
    target_link_libraries(fnm sps gl ${FFTW_LIBRARIES})
    if (FNM_PULSED_WAVE)
      target_link_libraries(fnm sofus)
    endif()
  else()
    link_libraries(fnm sps gl ${FFTW_LIBRARIES})
    if (FNM_PULSED_WAVE)
      target_link_libraries(fnm sofus)
    endif()
  endif()
else()
  link_libraries(fnm sps gl rt dl ${FFTW_LIBRARIES} pthread)
  if (FNM_PULSED_WAVE)
    target_link_libraries(fnm sofus)
  endif()
endif()

target_link_libraries(fnm sofus)

generate_export_header(fnm)# DEFINE_NO_DEPRECATED)

set_target_properties(fnm PROPERTIES
  PUBLIC_HEADER "${fnm_PUBLIC_HEADERS};${CMAKE_CURRENT_BINARY_DIR}/config.h;${CMAKE_CURRENT_BINARY_DIR}/fnm_export.h")

if (BUILD_DOCUMENTATION)
  add_dependencies(doc fnm)
endif()

# SWIG Target

if (BUILD_SWIG_INTERFACE)
  if (MSVC)
    if(_DEBUG_USING_PYTHON_RELEASE_RUNTIME)
      add_definitions(-DSWIG_PYTHON_INTERPRETER_NO_DEBUG)
    endif()
  endif()

  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  if (WIN32)
    find_package(Python3 COMPONENTS Interpreter Development NumPy)
  else()
    find_package(Python3 3.9 COMPONENTS Interpreter Development)
  endif()
  include_directories(${Python3_INCLUDE_DIRS})
  include_directories(${Python3_NumPy_INCLUDE_DIRS})

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}) # config.h

  if (BUILD_DOCUMENTATION)
    set(DOXYSWIG_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/documentation.i")

    # Generate SWIG documentation from Doxygen XML
    add_custom_target(swig_fnm_doxy
      COMMAND echo
      COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../bin/doxy2swig.py -c -a ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/xml/index.xml ${DOXYSWIG_OUTPUT}
      DEPENDS fnm doc
    )
  endif()

  if (MSVC)
    if (BUILD_SWIG_DOCUMENTATION)
      set_property(SOURCE swig_fnm.i PROPERTY SWIG_FLAGS "-D_SWIG_WIN32" "-DSWIG_INCLUDE_DOCUMENTATION")
    endif()
  endif()
  
  set_source_files_properties(swig_fnm.i PROPERTIES CPLUSPLUS ON)
  set_property(SOURCE swig_fnm.i PROPERTY SWIG_MODULE_NAME swig_fnm)
  set_property(SOURCE swig_fnm.i PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/..")

  set(swig_fnm_HEADERS
    fnm.hpp
  )

  if (BUILD_SWIG_DOCUMENTATION)
    set_source_files_properties(swig_fnm.i PROPERTIES OBJECT_DEPENDS "${DOXYSWIG_OUTPUT}")
  endif()
  execute_process(COMMAND swig -python -c++ -external-runtime ${CMAKE_CURRENT_BINARY_DIR}/swig_runtime.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  if (BUILD_SWIG_DOCUMENTATION AND NOT WIN32)
    set_property(SOURCE swig_fnm.i PROPERTY SWIG_FLAGS "-DSWIG_INCLUDE_DOCUMENTATION")
  endif()

  # By adding USE_SWIG_DEPENDENCIES TRUE in
  # set_source_files_properties, the swig wrapper is rebuilt each time
  # a C++ header is modified!
  # Available from CMake 3.20
  if(NOT ${CMAKE_VERSION} VERSION_LESS 3.20)
    set_property(SOURCE swig_fnm.i PROPERTY USE_SWIG_DEPENDENCIES TRUE)
  endif()

  swig_add_library(swig_fnm LANGUAGE python SOURCES swig_fnm.i ${swig_fnm_HEADERS} ${swig_deps_actual})

  # Removes generated file
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp
    COMMAND ${CMAKE_COMMAND} -E remove ${swig_generated_file_fullname}
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp
    DEPENDS ${swig_deps_actual}
    COMMENT "Removing old SWIG generated file" VERBATIM)

  # Custom target for establishing dependency
  add_custom_target(
    swigtrick
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/swig.stamp
    )
  if (BUILD_SWIG_DOCUMENTATION)
    # Depedency
    add_dependencies(swig_fnm swig_fnm_doxy)
  endif()

  if (MSVC)
    # Potential uninitialized variable in SWIG_AsVal_
    set_source_files_properties( ${swig_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/wd4701 /wd4127")

    # Allow to debug under windows, if debug versions of Python are missing
    if(_DEBUG_USING_PYTHON_RELEASE_RUNTIME)
      string(REPLACE "_d" "" PYTHON_LIBRARIES "${PYTHON_LIBRARIES}")
    endif()
  endif()
  
  swig_link_libraries(swig_fnm fnm ${Python3_LIBRARIES})

  add_custom_command(TARGET swig_fnm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/dicts.py"
    "$<TARGET_FILE_DIR:swig_fnm>")

  add_custom_command(TARGET swig_fnm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/pyutils.py"
    "$<TARGET_FILE_DIR:swig_fnm>")

  add_custom_command(TARGET swig_fnm POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/euler.py"
    "$<TARGET_FILE_DIR:swig_fnm>")

  # Add python files
  set(swig_python_files
    # Reference code
    dicts.py
    euler.py
    fnm.py
    # SWIG'ed code
    addpaths.py
    ${swig_extra_generated_files}
    fnm_arrays.py
    fnm_test.py
    # Examples (TODO: move to examples folder)
    test_cw_pressure.py
    test_element.py
    )
  
  install(FILES
    ${swig_python_files}
    DESTINATION python)

  if (WIN32)
    # pyconfig.h is not autogenerated on Windows. To avoid warnings, we
    # add a compiler directive
    get_directory_property(DirDefs COMPILE_DEFINITIONS )
    set_target_properties(swig_fnm PROPERTIES
      COMPILE_DEFINITIONS "${DirDefs};HAVE_ROUND")
  endif()
endif()

# Copy dependencies on WIN32
if (WIN32)
  if (BUILD_SHARED_LIBS)
    add_custom_command(TARGET fnm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "$<TARGET_FILE:gl>"
      "$<TARGET_FILE_DIR:fnm>")
    add_custom_command(TARGET fnm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "$<TARGET_FILE:sps>"
      "$<TARGET_FILE_DIR:fnm>")
  endif()

  if (BUILD_SWIG_INTERFACE)
    set(SwigFnmOutputFiles
      swig_fnm.py
      )
    foreach(SwigFnmOutputFile ${SwigFnmOutputFiles})
      add_custom_command(TARGET swig_fnm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/${SwigFnmOutputFile}"
        "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION:fnm>")
    endforeach()
  endif()

  # Copy FFTW dlls into target directory
  file(GLOB FFTW_DLLs "${FFTW_INCLUDES}/*.dll")

  # List of DLLs
  set(dll_copied_modules)

  foreach (input_dll_file IN LISTS FFTW_DLLs)
    # CMake mixes \ and /'s
    string(REPLACE "\\" "/" input_dll_file "${input_dll_file}")

    # Output file
    get_filename_component(dll_file "${input_dll_file}" NAME)

    set(output_dll_file
      "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION:fnm>/${dll_file}")

    add_custom_command(
      OUTPUT  "${output_dll_file}"
      DEPENDS "${input_dll_file}"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
      "${input_dll_file}"
      "${output_dll_file}"
      COMMENT "Copying ${dll_file} to the binary directory")

    list(APPEND dll_copied_modules
      "${output_dll_file}")
  endforeach()

  add_custom_target(sofus_fftw_copy ALL
    DEPENDS
    ${dll_copied_modules})
endif()

dump_cmake_variables("Python")
dump_cmake_variables("PYTHON")

install(TARGETS fnm swig_fnm
  EXPORT SofusTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/fnm" COMPONENT dev)
