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

# On POSIX compliant systems, we would like to support asynchronous execution
check_include_file(pthread.h HAVE_PTHREAD_H)
check_include_file(signal.h HAVE_SIGNAL_H)
check_include_file(mqueue.h HAVE_MQUEUE_H)
check_function_exists(round HAVE_ROUND)

# TODO: Fix why we exit before threads are finished
if (NOT USE_MQUEUE)
  set(HAVE_MQUEUE_H 0)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DSPS_DEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

if (MSVC OR CYGWIN)
  add_definitions(-DHAVE_CONFIG_H)
endif()

option(BUILD_SWIG_INTERFACE "Build SWIG interface" ON)

set(N_MAX_THREADS 8)
set(ACCURATE_TRIGONOMETRICS 0)
option (USE_Accurate_trigonometrics "Use scalar versions" OFF)
if (USE_Accurate_trigonometrics)
  set(ACCURATE_TRIGONOMETRICS 1)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

set(fnm_HEADERS
  ../sps/progress.hpp
  fnm_types.hpp
  FnmSIMD.hpp
  fnm_calc.hpp
  fnm_data.hpp
)

set(fnm_PUBLIC_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/fnm_export.h"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  fnm.hpp
)

set(fnm_SOURCES
  ../sps/progress.cpp
  fnm_types.cpp
  fnm_data.cpp
  fnm_calc.cpp
  fnm.cpp
)

if (UNIX AND NOT WIN32)
  set(fnm_HEADERS "${fnm_HEADERS}" "../sps/strace.hpp")
  set(fnm_SOURCES "${fnm_SOURCES}" "../sps/strace.cpp")
endif()

generate_product_version(
   VersionFilesOutputVariable
   NAME "${PROJECT_NAME}"
   ICON
   VERSION_MAJOR ${FNM_MAJOR_VERSION}
   VERSION_MINOR ${FNM_MINOR_VERSION}
   VERSION_PATCH ${FNM_PATHC_VERSION} #${BUILD_COUNTER}
   VERSION_REVISION ${FNM_TWEAK_VERSION} # Dev uses data format 20110222
)

add_library(fnm ${LIB_TYPE} ${fnm_SOURCES} ${fnm_HEADERS} ${fnm_PUBLIC_HEADERS} config.h.in ${VersionFilesOutputVariable})

if (WIN32)
  if (BUILD_SHARED_LIBS)
    # Works with CMake 2.4
    target_link_libraries(fnm sps gl)
  else()
    link_libraries(fnm sps gl)
  endif()
else()
  link_libraries(fnm sps gl rt)
endif()

generate_export_header(fnm)# DEFINE_NO_DEPRECATED)

set_target_properties(fnm PROPERTIES
  PUBLIC_HEADER "${fnm_PUBLIC_HEADERS};${CMAKE_CURRENT_BINARY_DIR}/config.h;${CMAKE_CURRENT_BINARY_DIR}/fnm_export.h")

if (BUILD_DOCUMENTATION)
  add_dependencies(doc fnm)
endif()

# SWIG Target

if (BUILD_SWIG_INTERFACE)
  if (MSVC)
    add_definitions(-DSWIG_PYTHON_INTERPRETER_NO_DEBUG)
  endif()

  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  
  find_package(PythonLibs REQUIRED)
  
  include_directories(${PYTHON_INCLUDE_PATH})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}) # config.h
  
  if (BUILD_DOCUMENTATION)
    set(DOXYSWIG_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/documentation.i")
    
    # Generate SWIG documentation from Doxygen XML
    add_custom_command(OUTPUT ${DOXYSWIG_OUTPUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND echo 
      COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../bin/doxy2swig.py ${PROJECT_SOURCE_DIR}/${CMAKE_DOXYGEN_OUT}/xml/index.xml ${CMAKE_CURRENT_BINARY_DIR}/documentation.i
      DEPENDS fnm doc
    )
    add_definitions(-DSWIG_INCLUDE_DOCUMENTATION)
  endif()

  if (MSVC)
    set(CMAKE_SWIG_FLAGS "-D_SWIG_WIN32")
  endif()
  
  set_source_files_properties(swig_fnm.i PROPERTIES CPLUSPLUS ON)
  
  set(swig_fnm_HEADERS
    fnm.hpp
    # Move elsewhere
#    block_match.hpp
  )
  
  set_source_files_properties(swig_fnm.i PROPERTIES OBJECT_DEPENDS "${DOXYSWIG_OUTPUT}")

  # HACK to make swig_fnm.i depend on input headers
  execute_process(COMMAND swig -M -python -c++ -I${CMAKE_CURRENT_BINARY_DIR}/.. -I${CMAKE_CURRENT_SOURCE_DIR}/.. swig_fnm.i
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE swig_deps
    INPUT_FILE swig_fnm.i)
  # Match all lines except the first one until " \"
  string(REGEX MATCHALL "\n  [^ ]+" temp ${swig_deps})
  set(swig_deps)
  foreach(t ${temp})
    string(STRIP "${t}" t)
    set(swig_deps ${swig_deps} "${t}")
  endforeach()
  set_property(DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${swig_deps} swig_fnm.i)
  
  swig_add_module(swig_fnm python swig_fnm.i ${swig_fnm_HEADERS})

  if (MSVC)
    # Potential uninitialized variable in SWIG_AsVal_
    set_source_files_properties( ${swig_generated_file_fullname} PROPERTIES COMPILE_FLAGS "/wd4701 /wd4127")

    # Allow to debug under windows, if debug versions of Python are missing
    string(REPLACE "_d" "" PYTHON_LIBRARIES "${PYTHON_LIBRARIES}")
  endif()
  
  swig_link_libraries(swig_fnm fnm ${PYTHON_LIBRARIES})

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
    set_target_properties(_swig_fnm PROPERTIES
      COMPILE_DEFINITIONS "${DirDefs};HAVE_ROUND")
  endif()
endif()

# Enable the test
find_package(PythonInterp)
         
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
      add_custom_command(TARGET _swig_fnm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/${SwigFnmOutputFile}"
        "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION:fnm>")
    endforeach()
  endif()
endif()

install(TARGETS fnm _swig_fnm
  EXPORT FnmTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/fnm" COMPONENT dev)
