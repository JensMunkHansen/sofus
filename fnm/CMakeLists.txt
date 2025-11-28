include(CheckIncludeFile)
include(CheckFunctionExists)
include(GenerateExportHeader)

# === Threading ===
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

if(Threads_FOUND)
  set(HAVE_THREAD 1)
endif()

# === Feature detection ===
check_include_file(pthread.h HAVE_PTHREAD_H)
check_include_file(signal.h HAVE_SIGNAL_H)
check_include_file(mqueue.h HAVE_MQUEUE_H)
check_function_exists(round HAVE_ROUND)

if(NOT USE_MQUEUE)
  set(HAVE_MQUEUE_H 0)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DSPS_DEBUG)
endif()

if(MSVC OR CYGWIN)
  add_definitions(-DHAVE_CONFIG_H)
endif()

# === Options ===
option(BUILD_SWIG_INTERFACE "Build SWIG interface" ON)
option(BUILD_SWIG_DOCUMENTATION "Build SWIG documentation" ON)

if(NOT BUILD_DOCUMENTATION)
  set(BUILD_SWIG_DOCUMENTATION OFF)
endif()

set(FNM_CLOSURE_FUNCTIONS 0)
option(FNM_ENABLE_Closures "FNM Enable closure functions" OFF)
if(FNM_ENABLE_Closures)
  set(FNM_CLOSURE_FUNCTIONS 1)
endif()

set(N_MAX_THREADS 16)

set(FNM_ENABLE_ATTENUATION 0)
option(FNM_ENABLE_Attenuation "FNM Enable attenuation" ON)
if(FNM_ENABLE_Attenuation)
  set(FNM_ENABLE_ATTENUATION 1)
endif()

set(FNM_DOUBLE_SUPPORT 0)
option(FNM_Double_support "Support double precision" OFF)
if(FNM_Double_support)
  set(FNM_DOUBLE_SUPPORT 1)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  @ONLY)

# === Dependencies ===
find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)

# === Source files ===
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
  fnm_transient.hpp
  fnm_transient_threaded.hpp
  fnm_response.hpp
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
  fnm_transient.cpp
  fnm_transient_threaded.cpp
  fnm_response.cpp
  if_fnm.cpp
  circular.cpp
  circular_data.cpp
  circular_calc.cpp
)

# === Library type ===
set(FNM_LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(FNM_LIB_TYPE SHARED)
endif()

# === Library target ===
add_library(fnm ${FNM_LIB_TYPE} ${fnm_SOURCES} ${fnm_HEADERS} ${fnm_PUBLIC_HEADERS} config.h.in)

target_include_directories(fnm PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

target_include_directories(fnm PRIVATE ${FFTW_INCLUDES})

# Link to build interface for compiler flags
if(TARGET build)
  target_link_libraries(fnm PRIVATE $<BUILD_INTERFACE:build>)
endif()

# Link dependencies
target_link_libraries(fnm PUBLIC
  sps
  gl
  sofus
  ${FFTW_LIBRARIES}
  Threads::Threads)

if(UNIX)
  target_link_libraries(fnm PUBLIC rt dl)
endif()

generate_export_header(fnm)

set_target_properties(fnm PROPERTIES
  PUBLIC_HEADER "${fnm_PUBLIC_HEADERS};${CMAKE_CURRENT_BINARY_DIR}/config.h;${CMAKE_CURRENT_BINARY_DIR}/fnm_export.h")

if(BUILD_DOCUMENTATION)
  add_dependencies(doc fnm)
endif()

# === SWIG Python bindings ===
if(BUILD_SWIG_INTERFACE)
  if(MSVC AND _DEBUG_USING_PYTHON_RELEASE_RUNTIME)
    add_definitions(-DSWIG_PYTHON_INTERPRETER_NO_DEBUG)
  endif()

  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})

  if(WIN32)
    find_package(Python3 COMPONENTS Interpreter Development NumPy)
  else()
    find_package(Python3 3.9 COMPONENTS Interpreter Development)
  endif()

  if(Python3_FOUND)
    if(BUILD_DOCUMENTATION)
      set(DOXYSWIG_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/documentation.i")

      add_custom_target(swig_fnm_doxy
        COMMAND echo
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../bin/doxy2swig.py -c -a
          ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/xml/index.xml ${DOXYSWIG_OUTPUT}
        DEPENDS fnm doc
      )
    endif()

    if(MSVC AND BUILD_SWIG_DOCUMENTATION)
      set_property(SOURCE swig_fnm.i PROPERTY SWIG_FLAGS "-D_SWIG_WIN32" "-DSWIG_INCLUDE_DOCUMENTATION")
    endif()

    set_source_files_properties(swig_fnm.i PROPERTIES CPLUSPLUS ON)
    set_property(SOURCE swig_fnm.i PROPERTY SWIG_MODULE_NAME swig_fnm)
    set_property(SOURCE swig_fnm.i PROPERTY INCLUDE_DIRECTORIES
      "${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")

    set(swig_fnm_HEADERS fnm.hpp)

    if(BUILD_SWIG_DOCUMENTATION)
      set_source_files_properties(swig_fnm.i PROPERTIES OBJECT_DEPENDS "${DOXYSWIG_OUTPUT}")
    endif()

    execute_process(
      COMMAND swig -python -c++ -external-runtime ${CMAKE_CURRENT_BINARY_DIR}/swig_runtime.h
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    if(BUILD_SWIG_DOCUMENTATION AND NOT WIN32)
      set_property(SOURCE swig_fnm.i PROPERTY SWIG_FLAGS "-DSWIG_INCLUDE_DOCUMENTATION")
    endif()

    if(NOT ${CMAKE_VERSION} VERSION_LESS 3.20)
      set_property(SOURCE swig_fnm.i PROPERTY USE_SWIG_DEPENDENCIES TRUE)
    endif()

    swig_add_library(swig_fnm
      LANGUAGE python
      SOURCES swig_fnm.i ${swig_fnm_HEADERS})

    target_include_directories(swig_fnm PRIVATE
      ${PROJECT_SOURCE_DIR}
      ${PROJECT_BINARY_DIR}
      ${Python3_INCLUDE_DIRS}
      ${Python3_NumPy_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR})

    if(BUILD_SWIG_DOCUMENTATION)
      add_dependencies(swig_fnm swig_fnm_doxy)
    endif()

    if(MSVC)
      set_source_files_properties(${swig_generated_file_fullname}
        PROPERTIES COMPILE_FLAGS "/wd4701 /wd4706 /wd4127 /wd4456")
      if(_DEBUG_USING_PYTHON_RELEASE_RUNTIME)
        string(REPLACE "_d" "" PYTHON_LIBRARIES "${PYTHON_LIBRARIES}")
      endif()
    endif()

    swig_link_libraries(swig_fnm fnm ${Python3_LIBRARIES})

    # Copy Python files to target directory
    set(fnm_python_support_files
      dicts.py
      grid.py
      euler.py
      ../python/pyutils.py)

    foreach(pyfile ${fnm_python_support_files})
      add_custom_command(TARGET swig_fnm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/${pyfile}"
        "$<TARGET_FILE_DIR:swig_fnm>")
    endforeach()

    # Install Python files
    set(swig_python_files
      dicts.py
      euler.py
      fnm.py
      grid.py
      addpaths.py
      ${CMAKE_CURRENT_SOURCE_DIR}/../python/pyutils.py
      ${swig_extra_generated_files}
      fnm_arrays.py
      fnm_test.py
      test_cw_pressure.py
      test_element.py
    )

    install(FILES ${swig_python_files}
      DESTINATION python/fnm)

    if(WIN32)
      get_directory_property(DirDefs COMPILE_DEFINITIONS)
      set_target_properties(swig_fnm PROPERTIES
        COMPILE_DEFINITIONS "${DirDefs};HAVE_ROUND")
    endif()
  endif()
endif()

# === Copy dependencies on WIN32 ===
if(WIN32)
  if(BUILD_SHARED_LIBS)
    add_custom_command(TARGET fnm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:gl>" "$<TARGET_FILE_DIR:fnm>")
    add_custom_command(TARGET fnm POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:sps>" "$<TARGET_FILE_DIR:fnm>")
  endif()

  if(BUILD_SWIG_INTERFACE AND Python3_FOUND)
    set(SwigFnmOutputFiles swig_fnm.py)
    foreach(SwigFnmOutputFile ${SwigFnmOutputFiles})
      add_custom_command(TARGET swig_fnm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/${SwigFnmOutputFile}"
        "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION:fnm>")
    endforeach()

    # Copy FFTW DLLs
    file(GLOB FFTW_DLLs "${FFTW_INCLUDES}/*.dll")
    set(dll_copied_modules)

    foreach(input_dll_file IN LISTS FFTW_DLLs)
      string(REPLACE "\\" "/" input_dll_file "${input_dll_file}")
      get_filename_component(dll_file "${input_dll_file}" NAME)
      set(output_dll_file "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION:fnm>/${dll_file}")

      add_custom_command(
        OUTPUT "${output_dll_file}"
        DEPENDS "${input_dll_file}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${input_dll_file}" "${output_dll_file}"
        COMMENT "Copying ${dll_file} to the binary directory")

      list(APPEND dll_copied_modules "${output_dll_file}")
    endforeach()

    add_custom_target(sofus_fftw_copy ALL DEPENDS ${dll_copied_modules})
  endif()
endif()

# === Installation ===
# Install to both FnmTargets and SofusTargets export sets
install(TARGETS fnm
  EXPORT FnmTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/fnm" COMPONENT dev)

install(TARGETS fnm
  EXPORT SofusTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/fnm" COMPONENT dev)

if(BUILD_SWIG_INTERFACE AND TARGET swig_fnm)
  install(TARGETS swig_fnm
    EXPORT FnmTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib)
  install(TARGETS swig_fnm
    EXPORT SofusTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib)
endif()
