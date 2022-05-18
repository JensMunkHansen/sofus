include(GenerateExportHeader)
include(GenerateProductVersion)

check_include_file(signal.h HAVE_SIGNAL_H)

option (SOFUS_SHARED_LIBS "Build Shared Libraries" OFF)
set(SOFUS_LIB_TYPE STATIC)
if (SOFUS_SHARED_LIBS)
  set(SOFUS_LIB_TYPE SHARED)
endif()

# Two-dimensional computation (Faster if ON)
set (USE_PROJECTIONS 0)
option (SOFUS_USE_Projections "Use projection for computing radii in planes" ON)
if (SOFUS_USE_Projections)
  set(USE_PROJECTIONS 1)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

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

add_library(sofus ${SOFUS_LIB_TYPE} ${sofus_SOURCES} ${sofus_HEADERS} config.h.in)

generate_export_header(sofus)

install(TARGETS sofus
  EXPORT FnmTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/fnm" COMPONENT dev)
