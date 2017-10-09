include(GenerateExportHeader)
include(GenerateProductVersion)

check_include_file(signal.h HAVE_SIGNAL_H)

# Two-dimensional computation (Faster if ON)
set (USE_PROJECTIONS 0)
option (SOFUS_USE_Projections "Use projection for computing radii in planes" ON)
if (SOFUS_USE_Projections)
  set(USE_PROJECTIONS 1)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)

set(sofus_HEADERS
  rect_int_limits.hpp
  )

set(sofus_SOURCES
  )


