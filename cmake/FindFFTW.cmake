# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

set(FFTW_LIBS_NAMES libfftw3f-3 fftw3f-3 libfftw3l-3 fftw3l-3 libfftw3-3 fftw3-3 libfftw3f fftw3f libfftw3l fftw3l libfftw3 fftw3)

foreach(lib ${FFTW_LIBS_NAMES})
  find_library(
    ${lib}_LIB
    NAMES ${lib}
    PATHS "C:/Program\ Files/fftw-3.3.4"
    PATH_SUFFIXES i386 x86_64
  )
  if (${lib}_LIB)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${${lib}_LIB})
  else()
    unset(${lib}_LIB CACHE)
  endif()
endforeach()

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
