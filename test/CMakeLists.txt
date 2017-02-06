find_package( GoogleTest REQUIRED )
find_package( Threads REQUIRED )
find_package( FFTW REQUIRED )

include_directories(${GTEST_INCLUDE_DIRS})

add_executable(fnm_test main.cpp)

# Targets are defined in GoogleTest.cmake
target_link_libraries(fnm_test fnm sps gl gtest ${FFTW_LIBRARIES})
if (FNM_PULSED_WAVE)
  target_link_libraries(fnm_test sofus)
endif()

if (BUILD_SHARED_LIBS)
  target_link_libraries(fnm_test dl)
endif()

# Link to threading library (general)
set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")

enable_testing()
add_test(NAME    fnm_test 
         COMMAND fnm_test)


# Not used

if (0)
  # Alternatively search for libraries
  
  # List of external libraries used
  set(FNM_TEST_LIBS_NAMES gtest gtest_main)
  
  # Search for libraries (binaries already present)
  foreach(lib ${FNM_TEST_LIBS_NAMES})
    find_library(
      ${lib}_LIB
      NAMES ${lib} 
      PATHS ${CMAKE_PREFIX_PATH}/lib
      PATH_SUFFIXES i386 x86_64
    )
    set(FNM_TEST_LIBS ${FNM_TEST_LIBS} ${${lib}_LIB})
  endforeach()
  
  target_link_libraries(fnm_test ${FNM_TEST_LIBS})
endif()

       
