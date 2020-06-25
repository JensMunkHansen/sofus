#find_package( GoogleTest REQUIRED )
find_package( GTest REQUIRED )
find_package( Threads REQUIRED )
find_package( FFTW REQUIRED )

include_directories(${GTEST_INCLUDE_DIRS})

if (MSVC)
  if (ENABLE_AVX)
    add_compile_flags("C;CXX" "/arch:AVX")
    add_definitions("/arch:AVX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  endif()
  if (ENABLE_AVX2)
    add_compile_flags("C;CXX" "/arch:AVX2")
    add_definitions("/arch:AVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
  endif()
endif()

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

add_executable(fnm_test main.cpp)

# Targets are defined in GoogleTest.cmake
target_link_libraries(fnm_test fnm sps gl ${GTEST_LIBRARIES} ${FFTW_LIBRARIES})
if (FNM_PULSED_WAVE)
  target_link_libraries(fnm_test sofus)
endif()

if (BUILD_SHARED_LIBS)
  target_link_libraries(fnm_test dl)
endif()

if (UNIX)
  target_link_libraries(fnm_test dl)
endif()

if (FNM_PULSED_WAVE)
  set(angle_sum_test_HEADERS
  )
  
  set(angle_sum_test_SOURCES
    angle_sum_test.cpp
  )
  
  if (UNIX AND NOT WIN32 AND NOT CYGWIN)
    set(angle_sum_test_HEADERS "${angle_sum_test_HEADERS}" "../sps/strace.hpp")
    set(angle_sum_test_SOURCES "${angle_sum_test_SOURCES}" "../sps/strace.cpp")
  endif()
  
  add_executable(angle_sum_test ${angle_sum_test_HEADERS} ${angle_sum_test_SOURCES})


  target_link_libraries(angle_sum_test ${GTEST_LIBRARIES})

  if (UNIX)
    target_link_libraries(angle_sum_test dl)
  endif()
endif()

if (NOT WIN32)
  if (NOT FNM_Double_support)
    add_executable(fnm_ctest ctest.c)
    target_link_libraries(fnm_ctest fnm sps gl ${FFTW_LIBRARIES} cmocka)
    if (UNIX)
      target_link_libraries(fnm_ctest dl)
    endif()
  endif()
endif()

enable_testing()
add_test(NAME    fnm_test 
         COMMAND fnm_test)

add_test(NAME    angle_sum_test 
         COMMAND angle_sum_test)

       
if (FNM_PULSED_WAVE)
  add_executable(sofus_test sofus_test.cpp)
  target_link_libraries(sofus_test fnm sofus sps gl ${GTEST_LIBRARIES} ${FFTW_LIBRARIES})
  if (UNIX)
    target_link_libraries(sofus_test dl)
  endif()

  add_test(NAME    sofus_test 
    COMMAND sofus_test)
endif()

if (NOT WIN32)
  add_executable(test_state ctest.c)
  add_dependencies(test_state test_state)
  add_test(test_state test_state)
  target_link_libraries(test_state fnm sps gl ${FFTW_LIBRARIES} cmocka dl)
endif()   

# Link to threading library (general)
set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")



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

       
