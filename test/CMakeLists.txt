# === Dependencies ===
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)
find_package(FFTW REQUIRED)

# === Helper function for test executables ===
function(add_fnm_test name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE
    fnm
    sps
    gl
    ${FFTW_LIBRARIES}
    GTest::GTest
    GTest::Main
    Threads::Threads)

  if(TARGET build)
    target_link_libraries(${name} PRIVATE build)
  endif()

  if(FNM_PULSED_WAVE)
    target_link_libraries(${name} PRIVATE sofus)
  endif()

  if(UNIX)
    target_link_libraries(${name} PRIVATE dl)
  endif()

  add_test(NAME ${name} COMMAND ${name})
endfunction()

# === Main FNM test ===
add_fnm_test(fnm_test main.cpp)

# === Pulsed wave tests ===
if(FNM_PULSED_WAVE)
  add_fnm_test(angle_sum_test angle_sum_test.cpp)
  add_fnm_test(sofus_test sofus_test.cpp)
endif()

# === C test (non-Windows, float only) ===
if(NOT WIN32 AND NOT FNM_Double_support)
  add_executable(fnm_ctest ctest.c)
  target_link_libraries(fnm_ctest PRIVATE
    fnm
    sps
    gl
    ${FFTW_LIBRARIES}
    cmocka)
  if(UNIX)
    target_link_libraries(fnm_ctest PRIVATE dl)
  endif()

  add_executable(test_state ctest.c)
  target_link_libraries(test_state PRIVATE
    fnm
    sps
    gl
    ${FFTW_LIBRARIES}
    cmocka
    dl)
  add_test(NAME test_state COMMAND test_state)
endif()

enable_testing()
