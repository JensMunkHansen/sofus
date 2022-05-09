cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

option (TEST_INCLUDE_PYTHON_TESTS "Include Python tests" ON)

# Copy DLL and Python wrapper to test directory if Python test are included
if (TEST_INCLUDE_PYTHON_TESTS)

  find_package(Python3 COMPONENTS Interpreter)

  add_custom_target(copy ALL)

  # Copy DLL/so to source directory for unit testing (SWIG)
  add_custom_command(TARGET copy
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:gl> ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS $<TARGET_FILE:gl>)

  add_custom_command(TARGET copy
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:fnm> ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS $<TARGET_FILE:fnm>)

  add_custom_command(TARGET copy
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sps> ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS $<TARGET_FILE:sps>)

  if (BUILD_SWIG_INTERFACE)
    add_custom_command(TARGET copy
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:swig_fnm> ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS $<TARGET_FILE:swig_fnm>)
  endif()

  enable_testing()
  
  message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")

  add_test(
    NAME    "${CMAKE_CURRENT_SOURCE_DIR}/fnm_attr_test.py"
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_attr_test.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(
    NAME    "${CMAKE_CURRENT_SOURCE_DIR}/fnm_reference_test.py"
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_reference_test.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(
    NAME    "${CMAKE_CURRENT_SOURCE_DIR}/fnm_simd_test.py"
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_simd_test.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  
endif(TEST_INCLUDE_PYTHON_TESTS)
  
