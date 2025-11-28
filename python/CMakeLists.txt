option(TEST_INCLUDE_PYTHON_TESTS "Include Python tests" ON)

# Copy DLL and Python wrapper to test directory if Python tests are included
if(TEST_INCLUDE_PYTHON_TESTS)
  find_package(Python3 COMPONENTS Interpreter)

  if(Python3_FOUND)
    add_custom_target(copy_python_test_files ALL
      COMMENT "Copying files for Python testing")

    # Copy libraries to source directory for unit testing (SWIG)
    add_custom_command(TARGET copy_python_test_files POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:gl> ${CMAKE_CURRENT_SOURCE_DIR})

    add_custom_command(TARGET copy_python_test_files POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:fnm> ${CMAKE_CURRENT_SOURCE_DIR})

    add_custom_command(TARGET copy_python_test_files POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:sps> ${CMAKE_CURRENT_SOURCE_DIR})

    if(BUILD_SWIG_INTERFACE AND TARGET swig_fnm)
      add_custom_command(TARGET copy_python_test_files POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:swig_fnm> ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    enable_testing()

    message(STATUS "Python executable: ${Python3_EXECUTABLE}")

    add_test(
      NAME fnm_attr_test
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_attr_test.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(
      NAME fnm_reference_test
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_reference_test.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(
      NAME fnm_simd_test
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fnm_simd_test.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
endif()
