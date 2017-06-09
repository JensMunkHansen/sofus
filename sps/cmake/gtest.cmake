########################### GTEST

# We need thread support
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rd_party)

if (NOT TARGET googletest)
  ExternalProject_Add(
      googletest
      DOWNLOAD_COMMAND ""
      SOURCE_DIR ${PROJECT_SOURCE_DIR}/3rd_party/gtest-1.7.0
  #    PREFIX ${CMAKE_BINARY_DIR}/3rd_party/gtest
  #    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
  #    SVN_REVISION -r "663"
      
      # Has effect on Windows only (when =OFF remember to force linking of dependencies)
      CMAKE_ARGS -Dgtest_force_shared_crt=ON
      # Alternatively use DGTEST_LINKED_AS_SHARED_LIBRARY
      
      # Disable install step
      INSTALL_COMMAND ""
      
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  
  # Get directories
  ExternalProject_Get_Property(googletest source_dir binary_dir)
  set(GTEST_INCLUDE_DIR ${source_dir}/include)
  
  # We cannot search for libraries, which has not been built yet 
  set(GTEST_TARGETS "gtest" "gtest_main")
  foreach (GTEST_TARGET ${GTEST_TARGETS})
    # Add imported library
    add_library(${GTEST_TARGET} UNKNOWN IMPORTED)
  
    # Set properties
    string(TOUPPER "${GTEST_TARGET}" _UPPER_GTEST_TARGET)
    if (MSVC)
      string(TOUPPER "$<CONFIGURATION:googletest>" _UPPER_CONF)
      # Library name (including path)
      set("${_UPPER_GTEST_TARGET}_PATH" ${binary_dir}/$<CONFIGURATION:googletest>/${CMAKE_FIND_LIBRARY_PREFIXES}${GTEST_TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX})
      set_property(TARGET ${GTEST_TARGET} PROPERTY IMPORTED_LOCATION_${UPPER_CONF} "${_UPPER_GTEST_TARGET}_PATH")
    else()
      # Library name (including path)
      set("${_UPPER_GTEST_TARGET}_PATH" ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}${GTEST_TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX})
      set_property(TARGET ${GTEST_TARGET} PROPERTY IMPORTED_LOCATION "${_UPPER_GTEST_TARGET}_PATH")
  #    set_property(TARGET ${GTEST_TARGET} PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
  #    set_property(TARGET ${GTEST_TARGET} PROPERTY INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}") # replaces above
      set_property(TARGET ${GTEST_TARGET} PROPERTY LINK_INTERFACE_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}") # overrides above
    endif()
  
    # Add dependencies
    add_dependencies(${GTEST_TARGET} googletest)
  endforeach(GTEST_TARGET)

  # Doesn't work, since they don't exist yet (see below)
  set_property(TARGET gtest PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include)

  # Macro for adding dependencies to gtest
  macro(add_gtest_dependencies test)
    add_dependencies(${test} gtest)
    include_directories(${source_dir}/include)
    target_link_libraries(${test} ${GTEST_PATH})
    target_link_libraries(${test} ${GTEST_MAIN_PATH})
    if (${UNIX})
      target_link_libraries(${test} pthread)
      if (BUILD_SHARED_LIBS)
	target_link_libraries(${test} dl)
      endif()
    endif()
  endmacro(add_gtest_dependencies)
endif()

