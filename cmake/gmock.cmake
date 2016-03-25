########################### GMOCK

# We need thread support
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Set default ExternalProject root directory
#set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rd_party)

if (NOT TARGET googlemock)
  ExternalProject_Add(
      googlemock
      URL https://googlemock.googlecode.com/files/gmock-1.7.0.zip
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/3rd_party/gmock

      # Has effect on Windows only (when =OFF remember to force linking of dependencies)
      # CMAKE_ARGS -Dgmock_force_shared_crt=ON
      
      # Disable install step
      INSTALL_COMMAND ""
      
      # Wrap download, configure and build steps in a script to log output
      LOG_DOWNLOAD ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
  )
  
  # Get directories
  ExternalProject_Get_Property(googlemock source_dir binary_dir)
  set(GMOCK_INCLUDE_DIR ${source_dir}/include)
  
  # We cannot search for libraries, which has not been built yet 
  set(GMOCK_TARGETS "gmock" "gmock_main")
  foreach (GMOCK_TARGET ${GMOCK_TARGETS})

    # Add imported library
    add_library(${GMOCK_TARGET} UNKNOWN IMPORTED)
  
    # Set properties
    string(TOUPPER "${GMOCK_TARGET}" _UPPER_GMOCK_TARGET)
    if (MSVC)
      string(TOUPPER "$<CONFIGURATION:googlemock>" _UPPER_CONF)
      # Library name (including path)
      set("${_UPPER_GMOCK_TARGET}_PATH" ${binary_dir}/$<CONFIGURATION:googlemock>/${CMAKE_FIND_LIBRARY_PREFIXES}${GMOCK_TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX})
      set_property(TARGET ${GMOCK_TARGET} PROPERTY IMPORTED_LOCATION_${UPPER_CONF} "${_UPPER_GMOCK_TARGET}_PATH")
    else()
      # Library name (including path)
      set("${_UPPER_GMOCK_TARGET}_PATH" ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}${GMOCK_TARGET}${CMAKE_STATIC_LIBRARY_SUFFIX})
      set_property(TARGET ${GMOCK_TARGET} PROPERTY IMPORTED_LOCATION "${_UPPER_GMOCK_TARGET}_PATH")
      set_property(TARGET ${GMOCK_TARGET} PROPERTY LINK_INTERFACE_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}") # overrides above
    endif()
  
    # Add dependencies
    add_dependencies(${GMOCK_TARGET} googlemock)
  endforeach(GMOCK_TARGET)

  # Doesn't work, since they don't exist yet (see below)
  set_property(TARGET gmock PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${source_dir}/include)

  # Macro for adding dependencies to gtest
  macro(add_gmock_dependencies test)
    add_dependencies(${test} gtest)
    include_directories(${source_dir}/include)
    target_link_libraries(${test} ${GMOCK_PATH})
    target_link_libraries(${test} ${GMOCK_MAIN_PATH})
    if (${UNIX})
      target_link_libraries(${test} pthread)
    endif()
  endmacro(add_gmock_dependencies)
endif()
