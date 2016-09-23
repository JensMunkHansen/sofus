# TODO: Generate GoogleTestConfig.cmake for the install tree.
cmake_minimum_required(VERSION 3.0)

project(Fnm_deps
  VERSION 0.1.0.0)

# We need thread support
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rd_party)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})

if (NOT TARGET googletest)
  ExternalProject_Add(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git

    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=Debug
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=Release
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL:PATH=MinSizeRel
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO:PATH=RelWithDebInfo
	       -Dgtest_force_shared_crt=ON # Has effect on Windows only (when =OFF remember to force linking of dependencies)
    
    # Alternatively use DGTEST_LINKED_AS_SHARED_LIBRARY
      
    INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
      
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    )

  # GoogleTest does not generate GoogleTestConfig.cmake (we do)
  
  ExternalProject_Get_Property(googletest source_dir binary_dir)

  set(CONF_INCLUDE_DIRS ${source_dir}/googletest/include)
  set(CONF_BINARY_DIR ${binary_dir})

  configure_file(GoogleTestConfig.cmake.in
    "${PROJECT_BINARY_DIR}/GoogleTestConfig.cmake" @ONLY)

endif()

