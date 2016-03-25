########################### NANOMSG

# TODO: Move to 3rd party

# We need thread support
find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Set default ExternalProject root directory
set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/3rd_party)

if (NOT TARGET nanomsg)
  ExternalProject_Add(nanomsg
    #GIT_REPOSITORY git@github.com:nanomsg/nanomsg.git
    URL https://github.com/nanomsg/nanomsg/releases/download/0.6-beta/nanomsg-0.6-beta.tar.gz
    #SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/nanomsg
    #CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/nanomsg/configure --prefix=${CMAKE_BINARY_DIR} 
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
    
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    )
  
  # Get directories
endif()

