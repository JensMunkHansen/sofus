# Setup installation directories and make relative paths absolute 
macro(setup_installation_directories)

  # Offer the user the choice of overriding the installation directories
  # ====================================================================
  set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
  set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
  set(INSTALL_INCLUDE_DIR include CACHE PATH
    "Installation directory for API header files")
  set(INSTALL_DOC_DIR doc CACHE PATH "Installation directory for documentation")
  if(WIN32 AND NOT CYGWIN)
    set(DEF_INSTALL_CMAKE_DIR CMake)
  else()
    set(DEF_INSTALL_CMAKE_DIR lib/CMake/${PROJECT_NAME})
  endif()
  set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
    "Installation directory for CMake files")
  
  # Make relative paths absolute
  # ============================
  foreach(p LIB BIN INCLUDE DOC CMAKE)
    set(var INSTALL_${p}_DIR)
    if(NOT IS_ABSOLUTE "${${var}}")
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
  endforeach()
endmacro()
