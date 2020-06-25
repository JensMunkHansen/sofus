macro (ADD_DOXYGEN_TARGET Target)
  # Figure this out
  set(CMAKE_DOXYGEN_OUT "${Target}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Figure this out to specify doxy file
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/.tmp")

  # Consider replacing DOXYGEN_EXECUTABLE with unbuffer sh -c "doxygen configfile 2>/proc/$$/fd/2" | animation 
  
  add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/doxygen.stamp  # Was ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} -b ${PROJECT_BINARY_DIR}/Doxyfile
    COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/doxygen.stamp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile # Was empty
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

  # DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html
  # ALL is removed to prevent it from building all the time
  add_custom_target( "${Target}" DEPENDS ${PROJECT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/doxygen.stamp)

  if(WIN32)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/html
      DESTINATION Documentation/${Target}
      COMPONENT Documentation)
  else()
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html
      DESTINATION "${INSTALL_DOC_DIR}/${Target}" COMPONENT dev)
  endif()
  SET(CPACK_NSIS_MENU_LINKS
    ${CPACK_NSIS_MENU_LINKS}
    "Documentation/${Target}html/index.html" "${Target} - (Doxygen doc)")

endmacro (ADD_DOXYGEN_TARGET)
