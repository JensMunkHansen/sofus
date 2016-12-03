macro (ADD_DOXYGEN_TARGET Target)
  set(CMAKE_DOXYGEN_OUT doc)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # New way
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/.tmp")

  add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/${CMAKE_DOXYGEN_OUT}/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

  add_custom_target( doc DEPENDS ${PROJECT_SOURCE_DIR}/doc/html/index.html )

  if(WIN32)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_DOXYGEN_OUT}/html
      DESTINATION Documentation/${Target}
      COMPONENT Documentation)
  else()
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html
      DESTINATION "${INSTALL_DOC_DIR}/${Target}" COMPONENT dev)
  endif()
  SET(CPACK_NSIS_MENU_LINKS
    ${CPACK_NSIS_MENU_LINKS}
    "Documentation/${Target}html/index.html" "${Target} - (Doxygen doc)")

endmacro (ADD_DOXYGEN_TARGET)
