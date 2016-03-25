macro (ADD_DOXYGEN_TARGET Target)
  set(CMAKE_DOXYGEN_OUT DoxygenDocs)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  add_custom_command(OUTPUT tmp_dir
    POST_BUILD
    COMMAND cmake -E make_directory "${PROJECT_BINARY_DIR}/.tmp")

  add_custom_target(doc
    DEPENDS tmp_dir
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

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
