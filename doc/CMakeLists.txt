
# Langage (fr or en)
set(LANG en)
set(BIL_DOC bil)


# Subdirectories of doc
set(BIL_EDITION_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/edition)
set(BIL_OXYGEN_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/oxygen)



file(WRITE ${BIL_EDITION_DIR}/versions.texi "@set BIL-VERSION   ${BIL_VERSION}")


#=======================================================================
# Target rules for the reference manual in different formats
if(UNIX AND NOT CYGWIN)
  install(FILES ${BIL_DOC}.1 DESTINATION ${BIL_INSTALL_FULL_MAN1DIR})
endif()

add_custom_target(doc DEPENDS pdf txt html info doxy)

set_directory_properties(PROPERTIES
  ADDITIONAL_CLEAN_FILES "${BIL_EDITION_DIR}/${BIL_DOC}.aux;${BIL_EDITION_DIR}/${BIL_DOC}.cp;${BIL_EDITION_DIR}/${BIL_DOC}.cps;${BIL_EDITION_DIR}/${BIL_DOC}.html;${BIL_EDITION_DIR}/${BIL_DOC}.info;${BIL_EDITION_DIR}/${BIL_DOC}.log;${BIL_EDITION_DIR}/${BIL_DOC}.pdf;${BIL_EDITION_DIR}/${BIL_DOC}.toc;${BIL_EDITION_DIR}/${BIL_DOC}.tp;${BIL_EDITION_DIR}/${BIL_DOC}.txt;${BIL_OXYGEN_DIR}/${BIL_DOC}.doxy;${CMAKE_CURRENT_SOURCE_DIR}/${BIL_DOC}.pdf;${CMAKE_CURRENT_SOURCE_DIR}/${BIL_DOC}.txt;${CMAKE_CURRENT_SOURCE_DIR}/${BIL_DOC}.html"
)



find_program(TEXI2PDF texi2pdf)
if(TEXI2PDF)
  add_custom_command(OUTPUT  ${BIL_DOC}.pdf 
    DEPENDS ${BIL_EDITION_DIR}/${BIL_DOC}.texi
    #COMMAND ${CMAKE_COMMAND} -E remove -f versions.texi
    #COMMAND ${CMAKE_COMMAND} -E echo "@set BIL-VERSION   ${BIL_VERSION}" > versions.texi
    COMMAND ${TEXI2PDF} ${BIL_DOC}.texi
    COMMAND ${CMAKE_COMMAND} -E copy ${BIL_DOC}.pdf ..
    WORKING_DIRECTORY ${BIL_EDITION_DIR})
  add_custom_target(pdf DEPENDS ${BIL_DOC}.pdf)
  install(FILES ${BIL_EDITION_DIR}/${BIL_DOC}.pdf DESTINATION ${BIL_INSTALL_FULL_DOCDIR} OPTIONAL)
endif()



find_program(MAKEINFO makeinfo)
if(MAKEINFO)
  add_custom_command(OUTPUT ${BIL_DOC}.info 
    DEPENDS ${BIL_EDITION_DIR}/${BIL_DOC}.texi 
    COMMAND ${MAKEINFO} --enable-encoding -o ${BIL_DOC}.info ${BIL_DOC}.texi
    WORKING_DIRECTORY ${BIL_EDITION_DIR})
  add_custom_target(info DEPENDS ${BIL_DOC}.info)
  install(FILES ${BIL_EDITION_DIR}/${BIL_DOC}.info DESTINATION ${BIL_INSTALL_FULL_INFODIR} OPTIONAL)
  
  add_custom_command(OUTPUT ${BIL_DOC}.txt 
    DEPENDS ${BIL_EDITION_DIR}/${BIL_DOC}.texi
    COMMAND ${MAKEINFO} --no-headers -o ${BIL_DOC}.txt ${BIL_DOC}.texi 
    COMMAND ${CMAKE_COMMAND} -E copy ${BIL_DOC}.txt ..
    WORKING_DIRECTORY ${BIL_EDITION_DIR})
  add_custom_target(txt DEPENDS ${BIL_DOC}.txt)
  
  add_custom_command(OUTPUT ${BIL_DOC}.html 
    DEPENDS ${BIL_EDITION_DIR}/${BIL_DOC}.texi 
    COMMAND ${MAKEINFO} --html --no-split -o ${BIL_DOC}.html ${BIL_DOC}.texi 
    COMMAND ${CMAKE_COMMAND} -E copy ${BIL_DOC}.html ..
    WORKING_DIRECTORY ${BIL_EDITION_DIR})
  add_custom_target(html DEPENDS ${BIL_DOC}.html)
  install(FILES ${BIL_EDITION_DIR}/${BIL_DOC}.html DESTINATION ${BIL_INSTALL_FULL_DOCDIR} OPTIONAL)
endif()




find_program(DOXYGEN doxygen)
if(DOXYGEN)
  set(BIL_OXYGEN_FILE bil-doxy)
  
  #[[
  file(COPY Doxyfile.def DESTINATION ${BIL_OXYGEN_DIR}/${BIL_OXYGEN_FILE})
  file(APPEND ${BIL_OXYGEN_FILE} 
    "PROJECT_NAME     = Bil\n
     PROJECT_NUMBER   = ${BIL_VERSION}\n
     OUTPUT_DIRECTORY = ${BIL_OXYGEN_DIR}\n
     INPUT            = ${BIL_PATH}/src\n
     HTML_OUTPUT      = ${BIL_OXYGEN_DIR}/${BIL_DOC}.doxy")
  ]]
  
  add_custom_command(OUTPUT ${BIL_DOC}.doxy
    COMMAND ${CMAKE_COMMAND} -E remove -f ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy Doxyfile.def ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "PROJECT_NAME     = Bil" >> ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "PROJECT_NUMBER   = ${BIL_VERSION}" >> ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "OUTPUT_DIRECTORY = ${BIL_OXYGEN_DIR}" >> ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "INPUT            = ${BIL_PATH}/src" >> ${BIL_OXYGEN_FILE}
    COMMAND ${CMAKE_COMMAND} -E echo "HTML_OUTPUT      = ${BIL_OXYGEN_DIR}/${BIL_DOC}.doxy" >> ${BIL_OXYGEN_FILE}
    COMMAND ${DOXYGEN} ${BIL_OXYGEN_FILE}
    WORKING_DIRECTORY ${BIL_OXYGEN_DIR})
  add_custom_target(doxy DEPENDS ${BIL_DOC}.doxy)
endif()

	
	

