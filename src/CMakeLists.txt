

include(DirsAndSourcesOfBil)
dirs_and_sources_of_bil(BIL_DIRSALL BIL_SRCALL)
#message(STATUS "directories:")
#message("${BIL_DIRSALL}")
#message(STATUS "source c and f files under src within two level depth:")
#message("${BIL_SRCALL}")




#[[
  Building the shared library target
 ]]
 
add_library(${BIL_LIB} SHARED ${BIL_SRCALL})

target_include_directories(${BIL_LIB}
    PRIVATE
    "${BIL_DIRSALL}"
    "${PROJECT_SOURCE_DIR}"
)

target_compile_features(${BIL_LIB}
    PUBLIC c_std_11 cxx_std_11
)

set_target_properties(${BIL_LIB}
  PROPERTIES 
  #INSTALL_RPATH "${BIL_INSTALL_LIBDIR}"
  LIBRARY_OUTPUT_DIRECTORY ${BIL_LIBDIR}
)

install(TARGETS ${BIL_LIB}
  DESTINATION ${BIL_INSTALL_LIBDIR}
  LIBRARY OPTIONAL)



#[[
  Building the executable target
 ]]

add_executable(${BIL_EXE} Main/Main.c)

target_include_directories(${BIL_EXE}
    PRIVATE
    "${BIL_DIRSALL}"
    "${PROJECT_SOURCE_DIR}"
)

target_compile_features(${BIL_EXE}
    PUBLIC c_std_11 cxx_std_11
)

#target_compile_options(${BIL_EXE}
#    PUBLIC
#    -pedantic ${WFLAGS} ${OFLAGS}
#)

set_target_properties(${BIL_EXE}
  PROPERTIES 
  INSTALL_RPATH "${BIL_INSTALL_LIBDIR}"
  RUNTIME_OUTPUT_DIRECTORY ${BIL_BINDIR}
)

target_link_libraries(${BIL_EXE} ${BIL_LIB} ${BIL_EXTRALIBS})

target_link_options(${BIL_EXE}
  PUBLIC "-O3 -g -gdwarf-3"
)

install(TARGETS ${BIL_EXE}
  DESTINATION ${BIL_INSTALL_BINDIR}
  RUNTIME OPTIONAL)



install(CODE
  "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${BIL_INSTALL_BINDIR}/${BIL_EXE} ${BIL_INSTALL_BINDIR}/bil)"
)
 




