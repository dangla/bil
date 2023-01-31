
include(DirsAndSourcesOfBil)
dirs_and_sources_of_bil(BIL_DIRSALL BIL_SRCALL)



#[[
  Building the shared library target
#]]
 
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
  #INSTALL_RPATH "${BIL_INSTALL_FULL_LIBDIR}"
  LIBRARY_OUTPUT_DIRECTORY ${BIL_LIBDIR}
)

install(TARGETS ${BIL_LIB}
  DESTINATION ${BIL_INSTALL_FULL_LIBDIR}
  LIBRARY OPTIONAL)



#[[
  Install the headers
#]]

include(AppendBilFiles)
append_bil_files("${BIL_DIRSALL}" "h" BIL_HDRALL)
install(FILES ${BIL_HDRALL} DESTINATION ${BIL_INSTALL_FULL_INCLUDEDIR})



#[[
  Building the executable target
#]]

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
  INSTALL_RPATH "${BILEXTRALIBS_INSTALL_RPATH}"
  RUNTIME_OUTPUT_DIRECTORY ${BIL_BINDIR}
)

target_link_libraries(${BIL_EXE} ${BIL_LIB} ${BIL_EXTRALIBS})

string(FIND "${BIL_EXTRALIBS}" "OPENMP" mtapi)
if(HAVE_OPENMP OR ${mtapi} GREATER -1)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    Message(STATUS "Compiler supports OpenMP.")
    target_link_libraries(${BIL_EXE} OpenMP::OpenMP_C)
  else()
    Message(STATUS "Compiler doesn't support OpenMP.")
  endif()
endif()


string(FIND "${BIL_EXTRALIBS}" "PTHREAD" mtapi)
if(HAVE_PTHREAD OR ${mtapi} GREATER -1)
  find_package(Threads)
  if(Threads_FOUND)
    Message(STATUS "Compiler supports Pthread.")
    target_link_libraries(${BIL_EXE} Threads::Threads)
  else()
    Message(STATUS "Compiler doesn't support Pthread.")
  endif()
endif()


string(FIND "${BIL_EXTRALIBS}" "dist" mtapi)
if(HAVE_MPI OR ${mtapi} GREATER -1)
  find_package(MPI)
  if(MPI_FOUND)
    Message(STATUS "Compiler supports MPI.")
    target_link_libraries(${BIL_EXE} MPI::MPI_C)
  else()
    Message(STATUS "Compiler doesn't support MPI.")
  endif()
endif()



target_link_options(${BIL_EXE}
  PUBLIC "-O3 -g -gdwarf-3"
)

install(TARGETS ${BIL_EXE}
  DESTINATION ${BIL_INSTALL_FULL_BINDIR}
  RUNTIME OPTIONAL)



install(CODE
  "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${BIL_INSTALL_FULL_BINDIR}/${BIL_EXE} ${BIL_INSTALL_FULL_BINDIR}/bil)"
)
 




