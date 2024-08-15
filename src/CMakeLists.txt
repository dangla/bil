
include(DirsAndSourcesOfBil)
dirs_and_sources_of_bil(BIL_DIRSALL BIL_SRCALL BIL_SRCMAIN)



#[[
  Building the shared library target
#]]
 
add_library(${BIL_LIB} SHARED ${BIL_SRCALL})

target_include_directories(${BIL_LIB}
    PRIVATE
    "${BIL_DIRSALL}"
    "${PROJECT_SOURCE_DIR}"
)

# autodiff requires a c++17 enabled compiler
target_compile_features(${BIL_LIB}
    PUBLIC c_std_11 cxx_std_17
)

set_target_properties(${BIL_LIB}
  PROPERTIES 
  INSTALL_RPATH "${BILEXTRALIBS_INSTALL_RPATH}"
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
include(CreateBilHeader)
create_bil_header("${BIL_HDRALL}" "${BIL_INCLUDEDIR}/bil.h")
install(FILES "${BIL_INCLUDEDIR}/bil.h" DESTINATION ${BIL_INSTALL_FULL_INCLUDEDIR})


#[[
  Building the executable target
#]]

#add_executable(${BIL_EXE} Main/Main.c Main/Main.cpp)
add_executable(${BIL_EXE} ${BIL_SRCMAIN})

target_include_directories(${BIL_EXE}
    PRIVATE
    "${BIL_DIRSALL}"
    "${PROJECT_SOURCE_DIR}"
)

target_compile_features(${BIL_EXE}
    PUBLIC c_std_11 cxx_std_17
)

#target_compile_options(${BIL_EXE}
#    PUBLIC
#    -pedantic ${WFLAGS} ${OFLAGS}
#)


set_target_properties(${BIL_EXE}
  PROPERTIES 
  INSTALL_RPATH "${BIL_INSTALL_FULL_LIBDIR};${BILEXTRALIBS_INSTALL_RPATH}"
  RUNTIME_OUTPUT_DIRECTORY ${BIL_BINDIR}
)

target_link_libraries(${BIL_EXE} ${BIL_LIB} ${BIL_EXTRALIBS})

if(HAVE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    Message(STATUS "Compiler supports OpenMP.")
    target_link_libraries(${BIL_EXE} OpenMP::OpenMP_C)
  else()
    Message(STATUS "Compiler doesn't support OpenMP.")
  endif()
else()
  string(FIND "${BIL_EXTRALIBS}" "OPENMP" openmp)
  if(${openmp} GREATER -1)
    Message(FATAL_ERROR "Please enable OpenMP.")
  endif()
endif()


if(HAVE_PTHREAD)
  find_package(Threads)
  if(Threads_FOUND)
    Message(STATUS "Compiler supports Pthread.")
    target_link_libraries(${BIL_EXE} Threads::Threads)
  else()
    Message(STATUS "Compiler doesn't support Pthread.")
  endif()
else()
  string(FIND "${BIL_EXTRALIBS}" "PTHREAD" pthread)
  if(${pthread} GREATER -1)
    Message(FATAL_ERROR "Please enable Pthread.")
  endif()
endif()


if(HAVE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    Message(STATUS "Compiler supports MPI.")
    #target_link_libraries(${BIL_EXE} MPI::MPI_C)
    target_link_libraries(${BIL_EXE} MPI::MPI_CXX)
  else()
    Message(STATUS "Compiler doesn't support MPI.")
  endif()
else()
  string(FIND "${BIL_EXTRALIBS}" "dist" dist)
  string(FIND "${BIL_EXTRALIBS}" "libpetsc" libpetsc)
  if(${dist} GREATER -1 OR ${libpetsc} GREATER -1)
    Message(FATAL_ERROR "Please enable MPI.")
  endif()
endif()


if(HAVE_AUTODIFF)
  UNSET(AUTODIFF_PATH CACHE)
  find_path(AUTODIFF_PATH forward/real.hpp PATHS ${CMAKE_INSTALL_FULL_INCLUDEDIR}/autodiff)
  if(AUTODIFF_PATH)
    Message(STATUS "Compiler supports autodiff.")
  else()
    Message(FATAL ERROR "Compiler doesn't support autodiff.")
  endif()
  
  find_package(Eigen3 REQUIRED)
  if(Eigen3_FOUND)
    Message(STATUS "Compiler supports Eigen3.")
  else()
    Message(FATAL ERROR "Compiler doesn't support Eigen3.")
  endif()
endif()


include(FeatureSummary)
feature_summary(INCLUDE_QUIET_PACKAGES DESCRIPTION "Found packages:" WHAT ALL)



target_link_options(${BIL_EXE}
  PUBLIC "-O3 -g -gdwarf-3"
)

install(TARGETS ${BIL_EXE}
  DESTINATION ${BIL_INSTALL_FULL_BINDIR}
  RUNTIME OPTIONAL)



install(CODE
  "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${BIL_INSTALL_FULL_BINDIR}/${BIL_EXE} ${BIL_INSTALL_FULL_BINDIR}/bil)"
)
 




