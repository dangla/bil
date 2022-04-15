cmake_minimum_required(VERSION 3.16)

# disable in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR “In-source build detected!”)
endif()

# if CMAKE_BUILD_TYPE is specified use it; otherwise set the default
# build type to "RelWithDebInfo" ("-O2 -g" with gcc) prior to calling
# project()
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else()
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose build type")
endif()


# set the project name
project(bil C CXX Fortran)


# Add the CMake directory for custon CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# Bil full path
# -------------
set(BIL_PATH ${CMAKE_CURRENT_SOURCE_DIR})


# Bil informations
# ----------------
include(SetBilInfos)


# Directories of Bil
# ------------------
set(BIL_LIBDIR ${BIL_PATH}/lib)
set(BIL_BINDIR ${BIL_PATH}/bin)
set(BIL_SRCDIR ${BIL_PATH}/src)


# Installation directories
# ------------------------
include(GNUInstallDirs)
set(BIL_INSTALL_FULL_LIBDIR  ${CMAKE_INSTALL_FULL_LIBDIR})
set(BIL_INSTALL_FULL_BINDIR  ${CMAKE_INSTALL_FULL_BINDIR})
set(BIL_INSTALL_FULL_MAN1DIR ${CMAKE_INSTALL_FULL_MANDIR}/man1)
set(BIL_INSTALL_FULL_DOCDIR  ${CMAKE_INSTALL_FULL_DOCDIR})
set(BIL_INSTALL_FULL_INFODIR ${CMAKE_INSTALL_FULL_INFODIR})


# Library of Bil
# ----------------
set(BIL_LIB bil-${BIL_VERSION}-${CMAKE_BUILD_TYPE})


# Bil executable
# --------------
set(BIL_EXE ${BIL_LIB}.exe)


# Extra-libraries
# ---------------
set(BIL_EXTRALIBS_FILE BilLib.h)
set(BIL_EXTRALIBS)
include(Extralibs)



#[[
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
]]



include(SetCompilerFlags)


# the executable and the library
add_subdirectory(src)

# the documentations
add_subdirectory(doc)



# Custom targets
# --------------


add_custom_target(githelp
    COMMAND echo "States of files: tracked or untracked"
    COMMAND echo "  tracked files: can be unmodified, modified or staged"
    COMMAND echo "  untracked files: anything else"
    COMMAND echo ""
    COMMAND echo "Checking"
    COMMAND echo "  checking the status of files: git status"
    COMMAND echo "  ...for a more compact output: git status -s"
    COMMAND echo ""
    COMMAND echo "Tracking"
    COMMAND echo "  tracking new files: git add [pathtonewfile.c]"
    COMMAND echo ""
    COMMAND echo "Staging"
    COMMAND echo "  staging modified files: git add [pathtofile.c]"
    COMMAND echo "  viewing staged and unstaged changes: git diff"
    COMMAND echo "  viewing what was staged that will go into the next commit: git diff --staged"
    COMMAND echo ""
    COMMAND echo "Committing"
    COMMAND echo "  committing the changes: git commit -m [message]"
    COMMAND echo "  skipping the staging:   git commit -a -m [message]"
    COMMAND echo ""
    COMMAND echo "Viewing"
    COMMAND echo "  viewing the commit history: git log"
    COMMAND echo "  showing the difference introduced in each commit: git log -p"
    COMMAND echo "  prints on a single line: git log --pretty=oneline"
    COMMAND echo ""
    COMMAND echo "Working with remotes"
    COMMAND echo "  showing your remotes: git remote"
    COMMAND echo "  with the URLs that Git has stored for the shortname to be used: git remote -v"
    COMMAND echo "  adding a remote: git remote add [shortname] [url]"
    COMMAND echo ""
    COMMAND echo "Pushing"
    COMMAND echo "  pushing to remote ifsttar: git push ifsttar master"
    COMMAND echo "  pushing to remote dangla: git push dangla master"
    COMMAND echo "  adding a remote: git remote add [shortname] [url]"
    COMMAND echo "  Note: Github will no longer accept account passwords"
    COMMAND echo "  for authenticating Git operations and will require "
    COMMAND echo "  a personal access token."
    COMMAND echo "  Such a personal access token can be generated"
    COMMAND echo "  in the github account by clicking:"
    COMMAND echo "  Settings, Developer settings, Personal access tokens."
    COMMAND echo "  Select the scope \"repo\"."
    COMMAND echo "  Use this personal access token as a password."
    COMMAND echo ""
    COMMAND echo "Tagging"
    COMMAND echo "  listing the available tags: git tag"
    COMMAND echo "  creating annotated tags: git tag -a [version] -m [message]"
    COMMAND echo "  tagging later: git tag -a [version] [checksum]"
    COMMENT "Reminders of git basic commands"
)



# Target rules for memcheck through valgrind.
# To pass the variable "foo" to make, use the command-line: make memcheck arg="foo"
find_program(VALGRIND valgrind)
if(VALGRIND)
  add_custom_target(memcheck
    COMMAND ${VALGRIND} --tool=memcheck --leak-check=full ${BIL_EXE} ${arg}
  )
endif()


# Epilogue
message("")
message("Bil has been configured for the OS: " ${BIL_OS})
message("with the following options: " ${BIL_CONFIG_OPTIONS})
message("")
message("C compiler:    ${CMAKE_C_COMPILER}")
message("C++ compiler:  ${CMAKE_CXX_COMPILER}")
message("Fortran compiler:    ${CMAKE_Fortran_COMPILER}")
message("C compiler flags:    ${CMAKE_C_FLAGS}")
message("C++ compiler flags:  ${CMAKE_CXX_FLAGS}")
message("Fortran compiler flags:  ${CMAKE_Fortran_FLAGS}")
if(CMAKE_BUILD_TYPE)
  message("Build type:  ${CMAKE_BUILD_TYPE}")
endif()
message("Install prefix:  ${CMAKE_INSTALL_PREFIX}")
message("")
message("Run 'ccmake ${CMAKE_CURRENT_SOURCE_DIR}' to fine-tune the configuration.")
message("")
