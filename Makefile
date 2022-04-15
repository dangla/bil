#
SHELL = /bin/sh

# Installation folders
PREFIX     := /usr/local
BINDIR     := ${PREFIX}/bin
DATADIR    := ${PREFIX}/share
INCLUDEDIR := ${PREFIX}/include
LIBDIR     := ${PREFIX}/lib


# Bil path
BIL_PATH_FILE := BilPath.h
BIL_PATH := ${shell pwd}


# Bil version
BIL_VERSION_FILE  := BilVersion.h
# BIL_VERSION       := ${word 2,${subst /bil-, ,${BIL_PATH}}}
BIL_VERSION       := ${shell cat ${BIL_PATH}/VERSION}
BIL_MAJOR_VERSION := ${word 1,${subst ., ,${BIL_VERSION}}}
BIL_MINOR_VERSION := ${word 2,${subst ., ,${BIL_VERSION}}}
BIL_PATCH_VERSION := ${word 3,${subst ., ,${BIL_VERSION}}}


# Extension for exec files
EXEEXT := .exe


# Executable
BIL_EXE := bil-${BIL_VERSION}${EXEEXT}


# General informations
BIL_INFO_FILE := BilInfo.h
BIL_SHORT_LICENSE := "GNU General Public License"
BIL_DATE      := ${shell date}
BIL_YEAR      := ${shell date +%Y}
BIL_HOST      := ${shell hostname}: ${shell hostname -I}
BIL_PACKAGER  := ${shell whoami}
BIL_OS        := ${shell uname -sr}
BIL_URL       := "http://bil.ifsttar.fr"
BIL_EMAIL     := "patrick.dangla@univ-eiffel.fr"
#BIL_COPYRIGHT := "Copyright \(C\) 2002-"${BIL_YEAR}" Patrick Dangla"
BIL_COPYRIGHT := "Copyright \(C\) 2002 Patrick Dangla"
BIL_PROGNAME  := "Bil, a modeling platform based on FEM/FVM"


# Source directories of Bil and HSL
#BIL_DIRS := src hsl
BIL_DIRS := src


# Directories of Bil
BIL_LIBDIR := ${BIL_PATH}/lib
BIL_BINDIR := ${BIL_PATH}/bin
BIL_SRCDIR := ${BIL_PATH}/src


# Exportations
export PREFIX
export BIL_VERSION
export BIL_PATH
export BIL_LIBDIR


# Libraries of Bil and HSL
#BIL_LIBS := -Llib -lbil -lhsl
#BIL_LIBS := -lbil-${BIL_VERSION} -lhsl
BIL_LIBS := -lbil-${BIL_VERSION}


# Directive file for the loading of extra-libraries
BIL_LIB_FILE  := BilLib.h


# Do we use extra-libraries?
include make.extralibs
BIL_LIBS += ${BIL_EXTRALIBS}


# Compiler commands
include make.inc


#=======================================================================
# Default target executed when no arguments are given to make.
all: bin doc


#=======================================================================
# Target rules for binaries

.PHONY: bin
bin: make.inc path_init version_init info_init lib compile link

compile: make.inc path_init version_init info_init lib
	@echo "\nCompilation"
	mkdir -p ${BIL_LIBDIR}
	for i in ${BIL_DIRS}; do (cd $$i && ${MAKE}); done

link: make.inc
	@echo "\nLinking with embedded location of shared libraries (rpath method), if any"
	mkdir -p ${BIL_BINDIR}
	${LINKER} ${OPTIM} ${BIL_SRCDIR}/Main/Main.o -L${BIL_LIBDIR} -Wl,-rpath=${BIL_LIBDIR} ${BIL_LIBS} -o ${BIL_BINDIR}/${BIL_EXE} ${LFLAGS}


#=======================================================================
# Target rules for installation

install: install-bin install-doc

install-bin: install-lib${LIBBILEXT}
	@echo "\nInstalling "${BIL_EXE}" in "${BINDIR}
	mkdir -p ${BINDIR}
	cp -f ${BIL_BINDIR}/${BIL_EXE} ${BINDIR}
	ln -sf ${BINDIR}/${BIL_EXE} ${BINDIR}/bil
	chmod 555 ${BINDIR}/${BIL_EXE}
	chmod 555 ${BINDIR}/bil
	@echo "\nChecking shared libraries required by Bil"
	ldd ${BINDIR}/bil
	@echo "\nChecking the shared library required by Bil (objdump)"
	#objdump -p ${BINDIR}/bil

install-libso:
	@echo "\nInstalling the library libbil-${BIL_VERSION}.so in "${LIBDIR}
	mkdir -p ${LIBDIR}
	cp -f ${BIL_LIBDIR}/libbil-${BIL_VERSION}.so ${LIBDIR}
	chmod 555 ${LIBDIR}/libbil-${BIL_VERSION}.so
	@echo "\nCreating symbolic links and setting up the cache file (ld.so.cache)"
	${LDCONFIG}
	@echo "\nLinking with default location of shared library "
	${LINKER} ${OPTIM} ${BIL_SRCDIR}/Main/Main.o -L${BIL_LIBDIR} ${BIL_LIBS} -o ${BIL_BINDIR}/${BIL_EXE} ${LFLAGS}
	@echo "\nChecking shared libraries required by libbil-${BIL_VERSION}.so"
	ldd ${LIBDIR}/libbil-${BIL_VERSION}.so
	@echo "\nChecking shared libraries required by libbil-${BIL_VERSION}.so (objdump)"
	#objdump -p ${LIBDIR}/libbil-${BIL_VERSION}.so

install-liba:
#	@echo "\nNo library libbil-${BIL_VERSION}.a to install in "${LIBDIR}"!"

install-doc:
	@echo "\nInstalling documentations"
	cd doc && ${MAKE} install-doc


#=======================================================================
# Target rules for cleaning

clean:
	@echo "\nCleaning directories (removing binaries,...)"
	for i in ${BIL_DIRS} doc examples base; do (cd $$i && ${MAKE} clean); done
	rm -f make.log
	@echo "\nRemoving libraries"
	rm -f ${BIL_LIBDIR}/*.a
	rm -f ${BIL_LIBDIR}/*.so
	find . -name "*~" -printf "\"%p\"\n" | xargs rm -f


clean-all: clean
	rm -f ${BIL_BINDIR}/${BIL_EXE}
	rm -f ${BIL_BINDIR}/*
	rm -f ${BIL_VERSION_FILE} 
	rm -f ${BIL_INFO_FILE} 
	rm -f ${BIL_LIB_FILE}
	rm -f ${BIL_PATH_FILE}
	cd doc && ${MAKE} clean-all


#=======================================================================
# Target rules for building path file

path:
	rm -f ${BIL_PATH_FILE}
	echo "#define BIL_PATH \"${BIL_PATH}\"" >  ${BIL_PATH_FILE}

path_init:
	@if [ ! -r ${BIL_PATH_FILE} ]; then ${MAKE} path ; fi


#=======================================================================
# Target rules for building version file

.PHONY: version
version:
	rm -f ${BIL_VERSION_FILE}
	echo "#define BIL_MAJOR_VERSION ${BIL_MAJOR_VERSION}" >  ${BIL_VERSION_FILE}
	echo "#define BIL_MINOR_VERSION ${BIL_MINOR_VERSION}" >> ${BIL_VERSION_FILE}
	echo "#define BIL_PATCH_VERSION ${BIL_PATCH_VERSION}" >> ${BIL_VERSION_FILE}
	echo "#define BIL_VERSION   \"${BIL_VERSION}\""       >> ${BIL_VERSION_FILE}

version_init:
	@if [ ! -r ${BIL_VERSION_FILE} ]; then ${MAKE} version ; fi


#=======================================================================
# Target rules for building info file

info:
	rm -f ${BIL_INFO_FILE}
	echo "#define BIL_PROGNAME  \"${BIL_PROGNAME}\""   > ${BIL_INFO_FILE}
	echo "#define BIL_COPYRIGHT \"${BIL_COPYRIGHT}\"" >> ${BIL_INFO_FILE}
	echo "#define BIL_DATE      \"${BIL_DATE}\""      >> ${BIL_INFO_FILE}
	echo "#define BIL_HOST      \"${BIL_HOST}\""      >> ${BIL_INFO_FILE}
	echo "#define BIL_PACKAGER  \"${BIL_PACKAGER}\""  >> ${BIL_INFO_FILE}
	echo "#define BIL_OS        \"${BIL_OS}\""        >> ${BIL_INFO_FILE}
	echo "#define BIL_SHORT_LICENSE \"${BIL_SHORT_LICENSE}\"" >> ${BIL_INFO_FILE}
	echo "#define BIL_URL       \"${BIL_URL}\""       >> ${BIL_INFO_FILE}
	echo "#define BIL_EMAIL     \"${BIL_EMAIL}\""     >> ${BIL_INFO_FILE}

info_init:
	@if [ ! -r ${BIL_INFO_FILE} ]; then ${MAKE} info ; fi


#=======================================================================
# Target rules for building extra library file

.PHONY: lib
lib:
	rm -f ${BIL_LIB_FILE}
	touch ${BIL_LIB_FILE} ;
	@if [ ${SLU_USE} = "YES" ] ; then \
		echo "#define SUPERLULIB  ${SLU_DIR}"    >>  ${BIL_LIB_FILE} ; \
	fi
	@if [ ${BLAS_USE} = "YES" ] ; then \
		echo "#define BLASLIB     ${BLAS_DIR}"   >> ${BIL_LIB_FILE} ; \
	fi
	@if [ ${LAPACK_USE} = "YES" ] ; then \
		echo "#define LAPACKLIB   ${LAPACK_DIR}" >> ${BIL_LIB_FILE} ; \
	fi


#=======================================================================
# Target rules for documentations

.PHONY: doc
doc:
	( ${MAKE} -C doc )


#=======================================================================
# Target rules for building some solution references

.PHONY: base
base:
	( ${MAKE} -C base )


#=======================================================================
# Target rules for archiving

targz: tar zip
	
tar: clean-all
	cd .. && tar cvf bil-${BIL_VERSION}-src.tar ${notdir ${BIL_PATH}}

zip:
	cd .. && gzip bil-${BIL_VERSION}-src.tar
	

#=======================================================================
# Target rules for githelp

githelp:
	@( echo "States of files: tracked or untracked" )
	@( echo "  tracked files: can be unmodified, modified or staged" )
	@( echo "  untracked files: anything else" )
	@( echo "" )
	@( echo "Checking" )
	@( echo "  checking the status of files: git status" )
	@( echo "  ...for a more compact output: git status -s" )
	@( echo "" )
	@( echo "Tracking" )
	@( echo "  tracking new files: git add <pathtonewfile.c>" )
	@( echo "" )
	@( echo "Staging" )
	@( echo "  staging modified files: git add <pathtofile.c>" )
	@( echo "  viewing staged and unstaged changes: git diff" )
	@( echo "  viewing what was staged that will go into the next commit: git diff --staged" )
	@( echo "" )
	@( echo "Committing" )
	@( echo "  committing the changes: git commit -m <message>" )
	@( echo "  skipping the staging:   git commit -a -m <message>" )
	@( echo "" )
	@( echo "Viewing" )
	@( echo "  viewing the commit history: git log" )
	@( echo "  showing the difference introduced in each commit: git log -p" )
	@( echo "  prints on a single line: git log --pretty=oneline" )
	@( echo "" )
	@( echo "Working with remotes" )
	@( echo "  showing your remotes: git remote" )
	@( echo "  with the URLs that Git has stored for the shortname to be used: git remote -v" )
	@( echo "  adding a remote: git remote add <shortname> <url>" )
	@( echo "" )
	@( echo "Pushing" )
	@( echo "  pushing to remote ifsttar: git push ifsttar master" )
	@( echo "  pushing to remote dangla: git push dangla master" )
	@( echo "  adding a remote: git remote add <shortname> <url>" )
	@( echo "  Note: Github will no longer accept account passwords" )
	@( echo "  for authenticating Git operations and will require " )
	@( echo "  a personal access token." )
	@( echo "  Such a personal access token can be generated" )
	@( echo "  in the github account by clicking:" )
	@( echo "  Settings, Developer settings, Personal access tokens.")
	@( echo "  Select the scope \"repo\".")
	@( echo "  Use this personal access token as a password.")
	@( echo "" )
	@( echo "Tagging" )
	@( echo "  listing the available tags: git tag" )
	@( echo "  creating annotated tags: git tag -a <version> -m <message>" )
	@( echo "  tagging later: git tag -a <version> <checksum>" )
	

#=======================================================================
# Target rules for tests

test:
	@( cd base; pwd )
	@( echo "BIL_PATH = ${BIL_PATH}" )
	@( echo "Last part of BIL_PATH = ${notdir ${BIL_PATH}}" )
	@( echo "BIL_COPYRIGHT = ${BIL_COPYRIGHT}" )
	@( echo "BIL_LIBS = ${BIL_LIBS}" )
	

#=======================================================================
# Target rules for memcheck through valgrind
# To pass the variable "foo" to make use: make memcheck arg="foo"

memcheck:
	@( valgrind --tool=memcheck --leak-check=full ${BIL_EXE} ${arg} )
