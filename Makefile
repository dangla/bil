#
SHELL = /bin/sh

# Installation folders
PREFIX     := /usr/local
BINDIR     := ${PREFIX}/bin
DATADIR    := ${PREFIX}/share
INCLUDEDIR := ${PREFIX}/include
LIBDIR     := ${PREFIX}/lib


# Bil path
BILPATH_H := BilPath.h
BIL_PATH := ${shell pwd}


# Bil version
BILVERSION_H  := BilVersion.h
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
BILINFO_H     := BilInfo.h
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
BILEXTRALIBS_H  := BilExtraLibs.h


# Do we use extra-libraries?
include make.extralibs
BIL_LIBS += ${BIL_EXTRALIBS}


# Compiler commands
include make.inc


# Options file
BILCONFIG_H := BilConfig.h


#=======================================================================
# Default target executed when no arguments are given to make.
all: bin doc


#=======================================================================
# Target rules for binaries

.PHONY: bin
bin: make.inc path_init version_init info_init config_init lib compile link

compile: make.inc path_init version_init info_init lib
	@echo "\nCompilation"
	mkdir -p ${BIL_LIBDIR}
	for i in ${BIL_DIRS}; do ( ${MAKE} -C $$i ); done

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
	( ${MAKE} -C doc install-doc )


#=======================================================================
# Target rules for cleaning

clean:
	@echo "\nCleaning directories (removing binaries,...)"
	for i in ${BIL_DIRS} doc examples base; do (${MAKE} -C $$i clean); done
	rm -f make.log
	@echo "\nRemoving libraries"
	rm -f ${BIL_LIBDIR}/*.a
	rm -f ${BIL_LIBDIR}/*.so
	find . -name "*~" -printf "\"%p\"\n" | xargs rm -f


clean-all: clean
	rm -f ${BIL_BINDIR}/${BIL_EXE}
	rm -f ${BIL_BINDIR}/*
	rm -f ${BILVERSION_H} 
	rm -f ${BILINFO_H} 
	rm -f ${BILEXTRALIBS_H}
	rm -f ${BILPATH_H}
	rm -f ${BILCONFIG_H}
	( ${MAKE} -C doc clean-all )


#=======================================================================
# Target rules for building path file

path:
	rm -f ${BILPATH_H}
	touch ${BILPATH_H}
	echo "#ifndef BILPATH_H"    >>  ${BILPATH_H}
	echo "#define BILPATH_H"    >>  ${BILPATH_H}
	echo "#define BIL_PATH \"${BIL_PATH}\"" >>  ${BILPATH_H}
	echo "#endif"    >>  ${BILPATH_H}

path_init:
	@if [ ! -r ${BILPATH_H} ]; then ${MAKE} path ; fi


#=======================================================================
# Target rules for building version file

.PHONY: version
version:
	rm -f ${BILVERSION_H}
	touch ${BILVERSION_H}
	echo "#ifndef BILVERSION_H"    >>  ${BILVERSION_H}
	echo "#define BILVERSION_H"    >>  ${BILVERSION_H}
	echo "#define BIL_MAJOR_VERSION ${BIL_MAJOR_VERSION}" >> ${BILVERSION_H}
	echo "#define BIL_MINOR_VERSION ${BIL_MINOR_VERSION}" >> ${BILVERSION_H}
	echo "#define BIL_PATCH_VERSION ${BIL_PATCH_VERSION}" >> ${BILVERSION_H}
	echo "#define BIL_VERSION   \"${BIL_VERSION}\""       >> ${BILVERSION_H}
	echo "#endif"    >>  ${BILVERSION_H}

version_init:
	@if [ ! -r ${BILVERSION_H} ]; then ${MAKE} version ; fi


#=======================================================================
# Target rules for building info file

info:
	rm -f ${BILINFO_H}
	touch ${BILINFO_H}
	echo "#ifndef BILINFO_H"    >>  ${BILINFO_H}
	echo "#define BILINFO_H"    >>  ${BILINFO_H}
	echo "#define BIL_PROGNAME  \"${BIL_PROGNAME}\""  >> ${BILINFO_H}
	echo "#define BIL_COPYRIGHT \"${BIL_COPYRIGHT}\"" >> ${BILINFO_H}
	echo "#define BIL_DATE      \"${BIL_DATE}\""      >> ${BILINFO_H}
	echo "#define BIL_HOST      \"${BIL_HOST}\""      >> ${BILINFO_H}
	echo "#define BIL_PACKAGER  \"${BIL_PACKAGER}\""  >> ${BILINFO_H}
	echo "#define BIL_OS        \"${BIL_OS}\""        >> ${BILINFO_H}
	echo "#define BIL_SHORT_LICENSE \"${BIL_SHORT_LICENSE}\"" >> ${BILINFO_H}
	echo "#define BIL_URL       \"${BIL_URL}\""       >> ${BILINFO_H}
	echo "#define BIL_EMAIL     \"${BIL_EMAIL}\""     >> ${BILINFO_H}
	echo "#endif"    >>  ${BILINFO_H}

info_init:
	@if [ ! -r ${BILINFO_H} ]; then ${MAKE} info ; fi


#=======================================================================
# Target rules for building extra library file

.PHONY: lib
lib:
	rm -f ${BILEXTRALIBS_H}
	touch ${BILEXTRALIBS_H}
	echo "#ifndef BILEXTRALIBS_H"    >>  ${BILEXTRALIBS_H}
	echo "#define BILEXTRALIBS_H"    >>  ${BILEXTRALIBS_H}
	if [ -n "${LAPACK_DIR}" ] ; then \
	  echo "#define LAPACKLIB   ${LAPACK_DIR}"  >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${BLAS_DIR}" ] ; then \
	  echo "#define BLASLIB     ${BLAS_DIR}"    >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${GSLCBLAS_DIR}" ] ; then \
	  echo "#define GSLCBLASLIB     ${GSLCBLAS_DIR}"    >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${SUPERLU_DIR}" ]; then \
	  echo "#define SUPERLULIB  ${SUPERLU_DIR}" >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${SUPERLUMT_DIR}" ]; then \
	  echo "#define SUPERLUMTLIB  ${SUPERLUMT_DIR}" >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${GSL_DIR}" ]; then \
	  echo "#define GSLLIB  ${GSL_DIR}" >> ${BILEXTRALIBS_H} ; \
	fi
	if [ -n "${PETSC_DIR}" ]; then \
	  echo "#define PETSCLIB  ${GSL_DIR}" >> ${BILEXTRALIBS_H} ; \
	fi
	echo "#endif"    >>  ${BILEXTRALIBS_H}


#=======================================================================
# Target rules for config file

config:
	rm -f ${BILCONFIG_H}
	touch ${BILCONFIG_H}
	echo "#ifndef BILCONFIG_H"    >>  ${BILCONFIG_H}
	echo "#define BILCONFIG_H"    >>  ${BILCONFIG_H}
	if [ ${ENABLE_OPENMP} = ON ]; then \
	  echo "#define HAVE_OPENMP" >> ${BILCONFIG_H} ; \
	fi
	if [ ${ENABLE_PTHREAD} = ON ]; then \
	  echo "#define HAVE_PTHREAD" >> ${BILCONFIG_H} ; \
	fi
	echo "#endif"    >>  ${BILCONFIG_H}

config_init:
	@if [ ! -r ${BILCONFIG_H} ]; then ${MAKE} config ; fi



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
	@( cat .githelp )
	

#=======================================================================
# Target rules for tests

test:
	@( echo "base:"; cd base; pwd )
	@( echo "BIL_PATH = ${BIL_PATH}" )
	@( echo "Last part of BIL_PATH = ${notdir ${BIL_PATH}}" )
	@( echo "BIL_COPYRIGHT = ${BIL_COPYRIGHT}" )
	@( echo "BIL_LIBS = ${BIL_LIBS}" )
	@( echo "SUPERLU_DIR = ${SUPERLU_DIR}" )
	@( echo "BLAS_DIR = ${BLAS_DIR}" )
	@( echo "LAPACK_DIR = ${LAPACK_DIR}" )
	@( echo "UEL_DIR = ${UEL_DIR}" )
	@( echo "CFLAGS = ${CFLAGS}" )
	

#=======================================================================
# Target rules for memcheck through valgrind
# To pass the variable "foo" to make use: make memcheck arg="foo"

memcheck:
	@( valgrind --tool=memcheck --leak-check=full bil ${arg} )
