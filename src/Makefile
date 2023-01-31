#
#
SHELL = /bin/sh

include ../make.inc
#${info in src/Makefile: CFLAGS = ${CFLAGS}}

LIBBILa     = ${BIL_LIBDIR}/libbil-${BIL_VERSION}.a

LIBBILso    = ${BIL_LIBDIR}/libbil-${BIL_VERSION}.so
LIBBILsoname= libbil-${BIL_VERSION}.so

# Library to be installed
LIBBIL      = ${BIL_LIBDIR}/libbil-${BIL_VERSION}.${LIBBILEXT}

sinclude Models/ListOfModels.inc
sinclude Modules/ListOfModules.inc


ifndef SELECTEDMODELS
step1: models
	@echo "Re-execution"
	${MAKE}
endif

ifndef SELECTEDMODULES
step2: modules
	@echo "Re-execution"
	${MAKE}
endif


# Directories
# -----------
DIRMODELFILES  = Models/ModelFiles/
DIRMODULEFILES = Modules/ModuleFiles/
DIRSALL = ${shell ls -d */ ./ */*/}
DIRSOTHER := ${filter-out ${DIRMODELFILES} ${DIRMODULEFILES},${DIRSALL}}


COMMA = ,


# Sources
# -------
# Main source
SRCMAIN = Main/Main.c

# Model sources
SELECTEDMODELS1 = ${addprefix ${DIRMODELFILES},${SELECTEDMODELS}}
ALLMODELFILES   = ${wildcard ${addsuffix *.c,${DIRMODELFILES}}}
ALLMODELFILES  += ${wildcard ${addsuffix *.cpp,${DIRMODELFILES}}}
SRCMODELSC   = ${filter ${addsuffix .c,${SELECTEDMODELS1}},${ALLMODELFILES}}
SRCMODELSCPP = ${filter ${addsuffix .cpp,${SELECTEDMODELS1}},${ALLMODELFILES}}
SRCMODELS    = ${SRCMODELSC} ${SRCMODELSCPP}

# Module sources
SELECTEDMODULES1 = ${addprefix ${DIRMODULEFILES},${SELECTEDMODULES}}
ALLMODULEFILES   = ${wildcard ${addsuffix *.c,${DIRMODULEFILES}}}
ALLMODULEFILES  += ${wildcard ${addsuffix *.cpp,${DIRMODULEFILES}}}
SRCMODULESC   = ${filter ${addsuffix .c,${SELECTEDMODULES1}},${ALLMODULEFILES}}
SRCMODULESCPP = ${filter ${addsuffix .cpp,${SELECTEDMODULES1}},${ALLMODULEFILES}}
SRCMODULES    = ${SRCMODULESC} ${SRCMODULESCPP}

# Other sources
SRCOTHERC   = ${wildcard ${addsuffix *.c,${DIRSOTHER}}}
SRCOTHERCPP = ${wildcard ${addsuffix *.cpp,${DIRSOTHER}}}
SRCOTHERF   = ${wildcard ${addsuffix *.f,${DIRSOTHER}}}
SRCOTHER    = ${SRCOTHERC} ${SRCOTHERCPP} ${SRCOTHERF}

# All sources
SRCALL = ${SRCOTHER} ${SRCMODELS} ${SRCMODULES}
# All sources except Main
SRCALL1 = ${filter-out ${SRCMAIN},${SRCALL}}


# Headers
# -------
HEADERSC   = ${wildcard ${addsuffix *.h,${DIRSALL}}}
HEADERSCPP = ${wildcard ${addsuffix *.hpp,${DIRSALL}}}
HEADERSALL = ${HEADERSC} ${HEADERSCPP}
HEADERSONLY = ${filter-out ${addsuffix .h,${basename ${SRCALL}}},${HEADERSALL}}


# Binaries
# --------
# Main binary
OBJMAIN    = Main/Main.o
OBJMODELS  = ${filter %.o,${SRCMODELS:.c=.o} ${SRCMODELS:.cpp=.o}}
OBJMODULES = ${filter %.o,${SRCMODULES:.c=.o} ${SRCMODULES:.cpp=.o}}
OBJOTHER   = ${filter %.o,${SRCOTHER:.c=.o} ${SRCOTHER:.cpp=.o} ${SRCOTHER:.f=.o}}
# All binaries
OBJALL     = ${OBJOTHER} ${OBJMODELS} ${OBJMODULES}
# All binaries except Main
OBJALL1 = ${filter-out ${OBJMAIN},${OBJALL}}



# Rules
# -----


#=======================================================================
# Default target executed when no arguments are given to make.
entry: parser models modules ${LIBBIL}



#=======================================================================
# Target rules for creating libraries

${LIBBILa}: ${SRCALL} ${OBJALL}
	@echo "\nCreating the static library "${LIBBILa}
	${AR} ${LIBBILa} ${OBJALL1}
	${RANLIB} ${LIBBILa}


${LIBBILso}: ${SRCALL} ${OBJALL}
	@echo "\nCreating the shared library "${LIBBILso}
	${CCSO} ${OPTIM} -Wl,-soname=${LIBBILsoname} -o ${LIBBILso} ${OBJALL1}
	


#=======================================================================
# We define pattern rules for sources older than headers

%.c %.cpp:: %.h
	@echo "\nTouch the source file: "$@
	touch $@

%.c %.cpp:: %.hpp
	@echo "\nTouch the source file: "$@
	touch $@



#=======================================================================
# We define pattern rules for compilation of objects in OBJALL

${DIRMODELFILES}%.o: ${DIRMODELFILES}%.c
	@echo "\nCompilation of the C model: "$<
	${CC} ${OPTIM} ${CFLAGS-MODELS} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODELFILES}%.o: ${DIRMODELFILES}%.cpp
	@echo "\nCompilation of the C++ model: "$<
	${CPP} ${OPTIM} ${CFLAGS-MODELS} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODULEFILES}%.o: ${DIRMODULEFILES}%.c
	@echo "\nCompilation of the C module: "$<
	${CC} ${OPTIM} ${CFLAGS-MODULES} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODULEFILES}%.o: ${DIRMODULEFILES}%.cpp
	@echo "\nCompilation of the C++ module: "$<
	${CPP} ${OPTIM} ${CFLAGS-MODULES} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

%.o: %.c
	@echo "\nCompilation of the C source: "$<
	${CC} ${OPTIM} ${CFLAGS} -c $< -o $@

%.o: %.cpp
	@echo "\nCompilation of the C++ source: "$<
	${CPP} ${OPTIM} ${CFLAGS} -c $< -o $@
	
%.o: %.f
	@echo "\nCompilation of the fortran source: "$<
	${F77} ${OPTIM} ${FFLAGS} -c $< -o $@



#=======================================================================
# Extra rules with no recipes (the recipes are given above).
# Force some objects to be recompiled if some prerequisites change.

${OBJALL}: ${HEADERSONLY} ../make.inc

Models/Models.o: Models/ListOfModels.inc

Modules/Modules.o: Modules/ListOfModules.inc



#=======================================================================
# Target rules for compiling models

.PHONY: models
models:
	@echo "\nLook for available models"
	( ${MAKE} -C ${DIRMODELFILES} )



#=======================================================================
# Target rules for compiling modules

.PHONY: modules
modules:
	@echo "\nLook for available modules"
	( ${MAKE} -C ${DIRMODULEFILES} )



#=======================================================================
# Target rules for cleaning

clean:
	rm -f ${wildcard ${addsuffix *.o,${DIRSOTHER}}}
	rm -f ${wildcard ${addsuffix *~,${DIRSOTHER}}}
	( ${MAKE} -C ${DIRMODELFILES} clean )
	( ${MAKE} -C ${DIRMODULEFILES} clean )



#=======================================================================
# Target rules for building the parser

parser:
	@echo "\nBuilding the parser"
	( ${MAKE} -C Parser )


EMPTY :=
SPACE := ${EMPTY} ${EMPTY}
# The two 'PRINT' methods work
#PRINT = echo -e ${addsuffix '\n',$(1)}
PRINT = echo -e ${foreach a,$(1),${a}'\n'}



#=======================================================================
# Target rules for listing sources

listsources:
	@echo '----------------------------------------------'
	@echo 'Sources'
	@echo '----------------------------------------------'
	@${call PRINT,${SRCALL}}

listsourcesdetailed:
	@echo '----------------------------------------------'
	@echo 'Sources'
	@echo '----------------------------------------------'
	@echo '1. All source files (SRCALL)'
	@${call PRINT,${SRCALL}}
	@echo '----------------------------------------------'
	@echo '2. Source Model Files (SRCMODELS)'
	@${call PRINT,${SRCMODELS}}
	@echo '----------------------------------------------'
	@echo '3. Source Module Files (SRCMODULES)'
	@${call PRINT,${SRCMODULES}}
	@echo '----------------------------------------------'
	@echo '4. Other source files .c (SRCOTHERC)'
	@${call PRINT,${SRCOTHERC}}
	@echo '----------------------------------------------'
	@echo '5. Other source files .cpp (SRCOTHERCPP)'
	@${call PRINT,${SRCOTHERCPP}}
	@echo '----------------------------------------------'




#=======================================================================
# Target rules for listing headers

listheaders:
	@echo '----------------------------------------------'
	@echo 'Headers'
	@echo '----------------------------------------------'
	@${call PRINT,${HEADERSALL}}

listheadersdetailed:
	@echo '----------------------------------------------'
	@echo 'Headers'
	@echo '----------------------------------------------'
	@echo '1. All Headers (HEADERSALL)'
	@${call PRINT,${HEADERSALL}}
	@echo '----------------------------------------------'
	@echo '2. Headers from sources .h (HEADERSC)'
	@${call PRINT,${HEADERSC}}
	@echo '----------------------------------------------'
	@echo '3. Headers from sources .hpp (HEADERSCPP)'
	@${call PRINT,${HEADERSCPP}}
	@echo '----------------------------------------------'
	@echo '4. Header files without associated source files (HEADERSONLY)'
	@${call PRINT,${HEADERSONLY}}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'




#=======================================================================
# Target rules for listing object files

listobjects:
	@echo '----------------------------------------------'
	@echo 'Objects (OBJALL)'
	@echo '----------------------------------------------'
	@${call PRINT,${OBJALL}}
	@echo '----------------------------------------------'




#=======================================================================
# Target rules for listing directories

listdirectories:
	@echo '----------------------------------------------'
	@echo 'Directories'
	@echo '----------------------------------------------'
	@echo '1. All directories (DIRSALL)'
	@${call PRINT,${DIRSALL}}
	@echo '----------------------------------------------'
	@echo '2. Directories except ModelFiles and ModuleFiles (DIRSOTHER)'
	@${call PRINT,${DIRSOTHER}}
	@echo '----------------------------------------------'




#=======================================================================
# Target rules for some operations

rename:
	for file in ${SRCALL} ${HEADERSALL}; do (\
	sed 's/^#include "..\/Tools\/Math.h/#include "Math_.h/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done




#=======================================================================
# Tests

DIRS = ${shell ls -d */ */*/}
LISTDIRS = ${shell ls -d ${1}/}
LISTN = ${foreach a,${1},${call LISTN,${a}}}

test:
	@echo '----------------------------------------------'
	@echo 'Directories'
	@echo '----------------------------------------------'
	@echo '1. All directories (DIRS)'
	@echo ${DIRS}
	@echo '----------------------------------------------'
	@${call PRINT,${DIRS}}



test1:
	@echo '----------------------------------------------'
	@echo 'Source files .cpp (SRCALL)'
	@${call PRINT,${wildcard ${addsuffix *.c,${DIRSOTHER}}} ${wildcard ${addsuffix *.cpp,${DIRSOTHER}}}}
	@echo '----------------------------------------------'
	@echo 'Object files (OBJCPP from SRCOTHERC)'
	@${call PRINT,${SRCOTHERC:%.cpp=%.o}}
	@${call PRINT,${filter %.o,${patsubst %.c,%.o,toto.cpp titi.c}}}

