#
#
SHELL = /bin/sh

include ../make.inc

LIBBILa     = ${BIL_LIBDIR}/libbil-${BIL_VERSION}.a

LIBBILso    = ${BIL_LIBDIR}/libbil-${BIL_VERSION}.so

#Library to be installed
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
SRCOTHER    = ${SRCOTHERC} ${SRCOTHERCPP}

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
OBJOTHER   = ${filter %.o,${SRCOTHER:.c=.o} ${SRCOTHER:.cpp=.o}}
# All binaries
OBJALL     = ${OBJOTHER} ${OBJMODELS} ${OBJMODULES}
# All binaries except Main
OBJALL1 = ${filter-out ${OBJMAIN},${OBJALL}}


# Rules
# -----
entry: models modules ${LIBBIL}


${LIBBILa}: ${SRCALL} ${OBJALL}
	@echo "\nCreating the static library "${LIBBILa}
	${AR} ${LIBBILa} ${OBJALL1}
	${RANLIB} ${LIBBILa}


${LIBBILso}: ${SRCALL} ${OBJALL}
	@echo "\nCreating the shared library "${LIBBILso}
	${CCSO} ${OPTIM} -o ${LIBBILso} ${OBJALL1}
	

# We define pattern rules for sources older than headers

%.c %.cpp:: %.h
	@echo "\nTouch the source file: "$@
	touch $@

%.c %.cpp:: %.hpp
	@echo "\nTouch the source file: "$@
	touch $@


# We define pattern rules for compilation of objects in OBJALL

${DIRMODELFILES}%.o: ${DIRMODELFILES}%.c
	@echo "\nCompilation of the C model: "$<
	${CC} ${OPTIM} ${CFLAGS-MODELS} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODELFILES}%.o: ${DIRMODELFILES}%.cpp
	@echo "\nCompilation of the C++ model: "$<
	${CC} ${OPTIM} ${CFLAGS-MODELS} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODULEFILES}%.o: ${DIRMODULEFILES}%.c
	@echo "\nCompilation of the C module: "$<
	${CC} ${OPTIM} ${CFLAGS-MODULES} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

${DIRMODULEFILES}%.o: ${DIRMODULEFILES}%.cpp
	@echo "\nCompilation of the C++ module: "$<
	${CC} ${OPTIM} ${CFLAGS-MODULES} -D 'BASENAME=${basename ${notdir $@}}' -c $< -o $@

%.o: %.c
	@echo "\nCompilation of the C source: "$<
	${CC} ${OPTIM} ${CFLAGS} -c $< -o $@

%.o: %.cpp
	@echo "\nCompilation of the C++ source: "$<
	${CC} ${OPTIM} ${CFLAGS} -c $< -o $@


# Extra rules with no recipes (the recipes are given above).
# Force some objects to be recompiled if some prerequisites change.
${OBJALL}: ${HEADERSONLY} ../make.inc

Models/Models.o: Models/ListOfModels.inc

Modules/Modules.o: Modules/ListOfModules.inc




.PHONY: models
models:
	@echo "\nLook for available models"
	cd ${DIRMODELFILES} && ${MAKE}

.PHONY: modules
modules:
	@echo "\nLook for available modules"
	cd ${DIRMODULEFILES} && ${MAKE}

clean:
	rm -f ${wildcard ${addsuffix *.o,${DIRSOTHER}}}
	rm -f ${wildcard ${addsuffix *~,${DIRSOTHER}}}
	cd ${DIRMODELFILES} && ${MAKE} clean
	cd ${DIRMODULEFILES} && ${MAKE} clean


test:
	@echo '----------------------------------------------'
	@echo 'Directories'
	@echo '----------------------------------------------'
	@echo '1. All directories'
	@echo 'DIRSALL = '${DIRSALL}
	@echo '----------------------------------------------'
	@echo '2. Directories except those of Models and Modules'
	@echo 'DIRSOTHER = '${DIRSOTHER}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Sources'
	@echo '----------------------------------------------'
	@echo '1. All source files'
	@echo 'SRCALL = '${SRCALL}
	@echo '----------------------------------------------'
	@echo '2. Source Model Files'
	@echo 'SRCMODELS = '${SRCMODELS}
	@echo '----------------------------------------------'
	@echo '3. Source Module Files'
	@echo 'SRCMODULES = '${SRCMODULES}
	@echo '----------------------------------------------'
	@echo '4. Other source files .c'
	@echo 'SRCOTHERC = '${SRCOTHERC}
	@echo '----------------------------------------------'
	@echo '5. Other source files .cpp'
	@echo 'SRCOTHERCPP = '${SRCOTHERCPP}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Headers'
	@echo '----------------------------------------------'
	@echo '1. All Headers'
	@echo 'HEADERSALL = '${HEADERSALL}
	@echo '----------------------------------------------'
	@echo '2. Headers from sources .h'
	@echo 'HEADERSC = '${HEADERSC}
	@echo '----------------------------------------------'
	@echo '3. Headers from sources .hpp'
	@echo 'HEADERSCPP = '${HEADERSCPP}
	@echo '----------------------------------------------'
	@echo '4. Header files without associated source files'
	@echo 'HEADERSONLY = '${HEADERSONLY}
	@echo '----------------------------------------------'
	@echo '----------------------------------------------'
	@echo 'Objects'
	@echo '----------------------------------------------'
	@echo 'OBJALL = '${OBJALL}
	@echo '----------------------------------------------'

test1:
	@echo '----------------------------------------------'
	@echo 'Source files .cpp'
	@echo 'SRCALL = '${wildcard ${addsuffix *.c,${DIRSOTHER}}} ${wildcard ${addsuffix *.cpp,${DIRSOTHER}}}
	@echo '----------------------------------------------'
	@echo 'Object files'
	@echo 'OBJCPP from SRCOTHERC = '${SRCOTHERC:%.cpp=%.o}
	@echo 'OBJCPP from SRCOTHERC = '${filter %.o,${patsubst %.c,%.o,toto.cpp titi.c}}

