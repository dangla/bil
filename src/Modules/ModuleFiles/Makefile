#
SHELL = /bin/sh


COMMA := ,
EMPTY :=
SPACE := ${EMPTY} ${EMPTY}

sinclude ../ListOfModules.inc

FILESC    = ${wildcard *.c}
FILESCPP  = ${wildcard *.cpp}
ALLFILES  = ${FILESC} ${FILESCPP}


ifndef SELECTEDMODULES

SELECTEDFILES   = ${shell ls ${ALLFILES}|sort -V}
SELECTEDMODULES = ${basename ${SELECTEDFILES}}

../ListOfModules.inc :
	echo "SELECTEDMODULES = "${SELECTEDMODULES} > ../ListOfModules.inc
	${MAKE}
	
else

SELECTEDFILES  = ${filter ${addsuffix .c,${SELECTEDMODULES}},${ALLFILES}}
SELECTEDFILES += ${filter ${addsuffix .cpp,${SELECTEDMODULES}},${ALLFILES}}

endif

OTHERFILES  = ${filter-out ${SELECTEDFILES},${ALLFILES}}


# Nb of selected modules
NBMODULES    = ${words ${SELECTEDMODULES}}
SELECTEDMODULESm  = ${addsuffix \#\#m,${SELECTEDMODULES}}
SELECTEDMODULENAMES  = ${addsuffix \",${addprefix \",${SELECTEDMODULES}}}


# To be printed in "ListOfModules.h"
XMODULES     = ${subst ${SPACE},${COMMA},${SELECTEDMODULESm}}
MODULENAMES  = ${subst ${SPACE},${COMMA},${SELECTEDMODULENAMES}}

../ListOfModules.h : ../ListOfModules.inc
	echo "#ifndef LISTOFMODULES_H"  > ../ListOfModules.h
	echo "#define LISTOFMODULES_H" >> ../ListOfModules.h
	echo "#define ListOfModules_Nb          "${NBMODULES}   >> ../ListOfModules.h
	echo "#define ListOfModules_Methods(m)  "${XMODULES}    >> ../ListOfModules.h
	echo "#define ListOfModules_Names       "${MODULENAMES} >> ../ListOfModules.h
	echo "#endif" >> ../ListOfModules.h


clean:
	rm -f *.o *~ 

clean-modules:
	rm -f ../ListOfModules.inc


test:
	@echo '----------------------------------------------'
	@echo 'All files C and C++ included here:'
	@echo 'ALLFILES = '${ALLFILES}
	@echo '----------------------------------------------'
	@echo 'C files:'
	@echo 'FILESC = '${FILESC}
	@echo '----------------------------------------------'
	@echo 'C++ files:'
	@echo 'FILESCPP = '${FILESCPP}
	@echo '----------------------------------------------'
	@echo 'Nb of selected modules:'
	@echo 'NBMODULES = '${NBMODULES}
	@echo '----------------------------------------------'
	@echo 'Selected files:'
	@echo 'SELECTEDFILES = '${SELECTEDFILES}
	@echo '----------------------------------------------'
	@echo 'Selected modules:'
	@echo 'SELECTEDMODULES = '${SELECTEDMODULES}
	@echo 'SELECTEDMODULESm = '${SELECTEDMODULESm}
	@echo 'XMODULES = '${XMODULES}
	@echo '----------------------------------------------'
	@echo 'Sorted modules = '${sort ${SELECTEDMODULES}}
	@echo '----------------------------------------------'
	@echo 'Other files:'
	@echo 'OTHERFILES = '${OTHERFILES}
	@echo '----------------------------------------------'
