#
SHELL = /bin/sh


COMMA := ,
EMPTY :=
SPACE := ${EMPTY} ${EMPTY}

sinclude ../ListOfModels.inc

FILESC    = ${wildcard *.c}
FILESCPP  = ${wildcard *.cpp}
ALLFILES  = ${FILESC} ${FILESCPP}


ifndef SELECTEDMODELS

SELECTEDFILES  = ${shell ls ${ALLFILES}|sort -V}
SELECTEDMODELS = ${basename ${SELECTEDFILES}}

../ListOfModels.inc:
	echo "SELECTEDMODELS = "${SELECTEDMODELS} > ../ListOfModels.inc
	${MAKE}

else

SELECTEDFILES  = ${filter ${addsuffix .c,${SELECTEDMODELS}},${ALLFILES}}
SELECTEDFILES += ${filter ${addsuffix .cpp,${SELECTEDMODELS}},${ALLFILES}}

endif

OTHERFILES  = ${filter-out ${SELECTEDFILES},${ALLFILES}}


# Nb of selected models
NBMODELS    = ${words ${SELECTEDMODELS}}
SELECTEDMODELSm  = ${addsuffix \#\#m,${SELECTEDMODELS}}
SELECTEDMODELNAMES  = ${addsuffix \",${addprefix \",${SELECTEDMODELS}}}


# To be printed in "ListOfModels.h"
XMODELS     = ${subst ${SPACE},${COMMA},${SELECTEDMODELSm}}
MODELNAMES  = ${subst ${SPACE},${COMMA},${SELECTEDMODELNAMES}}

../ListOfModels.h: ../ListOfModels.inc
	echo "#ifndef LISTOFMODELS_H"  > ../ListOfModels.h
	echo "#define LISTOFMODELS_H" >> ../ListOfModels.h
	echo "#define ListOfModels_NbOfModels  "${NBMODELS}   >> ../ListOfModels.h
	echo "#define ListOfModels_Methods(m)  "${XMODELS}    >> ../ListOfModels.h
	echo "#define ListOfModels_ModelNames  "${MODELNAMES} >> ../ListOfModels.h
	echo "#endif" >> ../ListOfModels.h


clean:
	rm -f *.o *~ 


clean-models:
	rm -f ../ListOfModels.inc


replace:
	for file in ${SELECTEDFILES}; do (\
	sed 's/^void ex/int ex/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_1:
	for file in ${SELECTEDFILES}; do (\
	sed 's/^#include "defs/#include "..\/defs/g' $$file > $$file.bak && \
	mv $$file.bak $$file &&\
	sed 's/^#include "lib/#include "..\/lib/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_2:
	for file in ${SELECTEDFILES}; do (\
	sed 's/^  printf(TITLE"\\n")/  printf(TITLE)/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_3:
	for file in ${ALLFILES}; do (\
	sed '/^void tb/ {s/,geom_t geom,int \*nv//}' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_4:
	for file in ${ALLFILES}; do (\
	sed '/chsurf/ {s/chsurf(x,r,dim,geom,dt,t,cg,el)/chsurf(x,r,dim,geom,dt,t,cg,el,el.fi)/}' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_5:
	for file in ${ALLFILES}; do (\
	sed 's/if(geom == 1)/if(geom == AXIS)/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_6:
	for file in ${ALLFILES}; do (\
	sed 's/define MODELINDEX/define MODELNAME/' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_7:
	for file in ${ALLFILES}; do (\
	sed 's/^#include "..\/defs.h"/#include "model.h"/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_8:
	for file in ${ALLFILES}; do (\
	sed '/^#include "..\/lib.h"/d' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_9:
	for file in ${ALLFILES}; do (\
	sed 's/mat->neq/Material_GetNbOfEquations(mat)/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_10:
	for file in ${ALLFILES}; do (\
	sed 's/^#include "Common.h/#include "CommonModel.h/g' $$file > $$file.bak && \
	mv $$file.bak $$file \
	); done


work_11:
	for file in ${ALLFILES}; do (\
	cp $$file $$file.bak && \
	sed 's/^#include "PredefinedMethods.h/#include "PredefinedModelMethods.h/g' $$file.bak > $$file \
	); done


work_test:
	for file in M52b.c; do (\
	cp $$file $$file.bak && \
	sed 's/^#include "PredefinedMethods.h/#include "PredefinedModelMethods.h/g' $$file.bak > $$file \
	); done


printfilenames:
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
	@echo 'Nb of selected models:'
	@echo 'NBMODELS = '${NBMODELS}
	@echo '----------------------------------------------'
	@echo 'Selected files:'
	@echo 'SELECTEDFILES = '${SELECTEDFILES}
	@echo '----------------------------------------------'
	@echo 'Selected models:'
	@echo 'SELECTEDMODELS = '${SELECTEDMODELS}
	@echo 'SELECTEDMODELSm = '${SELECTEDMODELSm}
	@echo 'XMODELS = '${XMODELS}
	@echo '----------------------------------------------'
	@echo 'Sorted models = '${sort ${SELECTEDMODELS}}
	@echo '----------------------------------------------'
	@echo 'Other files:'
	@echo 'OTHERFILES = '${OTHERFILES}
	@echo '----------------------------------------------'


showtitles:
	grep --include="*.c" --include="*.cpp" -T --color -e "#define TITLE" *.c*
