#
#
SHELL = /bin/sh

# Installation directories
MANDIR     = ${PREFIX}/man
INFODIR    = ${PREFIX}/info

# Langage (fr or en)
# LANG = en

# Documentation
#BIL_EDITION_DIR = ${shell ls -d edition-*}
BIL_EDITION_DIR = edition
BIL_DOC = bil

export BIL_DOC

# Subdirectories of doc
DOC_DIRS = ${BIL_EDITION_DIR} oxygen

all:
	cd ${BIL_EDITION_DIR} && rm -f versions.texi
	cd ${BIL_EDITION_DIR} && echo "@set BIL-VERSION   ${BIL_VERSION}" > versions.texi
	for i in ${DOC_DIRS}; do (cd $$i && ${MAKE} all); done
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.ps ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.pdf ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.txt ..
	rm -rf ${BIL_DOC}.html
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.html ..


#=======================================================================
# Target rules for the reference manual in different formats

doc:
	cd ${BIL_EDITION_DIR} && rm -f versions.texi
	cd ${BIL_EDITION_DIR} && echo "@set BIL-VERSION   ${BIL_VERSION}" > versions.texi
	cd ${BIL_EDITION_DIR} && ${MAKE} all
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.ps ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.pdf ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.txt ..
	rm -rf ${BIL_DOC}.html
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.html ..


#=======================================================================
# Target rules for oxygen

.PHONY: oxygen
oxygen:
	for i in ${DOC_DIRS}; do (cd $$i && ${MAKE} all); done
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.ps ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.pdf ..
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.txt ..
	rm -rf ${BIL_DOC}.html
	cd ${BIL_EDITION_DIR} && mv ${BIL_DOC}.html ..

clean:
	for i in ${DOC_DIRS}; do (cd $$i && ${MAKE} clean); done

clean-all: clean
	rm -f ${BIL_DOC}.ps ${BIL_DOC}.pdf ${BIL_DOC}.txt
	rm -rf ${BIL_DOC}.html 
	rm -rf ${BIL_EDITION_DIR}/versions.texi
	for i in ${DOC_DIRS}; do (cd $$i && ${MAKE} clean-all); done

install-doc:
	mkdir -p ${MANDIR}/man1
	cp -f ${BIL_DOC}.1 ${MANDIR}/man1
	#ln -sf ${MANDIR}/man1/${BIL_DOC}.1 ${MANDIR}/man1/bil.1
	mkdir -p ${INFODIR}
	cp ${BIL_EDITION_DIR}/${BIL_DOC}.info ${INFODIR}
	#ln -sf ${INFODIR}/${BIL_DOC}.info ${INFODIR}/bil.info
	cd ${BIL_EDITION_DIR} && ginstall-info --info-file=${BIL_DOC}.info --info-dir=${INFODIR}
