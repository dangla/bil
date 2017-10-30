#
#
SHELL = /bin/sh

# Edition de Bil
BIL_EDITION_FILE = edition.texi
BIL_EDITION_PATH = ${shell pwd}
BIL_EDITION       = ${word 2,${subst /edition-, ,${BIL_EDITION_PATH}}}
BIL_MAJOR_EDITION = ${word 1,${subst ., ,${BIL_EDITION}}}
BIL_MINOR_EDITION = ${word 2,${subst ., ,${BIL_EDITION}}}


dvi: edition_init
	texi2dvi ${BIL_DOC}.texi

ps: dvi
	dvips ${BIL_DOC} -o

pdf: edition_init
	texi2pdf ${BIL_DOC}.texi 

html: edition_init
	perl texi2html -init_file ${BIL_DOC}.t2h ${BIL_DOC}.texi
	rm -f index.html
	ln -s ${BIL_DOC}.html index.html

html1: edition_init
	perl texi2html -prefix ${BIL_DOC}-full -init_file ${BIL_DOC}-full.t2h ${BIL_DOC}.texi

html2: edition_init
	makeinfo --html -o ${BIL_DOC}.html ${BIL_DOC}.texi

info: edition_init
	makeinfo --enable-encoding -o ${BIL_DOC}.info ${BIL_DOC}.texi

install-info:
	mkdir -p /usr/local/info
	cp ${BIL_DOC}.info* /usr/local/info/
	install-info --info-dir=/usr/local/info --info-file=bil.info

infoz: info
	tar zcvf ${BIL_DOC}-info.tgz ${BIL_DOC}.info*

txt: edition_init
	makeinfo --no-headers -o ${BIL_DOC}.txt ${BIL_DOC}.texi

all: ps pdf html2 info txt

clean:
	rm -f *~ *.cp* *.fn* *.ky* *.pg* *.tp* *.vr* *.mv*\
          *.log *.toc *.aux *.dvi *.html *.txt\
          ${BIL_DOC}.ps ${BIL_DOC}.pdf ${BIL_DOC}.info* ${BIL_DOC}-info.tgz
	rm -rf ${BIL_DOC}.html

clean-all: clean
	rm -rf ${BIL_EDITION_FILE}

edition:
	rm -f ${BIL_EDITION_FILE}
	echo "@set MAJOR_EDITION ${BIL_MAJOR_EDITION}" >  ${BIL_EDITION_FILE}
	echo "@set MINOR_EDITION ${BIL_MINOR_EDITION}" >> ${BIL_EDITION_FILE}
	echo "@set EDITION   \"${BIL_EDITION}\""       >> ${BIL_EDITION_FILE}

edition_init:
	@if [ ! -r ${BIL_EDITION_FILE} ]; then ${MAKE} edition ; fi
