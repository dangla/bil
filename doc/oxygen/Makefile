#
#
SHELL = /bin/sh

BIL_OXYGEN_FILE = bil-doxy

all: html

.PHONY: html
html: doxy
	doxygen ${BIL_OXYGEN_FILE}

doxy:
	rm -f ${BIL_OXYGEN_FILE}
	cp Doxyfile.def ${BIL_OXYGEN_FILE}
	echo "PROJECT_NAME           = Bil" >> ${BIL_OXYGEN_FILE}
	echo "PROJECT_NUMBER         = ${BIL_VERSION}" >> ${BIL_OXYGEN_FILE}
	echo "OUTPUT_DIRECTORY       = ${BIL_PATH}/doc/oxygen" >> ${BIL_OXYGEN_FILE}
	echo "INPUT                  = ${BIL_PATH}/src" >> ${BIL_OXYGEN_FILE}
	echo "HTML_OUPUT             = ${BIL_DOC}.doxy" >> ${BIL_OXYGEN_FILE}

clean:
	rm -f *~

clean-all: clean
	rm -rf ${BIL_DOC}.doxy
	rm -rf html
