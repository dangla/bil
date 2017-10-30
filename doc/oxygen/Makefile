#
#
SHELL = /bin/sh

all: html

.PHONY: html
html: doxy
	doxygen Doxyfile

doxy:
	rm -f Doxyfile
	cp Doxyfile.def Doxyfile
	echo "PROJECT_NAME           = Bil" >> Doxyfile
	echo "PROJECT_NUMBER         = ${BIL_VERSION}" >> Doxyfile
	echo "OUTPUT_DIRECTORY       = ${BIL_PATH}/doc/oxygen" >> Doxyfile
	echo "INPUT                  = ${BIL_PATH}" >> Doxyfile
	echo "EXCLUDE                = ${BIL_PATH}/base \\ " >> Doxyfile
	echo "                         ${BIL_PATH}/examples \\ " >> Doxyfile
	echo "                         ${BIL_PATH}/doc" >> Doxyfile

clean:
	rm -f *~

clean-all: clean
	rm -rf html
