#
#
SHELL = /bin/sh

DIRSEXAMPLES = $(shell ls -d */)

clean:
	for i in ${DIRSEXAMPLES}; do (cd $$i \
	&& rm -f *~ *.t* *.p* *.sto *.graph* *.pos*\
	); done


EMPTY :=
SPACE := ${EMPTY} ${EMPTY}
# The two 'PRINT' methods work
#PRINT = echo -e ${addsuffix '\n',$(1)}
PRINT = echo -e ${foreach a,$(1),${a}'\n'}


listdirectories:
	@echo '----------------------------------------------'
	@echo 'Directories'
	@echo '----------------------------------------------'
	@echo '1. Example directories'
	@${call PRINT,${DIRSEXAMPLES}}
	@echo '----------------------------------------------'

