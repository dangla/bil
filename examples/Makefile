#
#
SHELL = /bin/sh

DIRSEXAMPLES = $(wildcard ex*)

clean:
	for i in ${DIRSEXAMPLES}; do (cd $$i \
	&& rm -f *~ *.t* *.p* *.sto *.graph* *.pos*\
	); done
