#
#
SHELL = /bin/sh

include ../make.inc

LIBa     = ${BIL_LIBDIR}/libhsl.a

SRC = $(wildcard *.f)

OBJ = ${SRC:.f=.o}

#.SUFFIXES: .o .f


# Rules
# -----
entry: ${LIBa}


${LIBa}: ${OBJ}
	@echo "\nCreating the static library "${LIBa}
	${AR} ${LIBa} ${OBJ}
	${RANLIB} ${LIBa}

#.f.o:
#	${F77} ${OPTIM} ${FFLAGS} -c $<


%.o: %.f ../make.inc
	${F77} ${OPTIM} ${FFLAGS} -c $<


clean:
	rm -f *.o
