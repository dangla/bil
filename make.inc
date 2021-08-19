# OS
#----
UNAME = Linux


# Compilers (C++11 for variadic macros)
#----------
CC     := gcc -x c   -std=c11
CPP    := gcc -x c++ -std=c++11
F77    := gfortran -cpp
LEX    := flex #-d
YACC   := bison #-t


# Compiler and preprocessor flags (recommended by GSL)
# see www.gnu.org/software/gsl/manual/html_node/GCC-warning-options-for-numerical-programs.html)


# Warning flags
#--------------
WFLAGS := -Wall -Wcast-align -Wcast-qual -Wnested-externs -Wpointer-arith -Wmissing-prototypes -Wstrict-prototypes -Wshadow -fno-common -Dinline=  -Wvariadic-macros #-W #-Wtraditional #-Wconversion #-Wwrite-strings #-Werror #-fshort-enums


# Directory flags (-I replaced by -iquote)
#----------------
# DFLAGS = -I${BIL_PATH}
DFLAGS = -iquote${BIL_PATH} -D '__BASEFILE__=${basename ${notdir $@}}'


# Preprocessor flags
#-------------------
IFLAGS-SRC = -iprefix ${BIL_PATH}/src/ -iwithprefix . -iwithprefix Common -iwithprefix DataSet -iwithprefix Help -iwithprefix Libraries -iwithprefix Macros -iwithprefix Main -iwithprefix Models -iwithprefix Modules -iwithprefix Outputs -iwithprefix Solver -iwithprefix Tools -iwithprefix Parser

IFLAGS-SOLVER = -iprefix ${BIL_PATH}/src/Solver/ -iwithprefix MatrixStorageFormats -iwithprefix ResolutionMethods

# Additional flags for models
IFLAGS-MODELS = -iprefix ${BIL_PATH}/src/Models/ -iwithprefix .  -iwithprefix DataBases  -iwithprefix Libraries  -iwithprefix Methods  -iwithprefix ConstitutiveLaws

# Additional flags for modules
IFLAGS-MODULES =  -iprefix ${BIL_PATH}/src/Modules/ -iwithprefix .


IFLAGS = ${IFLAGS-SRC} ${IFLAGS-SOLVER} ${IFLAGS-MODELS} ${IFLAGS-MODULES}


# C/CPP compiler flags
#---------------------
CFLAGS = -pedantic ${DFLAGS} ${WFLAGS} ${IFLAGS} -fPIC

# C compiler flags for models
CFLAGS-MODELS = ${CFLAGS} #${IFLAGS-MODELS}

# C compiler flags for modules
CFLAGS-MODULES = ${CFLAGS} #${IFLAGS-MODULES}


# Fortran compiler flags
#-----------------------
FFLAGS = -fPIC


# Optimization flags
#------------------
OPTIM  = -gdwarf-2 -g3 #-O2




# Linker
#-------
LINKER = g++


# Linker flags
#-------------
LFLAGS = -lgfortran




# Libraries
#----------

# Library extensions. Depends on which kind of library you want to install
# Choose either "so" (shared object) or "a" (static)
LIBBILEXT   = so

# How you create a static library on this machine
AR     = /usr/bin/ar ruvs
RANLIB = ranlib


# To create and install a shared library
CCSO = gcc -shared
LDCONFIG  = /sbin/ldconfig 

