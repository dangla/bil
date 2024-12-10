# OS
#----
UNAME = Linux


# Compilers (C++11 for variadic macros, C++17 for autodiff)
#----------
CC     := gcc -x c   -std=c11
CPP    := gcc -x c++ -std=c++17
# Force c++ compiler
CC     := ${CPP}
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
IFLAGS-SRC = -iprefix ${BIL_PATH}/src/ -iwithprefix . -iwithprefix Common -iwithprefix DataSet -iwithprefix Help -iwithprefix Libraries -iwithprefix Macros -iwithprefix Main -iwithprefix Models -iwithprefix Modules -iwithprefix Outputs -iwithprefix Solver -iwithprefix Tools -iwithprefix Parser -iwithprefix Parallelization

IFLAGS-SOLVER = -iprefix ${BIL_PATH}/src/Solver/ -iwithprefix MatrixStorageFormats -iwithprefix ResolutionMethods -iwithprefix VectorStorageFormats

# Additional flags for models
IFLAGS-MODELS = -iprefix ${BIL_PATH}/src/Models/ -iwithprefix .  -iwithprefix DataBases  -iwithprefix Libraries  -iprefix ${BIL_PATH}/src/Models/Methods/ -iwithprefix .  -iwithprefix ConstitutiveIntegrator -iwithprefix FiniteElementMethod -iwithprefix FiniteVolumeMethod -iprefix ${BIL_PATH}/src/Models/ConstitutiveLaws/ -iwithprefix . -iwithprefix PlasticityModels -iwithprefix DamageModels

# Additional flags for modules
IFLAGS-MODULES =  -iprefix ${BIL_PATH}/src/Modules/ -iwithprefix .


IFLAGS = ${IFLAGS-SRC} ${IFLAGS-SOLVER} ${IFLAGS-MODELS} ${IFLAGS-MODULES}


# Option flags
#-------------

sinclude ${BIL_PATH}/make.config
OFLAGS :=
#${info ENABLE_OPENMP = ${ENABLE_OPENMP}}
ifdef ENABLE_OPENMP
  ifeq (${ENABLE_OPENMP},ON)
    OFLAGS += -fopenmp
    #${info OFLAGS = ${OFLAGS}}
  endif
endif


# C/CPP compiler flags
#---------------------
CFLAGS = -pedantic ${DFLAGS} ${WFLAGS} ${IFLAGS} -fPIC ${OFLAGS}

# C compiler flags for models
CFLAGS-MODELS = ${CFLAGS} #${IFLAGS-MODELS}

# C compiler flags for modules
CFLAGS-MODULES = ${CFLAGS} #${IFLAGS-MODULES}


# Fortran compiler flags
#-----------------------
FFLAGS = -fPIC


# Optimization flags
#------------------
OPTIM  = -gdwarf-2 -g3     # Debug
#OPTIM  = -O3 -DNDEBUG      # Release
#OPTIM  = -O3 -g -gdwarf-3  # RelWithDebInfo




# Linker
#-------
LINKER = g++


# Linker flags
#-------------
LFLAGS = -lgfortran ${OFLAGS}




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

