
# Re-build the library file
message("Updating the library file: ${BIL_EXTRALIBS_FILE}")
file(REMOVE ${BIL_EXTRALIBS_FILE})
file(TOUCH  ${BIL_EXTRALIBS_FILE})

file(APPEND ${BIL_EXTRALIBS_FILE} "#ifndef BILLIB_H\n")
file(APPEND ${BIL_EXTRALIBS_FILE} "#define BILLIB_H\n")


# Do you use the libraries: BLAS, LAPACK, SuperLU (YES=ON/NO=OFF)?
if(NOT DEFINED BLAS_USE)
  set(BLAS_USE   OFF CACHE BOOL "Choose BLAS library")
endif()
if(NOT DEFINED LAPACK_USE)
  set(LAPACK_USE OFF CACHE BOOL "Choose LAPACK library")
endif()
if(NOT DEFINED SLU_USE)
  set(SLU_USE    OFF CACHE BOOL "Choose SuperLU library")
endif()
# Do you use the library: UEL (YES=ON/NO=OFF)?
if(NOT DEFINED UEL_USE)
  set(UEL_USE    OFF CACHE BOOL "Choose UEL library")
endif()


# Where are these libraries?
if(${SLU_USE})
  set(SLU_DIR  /usr/lib/x86_64-linux-gnu)
  set(SLU_LIB  ${SLU_DIR}/libsuperlu.so.5)
  set(BIL_EXTRALIBS ${BIL_EXTRALIBS} ${SLU_LIB})
  file(APPEND ${BIL_EXTRALIBS_FILE} "#define SUPERLULIB ${SLU_DIR}\n")
  message("Use the SuperLU library: " ${SLU_LIB})
endif()

if(${BLAS_USE})
  set(BLAS_DIR  /usr/lib)
  set(BLAS_LIB  ${BLAS_DIR}/libblas.so.3)
  set(BIL_EXTRALIBS ${BIL_EXTRALIBS} ${BLAS_LIB})
  file(APPEND ${BIL_EXTRALIBS_FILE} "#define BLASLIB ${BLAS_DIR}\n")
  message("Use the BLAS library: " ${BLAS_LIB})
endif()

if(${LAPACK_USE})
  set(LAPACK_DIR  /usr/lib)
  set(LAPACK_LIB  ${LAPACK_DIR}/liblapack.so.3)
  set(BIL_EXTRALIBS ${BIL_EXTRALIBS}  ${LAPACK_LIB})
  file(APPEND ${BIL_EXTRALIBS_FILE} "#define LAPACKLIB ${LAPACK_DIR}\n")
  message("Use the LAPACK library: " ${LAPACK_LIB})
endif()

if(${UEL_USE})
  set(UEL_DIR /home/dangla/Documents/Softwares/bil/Projects/IoannisStefanou)
  set(UEL_LIB ${UEL_DIR}/libUEL_64.so)
  set(BIL_EXTRALIBS ${BIL_EXTRALIBS} ${UEL_LIB})
  file(APPEND ${BIL_EXTRALIBS_FILE} "#define UELLIB ${UEL_DIR}\n")
  message("Use the UEL library: " ${UEL_LIB})
endif()


file(APPEND ${BIL_EXTRALIBS_FILE} "#endif")


if(${BIL_EXTRALIBS})
  message("Some extra-libraries are used.")
  message("Full paths to these extra-libraries:")
  message("${BIL_EXTRALIBS}")
else()
  message("No extra-libraries are used")
endif()


#welcome:
#       read -p "Do you use BLAS library (Y/N)?:" BLAS_USE
