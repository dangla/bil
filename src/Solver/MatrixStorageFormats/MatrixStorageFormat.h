#ifndef MATRIXSTORAGEFORMAT_H
#define MATRIXSTORAGEFORMAT_H



enum MatrixStorageFormat_e {        /* format of the matrix to be stored */
  MatrixStorageFormat_LDUSKL,       /* LDU Skyline format */
  MatrixStorageFormat_SKL,          /* Skyline format */
  MatrixStorageFormat_SuperLU,      /* Harwell-Boeing SuperLU format */
  MatrixStorageFormat_CCS,          /* Compressed column storage format */
  MatrixStorageFormat_NC = MatrixStorageFormat_CCS,
  MatrixStorageFormat_CRS,          /* Compressed row storage format */
  MatrixStorageFormat_Coordinate,   /* Coordinate format */
  MatrixStorageFormat_CSR,          /* Compressed sparse row format */
  MatrixStorageFormat_PetscAIJ,     /* PETSc AIJ format (same as CSR) */
  MatrixStorageFormat_NULL
} ;


/* class-like structure "MatrixStorageFormat_t" and attributes */

/* vacuous declarations and typedef names */
typedef enum MatrixStorageFormat_e MatrixStorageFormat_e ;
struct MatrixStorageFormat_s   ; typedef struct MatrixStorageFormat_s  MatrixStorageFormat_t ;



#include "Options.h"

extern MatrixStorageFormat_t* (MatrixStorageFormat_Create)(Options_t*) ;
extern void                   (MatrixStorageFormat_Delete)(void*) ;


/** The getters */
#define MatrixStorageFormat_GetType(MSF)              ((MSF)->type)
#define MatrixStorageFormat_GetOptions(MSF)           ((MSF)->options)


#include "Utils.h"

#define MatrixStorageFormat_Is(MSF,KEY) \
        (MatrixStorageFormat_GetType(MSF) == Utils_CAT(MatrixStorageFormat_,KEY))

#define MatrixStorageFormat_Type(KEY) \
        ((MatrixStorageFormat_e) Utils_CAT(MatrixStorageFormat_,KEY))



/* complete the structure types by using the typedef */

struct MatrixStorageFormat_s {
  MatrixStorageFormat_e type ;
  Options_t* options ;
} ;

#endif
