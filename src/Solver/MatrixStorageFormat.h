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
  MatrixStorageFormat_NULL
} ;


/* class-like structure "MatrixStorageFormat_t" and attributes */

/* vacuous declarations and typedef names */
typedef enum MatrixStorageFormat_e MatrixStorageFormat_t ;



#include "Utils.h"

#define MatrixStorageFormat_Is(MSF,KEY) \
        (MSF == Utils_CAT(MatrixStorageFormat_,KEY))
        
        
#endif
