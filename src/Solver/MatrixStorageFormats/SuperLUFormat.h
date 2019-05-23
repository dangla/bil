#ifndef SUPERLUFORMAT_H
#define SUPERLUFORMAT_H


/* class-like structure "SuperLUFormat_t" and attributes */

/* vacuous declarations and typedef names */
//struct SuperLUFormat_s    ; typedef struct SuperLUFormat_s    SuperLUFormat_t ;
#define SuperLUFormat_t     SuperMatrix  /* Mimic the SuperMatrix struct defined in SuperLu */



#include "Mesh.h"

#include "superlu.h"

  


extern SuperLUFormat_t* (SuperLUFormat_Create)(Mesh_t*) ;
extern void             (SuperLUFormat_Delete)(void*) ;







/** The getters */
#define SuperLUFormat_GetNbOfRows(a)                     ((a)->nrow)
#define SuperLUFormat_GetNbOfColumns(a)                  ((a)->ncol)
#define SuperLUFormat_GetStorage(a)                      ((a)->Store)
#define SuperLUFormat_GetStorageType(a)                  ((a)->Stype)
#define SuperLUFormat_GetDataType(a)                     ((a)->Dtype)
#define SuperLUFormat_GetMatrixType(a)                   ((a)->Mtype)




/* complete the structure types by using the typedef */

#if 0
/* Copy of the SuperMatrix struc as defined in SuperLU */
enum Stype_e {
    SLU_NC,    /* column-wise, no supernode */
    SLU_NR,    /* row-wize, no supernode */
    SLU_SC,    /* column-wise, supernode */
    SLU_SR,    /* row-wise, supernode */
    SLU_NCP,   /* column-wise, column-permuted, no supernode 
                  (The consecutive columns of nonzeros, after permutation,
                   may not be stored  contiguously.) */
    SLU_DN     /* Fortran style column-wise storage for dense matrix */
} ;

typedef enum Stype_e Stype_t ;

enum Dtype_e {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} ;

typedef enum Dtype_e Dtype_t ;

enum Mtype_e {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} ;

typedef enum Mtype_e Mtype_t ;
#endif


struct SuperLUFormat_s {
  Stype_t Stype;          /* Storage format (SLU_[NC,NR,SC,SR,NCP,DN]) */
  Dtype_t Dtype;          /* Data type (SLU_[S,D,C,Z])*/
  Mtype_t Mtype;          /* Mathematical property of the matrix (SLU_[GE,..]) */
  int    nrow;            /* number of rows */
  int    ncol;            /* number of columns */
  void*  Store;           /* pointer to the actual storage of the matrix */
} ;



#endif
