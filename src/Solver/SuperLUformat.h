#ifndef SUPERLUFORMAT_H
#define SUPERLUFORMAT_H


/* class-like structure "SuperLUformat_t" and attributes */

/* vacuous declarations and typedef names */
/* struct SuperLUformat_s    ; typedef struct SuperLUformat_s    SuperLUformat_t ; */
#define SuperLUformat_t     SuperMatrix  /* Mimic the SuperMatrix struct defined in SuperLu */


#include "Mesh.h"

#ifdef SLU_DIR
  #define val(x) #x
  #define xval(x) val(x)
  #include xval(SLU_DIR/SRC/dsp_defs.h)
  #undef val
  #undef xval
  
  
  extern SuperLUformat_t* (SuperLUformat_Create)(Mesh_t*) ;
  extern void           (SuperLUformat_Delete)(void*) ;
#endif





/** The getters */
#define SuperLUformat_GetNbOfRows(a)                     ((a)->nrow)
#define SuperLUformat_GetNbOfColumns(a)                  ((a)->ncol)
#define SuperLUformat_GetStorage(a)                      ((a)->Store)
#define SuperLUformat_GetStorageType(a)                  ((a)->Stype)
#define SuperLUformat_GetDataType(a)                     ((a)->Dtype)
#define SuperLUformat_GetMatrixType(a)                   ((a)->Mtype)




/* complete the structure types by using the typedef */

/* Copy of the SuperLUformat struc as defined in SuperLU */
#ifdef SLU_DIR
struct SuperLUformat_s {
  Stype_t Stype;          /* Storage format (SLU_[NC,NR,SC,SR,NCP,DN]) */
  Dtype_t Dtype;          /* Data type (SLU_[S,D,C,Z])*/
  Mtype_t Mtype;          /* Mathematical property of the matrix (SLU_[GE,..]) */
  int    nrow;            /* number of rows */
  int    ncol;            /* number of columns */
  void*  Store;           /* pointer to the actual storage of the matrix */
} ;
#endif


#endif
