#ifndef MATRIX_H
#define MATRIX_H

/* class-like structures "Matrix_t" and attributes */

/* vacuous declarations and typedef names */
struct Matrix_s       ; typedef struct Matrix_s       Matrix_t ;


#include "Mesh.h"
#include "Options.h"
#include "Element.h"
 
//extern Matrix_t*   (Matrix_Create)                (Mesh_t*,Options_t*) ;
extern Matrix_t*   (Matrix_Create)                (Mesh_t*,Options_t*,const int) ;
extern void        (Matrix_Delete)                (void*) ;
extern void        (Matrix_AssembleElementMatrix) (Matrix_t*,Element_t*,double*) ;
extern void        (Matrix_PrintMatrix)           (Matrix_t*,const char* keyword) ;


#define Matrix_GetMatrixIndex(MAT)               ((MAT)->index)
#define Matrix_GetMatrixStorageFormat(MAT)       ((MAT)->fmt)
#define Matrix_GetNbOfRows(MAT)                  ((MAT)->n)
#define Matrix_GetNbOfColumns(MAT)               ((MAT)->n)
#define Matrix_GetNbOfEntries(MAT)               ((MAT)->len)
#define Matrix_GetNbOfNonZeroValues(MAT)         ((MAT)->nnz)
#define Matrix_GetNonZeroValue(MAT)              ((MAT)->nzval)
#define Matrix_GetWorkSpace(MAT)                 ((MAT)->work)
#define Matrix_GetStorage(MAT)                   ((MAT)->store)
#define Matrix_GetState(MAT)                     ((MAT)->state)



/* States of the matrix */
#define Matrix_InitialState    (0)
#define Matrix_ModifiedState   (1)

#define Matrix_SetToInitialState(MAT) \
        do {Matrix_GetState(MAT) = Matrix_InitialState ;} while(0)
        
        
#define Matrix_WasNotModified(MAT) \
        (Matrix_GetState(MAT) == Matrix_InitialState)
        
        
#define Matrix_SetToModifiedState(MAT) \
        do {Matrix_GetState(MAT) = Matrix_ModifiedState ;} while(0)



/* Initialize the matrix */
#define Matrix_SetValuesToZero(MAT) \
        do { \
          unsigned int Matrix_k ; \
          for(Matrix_k = 0 ; Matrix_k < Matrix_GetNbOfNonZeroValues(MAT) ; Matrix_k++) { \
            Matrix_GetNonZeroValue(MAT)[Matrix_k] = 0. ; \
          } \
          Matrix_GetNbOfEntries(MAT) = 0 ; \
          Matrix_SetToInitialState(MAT) ; \
        } while(0)




#include "MatrixStorageFormat.h"

#define Matrix_StorageFormatIs(MAT,KEY) \
        MatrixStorageFormat_Is(Matrix_GetMatrixStorageFormat(MAT),KEY)
        
        

struct Matrix_s {             /* Matrix */
  unsigned int index ;        /* Matrix index */
  MatrixStorageFormat_t fmt ; /* Storage format */
  unsigned int    n ;         /* Nb of rows/columns */
  unsigned int    nnz ;       /* Nb of non zero values */
  unsigned int    len ;       /* Nb of entries in nzval */
  double*         nzval ;     /* Pointer to the non zero values */
  void*   work ;              /* Pointer to a working memory space */
  void*   store ;             /* Pointer to the actual storage of the matrix */
  char    state ;             /* State of the matrix */
} ;

#endif
