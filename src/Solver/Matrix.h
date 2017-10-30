#ifndef MATRIX_H
#define MATRIX_H

/* class-like structures "Matrix_t" and attributes */

/* vacuous declarations and typedef names */
struct Matrix_s       ; typedef struct Matrix_s       Matrix_t ;


#include "Mesh.h"
#include "Options.h"

extern Matrix_t*   Matrix_Create(Mesh_t*,Options_t*) ;
extern void Matrix_AssembleElementMatrix(Matrix_t*,double*,int*,int*,int) ;
extern void Matrix_PrintMatrix(Matrix_t*,const char* keyword) ;


#define Matrix_GetMatrixFormat(matrix)              ((matrix)->fmt)
#define Matrix_GetNbOfRows(matrix)                  ((matrix)->n)
#define Matrix_GetNbOfColumns(matrix)               ((matrix)->n)
#define Matrix_GetNbOfNonZeroValues(matrix)         ((matrix)->nnz)
#define Matrix_GetNonZeroValue(matrix)              ((matrix)->s)
#define Matrix_GetWorkSpace(matrix)                 ((matrix)->work)
#define Matrix_GetStorage(matrix)                   ((matrix)->store)
#define Matrix_GetState(matrix)                     ((matrix)->state)



#include "MatrixStorageFormat.h"

struct Matrix_s {             /* Matrix */
  MatrixFormat_t fmt ;        /* Storage format */
  unsigned int    n ;         /* Nb of rows/columns */
  unsigned int    nnz ;       /* Nb of non zero terms */
  double* s ;                 /* Pointer to the non zero terms */
  void*   work ;              /* Pointer to a working memory space */
  void*   store ;             /* Pointer to the actual storage of the matrix */
  char    state ;             /* State of the matrix */
} ;

#endif
