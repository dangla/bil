#ifndef MATRIX_H
#define MATRIX_H

/* class-like structures "Matrix_t" and attributes */

/* vacuous declarations and typedef names */
struct Matrix_s       ; typedef struct Matrix_s       Matrix_t ;


#include "Mesh.h"
#include "Options.h"

extern Matrix_t*   (Matrix_Create)(Mesh_t*,Options_t*) ;
extern void        (Matrix_Delete)(void*) ;
extern void        (Matrix_AssembleElementMatrix_)(Matrix_t*,double*,int*,int*,int) ;
extern void        (Matrix_PrintMatrix)(Matrix_t*,const char* keyword) ;


#define Matrix_GetMatrixStorageFormat(MAT)       ((MAT)->fmt)
#define Matrix_GetNbOfRows(MAT)                  ((MAT)->n)
#define Matrix_GetNbOfColumns(MAT)               ((MAT)->n)
#define Matrix_GetNbOfNonZeroValues(MAT)         ((MAT)->nnz)
#define Matrix_GetNonZeroValue(MAT)              ((MAT)->s)
#define Matrix_GetWorkSpace(MAT)                 ((MAT)->work)
#define Matrix_GetStorage(MAT)                   ((MAT)->store)
#define Matrix_GetState(MAT)                     ((MAT)->state)



#define Matrix_SetToInitialState(MAT) \
        do {Matrix_GetState(MAT) = 0 ;} while(0)
        
        
#define Matrix_WasNotModified(MAT) \
        (Matrix_GetState(MAT) == 0)
        
        
#define Matrix_SetToModifiedState(MAT) \
        do {Matrix_GetState(MAT) = 1 ;} while(0)


#define Matrix_SetValuesToZero(MAT) \
        do { \
          unsigned int k ; \
          for(k = 0 ; k < Matrix_GetNbOfNonZeroValues(MAT) ; k++) { \
            Matrix_GetNonZeroValue(MAT)[k] = 0. ; \
            Matrix_SetToInitialState(MAT) ; \
          } \
        } while(0)

        


#include "Element.h"
#include "Node.h"

#define Matrix_AssembleElementMatrix(MAT,EL,KE) \
  do { \
    int  cole[Element_MaxNbOfNodes*Model_MaxNbOfEquations] ; \
    int  lige[Element_MaxNbOfNodes*Model_MaxNbOfEquations] ; \
    int  nn = Element_GetNbOfNodes(EL) ; \
    int  neq = Element_GetNbOfEquations(EL) ; \
    int  i ; \
    for(i = 0 ; i < nn ; i++) { \
      Node_t* node_i = Element_GetNode(EL,i) ; \
      int    j ; \
      for(j = 0 ; j < neq ; j++) { \
        int ij = i*neq + j ; \
        int jj_col = Element_GetUnknownPosition(EL)[ij] ; \
        int jj_row = Element_GetEquationPosition(EL)[ij] ; \
        cole[ij] = (jj_col >= 0) ? Node_GetMatrixColumnIndex(node_i)[jj_col] : -1 ; \
        lige[ij] = (jj_row >= 0) ? Node_GetMatrixRowIndex(node_i)[jj_row] : -1 ; \
      } \
    } \
    Matrix_AssembleElementMatrix_(MAT,KE,cole,lige,nn*neq) ; \
  } while(0)





#include "MatrixStorageFormat.h"

#define Matrix_StorageFormatIs(MAT,KEY) \
        MatrixStorageFormat_Is(Matrix_GetMatrixStorageFormat(MAT),KEY)
        
        

struct Matrix_s {             /* Matrix */
  MatrixStorageFormat_t fmt ; /* Storage format */
  unsigned int    n ;         /* Nb of rows/columns */
  unsigned int    nnz ;       /* Nb of non zero terms */
  double* s ;                 /* Pointer to the non zero terms */
  void*   work ;              /* Pointer to a working memory space */
  void*   store ;             /* Pointer to the actual storage of the matrix */
  char    state ;             /* State of the matrix */
} ;

#endif
