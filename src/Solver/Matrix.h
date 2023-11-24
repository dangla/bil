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
extern int         (Matrix_AssembleElementMatrix) (Matrix_t*,Element_t*,double*) ;
extern void        (Matrix_PrintMatrix)           (Matrix_t*,const char* keyword) ;
extern void        (Matrix_SetValuesToZero)       (Matrix_t*) ;


#define Matrix_GetMatrixIndex(MAT)               ((MAT)->index)
#define Matrix_GetMatrixStorageFormat(MAT)       ((MAT)->fmt)
#define Matrix_GetNbOfRows(MAT)                  ((MAT)->n)
#define Matrix_GetNbOfColumns(MAT)               ((MAT)->n)
#define Matrix_GetNbOfEntries(MAT)               ((MAT)->len)
#define Matrix_GetNbOfNonZeroValues(MAT)         ((MAT)->nnz)
#define Matrix_GetNonZeroValue(MAT)              ((MAT)->nzval)
#define Matrix_GetStorage(MAT)                   ((MAT)->store)
#define Matrix_GetState(MAT)                     ((MAT)->state)
#define Matrix_GetSparsityPattern(MAT)           ((MAT)->sparsitypattern)
#define Matrix_GetRowPermutation(MAT)            ((MAT)->rowpermutation)
#define Matrix_GetColumnPermutation(MAT)         ((MAT)->columnpermutation)
#define Matrix_GetGenericWorkSpace(MAT)          ((MAT)->genericwork)

#define Matrix_GetOptions(MAT) \
        MatrixStorageFormat_GetOptions(Matrix_GetMatrixStorageFormat(MAT))



/* Entry state of the matrix */
#define Matrix_UnfactorizedState      (0)
#define Matrix_FactorizedState        (1)

#define Matrix_SetToUnfactorizedState(MAT) \
        do {Matrix_GetState(MAT) = Matrix_UnfactorizedState ;} while(0)
        
        
#define Matrix_SetToFactorizedState(MAT) \
        do {Matrix_GetState(MAT) = Matrix_FactorizedState ;} while(0)
        
        
#define Matrix_IsFactorized(MAT) \
        (Matrix_GetState(MAT) == Matrix_FactorizedState)
        
        
#define Matrix_IsNotFactorized(MAT) \
        (Matrix_GetState(MAT) == Matrix_UnfactorizedState)



/* Sparsity pattern of the matrix */
#define Matrix_SparsityPatternOFF    (0)
#define Matrix_SparsityPatternON     (1)

#define Matrix_SetSameSparsityPattern(MAT) \
        do {Matrix_GetSparsityPattern(MAT) = Matrix_SparsityPatternON ;} while(0)

#define Matrix_UnsetSameSparsityPattern(MAT) \
        do {Matrix_GetSparsityPattern(MAT) = Matrix_SparsityPatternOFF ;} while(0)

#define Matrix_HasSameSparsityPattern(MAT) \
        (Matrix_GetSparsityPattern(MAT) == Matrix_SparsityPatternON)

#define Matrix_HasNotSameSparsityPattern(MAT) \
        (Matrix_GetSparsityPattern(MAT) == Matrix_SparsityPatternOFF)
        


#if 0
/* Initialize the matrix */
#define Matrix_SetValuesToZero(MAT) \
        do { \
          if(0) { \
          } else if(Matrix_StorageFormatIs(MAT,PetscAIJ)) { \
            PetscAIJFormat_t* petscaij = (PetscAIJFormat_t*) Matrix_GetStorage(MAT) ;
            Mat* A = (Mat*) PetscAIJFormat_GetStorage(petscaij) ;
            MatZeroEntries(*A) ; \
          } else {
            unsigned int Matrix_k ; \
            for(Matrix_k = 0 ; Matrix_k < Matrix_GetNbOfNonZeroValues(MAT) ; Matrix_k++) { \
              Matrix_GetNonZeroValue(MAT)[Matrix_k] = 0. ; \
            } \
          } \
          Matrix_GetNbOfEntries(MAT) = 0 ; \
          Matrix_SetToUnfactorizedState(MAT) ; \
        } while(0)
#endif




#include "MatrixStorageFormat.h"
#include "GenericData.h"

#define Matrix_GetStorageFormatType(MAT) \
        MatrixStorageFormat_GetType(Matrix_GetMatrixStorageFormat(MAT))

#define Matrix_StorageFormatIs(MAT,KEY) \
        MatrixStorageFormat_Is(Matrix_GetMatrixStorageFormat(MAT),KEY)



#define Matrix_AppendGenericWorkSpace(MAT,GD) \
        Matrix_GetGenericWorkSpace(MAT) = GenericData_Append(Matrix_GetGenericWorkSpace(MAT),GD)
        
        

struct Matrix_s {             /* Matrix */
  unsigned int index ;        /* Matrix index */
  MatrixStorageFormat_t* fmt ; /* Storage format */
  unsigned int    n ;         /* Nb of rows/columns */
  unsigned int    nnz ;       /* Nb of non zero values */
  unsigned int    len ;       /* Nb of entries in nzval */
  double*         nzval ;     /* Pointer to the non zero values */
  void*   store ;             /* Pointer to the actual storage of the matrix */
  GenericData_t* genericwork ; /* Working space */
  char    state ;             /* Entry state of the matrix */
  char    sparsitypattern ;   /* Sparsity pattern of the matrix */
  int*    rowpermutation ;
  int*    columnpermutation ;
} ;

#endif
