#ifndef PETSCAIJFORMAT_H
#define PETSCAIJFORMAT_H


/* class-like structure "PetscAIJFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct PetscAIJFormat_s ; typedef struct PetscAIJFormat_s PetscAIJFormat_t ;


#include "Mesh.h"


extern PetscAIJFormat_t* (PetscAIJFormat_Create)(Mesh_t*,const int) ;
extern void              (PetscAIJFormat_Delete)(void*) ;



#if 1
/** The getters */
#define PetscAIJFormat_GetNbOfNonZeroValues(a)              ((a)->nnz)
#define PetscAIJFormat_GetNonZeroValue(a)                   ((a)->nzval)
#define PetscAIJFormat_GetIndexOfRow(a)                     ((a)->idwrow)
#define PetscAIJFormat_GetIndexOfColumn(a)                  ((a)->idxcol)
#define PetscAIJFormat_GetStorage(a)                        ((a)->store)
#endif



/* complete the structure types by using the typedef */
#if 1
struct PetscAIJFormat_s {
  int    nnz ;          /* nb of non zero values */
  double* nzval ;       /* Non zero values */
  int* idxrow ;         /* Indices of the rows */
  int* idxcol ;         /* Indices of the columns */
  void* store;          /* pointer to the actual storage of the matrix */
} ;
#endif



#endif
