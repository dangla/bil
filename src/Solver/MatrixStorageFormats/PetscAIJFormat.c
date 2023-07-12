#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "BilExtraLibs.h"
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "Mry.h"

#include "PetscAIJFormat.h"

#ifdef PETSCLIB
#include <petsc.h>
#endif




/* Extern functions */


#ifdef PETSCLIB

PetscAIJFormat_t* (PetscAIJFormat_Create)(Mesh_t* mesh,const int imatrix)
/* Create a matrix in PetscAIJFormat format */
{
  PetscAIJFormat_t* petscaij = (PetscAIJFormat_t*) Mry_New(PetscAIJFormat_t) ;
  Mat* aij = (Mat*) Mry_New(Mat) ;
  PetscInt ierror = MatCreate(PETSC_COMM_WORLD,aij)  ;
  
  /* undefined ref to PetscCall ??? */
  //PetscCall(ierror)  ;
  
  PetscAIJFormat_GetStorage(petscaij) = aij ;
  
  {
    int n = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
    
    MatSetSizes(*aij,PETSC_DECIDE,PETSC_DECIDE,n,n) ;
  }
  
  /* Nb of entries */
  {
    int nnz = Mesh_ComputeNbOfSelectedMatrixEntries(mesh,imatrix) ;
      
    PetscAIJFormat_GetNbOfNonZeroValues(petscaij) = nnz ;
  }
  
  return(aij) ;
}



void (PetscAIJFormat_Delete)(void* self)
{
  PetscAIJFormat_t* petscaij = (PetscAIJFormat_t*) self ;
  
  {
    Mat* aij = PetscAIJFormat_GetStorage(petscaij) ;
    
    MatDestroy(aij) ;
    free(aij) ;
  }
}

#endif
