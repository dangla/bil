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
  
  PetscAIJFormat_GetStorage(petscaij) = aij ;
  
  /* Create the matrix */
  {
    int n = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
    
    MatCreate(PETSC_COMM_WORLD,aij)  ;
    MatSetSizes(*aij,PETSC_DECIDE,PETSC_DECIDE,n,n) ;
    MatSetFromOptions(*aij) ;
    MatSetType(*aij,MATAIJ) ;
  }

  /* Preallocate the seq matrix aij */
  {
    int* nnzrow = Mesh_ComputeNbOfMatrixNonzerosPerRowAndColumn(mesh,imatrix) ;
    
    MatSeqAIJSetPreallocation(*aij,0,nnzrow) ;
    free(nnzrow) ;
  }
    
  /* The nb of processes and the rank of the calling process */
  {
    PetscMPIInt rank ;
    PetscMPIInt size ;
      
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n", size, rank);
  }
    
  /* Preallocate the MPI matrix aij */
  {
    PetscInt Istart ;
    PetscInt Iend ;
      
    MatGetOwnershipRange(*aij,&Istart,&Iend) ;
    
    {
      int nlocalrows = Iend - Istart ;
      int* d_nnzrow = Mesh_ComputeNbOfSubmatrixNonzerosPerRow(mesh,imatrix,Istart,Iend) ;
      int* o_nnzrow = d_nnzrow + nlocalrows ;
        
      MatMPIAIJSetPreallocation(*aij,0,d_nnzrow,0,o_nnzrow) ;
      free(d_nnzrow) ;
    }
  }
   
  /* Nb of entries */
  {
    int nnz = Mesh_ComputeNbOfSelectedMatrixEntries(mesh,imatrix) ;
      
    PetscAIJFormat_GetNbOfNonZeroValues(petscaij) = nnz ;
  }
  
  return(petscaij) ;
}



void (PetscAIJFormat_Delete)(void* self)
{
  PetscAIJFormat_t* petscaij = (PetscAIJFormat_t*) self ;
  
  {
    Mat* aij = PetscAIJFormat_GetStorage(petscaij) ;
    
    if(aij) {
      MatDestroy(aij) ;
      free(aij) ;
    }
  }
}

#endif
