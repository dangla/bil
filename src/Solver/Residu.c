#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "BilExtraLibs.h"
#include "Residu.h"


#if defined (PETSCLIB)
  #include <petsc.h>
#endif



/* Extern functions */

Residu_t*   (Residu_Create)(Mesh_t* mesh,Options_t* options,const int n_res,const int imatrix)
{
  Residu_t* residu = (Residu_t*) Mry_New(Residu_t) ;
  int n_col = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;

  Residu_GetLengthOfRHS(residu) = n_col ;
  Residu_GetNbOfRHS(residu) = n_res ;
  Residu_GetResiduIndex(residu) = imatrix ;
  
  
  /* Vector storage format */
  {
    Residu_GetVectorStorageFormat(residu) = VectorStorageFormat_Create(options) ;
  }


  /* Allocation of space for the right hand sides */
  {
    void* rhs = Mry_New(double[n_res*n_col]) ;
    
    Residu_GetRHS(residu) = rhs ;
  }
  
  
  /* Allocation of space for the solutions */
  {
    void* sol = Mry_New(double[n_res*n_col]) ;
    
    Residu_GetSolution(residu) = sol ;
  }


  #if defined (PETSCLIB)
  {
    if(Residu_StorageFormatIs(residu,PetscVec)) {
      
      /* Initialization */
      {
        Context_t* ctx = Options_GetContext(options) ;
        CommandLine_t* cmd = Context_GetCommandLine(ctx) ;
        int argc = CommandLine_GetNbOfArg(cmd) ;
        char** argv = CommandLine_GetArg(cmd) ;
        
        PetscInitialize(&argc,&argv,NULL,NULL) ;
      }

      /* The rhs */
      {
        double* rhs = (double*) Residu_GetRHS(residu) ;
        Vec* B = (Vec*) Mry_New(Vec) ;
        PetscInt n = n_col ;
        
        VecCreate(PETSC_COMM_WORLD,B) ;
        VecSetType(*B,VECSTANDARD) ;
        VecSetSizes(*B,PETSC_DECIDE,n) ;
        VecSetOption(*B,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE) ;
        VecSetFromOptions(*B) ;
        
        /* Create a preallocated seq vector */
        //VecCreateSeqWithArray(PETSC_COMM_WORLD,1,n,rhs,B) ;
        /* Create a preallocated MPI vector */
        //VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,n,rhs,B) ;
        
        VecPlaceArray(*B,rhs) ;
        
        Residu_GetStoragOfRHS(residu) = B ;
      }

      /* The solutions */
      {
        double* sol = (double*) Residu_GetSolution(residu) ;
        Vec* B = Residu_GetStoragOfRHS(residu) ;
        Vec* X = (Vec*) Mry_New(Vec) ;
        PetscInt n = n_col ;
        
        VecDuplicate(*B,X);
        //VecCreateSeqWithArray(PETSC_COMM_WORLD,1,n,sol,X) ;
        
        VecPlaceArray(*X,sol) ;
        
        Residu_GetStoragOfSolution(residu) = X ;
      }
    }
  }
  #endif
  
  
  return(residu) ;
}



void (Residu_Delete)(void* self)
{
  Residu_t* residu = (Residu_t*) self ;
  
  #if defined (PETSCLIB)
  {
    if(Residu_StorageFormatIs(residu,PetscVec)) {
      
      {
        Vec* B = Residu_GetStoragOfRHS(residu) ;
        
        if(B) {
          VecDestroy(B) ;
          free(B) ;
        }
      }

      {
        Vec* X = Residu_GetStoragOfSolution(residu) ;
        
        if(X) {
          VecDestroy(X) ;
          free(X) ;
        }
      }
    }
  }
  #endif
  
  {
    void* rhs = Residu_GetRHS(residu) ;
    
    if(rhs) {
      free(rhs) ;
      Residu_GetRHS(residu) = NULL ;
    }
  }
  
  {
    void* sol = Residu_GetSolution(residu) ;
    
    if(sol) {
      free(sol) ;
      Residu_GetSolution(residu) = NULL ;
    }
  }
}





void (Residu_AssembleElementResidu)(Residu_t* residu,Element_t* el,double* re)
{
  int imatrix = Residu_GetResiduIndex(residu) ;
  int  ndof = Element_GetNbOfDOF(el) ;
  int* row = Element_ComputeSelectedMatrixRowAndColumnIndices(el,imatrix) ;
  int* col = row + ndof ;
            
  if(Element_GetMaterial(el)) {
    if(!Residu_StorageFormatIs(residu,PetscVec)) {
      double* r = (double*) Residu_GetRHS(residu) ;
      int i ;
      
      for(i = 0 ; i < ndof ; i++) {
        int k = col[i] ;
        
        if(k >= 0) r[k] += re[i] ;
      }
      
      #if 0
      for(i = 0 ; i < nn ; i++) {
        Node_t* node_i = Element_GetNode(el,i) ;
        int j ;
              
        for(j = 0 ; j < neq ; j++) {
          int ij = i*neq + j ;
          int ii = Element_GetUnknownPosition(el)[ij] ;
                
          if(ii >= 0) {
            int k = Node_GetSelectedMatrixColumnIndexOf(node_i,ii,imatrix) ;
            if(k >= 0) r[k] += re[ij] ;
          }
        }
      }
      #endif
    
    #ifdef PETSCLIB
    /* format used in Petsc */
    } else if(Residu_StorageFormatIs(residu,PetscVec)) {
      Vec* B = Residu_GetStoragOfRHS(residu) ;
        
      VecSetValues(*B,ndof,col,re,ADD_VALUES) ;
    #endif
    
    } else {
      arret("Residu_AssembleElementResidu: unknown format") ;
    }
  }
  
  Element_FreeBufferFrom(el,row) ;
}



void Residu_PrintResidu(Residu_t* residu,const char* keyword)
{
  if(!Residu_StorageFormatIs(residu,PetscVec)) {
    double*  rhs = (double*) Residu_GetRHS(residu) ;
    int n_col = Residu_GetLengthOfRHS(residu) ;
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"residu:\n") ;
    fprintf(stdout,"n = %d\n",n_col) ;
            
    for(i = 0 ; i < n_col ; i++) {
      fprintf(stdout,"res %d: % e\n",i,rhs[i]) ;
    }
    
  #ifdef PETSCLIB
  /* format used in Petsc */
  } else if(Residu_StorageFormatIs(residu,PetscVec)) {
    Vec* b = Residu_GetStoragOfRHS(residu) ;
    
    VecView(*b,PETSC_VIEWER_STDOUT_WORLD) ;
  #endif

  } else {
    arret("Residu_PrintResidu: unknown format") ;
  }

  fflush(stdout) ;
}



void (Residu_SetValuesToZero)(Residu_t* residu)
/** zeros each element of the residu */
{
  if(0) {
  #ifdef PETSCLIB
  } else if(Residu_StorageFormatIs(residu,PetscVec)) {
    Vec* b = Residu_GetStoragOfRHS(residu) ;
    
    VecZeroEntries(*b) ;
  #endif
  } else {
    unsigned int n = Residu_GetNbOfRHS(residu)*Residu_GetLengthOfRHS(residu) ;
    double* rhs = (double*) Residu_GetRHS(residu) ;
    unsigned int k ;
          
    for(k = 0 ; k < n ; k++) {
      rhs[k] = 0. ;
    }
  }
}
