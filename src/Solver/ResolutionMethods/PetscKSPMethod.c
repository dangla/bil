#include "BilExtraLibs.h"

#ifdef PETSCLIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Solver.h"
#include "Message.h"
#include "Matrix.h"

#include "PetscKSPMethod.h"
#include "PetscAIJFormat.h"
#include "DistributedMS.h"

#include <petsc.h>


int   PetscKSPMethod_Solve(Solver_t* solver)
/** Resolution of a.x = b by Petsc KSP's method */
{
  Matrix_t* matrix = Solver_GetMatrix(solver) ;
  Residu_t* residu = Solver_GetResidu(solver) ;
  PetscAIJFormat_t* petscaij = (PetscAIJFormat_t*) Matrix_GetStorage(matrix) ;
  Mat* A = (Mat*) PetscAIJFormat_GetStorage(petscaij) ;
  Vec* B = (Vec*) Residu_GetStoragOfRHS(residu) ;
  Vec* X = (Vec*) Residu_GetStoragOfSolution(residu) ;


  /* Solve the linear system with KSP method */
  {
    GenericData_t* gw = Solver_GetGenericWorkSpace(solver) ;
    KSP* ksp = GenericData_FindData(gw,KSP,"ksp") ;
    
    MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY) ;
    MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY) ;
    
    VecAssemblyBegin(*B) ;
    VecAssemblyEnd(*B) ;
    
    /* To solve successive linear systems that have different 
     * preconditioner matrices we must call KSPSetOperators().
     **/
    KSPSetOperators(*ksp,*A,*A) ;
    KSPSetUp(*ksp) ;
    KSPSolve(*ksp,*B,*X) ;
    
    /* Store the solution */
    {
      double* sol = Residu_GetSolution(residu) ;
      int n = Residu_GetLengthOfRHS(residu) ;
      double* array ;
      int i ;
      
      VecGetArray(*X,&array) ;
      
      {
        PetscInt low ;
        PetscInt high ;
        int i ;
      
        VecGetOwnershipRange(*X,&low,&high) ;
        
        for(i = 0 ; i < low ; i++) {
          sol[i] = 0 ;
        }
        for(i = low ; i < high ; i++) {
          int k = i - low ;
          
          sol[i] = array[k] ;
        }
        for(i = high ; i < n ; i++) {
          sol[i] = 0 ;
        }
      }
      
      VecRestoreArray(*X,&array) ;
        
  
      /* Broadcast to other processors */
      {
        int size = DistributedMS_NbOfProcessors ;
        
        if(size > 1) {
          //int rank = DistributedMS_RankOfCallingProcess ;
          
          //MPI_Bcast(sol0,narray,MPI_DOUBLE,rank,MPI_COMM_WORLD) ;
          MPI_Allreduce(MPI_IN_PLACE,sol,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ;
        }
      }
    }
    
    {
      Options_t* options = Solver_GetOptions(solver) ;
      char* info = Options_GetPrintLevel(options) ;
      
      if(String_Is(info,"iter")) {
        PetscInt its ;
        PetscReal norm ;
      
        VecNorm(*X,NORM_INFINITY,&norm) ;
        KSPGetIterationNumber(*ksp,&its) ;
        PetscPrintf(PETSC_COMM_WORLD,"Norm of solution %g iterations %" PetscInt_FMT "\n",(double) norm,its) ;
      }
    }
  }
  
  return(0) ;
}

#endif

