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
#if 1
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
    KSPSolve(*ksp,*B,*X) ;
    
    {
      Options_t* options = Solver_GetOptions(solver) ;
      char* info = Options_GetPrintLevel(options) ;
      //int nprocs = Options_NbOfThreadsInSolver(options) ;
      
      if(String_Is(info,"iter")) {
        PetscInt its ;
        PetscReal norm ;
      
        VecNorm(*X,NORM_2,&norm) ;
        KSPGetIterationNumber(*ksp,&its) ;
        PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %" PetscInt_FMT "\n",(double) norm,its) ;
      }
    }
  }
#endif
  
  return(0) ;
}

#endif

