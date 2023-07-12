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

#include <petsc.h>


int   PetscKSPMethod_Solve(Solver_t* solver)
/** Resolution of a.x = b by Petsc KSP's method */
{
  Matrix_t* matrix = Solver_GetMatrix(solver) ;
  Residu_t* residu = Solver_GetResidu(solver) ;
  
  Mat* A = (Mat*) Matrix_GetStorage(matrix) ;
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

    KSPSolve(*ksp,*B,*X) ;
  }
#endif
  
  return(0) ;
}

#endif

