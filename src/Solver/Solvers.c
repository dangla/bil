#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "Mry.h"
#include "Solvers.h"
#include "Message.h"
#include "BilExtraLibs.h"



#if defined (PETSCLIB)
  #include <petsc.h>
#endif

/*
  Extern functions
*/


Solvers_t*  (Solvers_Create)(Mesh_t* mesh,Options_t* options,const int n)
{
  int NbOfMatrices = Mesh_GetNbOfMatrices(mesh) ;
  Solvers_t* solvers = (Solvers_t*) Mry_New(Solvers_t) ;
  
  Solvers_GetNbOfSolvers(solvers) = NbOfMatrices ;
  
  {
    Solver_t* solver = (Solver_t*) Mry_New(Solver_t[NbOfMatrices]) ;
    int i ;
    
    for(i = 0 ; i < NbOfMatrices ; i++) {
      Solver_t* solver_i = Solver_Create(mesh,options,n,i) ;
      
      solver[i] = solver_i[0] ;
      free(solver_i) ;
    }
    
    Solvers_GetSolver(solvers) = solver ;
  }
  
  return(solvers) ;
}



void  (Solvers_Delete)(void* self)
{
  Solvers_t* solvers = (Solvers_t*) self ;
  
  {
    int n = Solvers_GetNbOfSolvers(solvers) ;
    Solver_t* solver = Solvers_GetSolver(solvers) ;
    
    if(solver) {
      Mry_Delete(solver,n,Solver_Delete) ;
      free(solver) ;
      Solvers_GetSolver(solvers) = NULL ;
    }
  }
}
