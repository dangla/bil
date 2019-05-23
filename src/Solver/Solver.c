#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "Solver.h"
#include "Message.h"
#include "BilLib.h"
#include "ResolutionMethod.h"
#include "CroutMethod.h"
#ifdef SUPERLULIB
  #include "SuperLUMethod.h"
#endif
#include "MA38Method.h"




/*
  Extern functions
*/

Solver_t*  Solver_Create(Mesh_t* mesh,Options_t* options,const int n)
{
  Solver_t* solver = (Solver_t*) malloc(sizeof(Solver_t)) ;
  
  assert(solver) ;
  
  
  /*  Method */
  {
    char* method = Options_GetResolutionMethod(options) ;
    
    if(!strcmp(method,"crout")) {
    
      Solver_GetResolutionMethod(solver) = ResolutionMethod_Type(CROUT) ;
      Solver_GetSolve(solver) = CroutMethod_Solve ;
    
    #ifdef SUPERLULIB
    } else if(!strcmp(method,"slu")) {
    
      Solver_GetResolutionMethod(solver) = ResolutionMethod_Type(SLU) ;
      Solver_GetSolve(solver) = SuperLUMethod_Solve ;
    #endif

    } else if(!strcmp(method,"ma38")) {
    
      Solver_GetResolutionMethod(solver) = ResolutionMethod_Type(MA38) ;
      Solver_GetSolve(solver) = MA38Method_Solve ;
      
    } else {
      arret("Solver_Create(1): unknown method") ;
    }
  }
  
  
  /* Nb of rows/columns */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
    
    Solver_GetNbOfColumns(solver) = n_col ;
  }
  
  
  /* Allocation of space for the matrix */
  {
    Solver_GetMatrix(solver) = Matrix_Create(mesh,options) ;
  }
  
  
  /* Allocation of space for the right hand side */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    size_t sz = n*n_col*sizeof(double) ;
    double* rhs = (double*) malloc(sz) ;
    
    assert(rhs) ;
    
    Solver_GetRHS(solver) = rhs ;
  }
  
  
  /* Allocation of space for the solution */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    size_t sz = n*n_col*sizeof(double) ;
    double* sol = (double*) malloc(sz) ;
    
    assert(sol) ;
    
    Solver_GetSolution(solver) = sol ;
  }
  
  return(solver) ;
}




void  Solver_Delete(void* self)
{
  Solver_t** psolver = (Solver_t**) self ;
  Solver_t*  solver  = *psolver ;
  Matrix_t* a = Solver_GetMatrix(solver) ;
  
  Matrix_Delete(&a) ;
  free(Solver_GetRHS(solver)) ;
  free(Solver_GetSolution(solver)) ;
  free(solver) ;
  *psolver = NULL ;
}



void Solver_Print(Solver_t* solver,char* keyword)
{
  static int i_debug=0 ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"debug(%d)\n",i_debug++) ;
  fprintf(stdout,"-----\n") ;
            
  if(!strcmp(keyword,"residu")) {
    double*  rhs = Solver_GetRHS(solver) ;
    int n_col = Solver_GetNbOfColumns(solver) ;
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"residu:\n") ;
    fprintf(stdout,"n = %d\n",n_col) ;
            
    for(i = 0 ; i < n_col ; i++) {
      fprintf(stdout,"res %d: % e\n",i,rhs[i]) ;
    }
  }
  
  if(!strncmp(keyword,"matrix",6)) {
    Matrix_t*  a = Solver_GetMatrix(solver) ;
    
    Matrix_PrintMatrix(a,keyword) ;
  }

}

