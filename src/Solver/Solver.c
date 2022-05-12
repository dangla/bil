#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "Mry.h"
#include "Solver.h"
#include "Message.h"
#include "BilExtraLibs.h"
#include "ResolutionMethod.h"
#include "CroutMethod.h"
#ifdef SUPERLULIB
  #include "SuperLUMethod.h"
#endif
#ifdef BLASLIB
#include "MA38Method.h"
#endif




/*
  Extern functions
*/

#if 0
Solver_t*  (Solver_Create)(Mesh_t* mesh,Options_t* options,const int n)
{
  int NbOfMatrices = Mesh_GetNbOfMatrices(mesh) ;
  Solver_t* solver = (Solver_t*) Mry_New(Solver_t[NbOfMatrices]) ;
  
  {
    int i ;
    
    for(i = 0 ; i < NbOfMatrices ; i++) {
      Solver_t* solver_i = Solver_CreateSelectedMatrix(mesh,options,n,i) ;
      
      solver[i] = solver_i[0] ;
      free(solver_i) ;
    }
  }
  
  return(solver) ;
}
#endif



Solver_t*  (Solver_Create)(Mesh_t* mesh,Options_t* options,const int n,const int imatrix)
{
  Solver_t* solver = (Solver_t*) Mry_New(Solver_t) ;
  
  
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

    #ifdef BLASLIB
    } else if(!strcmp(method,"ma38")) {
    
      Solver_GetResolutionMethod(solver) = ResolutionMethod_Type(MA38) ;
      Solver_GetSolve(solver) = MA38Method_Solve ;
    #endif
      
    } else {
      arret("Solver_Create(1): unknown method") ;
    }
  }
  
  
  /* Nb of rows/columns */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
    
    Solver_GetNbOfColumns(solver) = n_col ;
  }
  
  
  /* Allocation of space for the residu */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    Residu_t* residu = Residu_Create(n_col,n) ;
    
    Residu_GetResiduIndex(residu) = imatrix ;
    
    Solver_GetResidu(solver) = residu ;
  }
  
  
  /* Allocation of space for the matrix */
  {
    Matrix_t* matrix = Matrix_Create(mesh,options,imatrix) ;
    
    Solver_GetMatrix(solver) = matrix ;
  }
  
  return(solver) ;
}



void  (Solver_Delete)(void* self)
{
  Solver_t* solver = (Solver_t*) self ;
  
  {
    Matrix_t* a = Solver_GetMatrix(solver) ;
    
    if(a) {
      Matrix_Delete(a) ;
      free(a) ;
      Solver_GetMatrix(solver) = NULL ;
    }
  }
  
  {
    Residu_t* rs = Solver_GetResidu(solver) ;

    if(rs) {
      Residu_Delete(rs) ;
      free(rs) ;
      Solver_GetResidu(solver) = NULL ;
    }
  }
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

