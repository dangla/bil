#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "Mry.h"
#include "Solver.h"
#include "Message.h"
#include "BilLib.h"
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


#if 1
Solver_t*  Solver_Create(Mesh_t* mesh,Options_t* options,const int n)
{
  int NbOfMatrices = Mesh_GetNbOfMatrices(mesh) ;
  Solver_t* solver = (Solver_t*) Mry_New(Solver_t[NbOfMatrices]) ;
  
  {
    int i ;
    
    for(i = 0 ; i < NbOfMatrices ; i++) {
      Solver_t* solver_i = Solver_CreateSelectedMatrix(mesh,options,n,i) ;
      
      solver[i] = solver_i[0] ;
    }
  }
  
  return(solver) ;
}
#endif



Solver_t*  Solver_CreateSelectedMatrix(Mesh_t* mesh,Options_t* options,const int n,const int imatrix)
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
  
  #if 0
  /* Allocation of space for the right hand side */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* rhs = (double*) Mry_New(double[n*n_col]) ;
    
    Solver_GetRHS(solver) = rhs ;
  }
  
  
  /* Allocation of space for the solution */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* sol = (double*) Mry_New(double[n*n_col]) ;
    
    Solver_GetSolution(solver) = sol ;
  }
  #endif
  
  
  /* Allocation of space for the matrix */
  {
    Solver_GetMatrix(solver) = Matrix_CreateSelectedMatrix(mesh,options,imatrix) ;
  }
  
  return(solver) ;
}



#if 0
Solver_t*  Solver_Create(Mesh_t* mesh,Options_t* options,const int n)
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
    int n_col = Mesh_GetNbOfMatrixColumns(mesh)[0] ;
    
    Solver_GetNbOfColumns(solver) = n_col ;
  }
  
  
  /* Allocation of space for the residu */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    Residu_t* residu = Residu_Create(n_col,n) ;
    
    Residu_GetResiduIndex(residu) = imatrix ;
    
    Solver_GetResidu(solver) = residu ;
  }
  
  
  #if 0
  /* Allocation of space for the right hand side */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* rhs = (double*) Mry_New(double[n*n_col]) ;
    
    Solver_GetRHS(solver) = rhs ;
  }
  
  
  /* Allocation of space for the solution */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* sol = (double*) Mry_New(double[n*n_col]) ;
    
    Solver_GetSolution(solver) = sol ;
  }
  #endif
  
  
  #if 0 // Suppress 2021/07/26
  /* Update the system */
  {
    Mesh_UpdateMatrixRowColumnIndexes(mesh) ;
  }
  
  
  /* Print */
  {
    char*   debug  = Options_GetPrintedInfos(options) ;
    
    if(!strcmp(debug,"numbering")) Mesh_PrintData(mesh,debug) ;
  }
  #endif
  
  
  /* Allocation of space for the matrix */
  {
    Solver_GetMatrix(solver) = Matrix_Create(mesh,options) ;
  }
  
  return(solver) ;
}
#endif




void  Solver_Delete(void* self)
{
  Solver_t** psolver = (Solver_t**) self ;
  Solver_t*  solver  = *psolver ;
  Matrix_t* a = Solver_GetMatrix(solver) ;
  Residu_t* rs = Solver_GetResidu(solver) ;
  
  Matrix_Delete(&a) ;
  Residu_Delete(&rs) ;
  //free(Solver_GetRHS(solver)) ;
  //free(Solver_GetSolution(solver)) ;
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

