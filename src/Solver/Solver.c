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

Solver_t*  (Solver_Create)(Mesh_t* mesh,Options_t* options,const int n_res,const int imatrix)
{
  Solver_t* solver = (Solver_t*) Mry_New(Solver_t) ;
  
  
  /* Resolution method */
  {
    Solver_GetResolutionMethod(solver) = ResolutionMethod_Create(options) ;
  }
  
  
  /* Nb of rows/columns */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
    
    Solver_GetNbOfColumns(solver) = n_col ;
  }
  
  
  /* Allocation of space for the residu */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    Residu_t* residu = Residu_Create(n_col,n_res) ;
    
    Residu_GetResiduIndex(residu) = imatrix ;
    
    Solver_GetResidu(solver) = residu ;
  }
  
  
  /* Allocation of space for the matrix */
  {
    Matrix_t* matrix = Matrix_Create(mesh,options,imatrix) ;
    
    Solver_GetMatrix(solver) = matrix ;
  }
  
  
  /* The solver */
  {
    if(Solver_ResolutionMethodIs(solver,CROUT)) {
      Solver_GetSolve(solver) = CroutMethod_Solve ;
      
    #ifdef SUPERLULIB
    } else if(Solver_ResolutionMethodIs(solver,SuperLU)) {
      Solver_GetSolve(solver) = SuperLUMethod_Solve ;
      
      /* Allocate work spaces for dgssvx */
      {
        int n = Solver_GetNbOfColumns(solver) ;
        void* etree  = Mry_New(int[n]) ;
        void* R      = Mry_New(double[n]) ;
        void* C      = Mry_New(double[n]) ;
        void* err    = Mry_New(double[2]) ;
        GenericData_t* getree = GenericData_Create(n,etree,int,"etree") ;
        GenericData_t* gR     = GenericData_Create(n,R,double,"R") ;
        GenericData_t* gC     = GenericData_Create(n,C,double,"C") ;
        GenericData_t* gerr   = GenericData_Create(2,err,double,"err") ;
        /* The FILL factor: FILL = 30 set by sp_ienv of superlu. */
        /* If lwork = 0 allocate space internally by system malloc */
        double fill = Options_GetFillFactor(options) ;
        Matrix_t* matrix = Solver_GetMatrix(solver) ;
        int nnz = Matrix_GetNbOfNonZeroValues(matrix) ;
        int lwork = floor(fill*nnz) ;

        Solver_AppendGenericWorkSpace(solver,getree) ;
        Solver_AppendGenericWorkSpace(solver,gR) ;
        Solver_AppendGenericWorkSpace(solver,gC) ;
        Solver_AppendGenericWorkSpace(solver,gerr) ;
      
        if(lwork) {
          void* work = Mry_New(double[lwork]) ;
          GenericData_t* gwork = GenericData_Create(lwork,work,double,"work") ;
        
          Solver_AppendGenericWorkSpace(solver,gwork) ;
        }
      }
    #endif
    
    #ifdef BLASLIB
    } else if(Solver_ResolutionMethodIs(solver,MA38)) {
      Solver_GetSolve(solver) = MA38Method_Solve ;

      /*  Allocate work space for ma38cd */
      {
        int n_col = Solver_GetNbOfColumns(solver) ;
        int lwork = 4 * n_col ;
        void* work = Mry_New(double[lwork]) ;
        GenericData_t* gwork = GenericData_Create(lwork,work,double,"work") ;
      
        Solver_AppendGenericWorkSpace(solver,gwork) ;
      }
    #endif
    
    } else {
      arret("Solver_Create: method not available") ;
    }
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
  
  {
    ResolutionMethod_t* rm = Solver_GetResolutionMethod(solver) ;
    
    if(rm) {
      ResolutionMethod_Delete(rm) ;
      free(rm) ;
      Solver_GetResolutionMethod(solver) = NULL ;
    }
  }
  
  {
    GenericData_t* genericwork = Solver_GetGenericWorkSpace(solver) ;
    
    if(genericwork) {
      GenericData_Delete(genericwork) ;
      free(genericwork) ;
      Solver_GetGenericWorkSpace(solver) = NULL ;
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

