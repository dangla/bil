#include "BilExtraLibs.h"

#ifdef SUPERLULIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Solver.h"
#include "Message.h"
#include "Matrix.h"
#include "DistributedMS.h"

#include "SuperLUMethod.h"
#include "SuperLUFormat.h"

#include "superlu.h"


int   SuperLUMethod_Solve(Solver_t* solver)
/** Resolution of a.x = b by SuperLU's method */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  double*   b = Solver_GetRHS(solver) ;
  double*   x = Solver_GetSolution(solver) ;
  int       n = Solver_GetNbOfColumns(solver) ;
  
  /* for dgssv */
  SuperLUFormat_t* A = (SuperLUFormat_t*) Matrix_GetStorage(a) ;
  SuperLUFormat_t  L ;
  SuperLUFormat_t  U ;
  SuperLUFormat_t  B ;
  DNformat         Bstore ;
  SuperLUFormat_t  X ;
  DNformat         Xstore ;
  superlu_options_t slu_options ;
  SuperLUStat_t stat ;
  int   info ;
  int*  perm_r = Matrix_GetRowPermutation(a) ;
  int*  perm_c = Matrix_GetColumnPermutation(a) ;
  int rank = DistributedMS_RankOfCallingProcess ;
  int size = DistributedMS_NbOfProcessors ;
  
  
  if(rank == 0) {
    /* Options */
    /* Set the default input slu_options:
       slu_options->Fact = DOFACT;
       slu_options->Equil = YES;
       slu_options->ColPerm = COLAMD;
       slu_options->DiagPivotThresh = 1.0;
       slu_options->Trans = NOTRANS;
       slu_options->IterRefine = NOREFINE;
       slu_options->SymmetricMode = NO;
       slu_options->PivotGrowth = NO;
       slu_options->ConditionNumber = NO;
       slu_options->PrintStat = YES;
    */
    set_default_options(&slu_options) ;
    slu_options.PrintStat = NO ;
    //slu_options.ColPerm = NATURAL ;
    //slu_options.DiagPivotThresh = 0. ;
    //slu_options.IterRefine = EXTRA ;
  
  
    /* Factorization */
    slu_options.Fact = DOFACT ;
    if(Matrix_IsFactorized(a)) {
      slu_options.Fact = FACTORED ;
    } else if(Matrix_HasSameSparsityPattern(a)) {
      /* Same sparsity pattern and different numerical entries */
      slu_options.Fact = SamePattern ;
      /* Same sparsity pattern and similar numerical */
      //slu_options.Fact = SamePattern_SameRowPerm ;
    }


    /* Initialize the statistics variables */
    StatInit(&stat) ;


    /* The rhs */
    {
      B.Stype = SLU_DN ;
      B.Dtype = SLU_D ;
      B.Mtype = SLU_GE ;
      B.nrow  = n ;
      B.ncol  = 1 ;
      B.Store = (DNformat *) &Bstore ;
      Bstore.lda = n ;
      Bstore.nzval = (double*) b ;
    }

    /* The solution */
    {
      X.Stype = SLU_DN ;
      X.Dtype = SLU_D ;
      X.Mtype = SLU_GE ;
      X.nrow  = n ;
      X.ncol  = 1 ;
      X.Store = (DNformat *) &Xstore ;
      Xstore.lda = n ;
      Xstore.nzval = (double*) x ;
    }
  
    #if 0
    {
      dgssv(&slu_options,A,perm_c,perm_r,&L,&U,&B,&stat,&info) ;
    
      {
        int i ;
      
        for(i = 0 ; i < n ; i++) {
          x[i] = b[i] ;
        }
      }
    }
    #endif


    #if 1
    {
      double         rpg ;
      double         rcond ;
      mem_usage_t    mem_usage ;
      GlobalLU_t     Glu ;
      char           equed[1] ;
    
      {
        GenericData_t* gw = Solver_GetGenericWorkSpace(solver) ;
        double* ferr  = GenericData_FindData(gw,double,"err") ;
        double* berr  = ferr + 1 ;
        int*    etree = GenericData_FindData(gw,int,"etree") ;
        double* R     = GenericData_FindData(gw,double,"R") ;
        double* C     = GenericData_FindData(gw,double,"C") ;
        GenericData_t* gwork  = GenericData_Find(gw,double,"work") ;
        void* work   = (gwork) ? GenericData_GetData(gwork) : NULL ;
        size_t lwork = (gwork) ? GenericData_GetSize(gwork)*GenericData_GetNbOfData(gwork) : 0 ;
      
        if(Matrix_IsFactorized(a)) {
          equed[0] = 'B' ;
        }
      
        dgssvx(&slu_options,A,perm_c,perm_r,etree,equed,R,C,&L,&U,work,lwork,&B,&X,&rpg,&rcond, ferr,berr,&Glu,&mem_usage,&stat,&info) ;
      
        Matrix_SetToFactorizedState(a) ;
        Matrix_SetSameSparsityPattern(a) ;
      }
    }
    #endif


    if(slu_options.PrintStat) StatPrint(&stat) ;
  
    StatFree(&stat);
  }
  
  /* Broadcast to other processors */
  #if DistributedMS_APIis(MPI)
  if(size > 1) {
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
  }
  #endif
  
  return(0) ;
}

#endif
