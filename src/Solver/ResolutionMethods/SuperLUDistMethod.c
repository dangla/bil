#include "BilExtraLibs.h"

#ifdef SUPERLUDISTLIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Solver.h"
#include "Message.h"
#include "Matrix.h"

#include "SuperLUDistMethod.h"
#include "SuperLUFormat.h"

#include "superlu.h"


int   SuperLUDistMethod_Solve(Solver_t* solver)
/** Resolution of a.x = b by SuperLU's method */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  double*   b = Solver_GetRHS(solver) ;
  double*   x = Solver_GetSolution(solver) ;
  int       n = Solver_GetNbOfColumns(solver) ;

  SuperLUFormat_t* A = (SuperLUFormat_t*) Matrix_GetStorage(a) ;
  superlu_dist_options_t sludist_options ;
  SuperLUStat_t stat;
  int   info ;


  {
    
    {
      GenericData_t* gw = Solver_GetGenericWorkSpace(solver) ;
      double berr[1] ;
      dScalePermstruct_t* ScalePermstruct = GenericData_FindData(gw,dScalePermstruct_t,"ScalePermstruct") ;
      dLUstruct_t* LUstruct = GenericData_FindData(gw,dLUstruct_t,"LUstruct") ;
      gridinfo_t* grid = GenericData_FindData(gw,gridinfo_t,"grid") ;
      
  
      /* Set the default input sludist_options: */
      /*
      sludist_options.Fact = DOFACT;
      sludist_options.Equil = YES;
      sludist_options.ParSymbFact = NO;
      #ifdef HAVE_PARMETIS
      sludist_options.ColPerm = METIS_AT_PLUS_A;
      #else
      sludist_options.ColPerm = MMD_AT_PLUS_A;
      #endif
      sludist_options.RowPerm = LargeDiag_MC64;
      sludist_options.ReplaceTinyPivot = NO;
      sludist_options.IterRefine = SLU_DOUBLE;
      sludist_options.Trans = NOTRANS;
      sludist_options.SolveInitialized = NO;
      sludist_options.RefineInitialized = NO;
      sludist_options.PrintStat = YES;
      sludist_options.lookahead_etree = NO;
      sludist_options.num_lookaheads = 10;
      sludist_options.superlu_maxsup = 256;
      sludist_options.superlu_relax = 60;
      strcpy(sludist_options.superlu_rankorder,"Z"); 
      strcpy(sludist_options.superlu_lbs,"GD");
      sludist_options.superlu_acc_offload = 1;
      sludist_options.superlu_n_gemm = 5000;
      sludist_options.superlu_max_buffer_size = 256000000;
      sludist_options.superlu_num_gpu_streams = 8;
      sludist_options.SymPattern = NO;
      sludist_options.Algo3d = NO;
      #ifdef SLU_HAVE_LAPACK
      sludist_options.DiagInv = YES;
      #else
      sludist_options.DiagInv = NO;
      #endif
      sludist_options.Use_TensorCore = NO;
      */
      set_default_options_dist(&sludist_options);
      sludist_options.ColPerm = NATURAL;
      sludist_options.RowPerm = NATURAL;
      sludist_options.PrintStat = NO;

      /* Initialize the statistics variables. */
      PStatInit(&stat);
  
      /* Factorization */
      if(Matrix_IsFactorized(a)) {
        sludist_options.Fact = FACTORED ;
      } else if(Matrix_HasSameSparsityPattern(a)) {
        /* Same sparsity pattern and different numerical entries */
        sludist_options.Fact = SamePattern ;
        /* Same sparsity pattern and similar numerical */
        //sludist_options.Fact = SamePattern_SameRowPerm ;
        //sludist_options.refact = YES ;
      } else {
      }

      pdgssvx_ABglobal(&sludist_options,A,ScalePermstruct,b,n,1,grid,LUstruct,berr,&stat,&info);
      
      Matrix_SetToFactorizedState(a) ;
      Matrix_SetSameSparsityPattern(a) ;
      
      if(info == 0) {
        int i ;
        
        for(i = 0 ; i < n ; i++) {
          x[i] = b[i] ;
        }
      }
      
      if(info < 0) {
        printf("SuperLUDistMethod_Solve: illegal %d-th argument in pdgssvx_ABglobal()\n",-info) ;
        return(-1) ;
      }
      
      
      #if 1
      if (info > 0 && info <= n) {
        
        printf("The matrix is singular to working precision.\n");
        printf("U(%d,%d) is exactly zero.\n",info);
        
        fflush(stdout);
        return(-1) ;
      }
      
      if (info > n) {
        
        printf("Number of bytes allocated when memory allocation failure occurred.\n",info - n) ;
        
        fflush(stdout);
        return(-1) ;
      }
      #endif
    }
  }
  
  return(0) ;
}

#endif

