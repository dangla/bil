#include "BilExtraLibs.h"

#ifdef SUPERLUMTLIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Solver.h"
#include "Message.h"
#include "Matrix.h"
#include "DistributedMS.h"

#include "SuperLUMTMethod.h"
#include "SuperLUFormat.h"

#include "superlu.h"



int   SuperLUMTMethod_Solve(Solver_t* solver)
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
  superlumt_options_t slumt_options ;
  int   info ;
  int*  perm_r = Matrix_GetRowPermutation(a) ;
  int*  perm_c = Matrix_GetColumnPermutation(a) ;
  Options_t* options = Solver_GetOptions(solver) ;
  //int nprocs = Options_NbOfThreadsInSolver(options) ;
  int nprocs = Options_NbOfThreads(options) ;
  int rank = DistributedMS_RankOfCallingProcess ;
  int size = DistributedMS_NbOfProcessors ;


  if(rank == 0) {
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


    #if 1
    {
      double         rpg ;
      double         rcond ;
      superlu_memusage_t    superlu_memusage ;
      //GlobalLU_t     Glu ;
      equed_t        equed[1] = {NOEQUIL} ;
    
      {
        GenericData_t* gw = Solver_GetGenericWorkSpace(solver) ;
        double* ferr  = GenericData_FindData(gw,double,"err") ;
        double* berr  = ferr + 1 ;
        int*    etree = GenericData_FindData(gw,int,"etree") ;
        int*    colcnt_h = GenericData_FindData(gw,int,"colcnt_h") ;
        int*    part_super_h = GenericData_FindData(gw,int,"part_super_h") ;
        double* R     = GenericData_FindData(gw,double,"R") ;
        double* C     = GenericData_FindData(gw,double,"C") ;
        GenericData_t* gwork  = GenericData_Find(gw,double,"work") ;
        void* work   = (gwork) ? GenericData_GetData(gwork) : NULL ;
        size_t lwork = (gwork) ? GenericData_GetSize(gwork)*GenericData_GetNbOfData(gwork) : 0 ;
      
  
        /* Options */
        slumt_options.nprocs = nprocs;
        slumt_options.fact = DOFACT;
        slumt_options.trans = NOTRANS;
        slumt_options.refact = NO;
        slumt_options.panel_size = sp_ienv(1);
        slumt_options.relax = sp_ienv(2);
        slumt_options.diag_pivot_thresh = 1;
        slumt_options.drop_tol = 0;
        slumt_options.ColPerm = COLAMD;
        slumt_options.usepr = NO;
        slumt_options.SymmetricMode = NO;
        slumt_options.PrintStat = NO;
        slumt_options.perm_c = perm_c;
        slumt_options.perm_r = perm_r;
        slumt_options.work = work;
        slumt_options.lwork = lwork;
        slumt_options.etree = etree;
        slumt_options.colcnt_h = colcnt_h;
        slumt_options.part_super_h = part_super_h;
  
  
        /* Factorization */
        if(Matrix_IsFactorized(a)) {
          slumt_options.fact = FACTORED ;
        } else if(Matrix_HasSameSparsityPattern(a)) {
          /* Same sparsity pattern and different numerical entries */
          //slumt_options.fact = SamePattern ;
          /* Same sparsity pattern and similar numerical */
          //slumt_options.Fact = SamePattern_SameRowPerm ;
          //slumt_options.refact = YES ;
        } else {
        /*
         * Get column permutation vector perm_c[], according to permc_spec:
         *   permc_spec = 0: natural ordering 
         *   permc_spec = 1: minimum degree ordering on structure of A'*A
         *   permc_spec = 2: minimum degree ordering on structure of A'+A
         *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
         */ 
          get_perm_c(3,A,perm_c);
        }
      
        if(Matrix_IsFactorized(a)) {
          //equed[0] = 'B' ;
        }

        pdgssvx(nprocs,&slumt_options,A,perm_c,perm_r,equed,R,C,&L,&U,&B,&X,&rpg,&rcond,ferr,berr,&superlu_memusage,&info) ;
      
        Matrix_SetToFactorizedState(a) ;
        Matrix_SetSameSparsityPattern(a) ;
      
        if(info > n+1) {
          printf("SuperLUMTMethod_Solve: memory allocation failed in pdgssvx()\n") ;
          return(-1) ;
        }
      
        if(info < 0) {
          printf("SuperLUMTMethod_Solve: illegal %d-th argument in pdgssvx()\n",-info) ;
          return(-1) ;
        }
      
      
        #if 1
        if (info == n+1) {
          int i ;
          int nrhs = 1 ;
          SCPformat* Lstore = (SCPformat *) L.Store;
          NCPformat* Ustore = (NCPformat *) U.Store;
        
          printf("The matrix is singular to working precision.\n");

          printf("Recip. pivot growth = %e\n", rpg);
          printf("Recip. condition number = %e\n", rcond);
          printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
          for (i = 0; i < nrhs; ++i) {
            printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
          }

          printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
          printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
          printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
          printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,superlu_memusage.expansions);
        
          fflush(stdout);
          return(-1) ;
        }
        #endif
      
        #if 0
        if (info > 0 && lwork == -1 ) {
          printf("** Estimated memory: " IFMT " bytes\n", info - n);
          printf("pdgssvx(): info " IFMT "\n----\n", info);
        }
        #endif
      }
    }
    #endif
  }
  
  /* Broadcast to other processors */
  if(size > 1) {
    #if DistributedMS_APIis(MPI)
      MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
    #else
      #error "Distributed memory system not available"
    #endif
  }
  
  return(0) ;
}
#endif

