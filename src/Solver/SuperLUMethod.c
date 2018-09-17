#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "BilLib.h"
#include "SuperLUMethod.h"
#include "Solver.h"
#include "Message.h"
#include "Matrix.h"
#include "SuperLUFormat.h"



#ifdef SUPERLULIB

  #define val(x) #x
  #define xval(x) val(x)
  //#include xval(SUPERLULIB/SRC/dsp_defs.h)
  //#include "/home/dangla/Documents/Softwares/getfem-4.0.0/superlu/slu_ddefs.h"
  #include "superlu.h"
  #undef val
  #undef xval

#endif


#ifdef SUPERLULIB

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
  superlu_options_t options ;
  SuperLUStat_t stat ;
  int   info ;
  /* int   *perm_r = a.work,*perm_c = (int *) a.work + n ; */
  static int*   perm_r ;
  static int*   perm_c ;
  /* for dgssvx */
  static int    iresol = 0 ;
  SuperLUFormat_t X ;
  DNformat        Xstore ;
  
  /* Allocate memory space for the permutations of rows and columns */
  if(!iresol) {
    perm_r = (int*) malloc(n*sizeof(int)) ;
    
    assert(perm_r) ;
    
    perm_c = (int*) malloc(n*sizeof(int)) ;
    
    assert(perm_c) ;
  }
  
  
  /* Options */
  /* Set the default input options:
     options->Fact = DOFACT;
     options->Equil = YES;
     options->ColPerm = COLAMD;
     options->DiagPivotThresh = 1.0;
     options->Trans = NOTRANS;
     options->IterRefine = NOREFINE;
     options->SymmetricMode = NO;
     options->PivotGrowth = NO;
     options->ConditionNumber = NO;
     options->PrintStat = YES;
  */
  set_default_options(&options) ;
  options.PrintStat = NO ;
  //options.ColPerm = NATURAL ;
  //options.DiagPivotThresh = 0. ;
  //options.IterRefine = EXTRA ;
  if(iresol) {
    options.Fact = SamePattern ;
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
#if 0
  if(!iresol) {
    dCreate_Dense_Matrix(&B,n,1,b,n,SLU_DN,SLU_D,SLU_GE) ;
  }
#endif
  

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
  if(!iresol) {
    dCreate_Dense_Matrix(&X,n,1,x,n,SLU_DN,SLU_D,SLU_GE) ;
  }
#endif


#if 0
  {
    dgssv(&options,A,perm_c,perm_r,&L,&U,&B,&stat,&info) ;
    
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
    /* facilitate multiple factorizations with SamePattern_SameRowPerm */
    GlobalLU_t     Glu ;
    char           equed[1] ;
    static double  ferr[1] ;
    static double  berr[1] ;
    static int*    etree ;
    static double* R ;
    static double* C ;
    static void*   work ;
    /* The FILL factor */
    int            ff = 0 ;
    /* lwork = 0 allocate space internally by system malloc */
    int     lwork = ff * Matrix_GetNbOfNonZeroValues(a) ;
  
    /* Allocate memory space for the static variables */
    if(!iresol) {
      etree = (int*) malloc(n*sizeof(int)) ;
      
      assert(etree) ;
      
      R = (double*) malloc(n*sizeof(double)) ;
      
      assert(R) ;
      
      C = (double*) malloc(n*sizeof(double)) ;
      
      assert(C) ;
      
      if(lwork > 0) {
        work = (double*) malloc(lwork*sizeof(double)) ;
        
        assert(work) ;
      }
    }
    
    {
      dgssvx(&options,A,perm_c,perm_r,etree,equed,R,C,&L,&U,work,lwork,&B,&X,&rpg,&rcond,ferr,berr,&Glu,&mem_usage,&stat,&info) ;
    }
  }
#endif


  if(options.PrintStat) StatPrint(&stat) ;
  
  StatFree(&stat);

  iresol++ ;
  
  return(0) ;
}

#endif

