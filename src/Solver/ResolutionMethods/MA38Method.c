
#include "BilExtraLibs.h"

#if defined (BLASLIB) && defined (LAPACKLIB)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#include "Solver.h"
#include "Message.h"
#include "Matrix.h"
#include "Math_.h"

#include "MA38Method.h"
#include "CoordinateFormat.h"
#include "DistributedMS.h"


static void MA38Method_PrintErrorDiagnostics(double*,int*,int*) ;


/* 
 * HSL - MA38 method
 * -----------------
 */
int   MA38Method_Solve(Solver_t* solver)
/** Solve a sparse unsymmetric system a.x = b 
 *  by a multifrontal approach through ma38 from HSL 
 */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  int     n = Solver_GetNbOfColumns(solver) ;
  
  CoordinateFormat_t* ac = (CoordinateFormat_t*) Matrix_GetStorage(a) ;
  int    lvalue = CoordinateFormat_GetLengthOfArrayValue(ac) ;
  int    lindex = CoordinateFormat_GetLengthOfArrayIndex(ac) ;
  double* value = CoordinateFormat_GetNonZeroValue(ac) ;
  int*    index = CoordinateFormat_GetIndex(ac) ;
  
  int    keep[20] ;
  int    icntl[20] ;
  double cntl[10] ;
  int    info[40] ;
  double rinfo[20] ;
  int rank = DistributedMS_RankOfCallingProcess ;
  int size = DistributedMS_NbOfProcessors ;
  
  
  if(rank == 0) {
    /* Initialize controls */
    ma38id_(keep,cntl,icntl) ;
  
    /* Printing controls:
     * 1: Print only error messages, 
     * 2: + warnings, 
     * 3: + terse diagnostics 
     */
    {
      Options_t* options = Matrix_GetOptions(a) ;
    
      icntl[2] = 1 ;
      //icntl[2] = 5 ;
      icntl[7] = 1 ; // Max nb of steps of iterative refinement.
    
      #if 0
      if(Options_IsToPrintOutAtEachIteration(options)) {
        icntl[2] = 5 ;
      }
      #endif
    }
  
    /* Factorize the matrix */
    {
      int ne = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      int job = 1 ; // 1 if icntl[7] > 0
      bool transa = false ;
  
      ma38ad_(&n,&ne,&job,&transa,&lvalue,&lindex,value,index,keep,cntl,icntl,info,rinfo) ;
    
      if(info[0] < 0) {
        double ffvalue = ((double) lvalue)/(2*ne) ;
        double ffindex = ((double) lindex)/(3*ne + 2*n + 1) ;
      
        printf("MA38Method_Solve: error in the factorization\n") ;
      
        printf("Length of INDEX allocated: %d\n",lindex) ;
      
        printf("Length of VALUE allocated: %d\n",lvalue) ;
  
        //printf("Number of entries: %d\n",ne) ;
        //printf("Dimension of the matrix: %d\n",n) ;
      
        printf("Fill factors = %e ; %e\n",ffindex,ffvalue) ;
      
        MA38Method_PrintErrorDiagnostics(cntl,icntl,info) ;

        return(-1) ;
      }
    }
  
    /* Solve a * x = b */
    {
      double* b = Solver_GetRHS(solver) ;
      double* x = Solver_GetSolution(solver) ;
      GenericData_t* gw = Solver_GetGenericWorkSpace(solver) ;
      double* w = GenericData_FindData(gw,double,"work") ;
      int job = 0 ;
      bool transc = false ;
    
      ma38cd_(&n,&job,&transc,&lvalue,&lindex,value,index,keep,b,x,w,cntl,icntl,info,rinfo) ;
    
      if(info[0] < 0) {
      
        printf("MA38Method_Solve: error in the resolution\n") ;
      
        return(-1) ;
      }
    }
  }
  
  /* Broadcast to other processors */
  if(size > 1) {
    double* x = Solver_GetSolution(solver) ;
    
    #if DistributedMS_APIis(MPI)
      MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
    #else
      #error "Distributed memory system not available "
    #endif
  }
  
  return(0) ;
}


void MA38Method_PrintErrorDiagnostics(double* cntl,int* icntl,int* info)
{
  int ne     = info[4] + info[2] + info[1] ;
  int n      = info[11] - info[10] - info[9] - info[6] ;
  
  printf("Number of duplicates: %d\n",info[1]) ;
  
  printf("Number of entries dropped: %d\n",info[2]) ;
  
  printf("Number of entries in the matrix after removing duplicates and invalid entries: %d\n",info[4]) ;
  
  printf("Number of entries in total: %d\n",ne) ;
  printf("Dimension of the matrix: %d\n",n) ;
  
  printf("Value of info from MA38: %d\n",info[0]) ;
  
  
  if(info[0] == -3 || info[0] == -5) {
    printf("Insufficient integer workspace. Increase the length of INDEX , or decrease fill-in. LINDEX must be set to at least the value of INFO(19) to progress beyond the point of failure. However, note that setting LINDEX to INFO(19) does not guarantee a successful call. Fill-in can usually be decreased by increasing ICNTL(5) and decreasing CNTL(1) . Fill-in can sometimes be decreased by factorizing the transpose matrix then calling MA38C/CD with TRANSC=.TRUE. None of these suggestions for reducing fill-in are guaranteed to work, since the pivot search is a heuristic.\n") ;
  }
  if(info[0] == -4 || info[0] == -5) {
    printf("Insufficient real workspace. Increase the length of VALUE , or decrease fill-in. LVALUE must be set to at least the value of INFO(21) to progress beyond the point of failure. However, note that setting LVALUE to INFO(21) does not guarantee a successful call. See the comments for error –3.\n") ;
  }
  
  printf("Garbage collections are performed on both INDEX and VALUE when the available space in either array is exhausted.\n") ;
  
  printf("Number of garbage collections caused by insufficient space in INDEX: %d\n",info[13]) ;

  printf("Number of garbage collections caused by insufficient space in VALUE: %d\n",info[14]) ;
  
  printf("If one of these number is excessively high, performance can be degraded. Try increasing the lengths of INDEX or VALUE if that occurs. (or try reducing fill-in using the suggestions listed in Section 2.3.1 under error –3).\n") ;
      
  printf("Amount of space used in INDEX on a call to MA38: %d\n",info[17]) ;
  
  printf("Minimum length of INDEX for MA38AD to succeed: %d\n",info[18]) ;
  
  printf("However, many garbage collections might be required if the problem were rerun with LINDEX set to this value.\n") ;
      
  printf("Amount of space used in VALUE on a call to MA38: %d\n",info[19]) ;
  
  printf("Minimum length of VALUE for MA38AD to succeed: %d\n",info[20]) ;
  
  printf("However, many garbage collections might be required if the problem were rerun with LVALUE set to this value.\n") ;
  
  printf("Number of steps of iterative refinement performed: %d\n",info[23]) ;
  
  {
    double ffindex = ((double) info[18])/(3*ne + 2*n + 1) ;
    double ffvalue = ((double) info[20])/(2*ne) ;
    double ff = Math_Max(ffindex,ffvalue) ;
    
    printf("Set the fill factor to at least %lf\n",ff) ;
    
    printf("i.e. use the option:\n") ;
    
    printf("-solver ma38 -ff %lf\n",ff) ;
  }
  
}

#endif


