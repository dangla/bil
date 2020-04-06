#ifdef BLASLIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "Solver.h"
#include "Message.h"
#include "Matrix.h"


#include "MA38Method.h"
#include "CoordinateFormat.h"




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
  
  /* Initialize controls */
  ma38id_(keep,cntl,icntl) ;
  
  /* Printing controls:
   * 1: Print only error messages, 
   * 2: + warnings, 
   * 3: + terse diagnostics 
   */
  {
    Options_t* options = CoordinateFormat_GetOptions(ac) ;
    
    icntl[2] = 1 ;
    
    if(Options_IsToPrintOutAtEachIteration(options)) {
      icntl[2] = 2 ;
    }
  }
  
  /* Factorize the matrix */
  {
    int ne = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
    int job = 0 ;
    bool transa = false ;
  
    ma38ad_(&n,&ne,&job,&transa,&lvalue,&lindex,value,index,keep,cntl,icntl,info,rinfo) ;
    
    if(info[0] < 0) {
      return(-1) ;
    }
  }
  
  /* Solve a * x = b */
  {
    double* b = Solver_GetRHS(solver) ;
    double* x = Solver_GetSolution(solver) ;
    double* w = (double*) Matrix_GetWorkSpace(a) ;
    int job = 0 ;
    bool transc = false ;
    
    ma38cd_(&n,&job,&transc,&lvalue,&lindex,value,index,keep,b,x,w,cntl,icntl,info,rinfo) ;
    
    if(info[0] < 0) {
      return(-1) ;
    }
  }
  
  return(0) ;
}

#endif


