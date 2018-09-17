#ifndef MA38METHOD_H
#define MA38METHOD_H


#include "Solver.h"

extern int        MA38Method_Solve(Solver_t*) ;


#if defined(__cplusplus)
  extern "C" {
#endif

extern void   ma38id_(int*,double*,int*) ;
extern void   ma38ad_(int*,int*,int*,bool*,int*,int*,double*,int*,int*,double*,int*,int*,double*) ;
extern void   ma38cd_(int*,int*,bool*,int*,int*,double*,int*,int*,double*,double*,double*,double*,int*,int*,double*) ;

#if defined(__cplusplus)
  }
#endif

#endif
