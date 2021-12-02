#ifndef SNIA_H
#define SNIA_H


#include "Solutions.h"
#include "DataSet.h"
#include "Solver.h"


extern int (SNIA_Initialize)(DataSet_t*,Solutions_t*) ;
extern int (SNIA_Increment)(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
extern int (SNIA_StepForward)(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
extern int (SNIA_Iterate)(DataSet_t*,Solutions_t*,Solver_t*) ;


#endif
