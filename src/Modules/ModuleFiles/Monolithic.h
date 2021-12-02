#ifndef MONOLITHIC_H
#define MONOLITHIC_H


#include "Solutions.h"
#include "DataSet.h"
#include "Solvers.h"


extern int (Monolithic_Initialize)(DataSet_t*,Solutions_t*) ;
extern int (Monolithic_Increment)(DataSet_t*,Solutions_t*,Solver_t*,OutputFiles_t*,double,double) ;
extern int (Monolithic_StepForward)(DataSet_t*,Solutions_t*,Solver_t*,double,double) ;
extern int (Monolithic_Iterate)(DataSet_t*,Solutions_t*,Solver_t*) ;

#endif
