#ifndef SOLVER_H
#define SOLVER_H

/* class-like structures "Solver_t" and attributes */

/* vacuous declarations and typedef names */
struct Solver_s       ; typedef struct Solver_s       Solver_t ;


#include "Mesh.h"
#include "Options.h"

extern Solver_t*  (Solver_Create)(Mesh_t*,Options_t*,const int) ;
extern void       (Solver_Delete)(void*) ;
extern void       (Solver_Print)(Solver_t*,char*) ;


#define Solver_GetResolutionMethod(SV)  ((SV)->mth)
#define Solver_GetNbOfRows(SV)          ((SV)->n)
#define Solver_GetNbOfColumns(SV)       ((SV)->n)
#define Solver_GetMatrix(SV)            ((SV)->a)
#define Solver_GetRHS(SV)               ((SV)->b)
#define Solver_GetSolution(SV)          ((SV)->x)
#define Solver_GetSolve(SV)             ((SV)->solve)




#define Solver_Solve(SV) \
        (Solver_GetSolve(SV)(SV))



/*  Typedef names of Methods */
typedef int  Solver_Solve_t(Solver_t*) ;



#include "ResolutionMethod.h"
#include "Matrix.h"

/* complete the structure types by using the typedef */
struct Solver_s {             /* System solver */
  Solver_Solve_t* solve ;
  ResolutionMethod_t mth ;    /* Method */
  unsigned int    n ;         /* Nb of rows/columns */
  Matrix_t* a ;               /* Matrix */
  double* b ;                 /* RHS */
  double* x ;                 /* Solution */
} ;

#endif
