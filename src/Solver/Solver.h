#ifndef SOLVER_H
#define SOLVER_H


enum ResolMethod_e {          /* type de methode de resolution */
  CROUT,                      /* methode de Crout */
  SLU                         /* methode SuperLU */
} ;

/* class-like structures "Solver_t" and attributes */

/* vacuous declarations and typedef names */
struct Solver_s       ; typedef struct Solver_s       Solver_t ;
typedef enum ResolMethod_e  ResolMethod_t ;


#include "Mesh.h"
#include "Options.h"
#include "Matrix.h"

extern Solver_t*  Solver_Create(Mesh_t*,Options_t*) ;
extern void       Solver_Print(Solver_t*,char*) ;


#define Solver_GetResolutionMethod(solver)  ((solver)->mth)
#define Solver_GetNbOfRows(solver)          ((solver)->n)
#define Solver_GetNbOfColumns(solver)       ((solver)->n)
#define Solver_GetMatrix(solver)            ((solver)->a)
#define Solver_GetRHS(solver)               ((solver)->b)
#define Solver_GetSolution(solver)          ((solver)->x)

#define Solver_GetSolve(solver)             ((solver)->solve)

#define Solver_Solve(solver)                (Solver_GetSolve(solver)(solver))

/*  Typedef names of Methods */
typedef int  Solver_Solve_t(Solver_t*) ;


/* complete the structure types by using the typedef */
struct Solver_s {             /* System solver */
  Solver_Solve_t* solve ;
  ResolMethod_t mth ;         /* Method */
  unsigned int    n ;         /* Nb of rows/columns */
  Matrix_t* a ;               /* Matrix */
  double* b ;                 /* RHS */
  double* x ;                 /* Solution */
} ;

#endif
