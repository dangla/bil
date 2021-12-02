#ifndef SOLVERS_H
#define SOLVERS_H

/* class-like structures "Solver_t" and attributes */

/* vacuous declarations and typedef names */
struct Solvers_s       ; typedef struct Solvers_s       Solvers_t ;


#include "Mesh.h"
#include "Options.h"

extern Solvers_t*  (Solvers_Create)(Mesh_t*,Options_t*,const int) ;
extern void        (Solvers_Delete)(void*) ;


#define Solvers_GetNbOfSolvers(SV)        ((SV)->nbofsolvers)
#define Solvers_GetSolver(SV)             ((SV)->solver)



#include "Solver.h"

/* complete the structure types by using the typedef */
struct Solvers_s {
  int nbofsolvers ;
  Solver_t* solver ;
} ;

#endif
