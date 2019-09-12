#ifndef ITERPROCESS_H
#define ITERPROCESS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct IterProcess_s  ; typedef struct IterProcess_s  IterProcess_t ;


#include "DataFile.h"
#include "ObVals.h"
#include "Nodes.h"
#include "Solver.h"

extern IterProcess_t*  IterProcess_Create(DataFile_t*,ObVals_t*) ;
extern int             IterProcess_SetCurrentError(IterProcess_t*,Nodes_t*,Solver_t*) ;
extern void            IterProcess_PrintCurrentError(IterProcess_t*) ;


#define IterProcess_GetNbOfIterations(IPR)           ((IPR)->niter)
#define IterProcess_GetNbOfRepetitions(IPR)          ((IPR)->nrecom)
#define IterProcess_GetTolerance(IPR)                ((IPR)->tol)
#define IterProcess_GetRepetitionIndex(IPR)          ((IPR)->irecom)
#define IterProcess_GetIterationIndex(IPR)           ((IPR)->iter)
#define IterProcess_GetCurrentError(IPR)             ((IPR)->error)
#define IterProcess_GetObValIndexOfCurrentError(IPR) ((IPR)->obvalindex)
#define IterProcess_GetNodeIndexOfCurrentError(IPR)  ((IPR)->nodeindex)
#define IterProcess_GetObVals(IPR)                   ((IPR)->obvals)



#define IterProcess_GetObVal(IPR) \
        ObVals_GetObVal(IterProcess_GetObVals(IPR))
        

/* Operations on iterations */
#define IterProcess_IncrementIterationIndex(IPR) \
        (IterProcess_GetIterationIndex(IPR)++)

#define IterProcess_LastIterationIsNotReached(IPR) \
        (IterProcess_GetIterationIndex(IPR) < IterProcess_GetNbOfIterations(IPR))

#define IterProcess_InitializeIterations(IPR) \
        (IterProcess_GetIterationIndex(IPR) = 0)


/* Operations on repetitions */
#define IterProcess_IncrementRepetitionIndex(IPR) \
       (IterProcess_GetRepetitionIndex(IPR)++)

#define IterProcess_LastRepetitionIsNotReached(IPR) \
        (IterProcess_GetRepetitionIndex(IPR) < IterProcess_GetNbOfRepetitions(IPR))

#define IterProcess_InitializeRepetitions(IPR) \
        (IterProcess_GetRepetitionIndex(IPR) = 0)


/* Operations on convergence */
#define IterProcess_ConvergenceIsMet(IPR) \
        (IterProcess_GetCurrentError(IPR) < IterProcess_GetTolerance(IPR))

#define IterProcess_ConvergenceIsNotMet(IPR) \
        (!IterProcess_ConvergenceIsMet(IPR))


/* Error on which unknown? */
#define IterProcess_GetNameOfTheCurrentError(IPR) \
        (ObVal_GetNameOfUnknown(IterProcess_GetObVal(IPR) + IterProcess_GetObValIndexOfCurrentError(IPR)))


struct IterProcess_s {        /* Iterative process */
  int    niter ;              /* Max nb of iterations */
  int    iter ;               /* Current iteration index */
  int    nrecom ;             /* Max nb of repetitions */
  int    irecom ;             /* Current repetition index */
  double tol ;                /* Tolerance */
  double error ;              /* Current error */
  int    obvalindex ;         /* Objective value index pertaining to the greatest error */
  int    nodeindex ;          /* Node index pertaining to the greatest error */
  ObVals_t* obvals ;          /* Objective variations */
} ;

#endif
