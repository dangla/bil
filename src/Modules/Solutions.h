#ifndef SOLUTIONS_H
#define SOLUTIONS_H

/* class-like structures "Solution_t" and attributes */

/* vacuous declarations and typedef names */

/* class-like structure "Solutions_t" */
struct Solutions_s     ; typedef struct Solutions_s     Solutions_t ;
/*   1. Solutions_t attributes */
struct Solution_s     ; typedef struct Solution_s     Solution_t ;



/* Declaration of Macros, Methods and Structures */

/* 1. Solutions_t
 * --------------*/
#include "Mesh.h"
 
extern Solutions_t*   (Solutions_Create)(Mesh_t*,int) ;
extern void           (Solutions_MergeExplicitTerms)(Solutions_t*) ;
extern void           (Solutions_StepForward)(Solutions_t*) ;
extern void           (Solutions_StepBackward)(Solutions_t*) ;


#define Solutions_GetNbOfSolutions(SOLS)    ((SOLS)->n_sol)
#define Solutions_GetSolution(SOLS)         ((SOLS)->solution)



/* Access to solutions */
#define Solutions_GetCurrentSolution(SOLS) \
        Solutions_GetSolution(SOLS)

#define Solutions_GetPreviousSolution(SOLS) \
        Solution_GetPreviousSolution(Solutions_GetSolution(SOLS))


/* Access to times */
#define Solutions_GetTime(SOLS) \
        Solution_GetTime(Solutions_GetSolution(SOLS))
        
#define Solutions_GetCurrentTime(SOLS) \
        Solution_GetTime(Solutions_GetCurrentSolution(SOLS))
        
#define Solutions_GetPreviousTime(SOLS) \
        Solution_GetTime(Solutions_GetPreviousSolution(SOLS))



struct Solutions_s {              /* Solutions */
  unsigned int n_sol ;            /* Nb of solutions */
  Solution_t* solution ;          /* Solution */
} ;
 
 
 

/* 2. Solution_t
 * -------------*/
 
extern void (Solution_Copy)(Solution_t*,Solution_t*) ;


#define Solution_GetTime(SOL)                 ((SOL)->t)
#define Solution_GetTimeStep(SOL)             ((SOL)->dt)
#define Solution_GetStepIndex(SOL)            ((SOL)->step)
#define Solution_GetNodesSol(SOL)             ((SOL)->nodessol)
#define Solution_GetElementsSol(SOL)          ((SOL)->elementssol)
#define Solution_GetPreviousSolution(SOL)     ((SOL)->sol_p)
#define Solution_GetNextSolution(SOL)         ((SOL)->sol_n)



/* Access to node solutions */
#define Solution_GetNodeSol(SOL) \
        NodesSol_GetNodeSol(Solution_GetNodesSol(SOL))

#define Solution_GetNbOfNodes(SOL) \
        NodesSol_GetNbOfNodes(Solution_GetNodesSol(SOL))

#define Solution_GetNbOfDOF(SOL) \
        NodesSol_GetNbOfDOF(Solution_GetNodesSol(SOL))

#define Solution_GetDOF(SOL) \
        NodesSol_GetDOF(Solution_GetNodesSol(SOL))



/* Access to element solutions */
#define Solution_GetElementSol(SOL) \
        ElementsSol_GetElementSol(Solution_GetElementsSol(SOL))

#define Solution_GetNbOfElements(SOL) \
        ElementsSol_GetNbOfElements(Solution_GetElementsSol(SOL))




#include "NodesSol.h"
#include "ElementsSol.h"

struct Solution_s {           /* solution at a time t */
  double t ;
  double dt ;
  int    step ;
  NodesSol_t* nodessol ;
  ElementsSol_t* elementssol ;
  Solution_t* sol_p ;         /* previous solution */
  Solution_t* sol_n ;         /* next solution */
} ;

#endif
