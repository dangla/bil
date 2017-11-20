#ifndef SOLUTION_H
#define SOLUTION_H

/* class-like structures "Solution_t" and attributes */

/* vacuous declarations and typedef names */

/* class-like structure "Solutions_t" */
struct Solution_s     ; typedef struct Solution_s     Solution_t ;



/* Declaration of Macros, Methods and Structures */
 
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
