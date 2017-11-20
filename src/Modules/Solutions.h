#ifndef SOLUTIONS_H
#define SOLUTIONS_H

/* class-like structures "Solution_t" and attributes */

/* vacuous declarations and typedef names */

/* class-like structure "Solutions_t" */
struct Solutions_s     ; typedef struct Solutions_s     Solutions_t ;



/* Declaration of Macros, Methods and Structures */

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


#include "Solution.h"

struct Solutions_s {              /* Solutions */
  unsigned int n_sol ;            /* Nb of solutions */
  Solution_t* solution ;          /* Solution */
} ;
 
 
 


#endif
