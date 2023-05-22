#ifndef SOLUTIONS_H
#define SOLUTIONS_H

/* class-like structures "Solution_t" and attributes */

/* vacuous declarations and typedef names */

/* class-like structure "Solutions_t" */
struct Solutions_s     ; typedef struct Solutions_s     Solutions_t ;



/* Declaration of Macros, Methods and Structures */

#include "Mesh.h"
 
extern Solutions_t*   (Solutions_Create)(Mesh_t*,const int) ;
extern void           (Solutions_Delete)(void*) ;
extern void           (Solutions_MergeExplicitTerms)(Solutions_t*) ;
extern void           (Solutions_MergeConstantTerms)(Solutions_t*) ;
extern void           (Solutions_StepForward)(Solutions_t*) ;
extern void           (Solutions_StepBackward)(Solutions_t*) ;
extern void           (Solutions_InitializeMeshPointers)(Solutions_t*,Mesh_t*) ;


#define Solutions_GetNbOfSolutions(SOLS)    ((SOLS)->n_sol)
#define Solutions_GetSolution(SOLS)         ((SOLS)->solution)
#define Solutions_GetSolutionHead(SOLS)     ((SOLS)->head)
#define Solutions_GetMergeIndex(SOLS)       ((SOLS)->mergeindex)



/* Access to solutions */
#define Solutions_GetCurrentSolution(SOLS) \
        Solutions_GetSolution(SOLS)

#define Solutions_GetPreviousSolution(SOLS) \
        Solution_GetPreviousSolution(Solutions_GetSolution(SOLS))

#define Solutions_GetNextSolution(SOLS) \
        Solution_GetNextSolution(Solutions_GetSolution(SOLS))

#define Solutions_GetSolutionInDistantPast(SOLS,i) \
        Solution_GetSolutionInDistantPast(Solutions_GetSolution(SOLS),i)


/* Access to times */
#define Solutions_GetTime(SOLS) \
        Solution_GetTime(Solutions_GetSolution(SOLS))
        
#define Solutions_GetCurrentTime(SOLS) \
        Solution_GetTime(Solutions_GetCurrentSolution(SOLS))
        
#define Solutions_GetPreviousTime(SOLS) \
        Solution_GetTime(Solutions_GetPreviousSolution(SOLS))



/* Step forward ; Step backward */
#if 1
#define Solutions_StepForward(SOLS) \
        do { \
          Solutions_GetSolution(SOLS) = Solution_GetNextSolution(Solutions_GetSolution(SOLS)) ; \
        } while(0)

#define Solutions_StepBackward(SOLS) \
        do { \
          Solutions_GetSolution(SOLS) = Solution_GetPreviousSolution(Solutions_GetSolution(SOLS)) ; \
        } while(0)
#endif



/* On merging explicit terms */
#define Solutions_SetIndexForMergedExplicitTerms(SOLS) \
        do { \
          Solutions_GetMergeIndex(SOLS) = 1 ; \
        } while(0)
        
#define Solutions_SetIndexForNotMergedExplicitTerms(SOLS) \
        do { \
          Solutions_GetMergeIndex(SOLS) = 0 ; \
        } while(0)

#define Solutions_ExplicitTermsAreMerged(SOLS) \
        (Solutions_GetMergeIndex(SOLS) == 1)



#include "Solution.h"

struct Solutions_s {              /* Solutions */
  unsigned int n_sol ;            /* Nb of solutions */
  char mergeindex ;
  Solution_t* head ;              /* Solution head */
  Solution_t* solution ;          /* Solution */
} ;
 
 
 


#endif
