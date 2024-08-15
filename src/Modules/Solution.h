#ifndef SOLUTION_H
#define SOLUTION_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "Solution_t" and attributes */

/* vacuous declarations and typedef names */

/* class-like structure "Solutions_t" */
struct Solution_s     ; typedef struct Solution_s     Solution_t ;



/* Declaration of Macros, Methods and Structures */

#include "Mesh.h"
 
extern Solution_t*  (Solution_Create)  (Mesh_t*) ;
extern void         (Solution_Delete)  (void*) ;
extern void         (Solution_Copy)    (Solution_t*,Solution_t*) ;
extern void         (Solution_InterpolateCurrentUnknowns)(Solution_t*,const int) ;
extern void         (Solution_CopySelectedSequentialUnknowns)(Solution_t*,Solution_t*,const int) ;
//extern Solution_t*  (Solution_GetSolutionInDistantFuture)(Solution_t*,unsigned int) ;
extern Solution_t*  (Solution_GetSolutionInDistantPast)(Solution_t*,unsigned int) ;


//#define Solution_GetTime(SOL)                 ((SOL)->t)
#define Solution_GetNbOfSequences(SOL)        ((SOL)->nbofsequences)
#define Solution_GetSequentialTime(SOL)       ((SOL)->sequentialtime)
#define Solution_GetSequentialTimeStep(SOL)   ((SOL)->sequentialtimestep)
#define Solution_GetSequentialStepIndex(SOL)  ((SOL)->sequentialstepindex)
#define Solution_GetNodesSol(SOL)             ((SOL)->nodessol)
#define Solution_GetElementsSol(SOL)          ((SOL)->elementssol)
#define Solution_GetPreviousSolution(SOL)     ((SOL)->sol_p)
#define Solution_GetNextSolution(SOL)         ((SOL)->sol_n)



/* The time, time step and step index */
#define Solution_GetTime(SOL) \
        Solution_GetSequentialTime(SOL)[Solution_GetNbOfSequences(SOL)-1]
        
#define Solution_GetTimeStep(SOL) \
        Solution_GetSequentialTimeStep(SOL)[Solution_GetNbOfSequences(SOL)-1]
        
#define Solution_GetStepIndex(SOL) \
        Solution_GetSequentialStepIndex(SOL)[Solution_GetNbOfSequences(SOL)-1]
        
#define Solution_GetPreviousTime(SOL) \
        Solution_GetSequentialTime(Solution_GetPreviousSolution(SOL))
        
#define Solution_GetPreviousTimeStep(SOL) \
        Solution_GetSequentialTimeStep(Solution_GetPreviousSolution(SOL))
        
#define Solution_GetPreviousStepIndex(SOL) \
        Solution_GetSequentialStepIndex(Solution_GetPreviousSolution(SOL))



/* Access to Nodes */
#define Solution_GetNodes(SOL) \
        NodesSol_GetNodes(Solution_GetNodesSol(SOL))


/* Access to node solutions */
#define Solution_GetNodeSol(SOL) \
        NodesSol_GetNodeSol(Solution_GetNodesSol(SOL))

#define Solution_GetNbOfNodes(SOL) \
        NodesSol_GetNbOfNodes(Solution_GetNodesSol(SOL))

//#define Solution_GetNbOfDOF(SOL) \
//        NodesSol_GetNbOfDOF(Solution_GetNodesSol(SOL))

//#define Solution_GetNodalValue(SOL) \
//        NodesSol_GetNodalValue(Solution_GetNodesSol(SOL))



/* Access to element solutions */
#define Solution_GetElementSol(SOL) \
        ElementsSol_GetElementSol(Solution_GetElementsSol(SOL))

#define Solution_GetNbOfElements(SOL) \
        ElementsSol_GetNbOfElements(Solution_GetElementsSol(SOL))



/* Copy the N solutions SOL+[0:N-1] in SOL+[1:N] */
/* #define Solution_CopyInDistantFuture(SOL,N) \
        do { \
          int dist = (N) - 1 ; \
          while(dist >= 0) { \
            Solution_t* sol_src  = Solution_GetSolutionInDistantFuture(SOL,dist) ; \
            Solution_t* sol_dest = Solution_GetNextSolution(sol_src) ; \
            Solution_Copy(sol_dest,sol_src) ; \
            dist-- ; \
          } \
        } while(0) */
        
        
#define Solution_InitializeSequentialTimes(SOL) \
        do { \
          int i ; \
          for(i = 0 ; i < Solution_GetNbOfSequences(SOL) ; i++) { \
            Solution_GetSequentialTime(SOL)[i] = Solution_GetTime(SOL) ; \
          } \
        } while(0)



/* Allocate memory */
#define Solution_AllocateMemoryForImplicitTerms(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_AllocateMemoryForImplicitTerms(elementssol) ; \
        } while(0)


#define Solution_AllocateMemoryForExplicitTerms(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_AllocateMemoryForExplicitTerms(elementssol) ; \
        } while(0)


#define Solution_AllocateMemoryForConstantTerms(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_AllocateMemoryForConstantTerms(elementssol) ; \
        } while(0)



/* Share data */
#define Solution_ShareConstantTerms(SOL_DEST,SOL_SRC) \
        do { \
          ElementsSol_t* elementssol_s = Solution_GetElementsSol(SOL_SRC) ; \
          ElementsSol_t* elementssol_d = Solution_GetElementsSol(SOL_DEST) ; \
          ElementsSol_ShareConstantGenericData(elementssol_d,elementssol_s) ; \
        } while(0)


#define Solution_ShareExplicitTerms(SOL_DEST,SOL_SRC) \
        do { \
          ElementsSol_t* elementssol_s = Solution_GetElementsSol(SOL_SRC) ; \
          ElementsSol_t* elementssol_d = Solution_GetElementsSol(SOL_DEST) ; \
          ElementsSol_ShareExplicitGenericData(elementssol_d,elementssol_s) ; \
        } while(0)



/* Delete data */
#define Solution_DeleteConstantTerms(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_DeleteConstantGenericData(elementssol) ; \
        } while(0)


#define Solution_DeleteExplicitTerms(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_DeleteExplicitGenericData(elementssol) ; \
        } while(0)



/* Set pointers to null */
#define Solution_DiscardExplicitGenericData(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_DiscardExplicitGenericData(elementssol) ; \
        } while(0)


#define Solution_DiscardConstantGenericData(SOL) \
        do { \
          ElementsSol_t* elementssol = Solution_GetElementsSol(SOL) ; \
          ElementsSol_DiscardConstantGenericData(elementssol) ; \
        } while(0)



#include "NodesSol.h"
#include "ElementsSol.h"

struct Solution_s {           /* solution at a time t */
  int   nbofsequences ;
  //double  t ;
  //double  dt ;
  //int     index ;
  double* sequentialtime ;
  double* sequentialtimestep ;
  int*    sequentialstepindex ;
  NodesSol_t* nodessol ;
  ElementsSol_t* elementssol ;
  Solution_t* sol_p ;         /* previous solution */
  Solution_t* sol_n ;         /* next solution */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
