#include <stdio.h>
#include "Solutions.h"
#include "Message.h"



static Solution_t*  (Solution_Create)(Mesh_t*,int) ;
static void         (Solutions_Initialize)(Mesh_t*,Solutions_t*) ;


/* Extern functions */

Solutions_t*   (Solutions_Create)(Mesh_t* mesh,int n_sol)
{
  Solutions_t* sols = (Solutions_t*) malloc(sizeof(Solutions_t)) ;
  
  if(!sols) arret("Solutions_Create") ;
  
  Solutions_GetNbOfSolutions(sols) = n_sol ;
  Solutions_GetSolution(sols) = Solution_Create(mesh,n_sol) ;
  
  Solutions_Initialize(mesh,sols) ;
  
  return(sols) ;
}



void   (Solutions_Initialize)(Mesh_t* mesh,Solutions_t* sols)
/** Define properties of elements including:
  *  1. shape and interpolation functions
  *  2. allocation memory of internal data (abstract data type)
  *  3. size of tables for (im/ex)plicit and constant terms 
  *     for default allocation (see below) */
{
  
  {
    Mesh_InitializeSolutionPointers(mesh,sols) ;
  
    Elements_DefineProperties(Mesh_GetElements(mesh)) ;
  }
  
  /* Default allocation memory of internal data
   * i.e. for (im/ex)plicit and constant terms */
  {
    int n_sol = Solutions_GetNbOfSolutions(sols) ;
    Solution_t* sol = Solutions_GetSolution(sols) ;
    int i ;
    
    /* Allocation for (im/ex)plicit terms */
    {
      for(i = 0 ; i < n_sol ; i++) {
        ElementsSol_t* elementssol = Solution_GetElementsSol(sol + i) ;
      
        ElementsSol_AllocateMemoryForImplicitTerms(elementssol) ;
        ElementsSol_AllocateMemoryForExplicitTerms(elementssol) ;
      }
    }
    
    /* Allocation for constant terms */
    {
      ElementsSol_t* elementssol0 = Solution_GetElementsSol(sol) ;
      
      ElementsSol_AllocateMemoryForConstantTerms(elementssol0) ;
    
      /* All the other elementssol share the same constant terms */
      for(i = 1 ; i < n_sol ; i++) {
        ElementsSol_t* elementssol = Solution_GetElementsSol(sol + i) ;
      
        ElementsSol_ShareConstantTermsFrom(elementssol0,elementssol) ;
      }
    }
  }
}



void  (Solutions_MergeExplicitTerms)(Solutions_t* sols)
{
  int n_sol = Solutions_GetNbOfSolutions(sols) ;
  Solution_t* sol = Solutions_GetSolution(sols) ;
  int   NbOfElements = Solution_GetNbOfElements(sol) ;
  ElementSol_t* elemtsol_0 = Solution_GetElementSol(sol) ;
  int   i ;
  
  for(i = 1 ; i < n_sol ; i++) {
    ElementSol_t* elemtsol_i = Solution_GetElementSol(sol + i) ;
    int j ;
    
    free(ElementSol_GetExplicitTerm(elemtsol_i)) ;
  
    for(j = 0 ; j < NbOfElements ; j++) {
      void* ve0 = ElementSol_GetExplicitTerm(elemtsol_0 + j) ;
      
      ElementSol_GetExplicitTerm(elemtsol_i + j) = ve0 ;
    }
  }
  
}



void Solutions_StepForward(Solutions_t* sols)
/** Step forward in the loop */
{
  Solution_t* sol = Solutions_GetSolution(sols) ;
  Solution_t* next = Solution_GetNextSolution(sol) ;
  
  Solutions_GetSolution(sols) = next ;
}



void Solutions_StepBackward(Solutions_t* sols)
/** Step backward in the loop */
{
  Solution_t* sol = Solutions_GetSolution(sols) ;
  Solution_t* prev = Solution_GetPreviousSolution(sol) ;
  
  Solutions_GetSolution(sols) = prev ;
}



Solution_t*   (Solution_Create)(Mesh_t* mesh,int n_sol)
{
  Solution_t* sol = (Solution_t*) malloc(n_sol*sizeof(Solution_t)) ;
  int    i ;
  
  if(!sol) arret("Solution_Create") ;

  for(i = 0 ; i < n_sol ; i++) {
    Solution_t* soli = sol + i ;
    
    Solution_GetNodesSol(soli)    = NodesSol_Create(mesh) ;
    Solution_GetElementsSol(soli) = ElementsSol_Create(mesh) ;
    Solution_GetTime(soli)      = 0 ;
    Solution_GetTimeStep(soli)  = 0 ;
    Solution_GetStepIndex(soli) = 0 ;
    Solution_GetPreviousSolution(soli) = NULL ;
    Solution_GetNextSolution(soli)     = NULL ;
  }
  

  /* Initialized as a circularly linked list */
  {
    unsigned int n_no = Mesh_GetNbOfNodes(mesh) ;
    unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
    
    for(i = 0 ; i < n_sol ; i++) {
      int i_1 = i ;
      int i_n = (i) ? i - 1 : n_sol - 1 ;
    
      Solution_GetPreviousSolution(sol + i_1) = sol + i_n ;
      Solution_GetNextSolution(sol + i_n)     = sol + i_1 ;
      
      {
        NodeSol_t* nodesol_1 = Solution_GetNodeSol(sol + i_1) ;
        NodeSol_t* nodesol_n = Solution_GetNodeSol(sol + i_n) ;
        unsigned int j ;
        
        for(j = 0 ; j < n_no ; j++) {
          NodeSol_GetPreviousNodeSol(nodesol_1 + j) = nodesol_n + j ;
        }
      }

      {
        ElementSol_t* elementsol_1 = Solution_GetElementSol(sol + i_1) ;
        ElementSol_t* elementsol_n = Solution_GetElementSol(sol + i_n) ;
        unsigned int j ;
        
        for(j = 0 ; j < n_el ; j++) {
          ElementSol_GetPreviousElementSol(elementsol_1 + j) = elementsol_n + j ;
        }
      }
    }
  }
  
  return(sol) ;
}
