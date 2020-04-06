#include <stdio.h>
#include <assert.h>
#include "Solutions.h"
#include "Mry.h"
#include "Message.h"



static void   (Solutions_Initialize)(Solutions_t*) ;
static void   (Solutions_AllocateMemory)(Solutions_t*) ;


/* Extern functions */

Solutions_t*   (Solutions_Create)(Mesh_t* mesh,const int n_sol)
{
  Solutions_t* sols = (Solutions_t*) Mry_New(Solutions_t) ;
  
  Solutions_GetNbOfSolutions(sols) = n_sol ;
  
  {
    Solution_t* sol = (Solution_t*) Mry_New(Solution_t[n_sol]) ;
    int i ;
      
    for(i = 0 ; i < n_sol ; i++) {
      Solution_t* soli = Solution_Create(mesh) ;
      
      sol[i] = soli[0] ;
    }
      
    Solutions_GetSolution(sols) = sol ;
  }
  
  Solutions_Initialize(sols) ;
  
  
  /* Define properties of elements including:
   *  1. shape and interpolation functions
   *  2. allocation memory of internal data (abstract data type)
   *  3. size of tables for (im/ex)plicit and constant terms 
   *     for default allocation.
   *  4. eliminate possible degrees of freedom at some nodes
   *     (Element_GetUnknownPosition/Element_GetEquationPosition set to < 0)
   *  5. merge unknowns of slave nodes
   *     (Node_GetMatrixRowIndex/Node_GetMatrixColumnIndex set to < 0)
   */
  {
    Mesh_InitializeSolutionPointers(mesh,sols) ;
    Elements_DefineProperties(Mesh_GetElements(mesh)) ;
  }
  
  /* In case of 4. we re-set equation continuity ie re-initialize
   * Element_GetUnknownPosition and Element_GetEquationPosition */
  Mesh_SetEquationContinuity(mesh) ;
  
  /* In case of 5. we update the matrix rows/columns:
   * Node_GetMatrixRowIndex and Node_GetMatrixColumnIndex */
  Mesh_UpdateMatrixRowColumnIndexes(mesh) ;

  Solutions_AllocateMemory(sols) ;
  
  return(sols) ;
}



void   (Solutions_Initialize)(Solutions_t* sols)
/** Initialized as a circularly linked list */
{
  {
    int n_sol = Solutions_GetNbOfSolutions(sols) ;
    Solution_t* sol = Solutions_GetSolution(sols) ;
    unsigned int n_no = Solution_GetNbOfNodes(sol) ;
    unsigned int n_el = Solution_GetNbOfElements(sol) ;
    int    i ;
    
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
}



void   (Solutions_AllocateMemory)(Solutions_t* sols)
 /**  Allocation memory of internal data (abstract data type)
  *   The size of tables for (im/ex)plicit and constant terms 
  *   must have been defined through the call to 
  *   Elements_DefineProperties.
  */
{
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
