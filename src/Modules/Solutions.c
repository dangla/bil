#include <stdio.h>
#include <assert.h>
#include "Mesh.h"
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
  Solutions_GetMergeIndex(sols) = 0 ;
  
  {
    Solution_t* sol = (Solution_t*) Mry_New(Solution_t[n_sol]) ;
    int i ;
      
    for(i = 0 ; i < n_sol ; i++) {
      Solution_t* soli = Solution_Create(mesh) ;
      
      sol[i] = soli[0] ;
      free(soli) ;
    }
      
    Solutions_GetSolution(sols) = sol ;
    Solutions_GetSolutionHead(sols) = sol ;
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
    Solutions_InitializeMeshPointers(sols,mesh) ;
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



void (Solutions_Delete)(void* self)
{
  Solutions_t* sols = (Solutions_t*) self ;
  
  {
    int n_sol = Solutions_GetNbOfSolutions(sols) ;
    Solution_t* sol = Solutions_GetSolutionHead(sols) ;
    
    if(sol) {
      int i ;
    
      for(i = 0 ; i < n_sol ; i++) {
        Solution_t* soli = sol + i ;
        
        /* Since the constant terms are shared we delete them only once.
         * So for i > 0 we set the pointers to null */
        if(i > 0) {
          Solution_DiscardConstantGenericData(soli) ;
          if(Solutions_ExplicitTermsAreMerged(sols)) {
            Solution_DiscardExplicitGenericData(soli) ;
          }
        }
      
        Solution_Delete(soli) ;
      }
      
      free(sol) ;
    }
  }
}



void   (Solutions_Initialize)(Solutions_t* sols)
/** Initialized as a circularly linked list */
{
  {
    int n_sol = Solutions_GetNbOfSolutions(sols) ;
    Solution_t* sol = Solutions_GetSolution(sols) ;
    int    i ;
    
    for(i = 0 ; i < n_sol ; i++) {
      int i_1 = i ;
      int i_n = (i) ? i - 1 : n_sol - 1 ;
    
      Solution_GetPreviousSolution(sol + i_1) = sol + i_n ;
      Solution_GetNextSolution(sol + i_n)     = sol + i_1 ;
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
  int n_sol = Solutions_GetNbOfSolutions(sols) ;
  Solution_t* sol = Solutions_GetSolution(sols) ;
  
  /* Default allocation memory of internal data
   * i.e. for (im/ex)plicit and constant terms */
  {
    int i ;

    /* Allocation for (im/ex)plicit terms */
    for(i = 0 ; i < n_sol ; i++) {
      Solution_AllocateMemoryForImplicitTerms(sol + i) ;
      Solution_AllocateMemoryForExplicitTerms(sol + i) ;
    }
  }

    /* Allocation for constant terms */
  {
    Solution_AllocateMemoryForConstantTerms(sol) ;
    Solutions_MergeConstantTerms(sols) ;
  }
}



void  (Solutions_MergeConstantTerms)(Solutions_t* sols)
/** Share the constant terms between all the solutions.
 */
{
  int n_sol = Solutions_GetNbOfSolutions(sols) ;
  Solution_t* sol = Solutions_GetSolution(sols) ;
  int   i ;
  
  for(i = 1 ; i < n_sol ; i++) {
    Solution_DeleteConstantTerms(sol + i) ;
    Solution_ShareConstantTerms(sol + i,sol) ;
  }
}



void  (Solutions_MergeExplicitTerms)(Solutions_t* sols)
/** Share the explicit terms between all the solutions.
 */
{
  int n_sol = Solutions_GetNbOfSolutions(sols) ;
  Solution_t* sol = Solutions_GetSolution(sols) ;
  int   i ;
  
  Solutions_SetIndexForMergedExplicitTerms(sols) ;
  
  for(i = 1 ; i < n_sol ; i++) {
    Solution_DeleteExplicitTerms(sol + i) ;
    Solution_ShareExplicitTerms(sol + i,sol) ;
  }
}



void (Solutions_StepForward)(Solutions_t* sols)
/** Step forward in the loop */
{
  Solution_t* sol = Solutions_GetSolution(sols) ;
  Solution_t* next = Solution_GetNextSolution(sol) ;
  
  Solutions_GetSolution(sols) = next ;
}



void (Solutions_StepBackward)(Solutions_t* sols)
/** Step backward in the loop */
{
  Solution_t* sol = Solutions_GetSolution(sols) ;
  Solution_t* prev = Solution_GetPreviousSolution(sol) ;
  
  Solutions_GetSolution(sols) = prev ;
}



void (Solutions_InitializeMeshPointers)(Solutions_t* sols,Mesh_t* mesh)
/** Initialize the pointers of nodes and elements to sols. */
{
  int n_no = Mesh_GetNbOfNodes(mesh) ;
  Node_t* no = Mesh_GetNode(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int    i ;

  for(i = 0 ; i < n_el ; i++) {
    Element_GetSolutions(el + i) = sols ;
  }
  
  for(i = 0 ; i < n_no ; i++) {
    Node_GetSolutions(no + i) = sols ;
  }

}
