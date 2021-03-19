#include <stdio.h>
#include "Solution.h"
#include "Message.h"
#include "Mry.h"




/* Extern functions */



Solution_t*   (Solution_Create)(Mesh_t* mesh)
{
  Solution_t* sol = (Solution_t*) Mry_New(Solution_t) ;
  
  {
      Solution_GetNodesSol(sol)    = NodesSol_Create(mesh) ;
      Solution_GetElementsSol(sol) = ElementsSol_Create(mesh) ;
      Solution_GetTime(sol)      = 0 ;
      Solution_GetTimeStep(sol)  = 0 ;
      Solution_GetStepIndex(sol) = 0 ;
      Solution_GetPreviousSolution(sol) = NULL ;
      Solution_GetNextSolution(sol)     = NULL ;
  }
  
  return(sol) ;
}



void Solution_Copy(Solution_t* sol_dest,Solution_t* sol_src)
/** Copy unknowns, (im/ex)plicit and constant terms 
 *  from sol_src to sol_dest */
{
  Solution_GetTime(sol_dest) = Solution_GetTime(sol_src) ;
  Solution_GetTimeStep(sol_dest) = Solution_GetTimeStep(sol_src) ;
  Solution_GetStepIndex(sol_dest) = Solution_GetStepIndex(sol_src) ;
  
  /* Nodal values */
  {
    NodesSol_t* nodessol_d = Solution_GetNodesSol(sol_dest) ;
    NodesSol_t* nodessol_s = Solution_GetNodesSol(sol_src) ;
    
    NodesSol_Copy(nodessol_d,nodessol_s) ;
  }

  /* (Im/Ex)plicit and constant terms */
  {
    ElementsSol_t* elementssol_d = Solution_GetElementsSol(sol_dest) ;
    ElementsSol_t* elementssol_s = Solution_GetElementsSol(sol_src) ;
    
    ElementsSol_Copy(elementssol_d,elementssol_s) ;
  }
}



/* Not used */
#if 0
void Solution_Move(Solution_t* dest,Solution_t* src)
{
  Solution_t* prev = Solution_GetPreviousSolution(src) ;
  Solution_t* next = Solution_GetNextSolution(src) ;
          
  *dest = *src ;
          
  Solution_GetPreviousSolution(next) = dest ;
  Solution_GetNextSolution(prev) = dest ;
}



void Solution_CopyPreviousToCurrentValues(Solution_t* sol)
{
  Solution_t* sol_p = Solution_GetPreviousSolution(sol) ;
  unsigned int n_u  = Solution_GetNbOfDOF(sol) ;
  unsigned int n_vi = Solution_GetNbOfImplicitTerms(sol) ;
  unsigned int n_ve = Solution_GetNbOfExplicitTerms(sol) ;
  double* u_1  = Solution_GetNodalValue(sol) ;
  double* vi_1 = Solution_GetImplicitTerm(sol) ;
  double* ve_1 = Solution_GetExplicitTerm(sol) ;
  double* u_n  = Solution_GetNodalValue(sol_p) ;
  double* vi_n = Solution_GetImplicitTerm(sol_p) ;
  double* ve_n = Solution_GetExplicitTerm(sol_p) ;
  unsigned int    i ;
  
  Solution_GetTime(sol) = Solution_GetTime(sol_p) ;
  Solution_GetTimeStep(sol) = Solution_GetTimeStep(sol_p) ;
  Solution_GetStepIndex(sol) = Solution_GetStepIndex(sol_p) ;

  for(i = 0 ; i < n_ve ; i++) {
    ve_1[i] = ve_n[i] ;
  }
    
  for(i = 0 ; i < n_vi ; i++) {
    vi_1[i] = vi_n[i] ;
  }
  
  for(i = 0 ; i < n_u ; i++) {
    u_1[i] = u_n[i] ;
  }
}
#endif
