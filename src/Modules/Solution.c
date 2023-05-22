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
  }
  
  Solution_GetNbOfSequences(sol) = Mesh_GetNbOfMatrices(mesh) ;
  
  /* Allocation of space for the times and the time steps */
  {
    int n = Solution_GetNbOfSequences(sol) ;
    double* t = (double*) Mry_New(double[2*n]) ;
    double* dt = t + n ;

    Solution_GetSequentialTime(sol)      = t ;
    Solution_GetSequentialTimeStep(sol)  = dt ;
  }
  
  /* Allocation of space for the step indexes */
  {
    int n = Solution_GetNbOfSequences(sol) ;
    int* index = (int*) Mry_New(int[n]) ;
    
    Solution_GetSequentialStepIndex(sol) = index ;
  }
  
  {
    Solution_GetPreviousSolution(sol) = NULL ;
    Solution_GetNextSolution(sol)     = NULL ;
  }
  
  return(sol) ;
}




void   (Solution_Delete)(void* self)
{
  Solution_t* sol = (Solution_t*) self ;
  
  {
    NodesSol_t* nodessol = Solution_GetNodesSol(sol) ;
    
    if(nodessol) {
      NodesSol_Delete(nodessol) ;
      free(nodessol) ;
      Solution_GetNodesSol(sol) = NULL ;
    }
  }
  
  {
    ElementsSol_t* elementssol = Solution_GetElementsSol(sol) ;
    
    if(elementssol) {
      ElementsSol_Delete(elementssol) ;
      free(elementssol) ;
      Solution_GetElementsSol(sol) = NULL ;
    }
  }

  {
    double* t = Solution_GetSequentialTime(sol) ;
    
    if(t) {
      free(t) ;
      Solution_GetSequentialTime(sol) = NULL ;
    }
  }
  
  {
    int* index = Solution_GetSequentialStepIndex(sol) ;
    
    if(index) {
      free(index) ;
      Solution_GetSequentialStepIndex(sol) = NULL ;
    }
  }
}



void (Solution_Copy)(Solution_t* sol_dest,Solution_t* sol_src)
/** Copy unknowns, (im/ex)plicit and constant terms 
 *  from sol_src to sol_dest */
{
  Solution_GetTime(sol_dest) = Solution_GetTime(sol_src) ;
  
  /* Times, time steps and step indexes */
  {
    int n = Solution_GetNbOfSequences(sol_src) ;
    double* t_src  = Solution_GetSequentialTime(sol_src) ;
    double* dt_src = Solution_GetSequentialTimeStep(sol_src) ;
    int* index_src = Solution_GetSequentialStepIndex(sol_src) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Solution_GetSequentialTime(sol_dest)[i] = t_src[i] ;
      Solution_GetSequentialTimeStep(sol_dest)[i] = dt_src[i] ;
      Solution_GetSequentialStepIndex(sol_dest)[i] = index_src[i] ;
    }
  }
  
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



void (Solution_CopySelectedSequentialUnknowns)(Solution_t* sol_dest,Solution_t* sol_src,const int sequentialindex)
/** Copy unknowns of sequential index "sequentialindex"
 *  from sol_src to sol_dest */
{
  /* Time, time step and step index */
  {
    double* t_src  = Solution_GetSequentialTime(sol_src) ;
    double* dt_src = Solution_GetSequentialTimeStep(sol_src) ;
    int* index_src = Solution_GetSequentialStepIndex(sol_src) ;
    
    {
      int i = sequentialindex ;
      
      Solution_GetSequentialTime(sol_dest)[i] = t_src[i] ;
      Solution_GetSequentialTimeStep(sol_dest)[i] = dt_src[i] ;
      Solution_GetSequentialStepIndex(sol_dest)[i] = index_src[i] ;
    }
  }
  
  /* Nodal values */
  {
    NodesSol_t* nodessol_d = Solution_GetNodesSol(sol_dest) ;
    NodesSol_t* nodessol_s = Solution_GetNodesSol(sol_src) ;
    
    {
      int i = sequentialindex ;
      
      NodesSol_CopySelectedSequentialUnknowns(nodessol_d,nodessol_s,i) ;
    }
  }
}


#if 1
Solution_t* (Solution_GetSolutionInDistantPast)(Solution_t* sol,unsigned int dist)
{
  while(dist--) sol = Solution_GetPreviousSolution(sol) ;
  return(sol) ;
}
#endif


#if 0
Solution_t* (Solution_GetSolutionInDistantFuture)(Solution_t* sol,unsigned int dist)
{
  while(dist--) sol = Solution_GetNextSolution(sol) ;
  return(sol) ;
}
#endif





#if 1
void (Solution_InterpolateCurrentUnknowns)(Solution_t* sol_1,const int sequentialindex)
/** Interpolate the current nodal values of unknowns for which
 *  the sequential indexes are lower than "sequentialindex", the
 *  latter being the sequential index of the unknowns which are
 *  being solved. */
{
  
  if(sequentialindex <= 0) return ;
  
  {
    Solution_t* sol_n = Solution_GetPreviousSolution(sol_1) ;
    Solution_t* sol_2 = Solution_GetNextSolution(sol_1) ;
    double* t_n = Solution_GetSequentialTime(sol_n) ;
    double* t_1 = Solution_GetSequentialTime(sol_1) ;
    double* t_2 = Solution_GetSequentialTime(sol_2) ;
    
    /* Set the sequential time of the interpolated unknowns to the
     * sequential time of the unknowns which are being solved, i.e.,
     * those of index "sequentialindex". */
    {
      int k ;

      for(k = 0 ; k < sequentialindex ; k++) {
        t_1[k] = t_1[sequentialindex] ;
      }
    }

    /* Interpolation of the specified unknowns */
    {
      Nodes_t* nodes = Solution_GetNodes(sol_1) ;
      unsigned int nb_nodes = Nodes_GetNbOfNodes(nodes) ;
      Node_t* node = Nodes_GetNode(nodes) ;
      int   i ;
  
      for(i = 0 ; i < nb_nodes ; i++) {
        Node_t* nodi = node + i ;
        int*    node_seq_ind = Node_GetSequentialIndexOfUnknown(nodi) ;
        int  nb_unk = Node_GetNbOfUnknowns(nodi) ;
        double* u_n = Node_GetPreviousUnknown(nodi) ;
        double* u_1 = Node_GetCurrentUnknown(nodi) ;
        double* u_2 = Node_GetNextUnknown(nodi) ;
        
        {
          int j ;

          for(j = 0 ; j < nb_unk ; j++) {
            int k = node_seq_ind[j] ;
          
            if(k < sequentialindex) {
              double dt_1 = t_1[k] - t_n[k] ;
              double dt_2 = t_2[k] - t_n[k] ;
              double a = (dt_2 > 0) ? dt_1 / dt_2 : 1 ;
              
              u_1[j] = u_n[j] + (u_2[j] - u_n[j]) * a ;
            }
          }
        }
      }
    }
  }
}
#endif
