#include <stdio.h>
#include "Mry.h"
#include "NodeSol.h"
#include "Message.h"


/* Extern functions */


NodeSol_t* (NodeSol_Create)(const int n)
{
  NodeSol_t* nodesol = (NodeSol_t*) Mry_New(NodeSol_t) ;
  
  {
    double* u = (double*) Mry_New(double[n]) ;
  
    NodeSol_GetNbOfUnknowns(nodesol) = n ;
  
    NodeSol_GetUnknown(nodesol) = u ;
  
    //NodeSol_GetPreviousNodeSol(nodesol) = NULL ;
    //NodeSol_GetNextNodeSol(nodesol)     = NULL ;
  }
  
  return(nodesol) ;
}



void (NodeSol_Delete)(void* self)
{
  NodeSol_t* nodesol = (NodeSol_t*) self ;
  
  {
    double* u = NodeSol_GetUnknown(nodesol) ;
    
    if(u) {
      free(u) ;
      NodeSol_GetUnknown(nodesol) = NULL ;
    }
  }
}



void (NodeSol_Copy)(NodeSol_t* nodesol_d,NodeSol_t* nodesol_s)
/** Copy the nodal unknowns from nodesol_src to nodesol_dest */
{
  /* Nodal values */
  {
    double* u_s = NodeSol_GetUnknown(nodesol_s) ;
    double* u_d = NodeSol_GetUnknown(nodesol_d) ;
        
    if(u_d != u_s) {
      unsigned int nu = NodeSol_GetNbOfUnknowns(nodesol_s) ;
      unsigned int i ;
      
      for(i = 0 ; i < nu ; i++) {
        u_d[i] = u_s[i] ;
      }
    }
  }
}
