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
  
    NodeSol_GetPreviousNodeSol(nodesol) = NULL ;
    NodeSol_GetNextNodeSol(nodesol)     = NULL ;
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



NodeSol_t* (NodeSol_GetNodeSolInDistantPast)(NodeSol_t* nodesol,unsigned int dist)
{
  while(dist--) nodesol = NodeSol_GetPreviousNodeSol(nodesol) ;
  return(nodesol) ;
}



NodeSol_t* (NodeSol_GetNodeSolInDistantFuture)(NodeSol_t* nodesol,unsigned int dist)
{
  while(dist--) nodesol = NodeSol_GetNextNodeSol(nodesol) ;
  return(nodesol) ;
}
