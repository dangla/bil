#include <stdio.h>
#include "NodesSol.h"
#include "Message.h"


/* Extern functions */

NodesSol_t* (NodesSol_Create)(Mesh_t* mesh)
{
  NodesSol_t* nodessol = (NodesSol_t*) malloc(sizeof(NodesSol_t)) ;
  unsigned int    i ;
  
  if(!nodessol) arret("NodesSol_Create") ;
  

  /* Verification */
  {
    unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
    Node_t* node = Mesh_GetNode(mesh) ;
    Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    unsigned int n_u = 0 ;
    
    for(i = 0 ; i < NbOfNodes ; i++) n_u += Node_GetNbOfUnknowns(node + i) ;
    
    if(n_u != Nodes_GetNbOfDOF(nodes)) {
      arret("NodesSol_Create") ;
    }
  }
  

  {
    /* Allocation of memory for the nodesol */
    {
      unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
      NodeSol_t* nodesol = (NodeSol_t*) malloc(NbOfNodes*sizeof(NodeSol_t)) ;
    
      if(!nodesol) {
        arret("NodesSol_Create(1): impossible d\'allouer la memoire") ;
      }
    
      NodesSol_GetNodeSol(nodessol) = nodesol ;
      NodesSol_GetNbOfNodes(nodessol) = NbOfNodes ;
      
      {
        Nodes_t* nodes = Mesh_GetNodes(mesh) ;
        Node_t* node = Mesh_GetNode(mesh) ;
        unsigned int    n_u = Nodes_GetNbOfDOF(nodes) ;
        
        /* The total nb of dof */
        NodesSol_GetNbOfDOF(nodessol) = n_u ;
        
        /* The nb of unknowns per nodesol */
        for(i = 0 ; i < NbOfNodes ; i++) {
          NodeSol_t* nodesol_i = NodesSol_GetNodeSol(nodessol) + i ;
        
          NodeSol_GetNbOfUnknowns(nodesol_i) = Node_GetNbOfUnknowns(node + i) ;
        }
      }
    }
  }
  
  NodesSol_AllocateMemory(nodessol) ;
  
  return(nodessol) ;
}



void (NodesSol_AllocateMemory)(NodesSol_t* nodessol)
{
  {
    /* Allocation of memory for the total nb of dof */
    unsigned int n_u = NodesSol_GetNbOfDOF(nodessol) ;
    double* u = (double*) calloc(n_u,sizeof(double)) ;
    
    if(!u) arret("NodesSol_AllocateMemory(1): allocation impossible") ;
    
    
    {
      unsigned int  NbOfNodes = NodesSol_GetNbOfNodes(nodessol) ;
      unsigned int  i ;
      
      for(i = 0 ; i < NbOfNodes ; i++) {
        NodeSol_t* nodesol_i = NodesSol_GetNodeSol(nodessol) + i ;
    
        if(i == 0) {
          NodeSol_GetUnknown(nodesol_i) = u ;
        } else {
          NodeSol_GetUnknown(nodesol_i) = NodeSol_GetUnknown(nodesol_i - 1) + NodeSol_GetNbOfUnknowns(nodesol_i - 1) ;
        }
      }
    }
  }
}



void NodesSol_Copy(NodesSol_t* nodessol_dest,NodesSol_t* nodessol_src)
/** Copy the nodal unknowns from nodessol_src to nodessol_dest */
{
  unsigned int nno = NodesSol_GetNbOfNodes(nodessol_src) ;
  
  NodesSol_GetNbOfNodes(nodessol_dest) = nno ;
  
  /* Nodal values */
  {
    unsigned int ndof = NodesSol_GetNbOfDOF(nodessol_src) ;
    double* u_s = NodesSol_GetDOF(nodessol_src) ;
    double* u_d = NodesSol_GetDOF(nodessol_dest) ;
    unsigned int i ;
    
    NodesSol_GetNbOfDOF(nodessol_dest) = ndof ;
    
    for(i = 0 ; i < ndof ; i++) {
      u_d[i] = u_s[i] ;
    }
  }
}



NodeSol_t* (NodeSol_GetDeepNodeSol)(NodeSol_t* nodesol,unsigned int depth)
{
  while(depth--) nodesol = NodeSol_GetPreviousNodeSol(nodesol) ;
  return(nodesol) ;
}



/* Not used from here */
#if 0
NodesSol_t* (NodesSol_Create1)(Mesh_t* mesh)
{
  NodesSol_t* nodessol = (NodesSol_t*) malloc(sizeof(NodesSol_t)) ;
  unsigned int    i ;
  
  if(!nodessol) arret("NodesSol_Create") ;
  

  /* Verification */
  {
    unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
    Node_t* node = Mesh_GetNode(mesh) ;
    Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    unsigned int n_u = 0 ;
    
    for(i = 0 ; i < NbOfNodes ; i++) n_u += Node_GetNbOfUnknowns(node + i) ;
    
    if(n_u != Nodes_GetNbOfDOF(nodes)) arret("NodesSol_Create") ;
  }
  
  
  /* Allocation of memory */
  {
    /* Allocation of memory for the nodesol */
    {
      unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
      NodeSol_t* nodesol = (NodeSol_t*) malloc(NbOfNodes*sizeof(NodeSol_t)) ;
    
      if(!nodesol) {
        arret("NodesSol_Create (1) : impossible d\'allouer la memoire") ;
      }
    
      NodesSol_GetNodeSol(nodessol) = nodesol ;
      NodesSol_GetNbOfNodes(nodessol) = NbOfNodes ;
      
      {
        Nodes_t* nodes = Mesh_GetNodes(mesh) ;
        Node_t* node = Mesh_GetNode(mesh) ;
        unsigned int    n_u = Nodes_GetNbOfDOF(nodes) ;
        
        /* The total nb of dof */
        NodesSol_GetNbOfDOF(nodessol) = n_u ;
        
        /* The nb of unknowns per nodesol */
        for(i = 0 ; i < NbOfNodes ; i++) {
          NodeSol_t* nodesol_i = NodesSol_GetNodeSol(nodessol) + i ;
        
          NodeSol_GetNbOfUnknowns(nodesol_i) = Node_GetNbOfUnknowns(node + i) ;
        }
      }
    }
  
    /* Allocation of memory for the total nb of dof */
    {
      unsigned int  NbOfNodes = NodesSol_GetNbOfNodes(nodessol) ;
      unsigned int n_u = NodesSol_GetNbOfDOF(nodessol) ;
      double* u = (double*) calloc(n_u,sizeof(double)) ;
    
      if(!u) arret("NodesSol_Create (2) : impossible d\'allouer la memoire") ;
      
  
      for(i = 0 ; i < NbOfNodes ; i++) {
        NodeSol_t* nodesol_i = NodesSol_GetNodeSol(nodessol) + i ;
    
        if(i == 0) {
          NodeSol_GetUnknown(nodesol_i) = u ;
        } else {
          NodeSol_GetUnknown(nodesol_i) = NodeSol_GetUnknown(nodesol_i - 1) + NodeSol_GetNbOfUnknowns(nodesol_i - 1) ;
        }
      }
    }
  }
  
  return(nodessol) ;
}
#endif
