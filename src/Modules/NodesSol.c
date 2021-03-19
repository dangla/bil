#include <stdio.h>
#include "Mry.h"
#include "NodesSol.h"
#include "Message.h"




/* Extern functions */

NodesSol_t* (NodesSol_Create)(Mesh_t* mesh)
{
  NodesSol_t* nodessol = (NodesSol_t*) Mry_New(NodesSol_t) ;
  

  /* Verification */
  {
    unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
    Node_t* node = Mesh_GetNode(mesh) ;
    Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    unsigned int n_u = 0 ;
    unsigned int    i ;
    
    for(i = 0 ; i < NbOfNodes ; i++) n_u += Node_GetNbOfUnknowns(node + i) ;
    
    if(n_u != Nodes_GetNbOfDOF(nodes)) {
      arret("NodesSol_Create") ;
    }
  }
  

  {
    /* Allocation of memory for the nodesol */
    {
      unsigned int  NbOfNodes = Mesh_GetNbOfNodes(mesh) ;
      Nodes_t* nodes = Mesh_GetNodes(mesh) ;
      unsigned int    n_u = Nodes_GetNbOfDOF(nodes) ;
      NodeSol_t* nodesol = (NodeSol_t*) Mry_New(NodeSol_t[NbOfNodes]) ;
      double* u = (double*) Mry_New(double[n_u]) ;
    
      NodesSol_GetNodeSol(nodessol) = nodesol ;
      NodesSol_GetNbOfNodes(nodessol) = NbOfNodes ;
      NodesSol_GetNbOfDOF(nodessol) = n_u ;
      NodesSol_GetNodalValue(nodessol) = u ;
      
      
      {
        Node_t* node = Mesh_GetNode(mesh) ;
        unsigned int    i ;
        
        /* The nb of unknowns per nodesol */
        for(i = 0 ; i < NbOfNodes ; i++) {
          NodeSol_t* nodesol_i = nodesol + i ;
          Node_t* node_i = node + i ;
        
          NodeSol_GetNbOfUnknowns(nodesol_i) = Node_GetNbOfUnknowns(node_i) ;
    
          NodeSol_GetUnknown(nodesol_i) = u ;
          
          u += NodeSol_GetNbOfUnknowns(nodesol_i) ;
        }
      }
    }
  }
  
  //NodesSol_AllocateMemory(nodessol) ;
  
  return(nodessol) ;
}



void NodesSol_Copy(NodesSol_t* nodessol_dest,NodesSol_t* nodessol_src)
/** Copy the nodal unknowns from nodessol_src to nodessol_dest */
{
  unsigned int nno = NodesSol_GetNbOfNodes(nodessol_src) ;
  
  NodesSol_GetNbOfNodes(nodessol_dest) = nno ;
  
  /* Nodal values */
  {
    unsigned int ndof = NodesSol_GetNbOfDOF(nodessol_src) ;
    double* u_s = NodesSol_GetNodalValue(nodessol_src) ;
    double* u_d = NodesSol_GetNodalValue(nodessol_dest) ;
    unsigned int i ;
    
    NodesSol_GetNbOfDOF(nodessol_dest) = ndof ;
    
    for(i = 0 ; i < ndof ; i++) {
      u_d[i] = u_s[i] ;
    }
  }
}



/* Intern functions */



/* Not used from here */
#if 0
static void           (NodesSol_AllocateMemory)(NodesSol_t*) ;

void (NodesSol_AllocateMemory)(NodesSol_t* nodessol)
{
  {
    /* Allocation of memory for the total nb of dof */
    //unsigned int n_u = NodesSol_GetNbOfDOF(nodessol) ;
    //double* u = (double*) Mry_New(double[n_u]) ;
    
    double* u = NodesSol_GetNodalValue(nodessol) ;
    
    
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
#endif

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
