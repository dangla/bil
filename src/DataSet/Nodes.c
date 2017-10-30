#include <stdlib.h>
#include "Message.h"
#include "Nodes.h"



void  Nodes_CreateMore(Nodes_t* nodes)
/** Matrix row and matrix column indexes */
{
  int    n_no = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;
  int    i ;
  int    n_inc ;
  int    n_equ ;



  /* Number of unknowns */
  n_inc = 0 ;
  for(i = 0 ; i < n_no ; i++) n_inc += Node_GetNbOfUnknowns(node + i) ;
  
  /* Number of equations */
  n_equ = 0 ;
  for(i = 0 ; i < n_no ; i++) n_equ += Node_GetNbOfEquations(node + i) ;
  
  if(n_inc != n_equ) {
    arret("Nodes_CreateMore(1): the numbers of unknowns and equations are different") ;
  }

  /* Nb of degrees of freedom */
  Nodes_GetNbOfDOF(nodes) = n_inc ;



  /* Allocation of space for the matrix column indexes */
  {
    size_t sz = n_inc*sizeof(int) ;
    int* colind = (int*) malloc(sz) ;
  
    if(!colind) arret("Nodes_CreateMore (2) : impossible d\'allouer la memoire") ;
  
    Node_GetMatrixColumnIndex(node) = colind ;
  
    for(i = 1 ; i < n_no ; i++) {
      Node_GetMatrixColumnIndex(node + i) = Node_GetMatrixColumnIndex(node + i - 1) + Node_GetNbOfUnknowns(node + i - 1) ;
    }
  }

  /* Allocation of space for the matrix row indexes */
  {
    size_t sz = n_equ*sizeof(int) ;
    int* rowind = (int*) malloc(sz) ;
  
    if(!rowind) arret("Nodes_CreateMore (3) : impossible d\'allouer la memoire") ;
  
    Node_GetMatrixRowIndex(node) = rowind ;
  
    for(i = 1 ; i < n_no ; i++) {
      Node_GetMatrixRowIndex(node + i) = Node_GetMatrixRowIndex(node + i - 1) + Node_GetNbOfEquations(node + i - 1) ;
    }
  }
}



/*
int   Nodes_ComputeNbOfMatrixRows(Nodes_t* nodes)
{
  int    i ;
  int    NbOfMatrixRows = 0 ;
  int    NbOfMatrixColumns = 0 ;

  NbOfMatrixRows = 0 ;
  NbOfMatrixColumns = 0 ;
  for(i = 0 ; i < (int) Nodes_GetNbOfNodes(nodes) ; i++) {
    Node_t* node_i = Nodes_GetNode(nodes) + i ;
    int   j ;

    for(j = 0 ; j < Node_GetNbOfUnknowns(node_i) ; j++) {
      if(Node_GetMatrixColumnIndex(node_i)[j] >= 0) {
        NbOfMatrixColumns += 1 ;
      }
    }

    for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
      if(Node_GetMatrixRowIndex(node_i)[j] >= 0) {
        NbOfMatrixRows += 1 ;
      }
    }
  }
  
  if(NbOfMatrixColumns != NbOfMatrixRows) arret("Mesh_ComputeNbOfMatrixRows") ;
  
  return(NbOfMatrixRows) ;
}
*/
