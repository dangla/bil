#include <stdlib.h>
#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "Model.h"
#include "Element.h"
#include "Nodes.h"




Nodes_t*  Nodes_New(const int nn,const int dim,const int nc)
{
  Nodes_t* nodes = (Nodes_t*) Mry_New(Nodes_t) ;
  
  Nodes_GetNbOfNodes(nodes) = nn ;
  Nodes_GetNbOfConnectivities(nodes) = nc ;
  
  /* Allocation of space for the nodes */
  {
    Node_t* node = (Node_t*) Mry_New(Node_t[nn]) ;
    
    Nodes_GetNode(nodes) = node ;
  }
  
  /* Allocation of space for the coordinates */
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    double* x = (double*) Mry_New(double[nn*dim]) ;
    
    Node_GetCoordinate(node) = x ;
  }
  
  /* Allocation of space for the pointers to "element" */
  {
    Element_t** pel = (Element_t**) Mry_New(Element_t*[nc]) ;
    Node_t* node = Nodes_GetNode(nodes) ;
        
    Node_GetPointerToElement(node) = pel ;
  }
  
  /* Initialization */
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    Element_t** pel = Node_GetPointerToElement(node) ;
    double* x = Node_GetCoordinate(node) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      Node_t* node_i = node + i ;
      
      Node_GetNodeIndex(node_i)         = i ;
      Node_GetCoordinate(node_i)        = x + i*dim ;
      Node_GetNbOfEquations(node_i)     = 0 ;
      Node_GetNbOfUnknowns(node_i)      = 0 ;
      Node_GetNameOfEquation(node_i)    = NULL ;
      Node_GetNameOfUnknown(node_i)     = NULL ;
      Node_GetObValIndex(node_i)        = 0 ;
      Node_GetMatrixColumnIndex(node_i) = NULL ;
      Node_GetMatrixRowIndex(node_i)    = NULL ;
      Node_GetNodeSol(node_i)           = NULL ;
      Node_GetNbOfElements(node_i)      = 0 ;  
      Node_GetPointerToElement(node_i)  = pel ;
    }
  }
  
  return(nodes) ;
}



void Nodes_Delete(void* self)
{
  Nodes_t** pnodes = (Nodes_t**) self ;
  Nodes_t*   nodes = *pnodes ;
  
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    double* x = Node_GetCoordinate(node) ;
    Element_t** pel = Node_GetPointerToElement(node) ;
    
    free(node) ;
    free(x) ;
    free(pel) ;
  }
  
  free(nodes) ;
}




void  Nodes_CreateMore(Nodes_t* nodes)
/** Allocate memory space at each node for:
 *  - the names of unknowns/equations
 *  - the indexes of matrix rows and matrix columns
 *  - the indexes of objective values
 */
{
  int    n_no = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;

  /*  The number of unknowns/equations per node:
   *  - Node_GetNbOfEquations(node_i)
   *  - Node_GetNbOfUnknowns(node_i)
   */
  {
    int in ;
  
    for(in = 0 ; in < n_no ; in++) {
      Node_t* node_i = node + in ;
      int nel = Node_GetNbOfElements(node_i) ;
      char* node_name_unk[Node_MaxNbOfEquations] ;
      char* node_name_eqn[Node_MaxNbOfEquations] ;
      int je ;
          
      Node_GetNbOfUnknowns(node_i)  = 0 ;
      Node_GetNbOfEquations(node_i) = 0 ;
          
      for(je = 0 ; je < nel ; je++) {
        Element_t* el = Node_GetElement(node_i,je) ;
        Material_t* mat = Element_GetMaterial(el) ;
    
        if(mat) {
          int neq = Element_GetNbOfEquations(el) ;
          char** mat_name_unk = Material_GetNameOfUnknown(mat) ;
          char** mat_name_eqn = Material_GetNameOfEquation(mat) ;
          int ieq ;
        
          for(ieq = 0 ; ieq < neq ; ieq++) {

            /* Continuity of unknowns: 
              * - number of unknowns: Node_GetNbOfUnknowns(node_i)
              */
            {
              int jun ;

              for(jun = 0 ; jun < Node_GetNbOfUnknowns(node_i) ; jun++) {
                if(!strcmp(node_name_unk[jun],mat_name_unk[ieq])) break ;
              }
          
              {
                /* Skip if this unknown has already been met */
                if(jun == Node_GetNbOfUnknowns(node_i)) {
                  /* Increment the number of unknowns at node i */
                  Node_GetNbOfUnknowns(node_i) += 1 ;
                  /* Set the name of the unknown ij */
                  if(Node_GetNbOfUnknowns(node_i) < Node_MaxNbOfEquations) {
                    node_name_unk[jun] = mat_name_unk[ieq] ;
                  } else {
                    arret("Too many equations!") ;
                  }
                }
              }
            }

            /* Continuity of equations:
              * - number of equations: Node_GetNbOfEquations(node_i)
              */
            {
              int jeq ;
            
              for(jeq = 0 ; jeq < Node_GetNbOfEquations(node_i) ; jeq++) {
                if(!strcmp(node_name_eqn[jeq],mat_name_eqn[ieq])) break ;
              }
          
              {
                /* Skip if this equation has already been met */
                if(jeq == Node_GetNbOfEquations(node_i)) {
                  /* Increment the number of equations */
                  Node_GetNbOfEquations(node_i) += 1 ;
                  /* Set the name of the equation ij */
                  if(Node_GetNbOfEquations(node_i) < Node_MaxNbOfEquations) {
                    node_name_eqn[jeq] = mat_name_eqn[ieq] ;
                  } else {
                    arret("Too many equations!") ;
                  }
                }
              }
            }
          }
        }
      }
      
      if(Node_GetNbOfUnknowns(node_i) != Node_GetNbOfEquations(node_i)) {
        int ind = Node_GetNodeSol(node_i) ;
        
        arret("Nodes_CreateMore: nb of unknowns and equations not equal at node %d",ind) ;
      }
    }
  }



  /* The total nb of degrees of freedom */
  {
    int    n_inc = 0 ;
    int    n_equ = 0 ;
    int    i ;
  
    /* Number of unknowns */
    for(i = 0 ; i < n_no ; i++) {
      n_inc += Node_GetNbOfUnknowns(node + i) ;
    }
  
    /* Number of equations */
    for(i = 0 ; i < n_no ; i++) {
      n_equ += Node_GetNbOfEquations(node + i) ;
    }
  
    if(n_inc != n_equ) {
      arret("Nodes_CreateMore(1): the numbers of unknowns and equations are different") ;
    }

    /* Nb of degrees of freedom */
    Nodes_GetNbOfDOF(nodes) = n_inc ;
  }
  
  
  /* Allocate memory space for names of equations and unknowns */
  {
    {
      int n_dof = Nodes_GetNbOfDOF(nodes) ;
      char** uname = (char**) Mry_New(char*[2*n_dof]) ;
      char** ename = uname + n_dof ;
      int in ;
  
      for(in = 0 ; in < n_no ; in++) {
        Node_GetNameOfUnknown(node + in)  = uname ;
        Node_GetNameOfEquation(node + in) = ename ;
        uname += Node_GetNbOfUnknowns(node + in) ;
        ename += Node_GetNbOfEquations(node + in) ;
      }
    }
  }


  /* Allocation of space for the matrix column indexes */
  {
    int n_dof = Nodes_GetNbOfDOF(nodes) ;
    int* colind = (int*) Mry_New(int[n_dof]) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_GetMatrixColumnIndex(node + i) = colind ;
      colind += Node_GetNbOfUnknowns(node + i) ;
    }
  }


  /* Allocation of space for the matrix row indexes */
  {
    int n_dof = Nodes_GetNbOfDOF(nodes) ;
    int* rowind = (int*) Mry_New(int[n_dof]) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_GetMatrixRowIndex(node + i) = rowind ;
      rowind += Node_GetNbOfEquations(node + i) ;
    }
  }
  
  
  
  /* Allocation of space for the index of objective values */
  {
    unsigned int n_dof = Nodes_GetNbOfDOF(nodes) ;
    unsigned short int* index = (unsigned short int*) Mry_New(unsigned short int[n_dof]) ;
    
    {
      int i ;
      
      for(i = 0 ; i < n_no ; i++) {
        Node_GetObValIndex(node + i) = index ;
        index += Node_GetNbOfEquations(node + i) ;
      }
    }
  }

  /* Space allocation for buffer */
  {
    Buffer_t* buf = Buffer_Create(Node_SizeOfBuffer) ;
    int i ;
  
    /* ATTENTION: same memory space (buffer) for all nodes */
    for(i = 0 ; i < n_no ; i++) {
      Node_GetBuffer(node + i) = buf ;
    }
  }
}



void Nodes_DeleteMore(void* self)
{
  Nodes_t** pnodes = (Nodes_t**) self ;
  Nodes_t*   nodes = *pnodes ;
  
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    char** uname = Node_GetNameOfUnknown(node) ;
    int* colind = Node_GetMatrixColumnIndex(node) ;
    int* rowind = Node_GetMatrixRowIndex(node) ;
    unsigned short int* index = Node_GetObValIndex(node) ;
    Buffer_t* buf = Node_GetBuffer(node) ;
    
    free(uname) ;
    free(colind) ;
    free(rowind) ;
    free(index) ;
    Buffer_Delete(&buf) ;
  }
}



void  Nodes_InitializeMatrixRowColumnIndexes(Nodes_t* nodes)
/** Initialization to arbitrarily negative value (-1) 
  * so as to eliminate dof of isolated nodes or
  * dof of negative position (see below).
  */
{
  {
    Node_t* no = Nodes_GetNode(nodes) ;
    int n_no = Nodes_GetNbOfNodes(nodes) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_t* node_i = no + i ;
      int   j ;
    
      for(j = 0 ; j < Node_GetNbOfUnknowns(node_i) ; j++)  {
        char* unk = Node_GetNameOfUnknown(node_i)[j] ;
        
        Node_EliminateMatrixColumnIndexForBCond(node_i,unk) ;
      }
    
      for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
        char* eqn = Node_GetNameOfEquation(node_i)[j] ;
              
        Node_EliminateMatrixRowIndexForBCond(node_i,eqn) ;
      }
    }
  }
}



void Nodes_SetMatrixRowColumnIndexes(Nodes_t* nodes,DataFile_t* datafile)
/** Set up the system */
{
  int   n_no = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;
  int NbOfMatrixRows = 0 ;
  int NbOfMatrixColumns = 0 ;
  
  {
    int*  perm = NULL ;
    int    i ;
    
    if(datafile) {
      perm = DataFile_ReadInversePermutationOfNodes(datafile,n_no) ;
    }
    
    
    /* In the order of the permuted nodes */
    for(i = 0 ; i < n_no ; i++) {
      int i_perm = (perm) ? perm[i] : i ;
      Node_t* node_i = node + i_perm ;
      int   j ;

      for(j = 0 ; j < Node_GetNbOfUnknowns(node_i) ; j++) {
        int icol = Node_GetMatrixColumnIndex(node_i)[j] ;
        
        if(icol >= 0) {
          Node_GetMatrixColumnIndex(node_i)[j] = NbOfMatrixColumns ;
          NbOfMatrixColumns += 1 ;
        }
      }

      for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
        int irow = Node_GetMatrixRowIndex(node_i)[j] ;
        
        if(irow >= 0) {
          Node_GetMatrixRowIndex(node_i)[j] = NbOfMatrixRows ;
          NbOfMatrixRows += 1 ;
        }
      }
    }
  
    if(perm) free(perm) ;
  }
  
  
  if(NbOfMatrixColumns != NbOfMatrixRows) {
    arret("Nodes_SetMatrixRowColumnIndexes") ;
  }
  
  Nodes_GetNbOfMatrixRows(nodes) = NbOfMatrixRows ;
  Nodes_GetNbOfMatrixColumns(nodes) = NbOfMatrixColumns ;
}



int Nodes_ComputeNbOfUnknownFields(Nodes_t* nodes)
/** Compute the nb of unknown (scalar) fields such as 
 *  displacement, temperature, concentration, ....
 **/
{
  int n_unknownnames = 0 ;
  int n_nodes = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;
    
    
  {
    char unknownname[Model_MaxNbOfEquations][Model_MaxLengthOfKeyWord] ;
    int i ;
    
    for(i = 0 ; i < n_nodes ; i++) {
      Node_t* node_i = node + i ;
      int neq = Node_GetNbOfEquations(node_i) ;
      int ieq ;
      
      for(ieq = 0 ; ieq < neq ; ieq++) {
        char* name_node  = Node_GetNameOfUnknown(node_i)[ieq] ;
        int    j ;
        
        for(j = 0 ; j < n_unknownnames ; j++) {
          char* name = unknownname[j] ;
          
          if(!strcmp(name,name_node)) break ;
        }
        
        if(j == n_unknownnames) {
          char* name = unknownname[j] ;
          
          strcpy(name,name_node) ;
          
          n_unknownnames += 1 ;
          
          if(n_unknownnames > Model_MaxNbOfEquations) {
            arret("Nodes_ComputeNbOfUnknownFields: too many equations") ;
          }
        }
      }
    }
  }
    
  return(n_unknownnames) ;
}


  
  

void Nodes_InitializeObValIndexes(Nodes_t* nodes)
/** Initialize the objective value indexes at the nodes 
 ** (Node_GetObValIndex).
 **/
{
  int n_nodes = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;
  ObVals_t* obvals = Nodes_GetObjectiveValues(nodes) ;
    
  {
    int i ;
    
    for(i = 0 ; i < n_nodes ; i++) {
      Node_t* node_i = node + i ;
      int neq = Node_GetNbOfEquations(node_i) ;
      int    ieq ;
      
      for(ieq = 0 ; ieq < neq ; ieq++) {
        char* name_node  = Node_GetNameOfUnknown(node_i)[ieq] ;
        int   j = ObVals_FindObValIndex(obvals,name_node) ;
        
        if(j >= 0) {
          
          Node_GetObValIndex(node_i)[ieq] = j ;
          
        } else {
          
          arret("Nodes_InitializeObValIndex: not enough objective values") ;
          
        }
      }
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
