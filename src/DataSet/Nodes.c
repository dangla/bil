#include <stdlib.h>
#include "Message.h"
#include "DataFile.h"
#include "String_.h"
#include "Mry.h"
#include "Model.h"
#include "Element.h"
#include "Nodes.h"
#include "Buffers.h"



static int  Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode(Nodes_t*) ;



Nodes_t*  Nodes_New(const int nn,const int dim,const int nc)
{
  Nodes_t* nodes = (Nodes_t*) Mry_New(Nodes_t) ;
  
  Nodes_GetNbOfNodes(nodes) = nn ;
  Nodes_GetNbOfConnectivities(nodes) = nc ;
  
  /* Initialize the number of matrices to 1 by default */
  Nodes_GetNbOfMatrices(nodes) = 1 ;
  
  /* Allocation of space for the nodes */
  {
    Node_t* node = (Node_t*) Mry_New(Node_t[nn]) ;
    
    Nodes_GetNode(nodes) = node ;
  }
  
  /* Allocation of space for the coordinates */
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    double* x = (double*) Mry_New(double[nn*dim]) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      Node_t* node_i = node + i ;
      
      Node_GetCoordinate(node_i)        = x + i*dim ;
    }
  }
  
  /* Allocation of space for the pointers to "element" */
  {
    Element_t** pel = (Element_t**) Mry_New(Element_t*[nc]) ;
    Node_t* node = Nodes_GetNode(nodes) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      Node_t* node_i = node + i ;
      
      Node_GetPointerToElement(node_i)  = pel ;
    }
  }
  
  /* Initialization */
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      Node_t* node_i = node + i ;
      
      Node_GetNodeIndex(node_i)         = i ;
      Node_GetNbOfEquations(node_i)     = 0 ;
      Node_GetNbOfUnknowns(node_i)      = 0 ;
      Node_GetNameOfEquation(node_i)    = NULL ;
      Node_GetNameOfUnknown(node_i)     = NULL ;
      Node_GetSequentialIndexOfUnknown(node_i) = NULL ;
      Node_GetObValIndex(node_i)        = 0 ;
      Node_GetMatrixColumnIndex(node_i) = NULL ;
      Node_GetMatrixRowIndex(node_i)    = NULL ;
      Node_GetNodeSol(node_i)           = NULL ;
      Node_GetNbOfElements(node_i)      = 0 ;
      Node_GetBuffers(node_i)           = NULL ;
    }
  }
  
  return(nodes) ;
}




void  Nodes_CreateMore(Nodes_t* nodes)
/** Allocate memory space at each node for:
 *  - the names of unknowns/equations
 *  - the sequential indexes of unknowns/equations
 *  - the indexes of matrix rows and matrix columns
 *  - the indexes of objective values
 */
{
  int    n_no = Nodes_GetNbOfNodes(nodes) ;
  Node_t* node = Nodes_GetNode(nodes) ;
  
  Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode(nodes) ;
  
  
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
  
  
  /* Allocate memory space for the sequential indexes of unknowns/equations */
  {
    {
      int n_dof = Nodes_GetNbOfDOF(nodes) ;
      int* index = (int*) Mry_New(int[n_dof]) ;
      int in ;
  
      for(in = 0 ; in < n_no ; in++) {
        Node_GetSequentialIndexOfUnknown(node + in)  = index ;
        index += Node_GetNbOfUnknowns(node + in) ;
      }
    }
  }
  
  
  /* Allocation of space for the nb of matrix rows and columns */
  {
    int  n = Nodes_MaxNbOfMatrices ;
    unsigned int* nb_rows = (unsigned int*) Mry_New(int[2*n]) ;
    unsigned int* nb_cols = nb_rows + n ;
    
    Nodes_GetNbOfMatrixRows(nodes) = nb_rows ;
    Nodes_GetNbOfMatrixColumns(nodes) = nb_cols ;
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
    Buffers_t* buf = Buffers_Create(Node_SizeOfBuffer) ;
    int i ;
  
    /* ATTENTION: same memory space (buffer) for all nodes */
    for(i = 0 ; i < n_no ; i++) {
      Node_GetBuffers(node + i) = buf ;
    }
  }
}



void (Nodes_Delete)(void* self)
{
  Nodes_t* nodes = (Nodes_t*) self ;
    
  {
    Node_t* node = Nodes_GetNode(nodes) ;
    
    {
      double* x = Node_GetCoordinate(node) ;
      
      if(x) {
        free(x) ;
      }
    }
  
    {
      Element_t** pel = Node_GetPointerToElement(node) ;
    
      if(pel) {
        free(pel) ;
      }
    }
    
    {
      char** uname = Node_GetNameOfUnknown(node) ;
      
      if(uname) {
        free(uname) ;
      }
    }
    
    {
      int* seq = Node_GetSequentialIndexOfUnknown(node) ;
      
      if(seq) {
        free(seq) ;
      }
    }
    
    {
      int* colind = Node_GetMatrixColumnIndex(node) ;
      
      if(colind) {
        free(colind) ;
      }
    }
    
    {
      int* rowind = Node_GetMatrixRowIndex(node) ;
      
      if(rowind) {
        free(rowind) ;
      }
    }
    
    {
      unsigned short int* index = Node_GetObValIndex(node) ;
      
      if(index) {
        free(index) ;
      }
    }
    
    {
      Buffers_t* buf = Node_GetBuffers(node) ;
      
      if(buf) {
        Buffers_Delete(buf) ;
        free(buf) ;
      }
    }
    
    free(node) ;
  }

  {
    unsigned int* nb_rows = Nodes_GetNbOfMatrixRows(nodes) ;
    
    if(nb_rows) {
      free(nb_rows) ;
    }
  }
}




int  Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode(Nodes_t* nodes)
/** Update 
 *  - the nb of unknowns and equations per node
 *  - the total number of degrees of freedom.
 *  without using informations from the nodes themselves, i.e., 
 *  using only informations from the elements connected to the nodes.
 *  This can be done before allocating memory space for the nodes.
 *  Return the total number of degrees of freedom.
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
                if(String_Is(node_name_unk[jun],mat_name_unk[ieq])) break ;
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
                    arret("Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode: too many equations!") ;
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
                if(String_Is(node_name_eqn[jeq],mat_name_eqn[ieq])) break ;
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
                    arret("Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode: too many equations!") ;
                  }
                }
              }
            }
          }
        }
      }
      
      if(Node_GetNbOfUnknowns(node_i) != Node_GetNbOfEquations(node_i)) {
        int ind = Node_GetNodeIndex(node_i) ;
        
        {
          int n_unk = Node_GetNbOfUnknowns(node_i) ;
          int i ;
          
          printf("\n") ;
          printf("name of unknowns:") ;
          for(i = 0 ; i < n_unk ; i++) {
            printf(" %s",node_name_unk[i]) ;
          }
        }
        {
          int n_equ = Node_GetNbOfEquations(node_i) ;
          int i ;
          
          printf("\n") ;
          printf("name of equations:") ;
          for(i = 0 ; i < n_equ ; i++) {
            printf(" %s",node_name_eqn[i]) ;
          }
        }
        
        arret("Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode: nb of unknowns and equations not equal at node %d",ind) ;
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
      arret("Nodes_UpdateTheNbOfUnknownsAndEquationsPerNode: the numbers of unknowns and equations are different") ;
    }

    /* Nb of degrees of freedom */
    Nodes_GetNbOfDOF(nodes) = n_inc ;
  }

  return(Nodes_GetNbOfDOF(nodes)) ;
}



void  Nodes_InitializeMatrixRowColumnIndexes(Nodes_t* nodes)
/** Initialization to (arbitrarily) negative value
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
  int n_matrices = Nodes_GetNbOfMatrices(nodes) ;
  unsigned int* nbofmatrixrows = Nodes_GetNbOfMatrixRows(nodes) ;
  unsigned int* nbofmatrixcols = Nodes_GetNbOfMatrixColumns(nodes) ;
  
  if(n_matrices > Nodes_MaxNbOfMatrices) {
    arret("Nodes_SetMatrixRowColumnIndexes: too many matrices") ;
  }
    
  
  {
    int i ;
    
    for(i = 0 ; i < n_matrices ; i++) {
      nbofmatrixrows[i] = 0 ;
      nbofmatrixcols[i] = 0 ;
    }
  }
  
  
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
          int seq_ind = Node_GetSequentialIndexOfUnknown(node_i)[j] ;
          
          if(seq_ind >= 0) {
            int imatrix = (seq_ind < n_matrices) ? seq_ind : n_matrices - 1 ;
            
            Node_GetMatrixColumnIndex(node_i)[j] = nbofmatrixcols[imatrix] ;
            nbofmatrixcols[imatrix] += 1 ;
            
          } else {
            Message_FatalError("Nodes_SetMatrixRowColumnIndexes:") ;
            // I'm questionning this possibility?
            //Node_GetMatrixColumnIndex(node_i)[j] = -1 ;
          }
        }
      }

      for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
        int irow = Node_GetMatrixRowIndex(node_i)[j] ;
        
        if(irow >= 0) {
          int seq_ind = Node_GetSequentialIndexOfUnknown(node_i)[j] ;
          
          if(seq_ind >= 0) {
            int imatrix = (seq_ind < n_matrices) ? seq_ind : n_matrices - 1 ;
            
            Node_GetMatrixRowIndex(node_i)[j] = nbofmatrixrows[imatrix] ;
            nbofmatrixrows[imatrix] += 1 ;
            
          } else {
            Message_FatalError("Nodes_SetMatrixRowColumnIndexes:") ;
            // I'm questionning this possibility?
            //Node_GetMatrixRowIndex(node_i)[j] = -1 ;
          }
        }
      }
    }
  
    if(perm) free(perm) ;
  }
  
  
  {
    int i ;
    
    for(i = 0 ; i < n_matrices ; i++) {
      if(nbofmatrixcols[i] != nbofmatrixrows[i]) {
        arret("Nodes_SetMatrixRowColumnIndexes") ;
      }
    }
  }
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




#if 0
void  Nodes_CreateMore(Nodes_t* nodes)
/** Allocate memory space at each node for:
 *  - the names of unknowns/equations
 *  - the sequential indexes of unknowns/equations
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
        int ind = Node_GetNodeIndex(node_i) ;
        
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
  
  
  /* Allocate memory space for the sequential indexes of unknowns/equations */
  {
    {
      int n_dof = Nodes_GetNbOfDOF(nodes) ;
      int* index = (int*) Mry_New(int[n_dof]) ;
      int in ;
  
      for(in = 0 ; in < n_no ; in++) {
        Node_GetSequentialIndexOfUnknown(node + in)  = index ;
        index += Node_GetNbOfUnknowns(node + in) ;
      }
    }
  }
  
  
  /* Allocation of space for the nb of matrix rows and columns */
  {
    int  n = Nodes_MaxNbOfMatrices ;
    unsigned int* nb_rows = (unsigned int*) Mry_New(int[2*n]) ;
    unsigned int* nb_cols = nb_rows + n ;
    
    Nodes_GetNbOfMatrixRows(nodes) = nb_rows ;
    Nodes_GetNbOfMatrixColumns(nodes) = nb_cols ;
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
    Buffers_t* buf = Buffers_Create(Node_SizeOfBuffer) ;
    int i ;
  
    /* ATTENTION: same memory space (buffer) for all nodes */
    for(i = 0 ; i < n_no ; i++) {
      Node_GetBuffers(node + i) = buf ;
    }
  }
}
#endif
