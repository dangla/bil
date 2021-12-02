#include <stdlib.h>
#include <ctype.h>
#include "String.h"
#include "Math.h"
#include "Message.h"
#include "Mry.h"
#include "Model.h"
#include "Node.h"


static int      (Node_CommonUnknownPositionIndexAtOverlappingNodes)(const Node_t*,const int,const char*) ;
static int      (Node_CommonEquationPositionIndexAtOverlappingNodes)(const Node_t*,const int,const char*) ;




#if 0
Node_t*  Node_Create(const int dim)
{
  Node_t* node = (Node_t*) Mry_New(Node_t) ;
  
  {
    double* x = (double*) Mry_New(double[dim]) ;
    
    Node_GetCoordinate(node) = x ;
  }
  
  return(node) ;
}
#endif



void (Node_Delete)(void* self)
{
  Node_t* node = (Node_t*) self ;

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
    Buffer_t* buf = Node_GetBuffer(node) ;
      
    if(buf) {
      Buffer_Delete(buf) ;
      free(buf) ;
      Node_GetBuffer(node) = NULL ;
    }
  }
}





int    (Node_FindUnknownPositionIndex)(const Node_t* node,const char* s)
/** Find the unknown position index whose name is pointed to by s */
{
  int n = Node_GetNbOfUnknowns(node) ;
  const char* const* ss = (const char* const*) Node_GetNameOfUnknown(node) ;
  int    i = String_FindPositionIndex(s,ss,n) ;
  
  return(i) ;

  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) {
      if(!strncmp(s,ss[i],strlen(s))) break ;
    }
    if(i == n) i = -1 ;
  }

  return(i) ;
}



int    (Node_FindEquationPositionIndex)(const Node_t* node,const char* s)
/** Find the equation position index whose name is pointed to by s */
{
  int n = Node_GetNbOfEquations(node) ;
  const char* const* ss = (const char* const*) Node_GetNameOfEquation(node) ;
  int    i = String_FindPositionIndex(s,ss,n) ;
  
  return(i) ;

  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) if(!strncmp(s,ss[i],strlen(s))) break ;
    if(i == n) i = -1 ;
  }

  return(i) ;
}



void Node_MakeUnknownContinuousAtOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
  int index = Node_CommonUnknownPositionIndexAtOverlappingNodes(node_k,n_node_k,name) ;
    
  {
    int j ;
    
    for(j = 0 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindUnknownPositionIndex(node_j,name) ;
      
      if(jj >= 0) {
        Node_GetMatrixColumnIndex(node_j)[jj] = index ;
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}



void Node_MakeEquationContinuousAtOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
  int index = Node_CommonEquationPositionIndexAtOverlappingNodes(node_k,n_node_k,name) ;
    
  {
    int j ;
    
    for(j = 0 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindEquationPositionIndex(node_j,name) ;
      
      if(jj >= 0) {
        Node_GetMatrixRowIndex(node_j)[jj] = index ;
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}




void Node_EliminateMatrixColumnIndexForOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
    
  {
    Node_t* node_i = node_k ;
    int ii = Node_FindUnknownPositionIndex(node_i,name) ;
    int j ;
      
    for(j = 1 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindUnknownPositionIndex(node_j,name) ;
        
      if(ii >= 0) {
        if(jj >= 0) {
          int ki = Node_GetMatrixColumnIndex(node_i)[ii] ;
          int kj = Node_GetMatrixColumnIndex(node_j)[jj] ;
            
          if(ki >= 0) {
              
            if(kj == ki) {
              Node_GetMatrixColumnIndex(node_j)[jj] = Node_IndexForOverlappingNode ;
            }
            
          } else if(kj >= 0) {
              
            node_i = node_j ;
            ii = jj ;
              
          }
        }
          
      } else {
          
        node_i = node_j ;
        ii = jj ;
          
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}


        
        
void Node_EliminateMatrixRowIndexForOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
    
  {
    Node_t* node_i = node_k ;
    int ii = Node_FindEquationPositionIndex(node_i,name) ;
    int j ;
      
    for(j = 1 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindEquationPositionIndex(node_j,name) ;
        
      if(ii >= 0) {
        if(jj >= 0) {
          int ki = Node_GetMatrixRowIndex(node_i)[ii] ;
          int kj = Node_GetMatrixRowIndex(node_j)[jj] ;
            
          if(ki >= 0) {
              
            if(kj == ki) {
              Node_GetMatrixRowIndex(node_j)[jj] = Node_IndexForOverlappingNode ;
            }
            
          } else if(kj >= 0) {
              
            node_i = node_j ;
            ii = jj ;
              
          }
        }
          
      } else {
          
        node_i = node_j ;
        ii = jj ;
          
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}




void Node_UpdateMatrixColumnIndexForOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
    
  {
    Node_t* node_i = node_k ;
    int ii = Node_FindUnknownPositionIndex(node_i,name) ;
    int j ;
      
    for(j = 1 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindUnknownPositionIndex(node_j,name) ;
        
      if(ii >= 0) {
        if(jj >= 0) {
          int ki = Node_GetMatrixColumnIndex(node_i)[ii] ;
          int kj = Node_GetMatrixColumnIndex(node_j)[jj] ;
            
          if(ki >= 0) {
              
            if(kj == Node_IndexForOverlappingNode) {
              Node_GetMatrixColumnIndex(node_j)[jj] = ki ;
            }
            
          } else if(kj >= 0) {
              
            node_i = node_j ;
            ii = jj ;
              
          }
        }
          
      } else {
          
        node_i = node_j ;
        ii = jj ;
          
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}




void Node_UpdateMatrixRowIndexForOverlappingNodes(const Node_t* node,const char* name)
{
  int n_node_k ;
  Node_t* node_k = Node_OverlappingNodes(node,&n_node_k) ;
    
  {
    Node_t* node_i = node_k ;
    int ii = Node_FindEquationPositionIndex(node_i,name) ;
    int j ;
      
    for(j = 1 ; j < n_node_k ; j++) {
      Node_t* node_j = node_k + j ;
      int jj = Node_FindEquationPositionIndex(node_j,name) ;
        
      if(ii >= 0) {
        if(jj >= 0) {
          int ki = Node_GetMatrixRowIndex(node_i)[ii] ;
          int kj = Node_GetMatrixRowIndex(node_j)[jj] ;
            
          if(ki >= 0) {
              
            if(kj == Node_IndexForOverlappingNode) {
              Node_GetMatrixRowIndex(node_j)[jj] = ki ;
            }
            
          } else if(kj >= 0) {
              
            node_i = node_j ;
            ii = jj ;
              
          }
        }
          
      } else {
          
        node_i = node_j ;
        ii = jj ;
          
      }
    }
  }
  
  Node_FreeBufferFrom(node,node_k) ;
}






/* Intern functions */

Node_t* Node_OverlappingNodes3(const Node_t* node,int* n_overlappingnodes,Node_t* overlappingnode)
/** Compute the overlapping nodes attached to the given node 
 *  belonging to zero-thickness element, including that node itself.
 * 
 *  Input:
 *  node: a pointer to Node_t for a node belonging to zero-thickness element
 *  overlappingnode: should be a null pointer at the first call. 
 *  A non-null pointer is used in subsequent recursive calls.
 * 
 *  Output:
 *  n_overlappingnodes: a pointer to int storing the nb of overlapping nodes.
 * 
 *  Return a pointer to Node_t with the overlapping nodes.
 **/
{
  Node_t* new_node ;
  int nb_new_node ;
  

  #define n_max  (20)
  if(overlappingnode) {
    
    new_node = overlappingnode ;
    nb_new_node = n_overlappingnodes[0] ;
    
    if(nb_new_node >= n_max) {
      arret("Node_OverlappingNodes3: too many nodes, increase n_max!") ;
    }
    
  } else {
    size_t SizeNeeded = n_max*sizeof(Node_t) ;
    
    new_node = (Node_t*) Node_AllocateInBuffer(node,SizeNeeded) ;
    
    /* We include the node itself in the list of overlapping nodes */
    new_node[0] = node[0] ;
    nb_new_node = 1 ;
  }
  #undef n_max
  
  
  {
    int nel = Node_GetNbOfElements(node) ;
    int ie ;

    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* element = Node_GetElement(node,ie) ;
      int in = Element_FindNodeIndex(element,node) ;

      if(in < 0) continue ;

      if(Element_HasZeroThickness(element)) {
        int jn = Element_OverlappingNode(element,in) ;

        /* No overlapping node! */
        if(jn == in) continue ;

        {
          Node_t* node_j = Element_GetNode(element,jn) ;
          int k ;

          /* Check if this node is already in the list */
          for(k = 0 ; k < nb_new_node ; k++) {
            if(Node_GetNodeIndex(new_node + k) == Node_GetNodeIndex(node_j)) {
              break ;
            }
          }

          /* New node to be added to the list (recursive case) */
          if(k == nb_new_node) {
            new_node[k] = node_j[0] ;
            nb_new_node += 1 ;
            /* Recursivity: we must do the same for this node */
            Node_OverlappingNodes3(node_j,&nb_new_node,new_node) ;
          }
                
          /* Jump to the next element until all nodes have been found (base case) */
        }
      }
    }
  }
  
  
  /* Only for the base case: order the new_node array 
   * from the smallest index to the largest index */
  if(overlappingnode) {
    int i ;
    
    for(i = 0 ; i < nb_new_node ; i++) {
      int j ;
        
      for(j = i + 1 ; j < nb_new_node ; j++) {
        int kj = Node_GetNodeIndex(new_node + j) ;
        /* ki should be re-calculated at each j in case of swapping */
        int ki = Node_GetNodeIndex(new_node + i) ;
          
        if(kj < ki) {
          Math_Swap(new_node[i],new_node[j],Node_t) ;
        }
      }
    }
  }
  
  
  n_overlappingnodes[0] = nb_new_node ;
  
  return(new_node) ;
}



int Node_CommonUnknownPositionIndexAtOverlappingNodes(const Node_t* node,const int n_node,const char* name)
{
  int index_min = 1 ;
  int index_max = 0 ;
  
  {
    int j ;
    
    for(j = 0 ; j < n_node ; j++) {
      const Node_t* node_j = node + j ;
      int jj = Node_FindUnknownPositionIndex(node_j,name) ;
      
      if(jj >= 0) {
        int kj = Node_GetMatrixColumnIndex(node_j)[jj] ;
        
        if(index_min > index_max) {
          index_max = kj ;
          index_min = kj ;
        } else {
          if(kj > index_max) index_max = kj ;
          if(kj < index_min) index_min = kj ;
        }
      }
    }
  }
  
  return(index_min) ;
}



int Node_CommonEquationPositionIndexAtOverlappingNodes(const Node_t* node,const int n_node,const char* name)
{
  int index_min = 1 ;
  int index_max = 0 ;
  
  {
    int j ;
    
    for(j = 0 ; j < n_node ; j++) {
      const Node_t* node_j = node + j ;
      int jj = Node_FindEquationPositionIndex(node_j,name) ;
      
      if(jj >= 0) {
        int kj = Node_GetMatrixRowIndex(node_j)[jj] ;
        
        if(index_min > index_max) {
          index_max = kj ;
          index_min = kj ;
        } else {
          if(kj > index_max) index_max = kj ;
          if(kj < index_min) index_min = kj ;
        }
      }
    }
  }
  
  return(index_min) ;
}



#if 0
int Node_UnknownIsPrescribedAtOverlappingNodes(const Node_t* node,const int n_node,const char* name)
{
  int j ;
    
  for(j = 0 ; j < n_node ; j++) {
    Node_t* node_j = node + j ;
    int jj = Node_FindUnknownPositionIndex(node_j,name) ;
      
    if(jj >= 0) {
      int kj = Node_GetMatrixColumnIndex(node_j)[jj] ;
          
      /* If one of the unknown is prescribed by the BC then
       * it must be prescribed at all overlapping nodes. */
      if(kj == Node_IndexForBCond) {
        return(1) ;
      }
    }
  }
  
  return(0) ;
}
#endif
