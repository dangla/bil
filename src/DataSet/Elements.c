#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "Buffers.h"
#include "Geometry.h"
#include "Nodes.h"
#include "IntFcts.h"
#include "Materials.h"
#include "Node.h"
#include "Elements.h"
#include "Message.h"
#include "Math_.h"
#include "Curves.h"
#include "Mry.h"


static double   (Elements_ComputeMaximumSizeOfElements)(Elements_t*) ;
static double   (Elements_ComputeMinimumSizeOfElements)(Elements_t*) ;



Elements_t*  (Elements_New)(const int n,const int nc)
{
  Elements_t* elts = (Elements_t*) Mry_New(Elements_t) ;
  
  Elements_GetNbOfElements(elts) = n ;
  Elements_GetNbOfConnectivities(elts) = nc ;
  
  
  /* Allocation of space for the elements */
  {
    Element_t* el = (Element_t*) Mry_New(Element_t[n]) ;
    
    Elements_GetElement(elts) = el ;
  }
  
  
  /* Allocation of space for the pointers to "node" */
  {
    Node_t** pnode = (Node_t**) Mry_New(Node_t*[nc]) ;
    Element_t* el = Elements_GetElement(elts) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Element_t* el_i = el + i ;
      
      Element_GetPointerToNode(el_i) = pnode ;
    }
  }
  
  
  /* Initialization */
  
  Elements_GetShapeFcts(elts) = NULL ;
  Elements_GetIntFcts(elts) = NULL ;
  
  {
    Element_t* el = Elements_GetElement(elts) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Element_t* el_i = el + i ;
      
      Element_GetElementIndex(el_i)      = i ;
      Element_GetDimension(el_i)         = -1 ;
      Element_GetNbOfNodes(el_i)         = 0 ;
      Element_GetRegionIndex(el_i)       = -1 ;
      Element_GetMaterial(el_i)          = NULL ;
      Element_GetMaterialIndex(el_i)     = -1 ;
      Element_GetShapeFct(el_i)          = NULL ;
      Element_GetIntFct(el_i)            = NULL ;
      Element_GetUnknownPosition(el_i)   = NULL ;
      Element_GetEquationPosition(el_i)  = NULL ;
      Element_GetBuffers(el_i)           = NULL ;
      //Element_GetElementSol(el_i)        = NULL ;
      Element_GetSolutions(el_i)         = NULL ;
    }
  }
  
  return(elts) ;
}



void (Elements_Delete)(void* self)
{
  Elements_t* elements = (Elements_t*) self ;
  
  {
    Element_t* el = Elements_GetElement(elements) ;
    
    {
      Node_t** pnode = Element_GetPointerToNode(el) ;
    
      if(pnode) {
        free(pnode) ;
      }
      
      Element_GetPointerToNode(el) = NULL ;
    }
    
    {
      short int* upos = Element_GetUnknownPosition(el) ;
    
      if(upos) {
        free(upos) ;
      }
      
      Element_GetUnknownPosition(el) = NULL ;
    }
    
    {
      Buffers_t* buf = Element_GetBuffers(el) ;
    
      if(buf) {
        Buffers_Delete(buf) ;
        free(buf) ;
      }
      
      Element_GetBuffers(el) = NULL ;
    }
    
    free(el) ;
  }
  
  {
    ShapeFcts_t* shapefcts = Elements_GetShapeFcts(elements) ;
    
    if(shapefcts) {
      ShapeFcts_Delete(shapefcts) ;
      free(shapefcts) ;
    }
  }
  
  {
    IntFcts_t* intfcts = Elements_GetIntFcts(elements) ;
    
    if(intfcts) {
      IntFcts_Delete(intfcts) ;
      free(intfcts) ;
    }
  }
}



void (Elements_LinkUp)(Elements_t* elements,Materials_t* materials)
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int n_mat = Materials_GetNbOfMaterials(materials) ;
  
  
  /* Link up element and material */
  {
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      int imat = Element_GetMaterialIndex(el + ie) ;
    
      if(imat >= 0 && imat < n_mat) {
        Element_GetMaterial(el + ie) = Materials_GetMaterial(materials) + imat ;
      } else {
        Element_GetMaterial(el + ie) = NULL ;
      }
    }
  }
}



void (Elements_CreateMore)(Elements_t* elements)
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  
  /* Pointers to unknowns and equations positions at nodes */
  {
    int n_pos = 0 ;
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      int imat = Element_GetMaterialIndex(el + ie) ;
      //Material_t* mat = Element_GetMaterial(el + ie) ;
    
      if(imat >= 0) {
        int nn = Element_GetNbOfNodes(el + ie) ;
        int neq = Element_GetNbOfEquations(el + ie) ;
        
        n_pos += nn*neq ;
      }
    }
  
    /* Memory space allocation with initialization to 0 */
    {
      //short int* pos = (short int* ) calloc(2*n_pos,sizeof(short int)) ;
      short int* upos = (short int* ) Mry_New(short int[2*n_pos]) ;
      short int* epos = upos + n_pos ;
      
      for(ie = 0 ; ie < 2*n_pos ; ie++) upos[ie] = 0 ;
      
      for(ie = 0 ; ie < n_el ; ie++) {
        int imat = Element_GetMaterialIndex(el + ie) ;
        //Material_t* mat = Element_GetMaterial(el + ie) ;
        
        Element_GetUnknownPosition(el + ie)  = upos ;
        Element_GetEquationPosition(el + ie) = epos ;
      
        if(imat >= 0) {
          int nn = Element_GetNbOfNodes(el + ie) ;
          int neq = Element_GetNbOfEquations(el + ie) ;
          
          upos += nn*neq ;
          epos += nn*neq ;
        }
      }
    }
  }

  /* Space allocation for buffer */
  {
    Buffers_t* buf = Buffers_Create(Element_SizeOfBuffer) ;
    int ie ;
  
    /* ATTENTION: same memory space (buffer) for all the elements */
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_GetBuffers(el + ie) = buf ;
    }
  }
  
  /* Create shape functions */
  {
    ShapeFcts_t* shapefcts = ShapeFcts_Create() ;
    int    ie ;
  
    Elements_GetShapeFcts(elements) = shapefcts ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* el_i = el + ie ;
      Material_t* mat = Element_GetMaterial(el_i) ;

      /* Find or create the shape functions */
      if(mat) {
        int  nn = Element_GetNbOfNodes(el_i) ;
        int  dim = Element_GetDimension(el_i) ;
        ShapeFct_t*  shapefct  = ShapeFcts_GetShapeFct(shapefcts) ;
        int  i = ShapeFcts_FindShapeFct(shapefcts,nn,dim) ;
      
        /* Element shape functions */
        Element_GetShapeFct(el_i) = shapefct + i ;
        
        {
          if(Element_HasZeroThickness(el_i)) {
            int nf = nn - Element_NbOfOverlappingNodes(el_i) ;
            int dim_h = dim - 1 ;
            int j  = ShapeFcts_FindShapeFct(shapefcts,nf,dim_h) ;

            Element_GetShapeFct(el_i) = shapefct + j ;
          }
        }
      }
    }
  }
  
  /* Create interpolation functions */
  {
    IntFcts_t* intfcts = IntFcts_Create() ;
    int    ie ;
    
    Elements_GetIntFcts(elements) = intfcts ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* el_i = el + ie ;
      Material_t* mat = Element_GetMaterial(el_i) ;

      /* Find or create default interpolation functions (Gauss type) */
      if(mat) {
        int  nn = Element_GetNbOfNodes(el_i) ;
        int  dim = Element_GetDimension(el_i) ;
        IntFct_t*  intfct  = IntFcts_GetIntFct(intfcts) ;
        int  i = IntFcts_FindIntFct(intfcts,nn,dim,"Gauss") ;

        /* Element interpolation functions */
        Element_GetIntFct(el_i)   = intfct + i ;
        
        {
          if(Element_HasZeroThickness(el_i)) {
            int nf = nn - Element_NbOfOverlappingNodes(el_i) ;
            int dim_h = dim - 1 ;
            int j  = IntFcts_FindIntFct(intfcts,nf,dim_h,"Gauss") ;

            Element_GetIntFct(el_i) = intfct + j ;
          }
        }
      }
    }
  }
  
  /* The max and min sizes of elements */
  Elements_GetMaximumSizeOfElements(elements) = Elements_ComputeMaximumSizeOfElements(elements) ;
  Elements_GetMinimumSizeOfElements(elements) = Elements_ComputeMinimumSizeOfElements(elements) ;
  
}



void  (Elements_DefineProperties)(Elements_t* elements)
/** Define the element properties including:
 **  1. the shape and interpolation functions
 **  2. the allocation memory of internal data (abstract data type)
 **  3. the size of tables for (im/ex)plicit and constant terms 
 **     for default allocation (see below) */
{
  int n_el = Elements_GetNbOfElements(elements) ;
  ShapeFcts_t* shapefcts = Elements_GetShapeFcts(elements) ;
  IntFcts_t* intfcts = Elements_GetIntFcts(elements) ;
  int    ie ;
  
  
  for(ie = 0 ; ie < n_el ; ie++) {
    Element_t* el = Elements_GetElement(elements) + ie ;
    Material_t* mat = Element_GetMaterial(el) ;

    /* Size of tables: initialization */
    Element_GetNbOfImplicitTerms(el) = 0 ;
    Element_GetNbOfExplicitTerms(el) = 0 ;
    Element_GetNbOfConstantTerms(el) = 0 ;

    if(mat) {
      /* Size of tables for (im/ex)plicit terms and constant terms
       * and possible other interpolation functions */
       
      /* To use interpolation functions with only part
       * of nodes (e.g. vertex nodes), cancel
       * equations and unknowns at the other nodes */
       
      /* To do so, proceed as follows in the model :
       * Element_GetUnknownPosition(el)  = -1
       * Element_GetEquationPosition(el) = -1 */
        
      Element_DefineProperties(el,intfcts,shapefcts) ;
      
      /* Update the size of tables in all the elementsol of the linked list */
      {
        int ni = Element_GetNbOfImplicitTerms(el) ;
        int ne = Element_GetNbOfExplicitTerms(el) ;
        int n0 = Element_GetNbOfConstantTerms(el) ;
        Solution_t* solution = Element_GetSolution(el) ;
  
        if(solution) {
          do {
            ElementSol_t* elementsol = Solution_GetElementSol(solution) + ie ;
            
            ElementSol_GetNbOfImplicitTerms(elementsol) = ni ;
            ElementSol_GetNbOfExplicitTerms(elementsol) = ne ;
            ElementSol_GetNbOfConstantTerms(elementsol) = n0 ;
            
            solution = Solution_GetPreviousSolution(solution) ;
          } while(solution != Element_GetSolution(el)) ;
        }
      }
    }
  }
  
  return ;
}



double (Elements_ComputeMaximumSizeOfElements)(Elements_t* elements)
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int ie ;
  double hmax = 0 ;
  
  for(ie = 0 ; ie < n_el ; ie++) {
    double h = Element_ComputeSize(el + ie) ;
    
    if(Element_IsSubmanifold(el + ie)) continue ;
    
    if(h > hmax) hmax = h ;
  }
  
  return(hmax) ;
}



double (Elements_ComputeMinimumSizeOfElements)(Elements_t* elements)
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int ie ;
  double hmin = -1 ;
  
  for(ie = 0 ; ie < n_el ; ie++) {
    double h = Element_ComputeSize(el + ie) ;
    
    if(Element_IsSubmanifold(el + ie)) continue ;
    
    if(hmin < 0 && h > 0) hmin = h ;
    
    if(h < hmin) hmin = h ;
  }
  
  return(hmin) ;
}




int (Elements_ComputeNbOfMatrixEntries)(Elements_t* elements)
{
  int nel = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int len = 0 ;
      
  {
    int ie ;
      
    for(ie = 0 ; ie < nel ; ie++) {
      Element_FreeBuffer(el + ie) ;
      len += Element_ComputeNbOfMatrixEntries(el + ie) ;
    }
  }
  
  return(len) ;
}




int (Elements_ComputeNbOfSelectedMatrixEntries)(Elements_t* elements,const int imatrix)
{
  int nel = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int len = 0 ;
      
  {
    int ie ;
      
    for(ie = 0 ; ie < nel ; ie++) {
      Element_FreeBuffer(el + ie) ;
      len += Element_ComputeNbOfSelectedMatrixEntries(el + ie,imatrix) ;
    }
  }
  
  return(len) ;
}



void  (Elements_InitializeMatrixRowColumnIndexes)(Elements_t* elements)
/** Initialize the Matrix Row/Column Indexes to 0 
 *  for nodes belonging to elements */
{
  {
    Element_t* el = Elements_GetElement(elements) ;
    int n_el = Elements_GetNbOfElements(elements) ;
    int    ie ;
  
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* elt_i = el + ie ;
      Material_t* mat = Element_GetMaterial(elt_i) ;
    
      if(mat) {
        int    nn  = Element_GetNbOfNodes(elt_i) ;
        int    neq = Element_GetNbOfEquations(elt_i) ;
        int i ;

        for(i = 0 ; i < nn ; i++) {
          Node_t* node_i = Element_GetNode(elt_i,i) ;
          int   j ;

          for(j = 0 ; j < neq ; j++) {
            int ij = i*neq + j ;
          
            /*  columns (unknowns) set to arbitrary >= 0 */
            {
              int ii = Element_GetUnknownPosition(elt_i)[ij] ;
              
              if(ii >= 0) Node_GetMatrixColumnIndex(node_i)[ii] = 0 ;
            }
          
            /*  rows (equations) set to arbitrary >= 0 */
            {
              int ii = Element_GetEquationPosition(elt_i)[ij] ;
              
              if(ii >= 0) Node_GetMatrixRowIndex(node_i)[ii] = 0 ;
            }
          }
        }
      }
    }
  }
}




void  (Elements_EliminateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t* elements)
/** Elimninate indexes of specific unknowns at the overlapping nodes 
 *  of zero-thickness interface element.
 **/
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int ie ;
    
  for(ie = 0 ; ie < n_el ; ie++) {
    Element_t* el_i = el + ie ;
      
    if(Element_HasZeroThickness(el_i)) {
      int neq = Element_GetNbOfEquations(el_i) ;
      ShapeFct_t* shapefct = Element_GetShapeFct(el_i) ;
      int nf = ShapeFct_GetNbOfFunctions(shapefct) ;
      int in ;
      
      for(in = 0 ; in < nf ; in++) {
        Node_t* node_in = Element_GetNode(el_i,in) ;
        int ieq ;
    
        for(ieq = 0 ; ieq < neq ; ieq++) {
          char* unk = Element_GetNameOfUnknown(el_i)[ieq] ;
          char* eqn = Element_GetNameOfEquation(el_i)[ieq] ;
          
          Node_EliminateMatrixColumnIndexForOverlappingNodes(node_in,unk) ;
          Node_EliminateMatrixRowIndexForOverlappingNodes(node_in,eqn) ;
        }
      }
    }
  }
}





void  (Elements_UpdateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t* elements)
/** Merge indexes of specific unknowns at the overlapping nodes 
 *  of zero-thickness interface element.
 **/
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int ie ;
    
  for(ie = 0 ; ie < n_el ; ie++) {
    Element_t* el_i = el + ie ;
      
    if(Element_HasZeroThickness(el_i)) {
      int neq = Element_GetNbOfEquations(el_i) ;
      ShapeFct_t* shapefct = Element_GetShapeFct(el_i) ;
      int nf = ShapeFct_GetNbOfFunctions(shapefct) ;
      int in ;
      
      for(in = 0 ; in < nf ; in++) {
        Node_t* node_in = Element_GetNode(el_i,in) ;
        int ieq ;
    
        for(ieq = 0 ; ieq < neq ; ieq++) {
          char* unk = Element_GetNameOfUnknown(el_i)[ieq] ;
          char* eqn = Element_GetNameOfEquation(el_i)[ieq] ;
          
          Node_UpdateMatrixColumnIndexForOverlappingNodes(node_in,unk) ;
          Node_UpdateMatrixRowIndexForOverlappingNodes(node_in,eqn) ;
        }
      }
    }
  }
}

