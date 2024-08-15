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


#if 1
Elements_t*  (Elements_New)(const int n,const int nc)
{
  Elements_t* elts = (Elements_t*) Mry_New(Elements_t) ;
  
  Elements_GetNbOfElements(elts) = n ;
  Elements_GetNbOfConnectivities(elts) = nc ;
  
  
  /* Allocation of space for the pointers to "node" */
  {
    Node_t** pnode = (Node_t**) Mry_New(Node_t*[nc]) ;
    
    Elements_GetPointerToNode(elts) = pnode ;
  }
  
  
  /* Allocation of space for each element */
  {
    Element_t* el = Mry_Create(Element_t,n,Element_New()) ;
    
    Elements_GetElement(elts) = el ;
  }
  
  
  /* Allocation of space for the regions */
  {
    Regions_t* regions = Regions_New() ;
    
    Elements_GetRegions(elts) = regions ;
  }
  
  
  /* Initialization */
  {
    Node_t** pnode = Elements_GetPointerToNode(elts) ;
    Element_t* el = Elements_GetElement(elts) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Element_t* el_i = el + i ;
      
      Element_GetElementIndex(el_i)  = i ;
      Element_GetPointerToNode(el_i) = pnode ;
    }
  }
  
  Elements_GetShapeFcts(elts) = NULL ;
  Elements_GetIntFcts(elts) = NULL ;
  
  return(elts) ;
}
#endif


#if 1
void (Elements_CreateMore)(Elements_t* elements)
{
  Buffers_t* buffers = Buffers_Create(Element_SizeOfBuffer) ;
  ShapeFcts_t* shapefcts = ShapeFcts_Create() ;
  IntFcts_t* intfcts = IntFcts_Create() ;
  
  {
    Elements_GetBuffers(elements) = buffers ;
    Elements_GetShapeFcts(elements) = shapefcts ;
    Elements_GetIntFcts(elements) = intfcts ;
  }
  
  /* Create additional memory space for each element */
  {
    int n_el = Elements_GetNbOfElements(elements) ;
    Element_t* el = Elements_GetElement(elements) ;
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_CreateMore(el + ie,buffers,shapefcts,intfcts) ;
    }
  }
  
  /* The max and min sizes of elements */
  Elements_GetMaximumSizeOfElements(elements) = Elements_ComputeMaximumSizeOfElements(elements) ;
  Elements_GetMinimumSizeOfElements(elements) = Elements_ComputeMinimumSizeOfElements(elements) ;
  
}
#endif



void (Elements_Delete)(void* self)
{
  Elements_t* elements = (Elements_t*) self ;
  
    
  {
    Node_t** pnode = Elements_GetPointerToNode(elements) ;
    
    if(pnode) {
      free(pnode) ;
    }
  }

  {
    int nel = Elements_GetNbOfElements(elements) ;
    Element_t* el = Elements_GetElement(elements) ;
      
    Mry_Delete(el,nel,Element_Delete) ;
    free(el) ;
  }

  {
    Regions_t* regions = Elements_GetRegions(elements) ;
    
    if(regions) {
      Regions_Delete(regions) ;
      free(regions) ;
    }
  }
    
  {
    Buffers_t* buf = Elements_GetBuffers(elements) ;
    
    if(buf) {
      Buffers_Delete(buf) ;
      free(buf) ;
    }
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
       
      /* To do so, proceed as follows in the model
       * (to cancel unknown/equation j at node i):
       * Element_GetUnknownPosition(el)[i*neq+j] = -1
       * Element_GetEquationPosition(el)[i*neq+j] = -1 */
        
      Element_DefineProperties(el,intfcts,shapefcts) ;
      
      /* Update the size of tables in all the elementsol of the circularly linked list */
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

