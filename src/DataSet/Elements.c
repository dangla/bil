#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "Buffer.h"
#include "Geometry.h"
#include "Nodes.h"
#include "IntFcts.h"
#include "Materials.h"
#include "Elements.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Curves.h"


static double   (Elements_ComputeMaximumSizeOfElements)(Elements_t*) ;
static double   (Elements_ComputeMinimumSizeOfElements)(Elements_t*) ;
static double   (Element_ComputeSize)(Element_t*) ;
static int      (Element_ComputeNbOfMatrixEntries)(Element_t*) ;



void Elements_CreateMore(Elements_t* elements,Materials_t* materials)
{
  int n_el = Elements_GetNbOfElements(elements) ;
  Element_t* el = Elements_GetElement(elements) ;
  int ie ;
  
  /* Pointers to material for each element */
  {
    for(ie = 0 ; ie < n_el ; ie++) {
      int imat = Element_GetMaterialIndex(el + ie) ;
    
      if(imat >= 0) {
        Element_GetMaterial(el + ie) = Materials_GetMaterial(materials) + imat ;
      } else {
        Element_GetMaterial(el + ie) = NULL ;
      }
    }
  }
  
  /* Pointers to unknowns and equations positions at nodes */
  {
    int n_pn = 0 ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Material_t* mat = Element_GetMaterial(el + ie) ;
      int nn = Element_GetNbOfNodes(el + ie) ;
      int neq = Element_GetNbOfEquations(el + ie) ;
    
      if(mat) n_pn += nn*neq ;
    }
  
  /* Memory space allocation with initialization to 0 */
    {
      short int* pntr_pn = (short int* ) calloc(2*n_pn,sizeof(short int)) ;
    
      if(!pntr_pn) {
        arret("Elements_CreateMore (1) : impossible d\'allouer la memoire") ;
      }
  
      Element_GetUnknownPosition(el)  = pntr_pn ;
      Element_GetEquationPosition(el) = pntr_pn + n_pn ;
    
      for(ie = 1 ; ie < n_el ; ie++) {
        Material_t* mat = Element_GetMaterial(el + ie) ;
        short int* pin = Element_GetUnknownPosition(el + ie - 1) ;
        short int* peq = Element_GetEquationPosition(el + ie - 1) ;
      
        if(mat) {
          int nn = Element_GetNbOfNodes(el + ie - 1) ;
          int neq = Element_GetNbOfEquations(el + ie - 1) ;
        
          Element_GetUnknownPosition(el + ie) = pin + nn*neq ;
          Element_GetEquationPosition(el + ie) = peq + nn*neq ;
        
        } else {
        
          Element_GetUnknownPosition(el + ie) = pin ;
          Element_GetEquationPosition(el + ie) = peq ;
        }
      }
    }
  }

  /* Space allocation for buffer */
  {
    Buffer_t* buf = Buffer_Create(Element_SizeOfBuffer) ;
  
    /* ATTENTION : same memory space (buffer) for all the elements */
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_GetBuffer(el + ie) = buf ;
    }
  }
  
  /* Create shape functions */
  Elements_GetShapeFcts(elements) = ShapeFcts_Create() ;
  
  /* Create interpolation functions */
  Elements_GetIntFcts(elements) = IntFcts_Create() ;
  
  /* The max and min sizes of elements */
  Elements_GetMaximumSizeOfElements(elements) = Elements_ComputeMaximumSizeOfElements(elements) ;
  Elements_GetMinimumSizeOfElements(elements) = Elements_ComputeMinimumSizeOfElements(elements) ;
  
}



void  Elements_DefineProperties(Elements_t* elements)
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


  /* Nb of shape functions: initialization */
  ShapeFcts_GetNbOfShapeFcts(shapefcts) = 0 ;

  /* Nb of interpolation functions: initialization */
  IntFcts_GetNbOfIntFcts(intfcts) = 0 ;
  
  
  for(ie = 0 ; ie < n_el ; ie++) {
    Element_t* el = Elements_GetElement(elements) + ie ;
    Material_t* mat = Element_GetMaterial(el) ;

    /* Size of tables: initialization */
    Element_GetNbOfImplicitTerms(el) = 0 ;
    Element_GetNbOfExplicitTerms(el) = 0 ;
    Element_GetNbOfConstantTerms(el) = 0 ;

    if(mat) {
      int  nn = Element_GetNbOfNodes(el) ;
      int  dim = Element_GetDimension(el) ;
      
      /* Find or create the shape functions */
      {
        ShapeFct_t*  shapefct  = ShapeFcts_GetShapeFct(shapefcts) ;
        int  i = ShapeFcts_FindShapeFct(shapefcts,nn,dim) ;
      
        /* Element shape functions */
        Element_GetShapeFct(el) = shapefct + i ;
      }
      
      /* Find or create default interpolation functions (Gauss type) */
      {
        IntFct_t*  intfct  = IntFcts_GetIntFct(intfcts) ;
        int  i = IntFcts_FindIntFct(intfcts,nn,dim,"Gauss") ;
      
        /* Element interpolation functions */
        Element_GetIntFct(el)   = intfct + i ;
      }
      
      /* Size of tables for (im/ex)plicit terms and constant terms
       * and possible other interpolation functions */
       
      /* To use interpolation functions with only part
       * of nodes (e.g. vertex nodes), cancel
       * equations and unknowns at the other nodes */
       
      /* To do so, proceed as follows in the model :
       * Element_GetUnknownPosition(el)  = -1
       * Element_GetEquationPosition(el) = -1 */
        
      Element_DefineProperties(el,intfcts) ;
      
      /* Update the size of tables in all the elementsol of the linked list */
      {
        int ni = Element_GetNbOfImplicitTerms(el) ;
        int ne = Element_GetNbOfExplicitTerms(el) ;
        int n0 = Element_GetNbOfConstantTerms(el) ;
        ElementSol_t* elementsol = Element_GetElementSol(el) ;
  
        if(elementsol) {
          do {
            ElementSol_GetNbOfImplicitTerms(elementsol) = ni ;
            ElementSol_GetNbOfExplicitTerms(elementsol) = ne ;
            ElementSol_GetNbOfConstantTerms(elementsol) = n0 ;
            elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
          } while(elementsol != Element_GetElementSol(el)) ;
        }
      }
    }
  }
  
  return ;
}



double Elements_ComputeMaximumSizeOfElements(Elements_t* elements)
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



double Elements_ComputeMinimumSizeOfElements(Elements_t* elements)
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



double Element_ComputeSize(Element_t* element)
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  double h = 0 ;
  double c[3] = {0,0,0} ;
  int i ;
  
  /* The center of element */
  for(i = 0 ; i < nn ; i++) {
    double* x = Element_GetNodeCoordinate(element,i) ;
    int j ;
    
    for(j = 0 ; j < dim ; j++) {
      c[j] += x[j]/nn ;
    }
  }
  
  /* The "radius" of element */
  for(i = 0 ; i < nn ; i++) {
    double* x = Element_GetNodeCoordinate(element,i) ;
    double r = 0 ;
    int j ;
    
    for(j = 0 ; j < dim ; j++) {
      double y = x[j] - c[j] ;
      
      r += y*y ;
    }
    
    r = sqrt(r) ;
    
    if(r > h) h = r ;
  }
  
  h *= 2 ;
  
  return(h) ;
}




int Elements_ComputeNbOfMatrixEntries(Elements_t* elements)
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




int Element_ComputeNbOfMatrixEntries(Element_t* element)
{
  int   ndof = Element_GetNbOfDOF(element) ;
  int*  row  = Element_ComputeMatrixRowAndColumnIndices(element) ;
  int*  col  = row + ndof ;
  int   len  = 0 ;
  
  {
    int   jdof ;

    for(jdof = 0 ; jdof < ndof ; jdof++) {
      int jcol = col[jdof] ;
    
      if(jcol < 0) continue ;

      {
        int idof ;
            
        for(idof = 0 ; idof < ndof ; idof++) {
          int irow = row[idof] ;
      
          if(irow < 0) continue ;
        
          len++ ;
        }
      }
    }
  }
  
  return(len) ;
}
