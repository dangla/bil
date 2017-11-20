#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "Node.h"
#include "Element.h"
#include "Tools/Math.h"
#include "Message.h"


#if 0
Element_AllocateGenericData(Element_t* el,int n,void* gdat,TypeId_t typeid,)
{
  

  /** The solutions from the microstructures */
  {
    /* Store "sols" for the whole history */
    {
      ElementSol_t* elementsol = Element_GetElementSol(el) ;

      if(elementsol) {
        do {
          Solutions_t* sols = (Solutions_t*) malloc(NbOfIntPoints*sizeof(Solutions_t)) ;
  
          if(!sols) arret("DefineElementProp") ;
    
          {
            int nsol_micro = 2 ;
            Mesh_t* mesh = DataSet_GetMesh(dataset) ;
            int i ;
    
            for(i = 0 ; i < NbOfIntPoints ; i++) {
              Solutions_t* solsi = Solutions_Create(mesh,nsol_micro) ;

              /* Not essential */
              Solutions_MergeExplicitTerms(solsi) ;
        
              sols[i] = *solsi ;
            }
          }
          
          {
            GenericData_t* gdat0 = ElementSol_GetImplicitGenericData(elementsol) ;
            GenericData_t* gdat  = GenericData_New() ;
          
            GenericData_Initialize(gdat,NbOfIntPoints,sols,Solutions_t) ;
          
            GenericData_InsertAfter(gdat0,gdat) ;
          }
        
          elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
        } while(elementsol != Element_GetElementSol(el)) ;
      }
    }
  }
}
#endif



double** Element_ComputePointerToNodalCoordinates(Element_t* element)
/** Compute the nodal coordinates */
{
  int nn = Element_GetNbOfNodes(element) ;
  size_t SizeNeeded = nn*sizeof(double*) ;
  double** x = (double**) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(element,i) ;
  }
  
  return(x) ;
}



double* Element_ComputeNodalCoordinates(Element_t* element)
/** Compute the nodal coordinates */
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  size_t SizeNeeded = nn*dim*sizeof(double) ;
  double* x = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    double* xi = Element_GetNodeCoordinate(element,i) ;
    int j ;
    
    for(j = 0 ; j < dim ; j++) {
      x[i*dim + j] = xi[j] ;
    }
  }
  
  return(x) ;
}



double** Element_ComputePointerToCurrentNodalUnknowns(Element_t* element)
/** Compute a pointer to the nodal unknowns at the current time */
{
  int nn = Element_GetNbOfNodes(element) ;
  size_t SizeNeeded = nn*(sizeof(double*)) ;
  double** u = (double**) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetCurrentNodalUnknown(element,i) ;
  }

  return(u) ;
}



double** Element_ComputePointerToPreviousNodalUnknowns(Element_t* element)
/** Compute a pointer to the nodal unknowns at the current time */
{
  int nn = Element_GetNbOfNodes(element) ;
  size_t SizeNeeded = nn*(sizeof(double*)) ;
  double** u = (double**) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetPreviousNodalUnknown(element,i) ;
  }

  return(u) ;
}



double* Element_ComputeDeepNodalUnknowns(Element_t* element,unsigned int depth)
/** Compute the nodal unknowns at the some depth */
{
  int nn = Element_GetNbOfNodes(element) ;
  int neq = Element_GetNbOfEquations(element) ;
  size_t SizeNeeded = nn*(neq*sizeof(double)) ;
  double* u = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    Node_t* node = Element_GetNode(element,i) ;
    double* v = Node_GetDeepUnknown(node,depth) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(element,i,j) ;
      if(jj >= 0) {
        u[i*neq + j] = v[jj] ;
      } else {
        u[i*neq + j] = 0. ;
      }
    }
  }

  return(u) ;
}



double* Element_ComputeCurrentNodalUnknowns(Element_t* element)
/** Compute the nodal unknowns at the current time */
{
  int nn = Element_GetNbOfNodes(element) ;
  int neq = Element_GetNbOfEquations(element) ;
  size_t SizeNeeded = nn*neq*sizeof(double) ;
  double* u = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    double* v = Element_GetCurrentNodalUnknown(element,i) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(element,i,j) ;
      
      if(jj >= 0) {
        u[i*neq + j] = v[jj] ;
      } else {
        u[i*neq + j] = 0. ;
      }
    }
  }

  return(u) ;
}



double* Element_ComputePreviousNodalUnknowns(Element_t* element)
/** Compute the nodal unknowns at the previous time */
{
  int nn = Element_GetNbOfNodes(element) ;
  int neq = Element_GetNbOfEquations(element) ;
  size_t SizeNeeded = nn*neq*sizeof(double) ;
  double* u = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    double* v = Element_GetPreviousNodalUnknown(element,i) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(element,i,j) ;
      
      if(jj >= 0) {
        u[i*neq + j] = v[jj] ;
      } else {
        u[i*neq + j] = 0. ;
      }
    }
  }

  return(u) ;
}



double* Element_ComputeIncrementalNodalUnknowns(Element_t* element)
/** Compute the incremental values of nodal unknowns */
{
  int nn = Element_GetNbOfNodes(element) ;
  int neq = Element_GetNbOfEquations(element) ;
  size_t SizeNeeded = nn*neq*sizeof(double) ;
  double* u = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    double* v1 = Element_GetCurrentNodalUnknown(element,i) ;
    double* vn = Element_GetPreviousNodalUnknown(element,i) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(element,i,j) ;
      
      if(jj >= 0) {
        u[i*neq + j] = v1[jj] - vn[jj] ;
      } else {
        u[i*neq + j] = 0. ;
      }
    }
  }

  return(u) ;
}



double* Element_ComputeIncrementalImplicitTerms(Element_t* element)
{
  int n = Element_GetNbOfImplicitTerms(element) ;
  double* vi1 = Element_GetCurrentImplicitTerm(element) ;
  double* vin = Element_GetPreviousImplicitTerm(element) ;
  size_t SizeNeeded = n*sizeof(double) ;
  double* vii = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int i ;
  
  for(i = 0 ; i < n ; i++) vii[i] = vi1[i] - vin[i] ;
  
  return(vii) ;
}



int    Element_FindUnknownPositionIndex(Element_t* element,char* s)
/** Find the unknown position index whose name is pointed to by s */
{
  int n = Element_GetNbOfEquations(element) ;
  char** ss = Element_GetNameOfUnknown(element) ;
  int    i ;

  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) if(!strncmp(s,ss[i],strlen(s))) break ;
    if(i == n) i = -1 ;
  }

  if(i < 0) arret("Element_FindUnknownPositionIndex (1) : position not known") ;
  return(i) ;
}



int    Element_FindEquationPositionIndex(Element_t* element,char* s)
/** Find the equation position index whose name is pointed to by s */
{
  int n = Element_GetNbOfEquations(element) ;
  char** ss = Element_GetNameOfEquation(element) ;
  int    i ;

  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) if(!strncmp(s,ss[i],strlen(s))) break ;
    if(i == n) i = -1 ;
  }

  if(i < 0) arret("Element_FindEquationPositionIndex (1) : position not known") ;
  return(i) ;
}




double*  Element_ComputeNormalVector(Element_t* element,double* dh,int nn)
/** Compute the unit outward normal vector to a submanifold element
 *  of dimension dim-1 */
{
#define DH(n,i) (dh[(n)*dim_h+(i)])
  size_t SizeNeeded = 3*sizeof(double) ;
  double* norm = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    dim_h = Element_GetDimension(element) ;
  int    dim   = Element_GetDimensionOfSpace(element) ;
  double c[3][3] ;
  
  
  if(dim_h != dim - 1) {
    arret("Element_ComputeNormalVector") ;
  }
  
  
  /* The jacobian */
  {
    //int nn = Element_GetNbOfNodes(element) ;
    int    i,j ;
    
    for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
      int    k ;
      
      c[i][j] = 0. ;
      for(k = 0 ; k < nn ; k++) {
        double* x = Element_GetNodeCoordinate(element,k) ;
        
        c[i][j] += x[i]*DH(k,j) ;
      }
    }
  }
  
  /* 1. Surface : norm = c1^c2 */
  if(dim_h == 2) {
    double v ;
    
    c[0][2] = c[1][0]*c[2][1] - c[2][0]*c[1][1] ;
    c[1][2] = c[2][0]*c[0][1] - c[0][0]*c[2][1] ;
    c[2][2] = c[0][0]*c[1][1] - c[1][0]*c[0][1] ;
    v = sqrt(c[0][2]*c[0][2] + c[1][2]*c[1][2] + c[2][2]*c[2][2]) ;
    norm[0] = c[0][2]/v ;
    norm[1] = c[1][2]/v ;
    norm[2] = c[2][2]/v ;
    
    return(norm) ;

  /* 2. Ligne : norm = c1^e_z */
  } else if(dim_h == 1) {
    double v = sqrt(c[0][0]*c[0][0] + c[1][0]*c[1][0]) ;
    
    norm[0] =  c[1][0]/v ;
    norm[1] = -c[0][0]/v ;
    norm[2] = 0. ;
    
    return(norm) ;

  /* 3. Point : norm = 1 */
  } else if(dim_h == 0) {
    norm[0] = 1. ;
    norm[1] = 0. ;
    norm[2] = 0. ;
    
    return(norm) ;
  }

  arret("Element_ComputeNormalVector") ;
  return(norm) ;
#undef DH
}



double*  Element_ComputeCoordinateInReferenceFrame(Element_t* el,double* x)
/** Compute the local coordinates in the reference element which maps 
 *  into coordinates "x" in the actual space */
{
#define DH(n,i)  (dh[(n)*dim_e+(i)])
  unsigned short int dim   = Element_GetDimensionOfSpace(el) ;
  unsigned short int dim_e = Element_GetDimension(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double* x_e[Element_MaxNbOfNodes] ;
  double diameter ;
  
  
  /* Node coordinates */
  {
    int    i ;
    
    for(i = 0 ; i < nn ; i++) {
      x_e[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }


  /* Diameter of element */
  {
    double x_max[3],x_min[3] ;
    int    i ;
    
    diameter = 0. ;
    for(i = 0 ; i < dim ; i++) {
      int   in ;
    
      x_max[i] = (x_min[i] = x_e[0][i]) ;
    
      for(in = 1 ; in < nn ; in++) {
        x_max[i] = (x_e[in][i] > x_max[i]) ? x_e[in][i] : x_max[i] ;
        x_min[i] = (x_e[in][i] < x_min[i]) ? x_e[in][i] : x_min[i] ;
      }
    
      x_max[i] -= x_min[i] ;
      diameter = (x_max[i] > diameter) ? x_max[i] : diameter ;
    }
  }


  /* Compute the coordinates in the reference element */
  {
    ShapeFct_t* shapefct = Element_GetShapeFct(el) ;
    double* a = ShapeFct_GetCoordinate(shapefct) ;
    int    max_iter = 10 ;
    double tol = 1.e-6 ;
    int    iter ;
    double err ;
    int i ;
    
    for(i = 0 ; i < 3 ; i++) {
      a[i] = 0 ;
    }
    
    for(iter = 0 ; iter < max_iter ; iter++) {
      double* h = ShapeFct_GetFunction(shapefct) ;
      double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
      double r[3] ;
      double k[9] ;
      
      ShapeFct_ComputeValuesAtPoint(dim_e,nn,a,h,dh) ;
    
      for(i = 0 ; i < dim ; i++) {
        int   j,in ;
      
        r[i] = x[i] ;
        for(j = 0 ; j < dim ; j++) k[dim*i+j] = 0. ;
      
        for(in = 0 ; in < nn ; in++) {
          r[i] -= h[in]*x_e[in][i] ;
          for(j = 0 ; j < dim_e ; j++) k[dim*i+j] += DH(in,j)*x_e[in][i] ;
        }
      }
    
      if(dim_e < dim) {
        if(dim_e == dim-1) {
          int  j = dim - 1 ;
          
          double* unorm = Element_ComputeNormalVector(el,dh,nn) ;
        
          for(i = 0 ; i < dim ; i++) k[dim*i + j] = unorm[i] ;
          
        } else {
          arret("Element_ComputeCoordinateInReferenceFrame(2)") ; /* a traiter plus tard */
        }
      }
    
      for(i = 0 ; i < dim ; i++) {
        int  j ;
      
        r[i] /= diameter ;
        
        for(j = 0 ; j < dim ; j++) k[dim*i+j] /= diameter ;
      }
    
      Math_SolveByGaussElimination(k,r,dim) ;
    
      for(i = 0 ; i < dim ; i++) a[i] += r[i] ;
    
      err = 0 ;
      for(i = 0 ; i < dim ; i++) {
        if(fabs(r[i]) > err) err = fabs(r[i]) ;
      }
    
      if(err < tol) break ;
    }
  
    if(err > tol) {
      arret("Element_ComputeCoordinateInReferenceFrame(3)") ;
    }
    
    return(a) ;
  }
  
  return(NULL) ;
#undef DH
}



int Element_ComputeNbOfSolutions(Element_t* el)
{
  ElementSol_t* elementsol = Element_GetElementSol(el) ;
  int n = 0 ;
  
  if(elementsol) {
    do {
      n++ ;
      elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
    } while(elementsol != Element_GetElementSol(el)) ;
  }
  
  return(n) ;
}
