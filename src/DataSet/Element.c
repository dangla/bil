#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "String_.h"
#include "Mry.h"
#include "Node.h"
#include "Element.h"
#include "Math_.h"
#include "Message.h"





Element_t*  (Element_New)(void)
{
  Element_t* el = (Element_t*) Mry_New(Element_t) ;
  
  {
    Element_GetElementIndex(el)      = -1 ;
    Element_GetPointerToNode(el)     = NULL ;
    Element_GetDimension(el)         = -1 ;
    Element_GetNbOfNodes(el)         = 0 ;
    //Element_GetRegionTag(el)       = -1 ;
    Element_GetRegion(el)            = NULL ;
    Element_GetMaterial(el)          = NULL ;
    Element_GetMaterialIndex(el)     = -1 ;
    Element_GetShapeFct(el)          = NULL ;
    Element_GetIntFct(el)            = NULL ;
    Element_GetUnknownPosition(el)   = NULL ;
    Element_GetEquationPosition(el)  = NULL ;
    Element_GetBuffers(el)           = NULL ;
    Element_GetSolutions(el)         = NULL ;
  }
  
  return(el) ;
}



void (Element_CreateMore)(Element_t* el,Buffers_t* buffers,ShapeFcts_t* shapefcts,IntFcts_t* intfcts)
{
  int imat = Element_GetMaterialIndex(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int neq = Element_GetNbOfEquations(el) ;
  int dim = Element_GetDimension(el) ;
  int ndof = nn*neq ;
  

  /* Memory space allocation for the pointers to unknowns and 
   * equations positions at nodes with initialization to 0 */
  if(imat >= 0) {
    short int* upos = (short int* ) Mry_New(short int[2*ndof]) ;
    short int* epos = upos + ndof ;
    int i ;
      
    for(i = 0 ; i < 2*ndof ; i++) upos[i] = 0 ;

    Element_GetUnknownPosition(el)  = upos ;
    Element_GetEquationPosition(el) = epos ;
  } else {
    Element_GetUnknownPosition(el)  = NULL ;
    Element_GetEquationPosition(el) = NULL ;
  }
  
  
  /* Memory space allocation for the matrix */
  {
    double* matrix = (double*) Mry_New(double[ndof*ndof]) ;
    
    Element_GetMatrix(el) = matrix ;
  }
  
  
  /* Memory space allocation for the residu */
  {
    double* residu = (double*) Mry_New(double[ndof]) ;
    
    Element_GetResidu(el) = residu ;
  }


  /* The buffer (the same for all elements) */
  {
    Element_GetBuffers(el) = buffers ;
  }


  /* Find or create the shape functions */
  {
    Material_t* mat = Element_GetMaterial(el) ;

    if(mat) {
      ShapeFct_t*  shapefct  = ShapeFcts_GetShapeFct(shapefcts) ;
      int  i = ShapeFcts_FindShapeFct(shapefcts,nn,dim) ;
      
      /* Element shape functions */
      Element_GetShapeFct(el) = shapefct + i ;
        
      {
        if(Element_HasZeroThickness(el)) {
          int nf = nn - Element_NbOfOverlappingNodes(el) ;
          int dim_h = dim - 1 ;
          int j  = ShapeFcts_FindShapeFct(shapefcts,nf,dim_h) ;

          Element_GetShapeFct(el) = shapefct + j ;
        }
      }
    }
  }


  /* Find or create default interpolation functions (Gauss type) */
  {
    Material_t* mat = Element_GetMaterial(el) ;

    if(mat) {
      IntFct_t*  intfct  = IntFcts_GetIntFct(intfcts) ;
      int  i = IntFcts_FindIntFct(intfcts,nn,dim,"Gauss") ;

      /* Element interpolation functions */
      Element_GetIntFct(el)   = intfct + i ;
        
      {
        if(Element_HasZeroThickness(el)) {
          int nf = nn - Element_NbOfOverlappingNodes(el) ;
          int dim_h = dim - 1 ;
          int j  = IntFcts_FindIntFct(intfcts,nf,dim_h,"Gauss") ;

          Element_GetIntFct(el) = intfct + j ;
        }
      }
    }
  }
}



void (Element_Delete)(void* self)
{
  Element_t* el = (Element_t*) self ;
  
  {
    short int* upos = Element_GetUnknownPosition(el) ;
    
    if(upos) {
      free(upos) ;
    }
      
    Element_GetUnknownPosition(el) = NULL ;
  }
  
  {
    double* matrix = Element_GetMatrix(el) ;
    
    if(matrix) {
      free(matrix) ;
    }
    
    Element_GetMatrix(el) = NULL ;
  }
  
  {
    double* residu = Element_GetResidu(el) ;
    
    if(residu) {
      free(residu) ;
    }
    
    Element_GetResidu(el) = NULL ;
  }
}



void (Element_AllocateMicrostructureSolutions)(Element_t const* el,Mesh_t* mesh,const int nsol)
/** Allocate space as implicit generic data at the interpolation points 
 *  of the element "el" for the solutions of the microstructure defined by "mesh".
 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  {
    int ie = Element_GetElementIndex(el) ;
    Solution_t* solution = Element_GetSolution(el) ;
      
    /* Store "sols" for the whole history */
    if(solution) {
      do {
        Solutions_t* sols = (Solutions_t*) Mry_New(Solutions_t[NbOfIntPoints]) ;
          
        {
          int i ;
    
          for(i = 0 ; i < NbOfIntPoints ; i++) {
            Solutions_t* solsi = Solutions_Create(mesh,nsol) ;
              
            sols[i] = solsi[0] ;
            free(solsi) ;
          }
        }
          
        {
          ElementSol_t* elementsol = Solution_GetElementSol(solution) + ie ;
          GenericData_t* gdat  = GenericData_Create(NbOfIntPoints,sols,Solutions_t,"Solutions") ;
          
          ElementSol_AddImplicitGenericData(elementsol,gdat) ;
        }

        solution = Solution_GetPreviousSolution(solution) ;
      } while(solution != Element_GetSolution(el)) ;
    }
  }
}



double** (Element_ComputePointerToNodalCoordinates)(Element_t const* element)
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



double* (Element_ComputeNodalCoordinates)(Element_t const* element)
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



double** (Element_ComputePointerToCurrentNodalUnknowns)(Element_t const* element)
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



double** (Element_ComputePointerToPreviousNodalUnknowns)(Element_t const* element)
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



double* (Element_ComputeDeepNodalUnknowns)(Element_t const* element,unsigned int depth)
/** Compute the nodal unknowns at the some depth */
{
  int nn = Element_GetNbOfNodes(element) ;
  int neq = Element_GetNbOfEquations(element) ;
  size_t SizeNeeded = nn*(neq*sizeof(double)) ;
  double* u = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    Node_t* node = Element_GetNode(element,i) ;
    double* v = Node_GetUnknownInDistantPast(node,depth) ;
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



double* (Element_ComputeCurrentNodalUnknowns)(Element_t const* element)
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



double* (Element_ComputePreviousNodalUnknowns)(Element_t const* element)
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



double* (Element_ComputeIncrementalNodalUnknowns)(Element_t const* element)
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



double* (Element_ComputeIncrementalImplicitTerms)(Element_t const* element)
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



int    (Element_FindUnknownPositionIndex)(Element_t const* element,const char* s)
/** Find the unknown position index whose name is pointed to by s */
{
  int n = Element_GetNbOfEquations(element) ;
  const char* const* ss = (const char* const*) Element_GetNameOfUnknown(element) ;
  int    i = String_FindPositionIndex(s,ss,n) ;
  
  return(i) ;
#if 0
  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) {
      if(!strncmp(s,ss[i],strlen(s))) break ;
    }
    if(i == n) i = -1 ;
  }

  //if(i < 0) arret("Element_FindUnknownPositionIndex: unknown position of %s",s) ;
  return(i) ;
#endif
}



int    (Element_FindEquationPositionIndex)(Element_t const* element,const char* s)
/** Find the equation position index whose name is pointed to by s */
{
  int n = Element_GetNbOfEquations(element) ;
  const char* const* ss = (const char* const*) Element_GetNameOfEquation(element) ;
  int    i = String_FindPositionIndex(s,ss,n) ;
  
  return(i) ;
#if 0
  if(isdigit(s[0])) { /* donne sous forme numerique */
    i  = atoi(s) - 1 ;
  } else {            /* donne sous forme alphabetique */
    for(i = 0 ; i < n ; i++) {
      if(!strncmp(s,ss[i],strlen(s))) break ;
    }
    if(i == n) i = -1 ;
  }

  //if(i < 0) arret("Element_FindEquationPositionIndex: unknown position of %s",s) ;
  return(i) ;
#endif
}



#if 1
double*  (Element_ComputeNormalVector)(Element_t const* element,double* dh,int nn,const int dim_h)
/** Compute the unit outward normal vector to a submanifold element
 *  of dimension dim-1 */
{
  size_t SizeNeeded = 3*sizeof(double) ;
  double* norm = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  //int    dim_h = Element_GetDimension(element) ;
  int    dim   = Element_GetDimensionOfSpace(element) ;

  if(dim_h != dim - 1) {
    arret("Element_ComputeNormalVector") ;
  }
  
  {
    double* jac = Element_ComputeJacobianMatrix(element,dh,nn,dim_h) ;

    /* The outward normal to the surface is opposed to the vectors of
     * the jacobian matrix (the vectors G2 or G3) because the material
     * is assumed to be located on the side of the surface pointed to
     * by G2 or G3 (see below Element_ComputeJacobianMatrix), i.e.:
     *   for dim_h = 1 the outward normal is opposed to G2 
     *   for dim_h = 2 the outward normal is opposed to G3 
     * 
     * ATTENTION: in the deprecated implementation the normal was G3
     * in case of 3D (dim_h = 2), see below.
     */
    {
      int j ;

      for(j = 0 ; j < 3 ; j++) {
        norm[j] = - jac[(j)*3 + (dim_h)] ;
      }
    }

    Element_FreeBufferFrom(element,jac) ;
  }

  return(norm) ;
}
#endif



#if 0
double*  (Element_ComputeNormalVector)(Element_t const* element,double* dh,int nn,const int dim_h)
/** Compute the unit outward normal vector to a submanifold element
 *  of dimension dim-1 */
{
#define DH(n,i) (dh[(n)*3+(i)])
//#define DH(n,i) (dh[(n)*dim_h+(i)])
  size_t SizeNeeded = 3*sizeof(double) ;
  double* norm = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  //int    dim_h = Element_GetDimension(element) ;
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
#endif



double* (Element_ComputeJacobianMatrix)(Element_t const* el,double* dh,int nn,const int dim_h)
/** Compute the 3x3 jacobian matrix */
{
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define JC(i,j)  (jc[(i)*3 + (j)])
  //int  dim_h = Element_GetDimension(el) ;
  int  dim   = Element_GetDimensionOfSpace(el) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* jc = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  
  #if 0
  ShapeFct_t* shapefct = Element_GetShapeFct(el) ;
  int  dim_h = ShapeFct_GetDimension(shapefct) ;
  int nn = ShapeFct_GetNbOfNodes(shapefct) ;
  double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
      
  ShapeFct_ComputeValuesAtPoint(dim_h,nn,a,NULL,dh) ;
  #endif
  

  if(dim_h > dim) {
    arret("Element_ComputeJacobianMatrix: impossible") ;
  }
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < 9 ; i++) jc[i] = 0. ;
  }
  
  
  /* The reduced jacobian matrix, JC */
  {
    int i ;
    
    for(i = 0 ; i < dim ; i++) {
      int    j ;
    
      for(j = 0 ; j < dim_h ; j++) {
        int k ;
        
        for(k = 0 ; k < nn ; k++) {
          double* x = Element_GetNodeCoordinate(el,k) ;
      
          JC(i,j) += x[i]*DH(k,j) ;
        }
      }
      /* The values JC(i,j) for (0 < i < dim , dim_h <= j < dim) are supplemented below. */
    }
    
    /* Supplement with 1 the diagonal terms i >= dim. */
    //for(i = dim_h ; i < 3 ; i++) {
    for(i = dim ; i < 3 ; i++) {
      JC(i,i) = 1 ;
    }
  }
  
  
  if(dim_h == dim) return(jc) ;


  /* 
   * In case of submanifold i.e. dim_h < dim, 
   * we complement the jacobian matrix:
   * 
   * 1. For a point in 1D, 2D and 3D: [J] = (e1,e2,e3)
   *    where (e1,e2,e3) is an orthonormed vector system.
   * 
   * 2. For a line in 2D or 3D: [J] = (J1,G2,G3)
   *    where (G2,G3) are orthonormed vectors in the normal plane to J1.
   *    That means that there is no deformation for any material element
   *    that is perpendicular to X axis in the reference frame.
   * 
   * 3. For a surface in 3D: [J] = (J1,J2,G3) 
   *    where G3 is the unit normal vector to (J1,J2).
   *    That means that here is no deformation for material element
   *    that is perpendicular to the XY plane in the reference frame.
   */
   
  

  /* In case of a point: [J] = (e1,e2,e3) */
  if(dim_h == 0) {
    int i ;
    
    for(i = 0 ; i < 3 ; i++) {
      JC(i,i) = 1 ;
    }
    
    return(jc) ;
  }


  #define J1(i)   JC((i),0)
  #define J2(i)   JC((i),1)
  #define J3(i)   JC((i),2)

  /* Compute G2 in case of line: [J] = (J1,G2,G3)  */
  if(dim_h == 1) {
    if(dim >= 2) {
      double c = sqrt(J1(0)*J1(0) + J1(1)*J1(1)) ;
      double b = sqrt(J1(0)*J1(0) + J1(2)*J1(2)) ;
      double a = sqrt(J1(1)*J1(1) + J1(2)*J1(2)) ;
      
      #define G2   J2
      /* G2 = (0,0,1) x J1 / |J1| */
      if(c > 0) { /* This is the case for dim = 2  (J1(2) = 0) */
        G2(0) = - J1(1)/c ;
        G2(1) =   J1(0)/c ;
        G2(2) =   0 ;
      /* G2 = (0,1,0) x J1 / |J1| */
      } else if(b > 0) {
        G2(0) =   J1(2)/b ;
        G2(1) =   0 ;
        G2(2) = - J1(0)/b ;
      /* G2 = (1,0,0) x J1 / |J1| */
      } else if(a > 0) {
        G2(0) =   0 ;
        G2(1) = - J1(2)/a ;
        G2(2) =   J1(1)/a ;
      /* J1 = (0,0,0) -> G2 = (0,1,0) */
      } else {
        G2(0) =   0 ;
        G2(1) =   1 ;
        G2(2) =   0 ;
        //arret("Element_ComputeJacobianMatrix: impossible") ;
      }
      #undef G2
    }
    
    if(dim == 2) { /* In this case G3 = (0,0,1) */
      return(jc) ;
    }
  }
    
    
    
  /* Compute G3 in case of a line or a surface: [J] = (J1,J2,G3) */
  if(dim_h >= 1) {
    if(dim == 3) {
      double v[3] ;
      
      v[0] = J1(1)*J2(2) - J1(2)*J2(1) ;
      v[1] = J1(2)*J2(0) - J1(0)*J2(2) ;
      v[2] = J1(0)*J2(1) - J1(1)*J2(0) ;
      
      {
        double a = 0 ;
        int i ;
        
        for(i = 0 ; i < 3 ; i++) {
          a += v[i]*v[i] ;
        }
      
        a = sqrt(a) ;
      
        #define G3   J3
        /* G3 = J1 x J2 / |J1 x J2| */
        if(a > 0) {
          for(i = 0 ; i < 3 ; i++) {
            G3(i) = v[i]/a ;
          }
        } else {
          double m1 = sqrt(J1(0)*J1(0) + J1(1)*J1(1) + J1(2)*J1(2)) ;
          double m2 = sqrt(J2(0)*J2(0) + J2(1)*J2(1) + J2(2)*J2(2)) ;
          
          arret("Element_ComputeJacobianMatrix: surface with zero width?") ;
          
          /* G3 =  */
          if(m1 > 0) {
            G3(0) =   0 ;
            G3(1) =   0 ;
            G3(2) =   1 ;
          }
        }
        #undef G3
      }
    }
  }
  
  
  #undef J3
  #undef J2
  #undef J1
  
  
  return(jc) ;
#undef DH
#undef JC
}



#if 1
double (Element_ComputeJacobianDeterminant)(Element_t const* el,double* dh,int nn,const int dim_h)
/** Compute the determinant of the jacobian matrix */
{
  double* jac  = Element_ComputeJacobianMatrix(el,dh,nn,dim_h) ;
  double det = Math_Compute3x3MatrixDeterminant(jac) ;
  
  Element_FreeBufferFrom(el,jac) ;
  
  if(det < 0) {
    int index = Element_GetElementIndex(el) ;
    
    arret("Element_ComputeJacobianDeterminant: negative determinant\n\
    (det = %e) at element %d",det,index) ;
  }
  
  return(det) ;
}
#endif


#if 1
double* (Element_ComputeInverseJacobianMatrix)(Element_t const* el,double* dh,int nn,const int dim_h)
/** Compute the inverse jacobian matrix */
{
  size_t SizeNeeded = 9*sizeof(double) ;
  double* cj = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  
  {
    double* jc = Element_ComputeJacobianMatrix(el,dh,nn,dim_h) ;
    double* b  = Math_Inverse3x3Matrix(jc) ;
    int i ;
    
    if(b) {
      for(i = 0 ; i < 9 ; i++) {
        cj[i] = b[i] ;
      }
    } else {
      arret("Element_ComputeInverseJacobianMatrix: not invertible") ;
    }
  }
    
  return(cj) ;
}
#endif



double* (Element_ComputeCoordinateVector)(Element_t const* el,double* h)
/** Compute the coordinate vector */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* coor = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;

  {
    int i ;
    
    for(i = 0 ; i < 3 ; i++) coor[i] = 0 ;
    
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(el,i) ;
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        coor[j] += h[i] * x[j] ;
      }
    }
    
  }

  return(coor) ;
}



#if 1
double*  (Element_ComputeCoordinateInReferenceFrame)(Element_t const* el,double* x)
/** Compute the local coordinates in the reference element frame 
 *  which map into coordinates "x" in the actual space */
{
  unsigned int dim   = Element_GetDimensionOfSpace(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double* x_e[Element_MaxNbOfNodes] ;
  double diameter ;
  
  
  /* Node coordinates */
  {
    int    i ;
    
    for(i = 0 ; i < nn ; i++) {
      x_e[i] = Element_GetNodeCoordinate(el,i) ;
      //Math_PrintVector(x_e[i],dim) ;
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
    int dim_h = ShapeFct_GetDimension(shapefct) ;
    int nf    = ShapeFct_GetNbOfNodes(shapefct) ;
    double* a = ShapeFct_GetCoordinate(shapefct) ;
    int    max_iter = 20 ;
    double tol = 1.e-6 ;
    int    iter ;
    double err ;
    int i ;
    
    for(i = 0 ; i < 3 ; i++) {
      a[i] = 0 ;
    }
    
    if(nn == 1) return(a) ;
    if(dim_h == 0) return(a) ;
    
    for(iter = 0 ; iter < max_iter ; iter++) {
      double* h = ShapeFct_GetFunction(shapefct) ;
      double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
      double r[3] ;
      double k[9] ;
      
      ShapeFct_ComputeValuesAtPoint(dim_h,nf,a,h,dh) ;
      
      //Math_PrintVector(a,3) ;
      
      {
        double* jac = Element_ComputeJacobianMatrix(el,dh,nf,dim_h) ;
        
        //Math_PrintMatrix(jac,3) ;
    
        for(i = 0 ; i < dim ; i++) {
          int   j,in ;
      
          r[i] = x[i] ;
      
          //for(in = 0 ; in < nn ; in++) {
          for(in = 0 ; in < nf ; in++) {
            r[i] -= h[in]*x_e[in][i] ;
          }
          
          for(j = 0 ; j < dim ; j++) {
            k[dim*i+j] = jac[i*3+j] ;
          }
        }
      }
    
      
      if(diameter > 0) {
       for(i = 0 ; i < dim ; i++) {
          int  j ;
      
          r[i] /= diameter ;
        
          for(j = 0 ; j < dim ; j++) k[dim*i+j] /= diameter ;
        }
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
}
#endif



int (Element_ComputeNbOfSolutions)(Element_t const* el)
{
  Solution_t* solution = Element_GetSolution(el) ;
  int n = 0 ;
      
  if(solution) {
    do {
      n++ ;
      solution = Solution_GetPreviousSolution(solution) ;
    } while(solution != Element_GetSolution(el)) ;
  }
  
  return(n) ;
}



int* (Element_ComputeMatrixRowAndColumnIndices)(Element_t const* el)
{
  int  nn  = Element_GetNbOfNodes(el) ;
  int  neq = Element_GetNbOfEquations(el) ;
  
  size_t SizeNeeded = 2*nn*neq*sizeof(int) ;
  int* rowind = (int*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int* colind = rowind + nn*neq ;
  int  i ;
    
  for(i = 0 ; i < nn ; i++) {
    Node_t* node_i = Element_GetNode(el,i) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int ij = i*neq + j ;
      int ii_col = Element_GetUnknownPosition(el)[ij] ;
      int ii_row = Element_GetEquationPosition(el)[ij] ;
      
      colind[ij] = (ii_col >= 0) ? Node_GetMatrixColumnIndex(node_i)[ii_col] : -1 ;
      rowind[ij] = (ii_row >= 0) ? Node_GetMatrixRowIndex(node_i)[ii_row] : -1 ;
    }
  }
  
  return(rowind) ;
}



int* (Element_ComputeSelectedMatrixRowAndColumnIndices)(Element_t const* el,const int imatrix)
{
  int nn   = Element_GetNbOfNodes(el) ;
  int neq  = Element_GetNbOfEquations(el) ;
  int ndof = nn*neq ;
  
  size_t SizeNeeded = 2*nn*neq*sizeof(int) ;
  int* rowind = (int*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int* colind = rowind + ndof ;
  int  i ;
    
  for(i = 0 ; i < nn ; i++) {
    Node_t* node_i = Element_GetNode(el,i) ;
    int    j ;
    
    for(j = 0 ; j < neq ; j++) {
      int ij = i*neq + j ;
      int ii_col = Element_GetUnknownPosition(el)[ij] ;
      int ii_row = Element_GetEquationPosition(el)[ij] ;
      
      colind[ij] = Node_GetSelectedMatrixColumnIndexOf(node_i,ii_col,imatrix) ;
      rowind[ij] = Node_GetSelectedMatrixRowIndexOf(node_i,ii_row,imatrix) ;
    }
  }
  
  return(rowind) ;
}



double (Element_ComputeSize)(Element_t const* element)
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



double* (Element_ComputeSizes)(Element_t const* element)
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* size = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;

  {
    ShapeFct_t* shapefct = Element_GetShapeFct(element) ;
    double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
    int dim_e = Element_GetDimension(element) ;
  
    {
      double  a[3] = {0,0,0} ;
      
      ShapeFct_ComputeValuesAtPoint(dim_e,nn,a,NULL,dh) ;
    }
    
    {
      double* jac = Element_ComputeJacobianMatrix(element,dh,nn,dim_e) ;
      int i ;
      
      /* Compute the length of each segment of the element
       * i.e. in the direction of the local frame. */
      #define J(i,j)  (jac[(i)*3 + (j)])
      for(i = 0 ; i < dim ; i++) {
        double d = 0 ;
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          d += J(j,i) * J(j,i) ;
        }
        
        size[i] = sqrt(d) ;
      }
      #undef J
      
      for(i = dim ; i < 3 ; i++) {
        size[i] = 0 ;
      }
      
      /* Order from largest to smallest */
      for(i = 0 ; i < 3 ; i++) {
        int j ;
        
        for(j = i + 1 ; j < 3 ; j++) {
          if(size[j] > size[i]) Math_Swap(size[i],size[j],double) ;
        }
      }
      
    }
  }
  
  return(size) ;
}




int (Element_HasZeroThickness)(Element_t const* element)
{
  int nn = Element_GetNbOfNodes(element) ;
  ShapeFct_t* shapefct = Element_GetShapeFct(element) ;
  int dim_e = Element_GetDimension(element) ;
  
  if(!shapefct) return(0) ;
  
  {
    double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
    double  a[3] = {0,0,0} ;
      
    ShapeFct_ComputeValuesAtPoint(dim_e,nn,a,NULL,dh) ;
  }
    
  {
    double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
    double* jac  = Element_ComputeJacobianMatrix(element,dh,nn,dim_e) ;
    double det = Math_Compute3x3MatrixDeterminant(jac) ;
    double h = Element_ComputeSize(element) ;
    double h3 = h*h*h ;
    double t = fabs(det) ;
    
    Element_FreeBufferFrom(element,jac) ;
    
    if(t < 1.e-10*h3) {
      return(Element_NbOfOverlappingNodes(element)) ;
    } else {
      return(0) ;
    }
  }
  
  return(0) ;
}



int (Element_NbOfOverlappingNodes)(Element_t const* element)
/** Return the number of overlapping nodes for zero-thickness element. */
{
  int nn  = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimension(element) ;
  
  if(dim == 0) {
    switch (nn) {
      default: break;
    }
  } else if(dim == 1) {
    switch (nn) {
      case 2 : return 1;              /* line 1 */
      default: break;
    }
  } else if(dim == 2) {
    switch (nn) {
      case 3 : return 1;              /* triangle 1 */
      case 4 : return 2;              /* quadrangle 1 */
      default: break;
    }
  } else if(dim == 3) {
    switch (nn) {
      case 4 : return 1;              /* tetrahedron 1 */
      case 8 : return 4;              /* hexahedron 1 */
      case 6 : return 3;              /* prism 1 */
      case 5 : return 1;              /* pyramid 1 */
      default: break;
    }
  }
  
  arret("Element_NbOfOverlappingNodes: not available") ;
  
  return(-1) ;
}




int (Element_OverlappingNode)(Element_t const* element,const int n)
/** Return the index of the first node other than n 
 *  overlapping node n or n if it fails. */
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  double* xn = Element_GetNodeCoordinate(element,n) ;
  double h = Element_ComputeSize(element) ;
  double tiny = h*1.e-10 ;

  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      double* xi = Element_GetNodeCoordinate(element,i) ;
      
      if(i == n) continue ;
      
      {
        double r = 0 ;
        int j ;
    
        for(j = 0 ; j < dim ; j++) {
          double y = xn[j] - xi[j] ;
      
          r += y*y ;
        }
    
        r = sqrt(r) ;
    
        if(r < tiny) return(i) ;
      }
    }
  }
  
  return(n) ;
}




#if 0
double (Element_ComputeThickness)(Element_t const* element)
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  double thickness = 0 ;
  
  /*  */
  {
    int nn1 = nn/2 ;
    int i ;
    
    for(i = 0 ; i < nn1 ; i++) {
      double* x1 = Element_GetNodeCoordinate(element,i) ;
      double* x2 = Element_GetNodeCoordinate(element,nn1 + i) ;
      double r = 0 ;
      int j ;
    
      for(j = 0 ; j < dim ; j++) {
        double y = x1[j] - x2[j] ;
      
        r += y*y ;
      }
    
      r = sqrt(r) ;
    
      if(r > thickness) thickness = r ;
    }
  }
  
  {
    int dim_e = Element_GetDimension(el) ;
    ShapeFct_t* shapefct = Element_GetShapeFct(el) ;
    //double* a = ShapeFct_GetCoordinate(shapefct) ;
    double* h = ShapeFct_GetFunction(shapefct) ;
    double* dh = ShapeFct_GetFunctionGradient(shapefct) ;
    double  a[3] = {0,0,0} ;
      
    ShapeFct_ComputeValuesAtPoint(dim_e,nn,a,h,dh) ;
  }
  
  return(thickness) ;
}
#endif



#if 0
double (Element_DistanceBetweenNodes)(Element_t const* element,const int i,const int j)
{
  int nn = Element_GetNbOfNodes(element) ;
  int dim = Element_GetDimensionOfSpace(element) ;
  double* xi = Element_GetNodeCoordinate(element,i) ;
  double* xj = Element_GetNodeCoordinate(element,j) ;
  double len = 0 ;
  
  {
    int k ;
    
    for(k = 0 ; k < dim ; k++) {
      double y = xj[k] - xi[k] ;
      
      len += y*y ;
    }
    
    len = sqrt(len) ;
  }
  
  return(len) ;
}
#endif




int (Element_ComputeNbOfMatrixEntries)(Element_t const* element)
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
  
  Element_FreeBufferFrom(element,row) ;
  
  return(len) ;
}




int (Element_ComputeNbOfSelectedMatrixEntries)(Element_t const* element,const int imatrix)
{
  int   ndof = Element_GetNbOfDOF(element) ;
  int*  row  = Element_ComputeSelectedMatrixRowAndColumnIndices(element,imatrix) ;
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
  
  Element_FreeBufferFrom(element,row) ;
  
  return(len) ;
}




void (Element_MakeUnknownContinuousAcrossZeroThicknessElement)(Element_t const* element,const char* name)
{
  ShapeFct_t* shapefct = Element_GetShapeFct(element) ;
  int nf = ShapeFct_GetNbOfNodes(shapefct) ;
  int in ;
  
  for(in = 0 ; in < nf ; in++) {
    Node_t* node = Element_GetNode(element,in) ;
          
    Node_MakeUnknownContinuousAtOverlappingNodes(node,name) ;
  }
}



void (Element_MakeEquationContinuousAcrossZeroThicknessElement)(Element_t const* element,const char* name)
{
  ShapeFct_t* shapefct = Element_GetShapeFct(element) ;
  int nf = ShapeFct_GetNbOfNodes(shapefct) ;
  int in ;
  
  for(in = 0 ; in < nf ; in++) {
    Node_t* node = Element_GetNode(element,in) ;
          
    Node_MakeEquationContinuousAtOverlappingNodes(node,name) ;
  }
}





int (Element_FindNodeIndex)(Element_t const* element,const Node_t* node)
/** Return the local node index matching that of node
 *  or -1 if it fails. */
{
  int nn  = Element_GetNbOfNodes(element) ;
  int in ;
          
  for(in = 0 ; in < nn ; in++) {
    if(Node_GetNodeIndex(node) == Element_GetNodeIndex(element,in)) {
      return(in) ;
    }
  }
          
  return(-1) ;
}



/* The two following functions have not been tested yet */
#if 0
double* (Element_ComputeDiscreteGradientOperator)(Element_t const* element,IntFct_t* intfct)
/** Compute the discrete gradient operator.
 *  On ouput:
 *  a pointer to double with allocated space of DIM*NN*NN
 *  DIM = dimension of space
 *  NN  = nb of nodes of the element
 *  The array dg contains the contribution of the element
 *  to the discrete gradient operator.
 *  Usage of dg: 
 *  the gradient of u at node i is calculated as follows
 *  gi(k) = 1/m_i sum_elts sum_j DG(i*dim+k,j) u_j  k=1,dim j=1,NN
 *  where DG(i,j) is defined as (dg[(i)*nn + (j)]) and
 *  m_i is lumped-mass at node i.
 *  This should be done on the set of elements containing the node i.
 *  Ref.
 *  D. Kuzmin, A guide to numerical methods for transport equations, 2010) */
{
  int dim = Element_GetDimensionOfSpace(element) ;
  int nn  = Element_GetNbOfNodes(element) ;
  int nf = IntFct_GetNbOfFunctions(intfct) ;
  size_t SizeNeeded = dim*nn*nn*sizeof(double) ;
  double* dg = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < dim*nn*nn ; i++) dg[i] = 0 ;
  }
  
  if(Element_IsSubmanifold(element)) {
    arret("Element_ComputeDiscreteGradientOperator") ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(element,i) ;
    }
  }
  
  
  if(nn > nf) {
    arret("Element_ComputeDiscreteGradientOperator") ;
  }
  

  {
    Symmetry_t sym = Element_GetSymmetry(element) ;
    int dim_h  = IntFct_GetDimension(intfct) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    double* weight = IntFct_GetWeight(intfct) ;
    double m[Element_MaxNbOfNodes] ;
    int p ;
    
    {
      int i ;
    
      for(i = 0 ; i < nf ; i++) m[i] = 0 ;
    }
    
    for(p = 0 ; p < np ; p++) {
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = Element_ComputeJacobianDeterminant(element,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = Element_ComputeInverseJacobianMatrix(element,dh,nf,dim_h) ;
    
      /* axisymetrical or spherical cases */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double radius = 0 ;
        int i ;
        
        for(i = 0 ; i < nf ; i++) radius += h[i]*x[i][0] ;
        a *= 2*M_PI*radius ;
        if(Symmetry_IsSpherical(sym)) a *= 2*radius ;
      }
    
  
      #define DG(i,j)  (dg[(i)*nn + (j)])
      #define CAJ(i,j) (caj[(i)*3 + (j)])
      #define DH(n,i)  (dh[(n)*3 + (i)])
      /* DG(i*dim+k,j) = int_V H(i)*DH(j,l)*CAJ(l,k)*J*dV */
      {
        int i ;
        
        for(i = 0 ; i < nf ; i++) {
          int k ;
          
          /* m is the surface/volume of the vertex-centered cell */
          m[i] += a*h[i] ;
        
          /* the discrete gradient operator: dg */
          for(k = 0 ; k < dim ; k++) {
            int j ;
          
            for(j = 0 ; j < nf ; j++) {
              int l ;
            
              for(l = 0 ; l < dim_h ; l++) {
                DG(i*dim+k,j) += a*h[i]*DH(j,l)*CAJ(l,k) ;
              }
            }
          }
        }
      }
      #undef CAJ
      #undef DH
      #undef DG
    }
  }
  
  return(dg) ;
}




double* (Element_ComputeLumpedMass)(Element_t const* element,IntFct_t* intfct)
/** Compute the lumped mass.
 *  On ouput:
 *  a pointer to double with allocated space of NN
 *  NN  = nb of nodes of the element
 *  Ref.
 *  D. Kuzmin, A guide to numerical methods for transport equations, 2010) */
{
  int nn  = Element_GetNbOfNodes(element) ;
  int nf = IntFct_GetNbOfFunctions(intfct) ;
  size_t SizeNeeded = nn*sizeof(double) ;
  double* lum = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) lum[i] = 0 ;
  }
  
  if(Element_IsSubmanifold(element)) {
    arret("Element_ComputeLumpedMass") ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(element,i) ;
    }
  }
  
  
  if(nn > nf) {
    arret("Element_ComputeLumpedMass") ;
  }
  

  {
    Symmetry_t sym = Element_GetSymmetry(element) ;
    int dim_h  = IntFct_GetDimension(intfct) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    double* weight = IntFct_GetWeight(intfct) ;
    int p ;
    
    for(p = 0 ; p < np ; p++) {
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = Element_ComputeJacobianDeterminant(element,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
    
      /* axisymetrical or spherical cases */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double radius = 0 ;
        int i ;
        
        for(i = 0 ; i < nf ; i++) radius += h[i]*x[i][0] ;
        a *= 2*M_PI*radius ;
        if(Symmetry_IsSpherical(sym)) a *= 2*radius ;
      }
    
  
      /* lum[i] = int_V H(i)*J*dV */
      {
        int i ;
        
        for(i = 0 ; i < nf ; i++) {
          lum[i] += a*h[i] ;
        }
      }
    }
  }
  
  return(lum) ;
}
#endif




#if 0
ElementSol_t* (Element_GetDeepElementSol)(Element_t const* el,unsigned int depth)
{
  Solution_t* solution = Element_GetSolution(el) ;
  ElementSol_t* elementsol ;
  
  while(depth--) {
    solution = Solution_GetPreviousSolution(solution) ;
  }
  
  elementsol = Solution_GetElementSol(solution) ;
  
  return(elementsol) ;
}
#endif

  

void  (Element_CopyCurrentSolutionIntoPreviousSolution)(Element_t const* el)
/** Copy the current solution into the previous solution */
{
  {
    ElementSol_t* elementsol = Element_GetElementSol(el) ;
    ElementSol_t* elementsol_n = Element_GetPreviousElementSol(el) ;
    
    ElementSol_Copy(elementsol_n,elementsol) ;
  }
  {
    int nn = Element_GetNbOfNodes(el) ;
    int i ;

    for(i = 0 ; i < nn ; i++) {
      Node_t* node = Element_GetNode(el,i) ;
      NodeSol_t* nodesol = Node_GetNodeSol(node) ;
      NodeSol_t* nodesol_n = Node_GetPreviousNodeSol(node) ;
      
      NodeSol_Copy(nodesol_n,nodesol) ;
    }
  }
}


double** (Element_ConvertToNodalUnknown)(Element_t const* el,double* f,int shift)
/** Convert the node quantity pointed to by "f" into a simulated nodal unknown "u"
 *  at position 0. */
{
  double uu[Element_MaxNbOfNodes][Node_MaxNbOfEquations];
  double* u[Element_MaxNbOfNodes];
  int nn = Element_GetNbOfNodes(el) ;
  
  {
    for(int i = 0; i < nn; i++) {
      u[i] = uu[i];
      uu[i][Element_GetNodalUnknownPosition(el,i,0)] = f[i*shift] ;
    }
  }
    
  return(u);
}



double   (Element_IntegrateOverElement)(Element_t* el,IntFct_t* intfct,double* f,int shift)
/** Integrate f over the element through the integration points. 
 *  Return the integration result. */
{
  int nn = IntFct_GetNbOfFunctions(intfct) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double* x[Element_MaxNbOfNodes] ;
  double sum = 0 ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }

  /* 0D */
  if(dim_h == 0) {
    if(nn == 1) {
      double radius = x[0][0] ;
      
      sum = f[0] ;
      
      if(Symmetry_IsCylindrical(sym)) sum *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) sum *= 4*M_PI*radius*radius ;
      
    } else {
      arret("FEM_IntegrateOverElement: impossible") ;
    }
    
    return(sum) ;
  }

  /* 1D, 2D, 3D */
  {
    int np = IntFct_GetNbOfPoints(intfct) ;
    double* weight = IntFct_GetWeight(intfct) ;
    int p ;
  
    for(p = 0 ; p < np ; p++) {
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = Element_ComputeJacobianDeterminant(el,dh,nn,dim_h) ;
      double a   = weight[p]*d ;
    
      /* Axisymmetrical or spherical case */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
        double radius = 0 ;
      
        for(i = 0 ; i < nn ; i++) radius += h[i]*x[i][0] ;
        
        a *= 2*M_PI*radius ;
        if(Symmetry_IsSpherical(sym)) a *= 2*radius ;
      }
    
      sum += a*f[p*shift] ;
    }
  }
  
  return(sum) ;
}



double*  (Element_ComputeIsoShapeFctInActualSpace)(Element_t* el,double* x)
/* fonction d'interpolation (h) et ses derivees (dh) en x */
{
  double* a = Element_ComputeCoordinateInReferenceFrame(el,x) ;
  
  {
    unsigned short int dim_e = Element_GetDimension(el) ;
    int nn = Element_GetNbOfNodes(el) ;
    size_t SizeNeeded = (nn*4)*sizeof(double) ;
    double* h = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
    double* dh = h + Element_GetNbOfNodes(el) ;
  
    ShapeFct_ComputeValuesAtPoint(dim_e,nn,a,h,dh) ;

    return(h) ;
  }
}



void (Element_CheckNumberingOfOverlappingNodes)(Element_t const* el,const int nf)
{
  int j ;

  for(j = 0 ; j < nf ; j++) {
    int jj = Element_OverlappingNode(el,j) ;

    if(jj != j && jj < nf) {
      arret("FEM_CheckNumberingOfOverlappingNodes: bad numbering") ;
    }
  }
}


double (Element_ComputeUnknown)(Element_t* el,double const* const* u,IntFct_t* intfct,int p,int inc)
/** Compute the unknown located at "inc" at the interpolation point "p" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
  int nn = IntFct_GetNbOfFunctions(intfct) ;
  double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
  double par = 0. ;
  int i ;
        
  for(i = 0 ; i < nn ; i++) {
    par += h[i]*U(i) ;
  }
        
  return(par) ;
#undef U
}



double* (Element_ComputeDisplacementVector)(Element_t* el,double const* const* u,IntFct_t* intfct,int p,int inc)
/** Compute the displacement vector located at 
 *  "inc,inc+1,inc+2" at the interpolation point "p" */
{
  size_t SizeNeeded = 3*sizeof(double) ;
  double* dis = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;

  {
    int dim = Element_GetDimensionOfSpace(el) ;
    int i ;
    
    /* In case of a zero-thickness element (fracture) 
     * we compute the discontinuous displacement vector */
    #if 0
    if(Element_HasZeroThickness(el)) {
      int nf = IntFct_GetNbOfFunctions(intfct) ;
      double* du[Element_MaxNbOfNodes] ;
      double  du1[Element_MaxNbOfNodes*Element_MaxNbOfDOF] ;
      
      for(i = 0 ; i < nf ; i++) {
        int ii = Element_OverlappingNode(el,i) ;
        int j ;
        
        du[i] = du1 + i*Element_MaxNbOfDOF ;
        
        for(j = 0 ; j < dim ; j++) {
          du[i][inc + j] = u[i][inc + j] - u[ii][inc + j] ;
        }
      }
    
      for(i = 0 ; i < dim ; i++) {
        dis[i] = Element_ComputeUnknown(el,du,intfct,p,inc+i) ;
      }
    } else 
    #endif
    {
      for(i = 0 ; i < dim ; i++) {
        dis[i] = Element_ComputeUnknown(el,u,intfct,p,inc+i) ;
      }
    }
    
    for(i = dim ; i < 3 ; i++) dis[i] = 0 ;
  }
        
  return(dis) ;
}



double* (Element_ComputeUnknownGradient)(Element_t* el,double const* const* u,IntFct_t* intfct,int p,int inc)
/** Compute the unknown gradient located at "inc" at the interpolation point "p" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  int dim = Element_GetDimensionOfSpace(el) ;
  int dim_e = Element_GetDimension(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int nn = IntFct_GetNbOfFunctions(intfct) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  
  
  /* One node mesh: for this special element the unknown stands for the gradient along x-axis */
  if(nn == 1 && dim_e == 3){
    int i ;
    
    for(i = 0 ; i < 3 ; i++) {
      grad[i] = 0 ;
    }

    grad[0] = U(0) ;
    
    return(grad) ;
  }
  

  {
    /* interpolation functions */
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    double* gu = grad ;
    double grf[3] = {0,0,0} ;
  
    /* the gradient in the reference frame */
    {
      int l ;
        
      for(l = 0 ; l < dim_h ; l++) {
        int k ;
      
        for(k = 0 ; k < nn ; k++) {
          grf[l] += U(k) * DH(k,l) ;
        }
      }
    }
  
    /* the parameter gradient (gu) */
    {
      /* inverse jacobian matrix (cj) */
      double* cj = Element_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
      for(i = 0 ; i < dim ; i++) {
        int l ;
        
        for(l = 0 ; l < dim_h ; l++) {
          gu[i] += grf[l] * CJ(l,i) ;
        }
      }
      
      Element_FreeBufferFrom(el,cj) ;
    }
  }
  
  return(grad) ;
  
#undef U
#undef DH
#undef CJ
}



double* (Element_ComputeLinearStrainTensor)(Element_t* el,double const* const* u,IntFct_t* intfct,int p,int inc)
/** Compute the 3D linearized strain tensor for the displacement vector 
 *  located at "inc" and at the interpolation point "p" */
{
#define U(n,i)   (u[n][Element_GetNodalUnknownPosition(el,n,inc + (i))])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define STRAIN(i,j) (strain[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int dim_e = Element_GetDimension(el) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* strain = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int nn = IntFct_GetNbOfFunctions(intfct) ;
  
  
  /* One node mesh: for this special element the displacement u_i stands for strain_ii */
  if(nn == 1 && dim_e == 3){
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      strain[i] = 0 ;
    }

    for(i = 0 ; i < dim ; i++) {
      STRAIN(i,i) = U(0,i) ;
    }
    
    return(strain) ;
  }


  {
    /* interpolation functions */
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    /* inverse jacobian matrix (cj) */
    double* cj = Element_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
    double  gu[3][3] ;
    int     i ;
  
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++) {
      int j ;
        
      for(j = 0 ; j < 3 ; j++)  {
        gu[i][j] = 0. ;
      }
    }
  
    /* displacement gradient (gu) */
    if(Element_HasZeroThickness(el)) {
      double* norm = Element_ComputeNormalVector(el,dh,nn,dim_h) ;
      double du[3] = {0,0,0} ;
      
      /* the displacement discontinuity */
      for(i = 0 ; i < dim ; i++) {
        int k ;
          
        for(k = 0 ; k < nn ; k++) {
          int kk = Element_OverlappingNode(el,k) ;
            
          du[i] += (U(k,i) - U(kk,i)) * h[k] ;
        }
      }
      
      #define NORM(i)  (norm[i])
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          gu[i][j] = du[i] * NORM(j) ;
        }
      }
      #undef NORM
      
    } else if(dim_h > 0) {
      double  grf[3][3] = {{0,0,0},{0,0,0},{0,0,0}} ;
      
      /* the gradient in the reference frame */
      for(i = 0 ; i < dim ; i++) {
        int l ;
          
        for(l = 0 ; l < dim_h ; l++) {
          int k ;
          
          for(k = 0 ; k < nn ; k++) {
            grf[i][l] += U(k,i) * DH(k,l) ;
          }
        }
      }
      
      /* the gradient in the current frame */
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          int l ;
          
          for(l = 0 ; l < dim_h ; l++) {
            gu[i][j] += grf[i][l] * CJ(l,j) ;
          }
        }
      }
    } else {
      Message_FatalError("FEM_ComputeLinearStrainTensor: dim_h <= 0!") ;
    }
      
    {
      //double* eps = strain ;
  
      /* Linearized strain tensor */
      for(i = 0 ; i < 3 ; i++) {
        int j ;
        
        for(j = 0 ; j < 3 ; j++) {
          STRAIN(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
        }
      }
  
      /* symmetric cases: axisymmetrical or spherical */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        
        if(Element_HasZeroThickness(el)) {
          /* No additif terms */
        } else {
          double radius = 0. ;
          double u_r   = 0. ;
          
          for(i = 0 ; i < nn ; i++) {
            double* x = Element_GetNodeCoordinate(el,i) ;
          
            radius += h[i]*x[0] ;
            u_r    += h[i]*U(i,0) ;
          }
          
          if(radius > 0) {
            STRAIN(2,2) += u_r/radius ;
          
            if(Symmetry_IsSpherical(sym)) STRAIN(1,1) += u_r/radius ;
          }
        }
      }
    }
      
    Element_FreeBufferFrom(el,cj) ;
  }
  
  return(strain) ;
#undef U
#undef DH
#undef STRAIN
#undef CJ
}



double* (Element_ComputeInternodeDistances)(Element_t* el)
/** Compute the intercell distances */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* dist = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < nn*nn ; i++)  dist[i] = 0. ;

  for(i = 0 ; i < nn ; i++) {
    double* xi = Element_GetNodeCoordinate(el,i) ;
    int j ;
    for(j = i + 1 ; j < nn ; j++) {
      double* xj = Element_GetNodeCoordinate(el,j) ;
      double dij = 0 ;
      int n ;
      for (n = 0 ; n < dim ; n++) {
        double xij = xj[n] - xi[n] ;
        dij += xij*xij ;
      }
      dij = sqrt(dij) ;
      dist[nn*i + j] = dij ;
      dist[nn*j + i] = dij ;
    }
  }
  
  return(dist) ;
}
