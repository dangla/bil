#include "FVM.h"
#include "Message.h"
#include "Math_.h"
#include "Geometry.h"
#include "Elements.h"
#include "Nodes.h"
#include "Session.h"
#include "GenericData.h"
#include "Mry.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <assert.h>


//static FVM_t* instancefvm = NULL ;

static FVM_t*  (FVM_Create)(void) ;
static int     (FVM_FindHalfSpace)(FVM_t*,int,int,double*) ;


FVM_t* (FVM_Create)(void)
{
  FVM_t* fvm = (FVM_t*) Mry_New(FVM_t) ;
  
  
  /* Space allocation for output */
  {
    size_t sz = FVM_MaxSizeOfOutput ;
    double* output = (double*) Mry_New(sz) ;
    
    FVM_GetOutput(fvm) = output ;
  }
  
  
  /* Space allocation for input */
  {
    size_t sz = FVM_MaxSizeOfInput ;
    double* input = (double*) Mry_New(sz) ;
    
    FVM_GetInput(fvm)= input ;
  }
  
  
  /* Space allocation for buffer */
  {
    Buffers_t* buf = Buffers_Create(FVM_SizeOfBuffer) ;
    
    FVM_GetBuffers(fvm) = buf ;
  }
  
  return(fvm) ;
}


void (FVM_Delete)(void* self)
{
  FVM_t* fvm = (FVM_t*) self ;
  
  free(FVM_GetOutput(fvm)) ;
  free(FVM_GetInput(fvm)) ;
  
  {
    Buffers_t* buf = FVM_GetBuffers(fvm) ;
    
    if(buf) {
      Buffers_Delete(buf) ;
      free(buf) ;
      FVM_GetBuffers(fvm) = NULL ;
    }
  }
}



#if 0
//FVM_t* FVM_GetInstance1(Element_t* el)
{
  if(!instancefvm) {
    instancefvm = FVM_Create() ;
  }
  
  FVM_GetElement(instancefvm) = el ;
  
  FVM_FreeBuffer(instancefvm) ;
  FVM_GetCellVolumes(instancefvm) = NULL ;
  FVM_GetCellSurfaceAreas(instancefvm) = NULL ;
  FVM_GetIntercellDistances(instancefvm) = NULL ;
  
  return(instancefvm) ;
}
#endif



FVM_t*  (FVM_GetInstance)(Element_t* el)
{
  GenericData_t* gdat = Session_FindGenericData(FVM_t,"FVM") ;
  
  if(!gdat) {
    FVM_t* fvm = FVM_Create() ;
    
    gdat = GenericData_Create(1,fvm,FVM_t,"FVM") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(FVM_t,"FVM")) ;
  }
  
  {
    FVM_t* fvm = (FVM_t*) GenericData_GetData(gdat) ;
  
    FVM_GetElement(fvm) = el ;
    FVM_FreeBuffer(fvm) ;
    FVM_GetCellVolumes(fvm) = NULL ;
    FVM_GetCellSurfaceAreas(fvm) = NULL ;
    FVM_GetIntercellDistances(fvm) = NULL ;
  
    return(fvm) ;
  }
}






double* (FVM_ComputeSurfaceLoadResidu)(FVM_t* fvm,Load_t* load,double t,double dt)
/* Compute the residu force due to surface loads (r) */
{
  Element_t* el = FVM_GetElement(fvm) ;
  double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  Geometry_t* geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  Node_t* *no = Element_GetPointerToNode(el) ;
  Field_t* field = Load_GetField(load) ;
  char    *load_eqn = Load_GetNameOfEquation(load) ;
  char    *load_type = Load_GetType(load) ;
  Function_t* function = Load_GetFunction(load) ;
  int    nn  = Element_GetNbOfNodes(el) ;
  int    neq = Element_GetNbOfEquations(el) ;
  int    ndof = nn*neq ;
  int    ieq = Element_FindEquationPositionIndex(el,load_eqn) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  size_t SizeNeeded = ndof*(sizeof(double)) ;
  double* r = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  double zero = 0. ;
  int    i ;

  /* initialization */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(field == NULL) return(r) ;

  /* Position index of the equation */
  if(ieq < 0) arret("FVM_ComputeSurfaceLoadResidu (1) : unknown equation") ;


  /* flux or cumulative flux*/
  if(strncmp(load_type,"flux",4) == 0 || strncmp(load_type,"cumulflux",4) == 0) {
    double ft = 1 ;
    
    
    if(strncmp(load_type,"flux",4) == 0) {
      ft = 1 ;
      
      if(function) {
        ft = Function_ComputeValue(function,t) ;
      }
    }
    
    if(strncmp(load_type,"cumulflux",4) == 0) {
      ft = dt ;
    
      if(function) {
        ft = Function_ComputeValue(function,t) - Function_ComputeValue(function,t - dt) ;
      }
    }
    
    if(dim == 1 && nn == 1) {
      double* x = Node_GetCoordinate(no[0]) ;
      double radius = x[0] ;
      
      r[ieq] = ft*Field_ComputeValueAtPoint(field,x,dim) ;
      
      if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
      return(r) ;
    }
    
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      
      for(i = 0 ; i < nn ; i++) {
        double* y = Node_GetCoordinate(no[i]) ;
        rb[i] = ft*Field_ComputeValueAtPoint(field,y,dim)*volume[i] ;
      }
      
      for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
      return(r) ;
    }
  }


  /* linear dependent flux: F = A*U */
  if(strncmp(load_type,"linearflux",4) == 0) {
    double ft = 1 ;
      
    if(function) {
      ft = Function_ComputeValue(function,t) ;
    }
    
    if(dim == 1 && nn == 1) {
      double v = u[0][ieq] ;
      double* x = Node_GetCoordinate(no[0]) ;
      double radius = x[0] ;
      
      r[ieq] = ft*Field_ComputeValueAtPoint(field,x,dim)*v ;
      
      if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
      return(r) ;
    }
    
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      
      for(i = 0 ; i < nn ; i++) {
        double v = u[i][ieq] ;
        double* y = Node_GetCoordinate(no[i]) ;
        rb[i] = ft*Field_ComputeValueAtPoint(field,y,dim)*volume[i]*v ;
      }
      
      for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
      return(r) ;
    }
  }
  
  
  /* force */
  if(strncmp(load_type,"force",5) == 0) {
    double ft = 1. ;
    
    if(function != NULL) {
      ft = Function_ComputeValue(function,t) ;
    }
    
    if(dim == 1 && nn == 1) {
      double* x = Node_GetCoordinate(no[0]) ;
      double radius = x[0] ;
      
      r[ieq] = ft*Field_ComputeValueAtPoint(field,x,dim) ;
      
      if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
      
      return(r) ;
    }
    
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      
      for(i = 0 ; i < nn ; i++) {
        double* y = Node_GetCoordinate(no[i]) ;
        rb[i] = ft*Field_ComputeValueAtPoint(field,y,dim) ;
      }
      
      for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
      return(r) ;
    }
  }
  
  arret("FVM_ComputeSurfaceLoadResidu (2) : unknown load") ;
  return(NULL) ;
}





double* (FVM_ComputeBodyForceResidu)(FVM_t* fvm,double const* f,int const dec)
/** Return the body force residu (r) */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    /* f is the "force" density at node i */
    r[i] = volume[i]*f[i*dec] ;
  }
  
  return(r) ;
}



double*  (FVM_ComputeFluxResidu)(FVM_t* fvm,double const* w,int const dec)
/** Return the flux residu (r) */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* surf = FVM_ComputeCellSurfaceAreas(fvm) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < nn ; i++)  r[i] = 0. ;

  for(i = 0 ; i < nn ; i++) {
    int j ;
    for(j = 0 ; j < nn ; j++) {
      /* wij is the outflow from cell i to cell j
       * Note that wii is not used */
      if(j != i) {
        double wij = w[i*dec + j] ;
        double a = surf[i*nn + j] ;
        r[i] += a*wij ;
      }
    }
  }
  
  return(r) ;
}



double* (FVM_ComputeMassAndFluxResidu)(FVM_t* fvm,double const* c,int const dec)
/** Return the body force and flux residu */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  int    i ;
  
  double* r = FVM_ComputeFluxResidu(fvm,c,dec) ;
  
  for(i = 0 ; i < nn ; i++) {
    /* f is the "force" density at node i */
    double f = c[i*dec + i] ;
    r[i] += volume[i]*f ;
  }
  
  return(r) ;
}



double* (FVM_ComputeMassBalanceEquationResidu)(FVM_t* fvm,double const* f,double const* f_n,double const dt)
/** Return a mass balance equation residu */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double g[FVM_MaxNbOfNodes*FVM_MaxNbOfNodes] ;
    
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        int p = i*nn + j ;
        
        if(i == j) {
          g[p] = f[p] - f_n[p] ;
        } else {
          g[p] = dt * f[p] ;
        }
      }
    }
  }
    
  {
    double* r = FVM_ComputeMassAndFluxResidu(fvm,g,nn) ;
      
    return(r) ;
  }
}



double* (FVM_ComputeMassMatrix)(FVM_t* fvm,double* c,int neq)
/** Return a pointer on a FV mass matrix  (Ndof*Ndof) with
 *  Ndof = nb of d.o.f. (N*Neq), 
 *  N    = nb of nodes and 
 *  Neq  = nb of equations (conservation of mass).
 *  The inputs c is an array NxNeqxNeq i.e. containing the 
 *  coefficients of N matrices NeqxNeq and should be given 
 *  in the following order (A_i is a matrix NeqxNeq):
 *  the coefficients of A_1(NeqxNeq) should be given first then 
 *  those of A_2(NeqxNeq) etc.. up to A_N(NeqxNeq).
 *  Each coefficient (kl) of the matrix A_i should represent the 
 *  derivative of the mass content of equation k with respect to 
 *  the unknown l.
 */
{
#define K(i,j)     (k[(i)*nn*neq + (j)])
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int ndof = nn*neq ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* k = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;
  
  for(i = 0 ; i < nn ; i++) {
    /* cii is the derivative of the mass content at node i
     * with respect to unknown at node i: cii = M_i,i */
    double* cii = c + i*neq*neq ;
    int e ;
    for(e = 0 ; e < neq ; e++) {   /* Equations */
      int u ;
      for(u = 0 ; u < neq ; u++) { /* Unknowns */
        K(e + i*neq,u + i*neq) = volume[i]*cii[e*neq + u] ;
      }
    }
  }
  
  return(k) ;
  
#undef K
}



double*  (FVM_ComputeIsotropicConductionMatrix)(FVM_t* fvm,double* c,int neq)
/** Return a pointer on a FV conduction matrix (Ndof*Ndof) 
 *  for isotropic material, with
 *  Ndof = nb of d.o.f. (N*Neq), 
 *  N    = nb of nodes and 
 *  Neq  = nb of equations (conservation of mass).
 *  The inputs c is an array NxNxNeqxNeq i.e. containing the 
 *  coefficients of NxN matrices NeqxNeq and should be given 
 *  in the following order (A_ij is a matrix NeqxNeq):
 *  the coefficients of A_11(NeqxNeq) should be given first then 
 *  those of A_12(NeqxNeq) etc.. up to A_1N(NeqxNeq), then 
 *  those of A_21(NeqxNeq) etc.. up to A_NN(NeqxNeq).
 *  Each coefficient (kl) of the matrix A_ij should represent the 
 *  derivative of the outflow of equation k from cell i to cell j 
 *  with respect to the unknown l at cell i. Note that we don't need
 *  the derivative of this outflow with respect to unknown l at node j
 *  since it is the same (with opposite sign) as that located in (kl)
 *  of the matrix A_ji. Note also that coefficients of matrices A_ii
 *  are not used.
 */
{
#define K(i,j)     (k[(i)*ndof + (j)])
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int ndof = nn*neq ;
  double* surf = FVM_ComputeCellSurfaceAreas(fvm) ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* k = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < ndof*ndof ; i++)  k[i] = 0. ;

  for(i = 0 ; i < nn ; i++) {
    int j ;
    for(j = i + 1 ; j < nn ; j++) {
      /* cij[e*neq + u] is the derivative of the outflow wij
       * from node i to node j for equation (e),
       * with respect to unknown (u) at node i: wij,u[i].
       * Note that we don't need wij,u[j] since it is equal 
       * to - wji,u[j]. Note also that cii is not used */
      double* cij = c + (i*nn + j)*neq*neq ;
      double* cji = c + (j*nn + i)*neq*neq ;
      double a = surf[i*nn + j] ;
      int e ;
      for(e = 0 ; e < neq ; e++) {   /* Equations */
        int u ;
        for(u = 0 ; u < neq ; u++) { /* Unknowns */
          double trij = a*cij[e*neq + u] ;
          double trji = a*cji[e*neq + u] ;
          K(e + i*neq,u + i*neq) = + trij ;
          K(e + i*neq,u + j*neq) = - trji ;
          K(e + j*neq,u + j*neq) = + trji ;
          K(e + j*neq,u + i*neq) = - trij ;
        }
      }
    }
  }
  
  return(k) ;

#undef K
}



double*  (FVM_ComputeMassAndIsotropicConductionMatrix)(FVM_t* fvm,double* c,int neq)
/** Mass and Conduction Matrix for isotropic material (k) */
{
#define K(i,j)     (k[(i)*ndof + (j)])
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int ndof = nn*neq ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  double* k ;
  int    i ;

  /* Conduction matrix */
  k = FVM_ComputeIsotropicConductionMatrix(fvm,c,neq) ;
  
  /* Mass matrix */
  for(i = 0 ; i < nn ; i++) {
    double* cii = c + (i*nn + i)*neq*neq ;
    int e ;
    for(e = 0 ; e < neq ; e++) {   /* Equations */
      int u ;
      for(u = 0 ; u < neq ; u++) { /* Unknowns */
        K(e + i*neq,u + i*neq) += volume[i]*cii[e*neq + u] ;
      }
    }
  }
  
  return(k) ;

#undef K
}



double* (FVM_ComputeCellVolumes)(FVM_t* fvm)
{
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimension(el) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double* volume = FVM_GetCellVolumes(fvm) ;
  
  if(volume) {
    /* Previously computed */
    return(volume) ;
  } else {
    size_t SizeNeeded = nn*sizeof(double) ;
    volume = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
    FVM_GetCellVolumes(fvm) = volume ;
  }
   
  /* 0D */
  if(dim == 0) {
    if(nn == 1) {
      volume[0] = 1 ;
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double x = Element_GetNodeCoordinate(el,0)[0] ;
        volume[0] = 2*M_PI*x ;
        if(Symmetry_IsSpherical(sym)) volume[0] *= 4*x ;
      }
      return(volume) ;
    } else {
      arret("FVM_ComputeCellVolumes (1)") ;
    }
    
  /* 1D */
  } else if(dim == 1) {
    if(nn == 2) {
      double x1 = Element_GetNodeCoordinate(el,1)[0] ;
      double x0 = Element_GetNodeCoordinate(el,0)[0] ;
      double dx = x1 - x0 ;
      int i ;
      for(i = 0 ; i < 2 ; i++) {
        volume[i] = fabs(dx)*0.5 ; 
      }
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double xm = (x1 + x0)*0.5 ;
        for(i = 0 ; i < 2 ; i++) {
          double x  = Element_GetNodeCoordinate(el,i)[0] ;
          double dm = x + xm ;
          volume[i] *= M_PI*dm ; 
          if(Symmetry_IsSpherical(sym)) volume[i] *= dm ;
        }
      }
      return(volume) ;
    } else {
      arret("FVM_ComputeCellVolumes (1)") ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* Triangle */
    if(nn == 3) {
      double* x0 = Element_GetNodeCoordinate(el,0) ;
      double* x1 = Element_GetNodeCoordinate(el,1) ;
      double* x2 = Element_GetNodeCoordinate(el,2) ;
      double a = 0,b = 0,c = 0 ;
      int i ;
      
      for(i = 0 ; i < dim ; i++) {
        c += (x1[i] - x0[i])*(x1[i] - x0[i]) ;
        b += (x2[i] - x0[i])*(x2[i] - x0[i]) ;
        a += (x2[i] - x1[i])*(x2[i] - x1[i]) ;
      }
      a = 0.5*sqrt(a) ;
      b = 0.5*sqrt(b) ;
      c = 0.5*sqrt(c) ;
      
      { /* (a,b,c) are the half lengths of the sides of the triangle */
        double p = 0.5*(a + b + c) ; /* quarter of perimeter */
        double s = 4*sqrt(p*(p - a)*(p - b)*(p - c)) ; /* area */
        double r = 2*a*b*c/s ; /* radius of the circumcircle */
        double sab = (r > c) ? sqrt(r*r - c*c) : 0 ;
        double sbc = (r > a) ? sqrt(r*r - a*a) : 0 ;
        double sca = (r > b) ? sqrt(r*r - b*b) : 0 ;
        double va = 0.5*(b*sca + c*sab) ;
        double vb = 0.5*(c*sab + a*sbc) ;
        double vc = 0.5*(a*sbc + b*sca) ;
        
        volume[0] = va ;
        volume[1] = vb ;
        volume[2] = vc ;
        
        /* Does the circumcenter lie inside the triangle (acute triangle) ?
         * Only if the length of any median is greater than the circumradius
         * The median ma is given by ma^2 = 2(a^2 + b^2 + c^2) - 3a^2 */
        {
          double a2 = a*a ;
          double b2 = b*b ;
          double c2 = c*c ;
          double m2 = (2*(a2 + b2 + c2) - r*r) ;
          
          if(3*a2 > m2 || 3*b2 > m2 || 3*c2 > m2) {
            int iel = Element_GetElementIndex(el) ;
            arret("FVM_ComputeCellVolumes: circumcenter not inside the triangle %d !",iel) ;
          }
        }
      }
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        arret("FVM_ComputeCellVolumes (2)") ;
      }
      return(volume) ;
      
    /* Quadrangle */
    } else if(nn == 4) {
      
      arret("FVM_ComputeCellVolumes: not yet done") ;
    } else {
      arret("FVM_ComputeCellVolumes: ") ;
    }
    
  /* 3D */
  } else if(dim == 3) {
    arret("FVM_ComputeCellVolumes (4)") ;
  }
  
  arret("FVM_ComputeCellVolumes") ;
  return(NULL) ;
}


double* (FVM_ComputeCellSurfaceAreas)(FVM_t* fvm)
{
#define AREA(i,j)     area[nn*(i) + (j)]
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimension(el) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double* area = FVM_GetCellSurfaceAreas(fvm) ;
  int i ;
  
  if(area) {
    /* Previously computed */
    return(area) ;
  } else {
    size_t SizeNeeded = nn*nn*sizeof(double) ;
    area = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
    FVM_GetCellSurfaceAreas(fvm) = area ;
  }
  
  /* The diagonal terms are not used */
  for(i = 0 ; i < nn ; i++) {
    area[nn*i + i] = 0. ;
  }
  
  
  /* 0D */
  if(dim == 0) {
    if(nn == 1) {
      /* No areas */
      return(area) ;
    } else {
      arret("FVM_ComputeCellSurfaceAreas (0)") ;
    }
    
  /* 1D */
  } else if(dim == 1) {
    if(nn == 2) {
      
      AREA(0,1) = 1 ;
      
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double x1 = Element_GetNodeCoordinate(el,1)[0] ;
        double x0 = Element_GetNodeCoordinate(el,0)[0] ;
        double dm = x1 + x0 ;
        AREA(0,1) = M_PI*dm ;
        if(Symmetry_IsSpherical(sym)) AREA(0,1) *= dm ;
      }
      
      AREA(1,0) = AREA(0,1) ;
  
      return(area) ;
    } else {
      arret("FVM_ComputeCellSurfaceAreas (1)") ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* Triangle */
    if(nn == 3) {
      double* x0 = Element_GetNodeCoordinate(el,0) ;
      double* x1 = Element_GetNodeCoordinate(el,1) ;
      double* x2 = Element_GetNodeCoordinate(el,2) ;
      double a = 0,b = 0,c = 0 ;
      
      for(i = 0 ; i < dim ; i++) {
        c += (x1[i] - x0[i])*(x1[i] - x0[i]) ;
        b += (x2[i] - x0[i])*(x2[i] - x0[i]) ;
        a += (x2[i] - x1[i])*(x2[i] - x1[i]) ;
      }
      a = 0.5*sqrt(a) ;
      b = 0.5*sqrt(b) ;
      c = 0.5*sqrt(c) ;
      
      { /* (a,b,c) are the half lengths of the sides of the triangle */
        double p = 0.5*(a + b + c) ; /* quarter of perimeter */
        double s = 4*sqrt(p*(p - a)*(p - b)*(p - c)) ; /* area */
        double r = 2*a*b*c/s ; /* radius of the circumcircle */
        double sab = (r > c) ? sqrt(r*r - c*c) : 0 ;
        double sbc = (r > a) ? sqrt(r*r - a*a) : 0 ;
        double sca = (r > b) ? sqrt(r*r - b*b) : 0 ;
        
        AREA(0,1) = (AREA(1,0) = sab) ;
        AREA(1,2) = (AREA(2,1) = sbc) ;
        AREA(2,0) = (AREA(0,2) = sca) ;
        
        /* Does the circumcenter lie inside the triangle (acute triangle) ?
         * Only if the length of any median is greater than the circumradius.
         * The median ma is given by ma^2 = 2(a^2 + b^2 + c^2) - 3a^2 */
        {
          double a2 = a*a ;
          double b2 = b*b ;
          double c2 = c*c ;
          double m2 = (2*(a2 + b2 + c2) - r*r) ;
          
          if(3*a2 > m2 || 3*b2 > m2 || 3*c2 > m2) {
            arret("FVM_ComputeCellSurfaceAreas: circumcenter not inside the triangle !") ;
          }
        }
      }
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        arret("FVM_ComputeCellSurfaceAreas (2)") ;
      }
      return(area) ;
    } else {
      arret("FVM_ComputeCellSurfaceAreas (3)") ;
    }
    
  /* 3D */
  } else if(dim == 3) {
    arret("FVM_ComputeCellSurfaceAreas (4)") ;
  }
  
  arret("FVM_ComputeCellSurfaceAreas") ;
  return(NULL) ;
#undef AREA
}


double* (FVM_ComputeCellVolumesAndSurfaceAreas)(FVM_t* fvm)
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn = Element_GetNbOfNodes(el) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  double* area = FVM_ComputeCellSurfaceAreas(fvm) ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    area[nn*i + i] = volume[i] ;
  }
  
  return(area) ;
}


short int   (FVM_FindLocalCellIndex)(FVM_t* fvm,double* s)
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn = Element_GetNbOfNodes(el) ;
  
  if(nn == 0) {
    return(0) ;
  } else {
    int i,j = 0 ;
    
    for(i = 1 ; i < nn ; i++) {
      j = FVM_FindHalfSpace(fvm,i,j,s) ;
    }
    return(j) ;
  }
  
  arret("FVM_FindLocalCellIndex") ;
  return(-1) ;
}



double* (FVM_ComputeIntercellDistances)(FVM_t* fvm)
/** Compute the intercell distances */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* dist = FVM_GetIntercellDistances(fvm) ;
  int    i ;
  
  if(dist) {
    /* Previously computed */
    return(dist) ;
  } else {
    size_t SizeNeeded = nn*nn*(sizeof(double)) ;
    dist = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
    FVM_GetIntercellDistances(fvm) = dist ;
  }
  
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



double* (FVM_ComputeTheNodalFluxVector)(FVM_t* fvm,double* w)
/** Return the flux vectors at nodes. "w" points to an
 *  array of NxN doubles where N is the number of nodes and where
 *  w[i*N + j] denotes the outflow rate from cell i to cell j. */
{
  #define W(i,j)    w[(i)*nn + (j)]
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimension(el) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  size_t SizeNeeded = 3*nn*sizeof(double) ;
  double* wn = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  
  {
    int i ;
    
    for(i = 0 ; i < 3*nn ; i++) wn[i] = 0. ;
  }
  
  
  /* 0D */
  if(dim == 0) {
    return(wn) ;
    
  /* 1D */
  } else if(dim == 1) {
    /* Line */
    if(nn == 2) {
      double* x0 = Element_GetNodeCoordinate(el,0) ;
      double* x1 = Element_GetNodeCoordinate(el,1) ;
      
      double a = 0 ;
      int i ;
      
      for(i = 0 ; i < dim ; i++) {
        a += (x1[i] - x0[i])*(x1[i] - x0[i]) ;
      }
      a = sqrt(a) ;
      
      wn[0] = W(0,1)*(x1[0] - x0[0])/(a) ;
      wn[3] = W(1,0)*(x0[0] - x1[0])/(a) ;

      return(wn) ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* Triangle */
    if(nn == 3) {
      double* x0 = Element_GetNodeCoordinate(el,0) ;
      double* x1 = Element_GetNodeCoordinate(el,1) ;
      double* x2 = Element_GetNodeCoordinate(el,2) ;
      
      double a = 0,b = 0,c = 0 ;
      int i ;
      
      for(i = 0 ; i < dim ; i++) {
        c += (x1[i] - x0[i])*(x1[i] - x0[i]) ;
        b += (x2[i] - x0[i])*(x2[i] - x0[i]) ;
        a += (x2[i] - x1[i])*(x2[i] - x1[i]) ;
      }
      a = sqrt(a) ;
      b = sqrt(b) ;
      c = sqrt(c) ;
      
      { /* (a,b,c) are the lengths of the sides of the triangle */
        double a2 = a*a ;
        double b2 = b*b ;
        double c2 = c*c ;
        double cosA = 0.5*(b2 + c2 - a2)/(b*c) ;
        double cosB = 0.5*(c2 + a2 - b2)/(c*a) ;
        double cosC = 0.5*(a2 + b2 - c2)/(a*b) ;
        double** x = Element_ComputePointerToNodalCoordinates(el) ;
        int k ;

        for(k = 0 ; k < dim ; k++) {
          wn[k]  = (W(0,1) - W(0,2)*cosA)*(x[1][k] - x[0][k])/(c) \
                 + (W(0,2) - W(0,1)*cosA)*(x[2][k] - x[0][k])/(b) ;
          wn[k] /= (1 - cosA*cosA) ;
        }
        for(k = 0 ; k < dim ; k++) {
          wn[3 + k]  = (W(1,0) - W(1,2)*cosB)*(x[0][k] - x[1][k])/(c) \
                     + (W(1,2) - W(1,0)*cosB)*(x[2][k] - x[1][k])/(a) ;
          wn[3 + k] /= (1 - cosB*cosB) ;
        }
        for(k = 0 ; k < dim ; k++) {
          wn[6 + k]  = (W(2,0) - W(2,1)*cosC)*(x[0][k] - x[2][k])/(b) \
                     + (W(2,1) - W(2,0)*cosC)*(x[1][k] - x[2][k])/(a) ;
          wn[6 + k] /= (1 - cosC*cosC) ;
        }
      }
      
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        arret("FVM_ComputeTheNodalFluxVector(2)") ;
      }

      return(wn) ;
    }
    
  /* 3D */
  } else if(dim == 3) {
    /* Tetrahedron */
    if(nn == 4) {
     
    arret("FVM_ComputeTheNodalFluxVector(4)") ;
    }
  }
  
  arret("FVM_ComputeTheNodalFluxVector") ;
  return(NULL) ;
  #undef W
}



double   (FVM_AverageCurrentImplicitTerm)(Mesh_t* mesh,const char* modelname,const int index,const int shift)
/** Volume averaging, over the mesh, the value stored at node i and at the position
 *  "index+i*shift" in the implicit terms of the model "modelname".
 */
{
  unsigned int nel = Mesh_GetNbOfElements(mesh) ;
  Element_t* el0 = Mesh_GetElement(mesh) ;
  double sum = 0 ;
  double vol = 0 ;
  
    /* Integration */
  {
    unsigned int ie ;
    
    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* el = el0 + ie ;
      Model_t* model = Element_GetModel(el) ;
      char* codename = Model_GetCodeNameOfModel(model) ;
    
      if(Element_IsSubmanifold(el)) continue ;
      
      if(String_Is(codename,modelname)) {
        double* vim = Element_GetCurrentImplicitTerm(el) ;
        int ni = Element_GetNbOfImplicitTerms(el) ;
        int nn = Element_GetNbOfNodes(el) ;
        FVM_t*    fvm    = FVM_GetInstance(el) ;
        double* volume = FVM_ComputeCellVolumes(fvm) ;
        double* vi = vim + index ;
        int i ; ;
        
        if(index + (nn-1)*shift >= ni) {
          arret("FVM_AverageCurrentImplicitTerm:") ;
        }
          
        for(i = 0 ; i < nn ; i++) {
          vol +=  volume[i] ;
          sum +=  vi[i*shift] ;
        }
      }
    }
    
    if(vol <= 0) {
      arret("FVM_AverageCurrentImplicitTerm:") ;
    }
  }

  return(sum/vol) ;
}



double   (FVM_AveragePreviousImplicitTerm)(Mesh_t* mesh,const char* modelname,const int index,const int shift)
/** Volume averaging, over the mesh, the value stored at node i and at the position
 *  "index+i*shift" in the implicit terms of the model "modelname".
 */
{
  unsigned int nel = Mesh_GetNbOfElements(mesh) ;
  Element_t* el0 = Mesh_GetElement(mesh) ;
  double sum = 0 ;
  double vol = 0 ;
  
    /* Integration */
  {
    unsigned int ie ;
    
    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* el = el0 + ie ;
      Model_t* model = Element_GetModel(el) ;
      char* codename = Model_GetCodeNameOfModel(model) ;
    
      if(Element_IsSubmanifold(el)) continue ;
      
      if(String_Is(codename,modelname)) {
        double* vim = Element_GetPreviousImplicitTerm(el) ;
        int ni = Element_GetNbOfImplicitTerms(el) ;
        int nn = Element_GetNbOfNodes(el) ;
        FVM_t*    fvm    = FVM_GetInstance(el) ;
        double* volume = FVM_ComputeCellVolumes(fvm) ;
        double* vi = vim + index ;
        int i ; ;
        
        if(index + (nn-1)*shift >= ni) {
          arret("FVM_AveragePreviousImplicitTerm:") ;
        }
          
        for(i = 0 ; i < nn ; i++) {
          vol +=  volume[i] ;
          sum +=  vi[i*shift] ;
        }
      }
    }
    
    if(vol <= 0) {
      arret("FVM_AveragePreviousImplicitTerm:") ;
    }
  }

  return(sum/vol) ;
}



double* (FVM_ComputeGradient)(FVM_t* fvm,double* u,IntFct_t* intfct,int p,int shift)
/** Compute the gradient for a nodal quantity u at the node point "p".
 *  Input:
 *    "shift": the nodal quantity is saved at i*shift
 *    "intfct" is expected to be interpolation functions at nodes.
 *  Output:
 *    a pointer to double[3] for the u gradient.
*/
{
#define U(n)   (u[n*shift])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int dim_e = Element_GetDimension(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
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




/*
 * Intern Functions
 */

int (FVM_FindHalfSpace)(FVM_t* fvm,int i1,int i2,double* s)
{
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double* x1 = Element_ComputePointerToNodalCoordinates(el)[i1] ;
  double* x2 = Element_ComputePointerToNodalCoordinates(el)[i2] ;
  double cos1 = 0 ;
  double cos2 = 0 ;
  int i ;
  
  for(i = 0 ; i < dim ; i++) {
    cos1 += (x2[i] - x1[i])*(s[i] - x1[i]) ;
    cos2 += (x1[i] - x2[i])*(s[i] - x2[i]) ;
  }
  
  return((cos1 < cos2) ? i1 : i2) ;
}



/*
 * Not used Functions
 */

#if 0
static double*    (FVM_ComputeRelativeCellSurfaceAreas)(FVM_t*) ;
static double*    (FVM_ComputeNormalGradientMatrix)(FVM_t*,int,double*) ;
static double*    (FVM_ComputeOutFlowMatrix)(FVM_t*,int,double*) ;



double* (FVM_ComputeRelativeCellSurfaceAreas)(FVM_t* fvm)
/** Compute the surface area gradients */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* surface = FVM_ComputeCellSurfaceAreas(fvm) ;
  double* flow = surface ;
  int    i ;

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
      flow[nn*i + j] /= dij ;
      flow[nn*j + i] /= dij ;
    }
  }
  return(flow) ;
}



double* (FVM_ComputeNormalGradientMatrix)(FVM_t* fvm,int shift,double* c)
/** Compute the c gradients normal to the cell surfaces */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* grad = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < nn*nn ; i++)  grad[i] = 0. ;

  for(i = 0 ; i < nn ; i++) {
    double* xi = Element_GetNodeCoordinate(el,i) ;
    double ci = c[shift*i] ;
    int j ;
    for(j = i + 1 ; j < nn ; j++) {
      double* xj = Element_GetNodeCoordinate(el,j) ;
      double cj = c[shift*j] ;
      double dij = 0 ;
      int n ;
      for (n = 0 ; n < dim ; n++) {
        double xij = xj[n] - xi[n] ;
        dij += xij*xij ;
      }
      dij = sqrt(dij) ;
      grad[nn*i + j] = (cj - ci)/dij ;
      grad[nn*j + i] = - grad[nn*i + j] ;
    }
  }
  
  return(grad) ;
}



double* (FVM_ComputeOutFlowMatrix)(FVM_t* fvm,int shift,double* c)
/** Compute the flow matrix from the cell values in c */
{
  Element_t* el = FVM_GetElement(fvm) ;
  int nn  = Element_GetNbOfNodes(el) ;
  double* surf = FVM_ComputeCellSurfaceAreas(fvm) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* flow = (double*) FVM_AllocateInBuffer(fvm,SizeNeeded) ;
  int    i ;
  
  /* initialization */
  for(i = 0 ; i < nn*nn ; i++)  flow[i] = 0. ;

  for(i = 0 ; i < nn ; i++) {
    double ci = c[shift*i] ;
    int j ;
    for(j = i + 1 ; j < nn ; j++) {
      double cj = c[shift*j] ;
      double a = surf[i*nn + j] ;
      flow[nn*i + j] = a*(cj - ci) ;
      flow[nn*j + i] = - flow[nn*i + j] ;
    }
  }
  
  return(flow) ;
}
#endif
