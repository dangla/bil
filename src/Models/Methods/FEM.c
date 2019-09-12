#include "FEM.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Geometry.h"
#include "Elements.h"
#include "ShapeFcts.h"
#include "IntFcts.h"
#include "Nodes.h"
#include "Models.h"
#include "Buffer.h"
#include "Session.h"
#include "GenericData.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <assert.h>


//static FEM_t* instancefem = NULL ;

static FEM_t*  FEM_Create(void) ;
static double* FEM_ComputeJacobianMatrix(Element_t*,double*,int) ;
static double* FEM_ComputeInverseJacobianMatrix(Element_t*,double*,int) ;
static double  FEM_ComputeJacobianDeterminant(Element_t*,double*,int) ;
//static double* FEM_ComputeNormalVector(Element_t*,double*,int) ;

/* 
   Extern Functions 
*/
FEM_t* FEM_Create(void)
{
  FEM_t*  fem    = (FEM_t*) malloc(sizeof(FEM_t)) ;
  
  if(!fem) assert(fem) ;
  
  
  /* Space allocation for output */
  {
    int nn = Element_MaxNbOfNodes ;
    int neq = Model_MaxNbOfEquations ;
    int ndof = nn*neq ;
    size_t sz = ndof*ndof*sizeof(double) ;
    double* output = (double*) malloc(sz) ;
    
    if(!output) arret("FEM_Create") ;
    
    FEM_GetOutput(fem) = output ;
  }
  
  
  /* Space allocation for input */
  {
    int np = IntFct_MaxNbOfIntPoints ;
    size_t sz = np*FEM_MaxShift*sizeof(double) ;
    double* input = (double*) malloc(sz) ;
    
    if(!input) arret("FEM_Create") ;
    
    FEM_GetInput(fem) = input ;
  }
  
  
  /* Space allocation for pintfct */
  {
    size_t sz = IntFcts_MaxNbOfIntFcts*sizeof(IntFct_t*) ;
    IntFct_t** pintfct = (IntFct_t**) malloc(sz) ;
    
    if(!pintfct) arret("FEM_Create") ;
    
    FEM_GetPointerToIntFct(fem) = pintfct ;
  }
  
  
  /* Space allocation for buffer */
  {
    Buffer_t* buf = Buffer_Create(FEM_SizeOfBuffer) ;
    
    FEM_GetBuffer(fem) = buf ;
  }
  
  return(fem) ;
}



void FEM_Delete(void* self)
{
  FEM_t** pfem = (FEM_t**) self ;
  FEM_t*   fem = *pfem ;
  
  free(FEM_GetOutput(fem)) ;
  free(FEM_GetInput(fem)) ;
  free(FEM_GetPointerToIntFct(fem)) ;
  Buffer_Delete(&FEM_GetBuffer(fem))  ;
  free(fem) ;
  *pfem = NULL ;
}


#if 0
FEM_t* FEM_GetInstance0(Element_t* el)
{
  if(!instancefem) {
    instancefem = FEM_Create() ;
  }
  
  FEM_GetElement(instancefem) = el ;
  
  FEM_FreeBuffer(instancefem) ;
  
  return(instancefem) ;
}
#endif



FEM_t*  (FEM_GetInstance)(Element_t* el)
{
  GenericData_t* gdat = Session_FindGenericData(FEM_t,"FEM") ;
  
  if(!gdat) {
    FEM_t* fem = FEM_Create() ;
    
    gdat = GenericData_Create(1,fem,FEM_t,"FEM") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(FEM_t,"FEM")) ;
  }
  
  {
    FEM_t* fem = (FEM_t*) GenericData_GetData(gdat) ;
  
    FEM_GetElement(fem) = el ;
  
    FEM_FreeBuffer(fem) ;
  
    return(fem) ;
  }
}



double*  FEM_ComputeElasticMatrix(FEM_t* fem,IntFct_t* fi,const double* c,const int dec)
/** Return a pointer on a FE elastic matrix (Ndof*Ndof)
 *  with Ndof = N * dim and N = nb of nodes */
{
#define KR(i,j)     (kr[(i)*ndof + (j)])
#define DH(n,i)     (dh[(n)*dim + (i)])
#define C(i,j,k,l)  (c[(((i)*3 + (j))*3 + (k))*3 + (l)])
#define CAJ(i,j)    (caj[(i)*dim + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = IntFct_GetNbOfNodes(fi) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int ndof = nn*dim ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* kr = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0., deux = 2. ;
  
  if(Element_IsSubmanifold(el)) arret("FEM_ComputeElasticMatrix") ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) kr[i] = zero ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , c += dec) {
    double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nn) ;
    double jcj[3][3][3][3],jc[3][3],cj[3][3] ;
    double rayon = zero ;
    int    j,k,l,r,s ;
    
    /* cas axisymetrique ou spherique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= deux*rayon ;
    }
    
    /* JCJ(r,i,k,s) = J(r,j)*C(i,j,k,l)*J(s,l) */
    for(r = 0 ; r < dim ; r++) for(i = 0 ; i < dim ; i++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim ; s++) {
      jcj[r][i][k][s] = zero ;
      for(j = 0 ; j < dim ; j++) for(l = 0 ; l < dim ; l++) {
        jcj[r][i][k][s] += CAJ(r,j)*C(i,j,k,l)*CAJ(s,l) ;
      }
    }
    
    /* KR(j,i,l,k) = DH(j,r)*JCJ(r,i,k,s)*DH(l,s) */
    for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(k = 0 ; k < dim ; k++) {
      for(r = 0 ; r < dim ; r++) for(s = 0 ; s < dim ; s++) {
        KR(j*dim+i,l*dim+k) += a*DH(j,r)*jcj[r][i][k][s]*DH(l,s) ;
      }
    }
    
    /* cas axisymetrique: 3 termes */
    if(Symmetry_IsCylindrical(sym)) {
      /* 1.a JC(r,i) = J(r,j)*C(i,j,theta,theta) */
      for(r = 0 ; r < dim ; r++) for(i = 0 ; i < dim ; i++) {
        jc[r][i] = zero ;
        for(j = 0 ; j < dim ; j++) jc[r][i] += CAJ(r,j)*C(i,j,2,2) ;
      }
      /* 1.b KR(j,i,l,0) = DH(j,r)*JC(r,i)*H(l)/r */
      for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(r = 0 ; r < dim ; r++) {
        KR(j*dim+i,l*dim) += a*DH(j,r)*jc[r][i]*h[l]/rayon ;
      }
      /* 2.a CJ(k,s) = C(theta,theta,k,l)*J(s,l) */
      for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim ; s++) {
        cj[k][s] = zero ; 
        for(l = 0 ; l < dim ; l++) cj[k][s] += C(2,2,k,l)*CAJ(s,l) ;
      }
      /* 2.b KR(j,0,l,k) = H(j)/r*CJ(k,s)*DH(l,s) */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn;l++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim ; s++) {
        KR(j*dim,l*dim+k) += a*h[j]/rayon*cj[k][s]*DH(l,s) ;
      }
      /* 3.  KR(j,0,l,0) = H(j)/r*C(theta,theta,theta,theta)*H(l)/r */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
        KR(j*dim,l*dim) += a*h[j]/rayon*C(2,2,2,2)*h[l]/rayon ;
      }
      
    /* cas spherique: 3 termes */
    } else if(Symmetry_IsSpherical(sym)) {
      double cc = C(1,1,1,1) + C(1,1,2,2) + C(2,2,1,1) + C(2,2,2,2) ;
      /* 1.a JC(r,i) = J(r,j)*(C(i,j,theta,theta)+C(i,j,phi,phi)) */
      for(r = 0 ; r < dim ; r++) for(i = 0 ; i < dim ; i++) {
        jc[r][i] = zero ;
        for(j = 0 ; j < dim ; j++) jc[r][i] += CAJ(r,j)*(C(i,j,1,1) + C(i,j,2,2)) ;
      }
      /* 1.b KR(j,i,l,0) = DH(j,r)*JC(r,i)*H(l)/r */
      for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(r = 0 ; r < dim ; r++) {
        KR(j*dim+i,l*dim) += a*DH(j,r)*jc[r][i]*h[l]/rayon ;
      }
      /* 2.a CJ(k,s) = (C(phi,phi,k,l)+C(theta,theta,k,l))*J(s,l) */
      for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim ; s++) {
        cj[k][s] = zero ; 
        for(l = 0 ; l < dim ; l++) cj[k][s] += (C(1,1,k,l) + C(2,2,k,l))*CAJ(s,l) ;
      }
      /* 2.b KR(j,0,l,k) = H(j)/r*CJ(k,s)*DH(l,s) */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim ; s++) {
        KR(j*dim,l*dim+k) += a*h[j]/rayon*cj[k][s]*DH(l,s) ;
      }
      /* 3.  KR(j,0,l,0) = H(j)/r*(C(phi,phi,phi,phi)+C(theta,theta,phi,phi)+C(phi,phi,theta,theta)+C(theta,theta,theta,theta))*H(l)/r */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
        KR(j*dim,l*dim) += a*h[j]/rayon*cc*h[l]/rayon ;
      }
    }
  }
  
  return(kr) ;
#undef KR
#undef DH
#undef C
#undef CAJ
}



double*  FEM_ComputeBiotMatrix(FEM_t* fem,IntFct_t* fi,const double* c,const int dec)
/** Return a pointer on a Biot-like coupling matrix (Ndof*N) 
 *  with Ndof = N * dim and N = nb of nodes */
{
#define KC(i,j)     (kc[(i)*nn + (j)])
#define DH(n,i)     (dh[(n)*dim + (i)])
#define C(i,j)      (c[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*dim + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nn  = IntFct_GetNbOfNodes(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int nrow = nn*dim ;
  int ncol = nn ;
  size_t SizeNeeded = nrow*ncol*(sizeof(double)) ;
  double* kc = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,j,p ;
  double zero = 0. ;
  
  if(Element_IsSubmanifold(el)) arret("FEM_ComputeBiotMatrix") ;
  
  /* initialisation */
  for(i = 0 ; i < nrow*ncol ; i++) kc[i] = zero ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , c += dec) {
    double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nn) ;
    double jc[3][3] ;
    double rayon ;
    int    k,l ;
    
    /* cas axisymetrique ou spherique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      rayon = zero ;
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    /* JC(k,i) = J(k,j)*C(j,i) */
    for(k = 0 ; k < dim ; k++) for(i = 0 ; i < dim ; i++) {
      jc[k][i] = zero ;
      for(j = 0 ; j < dim ; j++) jc[k][i] += CAJ(k,j)*C(j,i) ;
    }
    /* KC(j,i,l) = DH(j,k)*JC(k,i)*H(l) */
    for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(k = 0 ; k < dim ; k++) {
      KC(j*dim+i,l) += a*DH(j,k)*jc[k][i]*h[l] ;
    }
    /* cas axisymetrique: (r,z,theta) */
    if(Symmetry_IsCylindrical(sym)) {
      /* KC(j,0,l) = H(j)/r*C(theta,theta)*H(l) */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
        KC(j*dim,l) += a*h[j]/rayon*C(2,2)*h[l] ;
      }
    /* cas spherique: (r,theta,phi) */
    } else if(Symmetry_IsSpherical(sym)) {
      /* KC(j,0,l) = H(j)/r*(C(theta,theta)+C(phi,phi))*H(l) */
      for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
        KC(j*dim,l) += a*h[j]/rayon*(C(1,1)+C(2,2))*h[l] ;
      }
    }
  }
  
  return(kc) ;
  
#undef KC
#undef DH
#undef C
#undef CAJ
}




double* FEM_ComputeMassMatrix(FEM_t* fem,IntFct_t* fi,const double* c,const int dec)
/** Return a pointer on a FE mass matrix (N*N) 
 *  with N = nb of nodes*/
{
#define KM(i,j)     (km[(i)*nn + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nn  = IntFct_GetNbOfNodes(fi) ;
  int dim_h = IntFct_GetDimension(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* km = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0. ;
  
  if(Element_IsSubmanifold(el)) arret("FEM_ComputeMassMatrix") ;
  
  /* initialisation */
  for(i = 0 ; i < nn*nn ; i++) {
    km[i] = zero ;
  }

  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  
  /* 0D */
  if(dim_h == 0) {
    if(nn == 1) {
      KM(0,0) = c[0] ;
      if(Symmetry_IsCylindrical(sym)) KM(0,0) *= 2*M_PI*x[0][0] ;
      else if(Symmetry_IsSpherical(sym)) KM(0,0) *= 4*M_PI*x[0][0]*x[0][0] ;
    } else arret("FEM_ComputeMassMatrix : impossible") ;
    return(km) ;
  }

  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , c += dec) {
    double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*c[0]*d ;
    int    j ;
    
    /* cas axisymetrique ou shperique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      double rayon = zero ;
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      KM(i,j) += a*h[i]*h[j] ;
    }
  }
  
  return(km) ;
  
#undef KM
}





double*  FEM_ComputeConductionMatrix(FEM_t* fem,IntFct_t* fi,const double* c,const int dec)
/** Return a pointer on a FE conduction matrix (N*N) 
 *  with N = nb of nodes */
{
#define KC(i,j)     (kc[(i)*nn + (j)])
#define DH(n,i)     (dh[(n)*dim + (i)])
#define C(i,j)      (c[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*dim + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nn  = IntFct_GetNbOfNodes(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* kc = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0. ;
  
  if(Element_IsSubmanifold(el)) arret("FEM_ComputeConductionMatrix") ;
  
  /* initialisation */
  for(i = 0 ; i < nn*nn ; i++) {
    kc[i] = zero ;
  }
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , c += dec) {
    double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nn) ;
    double jcj[3][3] ;
    int    j,k,l ;
    
    /* cas axisymetrique ou shperique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      double rayon = zero ;
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    /* jcj = J(i,k)*C(k,l)*J(j,l) */
    for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
      jcj[i][j] = zero ;
      for(k = 0 ; k < dim ; k++) for(l = 0 ; l < dim ; l++)  {
        jcj[i][j] += CAJ(i,k)*C(k,l)*CAJ(j,l) ;
      }
    }
    /* KC(i,j) = DH(i,k)*JCJ(k,l)*DH(j,l) */
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      for(k = 0 ; k < dim ; k++) for(l = 0 ; l < dim ; l++) {
        KC(i,j) += a*DH(i,k)*jcj[k][l]*DH(j,l) ;
      }
    }
  }
  
  return(kc) ;

#undef KC
#undef DH
#undef C
#undef CAJ
}



double*  FEM_ComputePoroelasticMatrix(FEM_t* fem,IntFct_t* fi,const double* c,const int dec,const int n_dif)
/** Return a pointer on a FE poroelastic matrix  (Ndof*Ndof) 
 *  with Ndof = N*Neq, N = nb of nodes and Neq = dim + n_dif 
 *  The inputs c should be given in the following order 
 *  A0 to An, B0 to Bn etc.. with n = n_dif:
 *  | A0(9x9) A1(9x1) A2(9x1) ... | 
 *  | B0(9x9) B1(9x1) B2(9x1) ... | 
 *  | C0(1x9) C1(1)   C2(1)   ... |
 *  |  .......................... |
 */
{
#define K(i,j)      (k[(i)*ndof + (j)])
#define E_Mec       (0)
#define I_U         (0)
#define E_Hyd       (dim)
#define I_H         (dim)
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = IntFct_GetNbOfNodes(fi) ;
  int neq = dim + n_dif ;
  int ndof = nn*neq ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* k = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  int    i,n,m ;
  double zero = 0. ;
  
  if(Element_IsSubmanifold(el)) arret("FEM_ComputePoroelasticMatrix") ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* 
  ** 1.  Mechanics
  */
  /* 1.1 Elasticity */
  {
    double* kr = FEM_ComputeElasticMatrix(fem,fi,c,dec) ;
    
    for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
      int j ;
      
      for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
        K(E_Mec + i + n*neq,I_U + j + m*neq) = kr[(i+n*dim)*nn*dim + j + m*dim] ;
      }
    }
  }
  
  /* 1.2 Biot-like Coupling terms */
  for(i = 0 ; i < n_dif ; i++) {
    int I_h = I_H + i ;
    const double* c1 = c + 81 +i*9 ;
    double* kb = FEM_ComputeBiotMatrix(fem,fi,c1,dec) ;
    
    for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        K(E_Mec + j + n*neq,I_h + m*neq) = kb[(j+n*dim)*nn + m] ;
      }
    }
  }
  

  /* 2. Hydraulic equations */
  for(i = 0 ; i < n_dif ; i++) {
    int E_h = E_Hyd + i ;
    int j ;
    
    /* 2.1 Coupling terms */
    {
      const double* c1 = c + 81 + n_dif*9 + i*(9 + n_dif) ;
      double* kb = FEM_ComputeBiotMatrix(fem,fi,c1,dec) ;
    
      for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
        for(j = 0 ; j < dim ; j++) {
          K(E_h + m*neq,I_U + j + n*neq) = kb[(j + n*dim)*nn + m] ;
        }
      }
    }
    
    /* 2.2 Accumulations */
    for(j = 0 ; j < n_dif ; j++) {
      int I_h = I_H + j ;
      const double* c2 = c + 81 + n_dif*9 + i*(9 + n_dif) + 9 + j ;
      double* ka = FEM_ComputeMassMatrix(fem,fi,c2,dec) ;
      
      for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
        K(E_h + n*neq,I_h + m*neq) = ka[n*nn + m] ;
      }
    }
  }
  
  return(k) ;

#undef K
#undef E_Mec
#undef I_U
#undef E_Hyd
#undef I_H
}



void FEM_TransformMatrixFromDegree2IntoDegree1(FEM_t* fem,const int inc,const int equ,double* k)
{
#define K(i,j)     (k[(i)*nn*neq+(j)])
  Element_t* el = FEM_GetElement(fem) ;
  int    nn = Element_GetNbOfNodes(el) ;
  int    dim = Element_GetDimension(el) ;
  int    neq = Element_GetNbOfEquations(el) ;
  int    n_vertices = 0 ;
  int    edge_vertices_line[2]       = {0,1} ;
  int    edge_vertices_triangle[6]   = {0,1,1,2,2,0} ;
  int    edge_vertices_quadrangle[8] = {0,1,1,2,2,3,3,0} ;
  int    *edge_vertices[Element_MaxNbOfNodes] ;
  int    ilin,icol ;
  int    n ;

  if(nn > Element_MaxNbOfNodes) arret("FEM_TransformMatrixFromDegree2IntoDegree1") ;

  if(dim == 1) {
    if(nn == 3) {
      n_vertices = 2 ;
      edge_vertices[n_vertices] = edge_vertices_line ;
    } else return ;
  } else if(dim == 2) {
    if(nn == 6) {
      n_vertices = 3 ;
      edge_vertices[n_vertices] = edge_vertices_triangle ;
    } else if(nn == 8) {
      n_vertices = 4 ;
      edge_vertices[n_vertices] = edge_vertices_quadrangle ;
    } else return ;
  } else if(dim == 3) {
    arret("FEM_TransformMatrixFromDegree2IntoDegree1") ;
  } else {
    arret("FEM_TransformMatrixFromDegree2IntoDegree1") ;
  }

  for(n = n_vertices + 1 ; n < nn ; n++) {
    edge_vertices[n] = edge_vertices[n - 1] + 2 ;
  }

  /* Transformation de la matrice de degre 2 en degre 1 */
  for(ilin = 0 ; ilin < nn*neq ; ilin++) { /* boucle sur les lignes */
    for(n = n_vertices ; n < nn ; n++) {
      int i ;
      icol  = n*neq + inc ;
      for(i = 0 ; i < 2 ; i++) {
	      int    m    = edge_vertices[n][i] ;
	      int    jcol = m*neq + inc ;
	      K(ilin,jcol) += 0.5*K(ilin,icol) ;
      }
      K(ilin,icol) = 0. ;
    }
  }
  
  for(icol = 0 ; icol < nn*neq ; icol++) { /* boucle sur les colonnes */
    for(n = n_vertices ; n < nn ; n++) {
      int i ;
      ilin  = n*neq + equ ;
      for(i = 0 ; i < 2 ; i++) {
	      int    m     = edge_vertices[n][i] ;
	      int    jlin  = m*neq + equ ;
	      K(jlin,icol) += 0.5*K(ilin,icol) ;
      }
      K(ilin,icol) = 0. ;
    }
  }

  /* Relation lineaire aux noeuds milieux des aretes */
  for(n = n_vertices ; n < nn ; n++) {
    int i ;
    ilin  = n*neq + equ ;
    K(ilin,ilin) = -1 ;
    for(i = 0 ; i < 2 ; i++) {
      int    m     = edge_vertices[n][i] ;
      icol  = m*neq + inc ;
      K(ilin,icol) = 0.5 ;
    }
  }
}



double*   FEM_ComputeBodyForceResidu(FEM_t* fem,IntFct_t* intfct,const double* f,const int dec)
/* Compute the residu due to a force */
{
  Element_t* el = FEM_GetElement(fem) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = IntFct_GetDimension(intfct) ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0. ;

  /* initialisation */
  for(i = 0 ; i < nn ; i++) r[i] = zero ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }

  /* 0D */
  if(dim == 0) {
    if(nn == 1) {
      double radius = x[0][0] ;
      
      r[0] = f[0] ;
      
      if(Symmetry_IsCylindrical(sym)) r[0] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[0] *= 4*M_PI*radius*radius ;
      
    } else {
      arret("FEM_ComputeBodyForceResidu : impossible") ;
    }
    
    return(r) ;
  }

  /* 1D, 2D, 3D */
  for(p = 0 ; p < np ; p++) {
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    
    /* cas axisymetrique ou shperique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      double rayon = zero ;
      
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    /* R(i) = F*H(i) */
    for(i = 0 ; i < nn ; i++) r[i] += a*f[p*dec]*h[i] ;
  }
  
  return(r) ;
}



double*   FEM_ComputeStrainWorkResidu(FEM_t* fem,IntFct_t* intfct,const double* sig,const int dec)
/* Compute the residu due to strain work */
{
#define DH(n,i)     (dh[(n)*dim + (i)])
#define SIG(i,j)    (sig[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*dim + (j)])
#define R(n,i)      (r[(n)*dim + (i)])
  Element_t* el = FEM_GetElement(fem) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = IntFct_GetDimension(intfct) ;
  int ndof = nn*dim ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = ndof*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0. ;
  
  /* initialization */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , sig += dec) {
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nn) ;
    double rayon = 0. ;
    int    j ;
    
    /* cas axisymetrique ou shperique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    /* R(i,j) = DH(i,k)*J(k,l)*S(l,j) */
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < dim ; j++) {
      int    k,l ;
      for(k = 0 ; k < dim ; k++) for(l = 0 ; l < dim ; l++) {
        R(i,j) += a*DH(i,k)*CAJ(k,l)*SIG(l,j) ;
      }
    }
    
    /* cas axisymetrique ou spherique */
    if(Symmetry_IsCylindrical(sym)) {
      /* R(i,0) = H(i)/r*S(theta,theta) */
      for(i = 0 ; i < nn ; i++) {
        R(i,0) += a*h[i]/rayon*SIG(2,2) ;
      }
    } else if(Symmetry_IsSpherical(sym)) {
      /* R(i,0) = H(i)/r*(SIG(theta,theta)+SIG(phi,phi)) */
      for(i = 0 ; i < nn ; i++) {
        R(i,0) += a*h[i]/rayon*(SIG(1,1) + SIG(2,2)) ;
      }
    }
  }
  
  return(r) ;
  
#undef DH
#undef SIG
#undef CAJ
#undef R
}



double*   FEM_ComputeFluxResidu(FEM_t* fem,IntFct_t* intfct,const double* f,const int dec)
/* Compute the residu due to flux */
{
#define DH(n,i)     (dh[(n)*dim+(i)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
  Element_t* el = FEM_GetElement(fem) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = IntFct_GetDimension(intfct) ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  int    i,p ;
  double zero = 0. ;
  
  /* initialisation */
  for(i = 0 ; i < nn ; i++) r[i] = zero ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }
  
  /* boucle sur les points d'integration */
  for(p = 0 ; p < np ; p++ , f +=dec) {
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
    double a   = weight[p]*d ;
    double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nn) ;
    double rayon = 0. ;
    
    /* cas axisymetrique ou shperique */
    if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
      for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
      a *= 2*M_PI*rayon ;
      if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
    }
    
    /* R(i) = DH(i,k)*J(k,j)*F(j) */
    for(i = 0 ; i < nn ; i++) {
      int    j,k ;
      
      for(j = 0 ; j < dim ; j++) for(k = 0 ; k < dim ; k++) {
        r[i] += a*DH(i,k)*CAJ(k,j)*f[j] ;
      }
    }
  }
  
  return(r) ;
  
#undef DH
#undef CAJ
}


double* FEM_ComputeSurfaceLoadResidu(FEM_t* fem,IntFct_t* intfct,Load_t* load,const double t,const double dt)
/* Compute the residu force due to surface loads (r) */
{
  Element_t* el = FEM_GetElement(fem) ;
  Geometry_t* geom = Element_GetGeometry(el) ;
  unsigned short int dim = Geometry_GetDimension(geom) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  Node_t* *no = Element_GetPointerToNode(el) ;
  Field_t* field = Load_GetField(load) ;
  char    *load_eqn = Load_GetNameOfEquation(load) ;
  char    *load_type = Load_GetType(load) ;
  Function_t* function = Load_GetFunction(load) ;
  int    nn  = IntFct_GetNbOfNodes(intfct) ;
  int    np  = IntFct_GetNbOfPoints(intfct) ;
  int    neq = Element_GetNbOfEquations(el) ;
  int    ieq = Element_FindEquationPositionIndex(el,load_eqn) ;
  int    ndof = nn*neq ;
  size_t SizeNeeded = ndof*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Node_GetCoordinate(no[i]) ;
  }

  /* initialization */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(field == NULL) return(r) ;

  /* Position index of the equation */
  if(ieq < 0) arret("FEM_ComputeSurfaceLoadResidu (1) : unknown equation") ;

  /* flux */
  if(strncmp(load_type,"flux",4) == 0) {
    double ft = dt ;
    
    if(function != NULL) {
      ft = Function_ComputeValue(function,t) - Function_ComputeValue(function,t - dt) ;
    }
    
    if(dim == 1 && nn == 1) {
      double radius = x[0][0] ;
      r[ieq] = ft*Field_ComputeValueAtPoint(field,x[0],dim) ;
      if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
      return(r) ;
    }
    
    if(dim >= 2) {
      double* rb ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      
      for(p = 0 ; p < np ; p++) {
        double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < nn ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*Field_ComputeValueAtPoint(field,y,dim) ;
      }
      
      rb = FEM_ComputeBodyForceResidu(fem,intfct,f,1) ;
      
      for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
      FEM_FreeBufferFrom(fem,rb) ;
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
      double radius = x[0][0] ;
      r[ieq] = ft*Field_ComputeValueAtPoint(field,x[0],dim) ;
      if(Symmetry_IsCylindrical(sym)) r[ieq] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[ieq] *= 4*M_PI*radius*radius ;
      return(r) ;
    }
    
    if(dim >= 2) {
      double* rb ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      
      for(p = 0 ; p < np ; p++) {
        double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < nn ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*Field_ComputeValueAtPoint(field,y,dim) ;
      }
      
      rb = FEM_ComputeBodyForceResidu(fem,intfct,f,1) ;
      
      for(i = 0 ; i < nn ; i++) r[i*neq+ieq] = rb[i] ;
      
      FEM_FreeBufferFrom(fem,rb) ;
      return(r) ;
    }
  }
  
  
  /* pressure */
  if(strncmp(load_type,"press",5) == 0) {
    double ft = 1. ;
    
    if(function != NULL) {
      ft = Function_ComputeValue(function,t) ;
    }
    
    if(dim >= 2) {
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      
      for(p = 0 ; p < np ; p++) {
        double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
        double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
        double* n = Element_ComputeNormalVector(el,dh,nn) ;
        double y[3] = {0.,0.,0.,} ;
        double f0 ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < nn ; j++) y[i] += h[j]*x[j][i] ;
        }
        f0 = ft*Field_ComputeValueAtPoint(field,y,dim) ;
        for(i = 0 ; i < dim ; i++) f[p*dim+i] = -f0*n[i] ;
      }
      
      for(i = 0 ; i < dim ; i++) {
        int j ;
        double* rb = FEM_ComputeBodyForceResidu(fem,intfct,f+i,dim) ;
        
        for(j = 0 ; j < nn ; j++) r[j*neq+ieq+i] = rb[j] ;
        
        FEM_FreeBufferFrom(fem,rb) ;
      }
      
      return(r) ;
    }
  }
  
  
  /* tensor component */
  if(strncmp(load_type,"sig_",4) == 0) {
    int ii = load_type[4] - '1' ;
    int jj = load_type[5] - '1' ;
    double ft = 1. ;
    
    if(ii < 0 || ii >= dim) {
      arret("FEM_ComputeSurfaceLoadResidu (2) : unknown type") ;
    }
    
    if(jj < 0 || jj >= dim) {
      arret("FEM_ComputeSurfaceLoadResidu (2) : unknown type") ;
    }
    
    if(function != NULL) {
      ft = Function_ComputeValue(function,t) ;
    }
    
    if(dim >= 2) {
      double* rb ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      int j ;
      
      for(p = 0 ; p < np ; p++) {
        double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
        double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
        double* n = Element_ComputeNormalVector(el,dh,nn) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          for(j = 0 ; j < nn ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*Field_ComputeValueAtPoint(field,y,dim)*n[jj] ;
      }
      
      rb = FEM_ComputeBodyForceResidu(fem,intfct,f,1) ;
      
      for(j = 0 ; j < nn ; j++) r[j*neq+ieq+ii] = rb[j] ;
      
      FEM_FreeBufferFrom(fem,rb) ;
      return(r) ;
    }
  }
  
  arret("FEM_ComputeSurfaceLoadResidu (2) : unknown load") ;
  return(r) ;
}




void FEM_TransformResiduFromDegree2IntoDegree1(FEM_t* fem,const int equ,double* r)
{
  Element_t* el = FEM_GetElement(fem) ;
  int    nn = Element_GetNbOfNodes(el) ;
  int    dim = Element_GetDimension(el) ;
  int    neq = Element_GetNbOfEquations(el) ;
  int    n_vertices = 0 ;
  int    edge_vertices_line[2]       = {0,1} ;
  int    edge_vertices_triangle[6]   = {0,1,1,2,2,0} ;
  int    edge_vertices_quadrangle[8] = {0,1,1,2,2,3,3,0} ;
  int    *edge_vertices[Element_MaxNbOfNodes] ;
  int    n ;

  if(nn > Element_MaxNbOfNodes) arret("FEM_TransformResiduFromDegree2IntoDegree1") ;

  if(dim == 1) {
    if(nn == 3) {
      n_vertices = 2 ;
      edge_vertices[n_vertices] = edge_vertices_line ;
    } else return ;
  } else if(dim == 2) {
    if(nn == 6) {
      n_vertices = 3 ;
      edge_vertices[n_vertices] = edge_vertices_triangle ;
    } else if(nn == 8) {
      n_vertices = 4 ;
      edge_vertices[n_vertices] = edge_vertices_quadrangle ;
    } else return ;
  } else if(dim == 3) {
    arret("FEM_TransformResiduFromDegree2IntoDegree1") ;
  } else {
    arret("FEM_TransformResiduFromDegree2IntoDegree1") ;
  }

  for(n = n_vertices + 1 ; n < nn ; n++) {
    edge_vertices[n] = edge_vertices[n - 1] + 2 ;
  }

  for(n = n_vertices ; n < nn ; n++) {
    int ilin  = n*neq + equ ;
    int i ;
    for(i = 0 ; i < 2 ; i++) {
      int    m     = edge_vertices[n][i] ;
      int    jlin  = m*neq + equ ;
      r[jlin] += 0.5*r[ilin] ;
    }
    r[ilin] = 0. ;
  }
}



double*  FEM_ComputeIsoShapeFctInActualSpace(FEM_t* fem,double* x)
/* fonction d'interpolation (h) et ses derivees (dh) en x */
{
  Element_t* el = FEM_GetElement(fem) ;
  
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



void   FEM_AverageStresses(Mesh_t* mesh,double* stress)
{
  unsigned int nel = Mesh_GetNbOfElements(mesh) ;
  Element_t* el0 = Mesh_GetElement(mesh) ;
  double vol = FEM_ComputeVolume(mesh) ;
  
  /* Stress integration */
  {
    int n_models = 3 ;
    const char* modelswithstresses[3] = {"Elast","Plast","MechaMic"} ;
    const int   stressindex[3] = {0,0,0} ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double sig = 0 ;
      unsigned int ie ;
    
      for(ie = 0 ; ie < nel ; ie++) {
        Element_t* el = el0 + ie ;
        double* vim = Element_GetCurrentImplicitTerm(el) ;
        IntFct_t* intfct = Element_GetIntFct(el) ;
        int nbofintpoints = IntFct_GetNbOfPoints(intfct) ;
        int ni = Element_GetNbOfImplicitTerms(el) ;
        int nvi = ni/nbofintpoints ;
        FEM_t*    fem    = FEM_GetInstance(el) ;
    
        if(Element_IsSubmanifold(el)) continue ;
        
        {
          Model_t* model = Element_GetModel(el) ;
          char* codename = Model_GetCodeNameOfModel(model) ;
          int j = 0 ;
  
          while(j < n_models && strcmp(modelswithstresses[j],codename)) j++ ;
  
          if(j == n_models) {
            arret("FEM_AverageStresses: model not implemented") ;
          }
          
          vim += stressindex[j] ;
        }
    
        sig +=  FEM_IntegrateOverElement(fem,intfct,vim + i,nvi) ;
      }
    
      /* Stress average */
      stress[i] = sig/vol ;
    }
  }
}



double   FEM_ComputeVolume(Mesh_t* mesh)
{
  unsigned int nel = Mesh_GetNbOfElements(mesh) ;
  Element_t* el0 = Mesh_GetElement(mesh) ;
  double vol = 0 ;
  
  /* The volume */
  {
    unsigned int ie ;
    
    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* el = el0 + ie ;
      IntFct_t* intfct = Element_GetIntFct(el) ;
      FEM_t*    fem    = FEM_GetInstance(el) ;
      double one = 1 ;
    
      if(Element_IsSubmanifold(el)) continue ;
      
      vol +=  FEM_IntegrateOverElement(fem,intfct,&one,0) ;
    }
  }
  
  return(vol) ;
}



/* FEM_Compute functions in the form (FEM_t*,double**,IntFct_t*,int,int) */

double FEM_ComputeUnknown(FEM_t* fem,double** u,IntFct_t* intfct,int p,int inc)
/** Compute the unknown located at "inc" at the interpolation point "p" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(element,n,inc)])
  Element_t* element = FEM_GetElement(fem) ;
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



double* FEM_ComputeDisplacementVector(FEM_t* fem,double** u,IntFct_t* intfct,int p,int inc)
/** Compute the displacement vector located at 
 *  "inc,inc+1,inc+2" at the interpolation point "p" */
{
  size_t SizeNeeded = 3*sizeof(double) ;
  double* dis = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;

  {
    Element_t* element = FEM_GetElement(fem) ;
    int dim = Element_GetDimensionOfSpace(element) ;
    int i ;
    
    for(i = 0 ; i < dim ; i++) {
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,inc+i) ;
    }
    
    for(i = dim ; i < 3 ; i++) dis[i] = 0 ;
  }
        
  return(dis) ;
}



double* FEM_ComputeUnknownGradient(FEM_t* fem,double** u,IntFct_t* intfct,int p,int inc)
/** Compute the unknown gradient located at "inc" at the interpolation point "p" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(element,n,inc)])
#define DH(n,i)  (dh[(n)*dim+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  int    dim  = Element_GetDimension(element) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  

  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    /* interpolation functions */
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
    double* gu = grad ;
    int i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k,l ;
        
      for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim ; l++) {
        gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
      }
    }
      
    Element_FreeBufferFrom(element,cj) ;
  }
  
  return(grad) ;
  
#undef U
#undef DH
#undef CJ
}



double* FEM_ComputeLinearStrainTensor(FEM_t* fem,double** u,IntFct_t* intfct,int p,int inc)
/** Compute the 3D linearized strain tensor for the displacement vector 
 *  located at "inc" and at the interpolation point "p" */
{
#define U(n,i)   (u[n][Element_GetNodalUnknownPosition(element,n,inc + (i))])
#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  Symmetry_t sym = Element_GetSymmetry(element) ;
  int    dim  = Element_GetDimension(element) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* strain = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    
    /* interpolation functions */
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
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
    for(i = 0 ; i < dim ; i++) {
      int j ;
        
      for(j = 0 ; j < dim ; j++) {
        int k ;
          
        for(k = 0 ; k < nn ; k++) {
          int l ;
          
          for(l = 0 ; l < dim ; l++) {
            gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
          }
        }
      }
    }
      
    {
      double* eps = strain ;
  
      /* Linearized strain tensor */
      for(i = 0 ; i < 3 ; i++) {
        int j ;
        
        for(j = 0 ; j < 3 ; j++) {
          EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
        }
      }
  
      /* symmetric cases : axisymmetrical or spherical */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double rayon = 0. ;
        double u_r   = 0. ;
        
        for(i = 0 ; i < nn ; i++) {
          double* x = Element_GetNodeCoordinate(element,i) ;
          
          rayon += h[i]*x[0] ;
          u_r   += h[i]*U(i,0) ;
        }
        EPS(2,2) += u_r/rayon ;
        if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
      }
    }
      
    Element_FreeBufferFrom(element,cj) ;
  }
  
  return(strain) ;
#undef U
#undef DH
#undef EPS
#undef CJ
}

/* End */



/* FEM_Compute functions in the form (FEM_t*,double*,int,int) */

double FEM_ComputeCurrentUnknown(FEM_t* fem,double* h,int nn,int inc)
/** Compute the unknown at the current time located at "inc" */
{
#define U(n)   (Element_GetValueOfCurrentNodalUnknown(element,n,inc))
  Element_t* element = FEM_GetElement(fem) ;
  int    i ;
  double par = 0. ;
  
  for(i = 0 ; i < nn ; i++) {
    par += h[i]*U(i) ;
  }
  
  return (par) ;
#undef U
}


double FEM_ComputePreviousUnknown(FEM_t* fem,double* h,int nn,int inc)
/** Compute the unknown at the previous time located at "inc" */
{
#define U(n)   (Element_GetValueOfPreviousNodalUnknown(element,n,inc))
  Element_t* element = FEM_GetElement(fem) ;
  int    i ;
  double par = 0. ;
  
  for(i = 0 ; i < nn ; i++) {
    par += h[i]*U(i) ;
  }
  
  return (par) ;
#undef U
}



double* FEM_ComputeCurrentUnknownGradient(FEM_t* fem,double* dh,int nn,int inc)
/** Compute the current unknown gradient for the unknown located at "inc" */
{
#define U(n)   (Element_GetValueOfCurrentNodalUnknown(element,n,inc))
#define DH(n,i)  (dh[(n)*dim+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  int    dim  = Element_GetDimension(element) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
    double* gu = grad ;
    int    i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k ;
      
      for(k = 0 ; k < nn ; k++) {
        int l ;
        
        for(l = 0 ; l < dim ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
    }
  
    Element_FreeBufferFrom(element,cj) ;
  }
  
  return(grad) ;
  
#undef U
#undef DH
#undef CJ
}



double* FEM_ComputePreviousUnknownGradient(FEM_t* fem,double* dh,int nn,int inc)
/** Compute the previous unknown gradient for the unknown located at "inc" */
{
#define U(n)   (Element_GetValueOfPreviousNodalUnknown(element,n,inc))
#define DH(n,i)  (dh[(n)*dim+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  int    dim  = Element_GetDimension(element) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
    double* gu = grad ;
    int    i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k ;
      
      for(k = 0 ; k < nn ; k++) {
        int l ;
        
        for(l = 0 ; l < dim ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
    }
  
    Element_FreeBufferFrom(element,cj) ;
  }

  return(grad) ;
  
#undef U
#undef DH
#undef CJ
}
/* End */



/* FEM_Compute functions in the form (FEM_t*,double*,double*,int,int) */

double* FEM_ComputeCurrentLinearStrainTensor(FEM_t* fem,double* h,double* dh,int nn,int inc)
/** Compute the 3D linearized strain tensor for a displacement vector located at "inc" */
{
#define U(n,i)   (Element_GetValueOfCurrentNodalUnknown(element,n,inc + (i)))
#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(element) ;
  int    dim  = Element_GetDimension(element) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj) : Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim ; l++) {
      gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases : axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(element,i) ;
      rayon += h[i]*x[0] ;
      u_r   += h[i]*U(i,0) ;
    }
    EPS(2,2) += u_r/rayon ;
    if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
  }
  
  return(eps) ;
#undef U
#undef DH
#undef EPS
#undef CJ
}



double* FEM_ComputeIncrementalLinearStrainTensor(FEM_t* fem,double* h,double* dh,int nn,int inc)
/** Compute the 3D incremental linearized strain tensor for a displacement vector located at "inc" */
{
#define U(n,i)   (Element_GetValueOfCurrentNodalUnknown(element,n,inc + (i)))
#define U_n(n,i) (Element_GetValueOfPreviousNodalUnknown(element,n,inc + (i)))
#define DU(n,i)  (U(n,i) - U_n(n,i))
#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(element) ;
  int    dim  = Element_GetDimension(element) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj) : Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim ; l++) {
      gu[i][j] += DU(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases : axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(element,i) ;
      rayon += h[i]*x[0] ;
      u_r   += h[i]*DU(i,0) ;
    }
    EPS(2,2) += u_r/rayon ;
    if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
  }
  
  return(eps) ;
#undef U
#undef U_n
#undef DU
#undef DH
#undef EPS
#undef CJ
}



double* FEM_ComputePreviousLinearStrainTensor(FEM_t* fem,double* h,double* dh,int nn,int inc)
/** Compute the 3D linearized strain tensor for a displacement vector located at "inc" */
{
#define U(n,i)   (Element_GetValueOfPreviousNodalUnknown(element,n,inc + (i)))
#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(element) ;
  int    dim  = Element_GetDimension(element) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj) : Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim ; l++) {
      gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases : axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(element,i) ;
      rayon += h[i]*x[0] ;
      u_r   += h[i]*U(i,0) ;
    }
    EPS(2,2) += u_r/rayon ;
    if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
  }
  
  return(eps) ;
#undef U
#undef DH
#undef EPS
#undef CJ
}
/* End */



double   FEM_IntegrateOverElement(FEM_t* fem,IntFct_t* intfct,double* f,int shift)
/** Integrate f over the element through the integration points. 
 *  Return the integration result. */
{
  Element_t* el = FEM_GetElement(fem) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = IntFct_GetDimension(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double* x[Element_MaxNbOfNodes] ;
  double sum = 0 ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    x[i] = Element_GetNodeCoordinate(el,i) ;
  }

  /* 0D */
  if(dim == 0) {
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
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nn) ;
      double a   = weight[p]*d ;
    
      /* Axisymmetrical or spherical case */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
        double rayon = 0 ;
      
        for(i = 0 ; i < nn ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      sum += a*f[p*shift] ;
    }
  }
  
  return(sum) ;
}



/*
 * Intern functions
 */



#ifdef NOTDEFINED
double*  FEM_ComputeNormalVector(Element_t* element,double* dh,int nn)
/* Normale unitaire a un sous-espace de dimension dim-1 */
{
#define DH(n,i) (dh[(n)*dim_h+(i)])
  size_t SizeNeeded = 3*sizeof(double) ;
  double* norm = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    dim_h = Element_GetDimension(element) ;
  int    dim   = dim_h + 1 ;
  int    i,j ;
  double c[3][3] ;

  if(dim > 3) arret("FEM_ComputeNormalVector") ;
  
  /* le jacobien */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
    int    k ;
    c[i][j] = 0. ;
    for(k = 0 ; k < nn ; k++) {
      double* x = Element_GetNodeCoordinate(element,k) ;
      c[i][j] += x[i]*DH(k,j) ;
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

  arret("FEM_ComputeNormalVector") ;
  return(norm) ;
#undef DH
}
#endif


double* FEM_ComputeInverseJacobianMatrix(Element_t* element,double* dh,int nn)
/** Compute the inverse jacobian matrix */
{
#define DH(n,i)  (dh[(n)*dim_h+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  size_t SizeNeeded = 9*sizeof(double) ;
  double* cj = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    dim_h = Element_GetDimension(element) ;
  int    dim = Geometry_GetDimension(Element_GetGeometry(element));
  int    i,j ;
  double jc[3][3],dt,td ;
  double* jc1 = FEM_ComputeJacobianMatrix(element,dh,nn) ;
  
  /* The jacobian matrix */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
    jc[i][j] = jc1[i*dim_h + j] ;
  }
  
  Element_FreeBufferFrom(element,jc1) ;

  /* 1. Volume */
  if(dim_h == 3) {
    if(dim == 3) {
      dt  = jc[0][0]*jc[1][1]*jc[2][2] + jc[1][0]*jc[2][1]*jc[0][2]
          + jc[2][0]*jc[0][1]*jc[1][2] - jc[2][0]*jc[1][1]*jc[0][2]
          - jc[1][0]*jc[0][1]*jc[2][2] - jc[0][0]*jc[2][1]*jc[1][2] ;
      if(dt != 0.) {
        td = 1./dt ;
        CJ(0,0) = (jc[1][1]*jc[2][2] - jc[2][1]*jc[1][2])*td ;
        CJ(0,1) = (jc[2][1]*jc[0][2] - jc[0][1]*jc[2][2])*td ;
        CJ(0,2) = (jc[0][1]*jc[1][2] - jc[1][1]*jc[0][2])*td ;
        CJ(1,0) = (jc[1][2]*jc[2][0] - jc[2][2]*jc[1][0])*td ;
        CJ(1,1) = (jc[2][2]*jc[0][0] - jc[0][2]*jc[2][0])*td ;
        CJ(1,2) = (jc[0][2]*jc[1][0] - jc[1][2]*jc[0][0])*td ;
        CJ(2,0) = (jc[1][0]*jc[2][1] - jc[2][0]*jc[1][1])*td ;
        CJ(2,1) = (jc[2][0]*jc[0][1] - jc[0][0]*jc[2][1])*td ;
        CJ(2,2) = (jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1])*td ;
      } else arret("FEM_ComputeInverseJacobianMatrix (1)") ;
      return(cj) ;
    } else arret("FEM_ComputeInverseJacobianMatrix (2)") ;

  /* 2. Surface */
  } else if(dim_h == 2) {
    if(dim == 3) {
      double v ;
      /* on calcule la normale a la surface */
      jc[0][2] = jc[1][0]*jc[2][1] - jc[2][0]*jc[1][1] ;
      jc[1][2] = jc[2][0]*jc[0][1] - jc[0][0]*jc[2][1] ;
      jc[2][2] = jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1] ;
      v = sqrt(jc[0][2]*jc[0][2] + jc[1][2]*jc[1][2] + jc[2][2]*jc[2][2]) ;
      jc[0][2] /= v ;
      jc[1][2] /= v ;
      jc[2][2] /= v ;
      /* puis on inverse */
      dt  = jc[0][0]*jc[1][1]*jc[2][2] + jc[1][0]*jc[2][1]*jc[0][2]
          + jc[2][0]*jc[0][1]*jc[1][2] - jc[2][0]*jc[1][1]*jc[0][2]
          - jc[1][0]*jc[0][1]*jc[2][2] - jc[0][0]*jc[2][1]*jc[1][2] ;
      if(dt != 0.) {
        td = 1./dt ;
        CJ(0,0) = (jc[1][1]*jc[2][2] - jc[2][1]*jc[1][2])*td ;
        CJ(0,1) = (jc[2][1]*jc[0][2] - jc[0][1]*jc[2][2])*td ;
        CJ(0,2) = (jc[0][1]*jc[1][2] - jc[1][1]*jc[0][2])*td ;
        CJ(1,0) = (jc[1][2]*jc[2][0] - jc[2][2]*jc[1][0])*td ;
        CJ(1,1) = (jc[2][2]*jc[0][0] - jc[0][2]*jc[2][0])*td ;
        CJ(1,2) = (jc[0][2]*jc[1][0] - jc[1][2]*jc[0][0])*td ;
        CJ(2,0) = (jc[1][0]*jc[2][1] - jc[2][0]*jc[1][1])*td ;
        CJ(2,1) = (jc[2][0]*jc[0][1] - jc[0][0]*jc[2][1])*td ;
        CJ(2,2) = (jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1])*td ;
      } else arret("FEM_ComputeInverseJacobianMatrix (3)") ;
      return(cj) ;
    } else if(dim == 2) {
      dt = jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1] ;
      if(dt != 0.) {
        td = 1./dt ;
        CJ(0,0) =  jc[1][1]*td ;
        CJ(0,1) = -jc[0][1]*td ;
        CJ(1,0) = -jc[1][0]*td ;
        CJ(1,1) =  jc[0][0]*td ;
      } else arret("FEM_ComputeInverseJacobianMatrix (4)") ;
      return(cj) ;
    } else arret("FEM_ComputeInverseJacobianMatrix (5)") ;

  /* 3. Ligne */
  } else if(dim_h == 1) {
    if(dim == 3) {
      arret("FEM_ComputeInverseJacobianMatrix (6)") ;
    } else if(dim == 2) {
      /* on calcul la normale, c1^e_z, a la ligne */
      double v = sqrt(jc[0][0]*jc[0][0] + jc[1][0]*jc[1][0]) ;
      jc[0][1] =   jc[1][0]/v ;
      jc[1][1] = - jc[0][0]/v ;
      /* puis on inverse */
      dt = jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1] ;
      if(dt != 0.) {
        td = 1./dt ;
        CJ(0,0) =  jc[1][1]*td ;
        CJ(0,1) = -jc[0][1]*td ;
        CJ(1,0) = -jc[1][0]*td ;
        CJ(1,1) =  jc[0][0]*td ;
      } else arret("FEM_ComputeInverseJacobianMatrix (7)") ;
      return(cj) ;
    } else if(dim == 1) {
      dt = jc[0][0] ;
      if(dt != 0.) CJ(0,0) = 1./dt ;
      else arret("FEM_ComputeInverseJacobianMatrix (8)") ;
      return(cj) ;
    } else arret("FEM_ComputeInverseJacobianMatrix (9)") ;

  /* 4. Point */
  } else if(dim_h == 0) {
    return(cj) ;

  }

  arret("FEM_ComputeInverseJacobianMatrix (10)") ;
  return(NULL) ;
#undef CJ
#undef DH
}



double FEM_ComputeJacobianDeterminant(Element_t* element,double* dh,int nn)
/** Compute the determinant of the jacobian matrix */
{
#define DH(n,i)  (dh[(n)*dim_h+(i)])
#define J(i,j)   (jac[(i)*dim_h+(j)])
  int    dim_h = Element_GetDimension(element) ;
  int    dim   = Geometry_GetDimension(Element_GetGeometry(element));
  double* jac  = FEM_ComputeJacobianMatrix(element,dh,nn) ;
  double det ;
  
  /* 
     1. For a volume  : it's det(J1,J2,J3)
     2. For a surface : it's det(J1,J2,e3) 
        where e3 is the unit normal vector to (J1,J2)
     3. For a line    : it's det(J1,e2,e3) = |J1|
        where (e2,e3) are orthonormed vectors in the normal plane to J1
     4. For a point   : it's det(e1,e2,e3) = 1
        where (e1,e2,e3) is an orthonormed system.
  */

  /* 1. Volume : det(J1,J2,J3) */
  if(dim_h == 3) {
    if(dim == 3) {
      det  = J(0,0)*J(1,1)*J(2,2) - J(0,0)*J(2,1)*J(1,2)
           + J(1,0)*J(2,1)*J(0,2) - J(1,0)*J(0,1)*J(2,2)
           + J(2,0)*J(0,1)*J(1,2) - J(2,0)*J(1,1)*J(0,2) ;
      det  = fabs(det) ;
      
    } else {
      arret("FEM_ComputeJacobianDeterminant (1)") ;
    }

  /* 2. Surface : det(J1,J2,e3) = |c1^c2| */
  } else if(dim_h == 2) {
    if(dim == 3) {
      double v[3] ;
      
      v[0] = J(1,0)*J(2,1) - J(2,0)*J(1,1) ;
      v[1] = J(2,0)*J(0,1) - J(0,0)*J(2,1) ;
      v[2] = J(0,0)*J(1,1) - J(1,0)*J(0,1) ;
      
      det = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ;
      
    } else if(dim == 2) {
      det  = J(0,0)*J(1,1) - J(1,0)*J(0,1) ;
      det  = fabs(det) ;

    } else {
      arret("FEM_ComputeJacobianDeterminant (2)") ;
    }

  /* 3. Line : det(J1,e2,e3) = |J1|(e1,e2,e3) = |J1|   */
  } else if(dim_h == 1) {
    int i ;
    
    det = 0. ;
    for(i = 0 ; i < dim ; i++) {
      det += J(i,0)*J(i,0) ;
    }
    
    det = sqrt(det) ;

  /* 4. Point : det(e1,e2,e3) = 1 */
  } else if(dim_h == 0) {
    det = 1. ;
    
  } else {
  
    arret("FEM_ComputeJacobianDeterminant (3)") ;
  }
  
  Element_FreeBufferFrom(element,jac) ;
  
  return(det) ;
#undef DH
#undef J
}


double* FEM_ComputeJacobianMatrix(Element_t* element,double* dh,int nn)
/** Compute the jacobian matrix */
{
#define DH(n,i)  (dh[(n)*dim_h+(i)])
#define JC(i,j)  (jc[(i)*dim_h+(j)])
  size_t SizeNeeded = 9*sizeof(double) ;
  double* jc = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    dim_h = Element_GetDimension(element) ;
  int    dim = Geometry_GetDimension(Element_GetGeometry(element));
  int    i,j ;
  
  /* The Jacobian matrix */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
    int k ;
    JC(i,j) = 0. ;
    for(k = 0 ; k < nn ; k++) {
      double* x = Element_GetNodeCoordinate(element,k) ;
      JC(i,j) += x[i]*DH(k,j) ;
    }
  }
  
  return(jc) ;
#undef DH
#undef JC
}




/* FEM_Compute functions in the form (FEM_t*,IntFct_t*,double**,int) */

/* NOT USED */

#ifdef NOTDEFINED
double* FEM_ComputeUnknowns(FEM_t* fem,IntFct_t* intfct,double** u,int inc)
/** Compute the unknowns at the interpolation points located at "inc" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(element,n,inc)])
  Element_t* element = FEM_GetElement(fem) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  size_t SizeNeeded = np*sizeof(double) ;
  double* unk = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    int p ;
    
    for(p = 0 ; p < np ; p++) {
      /* interpolation functions */
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    
      /* Unknown */
      {
        double par = 0. ;
        int i ;
        
        for(i = 0 ; i < nn ; i++) {
          par += h[i]*U(i) ;
        }
        
        unk[p] = par ;
      }
    }
  }
  
  return(unk) ;
#undef U
}



double* FEM_ComputeUnknownGradients(FEM_t* fem,IntFct_t* intfct,double** u,int inc)
/** Compute the unknown gradients at the interpolation points located at "inc" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(element,n,inc)])
#define DH(n,i)  (dh[(n)*dim+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  int    dim  = Element_GetDimension(element) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  size_t SizeNeeded = np*3*sizeof(double) ;
  double* ugrd = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  

  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    int p ;
    
    for(p = 0 ; p < np ; p++) {
      /* interpolation functions */
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      /* inverse jacobian matrix (cj) */
      double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
      double gu[3] ;
      int i ;
  
      /* initialisation of gu */
      for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
      /* the parameter gradient (gu) */
      for(i = 0 ; i < dim ; i++) {
        int k,l ;
        
        for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
      
      {
        double* grad = ugrd + 3*p ;
        
        for(i = 0 ; i < 3 ; i++) grad[i] = gu[i] ;
      }
      
      Element_FreeBufferFrom(element,cj) ;
    }
  }
  
  return(ugrd) ;
  
#undef U
#undef DH
#undef CJ
}



double* FEM_ComputeLinearStrainTensors(FEM_t* fem,IntFct_t* intfct,double** u,int inc)
/** Compute the 3D linearized strain tensors for the displacement vectors 
 *  located at "inc" and at the interpolation points */
{
#define U(n,i)   (u[n][Element_GetNodalUnknownPosition(element,n,inc + (i))])
#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  Element_t* element = FEM_GetElement(fem) ;
  Symmetry_t sym = Element_GetSymmetry(element) ;
  int    dim  = Element_GetDimension(element) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  size_t SizeNeeded = np*9*sizeof(double) ;
  double* strains = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    int p ;
    
    for(p = 0 ; p < np ; p++) {
      /* interpolation functions */
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      /* inverse jacobian matrix (cj) */
      double* cj = FEM_ComputeInverseJacobianMatrix(element,dh,nn) ;
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
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          int k ;
          
          for(k = 0 ; k < nn ; k++) {
            int l ;
          
            for(l = 0 ; l < dim ; l++) {
              gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
            }
          }
        }
      }
      
      {
        double* eps = strains + 9*p ;
  
        /* Linearized strain tensor */
        for(i = 0 ; i < 3 ; i++) {
          int j ;
        
          for(j = 0 ; j < 3 ; j++) {
            EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
          }
        }
  
        /* symmetric cases : axisymmetrical or spherical */
        if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
          double rayon = 0. ;
          double u_r   = 0. ;
        
          for(i = 0 ; i < nn ; i++) {
            double* x = Element_GetNodeCoordinate(element,i) ;
          
            rayon += h[i]*x[0] ;
            u_r   += h[i]*U(i,0) ;
          }
          EPS(2,2) += u_r/rayon ;
          if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
        }
      }
      
      Element_FreeBufferFrom(element,cj) ;
    }
  }
  
  return(strains) ;
#undef U
#undef DH
#undef EPS
#undef CJ
}
#endif
/* End */
