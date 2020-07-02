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
#include "Mry.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <assert.h>


//static FEM_t* instancefem = NULL ;

static FEM_t*  FEM_Create(void) ;
//static double* FEM_ComputeJacobianMatrixOld(Element_t*,double*,int) ;
//static double* FEM_ComputeInverseJacobianMatrixOld(Element_t*,double*,int) ;
//static double  FEM_ComputeJacobianDeterminantOld(Element_t*,double*,int) ;
static double* FEM_ComputeJacobianMatrixNew(Element_t*,double*,int,const int) ;
static double* FEM_ComputeInverseJacobianMatrixNew(Element_t*,double*,int,const int) ;
static double  FEM_ComputeJacobianDeterminantNew(Element_t*,double*,int,const int) ;
//static double* FEM_ComputeNormalVector(Element_t*,double*,int) ;
static void    FEM_CheckNumberingOfOverlappingNodes(Element_t*,const int) ;

#define FEM_ComputeJacobianMatrix        FEM_ComputeJacobianMatrixNew
#define FEM_ComputeJacobianDeterminant   FEM_ComputeJacobianDeterminantNew
#define FEM_ComputeInverseJacobianMatrix FEM_ComputeInverseJacobianMatrixNew

/* 
   Extern Functions 
*/
FEM_t* FEM_Create(void)
{
  FEM_t*  fem    = (FEM_t*) Mry_New(FEM_t) ;
  
  
  /* Space allocation for output */
  {
    int nn = Element_MaxNbOfNodes ;
    int neq = Model_MaxNbOfEquations ;
    int ndof = nn*neq ;
    double* output = (double*) Mry_New(double[ndof*ndof]) ;
    
    FEM_GetOutput(fem) = output ;
  }
  
  
  /* Space allocation for input */
  {
    int np = IntFct_MaxNbOfIntPoints ;
    double* input = (double*) Mry_New(double[np*FEM_MaxShift]) ;
    
    FEM_GetInput(fem) = input ;
  }
  
  
  /* Space allocation for pintfct */
  {
    IntFct_t** pintfct = (IntFct_t**) Mry_New(IntFct_t*[IntFcts_MaxNbOfIntFcts]) ;
    
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
  //*pfem = NULL ;
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
#define DH(n,i)     (dh[(n)*3 + (i)])
//#define DH(n,i)     (dh[(n)*dim + (i)])
#define C(i,j,k,l)  (c[(((i)*3 + (j))*3 + (k))*3 + (l)])
#define CAJ(i,j)    (caj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int nf  = IntFct_GetNbOfFunctions(fi) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int dim_h = IntFct_GetDimension(fi) ;
  int ndof = nn*dim ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* kr = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < ndof*ndof ; i++) kr[i] = zero ;
  }
  
  if(Element_IsSubmanifold(el)) {
    arret("FEM_ComputeElasticMatrix") ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  


  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;
      
        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;

        for(p = 0 ; p < np ; p++ , c += dec) {
          double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          double* norm = Element_ComputeNormalVector(el,dh,nf,dim_h) ;
          double r[9] = {0,0,0,0,0,0,0,0,0} ;
          int i,j,k,l ;

          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            arret("FEM_ComputeElasticMatrix: To be improved") ;
          }
      
          /* */
          #define R(i,k)   r[(i)*3 + (k)]
          #define NORM(i)  (norm[i])
        
          for(i = 0 ; i < 3 ; i++) for(k = 0 ; k < 3 ; k++) {
          
            R(i,k) = 0 ;

            for(j = 0 ; j < 3 ; j++) for(l = 0 ; l < 3 ; l++) {
              double cijkl = C(i,j,k,l) + C(j,i,k,l) + C(i,j,l,k) + C(j,i,l,k) ;
            
              R(i,k) += NORM(j) * cijkl * NORM(l) ;
            }
          
            R(i,k) *= 0.25 ;
          }
        
          #undef NORM

          /* KR(j,i,l,k) = H(j) * R(i,k) * H(l) */
          for(j = 0 ; j < nf ; j++) for(i = 0 ; i < dim ; i++) {
            int jj = Element_OverlappingNode(el,j) ;
          
            for(l = 0 ; l < nf ; l++) for(k = 0 ; k < dim ; k++) {
              double k_jilk = a * h[j] * R(i,k) * h[l] ;
              int ll = Element_OverlappingNode(el,l) ;
          
              KR(j*dim  + i,l*dim  + k) +=   k_jilk ;
              KR(jj*dim + i,l*dim  + k) += - k_jilk ;
              KR(j*dim  + i,ll*dim + k) += - k_jilk ;
              KR(jj*dim + i,ll*dim + k) +=   k_jilk ;
            }
          }

          #undef R
        }
      
        return(kr) ;
      }
    }
  }

  
  /* Loop on integration points */
  {
    int    p ;
    
    for(p = 0 ; p < np ; p++ , c += dec) {
      double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
      double jcj[3][3][3][3],jc[3][3],cj[3][3] ;
      double radius = zero ;
      int    i,j,k,l,r,s ;
    
      /* The radius in axisymmetrical or spherical case */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        
        for(i = 0 ; i < nf ; i++) radius += h[i]*x[i][0] ;
        
        a *= 2*M_PI*radius ;
        
        if(Symmetry_IsSpherical(sym)) a *= 2*radius ;
      }
    
      /* JCJ(r,i,k,s) = J(r,j) * C(i,j,k,l) * J(s,l) */
      for(r = 0 ; r < dim_h ; r++) for(i = 0 ; i < dim ; i++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim_h ; s++) {
        jcj[r][i][k][s] = zero ;
        for(j = 0 ; j < dim ; j++) for(l = 0 ; l < dim ; l++) {
          jcj[r][i][k][s] += CAJ(r,j) * C(i,j,k,l) * CAJ(s,l) ;
        }
      }
    
      /* KR(j,i,l,k) = DH(j,r) * JCJ(r,i,k,s) * DH(l,s) */
      for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(k = 0 ; k < dim ; k++) {
        for(r = 0 ; r < dim_h ; r++) for(s = 0 ; s < dim_h ; s++) {
          KR(j*dim+i,l*dim+k) += a * DH(j,r) * jcj[r][i][k][s] * DH(l,s) ;
        }
      }
    
      /* Axisymmmetrical case: 3 terms */
      if(Symmetry_IsCylindrical(sym)) {
        /* 1.a JC(r,i) = J(r,j) * C(i,j,theta,theta) */
        for(r = 0 ; r < dim_h ; r++) for(i = 0 ; i < dim ; i++) {
          jc[r][i] = zero ;
          for(j = 0 ; j < dim ; j++) jc[r][i] += CAJ(r,j) * C(i,j,2,2) ;
        }
        
        /* 1.b KR(j,i,l,0) = DH(j,r) * JC(r,i) * H(l)/r */
        for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(r = 0 ; r < dim_h ; r++) {
          KR(j*dim+i,l*dim) += a * DH(j,r) * jc[r][i] * h[l]/radius ;
        }
        
        /* 2.a CJ(k,s) = C(theta,theta,k,l) * J(s,l) */
        for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim_h ; s++) {
          cj[k][s] = zero ; 
          for(l = 0 ; l < dim ; l++) cj[k][s] += C(2,2,k,l) * CAJ(s,l) ;
        }
        
        /* 2.b KR(j,0,l,k) = H(j)/r * CJ(k,s) * DH(l,s) */
        for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn;l++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim_h ; s++) {
          KR(j*dim,l*dim+k) += a * h[j]/radius * cj[k][s] * DH(l,s) ;
        }
        
        /* 3.  KR(j,0,l,0) = H(j)/r * C(theta,theta,theta,theta) * H(l)/r */
        for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
          KR(j*dim,l*dim) += a * h[j]/radius * C(2,2,2,2) * h[l]/radius ;
        }
      
      /* Spherical case: 3 terms */
      } else if(Symmetry_IsSpherical(sym)) {
        
        /* 1.a JC(r,i) = J(r,j) * (C(i,j,theta,theta) + C(i,j,phi,phi)) */
        for(r = 0 ; r < dim_h ; r++) for(i = 0 ; i < dim ; i++) {
          jc[r][i] = zero ;
          for(j = 0 ; j < dim ; j++) jc[r][i] += CAJ(r,j) * (C(i,j,1,1) + C(i,j,2,2)) ;
        }
        
        /* 1.b KR(j,i,l,0) = DH(j,r) * JC(r,i) * H(l)/r */
        for(j = 0 ; j < nn ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nn ; l++) for(r = 0 ; r < dim_h ; r++) {
          KR(j*dim+i,l*dim) += a * DH(j,r) * jc[r][i] * h[l]/radius ;
        }
        
        /* 2.a CJ(k,s) = (C(phi,phi,k,l) + C(theta,theta,k,l)) * J(s,l) */
        for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim_h ; s++) {
          cj[k][s] = zero ; 
          for(l = 0 ; l < dim ; l++) cj[k][s] += (C(1,1,k,l) + C(2,2,k,l)) * CAJ(s,l) ;
        }
        
        /* 2.b KR(j,0,l,k) = H(j)/r * CJ(k,s) * DH(l,s) */
        for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) for(k = 0 ; k < dim ; k++) for(s = 0 ; s < dim_h ; s++) {
          KR(j*dim,l*dim+k) += a * h[j]/radius * cj[k][s] * DH(l,s) ;
        }
        
        /* 3.  KR(j,0,l,0) = H(j)/r * (C(phi,phi,phi,phi) + C(theta,theta,phi,phi) + C(phi,phi,theta,theta) + C(theta,theta,theta,theta)) * H(l)/r */
        {
          double cc = C(1,1,1,1) + C(1,1,2,2) + C(2,2,1,1) + C(2,2,2,2) ;
          
          for(j = 0 ; j < nn ; j++) for(l = 0 ; l < nn ; l++) {
            KR(j*dim,l*dim) += a * h[j]/radius * cc * h[l]/radius ;
          }
        }
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
#define DH(n,i)     (dh[(n)*3 + (i)])
//#define DH(n,i)     (dh[(n)*dim + (i)])
#define C(i,j)      (c[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int dim_h = IntFct_GetDimension(fi) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nf  = IntFct_GetNbOfFunctions(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int nrow = nn*dim ;
  int ncol = nn ;
  size_t SizeNeeded = nrow*ncol*(sizeof(double)) ;
  double* kc = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* Initialisation */
  {
    int i ;
    
    for(i = 0 ; i < nrow*ncol ; i++) kc[i] = zero ;
  }
  
  if(Element_IsSubmanifold(el)) {
    arret("FEM_ComputeBiotMatrix") ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  

  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;
      
        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;

        for(p = 0 ; p < np ; p++ , c += dec) {
          double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          double* norm = Element_ComputeNormalVector(el,dh,nf,dim_h) ;
          double r[3] = {0,0,0} ;
          int i,j,l ;

          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            arret("FEM_ComputeBiotMatrix: To be improved") ;
          }
      
          /* 
           * Pay attention to the orientation of the normal:
           * the normal N is oriented from the opposite face
           * to the main face of the element, i.e. from ii to i */
          #define R(i)     (r[i])
          #define NORM(i)  (norm[i])
        
          for(i = 0 ; i < 3 ; i++) {
            R(i) = 0 ;
            for(j = 0 ; j < 3 ; j++) {
              R(i) += C(i,j) * NORM(j) ;
            }
          }
        
          #undef NORM

          /* KR(j,i,l) = H(j) * R(i) * H(l) */
          for(j = 0 ; j < nf ; j++) {
            int jj = Element_OverlappingNode(el,j) ;
            
            for(i = 0 ; i < dim ; i++) {
              for(l = 0 ; l < nf ; l++) {
                double k_jil = a * h[j] * R(i) * h[l] ;
                int ll = Element_OverlappingNode(el,l) ;
                
                KC(j*dim  + i,l ) +=   0.5 * k_jil ;
                KC(jj*dim + i,l ) += - 0.5 * k_jil ;
                KC(j*dim  + i,ll) +=   0.5 * k_jil ;
                KC(jj*dim + i,ll) += - 0.5 * k_jil ;
              }
            }
          }

          #undef R
        }
      
        return(kc) ;
      }
    }
  }
  
  /* boucle sur les points d'integration */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++ , c += dec) {
      double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
      double jc[3][3] ;
      double rayon ;
      int    i,j,k,l ;
    
      /* cas axisymetrique ou spherique */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        rayon = zero ;
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      /* JC(k,i) = J(k,j)*C(j,i) */
      for(k = 0 ; k < dim_h ; k++) for(i = 0 ; i < dim ; i++) {
        jc[k][i] = zero ;
        for(j = 0 ; j < dim ; j++) jc[k][i] += CAJ(k,j)*C(j,i) ;
      }
      /* KC(j,i,l) = DH(j,k)*JC(k,i)*H(l) */
      for(j = 0 ; j < nf ; j++) for(i = 0 ; i < dim ; i++) for(l = 0 ; l < nf ; l++) for(k = 0 ; k < dim_h ; k++) {
        KC(j*dim+i,l) += a*DH(j,k)*jc[k][i]*h[l] ;
      }
      /* cas axisymetrique: (r,z,theta) */
      if(Symmetry_IsCylindrical(sym)) {
        /* KC(j,0,l) = H(j)/r*C(theta,theta)*H(l) */
        for(j = 0 ; j < nf ; j++) for(l = 0 ; l < nf ; l++) {
          KC(j*dim,l) += a*h[j]/rayon*C(2,2)*h[l] ;
        }
      /* cas spherique: (r,theta,phi) */
      } else if(Symmetry_IsSpherical(sym)) {
        /* KC(j,0,l) = H(j)/r*(C(theta,theta)+C(phi,phi))*H(l) */
        for(j = 0 ; j < nf ; j++) for(l = 0 ; l < nf ; l++) {
          KC(j*dim,l) += a*h[j]/rayon*(C(1,1)+C(2,2))*h[l] ;
        }
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
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nf  = IntFct_GetNbOfFunctions(fi) ;
  int dim_h = IntFct_GetDimension(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* km = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn*nn ; i++) {
      km[i] = 0 ;
    }
  }
  
  //if(Element_IsSubmanifold(el)) arret("FEM_ComputeMassMatrix") ;

  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  


  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;
      
        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;
    
        for(p = 0 ; p < np ; p++ , c += dec) {
          double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*c[0]*d ;
          int    i ;
    
          /* axisymetrical or spherical cases */
          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            double rayon = 0 ;
            
            arret("FEM_ComputeMassMatrix: To be improved") ;
        
            for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
            a *= 2*M_PI*rayon ;
            if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
          }
    
          for(i = 0 ; i < nf ; i++) {
            int ii = Element_OverlappingNode(el,i) ;
            int j ;
        
            for(j = 0 ; j < nf ; j++) {
              int jj = Element_OverlappingNode(el,j) ;
              double kij = a * h[i] * h[j] ;
              
              KM(i ,j ) += 0.25 * kij ;
              KM(ii,j ) += 0.25 * kij ;
              KM(i ,jj) += 0.25 * kij ;
              KM(ii,jj) += 0.25 * kij ;
            }
          }
        }
      
        return(km) ;
      }
    }
  }
  
  
  /* 0D */
  if(dim_h == 0) {
    if(nn == 1) {
      KM(0,0) = c[0] ;
      if(Symmetry_IsCylindrical(sym)) KM(0,0) *= 2*M_PI*x[0][0] ;
      else if(Symmetry_IsSpherical(sym)) KM(0,0) *= 4*M_PI*x[0][0]*x[0][0] ;
    } else arret("FEM_ComputeMassMatrix: impossible") ;
    return(km) ;
  }

  /* Regular element */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++ , c += dec) {
      double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*c[0]*d ;
      int    i ;
    
      /* axisymetrical or spherical cases */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double rayon = 0 ;
        
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      for(i = 0 ; i < nf ; i++) {
        int j ;
        
        for(j = 0 ; j < nf ; j++) {
          KM(i,j) += a*h[i]*h[j] ;
        }
      }
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
#define DH(n,i)     (dh[(n)*3 + (i)])
//#define DH(n,i)     (dh[(n)*dim_h + (i)])
#define C(i,j)      (c[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int np  = IntFct_GetNbOfPoints(fi) ;
  int nf  = IntFct_GetNbOfFunctions(fi) ;
  int dim_h = IntFct_GetDimension(fi) ;
  double* weight = IntFct_GetWeight(fi) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*nn*(sizeof(double)) ;
  double* kc = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn*nn ; i++) {
      kc[i] = 0 ;
    }
  }
  
  //if(Element_IsSubmanifold(el)) arret("FEM_ComputeConductionMatrix") ;
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  


  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;
      
        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;
    
        for(p = 0 ; p < np ; p++ , c += dec) {
          double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
          double jcj[3][3] ;
          int    i,j,k,l ;
    
          /* axisymetrical or spherical cases */
          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            double rayon = 0 ;
            
            arret("FEM_ComputeConductionMatrix: To be improved") ;
        
            for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
            a *= 2*M_PI*rayon ;
            if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
          }
    
          /* jcj = J(i,k)*C(k,l)*J(j,l) */
          for(i = 0 ; i < dim_h ; i++) for(j = 0 ; j < dim_h ; j++) {
            jcj[i][j] = 0 ;
            for(k = 0 ; k < dim ; k++) for(l = 0 ; l < dim ; l++)  {
              jcj[i][j] += CAJ(i,k)*C(k,l)*CAJ(j,l) ;
            }
          }
      
          /* KC(i,j) = DH(i,k)*JCJ(k,l)*DH(j,l) */
          for(i = 0 ; i < nf ; i++) {
            int ii = Element_OverlappingNode(el,i) ;
            
            for(j = 0 ; j < nf ; j++) {
              int jj = Element_OverlappingNode(el,j) ;
            
              for(k = 0 ; k < dim_h ; k++) for(l = 0 ; l < dim_h ; l++) {
                double kcij = a * DH(i,k) * jcj[k][l] * DH(j,l) ;
                
                KC(i ,j ) += 0.25 * kcij ;
                KC(ii,j ) += 0.25 * kcij ;
                KC(i ,jj) += 0.25 * kcij ;
                KC(ii,jj) += 0.25 * kcij ;
              }
            }
          }
        }
      
        return(kc) ;
      }
    }
  }
  
  /* Loop on integration points */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++ , c += dec) {
      double* h  = IntFct_GetFunctionAtPoint(fi,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(fi,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
      double jcj[3][3] ;
      int    i,j,k,l ;
    
      /* axisymetrical or spherical cases */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double rayon = 0 ;
        
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      /* jcj = J(i,k)*C(k,l)*J(j,l) */
      for(i = 0 ; i < dim_h ; i++) for(j = 0 ; j < dim_h ; j++) {
        jcj[i][j] = 0 ;
        for(k = 0 ; k < dim ; k++) for(l = 0 ; l < dim ; l++)  {
          jcj[i][j] += CAJ(i,k)*C(k,l)*CAJ(j,l) ;
        }
      }
      
      /* KC(i,j) = DH(i,k)*JCJ(k,l)*DH(j,l) */
      for(i = 0 ; i < nf ; i++) for(j = 0 ; j < nf ; j++) {
        for(k = 0 ; k < dim_h ; k++) for(l = 0 ; l < dim_h ; l++) {
          KC(i,j) += a * DH(i,k) * jcj[k][l] * DH(j,l) ;
        }
      }
    }
  }
  
  return(kc) ;

#undef KC
#undef DH
#undef C
#undef CAJ
}



double*  FEM_ComputePoroelasticMatrix6(FEM_t* fem,IntFct_t* fi,const double* c,const int dec,const int n_dif,const int idis)
/** Return a pointer on a FE poroelastic matrix (Ndof x Ndof).
 * 
 *  Ndof = nb of degrees of freedom (= NN*Neq)
 *  NN   = nb of nodes 
 *  Neq  = nb of equations (= Dim + n_dif)
 *  Dim  = dimension of the problem
 * 
 *  The inputs are:
 * 
 *  n_dif = nb of Biot-like coupling terms (pressure, temperature, etc...)
 * 
 *  idis = position index of the first displacement in the unknown vector
 * 
 *  c = the entry matrix which should be given in the following order:
 * 
 *  K0 to Kn then A0 to An etc.. with n = n_dif:
 * 
 *  | K0(9x9) K1(9x1) K2(9x1) ... | 
 *  | A0(1x9) A1(1x1) A2(1x1) ... | 
 *  | B0(1x9) B1(1x1) B2(1x1) ... | 
 *  | ........................... |
 * 
 *   K0    = Stiffness matrix
 *   Kn    = Mechanic-hydraulic coupling terms
 *   A0,B0 = Hydraulic-mechanic coupling terms
 *   Ai,Bj = Hydraulic terms
 * 
 *  The outputs is provided in an order depending on "idis".
 *  The first displacement unknown U1 will be positionned at "idis".
 *  (1  ... idis ... Neq)
 *   |  ...  |   ...  |
 *   v       v        v
 *  (P1 ...  U1  ...  Pn)
 */
{
#define K(i,j)      (k[(i)*ndof + (j)])
#define E_Mec       indmec
#define I_U         E_Mec
#define E_Hyd(i)    (((i) < E_Mec) ? (i) : dim + (i))
#define I_H(i)      E_Hyd(i)
  int indmec = ((idis < 0) ? 0 : ((idis > n_dif) ? n_dif : idis)) ;
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int neq = dim + n_dif ;
  int ndof = nn*neq ;
  size_t SizeNeeded = ndof*ndof*(sizeof(double)) ;
  double* k = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  int    i,n,m ;
  double zero = 0. ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;
  
  
  if(Element_IsSubmanifold(el)) {
    arret("FEM_ComputePoroelasticMatrix") ;
  }

  /* 
  ** 1.  Mechanics
  */
  /* 1.1 Elasticity */
  {
    double* kr = FEM_ComputeElasticMatrix(fem,fi,c,dec) ;
    
    #define KR(i,j)     (kr[(i)*nn*dim + (j)])
    
    for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
      int j ;
      
      for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
        K(E_Mec + i + n*neq,I_U + j + m*neq) = KR(i + n*dim,j + m*dim) ;
      }
    }
    
    #undef KR
  }
  
  /* 1.2 Biot-like coupling terms */
  for(i = 0 ; i < n_dif ; i++) {
    int I_h = I_H(i) ;
    const double* c1 = c + 81 +i*9 ;
    double* kb = FEM_ComputeBiotMatrix(fem,fi,c1,dec) ;
    
    #define KB(i,j)     (kb[(i)*nn + (j)])
    
    for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        K(E_Mec + j + n*neq,I_h + m*neq) = KB(j + n*dim,m) ;
      }
    }
    
    #undef KB
  }
  

  /* 2. Hydraulic equations */
  for(i = 0 ; i < n_dif ; i++) {
    int E_h = E_Hyd(i) ;
    int j ;
    
    /* 2.1 Coupling terms */
    {
      const double* c1 = c + 81 + n_dif*9 + i*(9 + n_dif) ;
      double* kb = FEM_ComputeBiotMatrix(fem,fi,c1,dec) ;
    
      #define KB(i,j)     (kb[(i)*nn + (j)])
    
      for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
        for(j = 0 ; j < dim ; j++) {
          K(E_h + m*neq,I_U + j + n*neq) = KB(j + n*dim,m) ;
        }
      }
    
      #undef KB
    }
    
    /* 2.2 Accumulations */
    for(j = 0 ; j < n_dif ; j++) {
      int I_h = I_H(j) ;
      const double* c2 = c + 81 + n_dif*9 + i*(9 + n_dif) + 9 + j ;
      double* ka = FEM_ComputeMassMatrix(fem,fi,c2,dec) ;
    
      #define KA(i,j)     (ka[(i)*nn + (j)])
      
      for(n = 0 ; n < nn ; n++) for(m = 0 ; m < nn ; m++) {
        K(E_h + n*neq,I_h + m*neq) = KA(n,m) ;
      }
    
      #undef KA
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
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int nf = IntFct_GetNbOfFunctions(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;

  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) r[i] = zero ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  
  
  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;

        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;

        /* 1D, 2D, 3D */
        for(p = 0 ; p < np ; p++) {
          double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          int i ;
    
          /* cas axisymetrique ou shperique */
          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            double rayon = zero ;
            
            arret("FEM_ComputeBodyForceResidu: to be improved") ;
      
            for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
            a *= 2*M_PI*rayon ;
            if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
          }
    
          /* R(i) = F*H(i) */
          for(i = 0 ; i < nf ; i++) {
            int ii = Element_OverlappingNode(el,i) ;
            double ri = a * f[p*dec] * h[i] ;
            
            r[i ] += 0.5 * ri ;
            r[ii] += 0.5 * ri ;
          }
        }

        return(r) ;
      }
    }
  }


  /* 0D */
  if(dim_h == 0) {
    if(nn == 1) {
      double radius = x[0][0] ;
      
      r[0] = f[0] ;
      
      if(Symmetry_IsCylindrical(sym)) r[0] *= 2*M_PI*radius ;
      else if(Symmetry_IsSpherical(sym)) r[0] *= 4*M_PI*radius*radius ;
      
    } else {
      arret("FEM_ComputeBodyForceResidu: impossible") ;
    }
    
    return(r) ;
  }

  /* 1D, 2D, 3D */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++) {
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      int i ;
    
      /* cas axisymetrique ou shperique */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double rayon = zero ;
      
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      /* R(i) = F*H(i) */
      for(i = 0 ; i < nf ; i++) {
        r[i] += a * f[p*dec] * h[i] ;
      }
    }
  }
  
  return(r) ;
}



double*   FEM_ComputeStrainWorkResidu(FEM_t* fem,IntFct_t* intfct,const double* sig,const int dec)
/* Compute the residu due to strain work */
{
#define DH(n,i)     (dh[(n)*3 + (i)])
//#define DH(n,i)     (dh[(n)*dim_h + (i)])
#define SIG(i,j)    (sig[(i)*3 + (j)])
#define CAJ(i,j)    (caj[(i)*3 + (j)])
#define R(n,i)      (r[(n)*dim + (i)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int nf = IntFct_GetNbOfFunctions(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  int ndof = nn*dim ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = ndof*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = zero ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  
  
  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;

        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;

        for(p = 0 ; p < np ; p++ , sig += dec) {
          double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          double* norm = Element_ComputeNormalVector(el,dh,nf,dim_h) ;
          double sign[3] = {0,0,0} ;
          int i,j ;

          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            arret("FEM_ComputeStrainWorkResidu: to be improved") ;
          }

          /* Compute the vector stress: SIG.N 
           * Pay attention to the orientation of the normal:
           * the normal N is oriented from the opposite face
           * to the main face of the element, i.e. from ii to i */
          #define NORM(i)  (norm[i])
          for(i = 0 ; i < dim ; i++) {
            for(j = 0 ; j < dim ; j++) {
              sign[i] += SIG(i,j) * NORM(j) ;
            }
          }
          #undef NORM

          /* R(i,j) = H(i) * SIGN(j) */
          for(i = 0 ; i < nf ; i++) {
            int ii = Element_OverlappingNode(el,i) ;

            for(j = 0 ; j < dim ; j++) {
              R(i,j)  +=   a * h[i] * sign[j] ;
              R(ii,j) += - a * h[i] * sign[j] ;
            }
          }
        }

        return(r) ;
      }
    }
  }
  
  /* Regular element */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++ , sig += dec) {
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
      double rayon = 0. ;
      int i,j ;
    
      /* The radius in axisymmetrical or spherical case */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        
        a *= 2*M_PI*rayon ;
        
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      /* R(i,j) = DH(i,k) * J(k,l) * S(l,j) */
      for(i = 0 ; i < nf ; i++) for(j = 0 ; j < dim ; j++) {
        int    k,l ;
        
        for(k = 0 ; k < dim_h ; k++) for(l = 0 ; l < dim ; l++) {
          R(i,j) += a * DH(i,k) * CAJ(k,l) * SIG(l,j) ;
        }
      }
    
      /* Axisymmetrical or spherical case */
      if(Symmetry_IsCylindrical(sym)) {
        /* R(i,0) = H(i)/r * S(theta,theta) */
        for(i = 0 ; i < nf ; i++) {
          R(i,0) += a * h[i]/rayon * SIG(2,2) ;
        }
      } else if(Symmetry_IsSpherical(sym)) {
        /* R(i,0) = H(i)/r * (SIG(theta,theta) + SIG(phi,phi)) */
        for(i = 0 ; i < nf ; i++) {
          R(i,0) += a * h[i]/rayon * (SIG(1,1) + SIG(2,2)) ;
        }
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
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*dim_h+(i)])
#define CAJ(i,j)    (caj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int nn  = Element_GetNbOfNodes(el) ;
  int nf = IntFct_GetNbOfFunctions(intfct) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  double* weight = IntFct_GetWeight(intfct) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  size_t SizeNeeded = nn*(sizeof(double)) ;
  double* r = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  double* x[Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) r[i] = zero ;
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      x[i] = Element_GetNodeCoordinate(el,i) ;
    }
  }
  
  
  /* Interface elements */
  if(Element_HasZeroThickness(el)) {
    if(dim_h == dim - 1) {
      /* Assuming that the numbering of the element is formed with the nf first nodes  */
      if(nn > nf) {
        int p ;

        /* Check the numbering */
        FEM_CheckNumberingOfOverlappingNodes(el,nf) ;
    
        for(p = 0 ; p < np ; p++ , f +=dec) {
          double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
          double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
          double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
          double a   = weight[p]*d ;
          double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
          double rayon = 0. ;
          int i ;
    
          /* cas axisymetrique ou shperique */
          if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
            arret("FEM_ComputeFluxResidu: to be improved") ;
            for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
            a *= 2*M_PI*rayon ;
            if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
          }
    
          /* R(i) = DH(i,k)*J(k,j)*F(j) */
          for(i = 0 ; i < nf ; i++) {
            int ii = Element_OverlappingNode(el,i) ;
            int    j,k ;
      
            for(j = 0 ; j < dim ; j++) for(k = 0 ; k < dim_h ; k++) {
              double ri = a * DH(i,k) * CAJ(k,j) * f[j] ;
              
              r[i]  += 0.5 * ri ;
              r[ii] += 0.5 * ri ;
            }
          }
        }

        return(r) ;
      }
    }
  }
  
  /* Regular element */
  {
    int p ;
    
    for(p = 0 ; p < np ; p++ , f +=dec) {
      double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
      double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nf,dim_h) ;
      double a   = weight[p]*d ;
      double* caj = FEM_ComputeInverseJacobianMatrix(el,dh,nf,dim_h) ;
      double rayon = 0. ;
      int i ;
    
      /* cas axisymetrique ou shperique */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        for(i = 0 ; i < nf ; i++) rayon += h[i]*x[i][0] ;
        a *= 2*M_PI*rayon ;
        if(Symmetry_IsSpherical(sym)) a *= 2*rayon ;
      }
    
      /* R(i) = DH(i,k)*J(k,j)*F(j) */
      for(i = 0 ; i < nf ; i++) {
        int    j,k ;
      
        for(j = 0 ; j < dim ; j++) for(k = 0 ; k < dim_h ; k++) {
          r[i] += a * DH(i,k) * CAJ(k,j) * f[j] ;
        }
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
  int nn  = Element_GetNbOfNodes(el) ;
  unsigned short int dim = Geometry_GetDimension(geom) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  Node_t* *no = Element_GetPointerToNode(el) ;
  Field_t* field = Load_GetField(load) ;
  char    *load_eqn = Load_GetNameOfEquation(load) ;
  char    *load_type = Load_GetType(load) ;
  Function_t* function = Load_GetFunction(load) ;
  int dim_h  = IntFct_GetDimension(intfct) ;
  int    nf  = IntFct_GetNbOfFunctions(intfct) ;
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
  if(ieq < 0) arret("FEM_ComputeSurfaceLoadResidu (1): unknown equation") ;

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
          for(j = 0 ; j < nf ; j++) y[i] += h[j]*x[j][i] ;
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
          for(j = 0 ; j < nf ; j++) y[i] += h[j]*x[j][i] ;
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
        double* n = Element_ComputeNormalVector(el,dh,nf,dim_h) ;
        double y[3] = {0.,0.,0.,} ;
        double f0 ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < nf ; j++) y[i] += h[j]*x[j][i] ;
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
      arret("FEM_ComputeSurfaceLoadResidu (2): unknown type") ;
    }
    
    if(jj < 0 || jj >= dim) {
      arret("FEM_ComputeSurfaceLoadResidu (2): unknown type") ;
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
        double* n = Element_ComputeNormalVector(el,dh,nf,dim_h) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          for(j = 0 ; j < nf ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*Field_ComputeValueAtPoint(field,y,dim)*n[jj] ;
      }
      
      rb = FEM_ComputeBodyForceResidu(fem,intfct,f,1) ;
      
      for(j = 0 ; j < nn ; j++) r[j*neq+ieq+ii] = rb[j] ;
      
      FEM_FreeBufferFrom(fem,rb) ;
      return(r) ;
    }
  }
  
  arret("FEM_ComputeSurfaceLoadResidu (2): unknown load") ;
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
    #define N_MODELS (3)
    int n_models = N_MODELS ;
    const char* modelswithstresses[N_MODELS] = {"Elast","Plast","MechaMic"} ;
    const int   stressindex[N_MODELS] = {0,0,0} ;
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
    #undef N_MODELS
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
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
  Element_t* el = FEM_GetElement(fem) ;
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
    Element_t* el = FEM_GetElement(fem) ;
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
        dis[i] = FEM_ComputeUnknown(fem,du,intfct,p,inc+i) ;
      }
    } else 
    #endif
    {
      for(i = 0 ; i < dim ; i++) {
        dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,inc+i) ;
      }
    }
    
    for(i = dim ; i < 3 ; i++) dis[i] = 0 ;
  }
        
  return(dis) ;
}



double* FEM_ComputeUnknownGradient(FEM_t* fem,double** u,IntFct_t* intfct,int p,int inc)
/** Compute the unknown gradient located at "inc" at the interpolation point "p" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  

  {
    int    dim_h  = IntFct_GetDimension(intfct) ;
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    /* interpolation functions */
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
    double* gu = grad ;
    int i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k,l ;
        
      for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim_h ; l++) {
        gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
      }
    }
      
    Element_FreeBufferFrom(el,cj) ;
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
#define U(n,i)   (u[n][Element_GetNodalUnknownPosition(el,n,inc + (i))])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define EPS(i,j) (eps[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int dim_h = IntFct_GetDimension(intfct) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* strain = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    int nn = IntFct_GetNbOfFunctions(intfct) ;
    
    /* interpolation functions */
    double* h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
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
      
      #define NORM(i)  (norm[i])
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          int k ;
          
          for(k = 0 ; k < nn ; k++) {
            int kk = Element_OverlappingNode(el,k) ;
            
            gu[i][j] += (U(k,i) - U(kk,i)) * h[k] * NORM(j) ;
          }
        }
      }
      #undef NORM
      
    } else {
      
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          int k ;
          
          for(k = 0 ; k < nn ; k++) {
            int l ;
          
            for(l = 0 ; l < dim_h ; l++) {
              gu[i][j] += U(k,i) * DH(k,l) * CJ(l,j) ;
            }
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
  
      /* symmetric cases: axisymmetrical or spherical */
      if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
        double rayon = 0. ;
        double u_r   = 0. ;
        
        if(Element_HasZeroThickness(el)) {
          arret("FEM_ComputeLinearStrainTensor: not available yet") ;
        } else {
          for(i = 0 ; i < nn ; i++) {
            double* x = Element_GetNodeCoordinate(el,i) ;
          
            rayon += h[i]*x[0] ;
            u_r   += h[i]*U(i,0) ;
          }
          EPS(2,2) += u_r/rayon ;
          if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
        }
      }
    }
      
    Element_FreeBufferFrom(el,cj) ;
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
#define U(n)   (Element_GetValueOfCurrentNodalUnknown(el,n,inc))
  Element_t* el = FEM_GetElement(fem) ;
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
#define U(n)   (Element_GetValueOfPreviousNodalUnknown(el,n,inc))
  Element_t* el = FEM_GetElement(fem) ;
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
#define U(n)   (Element_GetValueOfCurrentNodalUnknown(el,n,inc))
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)]) 
  Element_t* el = FEM_GetElement(fem) ;
  int    dim_h  = Element_GetDimension(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
    double* gu = grad ;
    int    i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k ;
      
      for(k = 0 ; k < nn ; k++) {
        int l ;
        
        for(l = 0 ; l < dim_h ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
    }
  
    Element_FreeBufferFrom(el,cj) ;
  }
  
  return(grad) ;
  
#undef U
#undef DH
#undef CJ
}



double* FEM_ComputePreviousUnknownGradient(FEM_t* fem,double* dh,int nn,int inc)
/** Compute the previous unknown gradient for the unknown located at "inc" */
{
#define U(n)   (Element_GetValueOfPreviousNodalUnknown(el,n,inc))
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int    dim_h  = Element_GetDimension(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  size_t SizeNeeded = 3*sizeof(double) ;
  double* grad = (double*) FEM_AllocateInBuffer(fem,SizeNeeded) ;
  
  
  {
    /* inverse jacobian matrix (cj) */
    double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
    double* gu = grad ;
    int    i ;
  
    /* initialisation of gu */
    for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
    /* the parameter gradient (gu) */
    for(i = 0 ; i < dim ; i++) {
      int k ;
      
      for(k = 0 ; k < nn ; k++) {
        int l ;
        
        for(l = 0 ; l < dim_h ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
    }
  
    Element_FreeBufferFrom(el,cj) ;
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
#define U(n,i)   (Element_GetValueOfCurrentNodalUnknown(el,n,inc + (i)))
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define EPS(i,j) (eps[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int    dim_h  = Element_GetDimension(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj): Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim_h ; l++) {
      gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases: axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(el,i) ;
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
#define U(n,i)   (Element_GetValueOfCurrentNodalUnknown(el,n,inc + (i)))
#define U_n(n,i) (Element_GetValueOfPreviousNodalUnknown(el,n,inc + (i)))
#define DU(n,i)  (U(n,i) - U_n(n,i))
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define EPS(i,j) (eps[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int    dim_h  = Element_GetDimension(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj): Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim_h ; l++) {
      gu[i][j] += DU(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases: axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(el,i) ;
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
#define U(n,i)   (Element_GetValueOfPreviousNodalUnknown(el,n,inc + (i)))
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define EPS(i,j) (eps[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  size_t SizeNeeded = 9*sizeof(double) ;
  double* eps = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int    dim_h  = Element_GetDimension(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int    i,j ;
  double gu[3][3], *cj ;
  
  /* initialisation of gu */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++)  gu[i][j] = 0. ;
  
  /* inverse jacobian matrix (cj): Note that cj = eps ! */
  cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
  
  /* displacement gradient (gu) */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim ; j++) {
    int    k,l ;
    for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim_h ; l++) {
      gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
    }
  }
  
  /* Linearized strain tensor */
  for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
    EPS(i,j) = (gu[i][j] + gu[j][i])*0.5 ;
  }
  
  /* symmetric cases: axisymmetrical or spherical */
  if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
    double rayon = 0. ;
    double u_r   = 0. ;
    for(i = 0 ; i < nn ; i++) {
      double* x = Element_GetNodeCoordinate(el,i) ;
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
      double d   = FEM_ComputeJacobianDeterminant(el,dh,nn,dim_h) ;
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



double* FEM_ComputeJacobianMatrixNew(Element_t* el,double* dh,int nn,const int dim_h)
{
  return(Element_ComputeJacobianMatrix(el,dh,nn,dim_h)) ;
}




double FEM_ComputeJacobianDeterminantNew(Element_t* el,double* dh,int nn,const int dim_h)
/** Compute the determinant of the jacobian matrix */
{
  double* jac  = FEM_ComputeJacobianMatrixNew(el,dh,nn,dim_h) ;
  double det = Math_Compute3x3MatrixDeterminant(jac) ;
  
  if(det < 0) {
    arret("FEM_ComputeJacobianDeterminant: negative determinant (det = %f)",det) ;
  }
  
  Element_FreeBufferFrom(el,jac) ;
  
  return(det) ;
}



double* FEM_ComputeInverseJacobianMatrixNew(Element_t* el,double* dh,int nn,const int dim_h)
/** Compute the inverse jacobian matrix */
{
  size_t SizeNeeded = 9*sizeof(double) ;
  double* cj = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  
  {
    double* jc = FEM_ComputeJacobianMatrixNew(el,dh,nn,dim_h) ;
    double* b  = Math_Inverse3x3Matrix(jc) ;
    int i ;
    
    if(b) {
      for(i = 0 ; i < 9 ; i++) {
        cj[i] = b[i] ;
      }
    } else {
      arret("FEM_ComputeInverseJacobianMatrix: not invertible") ;
    }
  }
    
  return(cj) ;
}



void FEM_CheckNumberingOfOverlappingNodes(Element_t* el,const int nf)
{
  int j ;

  for(j = 0 ; j < nf ; j++) {
    int jj = Element_OverlappingNode(el,j) ;

    if(jj != j && jj < nf) {
      arret("FEM_CheckNumberingOfOverlappingNodes: bad numbering") ;
    }
  }
}





/* Not used from here */
#if 0
double* FEM_ComputeInverseJacobianMatrixOld(Element_t* el,double* dh,int nn)
/** Compute the inverse jacobian matrix */
{
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  size_t SizeNeeded = 9*sizeof(double) ;
  double* cj = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int    dim_h = Element_GetDimension(el) ;
  int    dim = Geometry_GetDimension(Element_GetGeometry(el));
  double jc[3][3],dt,td ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < 9 ; i++) cj[i] = 0. ;
  }
  
  
  /* The jacobian matrix */
  {
    double* jc1 = FEM_ComputeJacobianMatrixOld(el,dh,nn) ;
    int i ;
    
    for(i = 0 ; i < dim ; i++) {
      int j ;
      
      for(j = 0 ; j < dim_h ; j++) {
        jc[i][j] = jc1[i*3 + j] ;
      }
    }
  
    Element_FreeBufferFrom(el,jc1) ;
  }


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
#endif



#if 0
double* FEM_ComputeJacobianMatrixOld(Element_t* el,double* dh,int nn)
/** Compute the jacobian matrix */
{
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define JC(i,j)  (jc[(i)*3 + (j)])
  size_t SizeNeeded = 9*sizeof(double) ;
  double* jc = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int  dim_h = Element_GetDimension(el) ;
  int  dim   = Geometry_GetDimension(Element_GetGeometry(el));
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < 9 ; i++) jc[i] = 0. ;
  }
  
  
  /* The Jacobian matrix */
  {
    int i ;
    
    for(i = 0 ; i < dim ; i++) {
      int    j ;
    
      for(j = 0 ; j < dim_h ; j++) {
        int k ;
    
        JC(i,j) = 0. ;
        
        for(k = 0 ; k < nn ; k++) {
          double* x = Element_GetNodeCoordinate(el,k) ;
      
          JC(i,j) += x[i]*DH(k,j) ;
        }
      }
    }
  }
  
  return(jc) ;
#undef DH
#undef JC
}
#endif



#if 0
double FEM_ComputeJacobianDeterminantOld(Element_t* el,double* dh,int nn)
/** Compute the determinant of the jacobian matrix */
{
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define J(i,j)   (jac[(i)*3 + (j)])
  int    dim_h = Element_GetDimension(el) ;
  int    dim   = Geometry_GetDimension(Element_GetGeometry(el));
  double* jac  = FEM_ComputeJacobianMatrixOld(el,dh,nn) ;
  double det ;
  
  /* 
     1. For a volume: it's det(J1,J2,J3)
     2. For a surface: it's det(J1,J2,e3) 
        where e3 is the unit normal vector to (J1,J2)
     3. For a line   : it's det(J1,e2,e3) = |J1|
        where (e2,e3) are orthonormed vectors in the normal plane to J1
     4. For a point  : it's det(e1,e2,e3) = 1
        where (e1,e2,e3) is an orthonormed system.
  */

  /* 1. Volume: det(J1,J2,J3) */
  if(dim_h == 3) {
    if(dim == 3) {
      det  = J(0,0)*J(1,1)*J(2,2) - J(0,0)*J(2,1)*J(1,2)
           + J(1,0)*J(2,1)*J(0,2) - J(1,0)*J(0,1)*J(2,2)
           + J(2,0)*J(0,1)*J(1,2) - J(2,0)*J(1,1)*J(0,2) ;
      det  = fabs(det) ;
      
    } else {
      arret("FEM_ComputeJacobianDeterminant (1)") ;
    }

  /* 2. Surface: det(J1,J2,e3) = |c1^c2| */
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

  /* 3. Line: det(J1,e2,e3) = |J1|(e1,e2,e3) = |J1|   */
  } else if(dim_h == 1) {
    int i ;
    
    det = 0. ;
    for(i = 0 ; i < dim ; i++) {
      det += J(i,0)*J(i,0) ;
    }
    
    det = sqrt(det) ;

  /* 4. Point: det(e1,e2,e3) = 1 */
  } else if(dim_h == 0) {
    det = 1. ;
    
  } else {
  
    arret("FEM_ComputeJacobianDeterminant (3)") ;
  }
  
  Element_FreeBufferFrom(el,jac) ;
  
  return(det) ;
#undef DH
#undef J
}
#endif




/* FEM_Compute functions in the form (FEM_t*,IntFct_t*,double**,int) */

/* NOT USED */
#if 0
double*  FEM_ComputeNormalVector(Element_t* el,double* dh,int nn)
/* Normale unitaire a un sous-espace de dimension dim-1 */
{
#define DH(n,i) (dh[(n)*3+(i)])
//#define DH(n,i) (dh[(n)*dim_h+(i)])
  size_t SizeNeeded = 3*sizeof(double) ;
  double* norm = (double*) Element_AllocateInBuffer(el,SizeNeeded) ;
  int    dim_h = Element_GetDimension(el) ;
  int    dim   = dim_h + 1 ;
  int    i,j ;
  double c[3][3] ;

  if(dim > 3) arret("FEM_ComputeNormalVector") ;
  
  /* le jacobien */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
    int    k ;
    c[i][j] = 0. ;
    for(k = 0 ; k < nn ; k++) {
      double* x = Element_GetNodeCoordinate(el,k) ;
      c[i][j] += x[i]*DH(k,j) ;
    }
  }
  
  /* 1. Surface: norm = c1^c2 */
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

  /* 2. Ligne: norm = c1^e_z */
  } else if(dim_h == 1) {
    double v = sqrt(c[0][0]*c[0][0] + c[1][0]*c[1][0]) ;
    norm[0] =  c[1][0]/v ;
    norm[1] = -c[0][0]/v ;
    norm[2] = 0. ;
    return(norm) ;

  /* 3. Point: norm = 1 */
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



#if 0
double* FEM_ComputeUnknowns(FEM_t* fem,IntFct_t* intfct,double** u,int inc)
/** Compute the unknowns at the interpolation points located at "inc" */
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
  Element_t* el = FEM_GetElement(fem) ;
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
/** Compute the unknown gradients at the interpolation points located at "inc" 
 *  WRONG IMPLEMENTATION: check inverse jacobian matrix (dim/dim_h pb)*/
{
#define U(n)   (u[n][Element_GetNodalUnknownPosition(el,n,inc)])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  int    dim_h  = Element_GetDimension(el) ;
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
      double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
      double gu[3] ;
      int i ;
  
      /* initialisation of gu */
      for(i = 0 ; i < 3 ; i++)  gu[i] = 0. ;
  
      /* the parameter gradient (gu) */
      for(i = 0 ; i < dim_h ; i++) {
        int k,l ;
        
        for(k = 0 ; k < nn ; k++) for(l = 0 ; l < dim_h ; l++) {
          gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
        }
      }
      
      {
        double* grad = ugrd + 3*p ;
        
        for(i = 0 ; i < 3 ; i++) grad[i] = gu[i] ;
      }
      
      Element_FreeBufferFrom(el,cj) ;
    }
  }
  
  return(ugrd) ;
  
#undef U
#undef DH
#undef CJ
}



double* FEM_ComputeLinearStrainTensors(FEM_t* fem,IntFct_t* intfct,double** u,int inc)
/** Compute the 3D linearized strain tensors for the displacement vectors 
 *  located at "inc" and at the interpolation points.
 *  WRONG IMPLEMENTATION: check inverse jacobian matrix (dim/dim_h pb)*/
{
#define U(n,i)   (u[n][Element_GetNodalUnknownPosition(el,n,inc + (i))])
#define DH(n,i)  (dh[(n)*3 + (i)])
//#define DH(n,i)  (dh[(n)*dim_h + (i)])
#define EPS(i,j) (eps[(i)*3 + (j)])
#define CJ(i,j)  (cj[(i)*3 + (j)])
  Element_t* el = FEM_GetElement(fem) ;
  Symmetry_t sym = Element_GetSymmetry(el) ;
  int    dim_h  = Element_GetDimension(el) ;
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
      double* cj = FEM_ComputeInverseJacobianMatrix(el,dh,nn,dim_h) ;
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
      for(i = 0 ; i < dim_h ; i++) {
        int j ;
        
        for(j = 0 ; j < dim_h ; j++) {
          int k ;
          
          for(k = 0 ; k < nn ; k++) {
            int l ;
          
            for(l = 0 ; l < dim_h ; l++) {
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
  
        /* symmetric cases: axisymmetrical or spherical */
        if(Symmetry_IsCylindrical(sym) || Symmetry_IsSpherical(sym)) {
          double rayon = 0. ;
          double u_r   = 0. ;
        
          for(i = 0 ; i < nn ; i++) {
            double* x = Element_GetNodeCoordinate(el,i) ;
          
            rayon += h[i]*x[0] ;
            u_r   += h[i]*U(i,0) ;
          }
          EPS(2,2) += u_r/rayon ;
          if(Symmetry_IsSpherical(sym)) EPS(1,1) += u_r/rayon ;
        }
      }
      
      Element_FreeBufferFrom(el,cj) ;
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
