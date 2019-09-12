#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Message.h"
#include "Tools/Math.h"
#include "IntFct.h"

static void   IntFct_ComputeAtGaussPoints(IntFct_t*,int,int) ;
static void   IntFct_ComputeAtNodes(IntFct_t*,int,int) ;
static void   gaussp(int,double*,double*) ;
static void   hammer(int,double*,double*) ;
static void   gauss_tetraedre(int,double*,double*) ;
static void   (IntFct_ComputeIsoShapeFct)(int,int,double*,double*,double*) ;
static void   IntFct_ComputeAtMidSurfacePoints(IntFct_t*,int,int) ;
static void   midpoints(double*,double) ;
static void   normale(int,int,double**,double*,double*) ;
static void   (IntFct_AllocateMemory)(IntFct_t*) ;


/* Extern functions */

IntFct_t* (IntFct_Create)(int nn,int dim,const char* type)
{
  IntFct_t* intfct = (IntFct_t*) malloc(sizeof(IntFct_t)) ;
  
  if(!intfct) {
    arret("IntFct_Create(1)") ;
  }
  
  IntFct_GetDimension(intfct) = dim ;
  IntFct_GetNbOfNodes(intfct) = nn ;
  
  IntFct_AllocateMemory(intfct) ;
  
  strcpy(IntFct_GetType(intfct),type) ;
    
  if(IntFct_TypeIs(intfct,"Nodes")) {
    IntFct_ComputeAtNodes(intfct,nn,dim) ;
  } else if(IntFct_TypeIs(intfct,"Gauss")) {
    IntFct_ComputeAtGaussPoints(intfct,nn,dim) ;
  } else if(IntFct_TypeIs(intfct,"MidSurface")) {
    IntFct_ComputeAtMidSurfacePoints(intfct,nn,dim) ;
  } else {
    arret("IntFct_Create(3): type \"%s\" unknown",type) ;
  }
  
  return(intfct) ;
}




void (IntFct_AllocateMemory)(IntFct_t* intfct)
{
  
  {
    char* p = (char*) malloc(IntFct_MaxLengthOfKeyWord*sizeof(char)) ;
    
    if(!p) arret("IntFct_AllocateMemory(1)") ;
    
    IntFct_GetType(intfct) = p ;
  }
  
  {
    int dim = IntFct_GetDimension(intfct) ;
    int nn  = IntFct_GetNbOfNodes(intfct) ;
    int np  = IntFct_MaxNbOfIntPoints ;
    int k   = np*(1 + nn*(1 + dim) + dim) ;
    double* weight = (double*) malloc(k*sizeof(double)) ;
    
    if(!weight) arret("IntFct_AllocateMemory(2)") ;
    
    IntFct_GetWeight(intfct)           = weight ;
    IntFct_GetFunction(intfct)         = weight + np ;
    IntFct_GetFunctionGradient(intfct) = weight + np*(1 + nn) ;
    IntFct_GetPointCoordinates(intfct) = weight + np*(1 + nn*(1 + dim)) ;
  }
  
  return ;
}




double (IntFct_InterpolateAtPoint)(IntFct_t* intfct,double* var,int shift,int p)
/** Interpolate from nodal values stored in "var". 
 *  Return the interpolated value at the point p. */
{
#define U(n)   (var[(n)*shift])
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



void (IntFct_ComputeIsoShapeFct)(int dim,int nn,double* x,double* h,double* dh)
/* Compute Iso-shape functions (h) and their gradients (dh) at point x */
{
#define X         x[0]
#define Y         x[1]
#define Z         x[2]
#define DH(n,i)  (dh[(n)*dim+(i)])
#define DHx(n)    DH(n,0)
#define DHy(n)    DH(n,1)
#define DHz(n)    DH(n,2)
  
  /* 0D */
  if(dim == 0) {
    /* bord a 1 noeud */
    if(nn == 1) {
      h[0] = 1. ;
      return ;
    }
    
  /* 1D */
  } else if(dim == 1) {
    /* segment a 2 noeuds */
    if(nn == 2) {
      h[0]    = (1. - X)*0.5 ;
      h[1]    = (1. + X)*0.5 ;
      if(dh) {
        DHx(0)  = -0.5 ;
        DHx(1)  =  0.5 ;
      }
      return ;
    /* segment a 3 noeuds */
    } else if(nn == 3) {
      h[0]    = -X*(1. - X)*0.5 ;
      h[1]    =  X*(1. + X)*0.5 ;
      h[2]    = (1. - X*X) ;
      if(dh) {
        DHx(0)  = -0.5 + X ;
        DHx(1)  =  0.5 + X ;
        DHx(2)  = -2.*X ;
      }
      return ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* triangle a 3 noeuds */
    if(nn == 3) {
      h[0]    = 1. - X - Y ;
      h[1]    = X ;
      h[2]    = Y ;
      if(dh) {
        DHx(0)  = -1. ;
        DHy(0)  = -1. ;
        DHx(1)  =  1. ;
        DHy(1)  =  0. ;
        DHx(2)  =  0. ;
        DHy(2)  =  1. ;
      }
      return ;
    /* quadrilatere a 4 noeuds */
    } else if(nn == 4) {
      h[0]    =  (1. - X)*(1. - Y)*0.25 ;
      h[1]    =  (1. + X)*(1. - Y)*0.25 ;
      h[2]    =  (1. + X)*(1. + Y)*0.25 ;
      h[3]    =  (1. - X)*(1. + Y)*0.25 ;
      if(dh) {
        DHx(0)  = -(1. - Y)*0.25 ;
        DHy(0)  = -(1. - X)*0.25 ;
        DHx(1)  =  (1. - Y)*0.25 ;
        DHy(1)  = -(1. + X)*0.25 ;
        DHx(2)  =  (1. + Y)*0.25 ;
        DHy(2)  =  (1. + X)*0.25 ;
        DHx(3)  = -(1. + Y)*0.25 ;
        DHy(3)  =  (1. - X)*0.25 ;
      }
      return ;
    /* triangle a 6 noeuds */
    } else if(nn == 6) {
      h[0]    =  (1. - X - Y)*(1. - 2*X - 2*Y) ;
      h[1]    = -(1. - 2*X)*X ;
      h[2]    = -(1. - 2*Y)*Y ;
      h[3]    =  (1. - X - Y)*X*4 ;
      h[4]    =  X*Y*4 ;
      h[5]    =  (1. - X - Y)*Y*4 ;
      if(dh) {
        DHx(0)  = -3 + 4*X + 4*Y ;
        DHy(0)  = -3 + 4*X + 4*Y ;
        DHx(1)  = -1. + 4*X ;
        DHy(1)  =  0. ;
        DHx(2)  =  0. ;
        DHy(2)  = -1. + 4*Y ;
        DHx(3)  =  4*(1. - 2*X - Y) ;
        DHy(3)  = -4*X ;
        DHx(4)  =  4*Y ;
        DHy(4)  =  4*X ;
        DHx(5)  = -4*Y ;
        DHy(5)  =  4*(1. - X - 2*Y) ;
      }
      return ;
    /* quadrilatere a 8 noeuds */
    } else if(nn == 8) {
      double DX = 2*X,X2 = X*X ;
      double DY = 2*Y,Y2 = Y*Y ;
      h[0]    = -(1. - X)*(1. - Y)*(1. + X + Y)*0.25 ;
      h[1]    = -(1. + X)*(1. - Y)*(1. - X + Y)*0.25 ;
      h[2]    = -(1. + X)*(1. + Y)*(1. - X - Y)*0.25 ;
      h[3]    = -(1. - X)*(1. + Y)*(1. + X - Y)*0.25 ;
      h[4]    =  (1. - X2)*(1. - Y )*0.5 ;
      h[5]    =  (1. + X )*(1. - Y2)*0.5 ;
      h[6]    =  (1. - X2)*(1. + Y )*0.5 ;
      h[7]    =  (1. - X )*(1. - Y2)*0.5 ;
      if(dh) {
        DHx(0)  =  (1. - Y)*( DX + Y )*0.25 ;
        DHy(0)  =  (1. - X)*( X  + DY)*0.25 ;
        DHx(1)  = -(1. - Y)*(-DX + Y )*0.25 ;
        DHy(1)  =  (1. + X)*(-X  + DY)*0.25 ;
        DHx(2)  = -(1. + Y)*(-DX - Y )*0.25 ;
        DHy(2)  = -(1. + X)*(-X  - DY)*0.25 ;
        DHx(3)  =  (1. + Y)*( DX - Y )*0.25 ;
        DHy(3)  = -(1. - X)*( X  - DY)*0.25 ;
        DHx(4)  = -(1. - Y)*X ;
        DHy(4)  = -(1. - X2)*0.5 ;
        DHx(5)  =  (1. - Y2)*0.5 ;
        DHy(5)  = -(1. + X)*Y ;
        DHx(6)  = -(1. + Y)*X ;
        DHy(6)  =  (1. - X2)*0.5 ;
        DHx(7)  = -(1. - Y2)*0.5 ;
        DHy(7)  = -(1. - X)*Y ;
      }
      return ;
    }
    
    /* 3D */
  } else if(dim == 3) {
    /* tetraedre a 4 noeuds */
    if(nn == 4) {
      h[0]    = 1. - X - Y - Z;
      h[1]    = X ;
      h[2]    = Y ;
      h[3]    = Z ;
      if(dh) {
        DHx(0)  = -1. ;
        DHy(0)  = -1. ;
        DHz(0)  = -1. ;
        DHx(1)  =  1. ;
        DHy(1)  =  0. ;
        DHz(1)  =  0. ;
        DHx(2)  =  0. ;
        DHy(2)  =  1. ;
        DHz(2)  =  0. ;
        DHx(3)  =  0. ;
        DHy(3)  =  0. ;
        DHz(3)  =  1. ;
      }
      return ;
    /* hexaedre a 8 noeuds */
    } else if(nn == 8) {
      h[0]    =  (1. - X)*(1. - Y)*(1. - Z)*0.125 ;
      h[1]    =  (1. + X)*(1. - Y)*(1. - Z)*0.125 ;
      h[2]    =  (1. + X)*(1. + Y)*(1. - Z)*0.125 ;
      h[3]    =  (1. - X)*(1. + Y)*(1. - Z)*0.125 ;
      h[4]    =  (1. - X)*(1. - Y)*(1. + Z)*0.125 ;
      h[5]    =  (1. + X)*(1. - Y)*(1. + Z)*0.125 ;
      h[6]    =  (1. + X)*(1. + Y)*(1. + Z)*0.125 ;
      h[7]    =  (1. - X)*(1. + Y)*(1. + Z)*0.125 ;
      if(dh) {
        DHx(0)  = -(1. - Y)*(1. - Z)*0.125 ;
        DHy(0)  = -(1. - X)*(1. - Z)*0.125 ;
        DHz(0)  = -(1. - X)*(1. - Y)*0.125 ;
        DHx(1)  =  (1. - Y)*(1. - Z)*0.125 ;
        DHy(1)  = -(1. + X)*(1. - Z)*0.125 ;
        DHz(1)  = -(1. + X)*(1. - Y)*0.125 ;
        DHx(2)  =  (1. + Y)*(1. - Z)*0.125 ;
        DHy(2)  =  (1. + X)*(1. - Z)*0.125 ;
        DHz(2)  = -(1. + X)*(1. + Y)*0.125 ;
        DHx(3)  = -(1. + Y)*(1. - Z)*0.125 ;
        DHy(3)  =  (1. - X)*(1. - Z)*0.125 ;
        DHz(3)  = -(1. - X)*(1. + Y)*0.125 ;
        DHx(4)  = -(1. - Y)*(1. + Z)*0.125 ;
        DHy(4)  = -(1. - X)*(1. + Z)*0.125 ;
        DHz(4)  =  (1. - X)*(1. - Y)*0.125 ;
        DHx(5)  =  (1. - Y)*(1. + Z)*0.125 ;
        DHy(5)  = -(1. + X)*(1. + Z)*0.125 ;
        DHz(5)  =  (1. + X)*(1. - Y)*0.125 ;
        DHx(6)  =  (1. + Y)*(1. + Z)*0.125 ;
        DHy(6)  =  (1. + X)*(1. + Z)*0.125 ;
        DHz(6)  =  (1. + X)*(1. + Y)*0.125 ;
        DHx(7)  = -(1. + Y)*(1. + Z)*0.125 ;
        DHy(7)  =  (1. - X)*(1. + Z)*0.125 ;
        DHz(7)  =  (1. - X)*(1. + Y)*0.125 ;
      }
      return ;
    }
  }

  
  arret("IntFct_ComputeIsoShapeFct") ;
  
  return ;
#undef X
#undef Y
#undef Z
#undef DH
#undef DHx
#undef DHy
#undef DHz
}



int (IntFct_ComputeFunctionIndexAtPointOfReferenceFrame)(IntFct_t* intfct,double* x)
/* Compute function values at point x of the reference frame. 
 * Return the index of the point where values are recorded. */
{
  int np = IntFct_GetNbOfPoints(intfct) ;
  int p = np ;
  
  if(p + 1 > IntFct_MaxNbOfIntPoints) {
    arret("IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(1)") ;
  }
  
  {
    int     dim = IntFct_GetDimension(intfct) ;
    int     nn  = IntFct_GetNbOfFunctions(intfct) ;
    double* h   = IntFct_GetFunctionAtPoint(intfct,p) ;
    double* dh  = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    IntFct_ComputeIsoShapeFct(dim,nn,x,h,dh) ;
  }
  
  return(p) ;
}



/* Intern functions */


void   IntFct_ComputeAtGaussPoints(IntFct_t* f,int nn,int dim)
/** Interpolation functions at Gauss points */
{
  /* 0D */
  if(dim == 0) {
    /* 1 noeud */
    if(nn == 1) {
      int np = 1 ;
      double* w = IntFct_GetWeight(f) ;
      double* h = IntFct_GetFunction(f) ;
      
      IntFct_GetNbOfPoints(f) = np  ;
      h[0] = 1. ;
      w[0] = 1 ;
      
      return ;
    }
    
  /* 1D */
  } else if(dim == 1) {
    /* segment a 2 ou 3 noeuds */
    if(nn == 2 || nn == 3) {
      int np = 2 ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      int p ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) {
        arret("IntFct_ComputeAtGaussPoints (1)") ;
      }
      gaussp(np,a,w) ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* triangle a 3 ou 6 noeuds */
    if(nn == 3 || nn == 6) {
      int p,np = 3 ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) arret("IntFct_ComputeAtGaussPoints (2)") ;
      hammer(np,a,w) ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    /* quadrilatere a 4 ou 8 noeuds */
    } else if(nn == 4 || nn == 8) {
      int i,j,np = 4 ;
      double a[IntFct_MaxNbOfIntPoints] ;
      double w[IntFct_MaxNbOfIntPoints] ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) arret("IntFct_ComputeAtGaussPoints (3)") ;
      gaussp(2,a,w) ;
      
      for(i = 0 ; i < 2 ; i++) for(j = 0 ; j < 2 ; j++) {
        int p = 2*i + j ;
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        double* x = IntFct_GetPointCoordinates(f) + p*dim ;
        
        IntFct_GetWeight(f)[p]   = w[i]*w[j] ;
        x[0] = a[i] ;
        x[1] = a[j] ;
        IntFct_ComputeIsoShapeFct(dim,nn,x,h,dh) ;
      }
      
      return ;
    }
    
  /* 3D */
  } else if(dim == 3) {
    /* tetraedre a 4 ou 10 noeuds */
    if(nn == 4 || nn == 10) {
      int p,np = 4 ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) arret("IntFct_ComputeAtGaussPoints (4)") ;
      gauss_tetraedre(np,a,w) ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    /* hexaedre a 8 ou 27 noeuds */
    } else if(nn == 8 || nn == 27) {
      int i,j,k,np = 8 ;
      double a[IntFct_MaxNbOfIntPoints] ;
      double w[IntFct_MaxNbOfIntPoints] ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) arret("IntFct_ComputeAtGaussPoints (5)") ;
      gaussp(2,a,w) ;
      
      for(i = 0 ; i < 2 ; i++) for(j = 0 ; j < 2 ; j++) for(k = 0 ; k < 2 ; k++) {
        int p = 4*i + 2*j + k ;
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        double* x = IntFct_GetPointCoordinates(f) + p*dim ;
        
        IntFct_GetWeight(f)[p]   = w[i]*w[j]*w[k] ;
        
        x[0] = a[i] ;
        x[1] = a[j] ;
        x[2] = a[k] ;
        
        IntFct_ComputeIsoShapeFct(dim,nn,x,h,dh) ;
      }
      
      return ;
    }
  }
  
  arret("IntFct_ComputeAtGaussPoints (6)") ;
}



void   IntFct_ComputeAtNodes(IntFct_t* f,int nn,int dim)
/* Fonctions d'interpolation aux noeuds */
{
  int    np = nn ;

  IntFct_GetNbOfPoints(f) = np  ;
  if(np > IntFct_MaxNbOfIntPoints) arret("IntFct_ComputeAtNodes (1)") ;
  
  /* 0D */
  if(dim == 0) {
    /* 1 noeud */
    if(nn == 1) {
      double* h  = IntFct_GetFunction(f) ;
      
      h[0] = 1. ;
      
      return ;
    }
    
  /* 1D */
  } else if(dim == 1) {
    /* segment a 2 */
    if(nn == 2) {
      int p ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      a[0] = -1. ;
      a[1] =  1. ;
      w[0] = 1. ;
      w[1] = 1. ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p,h,dh) ;
      }
      
      return ;
    }
    
  /* 2D */
  } else if(dim == 2) {
    /* triangle a 3 */
    if(nn == 3) {
      int p ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      a[0]  = 0.  ; a[1]  = 0.  ;
      a[2]  = 1.  ; a[3]  = 0.  ;
      a[4]  = 0.  ; a[5]  = 1.  ;
      w[0] = 1./6.   ;
      w[1] = w[0] ;
      w[2] = w[0] ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    /* quadrilatere a 4 */
    } else if(nn == 4) {
      int p ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      a[0]  = -1.  ; a[1]  = -1.  ;
      a[2]  =  1.  ; a[3]  = -1.  ;
      a[4]  = -1.  ; a[5]  =  1.  ;
      a[6]  =  1.  ; a[7]  =  1.  ;
      w[0] = 1.      ;
      w[1] = w[0] ;
      w[2] = w[0] ;
      w[3] = w[0] ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    }
    
  /* 3D */
  } else if(dim == 3) {
    /* tetraedre a 4 */
    if(nn == 4) {
      int p ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      a[0]  = 0.  ; a[1]  = 0.  ; a[2]  = 0.  ;
      a[3]  = 1.  ; a[4]  = 0.  ; a[5]  = 0.  ;
      a[6]  = 0.  ; a[7]  = 1.  ; a[8]  = 0.  ;
      a[9]  = 0.  ; a[10] = 0.  ; a[11] = 1.  ;
      w[0] = 1./24.  ;
      w[1] = w[0] ;
      w[2] = w[0] ;
      w[3] = w[0] ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    /* hexaedre a 8 */
    } else if(nn == 8) {
      int p ;
      //double a[3*IntFct_MaxNbOfIntPoints] ;
      double* a = IntFct_GetPointCoordinates(f) ;
      double* w = IntFct_GetWeight(f) ;
      
      a[0]  = -1.  ; a[1]  = -1.  ; a[2]  = -1.  ;
      a[3]  =  1.  ; a[4]  = -1.  ; a[5]  = -1.  ;
      a[6]  =  1.  ; a[7]  =  1.  ; a[8]  = -1.  ;
      a[9]  = -1.  ; a[10] =  1.  ; a[11] = -1.  ;
      a[0]  = -1.  ; a[1]  = -1.  ; a[2]  =  1.  ;
      a[3]  =  1.  ; a[4]  = -1.  ; a[5]  =  1.  ;
      a[6]  =  1.  ; a[7]  =  1.  ; a[8]  =  1.  ;
      a[9]  = -1.  ; a[10] =  1.  ; a[11] =  1.  ;
      w[0] = 1.      ;
      w[1] = w[0] ;
      w[2] = w[0] ;
      w[3] = w[0] ;
      w[4] = w[0] ;
      w[5] = w[0] ;
      w[6] = w[0] ;
      w[7] = w[0] ;
      
      for(p = 0 ; p < np ; p++) {
        double* h  = IntFct_GetFunction(f)  + p*nn ;
        double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
      }
      
      return ;
    }
  }
  
  arret("IntFct_ComputeAtNodes (2)") ;
}




void   IntFct_ComputeAtMidSurfacePoints(IntFct_t* f,int nn,int dim)
/** Interpolation functions at mid points of surface */
{
  /* 0D/1D/2D/3D */
  if(dim >= 0 || dim <= 3) {
    if(nn == dim + 1) {
      int np = (dim > 1) ? dim + 1 : 1 ;
      double* a = IntFct_GetPointCoordinates(f) ;
      //double* w = IntFct_GetWeight(f) ;
      
      IntFct_GetNbOfPoints(f) = np ;
      if(np > IntFct_MaxNbOfIntPoints) {
        arret("IntFct_ComputeAtMidSurfacePoints(4)") ;
      }
      
      midpoints(a,dim) ;
      
      {
        int p ;
        
        for(p = 0 ; p < np ; p++) {
          double* h  = IntFct_GetFunction(f)  + p*nn ;
          double* dh = IntFct_GetFunctionGradient(f) + p*nn*dim ;
        
          IntFct_ComputeIsoShapeFct(dim,nn,a+p*dim,h,dh) ;
        }
      }
      
      return ;
    }
  }
  
  arret("IntFct_ComputeAtMidSurfacePoints") ;
}



void  gaussp(int n,double* a,double* w)
/* Points de Gauss */
{
  if(n == 1) {
    a[0] = 0. ;
    w[0] = 2. ;
    return ;
  } else if(n == 2) {
    a[0] = -0.5773502691896 ;
    a[1] = -a[0] ;
    w[0] = 1. ;
    w[1] = 1. ;
    return ;
  } else if(n == 3) {
    a[0] = -0.7745966692415 ;
    a[1] = 0. ;
    a[2] = -a[0] ;
    w[0] = 0.5555555555556 ;
    w[1] = 0.8888888888889 ;
    w[2] = w[0] ;
    return ;
  } else if(n == 4) {
    a[0] = -0.8611363115941 ;
    a[1] = -0.3399810435849 ;
    a[2] = -a[0] ;
    a[3] = -a[1] ;
    w[0] = 0.3478548451375 ;
    w[1] = 0.6521451548625 ;
    w[2] = w[0] ;
    w[3] = w[1] ;
    return ;
  } else if(n == 5) {
    a[0] = -0.9061798459387 ;
    a[1] = -0.5384693101057 ;
    a[2] = 0. ;
    a[3] = -a[0] ;
    a[4] = -a[1] ;
    w[0] = 0.2369268850562 ;
    w[1] = 0.4786286704994 ;
    w[2] = 0.5688888888889 ;
    w[3] = w[0] ;
    w[4] = w[1] ;
    return ;
  }
  arret("gaussp") ;
}

void hammer(int n,double* a,double* w)
/* Points d'integration sur un triangle */
{
  if(n == 1) {
    a[0] = 1./3. ;  a[1] = a[0] ;
    w[0] = 0.5 ;
    return ;
  } else if(n == 3) {
    a[0] = 1./6. ;  a[1] = a[0] ;
    a[2] = 2./3. ;  a[3] = a[0] ;
    a[4] = a[0]  ;  a[5] = a[2] ;
    w[0] = a[0] ;
    w[1] = w[0] ;
    w[2] = w[0] ;
    return ;
  } else if(n == 6) {
    a[0] = 0.445948490915965 ;  a[1] = a[0] ;
    a[2] = 1. - 2*a[0]       ;  a[3] = a[0] ;
    a[4] = a[0] ;            ;  a[5] = a[2] ;
    a[6] = 0.091576213509771 ;  a[7] = a[6] ;
    a[8] = 1. - 2*a[6]       ;  a[9] = a[6] ;
    a[10] = a[6]             ;  a[11] = a[8] ;
    w[0] = 0.111690794839005 ;
    w[1] = w[0] ;
    w[2] = w[0] ;
    w[3] = 0.054975871827661 ;
    w[4] = w[3] ;
    w[5] = w[3] ;
    return ;
  }
  arret("hammer") ;
}

void gauss_tetraedre(int n,double* a,double* w)
/* Points d'integration sur un tetraedre */
{
  if(n == 1) {
    a[0] = 1./4. ;  a[1] = a[0] ;  a[2] = a[0] ;
    w[0] = 1./6. ;
    return ;
  } else if(n == 4) {
    a[0] = (5. - sqrt(5.))/20. ;  a[1] = a[0] ;  a[2] = a[0] ;
    a[3] = (5. + 3.*sqrt(5.))/20. ;  a[4] = a[0] ;  a[5] = a[0] ;
    a[6] = a[0]  ;  a[7] = a[3] ;  a[8] = a[0] ;
    a[9] = a[0]  ;  a[10] = a[0] ;  a[11] = a[3] ;
    w[0] = 1./24. ;
    w[1] = w[0] ;
    w[2] = w[0] ;
    w[3] = w[0] ;
    return ;
  }
  arret("gauss_tetraedre") ;
}

void midpoints(double* a,double dim)
{
  if(dim == 0) {
    a[0] = 0 ;
  } else if(dim == 1) {
    a[0] = 0 ;
  } else if(dim == 2) {
    a[0] = 0.5 ; a[1] = 0 ;
    a[2] = 0   ; a[3] = 0.5 ;
    a[4] = 0.5 ; a[5] = 0.5 ;
  } else if(dim == 3) {
    a[0] = 0.5 ; a[1]  = 0   ; a[2]  = 0 ;
    a[3] = 0   ; a[4]  = 0.5 ; a[5]  = 0 ;
    a[6] = 0   ; a[7]  = 0   ; a[8]  = 0.5 ;
    a[9] = 0.5 ; a[10] = 0.5 ; a[11] = 0.5 ;
  } else {
    arret("midpoints") ;
  }
  
  return ;
}



void  IntFct_ComputeIsoShapeFctInActualSpace(int dim,int nn,double** x_e,double* x,int dim_e,double* h,double* dh)
/* fonction d'interpolation (h) et ses derivees (dh) en x */
{
#define DH(n,i)  (dh[(n)*dim_e+(i)])
  int    i,iter,max_iter = 10 ;
  double r[3],k[9],a[3] = {0.,0.,0},unorm[3] ;
  double x_max[3],x_min[3],d ;
  double err,tol = 1.e-6 ;


  /* diametre de l'element */
  d = 0. ;
  for(i = 0 ; i < dim ; i++) {
    int   in ;
    
    x_max[i] = (x_min[i] = x_e[0][i]) ;
    
    for(in = 1 ; in < nn ; in++) {
      x_max[i] = (x_e[in][i] > x_max[i]) ? x_e[in][i] : x_max[i] ;
      x_min[i] = (x_e[in][i] < x_min[i]) ? x_e[in][i] : x_min[i] ;
    }
    
    x_max[i] -= x_min[i] ;
    d = (x_max[i] > d) ? x_max[i] : d ;
  }


  /* calcul des coordonnees dans la configuration de reference */
  for(iter = 0 ; iter < max_iter ; iter++) {
    IntFct_ComputeIsoShapeFct(dim_e,nn,a,h,dh) ;
    
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
        
        normale(dim,nn,x_e,dh,unorm) ;
        
        for(i = 0 ; i < dim ; i++) k[dim*i + j] = unorm[i] ;
      } else {
        arret("IntFct_ComputeIsoShapeFctInActualSpace (1)") ; /* a traiter plus tard */
      }
    }
    
    for(i = 0 ; i < dim ; i++) {
      int  j ;
      
      r[i] /= d ;
      
      for(j = 0 ; j < dim ; j++) k[dim*i+j] /= d ;
    }
    
    Math_SolveByGaussElimination(k,r,dim) ;
    
    for(i = 0 ; i < dim ; i++) a[i] += r[i] ;
    
    err = 0 ;
    for(i = 0 ; i < dim ; i++) err = (fabs(r[i]) > err) ? fabs(r[i]) : err ;
    /* err /= d ; */
    if(err < tol) return ;
  }
  
  if(err > tol) arret("IntFct_ComputeIsoShapeFctInActualSpace (2)") ;
#undef DH
}




void normale(int dim,int nn,double** x,double* dh,double* norm)
/* Normale unitaire a un sous-espace de dimension dim-1 */
{
#define DH(n,i) (dh[(n)*dim_h+(i)])
  int    dim_h = dim - 1 ;
  int    i,j,k ;
  double c[3][3],v ;

  /* le jacobien */
  for(i=0;i<dim;i++) for(j=0;j<dim_h;j++) {
    c[i][j] = 0. ;
    for(k=0;k<nn;k++) c[i][j] += x[k][i]*DH(k,j) ; 
  }
  
  /* 1. Surface : norm = c1^c2 */
  if(dim_h == 2) {
    c[0][2] = c[1][0]*c[2][1] - c[2][0]*c[1][1] ;
    c[1][2] = c[2][0]*c[0][1] - c[0][0]*c[2][1] ;
    c[2][2] = c[0][0]*c[1][1] - c[1][0]*c[0][1] ;
    v = sqrt(c[0][2]*c[0][2] + c[1][2]*c[1][2] + c[2][2]*c[2][2]) ;
    norm[0] = c[0][2]/v ;
    norm[1] = c[1][2]/v ;
    norm[2] = c[2][2]/v ;
    return ;

  /* 2. Ligne : norm = c1^e_z */
  } else if(dim_h == 1) {
    v = sqrt(c[0][0]*c[0][0] + c[1][0]*c[1][0]) ;
    norm[0] =  c[1][0]/v ;
    norm[1] = -c[0][0]/v ;
    return ;

  /* 3. Point : norm = 1 */
  } else if(dim_h == 0) {
    norm[0] = 1. ;
    return ;
  }

  arret("normale") ;
#undef DH
}

