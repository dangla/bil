#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Message.h"
#include "Tools/Math.h"
#include "ShapeFct.h"


/* Extern functions */

void ShapeFct_ComputeValuesAtPoint(int dim,int nn,double* x,double* h,double* dh)
/* Compute shape functions (h) and their gradients (dh) at point x */
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
  
  arret("ShapeFct_ComputeValuesAtPoint") ;
  
  return ;
#undef X
#undef Y
#undef Z
#undef DH
#undef DHx
#undef DHy
#undef DHz
}
