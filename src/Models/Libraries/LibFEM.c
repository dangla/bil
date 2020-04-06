#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include "LibFEM.h"
#include "Message.h"
#include "Tools/Math.h"
#include "Geometry.h"
#include "Elements.h"
#include "IntFcts.h"
#include "Nodes.h"

/*
   Fonctions locales 
*/
static double det(int,int,double **,double *,int) ;
static double jac(int,int,double **,double *,double *) ;
static double*  ComputeNormalVector(Element_t*,double*,int) ;
static void   rssurfloads(Element_t*,IntFct_t*,double,double,Load_t*,double*) ;


/* 
   Fonctions globales 
*/
void mxcpl(double *kc,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym)
/* Matrice de couplage mecanique (kc) */
{
#define NN          (fi.nn)
#define KC(i,j)     (kc[(i)*NN+(j)])
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*fi.dim+(i)])
#define C(i,j)      (c[(i)*3+(j)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
  int    i,j,k,l,p ;
  double a,rayon = 0. ;
  double caj[9],jc[3][3] ;
  double *h,*dh ;
  double zero = 0., deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN*dim;i++) for(j=0;j<NN;j++) KC(i,j) = zero ;
  /* boucle sur les points d'integration */
  c -= dec ;
  for(p=0;p<fi.np;p++) {
    c += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*jac(dim,NN,x,dh,caj) ;
    /* cas axisymetrique ou spherique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* JC(k,i) = J(k,j)*C(j,i) */
    for(k=0;k<dim;k++) for(i=0;i<dim;i++) {
      jc[k][i] = zero ;
      for(j=0;j<dim;j++) jc[k][i] += CAJ(k,j)*C(j,i) ;
    }
    /* KC(j,i,l) = DH(j,k)*JC(k,i)*H(l) */
    for(j=0;j<NN;j++) for(i=0;i<dim;i++) for(l=0;l<NN;l++) for(k=0;k<dim;k++) {
      KC(j*dim+i,l) += a*DH(j,k)*jc[k][i]*h[l] ;
    }
    /* cas axisymetrique */
    if(sym == AXIS) {
      /* KC(j,0,l) = H(j)/r*C(theta,theta)*H(l) */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) {
        KC(j*dim,l) += a*h[j]/rayon*C(2,2)*h[l] ;
      }
    /* cas spherique */
    } else if(sym == SPHE) {
      /* KC(j,0,l) = H(j)/r*(C(theta,theta)+C(phi,phi))*H(l) */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) {
        KC(j*dim,l) += a*h[j]/rayon*(C(1,1)+C(2,2))*h[l] ;
      }
    }
  }
#undef NN
#undef KC
#undef DH
#undef C
#undef CAJ
}

void mxmass(double *km,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym)
/* Matrice de masse (km) */
{
#define NN          (fi.nn)
#define KM(i,j)     (km[(i)*NN+(j)])
  int    i,j,p ;
  double a,rayon ;
  double *h,*dh ;
  double zero = 0., deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN;i++) for(j=0;j<NN;j++) KM(i,j) = zero ;

  /* 0D */
  if(fi.dim == 0) {
    if(NN == 1) {
      KM(0,0) = c[0] ;
      if(sym == AXIS) KM(0,0) *= deux*M_PI*x[0][0] ;
      else if(sym == SPHE) KM(0,0) *= 4*M_PI*x[0][0]*x[0][0] ;
    } else arret("mxmass : impossible") ;
    return ;
  }

  /* boucle sur les points d'integration */
  c -= dec ;
  for(p=0;p<fi.np;p++) {
    c += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*c[0]*det(dim,NN,x,dh,fi.dim) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    for(i=0;i<NN;i++) for(j=0;j<NN;j++) KM(i,j) += a*h[i]*h[j] ;
  }
#undef NN
#undef KM
}



void mxcond(double *kc,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym)
/* Matrice de conduction (kc) */
{
#define NN          (fi.nn)
#define KC(i,j)     (kc[(i)*NN+(j)])
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*fi.dim+(i)])
#define C(i,j)      (c[(i)*3+(j)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
  int    i,j,k,l,p ;
  double a,rayon ;
  double caj[9],jcj[3][3] ;
  double *h,*dh ;
  double zero = 0., deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN;i++) for(j=0;j<NN;j++) KC(i,j) = zero ;

  /* boucle sur les points d'integration */
  c -= dec ;
  for(p=0;p<fi.np;p++) {
    c += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*jac(dim,NN,x,dh,caj) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* jcj = J(i,k)*C(k,l)*J(j,l) */
    for(i=0;i<dim;i++) for(j=0;j<dim;j++) {
      jcj[i][j] = zero ;
      for(k=0;k<dim;k++) for(l=0;l<dim;l++)  {
        jcj[i][j] += CAJ(i,k)*C(k,l)*CAJ(j,l) ;
      }
    }
    /* KC(i,j) = DH(i,k)*JCJ(k,l)*DH(j,l) */
    for(i=0;i<NN;i++) for(j=0;j<NN;j++) for(k=0;k<dim;k++) for(l=0;l<dim;l++) {
      KC(i,j) += a*DH(i,k)*jcj[k][l]*DH(j,l) ;
    }
  }
#undef NN
#undef KC
#undef DH
#undef C
#undef CAJ
}






void  mxrig(double *kr,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym)
/* Matrice de rigidite (kr) */
{
#define NN          (fi.nn)
#define KR(i,j)     (kr[(i)*NN*dim+(j)])
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*fi.dim+(i)])
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
  int    i,j,k,l,r,s,p ;
  double a,rayon = 0. ;
  double caj[9],jcj[3][3][3][3], jc[3][3], cj[3][3] ;
  double *h,*dh ;
  double zero = 0., deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN*dim;i++) for(j=0;j<NN*dim;j++) KR(i,j) = zero ;
  /* boucle sur les points d'integration */
  c -= dec ;
  for(p=0;p<fi.np;p++) {
    c += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*jac(dim,NN,x,dh,caj) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* JCJ(r,i,k,s) = J(r,j)*C(i,j,k,l)*J(s,l) */
    for(r=0;r<dim;r++) for(i=0;i<dim;i++) for(k=0;k<dim;k++) for(s=0;s<dim;s++) {
      jcj[r][i][k][s] = zero ;
      for(j=0;j<dim;j++) for(l=0;l<dim;l++) {
        jcj[r][i][k][s] += CAJ(r,j)*C(i,j,k,l)*CAJ(s,l) ;
      }
    }
    /* KR(j,i,l,k) = DH(j,r)*JCJ(r,i,k,s)*DH(l,s) */
    for(j=0;j<NN;j++) for(i=0;i<dim;i++) for(l=0;l<NN;l++) for(k=0;k<dim;k++) {
      for(r=0;r<dim;r++) for(s=0;s<dim;s++) {
        KR(j*dim+i,l*dim+k) += a*DH(j,r)*jcj[r][i][k][s]*DH(l,s) ;
      }
    }
    /* cas axisymetrique: 3 termes */
    if(sym == AXIS) {
      /* 1.a JC(r,i) = J(r,j)*C(i,j,theta,theta) */
      for(r=0;r<dim;r++) for(i=0;i<dim;i++) {
        jc[r][i] = zero ;
        for(j=0;j<dim;j++) jc[r][i] += CAJ(r,j)*C(i,j,2,2) ;
      }
      /* 1.b KR(j,i,l,0) = DH(j,r)*JC(r,i)*H(l)/r */
      for(j=0;j<NN;j++) for(i=0;i<dim;i++) for(l=0;l<NN;l++) for(r=0;r<dim;r++) {
        KR(j*dim+i,l*dim) += a*DH(j,r)*jc[r][i]*h[l]/rayon ;
      }
      /* 2.a CJ(k,s) = C(theta,theta,k,l)*J(s,l) */
      for(k=0;k<dim;k++) for(s=0;s<dim;s++) {
        cj[k][s] = zero ; 
        for(l=0;l<dim;l++) cj[k][s] += C(2,2,k,l)*CAJ(s,l) ;
      }
      /* 2.b KR(j,0,l,k) = H(j)/r*CJ(k,s)*DH(l,s) */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) for(k=0;k<dim;k++) for(s=0;s<dim;s++) {
        KR(j*dim,l*dim+k) += a*h[j]/rayon*cj[k][s]*DH(l,s) ;
      }
      /* 3.  KR(j,0,l,0) = H(j)/r*C(theta,theta,theta,theta)*H(l)/r */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) {
        KR(j*dim,l*dim) += a*h[j]/rayon*C(2,2,2,2)*h[l]/rayon ;
      }
    /* cas spherique: 3 termes */
    } else if(sym == SPHE) {
      /* 1.a JC(r,i) = J(r,j)*(C(i,j,theta,theta)+C(i,j,phi,phi)) */
      for(r=0;r<dim;r++) for(i=0;i<dim;i++) {
        jc[r][i] = zero ;
        for(j=0;j<dim;j++) jc[r][i] += CAJ(r,j)*(C(i,j,2,2)+C(i,j,1,1)) ;
      }
      /* 1.b KR(j,i,l,0) = DH(j,r)*JC(r,i)*H(l)/r */
      for(j=0;j<NN;j++) for(i=0;i<dim;i++) for(l=0;l<NN;l++) for(r=0;r<dim;r++) {
        KR(j*dim+i,l*dim) += a*DH(j,r)*jc[r][i]*h[l]/rayon ;
      }
      /* 2.a CJ(k,s) = (C(phi,phi,k,l)+C(theta,theta,k,l))*J(s,l) */
      for(k=0;k<dim;k++) for(s=0;s<dim;s++) {
        cj[k][s] = zero ; 
        for(l=0;l<dim;l++) cj[k][s] += (C(1,1,k,l)+C(2,2,k,l))*CAJ(s,l) ;
      }
      /* 2.b KR(j,0,l,k) = H(j)/r*CJ(k,s)*DH(l,s) */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) for(k=0;k<dim;k++) for(s=0;s<dim;s++) {
        KR(j*dim,l*dim+k) += a*h[j]/rayon*cj[k][s]*DH(l,s) ;
      }
      /* 3.  KR(j,0,l,0) = H(j)/r*(C(phi,phi,phi,phi)+C(theta,theta,phi,phi)+C(phi,phi,theta,theta)+C(theta,theta,theta,theta))*H(l)/r */
      for(j=0;j<NN;j++) for(l=0;l<NN;l++) {
        KR(j*dim,l*dim) += a*h[j]/rayon*(C(1,1,1,1)+C(1,1,2,2)+C(2,2,1,1)+C(2,2,2,2))*h[l]/rayon ;
      }
    }
  }
#undef NN
#undef KR
#undef DH
#undef C
#undef CAJ
}




void  mxbiot(double *k,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym,int neq,int e_mec,int i_u,int e_dif,int i_p)
/* Matrice de Biot Hydro-Mecanique (k) */
{
#define NN        (fi.nn)
#define K(i,j)    (k[(i)*NN*neq+(j)])
  int    i,j,n,m ;
  double kb[9*Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<NN*NN*neq*neq;i++) k[i] = zero ;

  /* 
  ** 1.  Mecanique
  */
  /* 1.1 rigidite du squelette */
  mxrig(kb,x,fi,c,dim,dec,sym) ;
  for(n=0;n<NN;n++) for(i=0;i<dim;i++) for(m=0;m<NN;m++) for(j=0;j<dim;j++) {
    K(e_mec+i+n*neq,i_u+j+m*neq) = kb[(i+n*dim)*NN*dim+j+m*dim] ;
  }
  c += 81 ;
  /* 1.2 couplage mecanique */
  mxcpl(kb,x,fi,c,dim,dec,sym) ;
  for(n=0;n<NN;n++) for(i=0;i<dim;i++) for(m=0;m<NN;m++) {
    K(e_mec+i+n*neq,i_p+m*neq) = kb[(i+n*dim)*NN+m] ;
  }
  c += 9 ;
  /*
  ** 2.  Hydraulique
  */
  /* 2.1 couplage diffusif */
  mxcpl(kb,x,fi,c,dim,dec,sym) ;
  for(n=0;n<NN;n++) for(i=0;i<dim;i++) for(m=0;m<NN;m++) {
    K(e_dif+m*neq,i_u+i+n*neq) = kb[(i+n*dim)*NN+m] ;
  }
  c += 9 ;
  /* 2.2 accumulation */
  mxmass(kb,x,fi,c,dim,dec,sym) ;
  for(n=0;n<NN;n++) for(m=0;m<NN;m++) {
    K(e_dif+n*neq,i_p+m*neq) = kb[n*NN+m] ;
  }
#undef NN
#undef K
}



void mxcmss(double *k,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym,int neq)
/* Matrices de masses couplees (k) */
{
#define NN        (fi.nn)
#define K(i,j)    (k[(i)*NN*neq+(j)])
  int    i,j,n,m ;
  double kb[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<NN*NN*neq*neq;i++) k[i] = zero ;

  for(i=0;i<neq;i++) for(j=0;j<neq;j++) {
    mxmass(kb,x,fi,c,dim,dec,sym) ;
    for(n=0;n<NN;n++) for(m=0;m<NN;m++) {
      K(i+n*neq,j+m*neq) = kb[n*NN+m] ;
    }
    c += 1 ;
  }
#undef NN
#undef K
}

void mxccnd(double *k,double **x,IntFct_t fi,double *c,int dim,int dec,Symmetry_t sym,int neq)
/* Matrices de conduction couplees (k) */
{
#define NN        (fi.nn)
#define K(i,j)    (k[(i)*NN*neq+(j)])
  int    i,j,n,m ;
  double kb[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<NN*NN*neq*neq;i++) k[i] = zero ;

  for(i=0;i<neq;i++) for(j=0;j<neq;j++) {
    mxcond(kb,x,fi,c,dim,dec,sym) ;
    for(n=0;n<NN;n++) for(m=0;m<NN;m++) {
      K(i+n*neq,j+m*neq) = kb[n*NN+m] ;
    }
    c += 9 ;
  }
#undef NN
#undef K
}

void   rscont(double *r,double **x,IntFct_t fi,double *s,int dim,int dec,Symmetry_t sym)
/* Forces dues aux contraintes */
{
#define NN          (fi.nn)
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*fi.dim+(i)])
#define S(i,j)      (s[(i)*3+(j)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
#define R(n,i)      (r[(n)*dim+(i)])
  int    i,j,k,l,p ;
  double a,rayon = 0. ;
  double caj[9] ;
  double *h,*dh ;
  double zero = 0.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN*dim;i++) r[i] = zero ;
  /* boucle sur les points d'integration */
  s -= dec ;
  for(p=0;p<fi.np;p++) {
    s += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*jac(dim,NN,x,dh,caj) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* R(i,j) = DH(i,k)*J(k,l)*S(l,j) */
    for(i=0;i<NN;i++) for(j=0;j<dim;j++) {
      for(k=0;k<dim;k++) for(l=0;l<dim;l++) {
        R(i,j) += a*DH(i,k)*CAJ(k,l)*S(l,j) ;
      }
    }
    /* cas axisymetrique ou spherique */
    if(sym == AXIS) {
      /* R(i,0) = H(i)/r*S(theta,theta) */
      for(i=0;i<NN;i++) R(i,0) += a*h[i]/rayon*S(2,2) ;
    } else if(sym == SPHE) {
      /* R(i,0) = H(i)/r*(S(theta,theta)+S(phi,phi)) */
      for(i=0;i<NN;i++) R(i,0) += a*h[i]/rayon*(S(1,1)+S(2,2)) ;
    }
  }
#undef NN
#undef DH
#undef S
#undef CAJ
#undef R
}

void   rsflux(double *r,double **x,IntFct_t fi,double *f,int dim,int dec,Symmetry_t sym)
/* Forces dues a un flux */
{
#define NN          (fi.nn)
#define DH(n,i)     (dh[(n)*3+(i)])
//#define DH(n,i)     (dh[(n)*fi.dim+(i)])
#define CAJ(i,j)    (caj[(i)*dim+(j)])
  int    i,j,k,p ;
  double a,rayon ;
  double caj[9] ;
  double *h,*dh ;
  double zero = 0.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<NN;i++) r[i] = zero ;
  /* boucle sur les points d'integration */
  f -= dec ;
  for(p=0;p<fi.np;p++) {
    f += dec ;
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*jac(dim,NN,x,dh,caj) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* R(i) = DH(i,k)*J(k,j)*F(j) */
    for(i=0;i<NN;i++) {
      for(j=0;j<dim;j++) for(k=0;k<dim;k++) r[i] += a*DH(i,k)*CAJ(k,j)*f[j] ;
    }
  }
#undef NN
#undef DH
#undef CAJ
}

void   rsmass(double *r,double **x,IntFct_t fi,double *f,int dim,int dec,Symmetry_t sym)
/* Forces de masse */
{
#define NN          (fi.nn)
  int    i,p ;
  double a,rayon ;
  double *h,*dh ;
  double zero = 0.,deux = 2. ;

  /* initialisation */
  for(i=0;i<NN;i++) r[i] = zero ;

  /* 0D */
  if(fi.dim == 0) {
    if(NN == 1) {
      r[0] = f[0] ;
      if(sym == AXIS) r[0] *= deux*M_PI*x[0][0] ;
      else if(sym == SPHE) r[0] *= 4*M_PI*x[0][0]*x[0][0] ;
    } else arret("rsmass : impossible") ;
    return ;
  }

  /* 1D, 2D, 3D */
  for(p=0;p<fi.np;p++) {
    h  = fi.h  + p*NN ;
    dh = fi.dh + p*NN*fi.dim ;
    a  = fi.w[p]*det(dim,NN,x,dh,fi.dim) ;
    /* cas axisymetrique ou shperique */
    if(sym == AXIS || sym == SPHE) {
      rayon = zero ;
      for(i=0;i<NN;i++) rayon += h[i]*x[i][0] ;
      a *= deux*M_PI*rayon ;
      if(sym == SPHE) a *= deux*rayon ;
    }
    /* R(i) = F*H(i) */
    for(i=0;i<NN;i++) r[i] += a*f[p*dec]*h[i] ;
  }
#undef NN
}



void chsurf(double **x,double *r,int dim,Symmetry_t sym,double dt,double t,Load_t cg,Element_t el,IntFct_t *fi)
/* Residu du aux chargements de surface (r) */
{
#define NEQ     Element_GetNbOfEquations(&el)
#define MAT     Element_GetMaterial(&el)
#define NN      Element_GetNbOfNodes(&el)
#define NP      IntFct_GetNbOfPoints(fi)
#define H(p)    IntFct_GetFunctionAtPoint(fi,p)
#define DH(p)   IntFct_GetFunctionGradientAtPoint(fi,p)
  int    ii,jj,ieq ;
  double zero = 0.,ft ;
  int    i ;

  /* initialisation */
  for(i = 0 ; i < NN*NEQ ; i++) r[i] = zero ;

  if(cg.ch == NULL) return ;

  /* Le numero de l'equation */
  ieq = Element_FindEquationPositionIndex(&el,cg.eqn) ;
  if(ieq < 0) arret("chsurf (1) : equation non connue") ;

  /* flux hydraulique */
  if(strncmp(cg.t,"flux",4) == 0) {
    if(cg.fn != NULL) ft = fonction(t,*cg.fn) - fonction(t-dt,*cg.fn) ;
    else ft = dt ;
    if(dim == 1 && NN == 1) {
      r[ieq] = ft*champ(x[0],dim,*cg.ch) ;
      if(sym == AXIS) r[ieq] *= 2.*M_PI*x[0][0] ;
      else if(sym == SPHE) r[ieq] *= 4.*M_PI*x[0][0]*x[0][0] ;
      return ;
    }
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      for(p = 0 ; p < NP ; p++) {
        double *h = H(p) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < NN ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*champ(y,dim,*cg.ch) ;
      }
      rsmass(rb,x,*fi,f,dim,1,sym) ;
      for(i = 0 ; i < NN ; i++) r[i*NEQ+ieq] = rb[i] ;
      return ;
    }
  }
  /* force mecanique */
  if(strncmp(cg.t,"force",5) == 0) {
    if(cg.fn != NULL)  ft = fonction(t,*cg.fn) ;
    else ft = 1. ;
    if(dim == 1 && NN == 1) {
      r[ieq] = ft*champ(x[0],dim,*cg.ch) ;
      if(sym == AXIS) r[ieq] *= 2.*M_PI*x[0][0] ;
      else if(sym == SPHE) r[ieq] *= 4.*M_PI*x[0][0]*x[0][0] ;
      return ;
    }
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      for(p = 0 ; p < NP ; p++) {
        double *h = H(p) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < NN ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*champ(y,dim,*cg.ch) ;
      }
      rsmass(rb,x,*fi,f,dim,1,sym) ;
      for(i = 0 ; i < NN ; i++) r[i*NEQ+ieq] = rb[i] ;
      return ;
    }
  }
  /* pression mecanique */
  if(strncmp(cg.t,"press",5) == 0) {
    if(cg.fn != NULL)  ft = fonction(t,*cg.fn) ;
    else ft = 1. ;
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      for(p = 0 ; p < NP ; p++) {
        double *h = H(p) ;
        double *dh = DH(p) ;
        double *n = ComputeNormalVector(&el,dh,NN) ;
        double y[3] = {0.,0.,0.,} ;
        double f0 ;
        for(i = 0 ; i < dim ; i++) {
          int j ;
          for(j = 0 ; j < NN ; j++) y[i] += h[j]*x[j][i] ;
        }
        f0 = ft*champ(y,dim,*cg.ch) ;
        for(i = 0 ; i < dim ; i++) f[p*dim+i] = -f0*n[i] ;
      }
      for(i = 0 ; i < dim ; i++) {
        int j ;
        rsmass(rb,x,*fi,f+i,dim,dim,sym) ;
        for(j = 0 ; j < NN ; j++) r[j*NEQ+ieq+i] = rb[j] ;
      }
      return ;
    }
  }
  /* composante de tenseur */
  if(strncmp(cg.t,"sig_",4) == 0 && (ii=cg.t[4]-'1') >= 0 && ii < dim \
  && (jj=cg.t[5]-'1') >= 0 && jj < dim) {
    if(cg.fn != NULL)  ft = fonction(t,*cg.fn) ;
    else ft = 1. ;
    if(dim >= 2) {
      double rb[Element_MaxNbOfNodes] ;
      double f[3*IntFct_MaxNbOfIntPoints] ;
      int p ;
      int j ;
      for(p = 0 ; p < NP ; p++) {
        double *h = H(p) ;
        double *dh = DH(p) ;
        double *n = ComputeNormalVector(&el,dh,NN) ;
        double y[3] = {0.,0.,0.,} ;
        for(i = 0 ; i < dim ; i++) {
          for(j = 0 ; j < NN ; j++) y[i] += h[j]*x[j][i] ;
        }
        f[p] = ft*champ(y,dim,*cg.ch)*n[jj] ;
      }
      rsmass(rb,x,*fi,f,dim,1,sym) ;
      for(j = 0 ; j < NN ; j++) r[j*NEQ+ieq+ii] = rb[j] ;
      return ;
    }
  }
  arret("chsurf (2) : chargement non prevu") ;
#undef NEQ
#undef NN
#undef MAT
#undef NP
#undef H
#undef DH
}


void rssurfloads(Element_t *el,IntFct_t *intfct,double t,double dt,Load_t *load,double *rs)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned short int dim = Geometry_GetDimension(geom) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  double *xx[Element_MaxNbOfNodes] ;
  int i ;
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    xx[i] = Node_GetCoordinate(no[i]) ;
  }
  chsurf(xx,rs,dim,sym,dt,t,*load,*el,intfct) ;
}


double param(double **u,double *h,int nn,int inc)
/* Calcul du parametre */
{
#define U(n)   (u[(n)][(inc)])
  int    i ;
  double par ;
  double zero = 0. ;
  
  par = zero ;
  for(i=0;i<nn;i++) par += h[i]*U(i) ;
  return (par) ;
  
#undef U
}

void def(double **x,double **u,double *h,double *dh,double *eps,int nn,int dim,Symmetry_t sym,int inc)
/* Tenseur des deformations (eps) */
{
#define U(n,i)   (u[(n)][inc+(i)])
#define DH(n,i)  (dh[(n)*3+(i)])
//#define DH(n,i)  (dh[(n)*dim+(i)])
#define EPS(i,j) (eps[(i)*3+(j)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  int    i,j,k,l ;
  double rayon,u_r ;
  double gu[3][3], cj[9] ;
  double zero = 0., deux = 2. ;
  
  /* initialisation de gu */
  for(i=0;i<3;i++) for(j=0;j<3;j++)  gu[i][j] = zero ;
  /* inverse de la matrice jacobienne (cj) */
  jac(dim,nn,x,dh,cj) ;
  /* gradient des deplacements (gu) */
  for(i=0;i<dim;i++) for(j=0;j<dim;j++) for(k=0;k<nn;k++) for(l=0;l<dim;l++) 
    gu[i][j] += U(k,i)*DH(k,l)*CJ(l,j) ;
  /* petites deformations */
  for(i=0;i<3;i++) for(j=0;j<3;j++) EPS(i,j) = (gu[i][j] + gu[j][i])/deux ;
  /* cas axisymetrique ou shperique */
  if(sym == AXIS || sym == SPHE) {
    rayon = zero ;
    u_r   = zero ;
    for(i=0;i<nn;i++) {
      rayon += h[i]*x[i][0] ;
      u_r   += h[i]*U(i,0) ;
    }
    EPS(2,2) += u_r/rayon ;
    if(sym == SPHE) EPS(1,1) += u_r/rayon ;
  }
#undef U
#undef DH
#undef EPS
#undef CJ
}

void grad(double **x,double **u,double *dh,double *gu,int nn,int dim,int inc)
/* Gradient du parametre i=inc (gu) */
{
#define U(n)     (u[(n)][(inc)])
#define DH(n,i)  (dh[(n)*3+(i)])
//#define DH(n,i)  (dh[(n)*dim+(i)])
#define CJ(i,j)  (cj[(i)*dim+(j)])
  int    i,k,l ;
  double cj[9] ;
  double zero = 0. ;
  
  /* initialisation de gu */
  for(i=0;i<3;i++)  gu[i] = zero ;
  /* inverse de la matrice jacobienne (cj) */
  jac(dim,nn,x,dh,cj) ;
  /* gradient (gu) */
  for(i=0;i<dim;i++) {
    for(k=0;k<nn;k++) for(l=0;l<dim;l++)
    gu[i] += U(k)*DH(k,l)*CJ(l,i) ;
  }
#undef U
#undef DH
#undef CJ
}



/*
   Fonctions locales 
*/

double det(int dim,int nn,double **x,double *dh,int dim_h)
/* Determinant de la matrice jacobienne */
{
#define DH(n,i) (dh[(n)*3+(i)])
//#define DH(n,i) (dh[(n)*dim_h+(i)])
  int    i,j,k ;
  double c[3][3],dt ;

  /* le jacobien */
  for(i=0;i<dim;i++) for(j=0;j<dim_h;j++) {
    c[i][j] = 0. ;
    for(k=0;k<nn;k++) c[i][j] += x[k][i]*DH(k,j) ; 
  }
  
  /* 
     1. Pour un volume   : c'est le det(c1,c2,c3)
     2. Pour une surface : c'est le det(c1,c2,e3) 
        ou e3 est la normale unitaire a (c1,c2)
     3. Pour une ligne   : c'est le det(c1,e2,e3) = |c1|
        ou (e2,e3) est orthonorme dans le plan normal a c1
     4. Pour un point    : c'est le det(e1,e2,e3) = 1
        ou (e1,e2,e3) est le triedre orthonorme.
  */

  /* 1. Volume : det(c1,c2,c3) */
  if(dim_h == 3) {
    if(dim == 3) {
      dt  = c[0][0]*c[1][1]*c[2][2] - c[0][0]*c[2][1]*c[1][2]
    + c[1][0]*c[2][1]*c[0][2] - c[1][0]*c[0][1]*c[2][2]
    + c[2][0]*c[0][1]*c[1][2] - c[2][0]*c[1][1]*c[0][2] ;
      return(fabs(dt)) ;
    } else arret("det (1)") ;

  /* 2. Surface : det(c1,c2,e3) = |c1^c2| */
  } else if(dim_h == 2) {
    if(dim == 3) {
      c[0][2] = c[1][0]*c[2][1] - c[2][0]*c[1][1] ;
      c[1][2] = c[2][0]*c[0][1] - c[0][0]*c[2][1] ;
      c[2][2] = c[0][0]*c[1][1] - c[1][0]*c[0][1] ;
      return(sqrt(c[0][2]*c[0][2] + c[1][2]*c[1][2] + c[2][2]*c[2][2])) ;
    } else if(dim == 2) {
      c[2][2] = c[0][0]*c[1][1] - c[1][0]*c[0][1] ;
      return(fabs(c[2][2])) ;
    } else arret("det (2)") ;

  /* 3. Ligne : det(c1,e2,e3) = |c1|(e1,e2,e3) = |c1|   */
  } else if(dim_h == 1) {
    dt = 0. ;
    for(i=0;i<dim;i++) dt += c[i][0]*c[i][0] ;
    return(sqrt(dt)) ;

  /* 4. Point : det(e1,e2,e3) = 1 */
  } else if(dim_h == 0) {
    return(1.) ;
  }
  
  arret("det (3)") ;
  return(0.) ;
#undef DH
}


double jac(int dim,int nn,double **x,double *dh,double *cj)
/* Determinant de la matrice jacobienne et son inverse (cj) */
{
#define CJ(i,j)  (cj[(i)*dim+(j)])
#define DH(n,i)  (dh[(n)*3+(i)])
//#define DH(n,i)  (dh[(n)*dim_h+(i)])
  int    dim_h = dim ; /* ne reste plus qu'a mettre dim_h en argument */
  int    i,j,k ;
  double jc[3][3],dt,td,v ;
  
  /* le jacobien */
  for(i=0;i<dim;i++) for(j=0;j<dim_h;j++) {
    jc[i][j] = 0. ;
    for(k=0;k<nn;k++) jc[i][j] += x[k][i]*DH(k,j) ; 
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
      } else arret("jac (1)") ;
      return(dt) ;
    } else arret("jac (2)") ;

  /* 2. Surface */
  } else if(dim_h == 2) {
    if(dim == 3) {
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
      } else arret("jac (3)") ;
      return(dt) ;
    } else if(dim == 2) {
      dt = jc[0][0]*jc[1][1] - jc[1][0]*jc[0][1] ;
      if(dt != 0.) {
        td = 1./dt ;
        CJ(0,0) =  jc[1][1]*td ;
        CJ(0,1) = -jc[0][1]*td ;
        CJ(1,0) = -jc[1][0]*td ;
        CJ(1,1) =  jc[0][0]*td ;
      } else arret("jac (4)") ;
      return(dt) ;
    } else arret("jac (5)") ;

  /* 3. Ligne */
  } else if(dim_h == 1) {
    if(dim == 3) {
      arret("jac (6)") ;
    } else if(dim == 2) {
      /* on calcul la normale, c1^e_z, a la ligne */
      v = sqrt(jc[0][0]*jc[0][0] + jc[1][0]*jc[1][0]) ;
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
      } else arret("jac (7)") ;
      return(dt) ;
    } else if(dim == 1) {
      dt = jc[0][0] ;
      if(dt != 0.) CJ(0,0) = 1./dt ;
      else arret("jac (8)") ;
      return(dt) ;
    } else arret("jac (9)") ;
  }

  arret("jac (10)") ;
  return(0.) ;
#undef CJ
#undef DH
}



double*  ComputeNormalVector(Element_t *element,double *dh,int nn)
/* Normale unitaire a un sous-espace de dimension dim-1 */
{
#define DH(n,i) (dh[(n)*3+(i)])
//#define DH(n,i) (dh[(n)*dim_h+(i)])
  size_t SizeNeeded = 3*sizeof(double) ;
  double *norm = (double*) Element_AllocateInBuffer(element,SizeNeeded) ;
  int    dim_h = Element_GetDimension(element) ;
  int    dim   = dim_h + 1 ;
  int    i,j ;
  double c[3][3] ;

  if(dim > 3) arret("ComputeNormalVector") ;
  
  /* le jacobien */
  for(i = 0 ; i < dim ; i++) for(j = 0 ; j < dim_h ; j++) {
    int    k ;
    c[i][j] = 0. ;
    for(k = 0 ; k < nn ; k++) {
      double *x = Element_GetNodeCoordinate(element,k) ;
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
    return(norm) ;

  /* 3. Point : norm = 1 */
  } else if(dim_h == 0) {
    norm[0] = 1. ;
    return(norm) ;
  }

  arret("ComputeNormalVector") ;
  return(norm) ;
#undef DH
}
