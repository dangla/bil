#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  20
#define TITLE "Equation de Richards (1D, inconnue S_l)"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ   (1)
#define E_liq (0)
#define I_s_l (0)
/* Fonctions */
static int    pm(const char *s) ;
/* Parametres */
static double gravite,phi,rho_l,k_int,mu_l,p_g,schema ;

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"p_g") == 0) return (5) ;
  else if(strcmp(s,"schema") == 0) return (6) ;
  else if(strcmp(s,"courbes") == 0) return (7) ;
  else return(-1) ;
}

int dm20(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 8 ;

  if(dim > 1) arret("dm20 : dimension > 1 non prevue") ;

  mat->neq        = NEQ ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->inc[I_s_l],"s_l") ;

  mat->pr[pm("schema")] = -1 ;
  dmat(mat,ficd,pm,n_donnees) ;
  return(n_donnees) ;
}

int qm20(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  fprintf(stdout,"\n\n\
L\'equation de conservation de la masse d\'eau, inconnue (s_l)\n") ;

  fprintf(stdout,"\n\
Exemple de donnees :\n\n") ;

  fprintf(ficd,"gravite = -9.81  # La gravite\n") ;
  fprintf(ficd,"phi = 0.38       # La porosite\n") ;
  fprintf(ficd,"rho_l = 1000     # La masse volumique du fluide\n") ;
  fprintf(ficd,"k_int = 8.9e-12  # La permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001     # La viscosite du fluide\n") ;
  fprintf(ficd,"p_g = 1.e5       # La pression du gaz\n") ;
  fprintf(ficd,"courbes = billes # Le nom du fichier p_c S_l k_rl\n") ;
  fprintf(ficd,"schema = 0       # schema centre, 1 = decentre en amont\n") ;

  return(1) ;
}

void tb20(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 3 ;
  el->n_ve = 2 ;
}


void ch20(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  r[0] = - r[0] ;
}

void in20(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define M_l(i)   (f[(i)])
#define W_l      (f[(2)])
#define K_l      (va[(0)])
#define KK_l     (va[(1)])
  double dx,sl ;
  int    i ;
  double deux = 2. ;
  
  if(el.dim < dim) return ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;

  dx = x[1][0] - x[0][0] ;
  /* MASSE LIQUIDE */
  for(i=0;i<2;i++) {
    M_l(i) = rho_l*phi*u[i][0] ;
  }
  /* COEFFICIENTS DE TRANSFERT */
  sl   = (u[0][0] + u[1][0])/deux ;
  K_l  = rho_l*k_int/mu_l*courbe(sl,el.mat->cb[1]) ;
  KK_l = - K_l*dcourbe(sl,el.mat->cb[0]) ;
  /* FLUX */
  W_l = - KK_l*(u[1][0] - u[0][0])/dx + K_l*rho_l*gravite ;

#undef M_l
#undef W_l
#undef K_l
#undef KK_l
}

int ex20(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define W_l      (f[(2)])
#define K_l      (va[(0)])
#define KK_l     (va[(1)])
  double sl,pl,z,h_0,h_1 ;
  double deux = 2. ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  schema  = el.mat->pr[pm("schema")] ;
  /*
    COEFFICIENTS DE TRANSFERT
  */
  sl  = (u[0][0] + u[1][0])/deux ;
  if(schema == 1) {
    z    = x[0][0] ;
    sl   = u[0][0] ;
    pl   = p_g - courbe(sl,el.mat->cb[0]) ;
    h_0  = pl - rho_l*gravite*z ;


    z    = x[1][0] ;
    sl   = u[1][0] ;
    pl   = p_g - courbe(sl,el.mat->cb[0]) ;
    h_1  = pl - rho_l*gravite*z ;

    sl   = (h_0 > h_1) ? u[0][0] : u[1][0] ;
  }
  K_l  = rho_l*k_int/mu_l*courbe(sl,el.mat->cb[1]) ;
  KK_l = - K_l*dcourbe(sl,el.mat->cb[0]) ;
  return(0) ;

#undef W_l
#undef K_l
#undef KK_l
}

int ct20(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define M_l1(i)   (f_1[(i)])
#define W_l1      (f_1[(2)])
#define K_l       (va[(0)])
#define KK_l      (va[(1)])
  double dx ;
  int    i ;
  
  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  /* masse de fluide */
  for(i=0;i<2;i++) {
    M_l1(i) = rho_l*phi*u_1[i][0] ;
  }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  W_l1 = - KK_l*(u_1[1][0] - u_1[0][0])/dx + K_l*rho_l*gravite ;

  return(0) ;
  
#undef M_l1
#undef W_l1
#undef K_l
}

int mx20(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2+(j)])
#define K_l       (va[(0)])
#define KK_l      (va[(1)])
  double dx,xm,tr_l ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    K(i,i) += volume[i]*rho_l*phi ;
  }
  /* termes d'ecoulement */
  tr_l = dt*surf*KK_l/dx ;
  K(0,0) += + tr_l ;
  K(0,1) += - tr_l ;
  K(1,0) += - tr_l ;
  K(1,1) += + tr_l ;

  return(0) ;

#undef K
#undef K_l
#undef KK_l
}

void rs20(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define M_l1(i)   (f_1[(i)])
#define W_l1      (f_1[(2)])
#define M_ln(i)   (f_n[(i)])
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;

  /* initialisation */
  for(i=0;i<el.nn;i++) r[i] = zero ;

  if(el.dim < dim) return ;
  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;

  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  r[0] -= volume[0]*(M_l1(0) - M_ln(0)) + dt*surf*W_l1 ;
  r[1] -= volume[1]*(M_l1(1) - M_ln(1)) - dt*surf*W_l1 ;

#undef M_l1
#undef W_l1
#undef M_ln
}

int so20(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define K_l     (va[(0)])
#define KK_l    (va[(1)])
  double dx,pl,sl,wl;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;

  /* initialisation */
  nso = 3 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* saturation */
  sl =  param(u,h_s,el.nn,0) ;
  /* pression */
  pl = p_g - courbe(sl,el.mat->cb[0]) ;
  /* flux */
  dx = x[1][0] - x[0][0] ;
  wl = - KK_l*(u[1][0] - u[0][0])/dx + K_l*rho_l*gravite ;

  /* quantites exploitees */
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"flux-liquide") ; r[1].n = 1 ;
  r[1].v[0] = wl ;
  strcpy(r[2].text,"saturation") ; r[2].n = 1 ;
  r[2].v[0] = sl ;
  return(nso) ;

#undef K_l
#undef KK_l
}
