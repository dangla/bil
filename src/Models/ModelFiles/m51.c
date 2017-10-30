#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  51
#define TITLE "Thermique avec changement de phase liquide-glace"
#define AUTHORS "FenChong-Dangla"

#include "OldMethods.h"

/* Fonctions */
static int    pm(const char *s) ;
/* Parametres */
static double phi0,rho_l,rho_g,C_sk,C_pl,C_pg,T_0,k_sk,k_l,k_g,L_0 ;

int pm(const char *s)
{
if(strcmp(s,"phi") == 0) return (0) ;
else if(strcmp(s,"rho_l") == 0) return (1) ;
else if(strcmp(s,"rho_g") == 0) return (2) ;
else if(strcmp(s,"C_sk") == 0) return (3) ;
else if(strcmp(s,"C_pl") == 0) return (4) ;
else if(strcmp(s,"C_pg") == 0) return (5) ;
else if(strcmp(s,"T_0") == 0) return (6) ;
else if(strcmp(s,"k_sk") == 0) return (7) ;
else if(strcmp(s,"k_l") == 0) return (8) ;
else if(strcmp(s,"k_g") == 0) return (9) ;
else if(strcmp(s,"L_0") == 0) return (10) ;
else if(strcmp(s,"n_p") == 0) return (11) ;
else if(strcmp(s,"T_1") == 0) return (12) ;
else if(strcmp(s,"T_2") == 0) return (13) ;
else if(strcmp(s,"courbe") == 0) return (14) ;
else
  { printf("donnee \"%s\" non connue (pm5)\n",s) ; exit(0) ; }
}

int dm51(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 15 ;
  
  mat->neq    = 1 ;
  strcpy(mat->eqn[0],"chaleur") ;
  strcpy(mat->inc[0],"tem") ;

  dmat(mat,ficd,pm,n_donnees) ;
  
  return(mat->n) ;
}

int qm51(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(1) ;

  printf("\n\n\
Une equation de la chaleur (tem)\n") ;
  
  printf("\n\n\
Exemple de donnees\n\n") ;
  
  fprintf(ficd,"phi = 0.1         # Porosite\n") ;
  fprintf(ficd,"rho_l = 0.998     # Masse volumique du liquide\n") ;
  fprintf(ficd,"rho_g = 0.91      # Masse volumique de la glace\n") ;
  fprintf(ficd,"C_sk = 2e+06      # Chaleur volumique du squeltte\n") ;
  fprintf(ficd,"C_pl = 4180       # Chaleur specifique du liquide\n") ;
  fprintf(ficd,"C_pg = 2000       # Chaleur specifique de la glace\n") ;
  fprintf(ficd,"T_0 = 273         # Temperature\n") ;
  fprintf(ficd,"k_sk = 1          # Conductivite thermique du squelette\n") ;
  fprintf(ficd,"k_l = 0.598       # Conductivite thermique du liquide\n") ;
  fprintf(ficd,"k_g = 2.2         # Conductivite thermique de la glace\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : T S_l\n") ;
  
  return(1) ;
}

void tb51(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 3 ;
  el->n_ve = 1 ;
}

void ch51(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  r[0] = - r[0] ;
}


void in51(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define S(i)   (f[(i)])
#define Q      (f[(2)])
#define Kth    (va[(0)])
  double phi, sl, sg, m0, m_l, m_g, k_e ;
  double dx, Tij ;
  int    i ;
  double un = 1., deux = 2. ;
  
  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  phi0   = el.mat->pr[pm("phi")] ;
  rho_l  = el.mat->pr[pm("rho_l")] ;
  rho_g  = el.mat->pr[pm("rho_g")] ;
  C_sk   = el.mat->pr[pm("C_sk")] ;
  C_pl   = el.mat->pr[pm("C_pl")] ;
  C_pg   = el.mat->pr[pm("C_pg")] ;
  T_0    = el.mat->pr[pm("T_0")] ;
  L_0    = el.mat->pr[pm("L_0")] ;
  k_sk   = el.mat->pr[pm("k_sk")] ;
  k_l    = el.mat->pr[pm("k_l")] ;
  k_g    = el.mat->pr[pm("k_g")] ;
  
  sl     = courbe(T_0,el.mat->cb[0]) ;
  sg     = un - sl ;
  m0     = phi0*(rho_l*sl+rho_g*sg) ;
  
  dx = x[1][0] - x[0][0] ;
  /* ENTROPIE */
  for(i=0;i<2;i++) {
    sl = courbe(u[i][0],el.mat->cb[0]) ;
    sg = un - sl ;
    phi = m0/(rho_l*sl+rho_g*sg) ;
    m_l = rho_l*phi*sl ;
    m_g = rho_g*phi*sg ;
    S(i) = (C_sk+m_l*C_pl+m_g*C_pg)*log(u[i][0]/T_0)+m_l*L_0/T_0 ;
  }
  /* COEFFICIENTS DE TRANSFERT */
  Tij = (u[0][0]+u[1][0])/deux ;
  sl  = courbe(Tij,el.mat->cb[0]) ;
  sg  = un - sl ;
  phi = m0/(rho_l*sl+rho_g*sg) ;
  k_e = sl*k_l+sg*k_g ;
  Kth = k_sk*(un - 3.*phi*(k_sk-k_e)/(3.*k_sk-(un-phi)*(k_sk-k_e)))/Tij ;
  /* FLUX */
  Q   = - Kth*(u[1][0] - u[0][0])/dx ;

#undef S
#undef Q
#undef Kth
}

int ex51(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define Kth      (va[(0)])
  double phi, sl, sg, m0, k_e ;
  double Tij ;
  double un = 1., deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi0   = el.mat->pr[pm("phi")] ;
  rho_l  = el.mat->pr[pm("rho_l")] ;
  rho_g  = el.mat->pr[pm("rho_g")] ;
  C_sk   = el.mat->pr[pm("C_sk")] ;
  C_pl   = el.mat->pr[pm("C_pl")] ;
  C_pg   = el.mat->pr[pm("C_pg")] ;
  T_0    = el.mat->pr[pm("T_0")] ;
  k_sk   = el.mat->pr[pm("k_sk")] ;
  k_l    = el.mat->pr[pm("k_l")] ;
  k_g    = el.mat->pr[pm("k_g")] ;
  /* COEFFICIENTS DE TRANSFERT */
  sl  = courbe(T_0,el.mat->cb[0]) ;
  sg  = un - sl ;
  m0  = phi0*(rho_l*sl+rho_g*sg) ;
  Tij = (u[0][0]+u[1][0])/deux ;
  sl  = courbe(Tij,el.mat->cb[0]) ;
  sg  = un - sl ;
  phi = m0/(rho_l*sl+rho_g*sg) ;
  k_e = sl*k_l+sg*k_g ;
  Kth = k_sk*(un - 3.*phi*(k_sk-k_e)/(3.*k_sk-(un-phi)*(k_sk-k_e)))/Tij ;
  return(0) ;
  
#undef Kth
}

int ct51(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define S_1(i)   (f_1[(i)])
#define Q_1      (f_1[(2)])
#define Kth      (va[(0)])
  double phi, sl, sg, m0, m_l, m_g ;
  double dx ;
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi0   = el.mat->pr[pm("phi")] ;
  rho_l  = el.mat->pr[pm("rho_l")] ;
  rho_g  = el.mat->pr[pm("rho_g")] ;
  C_sk   = el.mat->pr[pm("C_sk")] ;
  C_pl   = el.mat->pr[pm("C_pl")] ;
  C_pg   = el.mat->pr[pm("C_pg")] ;
  T_0    = el.mat->pr[pm("T_0")] ;
  L_0    = el.mat->pr[pm("L_0")] ;
  k_sk   = el.mat->pr[pm("k_sk")] ;
  k_l    = el.mat->pr[pm("k_l")] ;
  k_g    = el.mat->pr[pm("k_g")] ;
  
  sl     = courbe(T_0,el.mat->cb[0]) ;
  sg     = un - sl ;
  m0     = phi0*(rho_l*sl+rho_g*sg) ;
  
  dx = x[1][0] - x[0][0] ;
  /* ENTROPIE */
  for(i=0;i<2;i++) {
    sl  = courbe(u_1[i][0],el.mat->cb[0]) ;
    sg  = un - sl ;
    phi = m0/(rho_l*sl+rho_g*sg) ;
    m_l = rho_l*phi*sl ;
    m_g = rho_g*phi*sg ;
    S_1(i) = (C_sk+m_l*C_pl+m_g*C_pg)*log(u_1[i][0]/T_0)+m_l*L_0/T_0 ;
  }
  /* FLUX */
  Q_1 = - Kth*(u_1[1][0] - u_1[0][0])/dx ;

  return(0) ;
  
#undef S_1
#undef Q_1
#undef Kth
}

int  mx51(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)   (k[(i)*2+(j)])
#define Kth      (va[(0)])
  double phi, sl, sg, m0, m_l, m_g, dslsdT, trth, dm_l ;
  double dx, xm ;
  double volume[2], surf ;
  int    i ;
  double zero = 0., un = 1., deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn;i++) k[i] = zero ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi0   = el.mat->pr[pm("phi")] ;
  rho_l  = el.mat->pr[pm("rho_l")] ;
  rho_g  = el.mat->pr[pm("rho_g")] ;
  C_sk   = el.mat->pr[pm("C_sk")] ;
  C_pl   = el.mat->pr[pm("C_pl")] ;
  C_pg   = el.mat->pr[pm("C_pg")] ;
  T_0    = el.mat->pr[pm("T_0")] ;
  L_0    = el.mat->pr[pm("L_0")] ;
  k_sk   = el.mat->pr[pm("k_sk")] ;
  k_l    = el.mat->pr[pm("k_l")] ;
  k_g    = el.mat->pr[pm("k_g")] ;
  
  sl     = courbe(T_0,el.mat->cb[0]) ;
  sg     = un - sl ;
  m0     = phi0*(rho_l*sl+rho_g*sg) ;
  
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
    EQUATION DE LA CHALEUR : (S_1 - S_n) + dt * div(q/T) = 0
  */
  /* termes d'accumulation */
  for(i=0;i<2;i++)
    {
      sl  = courbe(u_1[i][0],el.mat->cb[0]) ;
      sg  = un - sl ;
      dslsdT = dcourbe(u_1[i][0],el.mat->cb[0]) ;
      phi = m0/(rho_l*sl+rho_g*sg) ;
      m_l = rho_l*phi*sl ;
      m_g = rho_g*phi*sg ;
      dm_l = rho_l*rho_g*phi*phi/m0*dslsdT ;
      K(i,i) += volume[i]*((C_sk+m_l*C_pl+m_g*C_pg)/u_1[i][0]
			   + ((C_pl-C_pg)*log(u_1[i][0]/T_0)+L_0/T_0)*dm_l) ;
  }
  /*
    termes d'ecoulement
  */
  trth    = dt*surf/dx*Kth ;
  K(0,0) += + trth ;
  K(0,1) += - trth ;
  K(1,0) += - trth ;
  K(1,1) += + trth ;

  return(0) ;
  
#undef K
#undef Kth
}

void rs51(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define S_1(i)   (f_1[(i)])
#define Q_1      (f_1[(2)])
#define S_n(i)   (f_n[(i)])
  double dx, xm ;
  double volume[2], surf ;
  int    i ;
  double zero = 0., un = 1., deux = 2. ;
  
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
    EQUATION DE LA CHALEUR : (S_1 - S_n) + dt * div(q/T) = 0
  */
  r[0] -= volume[0]*(S_1(0) - S_n(0)) + dt*surf*Q_1 ;
  r[1] -= volume[1]*(S_1(1) - S_n(1)) - dt*surf*Q_1 ;
  
#undef S_1
#undef Q_1
#undef S_n
}

int so51(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double sl,sg,m0,T,phi ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0., un = 1. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi0   = el.mat->pr[pm("phi")] ;
  rho_l  = el.mat->pr[pm("rho_l")] ;
  rho_g  = el.mat->pr[pm("rho_g")] ;
  C_sk   = el.mat->pr[pm("C_sk")] ;
  C_pl   = el.mat->pr[pm("C_pl")] ;
  C_pg   = el.mat->pr[pm("C_pg")] ;
  T_0    = el.mat->pr[pm("T_0")] ;
  k_sk   = el.mat->pr[pm("k_sk")] ;
  k_l    = el.mat->pr[pm("k_l")] ;
  k_g    = el.mat->pr[pm("k_g")] ;

  /* initialisation */
  nso = 3 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  /* Temperature */
  T   =  param(u,h_s,el.nn,0) ;
  /* Saturations */
  sl  = courbe(T_0,el.mat->cb[0]) ;
  sg  = un - sl ;
  m0  = rho_l*sl+rho_g*sg ;
  sl  = courbe(T,el.mat->cb[0]) ;
  sg  = un - sl ;
  /* Porosite */
  phi = phi0*m0/(rho_l*sl+rho_g*sg) ;
  /* quantites exploitees */
  strcpy(r[0].text,"T") ; r[0].n = 1 ;
  r[0].v[0] = T ;
  strcpy(r[1].text,"porosite") ; r[1].n = 1 ;
  r[1].v[0] = phi ;
  strcpy(r[2].text,"s_g") ; r[2].n = 1 ;
  r[2].v[0] = sg ;
  return(nso) ;
}
