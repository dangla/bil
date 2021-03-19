#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  14
#define TITLE "Ecoulement diphasique H2O-H2-air"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ      (3)

#define E_eau    (0)
#define E_hyd    (1)
#define E_air    (2)

#define I_p_l    (0)
#define I_p_h    (1)
#define I_p_a    (2)

#define P_l(n)   (u[(n)][I_p_l])
#define P_h(n)   (u[(n)][I_p_h])
#define P_a(n)   (u[(n)][I_p_a])

#define P_l1(n)  (u_1[(n)][I_p_l])
#define P_h1(n)  (u_1[(n)][I_p_h])
#define P_a1(n)  (u_1[(n)][I_p_a])

#define P_ln(n)  (u_n[(n)][I_p_l])
#define P_hn(n)  (u_n[(n)][I_p_h])
#define P_an(n)  (u_n[(n)][I_p_a])

#define M_l(n)   (f[(n)])
#define M_h(n)   (f[(n+2)])
#define M_a(n)   (f[(n+4)])
#define W_l      (f[(6)])
#define W_h      (f[(7)])
#define W_a      (f[(8)])

#define M_l1(n)  (f_1[(n)])
#define M_h1(n)  (f_1[(n+2)])
#define M_a1(n)  (f_1[(n+4)])
#define W_l1     (f_1[(6)])
#define W_h1     (f_1[(7)])
#define W_a1     (f_1[(8)])

#define M_ln(n)  (f_n[(n)])
#define M_hn(n)  (f_n[(n+2)])
#define M_an(n)  (f_n[(n+4)])
#define W_ln     (f_n[(6)])
#define W_hn     (f_n[(7)])
#define W_an     (f_n[(8)])

#define KD_l     (va[(0)])
#define KF_l     (va[(1)])
#define KD_d     (va[(2)])
#define KF_d     (va[(3)])
#define KD_h     (va[(4)])
#define KF_h     (va[(5)])
#define KD_a     (va[(6)])
#define KF_a     (va[(7)])
/* Fonctions */
static int    pm(char *s) ;
static double saturation(double,double,crbe_t) ;
static double dsaturation(double,double,crbe_t) ;
/* Parametres */
static double gravite,phi,rho_l,M_l,k_int,mu_l,mu_g,M_h2,M_a,p_l0,p_h0,p_a0,D_l,D_ah0,k_hen,RT_0,p_c3 ;

int pm(char *s)
{
if(strcmp(s,"gravite") == 0) return (0) ;
else if(strcmp(s,"phi") == 0) return (1) ;
else if(strcmp(s,"rho_l") == 0) return (2) ;
else if(strcmp(s,"M_l") == 0) return (3) ;
else if(strcmp(s,"k_int") == 0) return (4) ;
else if(strcmp(s,"mu_l") == 0) return (5) ;
else if(strcmp(s,"mu_g") == 0) return (6) ;
else if(strcmp(s,"p_l0") == 0) return (7) ;
else if(strcmp(s,"p_h0") == 0) return (8) ;
else if(strcmp(s,"p_a0") == 0) return (9) ;
else if(strcmp(s,"RT_0") == 0) return (10) ;
else if(strcmp(s,"M_h2") == 0) return (11) ;
else if(strcmp(s,"M_a") == 0) return (12) ;
else if(strcmp(s,"D_ah0") == 0) return (13) ;
else if(strcmp(s,"D_l") == 0) return (14) ;
else if(strcmp(s,"k_hen") == 0) return (15) ;
else if(strcmp(s,"n_p") == 0) return (16) ;
else if(strcmp(s,"p_c1") == 0) return (17) ;
else if(strcmp(s,"p_c2") == 0) return (18) ;
else if(strcmp(s,"p_c3") == 0) return (19) ;
else if(strcmp(s,"courbes") == 0) return (20) ;
else
{ printf("donnee \"%s\" non connue (pm14)\n",s) ; exit(0) ; }
}

int dm14(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 21 ;
  
  mat->neq        = NEQ ;
  strcpy(mat->eqn[E_eau],"liq") ;
  strcpy(mat->eqn[E_hyd],"hyd") ;
  strcpy(mat->eqn[E_air],"air") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_h],"p_h") ;
  strcpy(mat->inc[I_p_a],"p_a") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}

int qm14(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 3 equations :\n\
\t 1. Conservation de masse d\'eau   (p_l)\n\
\t 2. Conservation de masse d\'air   (p_a)\n\
\t 3. Conservation de la masse de H2 (gazeux et dissout) (p_h)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0       # La gravite\n") ;
  fprintf(ficd,"phi = 0.4         # La Porosite\n") ;
  fprintf(ficd,"rho_l = 1000      # Masse volumique du fluide\n") ;
  fprintf(ficd,"M_l = 18.e-3      # Masse molaire du liquide\n") ;
  fprintf(ficd,"M_a = 0.028       # Masse molaire de l\'air\n") ;
  fprintf(ficd,"M_h2 = 0.002      # Masse molaire de h2\n") ;
  fprintf(ficd,"k_int = 6.e-21    # Permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001      # Viscosite du liquide\n") ;
  fprintf(ficd,"mu_g = 1.e-05     # Viscosite du h2\n") ;
  fprintf(ficd,"k_hen = 7.6e-06   # Constante de Henry\n") ;
  fprintf(ficd,"p_c3 = 1e+06      # Pression capillaire limite\n") ;
  fprintf(ficd,"p_l0 = 100000.    # Pression de liquide de reference\n") ;
  fprintf(ficd,"p_h0 = 50000.     # Pression de l\'hydrogene gazeux de reference\n") ;
  fprintf(ficd,"p_a0 = 50000.     # Pression d\'air de reference\n") ;
  fprintf(ficd,"D_ah0 = 1.e-7     # Coefficient de diffusion air_hydrogene (D_ah0)\n") ;
  fprintf(ficd,"D_l = 1.e-11      # Coefficient de diffusion eau_hydrogene (D_l)\n") ;
  fprintf(ficd,"k_hen = 7.6e-06   # Constante de Henry pour l\'hydrogene\n") ;
  fprintf(ficd,"RT_0 = 2479.      # Constante des gaz parfaits fois la temperature\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}

void tb14(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 9 ;
  el->n_ve = 8 ;
}

void ch14(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;

  if(el.nn > 1) arret("ch14 : chargement non prevu") ;
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in14(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */
{
  double p_g[2],rho_h,rho_d,rho_a,c_h[2],c_a[2],D_ah,tau ;
  double p_lij,p_aij,p_hij;
  double pc,sl,sg,pg,ch,ca;
  double dx ;
  int    i ;
  double un = 1.,deux = 2. ;
  
  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  M_l     = el.mat->pr[pm("M_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  M_a     = el.mat->pr[pm("M_a")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_h0    = el.mat->pr[pm("p_h0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_ah0   = el.mat->pr[pm("D_ah0")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
  dx = x[1][0] - x[0][0] ;
  /* masses d'eau, d'hydrogène et d'air */
  for(i=0;i<2;i++) {
    p_g[i] = P_h(i) + P_a(i) ;
    c_h[i] = P_h(i)/p_g[i] ;
    c_a[i] = un - c_h[i] ;
    pc     = p_g[i] - P_l(i) ;
    sl     = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg     = un - sl ;
    rho_d  = P_h(i)*M_h2*k_hen;
    rho_h  = P_h(i)*M_h2/RT_0 ;
    rho_a  = P_a(i)*M_a/RT_0 ;
    M_l(i) = rho_l*phi*sl ;
    M_h(i) = phi*(rho_h*sg + rho_d*sl) ;
    M_a(i) = rho_a*phi*sg ;
  }
  /* coefficient de transfert */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_hij = (P_h(0)+P_h(1))/deux ;
  p_aij = (P_a(0)+P_a(1))/deux ;
  pg    = p_hij + p_aij ;
  pc    = pg - p_lij ;
  sl    = saturation(pc,p_c3,el.mat->cb[0]) ;
  sg    = un - sl ;
  rho_d = p_hij*M_h2*k_hen;
  rho_h = p_hij*M_h2/RT_0 ;
  rho_a = p_aij*M_a/RT_0 ;
  ch    = p_hij/pg ;
  ca    = p_aij/pg ;
  tau   = pow(phi,0.33)*pow(sg,2.33) ;
  D_ah  = D_ah0*(p_h0+p_a0)/pg;
  
  KD_l  = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;  /* Darcy eau liquide */
  KD_d  = rho_d*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;  /* Darcy hyd liquide */
  KF_l  = M_l*D_l*k_hen ;                                   /* Fick eau liquide */
  KF_d  = M_h2*D_l*k_hen ;                                   /* Fick hyd liquide */
  KD_h  = rho_h*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;/* Darcy hyd gazeux */
  KD_a  = rho_a*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;/* Darcy air */
  KF_h  = rho_h/ch*phi*sg*tau*D_ah ;                        /* Fick hyd gazeux */
  KF_a  = rho_a/ca*phi*sg*tau*D_ah ;                        /* Fick air */
  /* flux */
  W_l   = - KD_l*(P_l(1) - P_l(0))/dx + KF_l*(P_h(1) - P_h(0))/dx ;
  W_h   = - KD_d*(P_l(1) - P_l(0))/dx - KF_d*(P_h(1) - P_h(0))/dx- KD_h*(p_g[1] - p_g[0])/dx - KF_h*(c_h[1] - c_h[0])/dx ;
  W_a   = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
}

int ex14(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  double rho_d,rho_h,rho_a,D_ah,tau ;
  double p_lij,p_aij,p_hij ;
  double pc,sl,sg,pg,ch,ca ;
  double un = 1.,deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  M_l     = el.mat->pr[pm("M_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  M_a     = el.mat->pr[pm("M_a")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_h0    = el.mat->pr[pm("p_h0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_ah0   = el.mat->pr[pm("D_ah0")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /*
    COEFFICIENTS DE TRANSFERT
  */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_hij = (P_h(0)+P_h(1))/deux ;
  p_aij = (P_a(0)+P_a(1))/deux ;
  pg    = p_hij + p_aij ;
  pc    = pg - p_lij ;
  sl    = saturation(pc,p_c3,el.mat->cb[0]) ;
  sg    = un - sl ;
  rho_d = p_hij*M_h2*k_hen;
  rho_h = p_hij*M_h2/RT_0 ;
  rho_a = p_aij*M_a/RT_0 ;
  ch    = p_hij/pg ;
  ca    = p_aij/pg ;
  tau   = pow(phi,0.33)*pow(sg,2.33) ;
  D_ah  = D_ah0*(p_h0+p_a0)/pg;
  
  KD_l  = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;  /* Darcy eau liquide */
  KD_d  = rho_d*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;  /* Darcy hyd liquide */
  KF_l  = M_l*D_l*k_hen ;                                   /* Fick eau liquide */
  KF_d  = M_h2*D_l*k_hen ;                                   /* Fick hyd liquide */
  KD_h  = rho_h*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;/* Darcy hyd gazeux */
  KD_a  = rho_a*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;/* Darcy air */
  KF_h  = rho_h/ch*phi*sg*tau*D_ah ;                        /* Fick hyd gazeux */
  KF_a  = rho_a/ca*phi*sg*tau*D_ah ;                        /* Fick air */
  return(0) ;
}

int ct14(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double p_g[2],rho_d,rho_h,rho_a,c_h[2],c_a[2] ;
  double pc,sl,sg ;
  double dx ;
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  M_l     = el.mat->pr[pm("M_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  M_a     = el.mat->pr[pm("M_a")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_h0    = el.mat->pr[pm("p_h0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_ah0   = el.mat->pr[pm("D_ah0")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0   = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /* masses d'eau d'hydrogène et d'air*/
  for(i=0;i<2;i++) {
    /*if(P_a1(i) <= zero)  {
      printf("pression d\'air = %e\n",P_a1(i)) ;
      return(1) ;
      }*/
    p_g[i] = P_h1(i) + P_a1(i) ;
    c_h[i] = P_h1(i)/p_g[i] ;
    c_a[i] = un - c_h[i] ;
    pc     = p_g[i] - P_l1(i) ;
    sl    = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg     = un - sl ;
    rho_d  = P_h1(i)*M_h2*k_hen;
    rho_h  = P_h1(i)*M_h2/RT_0 ;
    rho_a  = P_a1(i)*M_a/RT_0 ;
    
    M_l1(i) = rho_l*phi*sl ;
    M_h1(i) = phi*(rho_h*sg + rho_d*sl) ;
    M_a1(i) = rho_a*phi*sg ;
  }
  /* flux */
  dx    = x[1][0] - x[0][0] ;
  W_l1  = - KD_l*(P_l1(1) - P_l1(0))/dx + KF_l*(P_h1(1) - P_h1(0))/dx ;
  W_h1  = - KD_d*(P_l1(1) - P_l1(0))/dx - KF_d*(P_h1(1) - P_h1(0))/dx- KD_h*(p_g[1] - p_g[0])/dx - KF_h*(c_h[1] - c_h[0])/dx ;
  W_a1   = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
  return(0) ;
}

int mx14(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double p_g[2];
  double rho_d,rho_h,rho_a,c_h[2],c_a[2] ;
  double pc,sl,sg,dslsdpc ;
  double trd_l,trf_l,trd_d,trf_d,trd_a,trf_a,trd_h,trf_h ;
  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  /* Numeros d'ordre des equations et des inconnues locales */
  int E_liq = E_eau;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;
  
  if(el.dim < dim) return(0) ;
  
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  M_l     = el.mat->pr[pm("M_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  M_a     = el.mat->pr[pm("M_a")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_h0    = el.mat->pr[pm("p_h0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_ah0   = el.mat->pr[pm("D_ah0")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0   = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
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
  
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    p_g[i] = P_h1(i) + P_a1(i) ;
    c_h[i] = P_h1(i)/p_g[i] ;
    c_a[i] = un - c_h[i];
    pc     = p_g[i] - P_l1(i) ;
    sl     = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg     = un - sl ;
    rho_d  = P_h1(i)*M_h2*k_hen;
    rho_h  = P_h1(i)*M_h2/RT_0 ;
    rho_a  = P_a1(i)*M_a/RT_0 ;
    dslsdpc = dsaturation(pc,p_c3,el.mat->cb[0]) ;
    /*
      CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
    */
    K(E_liq+i*NEQ,E_liq+i*NEQ)      += volume[i]*phi*rho_l*(-dslsdpc) ;
    K(E_liq+i*NEQ,E_hyd+i*NEQ)      += volume[i]*phi*rho_l*dslsdpc ;
    K(E_liq+i*NEQ,E_air+i*NEQ)      += volume[i]*phi*rho_l*dslsdpc ;
    /*
      CONSERVATION DE L'HYDROGENE  : (m_h1 - m_hn) + dt * div(w_h1) = 0
    */
    K(E_hyd+i*NEQ,E_liq+i*NEQ)      += volume[i]*phi*(rho_d*(-dslsdpc) + rho_h*dslsdpc) ;
    K(E_hyd+i*NEQ,E_hyd+i*NEQ)      += volume[i]*phi*(M_h2*sl*k_hen + rho_d*dslsdpc + M_h2/RT_0*sg - rho_h*dslsdpc) ;
    K(E_hyd+i*NEQ,E_air+i*NEQ)      += volume[i]*phi*(rho_d*dslsdpc + rho_h*(-dslsdpc)) ;
    /*
      CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
    */
    K(E_air+i*NEQ,E_liq+i*NEQ)      += volume[i]*phi*rho_a*dslsdpc ;
    K(E_air+i*NEQ,E_hyd+i*NEQ)      += volume[i]*phi*rho_a*(-dslsdpc) ;
    K(E_air+i*NEQ,E_air+i*NEQ)      += volume[i]*phi*(M_a/RT_0*sg - rho_a*dslsdpc) ;
    
  }
  /*
    termes d'ecoulement
  */
  trd_l = dt*surf/dx*KD_l ;
  trf_l = dt*surf/dx*KF_l ;
  trd_d = dt*surf/dx*KD_d ;
  trf_d = dt*surf/dx*KF_d ;
  trd_h = dt*surf/dx*KD_h ;
  trf_h = dt*surf/dx*KF_h ;
  trd_a = dt*surf/dx*KD_a ;
  trf_a = dt*surf/dx*KF_a ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  K(E_liq,E_liq)      += + trd_l ;
  K(E_liq,E_liq+NEQ)  += - trd_l ;
  K(E_liq,E_hyd)      += - trf_l ;
  K(E_liq,E_hyd+NEQ)  += + trf_l ;
  K(E_liq+NEQ,E_liq)      += - trd_l ;
  K(E_liq+NEQ,E_liq+NEQ)  += + trd_l ;
  K(E_liq+NEQ,E_hyd)      += - trf_l ;
  K(E_liq+NEQ,E_hyd+NEQ)  += + trf_l ;
  /*
    CONSERVATION DE L'HYDROGENE  : (m_h1 - m_hn) + dt * div(w_h1) = 0
  */
  K(E_hyd,E_liq)      += + trd_d ;
  K(E_hyd,E_liq+NEQ)  += - trd_d ;
  K(E_hyd,E_hyd)      += + trf_d + trd_h + trf_h*(c_a[0]/p_g[0]) ;
  K(E_hyd,E_hyd+NEQ)  += - trf_d - trd_h - trf_h*(c_a[1]/p_g[1]) ;
  K(E_hyd,E_air)      += + trd_h + trf_h*(-c_h[0]/p_g[0]) ;
  K(E_hyd,E_air+NEQ)  += - trd_h - trf_h*(-c_h[1]/p_g[1]);
  K(E_hyd+NEQ,E_liq)      += - trd_d ;
  K(E_hyd+NEQ,E_liq+NEQ)  += + trd_d ;
  K(E_hyd+NEQ,E_hyd)      += - trf_d - trd_h - trf_h*(c_a[0]/p_g[0]) ;
  K(E_hyd+NEQ,E_hyd+NEQ)  += + trf_d + trd_h + trf_h*(c_a[1]/p_g[1]) ;
  K(E_hyd+NEQ,E_air)      += - trd_h - trf_h*(-c_h[0]/p_g[0]) ;
  K(E_hyd+NEQ,E_air+NEQ)  += + trd_h + trf_h*(-c_h[1]/p_g[1]);
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
  */
  K(E_air,E_hyd)      += + trd_a + trf_a*(-c_a[0]/p_g[0]) ;
  K(E_air,E_hyd+NEQ)  += - trd_a - trf_a*(-c_a[1]/p_g[1]) ;
  K(E_air,E_air)      += + trd_a + trf_a*(c_h[0]/p_g[0]) ;
  K(E_air,E_air+NEQ)  += - trd_a - trf_a*(c_h[1]/p_g[1]);
  K(E_air+NEQ,E_hyd)      += - trd_a - trf_a*(-c_a[0]/p_g[0]) ;
  K(E_air+NEQ,E_hyd+NEQ)  += + trd_a + trf_a*(-c_a[1]/p_g[1]) ;
  K(E_air+NEQ,E_air)      += - trd_a - trf_a*(c_h[0]/p_g[0]) ;
  K(E_air+NEQ,E_air+NEQ)  += + trd_a + trf_a*(c_h[1]/p_g[1]);
  return(0) ;
  
#undef K
}

void rs14(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
double dx ,xm ;
double volume[2],surf ;
int    i ;
double zero = 0.,un = 1.,deux = 2. ;
/* Numeros d'ordre des equations et des inconnues locales */
int E_liq = E_eau;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;
  
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
R(0,E_liq) -= volume[0]*(M_l1(0) - M_ln(0)) + dt*surf*W_l1 ;
R(1,E_liq) -= volume[1]*(M_l1(1) - M_ln(1)) - dt*surf*W_l1 ;
/*
CONSERVATION DE L'HYDROGENE  : (m_h1 - m_hn) + dt * div(w_h1) = 0
*/
R(0,E_hyd) -= volume[0]*(M_h1(0) - M_hn(0)) + dt*surf*W_h1 ;
R(1,E_hyd) -= volume[1]*(M_h1(1) - M_hn(1)) - dt*surf*W_h1 ;
/*
CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a1) = 0
*/
R(0,E_air) -= volume[0]*(M_a1(0) - M_an(0)) + dt*surf*W_a1 ;
R(1,E_air) -= volume[1]*(M_a1(1) - M_an(1)) - dt*surf*W_a1 ;

#undef R
}

int so14(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_g[2] ;
  double c_h[2],c_a[2] ;
  double pl,ph,pa,pg,pc,sl ;
  double w_l,w_h,w_hld,w_hlf,w_hgd,w_hgf,w_a ;
  double dx ;
  int    i ;
  double un = 1.,deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  M_l     = el.mat->pr[pm("M_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  M_a     = el.mat->pr[pm("M_a")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_h0    = el.mat->pr[pm("p_h0")] ;
  p_a0    = el.mat->pr[pm("p_a0")] ;
  D_ah0   = el.mat->pr[pm("D_ah0")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0   = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /* Calcul de certaines quantites aux noeuds */
  for(i=0;i<2;i++) {
    p_g[i] = P_h(i) + P_a(i) ;
    c_h[i] = P_h(i)/p_g[i] ;
    c_a[i] = un - c_h[i] ;
  }
  /* pressions */
  pl = (P_l(0) + P_l(1))/deux ;
  ph = (P_h(0) + P_h(1))/deux ;
  pg = (p_g[0] + p_g[1])/deux ;
  pa = (P_a(0) + P_a(1))/deux ;
  pc = pg - pl ;
  /* saturation */
  sl     = saturation(pc,p_c3,el.mat->cb[0]) ;
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  w_l  = - KD_l*(P_l(1) - P_l(0))/dx + KF_l*(P_h(1) - P_h(0))/dx ;
  w_hld  = - KD_d*(P_l(1) - P_l(0))/dx ;
  w_hlf  = - KF_d*(P_h(1) - P_h(0))/dx ;
  w_hgd  = - KD_h*(p_g[1] - p_g[0])/dx ;
  w_hgf  = - KF_h*(c_h[1] - c_h[0])/dx ;
  w_h  = w_hld + w_hlf + w_hgd + w_hgf ;
  w_a  = - KD_a*(p_g[1] - p_g[0])/dx - KF_a*(c_a[1] - c_a[0])/dx ;
  
  /* quantites exploitees */
  strcpy(r[0].text,"pression_liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"pression_gaz") ; r[1].n = 1 ;
  r[1].v[0] = pg ;
  strcpy(r[2].text,"pression_hydrogene") ; r[2].n = 1 ;
  r[2].v[0] = ph ;
  strcpy(r[3].text,"pression_air") ; r[3].n = 1 ;
  r[3].v[0] = pa ;
  strcpy(r[4].text,"flux_liquide") ; r[4].n = 1 ;
  r[4].v[0] = w_l ;
  strcpy(r[5].text,"flux_hydrogene") ; r[5].n = 1 ;
  r[5].v[0] = w_h ;
  strcpy(r[6].text,"flux_air") ; r[6].n = 1 ;
  r[6].v[0] = w_a ;
  strcpy(r[7].text,"saturation") ; r[7].n = 1 ;
  r[7].v[0] = sl ;
  return (8) ;
}

double saturation(double pc,double pc3,crbe_t cb)
/* Degre de saturation regularise autour de 1 */
{
  double sl,sl3 ;
  double pc1 = cb.a[0],sl1 = cb.f[0] ;

  if(pc >= pc3 || pc1 >= pc3) sl = courbe(pc,cb) ;
  else {
    sl3 = courbe(pc3,cb) ;
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  return(sl) ;
}

double dsaturation(double pc,double pc3,crbe_t cb)
{
  int    n_i = cb.n - 1 ;
  double dpc = (cb.a[1] - cb.a[0])/n_i/10 ;

  return((saturation(pc + dpc,pc3,cb) - saturation(pc,pc3,cb))/dpc) ;
}

#undef P_l
#undef P_h
#undef P_a

#undef P_l1
#undef P_h1
#undef P_a1

#undef P_ln
#undef P_hn
#undef P_an

#undef M_l
#undef M_h
#undef M_a
#undef W_l
#undef W_h
#undef W_a

#undef M_l1
#undef M_h1
#undef M_a1
#undef W_l1
#undef W_h1
#undef W_a1

#undef M_ln
#undef M_hn
#undef M_an
#undef W_ln
#undef W_hn
#undef W_an

#undef KD_l
#undef KD_d
#undef KD_h
#undef KD_a
#undef KF_l
#undef KF_d
#undef KF_h
#undef KF_a
