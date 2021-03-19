#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  9
#define TITLE   "Ecoulement diphasique H2O-H2"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ   (2)
#define E_liq (0)
#define E_hyd (1)
#define I_p_l (0)
#define I_p_h (1)
/* Fonctions */
static int    pm(char *s) ;
static double saturation(double,double,crbe_t) ;
static double dsaturation(double,double,crbe_t) ;
/* Parametres */
static double gravite,phi,rho_l,k_int,k_hen,D_l,mu_l,mu_g,M_h2,RT_0,p_c3 ;

int pm(char *s)
{
if(strcmp(s,"gravite") == 0) return (0) ;
else if(strcmp(s,"phi") == 0) return (1) ;
else if(strcmp(s,"rho_l") == 0) return (2) ;
else if(strcmp(s,"D_l") == 0) return (3) ;
else if(strcmp(s,"k_int") == 0) return (4) ;
else if(strcmp(s,"mu_l") == 0) return (5) ;
else if(strcmp(s,"mu_g") == 0) return (6) ;
else if(strcmp(s,"RT_0") == 0) return (7) ;
else if(strcmp(s,"M_h2") == 0) return (8) ;
else if(strcmp(s,"k_hen") == 0) return (9) ;
else if(strcmp(s,"p_c3") == 0) return (10) ;
else if(strcmp(s,"courbes") == 0) return (11) ;
else
  { printf("donnee \"%s\" non connue (pm9)\n",s) ; exit(0) ; }
}

int dm9(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int   n_donnees = 12 ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->eqn[E_hyd],"hyd") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_h],"p_h") ;
  
  dmat(mat,ficd,pm,n_donnees) ;
  
  return(mat->n) ;
}

int qm9(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return NEQ ;
  
  printf("\n\n\
Le systeme est forme de 2 equations :\n\
\t 1. Conservation de masse d\'eau   (p_l)\n\
\t 2. Conservation de la masse de H2 (gazeux et dissout) (p_h)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0       # La gravite\n") ;
  fprintf(ficd,"phi = 0.4         # La Porosite\n") ;
  fprintf(ficd,"rho_l = 1000      # Masse volumique du fluide\n") ;
  fprintf(ficd,"M_h = 0.002       # Masse molaire du h2\n") ;
  fprintf(ficd,"k_int = 6.e-21    # Permeabilite intrinseque\n") ;
  fprintf(ficd,"D_l = 1.e-11      # Coeficient de diffusion\n") ;
  fprintf(ficd,"mu_l = 0.001      # Viscosite du liquide\n") ;
  fprintf(ficd,"mu_g = 1.e-05     # Viscosite du h2\n") ;
  fprintf(ficd,"RT_0 = 2479       # Temperature fois R (RT_0)\n") ;
  fprintf(ficd,"k_hen = 7.6e-06   # Constante de Henry\n") ;
  fprintf(ficd,"p_c3 = 1e+06      # Pression capillaire limite\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}

void tb9(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 6 ;
  el->n_ve = 3 ;
}

void ch9(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in9(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_h(i)   (u[(i)][I_p_h])
#define M_l(i)   (f[(i)])
#define M_h(i)   (f[(i+2)])
#define W_l      (f[(4)])
#define W_g      (f[(5)])
#define K_ll     (va[(0)])
#define K_ld     (va[(1)])
#define K_dg     (va[(2)])
  double pc,sl,dx,p_lij,p_hij,rho_d,K_g,rho_g,sg ;
  int    i ;
  double un = 1.,deux = 2. ;
  
  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
  dx = x[1][0] - x[0][0] ;
  /* MASSE FLUIDE */
  for(i=0;i<2;i++) {
    pc     = P_h(i) - P_l(i) ;
    sl     = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg     = un - sl ;
    rho_g  = P_h(i)*M_h2/RT_0 ;
    rho_d  = P_h(i)*M_h2*k_hen ;
    M_l(i) = rho_l*phi*sl ;
    M_h(i) = phi*(rho_g*sg + rho_d*sl) ;
  }
  /* COEFFICIENTS DE TRANSFERT */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_hij = (P_h(0)+P_h(1))/deux ; 
  pc    = p_hij - p_lij ;
  rho_g = M_h2*p_hij/RT_0 ;
  rho_d = p_hij*M_h2*k_hen ;
  K_ll  = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
  K_g   = rho_g*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;
  K_ld  = rho_d*K_ll/rho_l ;
  K_dg  = K_g + M_h2*k_hen*D_l ;
  /* FLUX */
  W_l   = - K_ll*(P_l(1) - P_l(0))/dx + K_ll*rho_l*gravite ;
  W_g   = - K_ld*(P_l(1) - P_l(0))/dx - K_dg*(P_h(1) - P_h(0))/dx ;

#undef P_l
#undef P_h
#undef M_l
#undef M_h
#undef W_l
#undef W_g
#undef K_ll
#undef K_ld
#undef K_dg
}

int ex9(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_h(i)   (u[(i)][I_p_h])
#define K_ll     (va[(0)])
#define K_ld     (va[(1)])
#define K_dg     (va[(2)])
  double pc,p_lij,p_hij,K_g,rho_d,rho_g ;
  double deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /*
    COEFFICIENTS DE TRANSFERT
  */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_hij = (P_h(0)+P_h(1))/deux ; 
  pc    = p_hij - p_lij ;
  rho_g = p_hij*M_h2/RT_0 ;
  rho_d = p_hij*M_h2*k_hen ;
  K_ll  = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
  K_g   = rho_g*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;
  K_ld  = rho_d*K_ll/rho_l ;
  K_dg  = K_g + M_h2*k_hen*D_l ;
  return(0) ;
  
#undef P_l
#undef P_h
#undef K_ll
#undef K_ld
#undef K_dg
}

int ct9(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define P_l1(i)   (u_1[(i)][I_p_l])
#define P_h1(i)   (u_1[(i)][I_p_h])
#define M_l1(i)   (f_1[(i)])
#define M_h1(i)   (f_1[(i+2)])
#define W_l1      (f_1[(4)])
#define W_g1      (f_1[(5)])
#define K_ll      (va[(0)])
#define K_ld      (va[(1)])
#define K_dg      (va[(2)])
  double pc,sl,dx,rho_g,rho_d,sg ;
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /* masse de fluide */
  for(i=0;i<2;i++) {
    pc      = P_h1(i) - P_l1(i) ;
    sl      = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg      = un - sl ;
    rho_g   = P_h1(i)*M_h2/RT_0 ;
    rho_d   = P_h1(i)*M_h2*k_hen ;
    if(P_h1(i) < 0.) {printf("  stop p_h<0\n") ; exit(0) ;}
    M_l1(i) = rho_l*phi*sl ;
    M_h1(i) = phi*(rho_g*sg + rho_d*sl) ;
  }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  W_l1 = - K_ll*(P_l1(1) - P_l1(0))/dx + K_ll*rho_l*gravite ;
  W_g1 = - K_ld*(P_l1(1) - P_l1(0))/dx - K_dg*(P_h1(1) - P_h1(0))/dx ;

  return(0) ;
  
#undef P_l1
#undef P_h1
#undef M_l1
#undef M_h1
#undef W_l1
#undef W_g1
#undef K_ll
#undef K_ld
#undef K_dg
}

int mx9(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
#define P_l1(i)   (u_1[(i)][I_p_l])
#define P_h1(i)   (u_1[(i)][I_p_h])
#define K_ll      (va[(0)])
#define K_ld      (va[(1)])
#define K_dg      (va[(2)])
  double pc,sl,dx,xm,tr_ll,tr_ld,tr_dg,dslsdpc,rho_g,rho_d,sg ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
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
  
  /* initialisation */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = zero ;
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    pc      = P_h1(i) - P_l1(i) ;
    sl      = saturation(pc,p_c3,el.mat->cb[0]) ;
    sg      = un - sl ;
    dslsdpc = dsaturation(pc,p_c3,el.mat->cb[0]) ;
    rho_g   = P_h1(i)*M_h2/RT_0 ;
    rho_d   = P_h1(i)*M_h2*k_hen ;
    /*
      CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
    */
    K(i*NEQ+E_liq,i*NEQ+I_p_l) += volume[i]*rho_l*phi*(-dslsdpc) ;
    K(i*NEQ+E_liq,i*NEQ+I_p_h) += volume[i]*rho_l*phi*dslsdpc ;
    /*
      CONSERVATION DE H2 : (m_h1 - m_hn) + dt * div(w_h1) = 0
    */
    K(i*NEQ+E_hyd,i*NEQ+I_p_l) += volume[i]*(rho_d-rho_g)*phi*(-dslsdpc) ;
    K(i*NEQ+E_hyd,i*NEQ+I_p_h) += volume[i]*phi*(M_h2/RT_0*sg + M_h2*k_hen*sl 
						 + (rho_d-rho_g)*dslsdpc) ;
  }
  /* termes d'ecoulement */
  tr_ll = dt*surf*K_ll/dx ;
  tr_ld = dt*surf*K_ld/dx ;
  tr_dg = dt*surf*K_dg/dx ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  K(E_liq,I_p_l)         += + tr_ll ;
  K(E_liq,NEQ+I_p_l)     += - tr_ll;
  K(NEQ+E_liq,I_p_l)     += - tr_ll ;
  K(NEQ+E_liq,NEQ+I_p_l) += + tr_ll ;
  /*
    CONSERVATION DE H2 : (m_g1 - m_gn) + dt * div(w_g1) = 0
  */ 
  K(E_hyd,I_p_l)         += + tr_ld ;
  K(E_hyd,I_p_h)         += + tr_dg ;
  K(E_hyd,NEQ+I_p_l)     += - tr_ld ;
  K(E_hyd,NEQ+I_p_h)     += - tr_dg ;
  K(NEQ+E_hyd,I_p_l)     += - tr_ld ;
  K(NEQ+E_hyd,I_p_h)     += - tr_dg ;
  K(NEQ+E_hyd,NEQ+I_p_l) += + tr_ld ;
  K(NEQ+E_hyd,NEQ+I_p_h) += + tr_dg ;

  return(0) ;

#undef K
#undef P_l1
#undef P_h1
#undef K_ll
#undef K_ld
#undef K_dg
}

void rs9(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

#define M_l1(i)   (f_1[(i)])
#define M_ln(i)   (f_n[(i)])
#define M_h1(i)   (f_1[(i+2)])
#define M_hn(i)   (f_n[(i+2)])
#define W_l1      (f_1[(4)])
#define W_g1      (f_1[(5)])
#define K_ll      (va[(0)])
#define K_ld      (va[(1)])
#define K_dg      (va[(2)])
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;
  
  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
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
  
  /* initialisation */
  for(i=0;i<2*NEQ;i++) r[i] = zero ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  R(0,E_liq) -= volume[0]*(M_l1(0) - M_ln(0)) + dt*surf*W_l1 ;
  R(1,E_liq) -= volume[1]*(M_l1(1) - M_ln(1)) - dt*surf*W_l1 ;
  /*
    CONSERVATION DE H2 : (m_g1 - m_gn) + dt * div(w_g1) = 0
  */ 
  R(0,E_hyd) -= volume[0]*(M_h1(0) - M_hn(0)) + dt*surf*W_g1 ;
  R(1,E_hyd) -= volume[1]*(M_h1(1) - M_hn(1)) - dt*surf*W_g1 ;
  
#undef R
#undef M_l1
#undef M_ln
#undef M_h1
#undef M_hn
#undef W_l1
#undef W_g1
}

int so9(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_h(i)   (u[(i)][I_p_h])
#define W_l      (f[(4)])
#define W_g      (f[(5)])
  double pc,sl,pl,pg ;
  double deux = 2. ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  D_l     = el.mat->pr[pm("D_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_h2    = el.mat->pr[pm("M_h2")] ;
  k_hen   = el.mat->pr[pm("k_hen")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /* Pressions */
  pl = (P_l(0) + P_l(1))/deux ;
  pg = (P_h(0) + P_h(1))/deux ;
  /* saturation */
  pc = pg - pl ;
  sl = saturation(pc,p_c3,el.mat->cb[0]) ;
  /* quantites exploitees */
  strcpy(r[0].text,"pression_liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"pression_h2") ; r[1].n = 1 ;
  r[1].v[0] = pg ;
  strcpy(r[2].text,"flux_liquide") ; r[2].n = 1 ;
  r[2].v[0] = W_l ;
  strcpy(r[3].text,"flux_h2") ; r[3].n = 1 ;
  r[3].v[0] = W_g ;
  strcpy(r[4].text,"saturation") ; r[4].n = 1 ;
  r[4].v[0] = sl ;
  return (5) ;
  
#undef P_l
#undef P_h
#undef W_l
#undef W_g
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
  double dpc = (cb.a[1] - cb.a[0])/n_i ;

  return((saturation(pc + dpc,pc3,cb) - saturation(pc,pc3,cb))/dpc) ;
}


#undef NEQ
#undef E_liq
#undef E_hyd
#undef I_p_l
#undef I_p_h
