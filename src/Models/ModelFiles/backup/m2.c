#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  2
#define TITLE "Ecoulement diphasique liquide-gaz"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ   (2)
#define E_liq (0)
#define E_gaz (1)
#define I_p_l (0)
#define I_p_g (1)
/* Fonctions */
static int    pm(char *s) ;
static double saturation(double,double,crbe_t) ;
static double dsaturation(double,double,crbe_t) ;
/* Parametres */
static double gravite,phi,rho_l,k_int,mu_l,mu_g,M_gsR,T_0,p_c3 ;

int pm(char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"mu_g") == 0) return (5) ;
  else if(strcmp(s,"T_0") == 0) return (6) ;
  else if(strcmp(s,"M_gsR") == 0) return (7) ;
  else if(strcmp(s,"p_c3") == 0) return (8) ;
  else if(strcmp(s,"courbes") == 0) return (9) ;
  else
    { printf("donnee \"%s\" non connue (pm2)\n",s) ; exit(0) ; }
}

int dm2(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 10 ;

  if(dim > 1) arret("dm2 : dimension > 1 non prevue") ;
  
  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->eqn[E_gaz],"gaz") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_g],"p_g") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm2(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
 Le systeme d\'equations est forme de 2 equations :\n\
\t 1. Conservation de la masse liquide (p_l) \n\
\t 2. Conservation de la masse gazeuse (p_g) \n\
Exemple de donnees :\n") ;

  fprintf(ficd,"gravite = -9.81    # Gravite\n") ;
  fprintf(ficd,"phi = 0.3          # Porosite\n") ;
  fprintf(ficd,"rho_l = 1000       # masse volumique d liquide\n") ;
  fprintf(ficd,"M_gsR = 3.46374e-3 # Masse molaire du gaz divisee par R\n") ;
  fprintf(ficd,"k_int = 4.4e-13    # Permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 1.e-3       # Viscosite du liquide\n") ;
  fprintf(ficd,"mu_g = 1.8e-5      # Viscosite du gaz\n") ;
  fprintf(ficd,"T_0 = 293          # Temperature\n") ;
  fprintf(ficd,"p_c3 = 300         # Pression capillaire limite\n") ;
  fprintf(ficd,"courbes = sol      # Nom du fichier p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}

void tb2(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 6 ;
  el->n_ve = 2 ;
}

void ch2(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;
  
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in2(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_g(i)   (u[(i)][I_p_g])
#define M_l(i)   (f[(i)])
#define M_g(i)   (f[(i+2)])
#define W_l      (f[(4)])
#define W_g      (f[(5)])
#define K_l      (va[(0)])
#define K_g      (va[(1)])
  double pc,sl,dx,p_lij,p_gij,rho_g ;
  int    i ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_gsR   = el.mat->pr[pm("M_gsR")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  
  dx = x[1][0] - x[0][0] ;
  /* MASSE FLUIDE */
  for(i=0;i<2;i++)
    {
      pc     = P_g(i) - P_l(i) ;
      sl     = saturation(pc,p_c3,el.mat->cb[0]) ;
      rho_g  = P_g(i)*M_gsR/T_0 ;
      M_l(i) = rho_l*phi*sl ;
      M_g(i) = rho_g*phi*(un - sl) ;
    }
  /* COEFFICIENTS DE TRANSFERT */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_gij = (P_g(0)+P_g(1))/deux ; 
  pc    = p_gij - p_lij ;
  rho_g = M_gsR*p_gij/T_0 ;
  K_l   = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
  K_g   = rho_g*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;
  /* FLUX */
  W_l   = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  W_g   = - K_g*(P_g(1) - P_g(0))/dx ;
  
#undef P_l
#undef P_g
#undef M_l
#undef M_g
#undef W_l
#undef W_g
#undef K_l
#undef K_g
}

int ex2(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_g(i)   (u[(i)][I_p_g])
#define K_l      (va[(0)])
#define K_g      (va[(1)])
  double pc,p_lij,p_gij,rho_g ;
  double deux = 2. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_gsR   = el.mat->pr[pm("M_gsR")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /*
    COEFFICIENTS DE TRANSFERT
  */
  p_lij = (P_l(0)+P_l(1))/deux ;
  p_gij = (P_g(0)+P_g(1))/deux ; 
  pc    = p_gij - p_lij ;
  rho_g = M_gsR*p_gij/T_0 ;
  K_l   = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
  K_g   = rho_g*k_int/mu_g*courbe(pc,el.mat->cb[2]) ;
  return(0) ;
  
#undef P_l
#undef P_g
#undef K_l
#undef K_g
}

int ct2(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define P_l1(i)   (u_1[(i)][I_p_l])
#define P_g1(i)   (u_1[(i)][I_p_g])
#define M_l1(i)   (f_1[(i)])
#define M_g1(i)   (f_1[(i+2)])
#define W_l1      (f_1[(4)])
#define W_g1      (f_1[(5)])
#define K_l       (va[(0)])
#define K_g       (va[(1)])
  double pc,sl,dx,rho_g ;
  int    i ;
  double un = 1. ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_gsR   = el.mat->pr[pm("M_gsR")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  /* masse de fluide */
  for(i=0;i<2;i++)
    {
      pc      = P_g1(i) - P_l1(i) ;
      sl      = saturation(pc,p_c3,el.mat->cb[0]) ;
      rho_g   = P_g1(i)*M_gsR/T_0 ;
      if(P_g1(i) < 0.) return(-1) ;
      M_l1(i) = rho_l*phi*sl ;
      M_g1(i) = rho_g*phi*(un - sl) ;
    }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  W_l1 = - K_l*(P_l1(1) - P_l1(0))/dx + K_l*rho_l*gravite ;
  W_g1 = - K_g*(P_g1(1) - P_g1(0))/dx ;

  return(0) ;
  
#undef P_l1
#undef P_g1
#undef M_l1
#undef M_g1
#undef W_l1
#undef W_g1
#undef K_l
#undef K_g
}

int  mx2(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
#define P_l1(i)   (u_1[(i)][I_p_l])
#define P_g1(i)   (u_1[(i)][I_p_g])
#define K_l       (va[(0)])
#define K_g       (va[(1)])
  double pc,sl,dx,xm,tr_l,tr_g,dslsdpc,dpc,rho_g ;
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
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_gsR   = el.mat->pr[pm("M_gsR")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;
  dpc     = (el.mat->cb[0].a[1] - el.mat->cb[0].a[0])/(el.mat->cb[0].n - 1) ;
  dpc     = dpc/100 ;
  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++)
    {
      volume[i] = fabs(dx)/deux ; 
      if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
    }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  
  /* termes d'accumulation */
  for(i=0;i<2;i++)
    {
      pc      = P_g1(i) - P_l1(i) ;
      sl      = saturation(pc,p_c3,el.mat->cb[0]) ;
      dslsdpc = dsaturation(pc,p_c3,el.mat->cb[0]) ;
      rho_g   = P_g1(i)*M_gsR/T_0 ;
      /*
	CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
      */
      K(i*NEQ+E_liq,i*NEQ+I_p_l) += volume[i]*rho_l*phi*(-dslsdpc) ;
      K(i*NEQ+E_liq,i*NEQ+I_p_g) += volume[i]*rho_l*phi*dslsdpc ;
      /*
	CONSERVATION DU GAZ : (m_g1 - m_gn) + dt * div(w_g1) = 0
      */
      K(i*NEQ+E_gaz,i*NEQ+I_p_l) += volume[i]*rho_g*phi*dslsdpc ;
      K(i*NEQ+E_gaz,i*NEQ+I_p_g) += volume[i]*(M_gsR/T_0*phi*(un - sl) 
					       - rho_g*phi*dslsdpc) ;
    }
  /* termes d'ecoulement */
  tr_l = dt*surf*K_l/dx ;
  tr_g = dt*surf*K_g/dx ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l1) = 0
  */
  K(E_liq,I_p_l)         += + tr_l ;
  K(E_liq,NEQ+I_p_l)     += - tr_l ;
  K(NEQ+E_liq,I_p_l)     += - tr_l ;
  K(NEQ+E_liq,NEQ+I_p_l) += + tr_l ;
  /*
    CONSERVATION DU GAZ : (m_g1 - m_gn) + dt * div(w_g1) = 0
  */ 
  K(E_gaz,I_p_g)         += + tr_g ;
  K(E_gaz,NEQ+I_p_g)     += - tr_g ;
  K(NEQ+E_gaz,I_p_g)     += - tr_g ;
  K(NEQ+E_gaz,NEQ+I_p_g) += + tr_g ;

  return(0) ;
  
#undef K
#undef P_l1
#undef P_g1
#undef K_l
#undef K_g
}

void rs2(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define M_l1(i)   (f_1[(i)])
#define M_ln(i)   (f_n[(i)])
#define M_g1(i)   (f_1[(i+2)])
#define M_gn(i)   (f_n[(i+2)])
#define W_l1      (f_1[(4)])
#define W_g1      (f_1[(5)])
#define K_l       (va[(0)])
#define K_g       (va[(1)])
#define R(n,i)    (r[(n)*NEQ+(i)])
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;
  
  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++)
    {
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
    CONSERVATION DU GAZ : (m_g1 - m_gn) + dt * div(w_g1) = 0
  */ 
  R(0,E_gaz) -= volume[0]*(M_g1(0) - M_gn(0)) + dt*surf*W_g1 ;
  R(1,E_gaz) -= volume[1]*(M_g1(1) - M_gn(1)) - dt*surf*W_g1 ;
  
#undef R
#undef M_l1
#undef M_ln
#undef M_g1
#undef M_gn
#undef W_l1
#undef W_g1
}

int so2(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define P_l(i)   (u[(i)][I_p_l])
#define P_g(i)   (u[(i)][I_p_g])
#define K_l      (va[(0)])
#define K_g      (va[(1)])
  double pl,pg,pc,sl,dx,wl,wg ;
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
  mu_g    = el.mat->pr[pm("mu_g")] ;
  M_gsR   = el.mat->pr[pm("M_gsR")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  p_c3    = el.mat->pr[pm("p_c3")] ;

  /* initialisation */
  nso = 5 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  /* Pressions */
  pl  = param(u,h_s,el.nn,I_p_l) ;
  pg  = param(u,h_s,el.nn,I_p_g) ;
  /* saturation */
  pc  = pg - pl ;
  sl  = saturation(pc,p_c3,el.mat->cb[0]) ;
  /* flux */
  dx  = x[1][0] - x[0][0] ;
  wl  = - K_l*(P_l(1) - P_l(0))/dx + K_l*rho_l*gravite ;
  wg  = - K_g*(P_g(1) - P_g(0))/dx ;
  /* quantites exploitees */
  strcpy(r[0].text,"p_l") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"p_g") ; r[1].n = 1 ;
  r[1].v[0] = pg ;
  strcpy(r[2].text,"w_l") ; r[2].n = 1 ;
  r[2].v[0] = wl ;
  strcpy(r[3].text,"w_g") ; r[3].n = 1 ;
  r[3].v[0] = wg ;
  strcpy(r[4].text,"saturation") ; r[4].n = 1 ;
  r[4].v[0] = sl ;
  return (nso) ;
  
#undef P_l
#undef P_g
#undef K_l
#undef K_g
}

double saturation(double pc,double pc3,crbe_t cb)
/* Degre de saturation regularise autour de 1 */
{
  double sl,sl1,sl3 ;
  double pc1 = cb.a[0] ;
  
  if(pc >= pc3 || pc1 >= pc3) sl = courbe(pc,cb) ;
  else {
    sl1 = cb.f[0] ;
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
#undef E_gaz
#undef I_p_l
#undef I_p_g
