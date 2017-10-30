/*
equations de richards
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  1
#include "OldMethods.h"

#define TITLE "Equation de Richards (1D)"
/* Macros */
#define NEQ (1)
#define NVE (1)
#define NVI (3)

#define E_liq (0)
#define I_p_l (0)

#define P_l(n)   (u[(n)][I_p_l])

#define M_l(n)   (f[(n)])
#define W_l      (f[(2)])
#define M_ln(n)  (f_n[(n)])

#define K_l      (va[(0)])

/* Fonctions */
static int    pm(char *s) ;
/* Parametres */
static double gravite,phi,rho_l,k_int,mu_l,p_g,schema ;

int pm(char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"p_g") == 0) return (5) ;
  else if(strcmp(s,"courbes") == 0) return (6) ;
  else if(strcmp(s,"schema") == 0) return (6) ;
  else {
    printf("donnee \"%s\" non connue (pm1)\n",s) ;
    exit(0) ;
  }
}

int dm1(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 7 ;

  if(dim > 1) arret("dm1 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->inc[I_p_l],"p_l") ;

  mat->pr[pm("schema")] = -1 ;
  dmat(mat,ficd,pm,n_donnees) ;
  schema  = mat->pr[pm("schema")] ;
  if(schema != -1) dmat(mat,ficd,pm,1) ;
  return(mat->n) ;
}

int qm1(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  fprintf(stdout,"\n\n\
L\'inconnue est la pression de liquide : p_l.\n\
Exemple de donnees :\n\n") ;

  fprintf(ficd,"gravite = -9.81  # La gravite\n") ;
  fprintf(ficd,"phi = 0.38       # La porosite\n") ;
  fprintf(ficd,"rho_l = 1000     # La masse volumique du fluide\n") ;
  fprintf(ficd,"k_int = 8.9e-12  # La permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001     # La viscosite du fluide\n") ;
  fprintf(ficd,"p_g = 1.e5       # La pression du gaz\n") ;
  fprintf(ficd,"courbes = billes # Le nom du fichier p_c S_l k_rl\n") ;

  return(NEQ) ;
}

void tb1(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch1(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in1(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  
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

  /* MASSE LIQUIDE */
  for(i=0;i<2;i++) {
    double pc  = p_g - P_l(i) ;
    double sl  = courbe(pc,el.mat->cb[0]) ;
    M_l(i) = rho_l*phi*sl ;
  }
  /* COEFFICIENTS DE TRANSFERT */
  {
    double pl   = (P_l(0) + P_l(1))*0.5 ;
    double pc   = p_g - pl ;
    double k_rl = courbe(pc,el.mat->cb[1]) ;
    K_l = rho_l*k_int/mu_l*k_rl ;
  }
  /* FLUX */
  {
    double dx = x[1][0] - x[0][0] ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    W_l = - K_l*grd_p_l + K_l*rho_l*gravite ;
  }

}

int ex1(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
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
  {
    double pl = (P_l(0) + P_l(1))*0.5 ;
    double pc ;
    double k_rl ;
    if(schema == 1) {
      double dx   = x[1][0] - x[0][0] ;
      double grdh = (P_l(1) - P_l(0))/dx - rho_l*gravite ;
      pl = (grdh > 0) ? P_l(1) : P_l(0) ;
    }
    pc   = p_g - pl ;
    k_rl = courbe(pc,el.mat->cb[1]) ;
    K_l  = rho_l*k_int/mu_l*k_rl ;
  }
  return(0) ;

}

int ct1(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
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
    double pc = p_g - P_l(i) ;
    double sl = courbe(pc,el.mat->cb[0]) ;
    M_l(i) = rho_l*phi*sl ;
  }

  /* flux */
  {
    double dx   = x[1][0] - x[0][0] ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    W_l = - K_l*grd_p_l + K_l*rho_l*gravite ;
  }
  
  return(0) ;
}

int mx1(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = zero ;

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
    double pc      = p_g - P_l(i) ;
    double dslsdpc = dcourbe(pc,el.mat->cb[0]) ;
    int    j = i*NEQ ;
    K(E_liq+j,I_p_l+j) += volume[i]*rho_l*phi*(-dslsdpc) ;
  }
  /* termes d'ecoulement */
  {
    double tr_l = dt*surf*K_l/dx ;
    K(E_liq,I_p_l)         += + tr_l ;
    K(E_liq,I_p_l+NEQ)     += - tr_l ;
    K(E_liq+NEQ,I_p_l)     += - tr_l ;
    K(E_liq+NEQ,I_p_l+NEQ) += + tr_l ;
  }

  return(0) ;

#undef K
}

void rs1(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
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
  r[0] -= volume[0]*(M_l(0) - M_ln(0)) + dt*surf*W_l ;
  r[1] -= volume[1]*(M_l(1) - M_ln(1)) - dt*surf*W_l ;

}

int so1(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
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

  /* quantites exploitees */
  {
    /* pression */
    double pl =  param(u,h_s,el.nn,I_p_l) ;
    /* saturation */
    double pc = p_g - pl ;
    double sl = courbe(pc,el.mat->cb[0]) ;
    /* flux */
    double dx = x[1][0] - x[0][0] ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    double wl = - K_l*grd_p_l + K_l*rho_l*gravite ;
    
    strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
    r[0].v[0] = pl ;
    strcpy(r[1].text,"flux-liquide") ; r[1].n = 3 ;
    r[1].v[0] = wl ;
    strcpy(r[2].text,"saturation") ; r[2].n = 1 ;
    r[2].v[0] = sl ;
  }

  return(nso) ;
}
