#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  11
#define TITLE "Electro-diffusion"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ      (3)
#define E_cat    (0)
#define E_ani    (1)
#define E_poi    (2)
#define I_c_c    (0)
#define I_c_a    (1)
#define I_psi    (2)
/* Fonctions */
static int    pm(const char *s) ;
/* Parametres */
static double phi,d_c,d_a,RT,eps,F,z ;

int pm(const char *s)
{
  if(!strcmp(s,"permittivite")) return (0) ;
  else if(!strcmp(s,"porosite")) return (1) ;
  else if(!strcmp(s,"D_c")) return (2) ;
  else if(!strcmp(s,"D_a")) return (3) ;
  else if(!strcmp(s,"RT")) return (4) ;
  else if(!strcmp(s,"F")) return (5) ;
  else if(!strcmp(s,"valence")) return (6) ;
  else
    { printf("donnee \"%s\" non connue (pm11)\n",s) ; exit(0) ; }
}

int dm11(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 7 ;

  if(dim > 1) arret("dm11 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_cat],"cations") ;
  strcpy(mat->eqn[E_ani],"anions") ;
  strcpy(mat->eqn[E_poi],"Poisson") ;
  strcpy(mat->inc[I_c_c],"c_c") ;
  strcpy(mat->inc[I_c_a],"c_a") ;
  strcpy(mat->inc[I_psi],"psi") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}


int qm11(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme d\'equation est forme de 3 equations:\n\
\t 1. l\'equation de Poisson (psi)\n\
\t 2. l\'equation de conservation de la masse des cations (c_c)\n\
\t 3. l\'equation de conservation de la masse des anions (c_a)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.3         # La porosite\n") ;
  fprintf(ficd,"permittivite = 8.85e11 # La permittivite\n") ;
  fprintf(ficd,"D_c = 0.6e-11          # Le coefficient de diffusion effectif du cation\n") ;
  fprintf(ficd,"D_a = 1.e-11           # Le coefficient de diffusion effectif de l'anion\n") ;
  fprintf(ficd,"RT = 2.43e3            # Le produit de la constante R avec la temperature T\n") ;
  fprintf(ficd,"F = 9.648e4            # La constante de Faraday\n") ;
  fprintf(ficd,"valence = 1            # La valence des ions\n") ;

  return(NEQ) ;
}

void tb11(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 9 ;
  el->n_ve = 2 ;
}

void ch11(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in11(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define C_c(n)   (u[(n)][I_c_c])
#define C_a(n)   (u[(n)][I_c_a])
#define PSI(n)   (u[(n)][I_psi])
#define N_c(n)   (f[(n)])
#define N_a(n)   (f[(n+2)])
#define Q(n)     (f[(n+4)])
#define J_c      (f[(6)])
#define J_a      (f[(7)])
#define E        (f[(8)])
#define K_c      (va[(0)])
#define K_a      (va[(1)])
  double dx,c_c,c_a ;
  int    i ;
  double deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  eps     = el.mat->pr[pm("permittivite")] ;
  d_c     = el.mat->pr[pm("D_c")] ;
  d_a     = el.mat->pr[pm("D_a")] ;
  RT      = el.mat->pr[pm("RT")] ;
  F       = el.mat->pr[pm("F")] ;
  z       = el.mat->pr[pm("valence")] ;
  
  dx = x[1][0] - x[0][0] ;
  /* nb de moles + et - */
  for(i=0;i<2;i++) {
    N_c(i) = phi*C_c(i) ;
    N_a(i) = phi*C_a(i) ;
    Q(i)   = F*z*(N_c(i) - N_a(i)) ;
  }
  /* coefficient de transfert */
  c_c  = (C_c(0) + C_c(1))/deux ;
  c_a  = (C_a(0) + C_a(1))/deux ;
  K_c  = d_c*z*F/RT*c_c ;
  K_a  = - d_a*z*F/RT*c_a ;
  /* flux */
  J_c  = - d_c*(C_c(1) - C_c(0))/dx - K_c*(PSI(1) - PSI(0))/dx ;
  J_a  = - d_a*(C_a(1) - C_a(0))/dx - K_a*(PSI(1) - PSI(0))/dx ;
  E    = - (PSI(1) - PSI(0))/dx ;

#undef C_c
#undef C_a
#undef PSI
#undef N_c
#undef N_a
#undef Q
#undef J_c
#undef J_a
#undef E
#undef K_c
#undef K_a
}

int ex11(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define C_c(n)   (u[(n)][I_c_c])
#define C_a(n)   (u[(n)][I_c_a])
#define K_c      (va[(0)])
#define K_a      (va[(1)])
  double c_c,c_a ;
  double deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  eps     = el.mat->pr[pm("permittivite")] ;
  d_c     = el.mat->pr[pm("D_c")] ;
  d_a     = el.mat->pr[pm("D_a")] ;
  RT      = el.mat->pr[pm("RT")] ;
  F       = el.mat->pr[pm("F")] ;
  z       = el.mat->pr[pm("valence")] ;
  /*
    COEFFICIENTS DE TRANSFERT
  */
  c_c  = (C_c(0) + C_c(1))/deux ;
  c_a  = (C_a(0) + C_a(1))/deux ;
  K_c  = d_c*z*F/RT*c_c ;
  K_a  = - d_a*z*F/RT*c_a ;
  return(0) ;
#undef C_c
#undef C_a
#undef K_c
#undef K_a
}

int ct11(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define C_c(n)   (u_1[(n)][I_c_c])
#define C_a(n)   (u_1[(n)][I_c_a])
#define PSI(n)   (u_1[(n)][I_psi])
#define N_c(n)   (f_1[(n)])
#define N_a(n)   (f_1[(n+2)])
#define Q(n)     (f_1[(n+4)])
#define J_c      (f_1[(6)])
#define J_a      (f_1[(7)])
#define E        (f_1[(8)])
#define K_c      (va[(0)])
#define K_a      (va[(1)])
  double dx ;
  int    i ;
  
  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  eps     = el.mat->pr[pm("permittivite")] ;
  d_c     = el.mat->pr[pm("D_c")] ;
  d_a     = el.mat->pr[pm("D_a")] ;
  RT      = el.mat->pr[pm("RT")] ;
  F       = el.mat->pr[pm("F")] ;
  z       = el.mat->pr[pm("valence")] ;
  /* masses d'eau et d'air, entropie */
  for(i=0;i<2;i++) {
    N_c(i) = phi*C_c(i) ;
    N_a(i) = phi*C_a(i) ;
    Q(i)   = F*z*(N_c(i) - N_a(i)) ;
  }
  /* flux */
  dx   = x[1][0] - x[0][0] ;
  J_c  = - d_c*(C_c(1) - C_c(0))/dx - K_c*(PSI(1) - PSI(0))/dx ;
  J_a  = - d_a*(C_a(1) - C_a(0))/dx - K_a*(PSI(1) - PSI(0))/dx ;
  E    = - (PSI(1) - PSI(0))/dx ;
  return(0) ;
#undef C_c
#undef C_a
#undef PSI
#undef N_c
#undef N_a
#undef Q
#undef J_c
#undef J_a
#undef E
#undef K_c
#undef K_a
}

int mx11(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K_c      (va[(0)])
#define K_a      (va[(1)])
#define K(i,j)   (k[(i)*2*NEQ+(j)])
  double trd_c,trd_a,trm_c,trm_a,tr_e ;
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  eps     = el.mat->pr[pm("permittivite")] ;
  d_c     = el.mat->pr[pm("D_c")] ;
  d_a     = el.mat->pr[pm("D_a")] ;
  RT      = el.mat->pr[pm("RT")] ;
  F       = el.mat->pr[pm("F")] ;
  z       = el.mat->pr[pm("valence")] ;
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
    /*
      CONSERVATION DE +   : (n_c1 - n_cn) + dt * div(J_c1) = 0
    */
    K(E_cat+i*NEQ,I_c_c+i*NEQ) += volume[i]*phi ;
    /*
      CONSERVATION DE -   : (n_a1 - n_an) + dt * div(J_a1) = 0
    */
    K(E_ani+i*NEQ,I_c_a+i*NEQ) += volume[i]*phi ;
    /*
      EQUATION DE POISSON : q - div(eps*E) = 0
    */
    K(E_poi+i*NEQ,I_c_c+i*NEQ) +=   volume[i]*F*z*phi ;
    K(E_poi+i*NEQ,I_c_a+i*NEQ) += - volume[i]*F*z*phi ;
  }
  /*
    termes d'ecoulement
  */
  trd_c  = dt*surf/dx*d_c ;
  trd_a  = dt*surf/dx*d_a ;
  trm_c  = dt*surf/dx*K_c ;
  trm_a  = dt*surf/dx*K_a ;
  tr_e   = surf/dx*eps ;
  /*
    CONSERVATION DE +   : (n_c1 - n_cn) + dt * div(J_c1) = 0
  */
  K(E_cat,I_c_c)          += + trd_c ;
  K(E_cat,I_c_c+NEQ)      += - trd_c ;
  K(E_cat+NEQ,I_c_c)      += - trd_c ;
  K(E_cat+NEQ,I_c_c+NEQ)  += + trd_c ;

  K(E_cat,I_psi)          += + trm_c ;
  K(E_cat,I_psi+NEQ)      += - trm_c ;
  K(E_cat+NEQ,I_psi)      += - trm_c ;
  K(E_cat+NEQ,I_psi+NEQ)  += + trm_c ;
  /*
    CONSERVATION DE -   : (n_a1 - n_an) + dt * div(J_a1) = 0
  */
  K(E_ani,I_c_a)          += + trd_a ;
  K(E_ani,I_c_a+NEQ)      += - trd_a ;
  K(E_ani+NEQ,I_c_a)      += - trd_a ;
  K(E_ani+NEQ,I_c_a+NEQ)  += + trd_a ;

  K(E_ani,I_psi)          += + trm_a ;
  K(E_ani,I_psi+NEQ)      += - trm_a ;
  K(E_ani+NEQ,I_psi)      += - trm_a ;
  K(E_ani+NEQ,I_psi+NEQ)  += + trm_a ;
  /*
    EQUATION DE POISSON : q - div(eps*E) = 0
  */
  K(E_poi,I_psi)          += - tr_e ;
  K(E_poi,I_psi+NEQ)      += + tr_e ;
  K(E_poi+NEQ,I_psi)      += + tr_e ;
  K(E_poi+NEQ,I_psi+NEQ)  += - tr_e ;

  return(0) ;

#undef K_c
#undef K_a
#undef K
}

void rs11(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define N_c1(n)   (f_1[(n)])
#define N_a1(n)   (f_1[(n+2)])
#define Q_1(n)    (f_1[(n+4)])
#define J_c1      (f_1[(6)])
#define J_a1      (f_1[(7)])
#define E_1       (f_1[(8)])
#define N_cn(n)   (f_n[(n)])
#define N_an(n)   (f_n[(n+2)])
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
    CONSERVATION DE +   : (n_c1 - n_cn) + dt * div(J_c1) = 0
  */
  R(0,E_cat) -= volume[0]*(N_c1(0) - N_cn(0)) + dt*surf*J_c1 ;
  R(1,E_cat) -= volume[1]*(N_c1(1) - N_cn(1)) - dt*surf*J_c1 ;
  /*
    CONSERVATION DE -   : (n_a1 - n_an) + dt * div(J_a1) = 0
  */
  R(0,E_ani) -= volume[0]*(N_a1(0) - N_an(0)) + dt*surf*J_a1 ;
  R(1,E_ani) -= volume[1]*(N_a1(1) - N_an(1)) - dt*surf*J_a1 ;
  /*
    EQUATION DE POISSON : q - div(eps*E) = 0
  */
  R(0,E_poi) -= volume[0]*Q_1(0) - surf*eps*E_1 ;
  R(1,E_poi) -= volume[1]*Q_1(1) + surf*eps*E_1 ;
  
#undef R
#undef N_c1
#undef N_cn
#undef N_a1
#undef N_an
#undef J_c1
#undef J_a1
#undef E_1
}

int so11(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define J_c      (f[(6)])
#define J_a      (f[(7)])
#define E        (f[(8)])
  double c_c,c_a,psi,q,j_e ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  eps     = el.mat->pr[pm("permittivite")] ;
  d_c     = el.mat->pr[pm("D_c")] ;
  d_a     = el.mat->pr[pm("D_a")] ;
  RT      = el.mat->pr[pm("RT")] ;
  F       = el.mat->pr[pm("F")] ;
  z       = el.mat->pr[pm("valence")] ;

  /* initialisation */
  nso = 8 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  /* concentrations */
  c_c  = param(u,h_s,el.nn,I_c_c) ;
  c_a  = param(u,h_s,el.nn,I_c_a) ;
  psi  = param(u,h_s,el.nn,I_psi) ;
  /* charge */
  q    = F*z*phi*(c_c - c_a) ;
  /* courant */
  j_e  = F*z*(J_c - J_a) ;
  /* quantites exploitees */
  strcpy(r[0].text,"concentration_c") ; r[0].n = 1 ;
  r[0].v[0] = c_c ;
  strcpy(r[1].text,"concentration_a") ; r[1].n = 1 ;
  r[1].v[0] = c_a ;
  strcpy(r[2].text,"potentiel") ; r[2].n = 1 ;
  r[2].v[0] = psi ;
  strcpy(r[3].text,"flux_c") ; r[3].n = 1 ;
  r[3].v[0] = J_c ;
  strcpy(r[4].text,"flux_a") ; r[4].n = 1 ;
  r[4].v[0] = J_a ;
  strcpy(r[5].text,"champ_E") ; r[5].n = 1 ;
  r[5].v[0] = E ;
  strcpy(r[6].text,"charge") ; r[6].n = 1 ;
  r[6].v[0] = q ;
  strcpy(r[7].text,"courant") ; r[7].n = 1 ;
  r[7].v[0] = j_e ;
  return (nso) ;
  
#undef J_c
#undef J_a
#undef E
}
