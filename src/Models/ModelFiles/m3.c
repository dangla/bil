#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  3
#define TITLE "Equation de Poisson-Boltzmann"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ (1)
#define NVE (2)
#define NVA (0)
/* Fonctions */
static int    pm(const char *s) ;
/* Parametres */
static double sigma,e,v,epsilon,T_0,n_0,k_B ;
static int    ne ;

int pm(const char *s)
{
  if(strcmp(s,"sigma") == 0) return (0) ;
  else if(strcmp(s,"e") == 0) return (1) ;
  else if(strcmp(s,"v") == 0) return (2) ;
  else if(strcmp(s,"epsilon") == 0) return (3) ;
  else if(strcmp(s,"T_0") == 0) return (4) ;
  else if(strcmp(s,"n_0") == 0) return (5) ;
  else if(strcmp(s,"k_B") == 0) return (6) ;
  else if(strcmp(s,"n_1") == 0) return (7) ;
  else if(strcmp(s,"ne") == 0) return (8) ;
  else
    { printf("donnee \"%s\" non connue (pm3)\n",s) ; exit(0) ; }
}

int dm3(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 9 ;
  
  mat->neq    = 1 ;
  strcpy(mat->eqn[0],"eq_PB") ;
  strcpy(mat->inc[0],"psi") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm3(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Equation de Poisson-Boltzmann, inconnue le potentiel electrique (psi)\n") ;

  fprintf(stdout,"\n\
Exemple de donnees :\n\n") ;


  fprintf(ficd,"sigma = -0.2         # Charge de surface\n") ;
  fprintf(ficd,"e = 1.6e-19          # Charge electrique\n") ;
  fprintf(ficd,"v = 2                # Valence\n") ;
  fprintf(ficd,"epsilon = 7.0832e-10 # Constante dielectrique\n") ;
  fprintf(ficd,"k_B = 1.38e-23       # Constante de Boltzmann\n") ;
  fprintf(ficd,"T_0 = 300            # Temperature\n") ;
  fprintf(ficd,"n_0 = 6.022e26       # Concentration ionique initiale\n") ;
  fprintf(ficd,"n_1 = 0              # Concentration ionique n(t) = n_0 + n_1*t\n") ;
  fprintf(ficd,"ne = 100             # Nb d'elements du maillage\n") ;
  
  return(1) ;
}

void tb3(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVE ;
  el->n_ve = NVA ;
}

void ch3(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  /*
    Donnees
  */
  sigma   = el.mat->pr[pm("sigma")] ;
  e       = el.mat->pr[pm("e")] ;
  v       = el.mat->pr[pm("v")] ;
  epsilon = el.mat->pr[pm("epsilon")] ;
  T_0     = el.mat->pr[pm("T_0")] ;
  k_B     = el.mat->pr[pm("k_B")] ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  r[0] *= v*e*sigma/(epsilon*k_B*T_0) ;
}

void in3(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define Psi(n)   (f[(n)])
  int    i ;
  
  if(el.dim < dim) return ;

  /* potentiels */
  for(i=0;i<2;i++) {
    Psi(i) = u[i][0] ;
  }
#undef Psi
}

int ex3(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  return(0) ;
}

int ct3(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define Psi(n)   (f_1[(n)])
  int    i ;
  
  if(el.dim < dim) return(0) ;

  /* potentiels */
  for(i=0;i<2;i++) {
    Psi(i) = u_1[i][0] ;
  }

  return(0) ;
#undef Psi
}

int mx3(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2+(j)])
  double ul_d2 ;
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
  sigma   = el.mat->pr[pm("sigma")] ;
  e       = el.mat->pr[pm("e")] ;
  v       = el.mat->pr[pm("v")] ;
  epsilon = el.mat->pr[pm("epsilon")] ;
  n_0     = el.mat->pr[pm("n_0")] + el.mat->pr[pm("n_1")]*t ;
  T_0     = el.mat->pr[pm("T_0")] ;
  k_B     = el.mat->pr[pm("k_B")] ;
  ul_d2   = 2.*v*v*e*e*n_0/(epsilon*k_B*T_0) ;
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
    POISSON-BOLTZMANN : 1/L^2sinh(u) - div(grad(u)) = 0
  */
  K(0,0) = volume[0]*ul_d2*cosh(u_1[0][0]) + surf/dx  ;
  K(0,1) = - surf/dx  ;
  K(1,0) = - surf/dx  ;
  K(1,1) = volume[1]*ul_d2*cosh(u_1[1][0]) + surf/dx  ;
  
  return(0) ;
#undef K
}

void rs3(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  double ul_d2 ;
  double dx, xm ;
  double volume[2], surf ;
  int    i ;
  double zero = 0., un = 1., deux = 2. ;                                          
                                                                                
  /* initialisation */
  for(i=0;i<el.nn;i++) r[i] = zero ;
                                                                                
  if(el.dim < dim) return ;

  /*
    Donnees
  */
  sigma   = el.mat->pr[pm("sigma")] ;
  e       = el.mat->pr[pm("e")] ;
  v       = el.mat->pr[pm("v")] ;
  epsilon = el.mat->pr[pm("epsilon")] ;
  n_0     = el.mat->pr[pm("n_0")] + el.mat->pr[pm("n_1")]*t ;
  T_0     = el.mat->pr[pm("T_0")] ;
  k_B     = el.mat->pr[pm("k_B")] ;
  ul_d2   = 2.*v*v*e*e*n_0/(epsilon*k_B*T_0) ;
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
  for(i=0;i<2;i++) r[i] = zero ;
  /*
    POISSON-BOLTZMANN : 1/L^2sinh(u) - div(grad(u)) = 0
   */
  r[0] -= volume[0]*ul_d2*sinh(u_1[0][0]) - surf*(u_1[1][0]-u_1[0][0])/dx ;
  r[1] -= volume[1]*ul_d2*sinh(u_1[1][0]) + surf*(u_1[1][0]-u_1[0][0])/dx ;
}

int so3(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define Psi(i,n)    (f[(i)*NVE+n])
  double phi,n_plus,n_moins,lam_plus,lam_moins,som ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  int    i,j,nso ; 
 
  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  e       = el.mat->pr[pm("e")] ;
  v       = el.mat->pr[pm("v")] ;
  n_0     = el.mat->pr[pm("n_0")] + el.mat->pr[pm("n_1")]*t ;
  ne      = floor(el.mat->pr[pm("ne")] + 0.5) ;

  /* initialisation */
  nso = 7 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = 0. ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* le potentiel electrique */
  phi = param(u,h_s,el.nn,0) ;
  /* concentrations ioniques */
  n_plus  = n_0*exp(-phi) ;
  n_moins = n_0*exp(+phi) ;
  /* inverse de l'activite */
  lam_moins = 0. ;
  lam_plus  = 0. ;
  som       = 0. ;
  if(x[0][0] == 0. && t > 0.) {
    for(i=0;i<ne;i++) {
      lam_moins += (exp(+Psi(i,0))+exp(+Psi(i,1)))/2.*(x[0][i+1] - x[0][i]) ;
      lam_plus  += (exp(-Psi(i,0))+exp(-Psi(i,1)))/2.*(x[0][i+1] - x[0][i]) ;
      som       += (Psi(i,0) + Psi(i,1))/2.*(x[0][i+1] - x[0][i]) ;
    }
    lam_moins /= (x[0][ne] - x[0][0]) ;
    lam_plus  /= (x[0][ne] - x[0][0]) ;
    som       -= Psi(ne-1,1)*(x[0][ne] - x[0][0]) ;
  }
  /* quantites exploitees */
  strcpy(r[0].text,"potentiel") ; r[0].n = 1 ;
  r[0].v[0] = phi ;
  strcpy(r[1].text,"n+") ; r[1].n = 1 ;
  r[1].v[0] = n_plus ;
  strcpy(r[2].text,"n-") ; r[2].n = 1 ;
  r[2].v[0] = n_moins ;
  strcpy(r[3].text,"lambda-") ; r[3].n = 1 ;
  r[3].v[0] = lam_moins ;
  strcpy(r[4].text,"somme") ; r[4].n = 1 ;
  r[4].v[0] = som ;
  strcpy(r[5].text,"lambda+") ; r[5].n = 1 ;
  r[5].v[0] = lam_plus ;
  strcpy(r[6].text,"charge") ; r[6].n = 1 ;
  r[6].v[0] = e*v*(n_plus - n_moins) ;
  return(nso) ;
#undef Psi
}
