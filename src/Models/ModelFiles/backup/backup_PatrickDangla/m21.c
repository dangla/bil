#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  21
#define TITLE "Electrostatique"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ (1)
#define NVE (dim)
#define NVA (0)
/* Fonctions */
static int    pm(const char *s) ;
/* Parametres */
static double eps ;

int pm(const char *s) {
  if(strcmp(s,"permittivite") == 0) return (0) ;
  else
    { printf("donnee \"%s\" non connue (pm21)\n",s) ; exit(0) ; }
}

int dm21(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 1 ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[0],"eq_Poisson") ;
  strcpy(mat->inc[0],"pot") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm21(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Equation fondamentale de l\'electrostatique (pot)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"permittivite = 7.0832e-10 # Permittivite du milieu\n") ;
  
  return(NEQ) ;
}

void tb21(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVE*el->fi->np ;
  el->n_ve = NVA*el->fi->np ;
}


void ch21(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;

}


void in21(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double grd[3] ;
  int    i,p ;
  double *dh ;

  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  eps = el.mat->pr[pm("permittivite")] ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    dh = el.fi->dh + p*dim*el.nn ;
    /* gradient */
    grad(x,u,dh,grd,el.nn,dim,0) ;
    /* rangement dans f */
    for(i=0;i<dim;i++) f[p*NVE+i]  = - eps*grd[i] ;
  }
}


int ex21(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  return(0) ;
}

int ct21(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double grd[3] ;
  int    i,p ;
  double *dh ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  eps = el.mat->pr[pm("permittivite")] ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    dh = el.fi->dh + p*dim*el.nn ;
    /* gradient */
    grad(x,u_1,dh,grd,el.nn,dim,0) ;
    /* rangement dans f_1 */
    for(i=0;i<dim;i++) f_1[p*NVE+i]  = - eps*grd[i] ;
  }

  return(0) ;
}


int  mx21(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
  int    i,dec ;
  double c[MAX_PGAUSS*9] ;

  /* initialisation */
  for(i=0;i<el.nn*el.nn;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  eps = el.mat->pr[pm("permittivite")] ;
  
  for(i=0;i<9;i++) c[i] = 0. ;
  c[0] = eps ;
  c[4] = eps ;
  c[8] = eps ;
  dec = 0 ;
  mxcond(k,x,*el.fi,c,dim,dec,geom) ;

  return(0) ;
}


void rs21(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;

  /* initialisation */
  for(i=0;i<el.nn;i++) r[i] = 0. ;

  if(el.dim < dim) return ;
  
  rsflux(r,x,*el.fi,f_1,dim,NVE,geom) ;
}


int so21(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double psi,grd[3] ;
  int    i,j,p,nso ;
  double *dh,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  eps = el.mat->pr[pm("permittivite")] ;
  
  /* initialisation */
  nso = 2 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* potentiel */
  psi = param(u,h_s,el.nn,0) ;

  /* quantites exploitees */
  strcpy(r[0].text,"potentiel") ; r[0].n = 1 ;
  r[0].v[0] = psi ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    dh = el.fi->dh + p*dim*el.nn ;
    /* gradient */
    grad(x,u,dh,grd,el.nn,dim,0) ;
    /* quantites exploitees par element */
    strcpy(r[1].text,"champ_electrique") ; r[1].n = 3 ;
    for(i=0;i<dim;i++) r[1].v[i] += - grd[i]/el.fi->np ;
  }
  return (nso) ;
}
#undef NEQ
#undef NVE
#undef NVA
