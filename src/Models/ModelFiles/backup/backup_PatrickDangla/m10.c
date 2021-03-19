#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  10
#define TITLE "Equation de Richards (3D)"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ (1)
#define NVE (1+dim)
#define NVA (1)
/* Fonctions */
static int    pm(char *s) ;
static int    c10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c) ;
static int    k10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c) ;
/* Parametres */
static double gravite,phi,rho_l,k_int,mu_l,p_g ;

int pm(char *s) {
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"p_g") == 0) return (5) ;
  else if(strcmp(s,"courbes") == 0) return (6) ;
  else
    { printf("donnee \"%s\" non connue (pm10)\n",s) ; exit(0) ; }
}

int dm10(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 7 ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[0],"liq") ;
  strcpy(mat->inc[0],"p_l") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}

int qm10(int dim,FILE *ficd)
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

void tb10(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = (1+dim)*el->fi->np ;
  el->n_ve = el->fi->np ;
}


void ch10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;

}


void in10(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double pc,pl,sl,gpl[3],w_l[3],k_l,m_l ;
  int    i,p ;
  double *h,*dh ;

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
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,0) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* masse liquide */
    m_l = rho_l*phi*sl ;
    /* coefficient de transfert */
    k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    /* flux */
    grad(x,u,dh,gpl,el.nn,dim,0) ;
    for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
    w_l[dim-1] += k_l*rho_l*gravite ;
    /* rangement dans f */
    f[p*NVE] = m_l ;
    for(i=0;i<dim;i++) f[p*NVE+1+i]  = w_l[i] ;
    /* rangement dans va */
    va[p*NVA] = k_l ;
  }
}


int ex10(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  int    p ;
  double pl,pc,k_l ;
  double *h ;

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
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,0) ;
    pc  = p_g - pl ;
    /* coefficient de transfert */
    k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    /* rangement dans va */
    va[p*NVA] = k_l ;
  }
  return(0) ;
}

int ct10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double pc,pl,sl,gpl[3],w_l[3],k_l,m_l ;
  int    i,p ;
  double *h,*dh ;

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
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    pl  = param(u_1,h,el.nn,0) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* masse liquide */
    m_l = rho_l*phi*sl ;
    /* coefficient de transfert */
    k_l = va[p*NVA] ;
    /* flux */
    grad(x,u_1,dh,gpl,el.nn,dim,0) ;
    for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
    w_l[dim-1] += k_l*rho_l*gravite ;
    /* rangement dans f_1 */
    f_1[p*NVE] = m_l ;
    for(i=0;i<dim;i++) f_1[p*NVE+1+i]  = w_l[i] ;
  }

  return(0) ;
}


int mx10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
  int    i,dec ;
  double c[MAX_PGAUSS*9] ;
  double kb[MAX_NOEUDS*MAX_NOEUDS] ;

  /* initialisation */
  for(i=0;i<el.nn*el.nn;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;
  
  dec = c10(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxmass(k,x,*el.fi,c,dim,dec,geom) ;
  dec = k10(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxcond(kb,x,*el.fi,c,dim,dec,geom) ;
  for(i=0;i<el.nn*el.nn;i++) k[i] += dt*kb[i] ;

  return(0) ;
}


void rs10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;
  double rb[MAX_NOEUDS],g1[MAX_NOEUDS] ;

  /* initialisation */
  for(i=0;i<el.nn;i++) r[i] = 0. ;

  if(el.dim < dim) return ;
  
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVE] - f_n[i*NVE] ;
  rsmass(r,x,*el.fi,g1,dim,1,geom) ;
  rsflux(rb,x,*el.fi,f_1+1,dim,NVE,geom) ;
  for(i=0;i<el.nn;i++) r[i] = -r[i] + dt*rb[i] ;
}


int so10(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double pc,pl,sl,gpl[3],w_l[3],k_l ;
  int    i,j,p,nso ;
  double *h,*dh,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
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

  /* pression */
  pl = param(u,h_s,el.nn,0) ;
  pc = p_g - pl ;
  /* saturation */
  sl  = courbe(pc,el.mat->cb[0]) ;

  /* quantites exploitees */
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[2].text,"saturation") ; r[2].n = 1 ;
  r[2].v[0] = sl ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* coefficient de transfert */
    k_l = va[p*NVA] ;
    /* flux */
    grad(x,u,dh,gpl,el.nn,dim,0) ;
    for(i=0;i<dim;i++) w_l[i] = - k_l*gpl[i] ;
    w_l[dim-1] += k_l*rho_l*gravite ;
    /* quantites exploitees par element */
    strcpy(r[1].text,"flux-liquide") ; r[1].n = 3 ;
    for(i=0;i<dim;i++) r[1].v[i] += w_l[i]/el.fi->np ;
  }
  return (nso) ;
}

int c10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  double pl,pc,sl,dslsdpc ;
  int    dec,p ;
  double *h ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  phi     = el.mat->pr[pm("phi")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  
  dec = 1 ;
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u_1,h,el.nn,0) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    dslsdpc = dcourbe(pc,el.mat->cb[0]) ;
    c[p] = -rho_l*phi*dslsdpc ;
  }
  return(dec) ;
}


int k10(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec ;
  double *c1 ;
  int    i,p ;
  double zero = 0. ;
  
  dec = 9 ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* permeabilite */
    c1[0] = va[p*NVA] ;
    c1[4] = va[p*NVA] ;
    c1[8] = va[p*NVA] ;
  }
  return(dec) ;
}
#undef NEQ
#undef NVE
#undef NVA
