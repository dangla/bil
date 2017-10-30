#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  30
#define TITLE "Ecoulement a double permeabilite"
#define AUTHORS "Lassabatere"

#include "OldMethods.h"

/* Macros */
#define NEQ   (2)
#define E_fra (0)
#define E_mat (1)
#define I_p_f (0)
#define I_p_m (1)
#define NVI   (3+2*dim)
#define NVE   (3)
/* Fonctions */
static int    pm(const char*) ;
static int    c30(double**,double**,double**,double*,double*,double*,elem_t,geom_t,double,double*) ;
static int    k30(double**,double**,double**,double*,double*,double*,elem_t,geom_t,double,double*) ;
/*
static double (*xcourbe[2])(double,crbe_t) = {courbe,courbe} ;
static double (*xdcourbe[2])(double,crbe_t) = {dcourbe,dcourbe} ;
*/
/* Parametres */
static double gravite,v_f,v_m,rho,eps_f0,eps_m0,K_f,K_m,K_a,mu,p_g,S_s,p_f0,p_m0 ;

/* Macros pour les fonctions courbe/courbelog et dcourbe/dcourbelog */
#define LOG    1      /* 1 = echelle log, 0 = echelle normale */
#define COURBE(a,b)   courbe(a,b)
#define DCOURBE(a,b)  dcourbe(a,b)

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"w_f") == 0) return (1) ;
  else if(strcmp(s,"rho") == 0) return (2) ;
  else if(strcmp(s,"eps_f") == 0) return (3) ;
  else if(strcmp(s,"eps_m") == 0) return (4) ;
  else if(strcmp(s,"K_f") == 0) return (5) ;
  else if(strcmp(s,"K_m") == 0) return (6) ;
  else if(strcmp(s,"K_a") == 0) return (7) ;
  else if(strcmp(s,"mu") == 0) return (8) ;
  else if(strcmp(s,"p_g") == 0) return (9) ;
  else if(strcmp(s,"S_s") == 0) return (10) ;
  else if(strcmp(s,"p_f0") == 0) return (11) ;
  else if(strcmp(s,"p_m0") == 0) return (12) ;
  else
    { printf("donnee \"%s\" non connue (pm30)\n",s) ; exit(0) ; }
}


int dm30(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 13,n_courbes = 5 ;
  
  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_fra],"fra") ;
  strcpy(mat->eqn[E_mat],"mat") ;
  strcpy(mat->inc[I_p_f],"p_f") ;
  strcpy(mat->inc[I_p_m],"p_m") ;

  lit_mate(mat,ficd,pm,n_donnees,n_courbes) ;
  return(mat->n) ;
}

int qm30(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 2 equations:\n\
\t- conservation de l\'eau du milieu fracture  (p_f)\n\
\t- conservation de l\'eau du milieu matrice   (p_m)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = -9.81.   # gravite\n") ;
  fprintf(ficd,"rho = 1000         # masse volumique\n") ;
  fprintf(ficd,"w_f = 0.05         # frac. volumique du milieu fracture\n") ;
  fprintf(ficd,"eps_f = 0.5        # porosite du milieu fracture\n") ;
  fprintf(ficd,"eps_m = 0.5        # porosite du milieu matrice\n") ;
  fprintf(ficd,"K_f = 2.36e-11     # permeabilite du milieu fracture\n") ;
  fprintf(ficd,"K_m = 1.24e-11     # permeabilite du milieu matrice\n") ;
  fprintf(ficd,"K_a = 1.416e-12    # permeabilite de l\'interface f/m\n") ;
  fprintf(ficd,"mu = 1.e-3         # viscosite\n") ;
  fprintf(ficd,"p_g = 0            # pression du gaz\n") ;
  fprintf(ficd,"p_f0 = -1.e-4      # pression de reference du m. f.\n") ;
  fprintf(ficd,"p_m0 = -1.e-4      # pression de reference du m. m.\n") ;
  fprintf(ficd,"S_s = 1.02e-9      # souplesse elastique des pores\n") ;
  fprintf(ficd,"Courbes = my_file  # Nom du fichier : p_c S_l_f k_rl_f S_l_m k_rl_m k_rl_int\n") ;
  return(NEQ) ;
}

void tb30(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}

void ch30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = - r[i] ;
}

void in30(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double p_f,p_m,s_f,s_m,m_f,m_m,dm_a ;
  double gp_f[3],gp_m[3],w_f[3],w_m[3],k_f,k_m,k_a ;
  double eps_f,eps_m ;
  int    i,p ;
  double *h,*dh ;
  crbe_t *cb = el.mat->cb ;

  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho     = el.mat->pr[pm("rho")] ;
  v_f     = el.mat->pr[pm("w_f")] ;
  eps_f0  = el.mat->pr[pm("eps_f")] ;
  eps_m0  = el.mat->pr[pm("eps_m")] ;
  K_f     = el.mat->pr[pm("K_f")] ;
  K_m     = el.mat->pr[pm("K_m")] ;
  K_a     = el.mat->pr[pm("K_a")] ;
  mu      = el.mat->pr[pm("mu")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  S_s     = el.mat->pr[pm("S_s")] ;
  p_f0    = el.mat->pr[pm("p_f0")] ;
  p_m0    = el.mat->pr[pm("p_m0")] ;

  v_m = 1. - v_f ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    p_f  = param(u,h,el.nn,I_p_f) ;
    p_m  = param(u,h,el.nn,I_p_m) ;
    /* saturations */
    s_f = COURBE(p_g - p_f,cb[0]) ;
    s_m = COURBE(p_g - p_m,cb[2]) ;
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* masses liquides */
    m_f = rho*v_f*eps_f*s_f ;
    m_m = rho*v_m*eps_m*s_m ;
    /* coefficient de transfert */
    k_f = rho*v_f*K_f/mu*COURBE(p_g - p_f,cb[1]) ;
    k_m = rho*v_m*K_m/mu*COURBE(p_g - p_m,cb[3]) ;
    k_a = rho*K_a/mu*0.5*(COURBE(p_g - p_f,cb[4]) + COURBE(p_g - p_m,cb[4])) ;
    /* flux de masse */
    grad(x,u,dh,gp_f,el.nn,dim,I_p_f) ;
    for(i=0;i<3;i++) w_f[i] = - k_f*gp_f[i] ;
    w_f[dim-1] += k_f*rho*gravite ;
    grad(x,u,dh,gp_m,el.nn,dim,I_p_m) ;
    for(i=0;i<3;i++) w_m[i] = - k_m*gp_m[i] ;
    w_m[dim-1] += k_m*rho*gravite ;
    /* terme d'echange */
    dm_a = k_a*(p_f - p_m) ;
    /* rangement dans f */
    f[p*NVI+E_fra] = m_f ;
    f[p*NVI+E_mat] = m_m ;
    for(i=0;i<dim;i++) {
      f[p*NVI+NEQ+E_fra*dim+i]  = w_f[i] ;
      f[p*NVI+NEQ+E_mat*dim+i]  = w_m[i] ;
    }
    f[p*NVI+NEQ+NEQ*dim] = dm_a ;
    /* rangement dans va */
    va[p*NVE]   = k_f ;
    va[p*NVE+1] = k_m ;
    va[p*NVE+2] = k_a ;
  }
}

int ex30(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double p_f,p_m,k_f,k_m,k_a ;
  int    p ;
  double *h ;
  crbe_t *cb = el.mat->cb ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho     = el.mat->pr[pm("rho")] ;
  v_f     = el.mat->pr[pm("w_f")] ;
  eps_f0  = el.mat->pr[pm("eps_f")] ;
  eps_m0  = el.mat->pr[pm("eps_m")] ;
  K_f     = el.mat->pr[pm("K_f")] ;
  K_m     = el.mat->pr[pm("K_m")] ;
  K_a     = el.mat->pr[pm("K_a")] ;
  mu      = el.mat->pr[pm("mu")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  S_s     = el.mat->pr[pm("S_s")] ;
  p_f0    = el.mat->pr[pm("p_f0")] ;
  p_m0    = el.mat->pr[pm("p_m0")] ;

  v_m = 1. - v_f ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pressions */
    p_f  = param(u,h,el.nn,I_p_f) ;
    p_m  = param(u,h,el.nn,I_p_m) ;
    /* coefficient de transfert */
    k_f = rho*v_f*K_f/mu*COURBE(p_g - p_f,cb[1]) ;
    k_m = rho*v_m*K_m/mu*COURBE(p_g - p_m,cb[3]) ;
    k_a = rho*K_a/mu*0.5*(COURBE(p_g - p_f,cb[4]) + COURBE(p_g - p_m,cb[4])) ;
    /* rangement dans va */
    va[p*NVE]   = k_f ;
    va[p*NVE+1] = k_m ;
    va[p*NVE+2] = k_a ;
  }
  return(0) ;
}

int ct30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double p_f,p_m,s_f,s_m,m_f,m_m,dm_a ;
  double gp_f[3],gp_m[3],w_f[3],w_m[3],k_f,k_m,k_a ;
  double eps_f,eps_m ;
  int    i,p ;
  double *h,*dh ;
  crbe_t *cb = el.mat->cb ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho     = el.mat->pr[pm("rho")] ;
  v_f     = el.mat->pr[pm("w_f")] ;
  eps_f0  = el.mat->pr[pm("eps_f")] ;
  eps_m0  = el.mat->pr[pm("eps_m")] ;
  K_f     = el.mat->pr[pm("K_f")] ;
  K_m     = el.mat->pr[pm("K_m")] ;
  K_a     = el.mat->pr[pm("K_a")] ;
  mu      = el.mat->pr[pm("mu")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  S_s     = el.mat->pr[pm("S_s")] ;
  p_f0    = el.mat->pr[pm("p_f0")] ;
  p_m0    = el.mat->pr[pm("p_m0")] ;

  v_m = 1. - v_f ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    p_f  = param(u_1,h,el.nn,I_p_f) ;
    p_m  = param(u_1,h,el.nn,I_p_m) ;
    /* saturations */
    s_f = COURBE(p_g - p_f,cb[0]) ;
    s_m = COURBE(p_g - p_m,cb[2]) ;
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* masses liquides */
    m_f = rho*v_f*eps_f*s_f ;
    m_m = rho*v_m*eps_m*s_m ;
    /* coefficient de transfert */
    k_f = va[p*NVE] ;
    k_m = va[p*NVE+1] ;
    k_a = va[p*NVE+2] ;
    /* flux de masse */
    grad(x,u_1,dh,gp_f,el.nn,dim,I_p_f) ;
    for(i=0;i<3;i++) w_f[i] = - k_f*gp_f[i] ;
    w_f[dim-1] += k_f*rho*gravite ;
    grad(x,u_1,dh,gp_m,el.nn,dim,I_p_m) ;
    for(i=0;i<3;i++) w_m[i] = - k_m*gp_m[i] ;
    w_m[dim-1] += k_m*rho*gravite ;
    /* terme d'echange */
    dm_a = k_a*(p_f - p_m) ;
    /* rangement dans f */
    f_1[p*NVI+E_fra] = m_f ;
    f_1[p*NVI+E_mat] = m_m ;
    for(i=0;i<dim;i++) {
      f_1[p*NVI+NEQ+E_fra*dim+i]  = w_f[i] ;
      f_1[p*NVI+NEQ+E_mat*dim+i]  = w_m[i] ;
    }
    f_1[p*NVI+NEQ+NEQ*dim] = dm_a ;
  }

  return(0) ;
}

int  mx30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
  int    i,dec ;
  double c[MAX_PGAUSS*9*NEQ*NEQ] ;
  double kb[NEQ*NEQ*MAX_NOEUDS*MAX_NOEUDS] ;

  /* initialisation */
  for(i=0;i<NEQ*NEQ*el.nn*el.nn;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    loi de comportement
  */
  dec = c30(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxcmss(k,x,*el.fi,c,dim,dec,geom,NEQ) ;
  /*
    conduction
  */
  dec = k30(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxccnd(kb,x,*el.fi,c,dim,dec,geom,NEQ) ;
  for(i=0;i<NEQ*NEQ*el.nn*el.nn;i++) k[i] += dt*kb[i] ;

  return(0) ;
}

void rs30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;
  double rb[MAX_NOEUDS],g1[MAX_PGAUSS],*f_2 ;

  /* initialisation */
  for(i=0;i<NEQ*el.nn;i++) r[i] = 0. ;

  if(el.dim < dim) return ;
  
  /* fracture */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_fra] - f_n[i*NVI+E_fra] + dt*f_1[i*NVI+NEQ+NEQ*dim] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fra] = - rb[i] ;
  f_2 = f_1 + NEQ + E_fra*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fra] +=  dt*rb[i] ;

  /* matrice */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_mat] - f_n[i*NVI+E_mat] - dt*f_1[i*NVI+NEQ+NEQ*dim] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_mat] = - rb[i] ;
  f_2 = f_1 + NEQ + E_mat*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_mat] +=  dt*rb[i] ;
}

int so30(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_f,p_m,s_f,s_m,gp_f[3],gp_m[3],w_f[3],w_m[3],k_f,k_m ;
  int    i,j,p,nso ;
  double *h,*dh,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  crbe_t *cb = el.mat->cb ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho     = el.mat->pr[pm("rho")] ;
  v_f     = el.mat->pr[pm("w_f")] ;
  eps_f0  = el.mat->pr[pm("eps_f")] ;
  eps_m0  = el.mat->pr[pm("eps_m")] ;
  K_f     = el.mat->pr[pm("K_f")] ;
  K_m     = el.mat->pr[pm("K_m")] ;
  K_a     = el.mat->pr[pm("K_a")] ;
  mu      = el.mat->pr[pm("mu")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  S_s     = el.mat->pr[pm("S_s")] ;
  p_f0    = el.mat->pr[pm("p_f0")] ;
  p_m0    = el.mat->pr[pm("p_m0")] ;

  v_m = 1. - v_f ;
  
  /* initialisation */
  nso = 6 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pressions */
  p_f  = param(u,h_s,el.nn,I_p_f) ;
  p_m  = param(u,h_s,el.nn,I_p_m) ;
  /* saturations */
  s_f = COURBE(p_g - p_f,cb[0]) ;
  s_m = COURBE(p_g - p_m,cb[2]) ;

  /* quantites exploitees */
  strcpy(r[0].text,"pression-fracture") ; r[0].n = 1 ;
  r[0].v[0] = p_f ;
  strcpy(r[1].text,"pression-matrice") ; r[1].n = 1 ;
  r[1].v[0] = p_m ;
  strcpy(r[2].text,"saturation-fracture") ; r[2].n = 1 ;
  r[2].v[0] = s_f ;
  strcpy(r[3].text,"saturation-matrice") ; r[3].n = 1 ;
  r[3].v[0] = s_m ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* coefficient de transfert */
    k_f = va[p*NVE] ;
    k_m = va[p*NVE+1] ;
    /* flux de masse */
    grad(x,u,dh,gp_f,el.nn,dim,I_p_f) ;
    for(i=0;i<3;i++) w_f[i] = - k_f*gp_f[i] ;
    w_f[dim-1] += k_f*rho*gravite ;
    grad(x,u,dh,gp_m,el.nn,dim,I_p_m) ;
    for(i=0;i<3;i++) w_m[i] = - k_m*gp_m[i] ;
    w_m[dim-1] += k_m*rho*gravite ;
    /* quantites exploitees par element */
    strcpy(r[4].text,"flux_fracture") ; r[4].n = 3 ;
    for(i=0;i<dim;i++) r[4].v[i] += w_f[i]/el.fi->np ;
    strcpy(r[5].text,"flux_matrice") ; r[5].n = 3 ;
    for(i=0;i<dim;i++) r[5].v[i] += w_m[i]/el.fi->np ;
  }
  return (nso) ;
}


int c30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  double p_f,p_m,k_a ;
  double s_f,s_m,ds_fsdpc,ds_msdpc ;
  double eps_f,eps_m ;
  crbe_t *cb = el.mat->cb ;
  int    dec,p ;
  double *h,*c1 ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho     = el.mat->pr[pm("rho")] ;
  v_f     = el.mat->pr[pm("w_f")] ;
  eps_f0  = el.mat->pr[pm("eps_f")] ;
  eps_m0  = el.mat->pr[pm("eps_m")] ;
  K_f     = el.mat->pr[pm("K_f")] ;
  K_m     = el.mat->pr[pm("K_m")] ;
  K_a     = el.mat->pr[pm("K_a")] ;
  mu      = el.mat->pr[pm("mu")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  S_s     = el.mat->pr[pm("S_s")] ;
  p_f0    = el.mat->pr[pm("p_f0")] ;
  p_m0    = el.mat->pr[pm("p_m0")] ;

  v_m = 1. - v_f ;

  dec = NEQ*NEQ ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pressions */
    p_f  = param(u_1,h,el.nn,I_p_f) ;
    p_m  = param(u_1,h,el.nn,I_p_m) ;
    /* saturations */
    s_f = COURBE(p_g - p_f,cb[0]) ;
    s_m = COURBE(p_g - p_m,cb[2]) ;
    /* derivees des saturations */
    ds_fsdpc = DCOURBE(p_g - p_f,cb[0]) ;
    ds_msdpc = DCOURBE(p_g - p_m,cb[2]) ;
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* coefficient de transfert */
    k_a = va[p*NVE+2] ;
    /* fracture */
    c1[E_fra*NEQ+I_p_f] =   rho*v_f*(eps_f*(-ds_fsdpc) + S_s*s_f) + dt*k_a ;
    c1[E_fra*NEQ+I_p_m] = - dt*k_a ;
    /* matrice */
    c1[E_mat*NEQ+I_p_f] = - dt*k_a ;
    c1[E_mat*NEQ+I_p_m] =   rho*v_m*(eps_m*(-ds_msdpc) + S_s*s_m) + dt*k_a ;
  }
  return(dec) ;
}


int k30(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec ;
  double *c1 ;
  int    i,p ;
  double zero = 0. ;
  
  dec = 9*NEQ*NEQ ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* permeabilite fracture */
    c1[(E_fra*NEQ+I_p_f)*9+0] = va[p*NVE] ;
    c1[(E_fra*NEQ+I_p_f)*9+4] = va[p*NVE] ;
    c1[(E_fra*NEQ+I_p_f)*9+8] = va[p*NVE] ;
    /* permeabilite matrice */
    c1[(E_mat*NEQ+I_p_m)*9+0] = va[p*NVE+1] ;
    c1[(E_mat*NEQ+I_p_m)*9+4] = va[p*NVE+1] ;
    c1[(E_mat*NEQ+I_p_m)*9+8] = va[p*NVE+1] ;
  }
  return(dec) ;
}

#undef NEQ
#undef NVI
#undef NVE

#undef COURBE
#undef DCOURBE
