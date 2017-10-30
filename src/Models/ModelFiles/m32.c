#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  32
#define TITLE "Ecoulement a double permeabilite MIM"
#define AUTHORS "Lassabatere"

#include "OldMethods.h"

/* Macros */
#define NEQ   (4)
#define E_fra (0)
#define E_mat (1)
#define E_sf  (2)
#define E_sm  (3)
#define I_p_f (0)
#define I_p_m (1)
#define I_c_f (2)
#define I_c_m (3)
#define NVI   (8+4*dim)
#define NVE   (25)
/* Fonctions */
static int    pm(const char*) ;
static int    c32(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static int    k32(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static double tortuosite_f(double,double) ;
static double tortuosite_m(double,double) ;
static double normv(double*) ;
/* Parametres */
static double gravite,v_f,v_m,rho,eps_f0,eps_m0,K_f,K_m,K_a,mu,p_g,S_s,p_f0,p_m0,d0_sf,d0_sm,lambda_sf,lambda_sm,alpha_s,fmobile_f,fmobile_m,beta_sf,beta_sm ;

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
  else if(strcmp(s,"D_sf") == 0) return (13) ;
  else if(strcmp(s,"Lambda_sf") == 0) return (14) ;
  else if(strcmp(s,"D_sm") == 0) return (15) ;
  else if(strcmp(s,"Lambda_sm") == 0) return (16) ;
  else if(strcmp(s,"alpha_s") == 0) return (17) ;
  else if(strcmp(s,"fmobile_f") == 0) return (18) ;
  else if(strcmp(s,"fmobile_m") == 0) return (19) ;
  else if(strcmp(s,"beta_sf") == 0) return (20) ;
  else if(strcmp(s,"beta_sm") == 0) return (21) ;
  else
    { printf("donnee \"%s\" non connue (pm32)\n",s) ; exit(0) ; }
}


int dm32(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 22,n_courbes = 5 ;
  
  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_fra],"fra") ;
  strcpy(mat->eqn[E_mat],"mat") ;
  strcpy(mat->eqn[E_sf],"solute_f") ;
  strcpy(mat->eqn[E_sm],"solute_m") ;
  strcpy(mat->inc[I_p_f],"p_f") ;
  strcpy(mat->inc[I_p_m],"p_m") ;
  strcpy(mat->inc[I_c_f],"c_f") ;
  strcpy(mat->inc[I_c_m],"c_m") ;

  lit_mate(mat,ficd,pm,n_donnees,n_courbes) ;
  return(mat->n) ;
}

int qm32(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 4 equations:\n\
\t- conservation de l\'eau du milieu fracture  (p_f)\n\
\t- conservation de l\'eau du milieu matrice   (p_m)\n\
\t- conservation du solute du milieu fracture (c_f)\n\
\t- conservation du solute du milieu matrice  (c_m)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0.       # gravite\n") ;
  fprintf(ficd,"rho = 1000         # masse volumique\n") ;
  fprintf(ficd,"w_f = 0.3          # frac. volumique du milieu fracture\n") ;
  fprintf(ficd,"eps_f = 0.395      # porosite du milieu fracture\n") ;
  fprintf(ficd,"eps_m = 0.395      # porosite du milieu matrice\n") ;
  fprintf(ficd,"K_f = 2.1e-8       # permeabilite du milieu fracture\n") ;
  fprintf(ficd,"K_m = 0.5e-8       # permeabilite du milieu matrice\n") ;
  fprintf(ficd,"K_a = 0            # permeabilite de l\'interface f/m\n") ;
  fprintf(ficd,"mu = 1.e-3         # viscosite\n") ;
  fprintf(ficd,"p_g = 0            # pression du gaz\n") ;
  fprintf(ficd,"p_f0 = -1.e-4      # pression de reference du m. f.\n") ;
  fprintf(ficd,"p_m0 = -1.e-4      # pression de reference du m. m.\n") ;
  fprintf(ficd,"S_s = 0            # souplesse elastique des pores\n") ;
  fprintf(ficd,"D_sf = 0           # coef. de diffusion eff. du m. f.\n") ;
  fprintf(ficd,"D_sm = 0           # coef. de diffusion eff. du m. m.\n") ;
  fprintf(ficd,"Lambda_sf = 1.3e-2 # coef. de dispersion du m. fracture\n") ;
  fprintf(ficd,"Lambda_sm = 0.3e-2 # coef. de dispersion du m. matrice\n") ;
  fprintf(ficd,"alpha_s = 0        # coef. d\'echange d\'eau f/m\n") ;
  fprintf(ficd,"fmobile_f = 0.6    # frac. vol. d\'eau mobile du m. f.\n") ;
  fprintf(ficd,"fmobile_m = 0.6    # frac. vol. d\'eau mobile du m. m.\n") ;
  fprintf(ficd,"beta_sf = 1.e-5    # coef. d\'echange de solute MIM du m. f.\n") ;
  fprintf(ficd,"beta_sm = 1.e-5    # coef. d\'echange de solute MIM du m. m.\n") ;
  fprintf(ficd,"Courbes = my_file  # Nom du fichier : p_c S_l_f k_rl_f S_l_m k_rl_m k_rl_int\n") ;

  return(NEQ) ;
}

void tb32(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}

void ch32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = - r[i] ;
}

void in32(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double p_f,p_m,s_f,s_m,c_f,c_m,c_f_im,c_m_im ;
  double m_f,m_m,dm_a,m_sf,m_sm,dm_sa,m_sf_im,m_sm_im ;
  double gp_f[3],gp_m[3],w_f[3],w_m[3],gc_f[3],gc_m[3],w_sf[3],w_sm[3] ;
  double k_f,k_m,k_a,k_sf,k_sm,d_sf[3][3],d_sm[3][3],k_sa,d_sa ;
  double eps_f,eps_m ;
  double ww_f,ww_m ;
  int    i,j,p ;
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
  d0_sf   = el.mat->pr[pm("D_sf")] ;
  d0_sm   = el.mat->pr[pm("D_sm")] ;
  lambda_sf = el.mat->pr[pm("Lambda_sf")] ;
  lambda_sm = el.mat->pr[pm("Lambda_sm")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;
  fmobile_f = el.mat->pr[pm("fmobile_f")] ;
  fmobile_m = el.mat->pr[pm("fmobile_m")] ;
  beta_sf = el.mat->pr[pm("beta_sf")] ;
  beta_sm = el.mat->pr[pm("beta_sm")] ;

  v_m = 1. - v_f ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    p_f  = param(u,h,el.nn,I_p_f) ;
    p_m  = param(u,h,el.nn,I_p_m) ;
    /* concentrations */
    c_f  = param(u,h,el.nn,I_c_f) ;
    c_m  = param(u,h,el.nn,I_c_m) ;
    c_f_im = c_f ;
    c_m_im = c_m ;
    /* saturations */
    s_f = COURBE(p_g - p_f,cb[0]) ;
    s_m = COURBE(p_g - p_m,cb[2]) ;
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* masses */
    m_f  = rho*v_f*eps_f*s_f ;
    m_m  = rho*v_m*eps_m*s_m ;
    m_sf_im = (1. - fmobile_f)*c_f_im*v_f*eps_f*s_f ;
    m_sm_im = (1. - fmobile_m)*c_m_im*v_m*eps_m*s_m ;
    m_sf = fmobile_f*c_f*v_f*eps_f*s_f + m_sf_im  ;
    m_sm = fmobile_m*c_m*v_m*eps_m*s_m + m_sm_im ;
    /* coefficient de transfert */
    k_f  = rho*v_f*K_f/mu*COURBE(p_g - p_f,cb[1]) ;
    k_m  = rho*v_m*K_m/mu*COURBE(p_g - p_m,cb[3]) ;
    k_a  = rho*K_a/mu*0.5*(COURBE(p_g - p_f,cb[4]) + COURBE(p_g - p_m,cb[4])) ;
    /* flux de masse d'eau */
    grad(x,u,dh,gp_f,el.nn,dim,I_p_f) ;
    grad(x,u,dh,gp_m,el.nn,dim,I_p_m) ;
    grad(x,u,dh,gc_f,el.nn,dim,I_c_f) ;
    grad(x,u,dh,gc_m,el.nn,dim,I_c_m) ;

    for(i=0;i<3;i++) w_f[i] = - k_f*gp_f[i] ;
    w_f[dim-1] += k_f*rho*gravite ;

    for(i=0;i<3;i++) w_m[i] = - k_m*gp_m[i] ;
    w_m[dim-1] += k_m*rho*gravite ;
    /* coefficient de diffusion et dispersion */
    k_sf = c_f/rho*k_f ;
    k_sm = c_m/rho*k_m ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sf[i][j] = 0. ;
    d_sf[0][0] = v_f*eps_f*s_f*tortuosite_f(eps_f,s_f)*d0_sf ;
    d_sf[1][1] = d_sf[0][0] ;
    d_sf[2][2] = d_sf[0][0] ;
    if((ww_f = normv(w_f)) > 0.) {
      ww_f = lambda_sf/(rho*ww_f) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	d_sf[i][j] += ww_f*w_f[i]*w_f[j] ;
      }
    }

    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sm[i][j] = 0. ;
    d_sm[0][0] = v_m*eps_m*s_m*tortuosite_m(eps_m,s_m)*d0_sm ;
    d_sm[1][1] = d_sm[0][0] ;
    d_sm[2][2] = d_sm[0][0] ;
    if((ww_m = normv(w_m)) > 0.) {
      ww_m = lambda_sm/(rho*ww_m) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	d_sm[i][j] += ww_m*w_m[i]*w_m[j] ;
      }
    }

    /* k_sa = 0.5*(c_f + c_m)*k_a/rho ; */
    k_sa = (p_f > p_m) ? c_f*k_a/rho : c_m*k_a/rho ;
    d_sa = alpha_s ;
    /* flux de masse du solute */

    for(i=0;i<3;i++) w_sf[i] = - k_sf*gp_f[i] ;
    w_sf[dim-1] += k_sf*rho*gravite ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) w_sf[i] += - d_sf[i][j]*gc_f[j] ;

    for(i=0;i<3;i++) w_sm[i] = - k_sm*gp_m[i] ;
    w_sm[dim-1] += k_sm*rho*gravite ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) w_sm[i] += - d_sm[i][j]*gc_m[j] ;
    /* termes d'echange */
    dm_a  = k_a*(p_f - p_m) ;
    dm_sa = k_sa*(p_f - p_m) + d_sa*(c_f - c_m) ;
    /* rangement dans f */
    f[p*NVI+E_fra] = m_f ;
    f[p*NVI+E_mat] = m_m ;
    f[p*NVI+E_sf]  = m_sf ;
    f[p*NVI+E_sm]  = m_sm ;
    for(i=0;i<dim;i++) {
      f[p*NVI+NEQ+E_fra*dim+i]  = w_f[i] ;
      f[p*NVI+NEQ+E_mat*dim+i]  = w_m[i] ;
      f[p*NVI+NEQ+E_sf*dim+i]   = w_sf[i] ;
      f[p*NVI+NEQ+E_sm*dim+i]   = w_sm[i] ;
    }
    f[p*NVI+NEQ+NEQ*dim+0] = dm_a ;
    f[p*NVI+NEQ+NEQ*dim+1] = dm_sa ;
    f[p*NVI+NEQ+NEQ*dim+2] = m_sf_im ;
    f[p*NVI+NEQ+NEQ*dim+3] = m_sm_im ;
    /* rangement dans va */
    va[p*NVE]   = k_f ;
    va[p*NVE+1] = k_m ;
    va[p*NVE+2] = k_a ;
    va[p*NVE+3] = k_sf ;
    va[p*NVE+4] = k_sm ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      va[p*NVE+5+i*3+j]  = d_sf[i][j] ;
      va[p*NVE+14+i*3+j] = d_sm[i][j] ;
    }

    va[p*NVE+23] = k_sa ;
    va[p*NVE+24] = d_sa ;
  }
}

int ex32(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double p_f,p_m,s_f,s_m,c_f,c_m ;
  double k_f,k_m,k_a,k_sf,k_sm,d_sf[3][3],d_sm[3][3],k_sa,d_sa ;
  double w_f[3],w_m[3] ;
  double ww_f,ww_m ;
  double eps_f,eps_m ;
  int    i,j,p ;
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
  d0_sf   = el.mat->pr[pm("D_sf")] ;
  d0_sm   = el.mat->pr[pm("D_sm")] ;
  lambda_sf = el.mat->pr[pm("Lambda_sf")] ;
  lambda_sm = el.mat->pr[pm("Lambda_sm")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;

  v_m = 1. - v_f ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pressions */
    p_f  = param(u,h,el.nn,I_p_f) ;
    p_m  = param(u,h,el.nn,I_p_m) ;
    /* concentrations */
    c_f  = param(u,h,el.nn,I_c_f) ;
    c_m  = param(u,h,el.nn,I_c_m) ;
    /* saturations */
    s_f = COURBE(p_g - p_f,cb[0]) ;
    s_m = COURBE(p_g - p_m,cb[2]) ;
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* flux */
    for(i=0;i<3;i++) w_f[i] = (i < dim) ? f[p*NVI+NEQ+E_fra*dim+i] : 0. ;
    for(i=0;i<3;i++) w_m[i] = (i < dim) ? f[p*NVI+NEQ+E_mat*dim+i] : 0. ;
    /* coefficient de transfert */
    k_f = rho*v_f*K_f/mu*COURBE(p_g - p_f,cb[1]) ;
    k_m = rho*v_m*K_m/mu*COURBE(p_g - p_m,cb[3]) ;
    k_a = rho*K_a/mu*0.5*(COURBE(p_g - p_f,cb[4]) + COURBE(p_g - p_m,cb[4])) ;
    k_sf = c_f/rho*k_f ;
    k_sm = c_m/rho*k_m ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sf[i][j] = 0. ;
    d_sf[0][0] = v_f*eps_f*s_f*tortuosite_f(eps_f,s_f)*d0_sf ;
    d_sf[1][1] = d_sf[0][0] ;
    d_sf[2][2] = d_sf[0][0] ;
    if((ww_f = normv(w_f)) > 0.) {
      ww_f = lambda_sf/(rho*ww_f) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	d_sf[i][j] += ww_f*w_f[i]*w_f[j] ;
      }
    }

    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sm[i][j] = 0. ;
    d_sm[0][0] = v_m*eps_m*s_m*tortuosite_m(eps_m,s_m)*d0_sm ;
    d_sm[1][1] = d_sm[0][0] ;
    d_sm[2][2] = d_sm[0][0] ;
    if((ww_m = normv(w_m)) > 0.) {
      ww_m = lambda_sm/(rho*ww_m) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	d_sm[i][j] += ww_m*w_m[i]*w_m[j] ;
      }
    }

    /* k_sa = 0.5*(c_f + c_m)*k_a/rho ; */
    k_sa = (p_f > p_m) ? c_f*k_a/rho : c_m*k_a/rho ;
    d_sa = alpha_s ;
    /* rangement dans va */
    va[p*NVE]   = k_f ;
    va[p*NVE+1] = k_m ;
    va[p*NVE+2] = k_a ;
    va[p*NVE+3] = k_sf ;
    va[p*NVE+4] = k_sm ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      va[p*NVE+5+i*3+j]  = d_sf[i][j] ;
      va[p*NVE+14+i*3+j] = d_sm[i][j] ;
    }

    va[p*NVE+23] = k_sa ;
    va[p*NVE+24] = d_sa ;
  }
  return(0) ;
}

int ct32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double p_f,p_m,s_f,s_m,c_f,c_m,c_f_im,c_m_im ;
  double m_f,m_m,dm_a,m_sf,m_sm,dm_sa,m_sf_im,m_sm_im ;
  double gp_f[3],gp_m[3],w_f[3],w_m[3],gc_f[3],gc_m[3],w_sf[3],w_sm[3] ;
  double k_f,k_m,k_a,k_sf,k_sm,k_sa,d_sa ;
  double d_sf[3][3],d_sm[3][3] ;
  double eps_f,eps_m ;
  int    i,j,p ;
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
  d0_sf   = el.mat->pr[pm("D_sf")] ;
  d0_sm   = el.mat->pr[pm("D_sm")] ;
  lambda_sf = el.mat->pr[pm("Lambda_sf")] ;
  lambda_sm = el.mat->pr[pm("Lambda_sm")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;
  fmobile_f = el.mat->pr[pm("fmobile_f")] ;
  fmobile_m = el.mat->pr[pm("fmobile_m")] ;
  beta_sf = el.mat->pr[pm("beta_sf")] ;
  beta_sm = el.mat->pr[pm("beta_sm")] ;

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
    /* concentrations */
    c_f  = param(u_1,h,el.nn,I_c_f) ;
    c_m  = param(u_1,h,el.nn,I_c_m) ;
    {
      double m_sf_im_n = f_n[p*NVI+NEQ+NEQ*dim+2] ;
      double m_sm_im_n = f_n[p*NVI+NEQ+NEQ*dim+3] ;
      c_f_im = (m_sf_im_n + beta_sf*dt*c_f)/((1. - fmobile_f)*v_f*eps_f*s_f + beta_sf*dt) ;
      c_m_im = (m_sm_im_n + beta_sm*dt*c_m)/((1. - fmobile_m)*v_m*eps_m*s_m + beta_sm*dt) ;
    }
    /* masses */
    m_f  = rho*v_f*eps_f*s_f ;
    m_m  = rho*v_m*eps_m*s_m ;
    m_sf_im = (1. - fmobile_f)*c_f_im*v_f*eps_f*s_f ;
    m_sm_im = (1. - fmobile_m)*c_m_im*v_m*eps_m*s_m ;
    m_sf = c_f*fmobile_f*v_f*eps_f*s_f + m_sf_im ;
    m_sm = c_m*fmobile_m*v_m*eps_m*s_m + m_sm_im ;
    /* coefficient de transfert */
    k_f  = va[p*NVE] ;
    k_m  = va[p*NVE+1] ;
    k_a  = va[p*NVE+2] ;
    k_sf = va[p*NVE+3] ;
    k_sm = va[p*NVE+4] ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      d_sf[i][j] = va[p*NVE+5+i*3+j] ;
      d_sm[i][j] = va[p*NVE+14+i*3+j] ;
    }

    k_sa = va[p*NVE+23] ;
    d_sa = va[p*NVE+24] ;
    /* flux de masse */
    grad(x,u_1,dh,gp_f,el.nn,dim,I_p_f) ;
    grad(x,u_1,dh,gp_m,el.nn,dim,I_p_m) ;
    grad(x,u_1,dh,gc_f,el.nn,dim,I_c_f) ;
    grad(x,u_1,dh,gc_m,el.nn,dim,I_c_m) ;

    for(i=0;i<3;i++) w_f[i] = - k_f*gp_f[i] ;
    w_f[dim-1] += k_f*rho*gravite ;

    for(i=0;i<3;i++) w_m[i] = - k_m*gp_m[i] ;
    w_m[dim-1] += k_m*rho*gravite ;

    for(i=0;i<3;i++) w_sf[i] = - k_sf*gp_f[i] ;
    w_sf[dim-1] += k_sf*rho*gravite ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) w_sf[i] += - d_sf[i][j]*gc_f[j] ;

    for(i=0;i<3;i++) w_sm[i] = - k_sm*gp_m[i] ;
    w_sm[dim-1] += k_sm*rho*gravite ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) w_sm[i] += - d_sm[i][j]*gc_m[j] ;
    /* terme d'echange */
    dm_a  = k_a*(p_f - p_m) ;
    dm_sa = k_sa*(p_f - p_m) + d_sa*(c_f - c_m) ;
    /* rangement dans f */
    f_1[p*NVI+E_fra] = m_f ;
    f_1[p*NVI+E_mat] = m_m ;
    f_1[p*NVI+E_sf]  = m_sf ;
    f_1[p*NVI+E_sm]  = m_sm ;
    for(i=0;i<dim;i++) {
      f_1[p*NVI+NEQ+E_fra*dim+i]  = w_f[i] ;
      f_1[p*NVI+NEQ+E_mat*dim+i]  = w_m[i] ;
      f_1[p*NVI+NEQ+E_sf*dim+i]   = w_sf[i] ;
      f_1[p*NVI+NEQ+E_sm*dim+i]   = w_sm[i] ;
    }
    f_1[p*NVI+NEQ+NEQ*dim+0] = dm_a ;
    f_1[p*NVI+NEQ+NEQ*dim+1] = dm_sa ;
    f_1[p*NVI+NEQ+NEQ*dim+2] = m_sf_im ;
    f_1[p*NVI+NEQ+NEQ*dim+3] = m_sm_im ;
  }

  return(0) ;
}

int  mx32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
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
  dec = c32(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  mxcmss(k,x,*el.fi,c,dim,dec,geom,NEQ) ;
  /*
    conduction
  */
  dec = k32(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  mxccnd(kb,x,*el.fi,c,dim,dec,geom,NEQ) ;
  for(i=0;i<NEQ*NEQ*el.nn*el.nn;i++) k[i] += dt*kb[i] ;

  return(0) ;
}

void rs32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;
  double rb[MAX_NOEUDS],g1[MAX_PGAUSS],*f_2 ;

  /* initialisation */
  for(i=0;i<NEQ*el.nn;i++) r[i] = 0. ;

  if(el.dim < dim) return ;
  
  /* fracture */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_fra] - f_n[i*NVI+E_fra] + dt*f_1[i*NVI+NEQ+NEQ*dim+0] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fra] = - rb[i] ;
  f_2 = f_1 + NEQ + E_fra*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fra] +=  dt*rb[i] ;

  /* matrice */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_mat] - f_n[i*NVI+E_mat] - dt*f_1[i*NVI+NEQ+NEQ*dim+0] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_mat] = - rb[i] ;
  f_2 = f_1 + NEQ + E_mat*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_mat] +=  dt*rb[i] ;
  
  /* solute fracture */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_sf] - f_n[i*NVI+E_sf] + dt*f_1[i*NVI+NEQ+NEQ*dim+1] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_sf] = - rb[i] ;
  f_2 = f_1 + NEQ + E_sf*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_sf] +=  dt*rb[i] ;
  
  /* solute matrice */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_sm] - f_n[i*NVI+E_sm] - dt*f_1[i*NVI+NEQ+NEQ*dim+1] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_sm] = - rb[i] ;
  f_2 = f_1 + NEQ + E_sm*dim ;
  rsflux(rb,x,*el.fi,f_2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_sm] +=  dt*rb[i] ;
}

int so32(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double p_f,p_m,s_f,s_m,c_f,c_m ;
  double w_f[3],w_m[3],w_sf[3],w_sm[3] ;
  int    i,j,p,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  crbe_t *cb = el.mat->cb ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  
  /* initialisation */
  nso = 10 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pressions */
  p_f  = param(u,h_s,el.nn,I_p_f) ;
  p_m  = param(u,h_s,el.nn,I_p_m) ;
  /* concentrations */
  c_f  = param(u,h_s,el.nn,I_c_f) ;
  c_m  = param(u,h_s,el.nn,I_c_m) ;
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
  strcpy(r[4].text,"concentration-fracture") ; r[4].n = 1 ;
  r[4].v[0] = c_f ;
  strcpy(r[5].text,"concentration-matrice") ; r[5].n = 1 ;
  r[5].v[0] = c_m ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* flux de masse */
    for(i=0;i<dim;i++) w_f[i]  = f[p*NVI+NEQ+E_fra*dim+i] ;
    for(i=0;i<dim;i++) w_m[i]  = f[p*NVI+NEQ+E_mat*dim+i] ;
    for(i=0;i<dim;i++) w_sf[i] = f[p*NVI+NEQ+E_sf*dim+i] ;
    for(i=0;i<dim;i++) w_sm[i] = f[p*NVI+NEQ+E_sm*dim+i] ;
    /* quantites exploitees par element */
    strcpy(r[6].text,"flux_fracture") ; r[6].n = 3 ;
    for(i=0;i<dim;i++) r[6].v[i] += w_f[i]/el.fi->np ;
    strcpy(r[7].text,"flux_matrice") ; r[7].n = 3 ;
    for(i=0;i<dim;i++) r[7].v[i] += w_m[i]/el.fi->np ;
    strcpy(r[8].text,"flux_solute-fracture") ; r[8].n = 3 ;
    for(i=0;i<dim;i++) r[8].v[i] += w_sf[i]/el.fi->np ;
    strcpy(r[9].text,"flux_solute-matrice") ; r[9].n = 3 ;
    for(i=0;i<dim;i++) r[9].v[i] += w_sm[i]/el.fi->np ;
  }
  return (nso) ;
}


int c32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  double p_f,p_m,c_f,c_m,c_f_im,c_m_im ;
  double s_f,s_m,ds_fsdpc,ds_msdpc ;
  double dc_f_imsdp_f,dc_f_imsdc_f,dc_m_imsdp_m,dc_m_imsdc_m ;
  double k_a,k_sa,d_sa ;
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
  d0_sf   = el.mat->pr[pm("D_sf")] ;
  d0_sm   = el.mat->pr[pm("D_sm")] ;
  lambda_sf = el.mat->pr[pm("Lambda_sf")] ;
  lambda_sm = el.mat->pr[pm("Lambda_sm")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;
  fmobile_f = el.mat->pr[pm("fmobile_f")] ;
  fmobile_m = el.mat->pr[pm("fmobile_m")] ;
  beta_sf = el.mat->pr[pm("beta_sf")] ;
  beta_sm = el.mat->pr[pm("beta_sm")] ;

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
    /* porosites */
    eps_f = eps_f0 + S_s*(p_f - p_f0) ;
    eps_m = eps_m0 + S_s*(p_m - p_m0) ;
    /* concentrations */
    c_f  = param(u_1,h,el.nn,I_c_f) ;
    c_m  = param(u_1,h,el.nn,I_c_m) ;
    {
      double m_sf_im_n = f_n[p*NVI+NEQ+NEQ*dim+2] ;
      double m_sm_im_n = f_n[p*NVI+NEQ+NEQ*dim+3] ;
      c_f_im = (m_sf_im_n + beta_sf*dt*c_f)/((1. - fmobile_f)*v_f*eps_f*s_f + beta_sf*dt) ;
      c_m_im = (m_sm_im_n + beta_sm*dt*c_m)/((1. - fmobile_m)*v_m*eps_m*s_m + beta_sm*dt) ;
    }
    /* derivees des saturations */
    ds_fsdpc = DCOURBE(p_g - p_f,cb[0]) ;
    ds_msdpc = DCOURBE(p_g - p_m,cb[2]) ;
    /* derivees des concentrations de solutes type im */
    dc_f_imsdp_f = - c_f_im*(1. - fmobile_f)*v_f*(-eps_f*ds_fsdpc + s_f*S_s)/((1. - fmobile_f)*v_f*eps_f*s_f + beta_sf*dt) ;
    dc_f_imsdc_f = beta_sf*dt/((1. - fmobile_f)*v_f*eps_f*s_f + beta_sf*dt) ;
    dc_m_imsdp_m = - c_m_im*(1. - fmobile_m)*v_m*(-eps_m*ds_msdpc + s_m*S_s)/((1. - fmobile_m)*v_m*eps_m*s_m + beta_sm*dt) ;
    dc_m_imsdc_m = beta_sm*dt/((1. - fmobile_m)*v_m*eps_m*s_m + beta_sm*dt) ;
    /* coefficient de transfert */
    k_a  = va[p*NVE+2] ;
    k_sa = va[p*NVE+23] ;
    d_sa = va[p*NVE+24] ;
    /* fracture */
    c1[E_fra*NEQ+I_p_f] =   rho*v_f*(eps_f*(-ds_fsdpc) + S_s*s_f) + dt*k_a ;
    c1[E_fra*NEQ+I_p_m] = - dt*k_a ;
    /* matrice */
    c1[E_mat*NEQ+I_p_f] = - dt*k_a ;
    c1[E_mat*NEQ+I_p_m] =   rho*v_m*(eps_m*(-ds_msdpc) + S_s*s_m) + dt*k_a ;
    /* solute fracture */
    c1[E_sf*NEQ+I_p_f]  =   (c_f*fmobile_f + c_f_im*(1. - fmobile_f))*v_f*(eps_f*(-ds_fsdpc) + S_s*s_f) + dc_f_imsdp_f*(1. - fmobile_f)*v_f*eps_f*s_f + dt*k_sa ;
    c1[E_sf*NEQ+I_p_m]  = - dt*k_sa ;
    c1[E_sf*NEQ+I_c_f]  =   fmobile_f*v_f*eps_f*s_f + dc_f_imsdc_f*(1. - fmobile_f)*v_f*eps_f*s_f + dt*d_sa ;
    c1[E_sf*NEQ+I_c_m]  = - dt*d_sa ;
    /* solute matrice */
    c1[E_sm*NEQ+I_p_f]  = - dt*k_sa ;
    c1[E_sm*NEQ+I_p_m]  =   (c_m*fmobile_m + c_m_im*(1. - fmobile_m))*v_m*(eps_m*(-ds_msdpc) + S_s*s_m) + dc_m_imsdp_m*(1. - fmobile_m)*v_m*eps_m*s_m + dt*k_sa ;
    c1[E_sm*NEQ+I_c_f]  = - dt*d_sa ;
    c1[E_sm*NEQ+I_c_m]  =   fmobile_m*v_m*eps_m*s_m + dc_m_imsdc_m*(1. - fmobile_m)*v_m*eps_m*s_m + dt*d_sa ;
  }
  return(dec) ;
}


int k32(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  double k_f,k_m,k_a,k_sf,k_sm,k_sa,d_sa ;
  double d_sf[3][3],d_sm[3][3] ;
  int    dec ;
  double *c1 ;
  int    i,j,p ;
  double zero = 0. ;
  
  dec = 9*NEQ*NEQ ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* coefficients de transfert */
    k_f  = va[p*NVE] ;
    k_m  = va[p*NVE+1] ;
    k_a  = va[p*NVE+2] ;
    k_sf = va[p*NVE+3] ;
    k_sm = va[p*NVE+4] ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sf[i][j] = va[p*NVE+5+i*3+j] ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) d_sm[i][j] = va[p*NVE+14+i*3+j] ;

    k_sa = va[p*NVE+23] ;
    d_sa = va[p*NVE+24] ;
    /* fracture */
    c1[(E_fra*NEQ+I_p_f)*9+0] = k_f ;
    c1[(E_fra*NEQ+I_p_f)*9+4] = k_f ;
    c1[(E_fra*NEQ+I_p_f)*9+8] = k_f ;
    /* matrice */
    c1[(E_mat*NEQ+I_p_m)*9+0] = k_m ;
    c1[(E_mat*NEQ+I_p_m)*9+4] = k_m ;
    c1[(E_mat*NEQ+I_p_m)*9+8] = k_m ;
    /* solute fracture */
    c1[(E_sf*NEQ+I_p_f)*9+0]  = k_sf ;
    c1[(E_sf*NEQ+I_p_f)*9+4]  = k_sf ;
    c1[(E_sf*NEQ+I_p_f)*9+8]  = k_sf ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) c1[(E_sf*NEQ+I_c_f)*9+i*3+j] = d_sf[i][j] ;
    /* solute matrice */
    c1[(E_sm*NEQ+I_p_m)*9+0]  = k_sm ;
    c1[(E_sm*NEQ+I_p_m)*9+4]  = k_sm ;
    c1[(E_sm*NEQ+I_p_m)*9+8]  = k_sm ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) c1[(E_sm*NEQ+I_c_m)*9+i*3+j] = d_sm[i][j] ;
  }
  return(dec) ;
}

double tortuosite_f(double phi,double s)
{
  return(1.) ;
}

double tortuosite_m(double phi,double s)
{
  return(1.) ;
}

double normv(double *v)
{
  double vv = 0. ;
  int i ;

  for(i=0;i<3;i++) vv += v[i]*v[i] ;

  return(sqrt(vv)) ;
}

#undef NEQ
#undef NVI
#undef NVE

#undef COURBE
#undef DCOURBE
