#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  16
#define TITLE "Poroplasticite (Drucker-Prager)"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ   (1+dim)
#define NVE   (26)
#define NVA   (1)
#define E_liq (dim)
#define E_mec (0)
#define I_p_l (dim)
#define I_u   (0)
/* Fonctions */
static int    pm(const char *s) ;
static int    c16(double **,double **,double **,double *,double *,double *,elem_t,int,geom_t,double,double *) ;
static int    k16(double **,double **,double **,double *,double *,double *,elem_t,int,geom_t,double,double *) ;
static double rn16(double *,double *,mate_t) ;
static double cr16(double *,double *,double *,mate_t) ;
/* Parametres */
static double gravite,young,poisson,phi0,k_int,mu_l,rho_l0,p_l0,b,N,rho_s,cohe,af,ad,beta,k_l,sig0_11,sig0_22,sig0_33,a_int ;

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"young") == 0) return (1) ;
  else if(strcmp(s,"poisson") == 0) return (2) ;
  else if(strcmp(s,"porosite") == 0) return (3) ;
  else if(strcmp(s,"rho_l") == 0) return (4) ;
  else if(strcmp(s,"k_int") == 0) return (5) ;
  else if(strcmp(s,"mu_l") == 0) return (6) ;
  else if(strcmp(s,"b") == 0) return (7) ;
  else if(strcmp(s,"N") == 0) return (8) ;
  else if(strcmp(s,"rho_s") == 0) return (9) ;
  else if(strcmp(s,"cohesion") == 0) return (10) ;
  else if(strcmp(s,"frottement") == 0) return (11) ;
  else if(strcmp(s,"dilatance") == 0) return (12) ;
  else if(strcmp(s,"beta") == 0) return (13) ;
  else if(strcmp(s,"p_l0") == 0) return (14) ;
  else if(strcmp(s,"k_l") == 0) return (15) ;
  else if(strcmp(s,"sig0_11") == 0) return (16) ;
  else if(strcmp(s,"sig0_22") == 0) return (17) ;
  else if(strcmp(s,"sig0_33") == 0) return (18) ;
  else if(strcmp(s,"a_int") == 0) return (19) ;
  else { printf("donnee \"%s\" non connue (pm16)\n",s) ; exit(0) ; }
}

int dm16(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int i,n_donnees = 20 ;
  
  mat->neq = NEQ ;
  for(i=0;i<dim;i++) sprintf(mat->eqn[E_mec+i],"meca_%d",i+1) ;
  for(i=0;i<dim;i++) sprintf(mat->inc[I_u+i],"u_%d",i+1) ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  
  dmat(mat,ficd,pm,n_donnees) ;

  /*
    Donnees
  */
  gravite = mat->pr[pm("gravite")] ;
  young   = mat->pr[pm("young")] ;
  poisson = mat->pr[pm("poisson")] ;
  phi0    = mat->pr[pm("porosite")] ;
  k_int   = mat->pr[pm("k_int")] ;
  a_int   = mat->pr[pm("a_int")] ;
  mu_l    = mat->pr[pm("mu_l")] ;
  rho_l0  = mat->pr[pm("rho_l")] ;
  k_l     = mat->pr[pm("k_l")] ;
  rho_s   = mat->pr[pm("rho_s")] ;
  p_l0    = mat->pr[pm("p_l0")] ;
  b       = mat->pr[pm("b")] ;
  N       = mat->pr[pm("N")] ;
  beta    = mat->pr[pm("beta")] ;
  sig0_11 = mat->pr[pm("sig0_11")] ;
  sig0_22 = mat->pr[pm("sig0_22")] ;
  sig0_33 = mat->pr[pm("sig0_33")] ;
  cohe    = mat->pr[pm("cohesion")] ;
  af      = mat->pr[pm("frottement")]*M_PI/180. ;
  ad      = mat->pr[pm("dilatance")]*M_PI/180. ;

  return(mat->n) ;
}

int qm16(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est formee de (1 + dim) equations :\n\
\t 1. l\'equation de conservation de la masse d\'eau (p_l)\n\
\t 2. les equations d\'equilibre mecanique (u_1,u_2,u_3)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;
  
  fprintf(ficd,"gravite = 0       # gravite\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 5.8e+09   # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.3     # coefficient de Poisson\n") ;
  fprintf(ficd,"porosite = 0.15   # porosite\n") ;
  fprintf(ficd,"rho_l = 1000      # masse volumique du fluide\n") ;
  fprintf(ficd,"p_l0 = 4.7e+06    # pression initiale du fluide\n") ;
  fprintf(ficd,"p_g = 0           # pression du gaz\n") ;
  fprintf(ficd,"k_l = 2e+09       # module de compression du fluide\n") ;
  fprintf(ficd,"k_int = 1e-19     # permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001      # viscosite du liquide\n") ;
  fprintf(ficd,"b = 0.8           # coefficient de Biot\n") ;
  fprintf(ficd,"N = 4.e-11        # compressibilite des pores\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion\n") ;
  fprintf(ficd,"frottement = 25   # frottement\n") ;
  fprintf(ficd,"dilatance = 25    # dilatance \n") ;
  fprintf(ficd,"beta = 0.8        # coefficient beta\n") ;
  fprintf(ficd,"sig0_11 = -11.5e6 # contrainte initiale sig0_11\n") ;
  fprintf(ficd,"sig0_22 = -11.5e6 # contrainte initiale sig0_22\n") ;
  fprintf(ficd,"sig0_33 = -11.5e6 # contrainte initiale sig0_33\n") ;
  fprintf(ficd,"a_int = 2.e10     # coefficient a_int\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;
  
  return(NEQ) ;
}

void tb16(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVE*el->fi->np ;
  el->n_ve = NVA*el->fi->np ;
}

void ch16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  /* hydraulique */
  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_liq) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  } else if(!strcmp(cg.eqn,el.mat->eqn[E_liq])) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }
}

void in16(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double dmu,lame,rho_l,k_h ;
  double eps[9],sig[9],phi,dphi,pl,dpl,tre,m_l,w_l[3],gpl[3],eps_p[9],x_p[3] ; 
  double dfsds[9],dgsds[9],crit ;
  int    i,p ;
  double *h,*dh ;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  a_int   = el.mat->pr[pm("a_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l0  = el.mat->pr[pm("rho_l")] ;
  k_l     = el.mat->pr[pm("k_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  beta    = el.mat->pr[pm("beta")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;
  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    tre  = eps[0] + eps[4] + eps[8] ;
    /* pressions */
    pl  = param(u,h,el.nn,I_p_l) ;
    dpl = pl - p_l0 ;
    /* porosite */
    dphi = b*tre + N*dpl ;
    phi  = phi0 + dphi ;
    /* coordonnees */
    for(i=0;i<3;i++) x_p[i] = param(x,h,el.nn,i) ;
    /* contraintes initiales */
    /*
    {
    Field_t *ch = Material_GetField(Element_GetMaterial(&el)) ;
    i = floor(el.mat->pr[pm("sig0_11")]+0.5) ;
    sig0_11 = (i < 1) ? zero : champ(x_p,dim,ch[i-1]) ;
    i = floor(el.mat->pr[pm("sig0_22")]+0.5) ;
    sig0_22 = (i < 1) ? zero : champ(x_p,dim,ch[i-1]) ;
    i = floor(el.mat->pr[pm("sig0_33")]+0.5) ;
    sig0_33 = (i < 1) ? zero : champ(x_p,dim,ch[i-1]) ;
    }
    */
    /* contraintes effectives */
    for(i=0;i<9;i++) sig[i] = dmu*eps[i] ;
    sig[0] += sig0_11 + lame*tre - b*dpl + beta*pl ;
    sig[4] += sig0_22 + lame*tre - b*dpl + beta*pl ;
    sig[8] += sig0_33 + lame*tre - b*dpl + beta*pl ;
    /* deformations plastiques */
    for(i=0;i<9;i++) eps_p[i] = 0. ;
    /* critere */
    crit = cr16(sig,dfsds,dgsds,*el.mat) ;
    /* contraintes totales */
    sig[0] += - beta*pl ;
    sig[4] += - beta*pl ;
    sig[8] += - beta*pl ; 
    /* masse volumique du fluide */
    rho_l = rho_l0*(1. + dpl/k_l) ;
    /* masse liquide */
    m_l = rho_l*phi ;
    /* coefficient de transfert */
    k_h = rho_l*k_int/mu_l ;
    if(dphi > 0) {
      if(dphi < 1.e-2) {
	k_h *= 1. + a_int*dphi*dphi*dphi ;
      } else k_h *= 1. + a_int*1.e-6 ;
    }
    /* flux */
    grad(x,u,dh,gpl,el.nn,dim,I_p_l) ;
    for(i=0;i<3;i++) w_l[i] = - k_h*gpl[i] ;
    w_l[dim-1] += k_h*rho_l*gravite ;
    /* rangement dans f et va */
    va[p*NVA] = k_h ;
    f[p*NVE] = m_l ;
    for(i=0;i<3;i++) f[p*NVE+1+i]  = w_l[i] ;
    for(i=0;i<9;i++) f[p*NVE+4+i]  = sig[i] ;
    for(i=0;i<3;i++) f[p*NVE+13+i] = zero ;
    f[p*NVE+13+dim-1] = (rho_s+m_l)*gravite ;
    for(i=0;i<9;i++) f[p*NVE+16+i]  = eps_p[i] ;
    f[p*NVE+25] = crit ;
  }
}

int ex16(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double rho_l,k_h ;
  double eps[9],eps_p[9],dphi,pl,tre,phi_p,tre_p ; 
  int    i,p ;
  double *h,*dh ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  k_int   = el.mat->pr[pm("k_int")] ;
  a_int   = el.mat->pr[pm("a_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l0  = el.mat->pr[pm("rho_l")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  k_l     = el.mat->pr[pm("k_l")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  beta    = el.mat->pr[pm("beta")] ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    tre  = eps[0] + eps[4] + eps[8] ;
    /* pressions */
    pl  = param(u,h,el.nn,I_p_l) ;
    /* deformations plastiques */
    for(i=0;i<9;i++) eps_p[i] = f[p*NVE+16+i] ;
    /* porosite */
    tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
    phi_p = beta*tre_p ;
    dphi  = b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
    /* masse volumique du fluide */
    rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;
    /* rangement dans va */
    k_h = rho_l*k_int/mu_l ;
    if(dphi > 0) {
      if(dphi < 1.e-2) {
	k_h *= 1. + a_int*dphi*dphi*dphi ;
      } else k_h *= 1. + a_int*1.e-6 ;
    }
    va[p*NVA] = k_h ;
  }
  return(0) ;
}

int ct16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double dmu,lame,mu,rho_l,k_h ;
  double eps[9],sig[9],phi,pl,tre,m_l,w_l[3],gpl[3],eps_p[9],phi_p,tre_p ;
  double eps_n[9],pl_n ;
  double deps[9],trde,dpl ;
  double crit ;
  int    i,p ;
  double *h,*dh ;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l0  = el.mat->pr[pm("rho_l")] ;
  k_l     = el.mat->pr[pm("k_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  beta    = el.mat->pr[pm("beta")] ;
  dmu     = young/(un+poisson) ;
  mu      = dmu/2. ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u_1,h,dh,eps,el.nn,dim,geom,I_u) ;
    def(x,u_n,h,dh,eps_n,el.nn,dim,geom,I_u) ;
    for(i=0;i<9;i++) deps[i] =  eps[i] - eps_n[i] ;
    tre  = eps[0] + eps[4] + eps[8] ;
    trde = deps[0] + deps[4] + deps[8] ;
    /* pressions */
    pl   = param(u_1,h,el.nn,I_p_l) ;
    pl_n = param(u_n,h,el.nn,I_p_l) ;
    dpl  = pl - pl_n ;
    /* contraintes effectives */
    for(i=0;i<9;i++) sig[i] = f_n[p*NVE+4+i] + dmu*deps[i] ;
    sig[0] += lame*trde - b*dpl + beta*pl ;
    sig[4] += lame*trde - b*dpl + beta*pl ;
    sig[8] += lame*trde - b*dpl + beta*pl ;
    /* deformations plastiques */
    for(i=0;i<9;i++) eps_p[i] = f_n[p*NVE+16+i] ;
    /* projection */
    crit = rn16(sig,eps_p,*el.mat) ;
    /* contraintes totales */
    sig[0] -= beta*pl ;
    sig[4] -= beta*pl ;
    sig[8] -= beta*pl ;
    /* porosite */
    tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
    phi_p = beta*tre_p ;
    phi   = phi0 + b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
    /* masse volumique du fluide */
    rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;
    /* masse liquide */
    m_l = rho_l*phi ;
    /* coefficient de transfert */
    k_h = va[p*NVA] ;
    /* flux */
    grad(x,u_1,dh,gpl,el.nn,dim,I_p_l) ;
    for(i=0;i<3;i++) w_l[i] = - k_h*gpl[i] ;
    w_l[dim-1] += k_h*rho_l*gravite ;
    /* rangement dans f_1 */
    f_1[p*NVE] = m_l ;
    for(i=0;i<3;i++) f_1[p*NVE+1+i] = w_l[i] ;
    for(i=0;i<9;i++) f_1[p*NVE+4+i] = sig[i] ;
    for(i=0;i<3;i++) f_1[p*NVE+13+i] = zero ;
    f_1[p*NVE+13+dim-1] = (rho_s+m_l)*gravite ;
    for(i=0;i<9;i++) f_1[p*NVE+16+i]  = eps_p[i] ; 
    f_1[p*NVE+25] = crit ;
  }
  return(0) ;
}

int mx16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  int    i,j,n,m,dec ;
  double c[MAX_PGAUSS*100] ;
  double kb[9*MAX_NOEUDS*MAX_NOEUDS] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
  ** Matrice de comportement
  */
  dec = c16(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  /* 
  ** 1.  Mecanique
  */
  /* 1.1 rigidite du squelette */
  mxrig(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) for(j=0;j<dim;j++) {
    K(E_mec+i+n*NEQ,I_u+j+m*NEQ) = kb[(i+n*dim)*el.nn*dim+j+m*dim] ;
  }
  /* 1.2 couplage mecanique */
  mxcpl(kb,x,*el.fi,c+81,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) {
    K(E_mec+i+n*NEQ,I_p_l+m*NEQ) = kb[(i+n*dim)*el.nn+m] ;
  }
  /*
  ** 2.  Hydraulique
  */
  /* 2.1 couplage diffusif */
  mxcpl(kb,x,*el.fi,c+90,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) {
    K(E_liq+m*NEQ,I_u+i+n*NEQ) = kb[(i+n*dim)*el.nn+m] ;
  }
  /* 2.2 accumulation */
  mxmass(kb,x,*el.fi,c+99,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_liq+n*NEQ,I_p_l+m*NEQ) = kb[n*el.nn+m] ;
  }

  /*
  ** Matrice de conduction
  */
  dec = k16(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  mxcond(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_liq+n*NEQ,I_p_l+m*NEQ) += dt*kb[n*el.nn+m] ;
  }

  return(0) ;
#undef K
}

void rs16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i,n ;
  double rb[3*MAX_NOEUDS],g1[MAX_PGAUSS] ;
  double zero = 0. ;
  
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /* 1. Mecanique */
  rscont(rb,x,*el.fi,f_1+4,dim,NVE,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) R(n,E_mec+i) -= rb[i+n*dim] ;
  rsmass(rb,x,*el.fi,f_1+13+dim-1,dim,NVE,geom) ;
  for(n=0;n<el.nn;n++) R(n,E_mec+dim-1) -= -rb[n] ;
  /* 2. Hydraulique */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVE] - f_n[i*NVE] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= rb[i] ;
  rsflux(rb,x,*el.fi,f_1+1,dim,NVE,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= -dt*rb[i] ;
#undef R
}

int so16(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,p,nso ;
  double *h,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  
  /* initialisation */
  nso = 6 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] =  param(u,h_s,el.nn,I_p_l) ;
  strcpy(r[1].text,"deplacements") ; r[1].n = 3 ;
  for(i=0;i<dim;i++) r[1].v[i] = param(u,h_s,el.nn,I_u+i) ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    h  = el.fi->h  + p*el.nn ;
    strcpy(r[2].text,"flux-de-masse") ; r[2].n = 3 ;
    for(i=0;i<dim;i++) r[2].v[i] += f[p*NVE+1+i]/el.fi->np ;
    strcpy(r[3].text,"contraintes") ; r[3].n = 9 ;
    for(i=0;i<9;i++) r[3].v[i] += f[p*NVE+4+i]/el.fi->np ;
    strcpy(r[4].text,"deformations-plastiques") ; r[4].n = 9 ;
    for(i=0;i<9;i++) r[4].v[i] += f[p*NVE+16+i]/el.fi->np ;
    strcpy(r[5].text,"k_h") ; r[5].n = 1 ;
    r[5].v[0] += va[p*NVA]/el.fi->np ;
  }
  return (nso) ;
}

int c16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  double dmu,lame,mu ;
  double pl,rho_l ;
  double sig[9],eps_p[9],phi,eps[9],tre,tre_p,phi_p ;
  double crit,dfsds[3][3],dgsds[3][3],fc[3][3],cg[3][3],fcg,trf,trg ;
  int    dec ;
  int    i,j,k,l,p ;
  double *h,*dh,*c1 ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("porosite")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l0  = el.mat->pr[pm("rho_l")] ;
  k_l     = el.mat->pr[pm("k_l")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  beta    = el.mat->pr[pm("beta")] ;
  dmu     = young/(un+poisson) ;
  mu      = dmu/deux ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  
  dec = 100 ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u_1,h,dh,eps,el.nn,dim,geom,I_u) ;
    tre  = eps[0] + eps[4] + eps[8] ;
    /* pressions */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    /* contraintes effectives */
    for(i=0;i<9;i++) sig[i] = f_1[p*NVE+4+i] ;
    sig[0] += beta*pl ;
    sig[4] += beta*pl ;
    sig[8] += beta*pl ;
    /* deformations plastiques */
    for(i=0;i<9;i++) eps_p[i] = f_1[p*NVE+16+i] ;
    /* porosite */
    tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
    phi_p = beta*tre_p ;
    phi   = phi0 + b*(tre - tre_p) + N*(pl - p_l0) + phi_p ;
    /* masse volumique du fluide */
    rho_l = rho_l0*(1. + (pl - p_l0)/k_l) ;

    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    
    /* mecanique */
    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      C1(i,i,j,j) += lame ;
      C1(i,j,i,j) += mu ;
      C1(i,j,j,i) += mu ;
    }
    /* critere */
    crit = f_1[p*NVE+25] ;
    if(crit >= 0.) {
      cr16(sig,dfsds[0],dgsds[0],*el.mat) ;
      fcg = 0. ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	fc[i][j] = 0. ;
	cg[i][j] = 0. ;
	for(k=0;k<3;k++) for(l=0;l<3;l++) {
	  fc[i][j] += dfsds[k][l]*C1(l,k,i,j) ;
	  cg[i][j] += C1(i,j,k,l)*dgsds[l][k] ;
	  fcg += dfsds[j][i]*C1(i,j,k,l)*dgsds[l][k] ;
	}
      }
      if(fcg > 0.) fcg = 1./fcg ;
      else {
	printf("\n") ;
	printf("dfsds = ") ;
	for(i=0;i<3;i++) for(j=0;j<3;j++) printf(" %e",dfsds[i][j]) ;
	printf("\n") ;
	printf("dgsds = ") ;
	for(i=0;i<3;i++) for(j=0;j<3;j++) printf(" %e",dgsds[i][j]) ;
	printf("\n") ;
	printf("fcg = %e\n",fcg) ;
	printf("x =") ;
	for(i=0;i<dim;i++) printf(" %e",param(x,h,el.nn,i)) ;
	printf("\n") ;
	exit(0) ;
      }
      for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) for(l=0;l<3;l++) {
	C1(i,j,k,l) -= cg[i][j]*fc[k][l]*fcg ;
      }
    }
    c1 += 81 ;
    for(i=0;i<3;i++) B1(i,i) = -b ;
    if(crit >= 0.) {
      trf = dfsds[0][0] + dfsds[1][1] + dfsds[2][2] ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	B1(i,j) -= cg[i][j]*(beta-b)*trf*fcg ;
      }
    }
    
    c1 += 9 ;
    /* hydraulique */
    for(i=0;i<3;i++) B1(i,i) = rho_l*b ;
    if(crit >= 0.) {
      trg = dgsds[0][0] + dgsds[1][1] + dgsds[2][2] ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	B1(i,j) += rho_l*fc[i][j]*(beta-b)*trg*fcg ;
      }
    }
    c1 += 9 ;
    c1[0] = rho_l*N + rho_l0*phi/k_l ;
    if(crit >= 0.) {
      c1[0] += rho_l*(beta-b)*(beta-b)*trf*trg*fcg ;
    }
  }
  return(dec) ;
  
#undef C1
#undef B1
}

int k16(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
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

double cr16(double *sig,double *dfsds,double *dgsds,mate_t mat)
/* Critere de Drucker-Prager */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9] ;
  double crit,lambda ;
  double ff,dd,cc,k,dmu,mu ;
  int    i ;
  /*
    Donnees
  */
  young   = mat.pr[pm("young")] ;
  poisson = mat.pr[pm("poisson")] ;
  cohe    = mat.pr[pm("cohesion")] ;
  af      = mat.pr[pm("frottement")]*M_PI/180. ;
  ad      = mat.pr[pm("dilatance")]*M_PI/180. ;
  dmu     = young/(1.+poisson) ;
  mu      = dmu/2. ;
  k       = young/(1. - 2.*poisson)/3. ;
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc      = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*cohe ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc      = 6.*cos(af)/(3. - sin(af))*cohe ;
  /*
    Le critere
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q + ff*p - cc ;
  /*
    Le deviateur
  */
  for(i=0;i<9;i++) dev[i] = sig[i] - p*id[i] ;
  /*
    Le gradient du critere
  */
  if(q > 0.) {
    for(i=0;i<9;i++) dfsds[i] = 1.5*dev[i]/q + id[i]*ff/3. ;
  } else {
    for(i=0;i<9;i++) dfsds[i] = id[i]*ff/3. ;
  }
  /*
    Le gradient du potentiel
  */
  /* cas elastique */
  if(crit <= 0.) {
    if(q > 0.) {
      for(i=0;i<9;i++) dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
    } else {
      for(i=0;i<9;i++) dgsds[i] = id[i]*dd/3. ;
    }
  }
  /* cas plastique */
  if(crit > 0.) {
    /* cas regulier */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      for(i=0;i<9;i++) dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
    /* regime de coin */
    } else {
      lambda = (ff*p - cc)/(k*ff*dd) ;
      for(i=0;i<9;i++) dgsds[i] = dev[i]/(dmu*lambda) + id[i]*dd/3. ;
    }
  }
  return(crit) ;
}

double rn16(double *sig,double *eps_p,mate_t mat)
/* Critere de Drucker-Prager : return mapping */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,p,q,dev_t[9] ;
  double deps_p[9] ;
  double crit,lambda ;
  double ff,dd,cc,k,dmu,mu ;
  int    i ;
  /*
    Donnees
  */
  young   = mat.pr[pm("young")] ;
  poisson = mat.pr[pm("poisson")] ;
  cohe    = mat.pr[pm("cohesion")] ;
  af      = mat.pr[pm("frottement")]*M_PI/180. ;
  ad      = mat.pr[pm("dilatance")]*M_PI/180. ;
  dmu     = young/(1.+poisson) ;
  mu      = dmu/2. ;
  k       = young/(1. - 2.*poisson)/3. ;
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc      = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*cohe ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc      = 6.*cos(af)/(3. - sin(af))*cohe ;
  /*
    Le critere
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*j2(sig)) ;
  crit = q_t + ff*p_t - cc ;
  /*
    Le deviateur
  */
  for(i=0;i<9;i++) dev_t[i] = sig[i] - p_t*id[i] ;
  /*
    Deformations plastiques et contraintes
  */
  if(crit >= 0.) {
    lambda = crit/(3*mu+k*ff*dd) ;
    q = q_t - lambda*3*mu ;
    /* cas regulier */
    if(q > 0.) {
      p = p_t - lambda*k*dd ;
      for(i=0;i<9;i++) {
	deps_p[i] = lambda*(1.5*dev_t[i]/q_t + id[i]*dd/3.) ;
	sig[i]    = dev_t[i]*q/q_t + p*id[i] ;
      }
    /* regime de coin */
    } else {
      lambda = (ff*p_t - cc)/(k*ff*dd) ;
      p = cc/ff ;
      q = 0. ;
      for(i=0;i<9;i++) {
	deps_p[i] = dev_t[i]/dmu + lambda*id[i]*dd/3. ;
	sig[i]    = p*id[i] ;
      }
    }
    for(i=0;i<9;i++) eps_p[i] += deps_p[i] ;
  }
  return(crit) ;
}

#undef NEQ
#undef NVE
#undef NVA
#undef E_liq
#undef E_mec
#undef I_p_l
#undef I_u
