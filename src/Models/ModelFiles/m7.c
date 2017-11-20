#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

#define MODELINDEX  7
#define TITLE   "Poroelasticite non sature"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ   (1+dim)
#define NVI   (17)
#define NVE   (1)

#define E_liq (dim)
#define E_mec (0)

#define I_p_l (dim)
#define I_u   (0)

#define M_l     (f[0])
#define W_l     (f + 1)
#define SIG     (f + 4)
#define F_MASS  (f + 13)
#define PHI     (f[16])

#define M_ln    (f_n[0])

#define K_l     (va[0])

/* Fonctions */
static int    pm(const char*) ;
static double pie(double,double,crbe_t) ;
static int    c7(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static int    k7(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static void mxnd(double*,elem_t,int,int) ;
static void rsnd(double*,elem_t,int) ;
/* Parametres */
static double gravite ;
static double young,poisson,b,N ;
static double k_int,mu_l ;
static double phi0,rho_l,rho_s ;
static double p_g,p_l0,sig0_11,sig0_22,sig0_33 ;

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"young") == 0) return (1) ;
  else if(strcmp(s,"poisson") == 0) return (2) ;
  else if(strcmp(s,"phi") == 0) return (3) ;
  else if(strcmp(s,"rho_l") == 0) return (4) ;
  else if(strcmp(s,"p_l0") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_l") == 0) return (7) ;
  else if(strcmp(s,"b") == 0) return (8) ;
  else if(strcmp(s,"N") == 0) return (9) ;
  else if(strcmp(s,"p_g") == 0) return (10) ;
  else if(strcmp(s,"rho_s") == 0) return (11) ;
  else if(strcmp(s,"sig0_11") == 0) return (12) ;
  else if(strcmp(s,"sig0_22") == 0) return (13) ;
  else if(strcmp(s,"sig0_33") == 0) return (14) ;
  else if(strcmp(s,"courbes") == 0) return (15) ;
  else return(-1) ;
}

int dm7(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 16 ;
  int i ;

  mat->neq = NEQ ;

  for(i=0;i<dim;i++) {
    sprintf(mat->eqn[E_mec + i],"meca_%d",i+1) ;
    sprintf(mat->inc[I_u + i],"u_%d",i+1) ;
  }
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->inc[I_p_l],"p_l") ;

  /* Par defaut tout a 0 */
  for(i=0;i<n_donnees;i++) mat->pr[i] = 0. ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}

int qm7(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{

  printf(TITLE) ;

  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 1+dim equations:\n\
\t 1. l\'equation de conservation de la masse d\'eau (p_l)\n\
\t 2. les equations d\'equilibre mecanique (u_1,u_2,u_3)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0       # gravite\n") ;
  fprintf(ficd,"rho_s = 0         # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 0.833333  # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.25    # coefficient de Poisson\n") ;
  fprintf(ficd,"phi = 0.3         # porosite\n") ;
  fprintf(ficd,"rho_l = 1         # masse volumique du fluide\n") ;
  fprintf(ficd,"p_l0 = 0          # pression initiale du fluid\n") ;
  fprintf(ficd,"k_int = 1         # permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 1          # viscosite du liquide\n") ;
  fprintf(ficd,"p_g = 0           # pression de gaz\n") ;
  fprintf(ficd,"b = 1             # coefficient de Biot\n") ;
  fprintf(ficd,"N = 0             # compressibilite des pores\n") ;
  fprintf(ficd,"sig0_11 = 0       # contrainte initiale 11\n") ;
  fprintf(ficd,"sig0_22 = 0       # contrainte initiale 22\n") ;
  fprintf(ficd,"sig0_33 = 0       # contrainte initiale 33\n") ;
  fprintf(ficd,"courbes = My_file # Nom du fichier : p_c S_l k_rl\n") ;

  return(NEQ) ;
}

void tb7(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}

void ch7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  /* hydraulique */
  if(Element_FindEquationPositionIndex(&el,cg.eqn) == E_liq) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }
}

void in7(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */
{
  double dmu,lame,pp0 ;
  int    i,p ;
  double un = 1.,deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;

  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  pp0     = pie(p_l0,p_g,el.mat->cb[0]) ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,f += NVI,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    double *dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    double pl  = param(u,h,el.nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = pie(pl,p_g,el.mat->cb[0]) ;
    /* deformations */
    double eps[9] ;
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    {
      double tre  = eps[0] + eps[4] + eps[8] ;
      /* porosite */
      double dphi = b*tre + N*(pp - pp0) ;
      double phi  = phi0 + dphi ;
      /* saturation */
      double sl  = courbe(pc,el.mat->cb[0]) ;
      /* contraintes */
      double sig[9] ;
      for(i=0;i<9;i++) sig[i] = dmu*eps[i] ;
      sig[0] += sig0_11 + lame*tre - b*(pp - pp0) ;
      sig[4] += sig0_22 + lame*tre - b*(pp - pp0) ;
      sig[8] += sig0_33 + lame*tre - b*(pp - pp0) ;
      {
	/* masse liquide */
	double m_l = rho_l*phi*sl ;
	/* coefficient de transfert */
	double k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
	/* flux */
	double w_l[3],gpl[3] ;
	grad(x,u,dh,gpl,el.nn,dim,I_p_l) ;
	for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
	w_l[dim-1] += k_l*rho_l*gravite ;
	/* rangement dans f */
	M_l = m_l ;
	for(i=0;i<3;i++) W_l[i]    = w_l[i] ;
	for(i=0;i<9;i++) SIG[i]    = sig[i] ;
	for(i=0;i<3;i++) F_MASS[i] = 0. ;
	F_MASS[dim-1] = (rho_s+m_l)*gravite ;
	PHI = phi ;
      }
    }
  }
}

int ex7(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  int    p ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    /* pressions */
    double pl  = param(u,h,el.nn,I_p_l) ;
    double pc  = p_g - pl ;
    /* coefficient de transfert */
    double k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    /* rangement dans va */
    K_l = k_l ;
  }
  return(0) ;
}

int ct7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double dmu,lame,pp0 ;
  int    i,p ;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;

  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  pp0     = pie(p_l0,p_g,el.mat->cb[0]) ;


  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,f += NVI,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    double *dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    double pl  = param(u,h,el.nn,I_p_l) ;
    double pc  = p_g - pl ;
    double pp  = pie(pl,p_g,el.mat->cb[0]) ;
    /* deformations */
    double eps[9] ;
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    {
      double tre  = eps[0] + eps[4] + eps[8] ;
      /* porosite */
      double dphi = b*tre + N*(pp - pp0) ;
      double phi  = phi0 + dphi ;
      /* saturation */
      double sl  = courbe(pc,el.mat->cb[0]) ;
      /* contraintes */
      double sig[9] ;
      for(i=0;i<9;i++) sig[i] = dmu*eps[i] ;
      sig[0] += sig0_11 + lame*tre - b*(pp - pp0) ;
      sig[4] += sig0_22 + lame*tre - b*(pp - pp0) ;
      sig[8] += sig0_33 + lame*tre - b*(pp - pp0) ;
      {
	/* masse liquide */
	double m_l = rho_l*phi*sl ;
	/* coefficient de transfert */
	double k_l = K_l ;
	/* flux */
	double w_l[3],gpl[3] ;
	grad(x,u,dh,gpl,el.nn,dim,I_p_l) ;
	for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
	w_l[dim-1] += k_l*rho_l*gravite ;
	/* rangement dans f */
	M_l = m_l ;
	for(i=0;i<3;i++) W_l[i] = w_l[i] ;
	for(i=0;i<9;i++) SIG[i] = sig[i] ;
	for(i=0;i<3;i++) F_MASS[i] = zero ;
	F_MASS[dim-1] = (rho_s+m_l)*gravite ;
	PHI = phi ;
      }
    }
  }

  return(0) ;
}

int mx7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  int    i,n,m,dec ;
  double c[MAX_PGAUSS*100] ;
  double kb[MAX_NOEUDS*MAX_NOEUDS] ;
  double zero = 0. ;

  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
  ** Matrice de comportement
  */
  dec = c7(x,u,u_n,f,f_n,va,el,dim,geom,dt,c) ;
  mxbiot(k,x,*el.fi,c,dim,dec,geom,NEQ,E_mec,I_u,E_liq,I_p_l) ;
  /*
  ** Matrice de conduction
  */
  dec = k7(x,u,u_n,f,f_n,va,el,dim,geom,dt,c) ;
  mxcond(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_liq+n*NEQ,I_p_l+m*NEQ) += dt*kb[n*el.nn+m] ;
  }

  if(strstr(el.mat->method,"P2P1")) mxnd(k,el,I_p_l,E_liq) ; /* elements P2P1 */
  return(0) ;
#undef K
}

void rs7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i,n ;
  double rb[3*MAX_NOEUDS],g1[MAX_PGAUSS] ;
  double *f_1 = f ;
  double zero = 0. ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /* 1. Mecanique */
  rscont(rb,x,*el.fi,SIG,dim,NVI,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) R(n,E_mec+i) -= rb[i+n*dim] ;
  rsmass(rb,x,*el.fi,F_MASS+dim-1,dim,NVI,geom) ;
  for(n=0;n<el.nn;n++) R(n,E_mec+dim-1) -= -rb[n] ;
  /* 2. Hydraulique */
  /* 2.1 Termes d'accumulation */
  for(i=0;i<el.fi->np;i++,f += NVI,f_n += NVI) g1[i] = M_l - M_ln ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= rb[i] ;
  /* 2.2 Termes de transport */
  f = f_1 ;
  rsflux(rb,x,*el.fi,W_l,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= -dt*rb[i] ;


  if(strstr(el.mat->method,"P2P1")) rsnd(r,el,E_liq) ; /* elements P2P1 */
#undef R
}

int so7(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define NSO  (7)
  double dmu,lame,pp0 ;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;

  dmu     = young/(un+poisson) ;
  lame    = dmu*poisson/(un-deux*poisson) ;
  pp0     = pie(p_l0,p_g,el.mat->cb[0]) ;

  /* initialisation */
  {
    int    i,j ;
    for(i=0;i<NSO;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;
  }

  {
    int    nso = 0 ;
    double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
    /* quantites exploitees en s */
    fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
    {
      int    i,p ;
      /* pressions */
      double pl  = param(u,h_s,el.nn,I_p_l) ;
      double pc  = p_g - pl ;
      double pp  = pie(pl,p_g,el.mat->cb[0]) ;
      /* saturation */
      double sl  = courbe(pc,el.mat->cb[0]) ;

      strcpy(r[nso].text,"pression-liquide") ; r[nso].n = 1 ;
      r[nso++].v[0] =  pl ;
      
      strcpy(r[nso].text,"deplacements") ; r[nso].n = 3 ;
      for(i=0;i<dim;i++) r[nso].v[i] = param(u,h_s,el.nn,I_u+i) ;
      nso++ ;
      
      /* boucle sur les points d'integration */
      for(p=0;p<el.fi->np;p++,f += NVI,va += NVE) {
	strcpy(r[nso].text,"flux-liquide") ; r[nso].n = 3 ;
	for(i=0;i<dim;i++) r[nso].v[i] += W_l[i]/el.fi->np ;
	
	strcpy(r[nso+1].text,"contraintes") ; r[nso+1].n = 9 ;
	for(i=0;i<9;i++) r[nso+1].v[i] += SIG[i]/el.fi->np ;
	
	strcpy(r[nso+2].text,"porosite") ; r[nso+2].n = 1 ;
	r[nso+2].v[0] += PHI/el.fi->np ;
      }
      nso += 3 ;
      
      strcpy(r[nso].text,"saturation") ; r[nso].n = 1 ;
      r[nso++].v[0] = sl ;
      
      strcpy(r[nso].text,"pression-pi") ; r[nso].n = 1 ;
      r[nso++].v[0] = pp ;
    }
    if(nso != NSO) arret("so7") ;
  }

  return (NSO) ;
}

double pie(double p_l,double p_g,crbe_t cb)
{
  int    j ;
  double pc,sl,sg,u ;
  int    np = 3 ;
  double a[3] = {0.93246951420,0.66120938646,0.23861918608} ;
  double w[3] = {0.17132449237,0.36076157304,0.46791393457} ;
  double zero = 0.,un = 1.,deux = 2. ;

  pc = p_g - p_l ;
  sl = courbe(pc,cb) ;
  sg = un - sl ;
  u  = zero ;
  for(j=0;j<np;j++) {
    u += w[j]*(courbe(pc/deux*(un+a[j]),cb)
       +       courbe(pc/deux*(un-a[j]),cb)) ;
  }
  u  = u*pc/deux - sl*pc ;
  return(sl*p_l + sg*p_g - 2./3.*u) ;
}

int c7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  double dmu,lame,mu ;
  int    dec = 100 ;
  int    i,j,p ;
  double zero = 0.,un = 1.,deux = 2. ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;

  dmu     = young/(un+poisson) ;
  mu      = dmu/deux ;
  lame    = dmu*poisson/(un-deux*poisson) ;

  for(p=0;p<el.fi->np;p++) {
    double *c1 = c + p*dec ;
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    /* pressions */
    double pl  = param(u,h,el.nn,I_p_l) ;
    double pc  = p_g - pl ;
    /* saturation */
    double sl  = courbe(pc,el.mat->cb[0]) ;
    double dslsdpc = dcourbe(pc,el.mat->cb[0]) ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;

    /* mecanique */
    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      C1(i,i,j,j) += lame ;
      C1(i,j,i,j) += mu ;
      C1(i,j,j,i) += mu ;
    }
    c1 += 81 ;
    for(i=0;i<3;i++) B1(i,i) = -b*(sl+pc*dslsdpc/3.) ;

    c1 += 9 ;
    /* hydraulique */
    for(i=0;i<3;i++) B1(i,i) = rho_l*sl*b ;
    c1 += 9 ;
    c1[0] = -rho_l*phi0*dslsdpc + rho_l*sl*N*(sl+pc*dslsdpc/3.) ;
  }
  return(dec) ;

#undef C1
#undef B1
}

int k7(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec = 9 ;
  int    i,p ;
  double zero = 0. ;

  for(p=0;p<el.fi->np;p++,va += NVE) {
    double *c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* permeabilite */
    c1[0] = K_l ;
    c1[4] = K_l ;
    c1[8] = K_l ;
  }
  return(dec) ;
}

#undef NEQ


void mxnd(double *k,elem_t el,int inc,int equ)
{
#define NN         (Element_GetNbOfNodes(&el))
#define DIM        (Element_GetDimension(&el))
#define NEQ        (Element_GetNbOfEquations(&el))
#define K(i,j)     (k[(i)*NN*NEQ+(j)])
  int    n_vertices = 0 ;
  int    edge_vertices_line[2]       = {0,1} ;
  int    edge_vertices_triangle[6]   = {0,1,1,2,2,0} ;
  int    edge_vertices_quadrangle[8] = {0,1,1,2,2,3,3,0} ;
  int    *edge_vertices[MAX_NOEUDS] ;
  int    ilin,icol ;
  int    i,n ;

  if(NN > MAX_NOEUDS) arret("mxnd") ;

  if(DIM == 1) {
    if(NN == 3) {
      n_vertices = 2 ;
      edge_vertices[n_vertices] = edge_vertices_line ;
    } else return ;
  } else if(DIM == 2) {
    if(NN == 6) {
      n_vertices = 3 ;
      edge_vertices[n_vertices] = edge_vertices_triangle ;
    } else if(NN == 8) {
      n_vertices = 4 ;
      edge_vertices[n_vertices] = edge_vertices_quadrangle ;
    } else return ;
  } else if(DIM == 3) {
    arret("mxnd") ;
  } else {
    arret("mxnd") ;
  }

  for(n=n_vertices+1;n<NN;n++) edge_vertices[n] = edge_vertices[n-1] + 2 ;

  /* Transformation de la matrice de degre 2 en degre 1 */
  for(ilin=0;ilin<NN*NEQ;ilin++) { /* boucle sur les lignes */
    for(n=n_vertices;n<NN;n++) {
      int    icol  = n*NEQ + inc ;
      for(i=0;i<2;i++) {
	int    m    = edge_vertices[n][i] ;
	int    jcol = m*NEQ + inc ;
	K(ilin,jcol) += 0.5*K(ilin,icol) ;
      }
      K(ilin,icol) = 0. ;
    }
  }
  for(icol=0;icol<NN*NEQ;icol++) { /* boucle sur les colonnes */
    for(n=n_vertices;n<NN;n++) {
      int    ilin  = n*NEQ + equ ;
      for(i=0;i<2;i++) {
	int    m     = edge_vertices[n][i] ;
	int    jlin  = m*NEQ + equ ;
	K(jlin,icol) += 0.5*K(ilin,icol) ;
      }
      K(ilin,icol) = 0. ;
    }
  }

  /* Relation lineaire aux noeuds milieux des aretes */
  for(n=n_vertices;n<NN;n++) {
    int    ilin  = n*NEQ + equ ;
    K(ilin,ilin) = -1 ;
    for(i=0;i<2;i++) {
      int    m     = edge_vertices[n][i] ;
      int    icol  = m*NEQ + inc ;
      K(ilin,icol) = 0.5 ;
    }
  }
#undef NN
#undef NEQ
#undef DIM
#undef K
}



void rsnd(double *r,elem_t el,int equ)
{
#define NN         (Element_GetNbOfNodes(&el))
#define DIM        (Element_GetDimension(&el))
#define NEQ        (Element_GetNbOfEquations(&el))
  int    n_vertices = 0 ;
  int    edge_vertices_line[2]       = {0,1} ;
  int    edge_vertices_triangle[6]   = {0,1,1,2,2,0} ;
  int    edge_vertices_quadrangle[8] = {0,1,1,2,2,3,3,0} ;
  int    *edge_vertices[MAX_NOEUDS] ;
  int    i,n ;

  if(NN > MAX_NOEUDS) arret("rsnd") ;

  if(DIM == 1) {
    if(NN == 3) {
      n_vertices = 2 ;
      edge_vertices[n_vertices] = edge_vertices_line ;
    } else return ;
  } else if(DIM == 2) {
    if(NN == 6) {
      n_vertices = 3 ;
      edge_vertices[n_vertices] = edge_vertices_triangle ;
    } else if(NN == 8) {
      n_vertices = 4 ;
      edge_vertices[n_vertices] = edge_vertices_quadrangle ;
    } else return ;
  } else if(DIM == 3) {
    arret("rsnd") ;
  } else {
    arret("rsnd") ;
  }

  for(n=n_vertices+1;n<NN;n++) edge_vertices[n] = edge_vertices[n-1] + 2 ;

  for(n=n_vertices;n<NN;n++) {
    int    ilin  = n*NEQ + equ ;
    for(i=0;i<2;i++) {
      int    m     = edge_vertices[n][i] ;
      int    jlin  = m*NEQ + equ ;
      r[jlin] += 0.5*r[ilin] ;
    }
    r[ilin] = 0. ;
  }

#undef NN
#undef NEQ
#undef DIM
}
