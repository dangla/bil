#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Common.h"

#define MODELINDEX  60
#define TITLE   "Poroelasticity with surface adsorption"
#define AUTHORS "Nikoosokhan"

#include "OldMethods.h"

/* Macros */
#define NEQ     (1+dim)
#define NVI     (17)
#define NVE     (1)

#define E_co2   (dim)
#define E_mec   (0)

#define I_p_co2 (dim)
#define I_u     (0)

#define N_CO2   (f[0])
#define W_CO2   (f + 1)
#define SIG     (f + 4)
#define F_MASS  (f + 13)
#define PHI     (f[16])

#define N_CO2n  (f_n[0])

#define K_CO2   (va[0])

/* constantes physiques */
#define RT        (2436.)   /* produit de R=8.3143 et T=293 (Pa.m3/mole) */
/* masse molaire (kg/mole) */
#define M_CO2     (44.e-3)

/* Fonctions */
static int    pm(const char*) ;
static int    c60(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static int    k60(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t,double,double*) ;
static double concentration_ads(double,double,double) ;
static double contrainte_sup(double,double,double) ;
/* Parametres */
static double gravite ;
static double young,poisson,b,N ;
static double k_int,mu_co2 ;
static double omega0,phi0,rho_s ;
static double p_co20,sig0_11,sig0_22,sig0_33 ;
static double a0,Gamma_max ;
static double c_eps,c_phi ;

#define CONCENTRATION_ADS(p)     concentration_ads(p,a0,Gamma_max)
#define CONTRAINTE_SUP(p)        contrainte_sup(p,a0,Gamma_max)

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"young") == 0) return (1) ;
  else if(strcmp(s,"poisson") == 0) return (2) ;
  else if(strcmp(s,"phi") == 0) return (3) ;
  else if(strcmp(s,"omega") == 0) return (4) ;
  else if(strcmp(s,"p_co20") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_co2") == 0) return (7) ;
  else if(strcmp(s,"b") == 0) return (8) ;
  else if(strcmp(s,"N") == 0) return (9) ;
  else if(strcmp(s,"rho_s") == 0) return (10) ;
  else if(strcmp(s,"sig0_11") == 0) return (11) ;
  else if(strcmp(s,"sig0_22") == 0) return (12) ;
  else if(strcmp(s,"sig0_33") == 0) return (13) ;
  else if(strcmp(s,"a0") == 0) return (14) ;
  else if(strcmp(s,"Gamma_max") == 0) return (15) ;
  else if(strcmp(s,"c_eps") == 0) return (16) ;
  else if(strcmp(s,"c_phi") == 0) return (17) ;
  else return(-1) ;
}

int dm60(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 18,n_courbes = 0 ;
  int i ;

  mat->neq = NEQ ;

  for(i=0;i<dim;i++) {
    sprintf(mat->eqn[E_mec + i],"meca_%d",i+1) ;
    sprintf(mat->inc[I_u + i],"u_%d",i+1) ;
  }
  strcpy(mat->eqn[E_co2],"co2") ;
  strcpy(mat->inc[I_p_co2],"p_co2") ;

  /* Par defaut tout a 0 */
  for(i=0;i<n_donnees;i++) mat->pr[i] = 0. ;
  { /* initialisation automatique (a completer) SAEID */
    mat->pr[pm("b")]       = 1. ;
    mat->pr[pm("mu_co2")]  = 1.e-3 ;
    mat->pr[pm("Gamma_max")]  = 2.e-7 ; /* M. Vandamme et al (2010) */
  }

  lit_mate(mat,ficd,pm,n_donnees,n_courbes) ;

  /* Pour creer une courbe dans un fichier temporaire */
  /*
  {
    char   *nomtmp = (char *) tempnam(".","curve") ;
    char   line[500] ;
    a0      = mat->pr[pm("a0")] ;
    Gamma_max  = mat->pr[pm("Gamma_max")] ;
    sprintf(line,"Courbe = %s p_co2 = Range{x1 = 0,x2 = %e,n = 101} Gamma = Langmuir(1){Gmax = %e,a = %e}",nomtmp,10*a0,Gamma_max,a0) ;
    free(nomtmp) ;
    ecrit_courbe(line) ;
    mat->nc += lit_courbe(mat,line) ;
  }
  */

  return(mat->n) ;
}

int qm60(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
The system consists in 1+dim equations:\n\
\t 1. Mass Balance Equation (p_co2)\n\
\t 2. Mechanical Equilibrium (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;

  fprintf(ficd,"gravite = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 0         # masse density of the dry material\n") ;
  fprintf(ficd,"young = 0.8e9     # Young's modulus\n") ;
  fprintf(ficd,"poisson = 0.25    # Poisson's ratio\n") ;
  fprintf(ficd,"phi = 0.3         # porosity\n") ;
  fprintf(ficd,"p_co20 = 0        # initial pressure of CO2\n") ;
  fprintf(ficd,"k_int = 1.e-12    # intrinsic permeability\n") ;
  fprintf(ficd,"mu_co2 = 1.e-3    # viscosity of CO2\n") ;
  fprintf(ficd,"b = 1             # Biot's coefficient\n") ;
  fprintf(ficd,"N = 0             # compressibility of pores\n") ;
  fprintf(ficd,"sig0_11 = 0       # initial stress 11\n") ;
  fprintf(ficd,"sig0_22 = 0       # initial stress 22\n") ;
  fprintf(ficd,"sig0_33 = 0       # initial stress 33\n") ;
  fprintf(ficd,"omega = 0         # solid surface per unit volume\n") ;
  fprintf(ficd,"a0 = 1.e6         # pressure at Gamma_max/2\n") ;
  fprintf(ficd,"Gamma_max = 2.e-7 # max adsorbed concentration (for p = inf)\n") ;
  fprintf(ficd,"c_eps = 0         # derivative d_omega/d_eps\n") ;
  fprintf(ficd,"c_phi = 0         # derivative d_omega/d_phi\n") ;
  
  return(NEQ) ;
}

void tb60(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}

void ch60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  /* hydraulique */
  if(Element_FindEquationPositionIndex(&el,cg.eqn) == E_co2) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }
}

void in60(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double dmu,lame ;
  int    p ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_co2  = el.mat->pr[pm("mu_co2")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_co20  = el.mat->pr[pm("p_co20")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;
  omega0  = el.mat->pr[pm("omega")] ;
  a0      = el.mat->pr[pm("a0")] ;
  Gamma_max     = el.mat->pr[pm("Gamma_max")] ;
  c_eps   = el.mat->pr[pm("c_eps")] ;
  c_phi   = el.mat->pr[pm("c_phi")] ;

  dmu     = young/(1 + poisson) ;
  lame    = dmu*poisson/(1 - 2*poisson) ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,f += NVI,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    double *dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    double p_co2  = param(u,h,el.nn,I_p_co2) ;
    /* contrainte superficielle */
    double sig_s   = CONTRAINTE_SUP(p_co2) ;
    double sig_s0  = CONTRAINTE_SUP(p_co20) ;
    double dsig_s  = sig_s - sig_s0 ;
    /* pression effective */
    double pp_co2  = p_co2  - c_phi*sig_s ;
    double pp_co20 = p_co20 - c_phi*sig_s0 ;
    double dpp_co2 = pp_co2 - pp_co20 ;
    /* concentration adsorbee */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    /* densite molaire en co2 */
    double rho_co2 = p_co2/RT ;
    /* deformations */
    double eps[9] ;
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    {
      double tre  = eps[0] + eps[4] + eps[8] ;
      /* porosite */
      double dphi = b*tre + N*dpp_co2 ;
      double phi  = phi0 + dphi ;
      /* contraintes */
      double sig[9] ;
      int    i ;
      for(i=0;i<9;i++) sig[i] = dmu*eps[i] ;
      sig[0] += sig0_11 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      sig[4] += sig0_22 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      sig[8] += sig0_33 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      {
	/* surface solide interne */
	double omega   = omega0 + c_eps*tre + c_phi*dphi ;
	/* contenu molaire en co2 */
	double n_co2   = rho_co2*phi + c_co2*omega ;
	/* coefficient de transfert */
	double k_co2   = rho_co2*k_int/mu_co2 ;
	/* flux */
	double w_co2[3],grd_co2[3] ;
	grad(x,u,dh,grd_co2,el.nn,dim,I_p_co2) ;
	for(i=0;i<3;i++) w_co2[i] = - k_co2*grd_co2[i] ;
	w_co2[dim-1] += k_co2*rho_co2*gravite ;
	/* rangement dans f */
	N_CO2 = n_co2 ;
	for(i=0;i<3;i++) W_CO2[i]  = w_co2[i] ;
	for(i=0;i<9;i++) SIG[i]    = sig[i] ;
	for(i=0;i<3;i++) F_MASS[i] = 0. ;
	F_MASS[dim-1] = (rho_s + M_CO2*n_co2)*gravite ;
	PHI = phi ;
	/* rangement dans va */
	K_CO2 = k_co2 ;
      }
    }
  }
}

int ex60(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  int    p ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_co2  = el.mat->pr[pm("mu_co2")] ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    /* pression */
    double p_co2  = param(u,h,el.nn,I_p_co2) ;
    /* densite molaire en co2 */
    double rho_co2 = p_co2/RT ;
    /* coefficient de transfert */
    double k_co2 = rho_co2*k_int/mu_co2 ;
    /* rangement dans va */
    K_CO2 = k_co2 ;
  }

  return(0) ;
}

int ct60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double dmu,lame ;
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
  mu_co2  = el.mat->pr[pm("mu_co2")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_co20  = el.mat->pr[pm("p_co20")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;
  omega0  = el.mat->pr[pm("omega")] ;
  a0      = el.mat->pr[pm("a0")] ;
  Gamma_max     = el.mat->pr[pm("Gamma_max")] ;
  c_eps   = el.mat->pr[pm("c_eps")] ;
  c_phi   = el.mat->pr[pm("c_phi")] ;

  dmu     = young/(1 + poisson) ;
  lame    = dmu*poisson/(1 - 2*poisson) ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++,f += NVI,va += NVE) {
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    double *dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    double p_co2  = param(u,h,el.nn,I_p_co2) ;
    /* contrainte superficielle */
    double sig_s   = CONTRAINTE_SUP(p_co2) ;
    double sig_s0  = CONTRAINTE_SUP(p_co20) ;
    double dsig_s  = sig_s - sig_s0 ;
    /* pression effective */
    double pp_co2  = p_co2  - c_phi*sig_s ;
    double pp_co20 = p_co20 - c_phi*sig_s0 ;
    double dpp_co2 = pp_co2 - pp_co20 ;
    /* concentration adsorbee */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    /* densite molaire en co2 */
    double rho_co2 = p_co2/RT ;
    /* deformations */
    double eps[9] ;
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    {
      double tre  = eps[0] + eps[4] + eps[8] ;
      /* porosite */
      double dphi = b*tre + N*dpp_co2 ;
      double phi  = phi0 + dphi ;
      /* contraintes */
      double sig[9] ;
      int    i ;
      for(i=0;i<9;i++) sig[i] = dmu*eps[i] ;
      sig[0] += sig0_11 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      sig[4] += sig0_22 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      sig[8] += sig0_33 + lame*tre - b*dpp_co2 + c_eps*dsig_s ;
      {
	/* surface solide interne */
	double omega   = omega0 + c_eps*tre + c_phi*dphi ;
	/* contenu molaire co2 */
	double n_co2   = rho_co2*phi + c_co2*omega ;
	/* coefficient de transfert */
	double k_co2   = K_CO2 ;
	/* flux */
	double w_co2[3],grd_co2[3] ;
	grad(x,u,dh,grd_co2,el.nn,dim,I_p_co2) ;
	for(i=0;i<3;i++) w_co2[i] = - k_co2*grd_co2[i] ;
	w_co2[dim-1] += k_co2*rho_co2*gravite ;
	/* rangement dans f */
	N_CO2 = n_co2 ;
	for(i=0;i<3;i++) W_CO2[i]  = w_co2[i] ;
	for(i=0;i<9;i++) SIG[i]    = sig[i] ;
	for(i=0;i<3;i++) F_MASS[i] = 0. ;
	F_MASS[dim-1] = (rho_s + M_CO2*n_co2)*gravite ;
	PHI = phi ;
      }
    }
  }

  return(0) ;
}

int mx60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  int    i,n,m,dec ;
  double c[MAX_PGAUSS*100] ;
  double kb[MAX_NOEUDS*MAX_NOEUDS] ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
  ** Matrice de comportement
  */
  dec = c60(x,u,u_n,f,f_n,va,el,dim,geom,dt,c) ;
  mxbiot(k,x,*el.fi,c,dim,dec,geom,NEQ,E_mec,I_u,E_co2,I_p_co2) ;
  /*
  ** Matrice de conduction
  */
  dec = k60(x,u,u_n,f,f_n,va,el,dim,geom,dt,c) ;
  mxcond(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_co2+n*NEQ,I_p_co2+m*NEQ) += dt*kb[n*el.nn+m] ;
  }

  return(0) ;
#undef K
}

void rs60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i,n ;
  double rb[3*MAX_NOEUDS],g1[MAX_PGAUSS] ;
  double *f_1 = f ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0. ;

  if(el.dim < dim) return ;

  /* 1. Mecanique */
  /* 1.1 Contraintes */
  rscont(rb,x,*el.fi,SIG,dim,NVI,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) R(n,E_mec+i) -= rb[i+n*dim] ;
  /* 1.2 Forces de masse */
  rsmass(rb,x,*el.fi,F_MASS+dim-1,dim,NVI,geom) ;
  for(n=0;n<el.nn;n++) R(n,E_mec+dim-1) -= -rb[n] ;

  /* 2. Hydraulique */
  /* 2.1 Termes d'accumulation */
  for(i=0;i<el.fi->np;i++,f += NVI,f_n += NVI) g1[i] = N_CO2 - N_CO2n ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_co2) -= rb[i] ;
  /* 2.2 Termes de transport */
  f = f_1 ;
  rsflux(rb,x,*el.fi,W_CO2,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_co2) -= -dt*rb[i] ;

#undef R
}

int so60(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define NSO   (8)

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_co2  = el.mat->pr[pm("mu_co2")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_co20  = el.mat->pr[pm("p_co20")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  sig0_11 = el.mat->pr[pm("sig0_11")] ;
  sig0_22 = el.mat->pr[pm("sig0_22")] ;
  sig0_33 = el.mat->pr[pm("sig0_33")] ;
  omega0  = el.mat->pr[pm("omega")] ;
  a0      = el.mat->pr[pm("a0")] ;
  Gamma_max     = el.mat->pr[pm("Gamma_max")] ;
  c_eps   = el.mat->pr[pm("c_eps")] ;
  c_phi   = el.mat->pr[pm("c_phi")] ;

  /* initialisation */
  {
    int    i,j ;
    for(i=0;i<NSO;i++) for(j=0;j<9;j++) r[i].v[j] = 0. ;
  }

  {
    int    nso = 0 ;
    double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
    /* quantites exploitees en s */
    fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
    {
      int    i,p ;
      /* pressions */
      double p_co2   = param(u,h_s,el.nn,I_p_co2) ;
      /* concentration adsorbe */
      double c_co2   = CONCENTRATION_ADS(p_co2) ;
      /* densite molaire en co2 */
      double rho_co2 = p_co2/RT ;
      /* contrainte superficielle */
      double sig_s   = CONTRAINTE_SUP(p_co2) ;
	
      strcpy(r[nso].text,"pression-co2") ; r[nso].n = 1 ;
      r[nso++].v[0] =  p_co2 ;
      
      strcpy(r[nso].text,"deplacements") ; r[nso].n = 3 ;
      for(i=0;i<dim;i++) r[nso].v[i] = param(u,h_s,el.nn,I_u+i) ;
      nso++ ;
      
      strcpy(r[nso].text,"concentration_adsorbee") ; r[nso].n = 1 ;
      r[nso++].v[0] =  c_co2 ;
      
      strcpy(r[nso].text,"contrainte_superficielle") ; r[nso].n = 1 ;
      r[nso++].v[0] =  sig_s ;
      
      strcpy(r[nso].text,"densite_molaire_co2") ; r[nso].n = 1 ;
      r[nso++].v[0] =  rho_co2 ;
      
      /* boucle sur les points d'integration */
      for(p=0;p<el.fi->np;p++,f += NVI) {
	strcpy(r[nso].text,"porosite") ; r[nso].n = 1 ;
	r[nso].v[0] += PHI/el.fi->np ;
	
	strcpy(r[nso+1].text,"contraintes") ; r[nso+1].n = 9 ;
	for(i=0;i<9;i++) r[nso+1].v[i] += SIG[i]/el.fi->np ;
	
	strcpy(r[nso+2].text,"flux") ; r[nso+2].n = 3 ;
	for(i=0;i<3;i++) r[nso+2].v[i] += W_CO2[i]/el.fi->np ;
      }
    }
    nso += 3 ;
    if(nso != NSO) arret("so60") ;
  }

  return (NSO) ;
}

int c60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  double dmu,lame,mu ;
  int    dec = 100 ;
  int    p ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  young   = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_co2  = el.mat->pr[pm("mu_co2")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_co20  = el.mat->pr[pm("p_co20")] ;
  b       = el.mat->pr[pm("b")] ;
  N       = el.mat->pr[pm("N")] ;
  omega0  = el.mat->pr[pm("omega")] ;
  a0      = el.mat->pr[pm("a0")] ;
  Gamma_max     = el.mat->pr[pm("Gamma_max")] ;
  c_eps   = el.mat->pr[pm("c_eps")] ;
  c_phi   = el.mat->pr[pm("c_phi")] ;

  dmu     = young/(1 + poisson) ;
  mu      = dmu/2 ;
  lame    = dmu*poisson/(1 - 2*poisson) ;
  
  for(p=0;p<el.fi->np;p++) {
    int    i,j ;
    double *c1 = c + p*dec ;
    /* fonctions d'interpolation */
    double *h  = el.fi->h  + p*el.nn ;
    double *dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    double p_co2  = param(u,h,el.nn,I_p_co2) ;
    /* contrainte superficielle */
    double sig_s   = CONTRAINTE_SUP(p_co2) ;
    double sig_s0  = CONTRAINTE_SUP(p_co20) ;
    /* pression effective */
    double pp_co2  = p_co2  - sig_s*c_phi ;
    double pp_co20 = p_co20 - sig_s0*c_phi ;
    /* densite molaire en co2 */
    double rho_co2 = p_co2/RT ;
    /* concentration adsorbe */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    /* deformations */
    double eps[9] ;
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    {
      double tre  = eps[0] + eps[4] + eps[8] ;
      /* porosite */
      double dphi = b*tre + N*(pp_co2 - pp_co20) ;
      double phi  = phi0 + dphi ;
      /* surface solide interne */
      double omega   = omega0 + c_eps*tre + c_phi*dphi ;

      /* derivees */
      /* ... par rapport a p_co2 */
      double dp_co2 = 10. ; /* en Pa */
      double p_co22 = p_co2 + dp_co2 ;
      double sig_s2 = CONTRAINTE_SUP(p_co22) ;
      double c_co22 = CONCENTRATION_ADS(p_co22) ;
      double dsig_ssdp_co2   = (sig_s2 - sig_s)/dp_co2 ;
      double dpp_co2sdp_co2  = 1. - c_phi*dsig_ssdp_co2 ;
      double drho_co2sdp_co2 = 1./RT ;
      double dphisdp_co2     = N*dpp_co2sdp_co2 ;
      double domegasdp_co2   = c_phi*dphisdp_co2 ;
      double dc_co2sdp_co2   = (c_co22 - c_co2)/dp_co2 ;
      
      /* initialisation */
      for(i=0;i<dec;i++) c1[i] = 0. ;
      
      /* mecanique */
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	  C1(i,i,j,j) += lame ;
	  C1(i,j,i,j) += mu ;
	  C1(i,j,j,i) += mu ;
	}
      c1 += 81 ;
      for(i=0;i<3;i++) B1(i,i) = -b*dpp_co2sdp_co2 + c_eps*dsig_ssdp_co2  ;
      
      c1 += 9 ;
      /* hydraulique */
      for(i=0;i<3;i++) B1(i,i) = rho_co2*b + c_co2*c_eps ;
      c1 += 9 ;
      c1[0] = rho_co2*dphisdp_co2 + drho_co2sdp_co2*phi + dc_co2sdp_co2*omega + c_co2*domegasdp_co2 ;
    }
  }
  return(dec) ;
  
#undef C1
#undef B1
}

int k60(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec = 9 ;
  int    p ;
  
  for(p=0;p<el.fi->np;p++,va += NVE) {
    double *c1 = c + p*dec ;
    int    i ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = 0. ;
    /* permeabilite */
    c1[0] = K_CO2 ;
    c1[4] = K_CO2 ;
    c1[8] = K_CO2 ;
  }
  return(dec) ;
}

double concentration_ads(double p,double a,double Gamma_max)
{
  return(Gamma_max*p/(a + p)) ;
}

double contrainte_sup(double p,double a,double Gamma_max)
{
  return(-Gamma_max*RT*log(1 + p/a)) ;
}
