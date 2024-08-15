/*
	Sel NaCl
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

#define MODELINDEX  27
#define TITLE "Sechage isotherme avec sels (Nov. 2008)"
#define AUTHORS "Nguyen"

#include "OldMethods.h"

/* Macros */
#define NEQ 	(2)
#define NVI     (10)
#define NVE     (12)
#define NVE0    (6)
#define NV0     (2)

#define E_eau	(0)
#define E_Cl	(1)

#define I_H_r   (0)
#define I_C_Cl  (1)

#define H_r(n)      (u[(n)][I_H_r])
#define C_Cl_t(n)   (u[(n)][I_C_Cl])

#define M_W(n)      (f[(n)])
#define N_Cl(n)     (f[(n+2)])

#define W_W         (f[4])
#define W_Cl        (f[5])

#define C_Cl_l(n)   (f[(6+n)])
#define P_l(n)      (f[(8+n)])

#define KD_W        (va[(0)])
#define KF_V        (va[(1)])

#define KF_Cl       (va[(2)])
#define KF_Na       (va[(3)])

#define KD_Cl       (va[(4)])
#define KD_Na       (va[(5)])

#define LNA_W(n)    (va[(6+n)])
#define PHI(n)      (va[(8+n)])


#define PHI0(n)     (v0[(0+n)])

/* valences */
#define z_cl      (-1.)
#define z_na      (1.)

/* volumes molaires partiels liquides (m3/mole) */
#define v_h2o     (1.80e-5)
#define v_na      (1.87e-5)
#define v_cl      (2.52e-6)

/* volumes molaires solides (m3/mole) */
#define v_nacl    (24.5e-6)

/* coefficients de diffusion moleculaire (m2/s) */
#define do_cl     (2.032e-9)
#define do_na     (1.334e-9)
#define do_va     (2.42e-5)

/* constante d'equilibre */
#define K_nacl    (6.e3)      /* Solubilite de NaCl (moles/m3) */

/* constante physique */
#define FARADAY   (9.64846e4) /* Faraday (C/mole) */
#define TEMPERATURE         (293.)      /* Temperature (K) */
#define RT        (2436.)     /* Produit R = 8.3143 et T = 293. (J/mole) */

/* viscosites (Pa.s) */
#define mu_g      (1.8e-5)
#define mu_l      (1.002e-3)

/* Masses molaires (kg/mole) */
#define M_h2o     (1.8e-2)
#define M_air     (2.896e-2)

/* autres */
#define p_atm     (1.01325e5) /* Pression atmospherique (Pa) */
#define p_g       (0.)        /* Pression du gaz (Pa) */

/* Pression de vapeur (Pa) */
#define P_VS(T)	  (609.14*pow(10.,7.45*(T - 273.)/(T - 38.)))

/* Fonctions */
static int    pm(const char *s) ;
static double activite(double,double,double) ;
static double lna_i(double,double,double,double,double,double) ;
extern double lng_LinLee(double,double,double,double,double,double) ;
extern double lng_TQN(double,double,double,double,double,double,double,double) ;


static double tortuosite_l(double) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
static double concentration_cl_l(double) ;
static double porosite(double,double,double) ;


/* Parametres */
static double phi0,r_d,k_int ;
static double d_cl ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1) ;  
  else if(strcmp(s,"r_d") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"courbes") == 0) return (4) ;
  else {
    printf("donnee \"%s\" non connue (pm27)\n",s) ; exit(0) ;
  }
}

int dm27(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    n_donnees = 5 ;

  if(dim > 1) arret("dm27 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;

  strcpy(mat->eqn[E_eau],"liq") ;
  strcpy(mat->eqn[E_Cl],"sel") ;

  strcpy(mat->inc[I_H_r],"h_r") ;
  strcpy(mat->inc[I_C_Cl],"c_cl") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}

int qm27(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 2 equations :\n\
\t 1. Conservation de la masse d\'eau  (h_r)\n\
\t 2. Conservation de la masse de NaCl (c_cl)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.12  # Porosite\n") ;
  fprintf(ficd,"k_int = 1.e-20   # Permeabilite intrinseque (m2)\n") ;
  fprintf(ficd,"D_Cl = 6.25e-12  # Diffusion effective de Cl (m2/s)\n") ;
  fprintf(ficd,"r_d = 1.         # Rapport des tortuosites anions/cations\n") ;
  fprintf(ficd,"courbes = myfile # nom du fichier : p_c S_l k_rl tau_l\n") ;

  return(NEQ) ;
}

void tb27(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ; /* implicite */
  el->n_ve = NVE ; /* explicite */
  Element_GetNbOfConstantTerms(el) = NV0 ;
}

void ch27(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define MAT     (el.mat)
  double* v0 = Element_GetConstantTerm(&el) ;

  int    i ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0. ;

  /*
    Donnees
  */
  phi0     = el.mat->pr[pm("porosite")] ;
  phi0     = PHI0(0) ;

  /* flux d'echange */
  if(!strncmp(cg.t,"echange",4)) {
    /* seulement sur l'equation de conservation de l'eau */
    if(isdigit(cg.eqn[0])) { /* donne sous forme numerique */
      double ieq  = atoi(cg.eqn) - 1 ;
      if(ieq != E_eau) arret("ch27 : non prevu") ;
    } else {                 /* donne sous forme alphabetique */
      if(strcmp(cg.eqn,MAT->eqn[E_eau])) arret("ch27 : non prevu") ;
    }
    {
      double c_cl_t = C_Cl_t(0) ;  
      double h_r    = H_r(0) ;
      double p_vs   = P_VS(TEMPERATURE) ;  
      double p_v    = h_r*p_vs ; 
      double c_cl,c_na,c_w ;
      double s_w ;
      double p_l,p_c ;
      double phi ;
      double lna_w ;
      double rho_v,rho_w ;
      double h_ext = champ(x[0],dim,*cg.ch)*fonction(t,*cg.fn) ;
      double dx    = 1.e-2 ; /* humidite h_ext a une distance dx = 1 cm */
      
      /* concentrations */
      c_cl    = concentration_cl_l(c_cl_t) ;
      c_na    = c_cl ;
      c_w     = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;
      
      /* activite de l'eau */
      lna_w = activite(c_cl,c_na,c_w) ;
      
      /* pressions */
      p_l = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
      p_c = p_g - p_l ;
      
      /* masses volumiques */
      rho_v = M_h2o/RT*p_v ;
      rho_w = M_h2o*c_w ;
      
      /* saturations */
      s_w  = courbe(p_c,el.mat->cb[0]) ;
      
      /* porosite */
      phi  = porosite(phi0,s_w,c_cl_t) ;
      
      r[E_eau] = - dt*do_va*M_h2o/RT*p_vs*(h_ext - h_r)/dx ;
      r[E_eau] *= pow(phi/phi0,3)*pow((1 - phi0)/(1 - phi),2) ;
    }
  } else {
    chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  }

  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;

#undef MAT
}

void in27(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double* v0 = Element_GetConstantTerm(&el) ;
  int    i ;

  /*
    Donnees
  */
  phi0     = el.mat->pr[pm("porosite")] ;
  if(phi0 < 0.) {
    int   fct_phi = - floor(phi0 + 0.5) - 1 ;
    Field_t *ch = Material_GetField(Element_GetMaterial(&el)) ;
    for(i=0;i<el.nn;i++) PHI0(i) = champ(x[i],dim,ch[fct_phi]) ;
  } else {
    for(i=0;i<el.nn;i++) PHI0(i) = phi0 ;
  }

  if(el.dim < dim) return ;

  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double c_cl_t = C_Cl_t(i) ;  
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ;  
    double p_v    = h_r*p_vs ;
    double c_cl,c_na,c_w ;
    double s_w,s_g ;
    double p_l,p_c ;
    double phi ;
    double lna_w ;
    double rho_v,rho_w ;

    /* concentrations */
    c_cl    = concentration_cl_l(c_cl_t) ;
    c_na    = c_cl ;
    c_w     = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;
  
    /* activite de l'eau */
    lna_w = activite(c_cl,c_na,c_w) ;

    /* pressions */
    if(h_r < 0.) {
      s_w = - h_r ;
      if(el.mat->nc > 2) {
	p_c = courbe(s_w,el.mat->cb[3]) ;
      } else {
	arret("in27 : nombre de courbes insuffisant") ;
	return ;
      }
      h_r = exp(lna_w - v_h2o/RT*p_c) ;  
      p_v = h_r*p_vs ;
    }
    p_l = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
    p_c = p_g - p_l ;

    /* masses volumiques */
    rho_v = M_h2o/RT*p_v ;
    rho_w = M_h2o*c_w ;

    /* saturations */
    s_w  = courbe(p_c,el.mat->cb[0]) ;
    s_g  = 1. - s_w ;
        
    /* porosite */
    phi0 = PHI0(i) ;
    phi  = porosite(phi0,s_w,c_cl_t) ;

    /* contenus molaires */
    N_Cl(i)  = c_cl_t*phi0*s_w ;
    M_W(i)   = rho_w*phi*s_w + rho_v*phi*s_g ;

    /* sauvegarde */
    C_Cl_l(i) = c_cl ;
    P_l(i)    = p_l ;
    H_r(i)    = h_r ;
  }

  /* Coefficients de transfert */
  {
    ex_t ex27 ;
    ex27(x,u,f,va,el,dim,geom,0.) ; 
  }

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex27(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double* v0 = Element_GetConstantTerm(&el) ;
  int i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi0     = el.mat->pr[pm("porosite")] ;
  d_cl     = el.mat->pr[pm("D_Cl")] ;
  r_d      = el.mat->pr[pm("r_d")] ;
  k_int	   = el.mat->pr[pm("k_int")] ;

  for(i=0;i<NVE0;i++) va[i] = 0. ;

  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double c_cl_t = C_Cl_t(i) ;  
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ;
    double c_cl,c_na,c_w ;
    double s_w,s_g ;
    double p_l,p_c ;
    double phi ;
    double lna_w ;
    double tau_g,tau_l,tau_ani,tau_cat ;
    double k_rl,k_h ;

    /* concentrations */
    c_cl = concentration_cl_l(c_cl_t) ;
    c_na = c_cl ;
    c_w  = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;

    /* activite de l'eau */ 
    lna_w = activite(c_cl,c_na,c_w) ;	
	 
    /* pressions */
    p_l = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
    p_c = p_g - p_l ;

    /* saturations */
    s_w  = courbe(p_c,el.mat->cb[0]) ;
    s_g  = 1. - s_w ;

    /* porosite */
    phi0 = PHI0(i) ;
    phi  = porosite(phi0,s_w,c_cl_t) ;

    /* permeabilite */
    k_rl    = courbe(p_c,el.mat->cb[1]) ;
    k_h     = k_int/mu_l*k_rl ;
    k_h    *= pow(phi/phi0,3)*pow((1. - phi0)/(1. - phi),2) ;

    /* tortuosites gaz et liquide*/
    tau_g   = pow(phi,4./3.)*pow(s_g,10./3.) ;
    tau_l   = tortuosite_l(phi)*courbe(p_c,el.mat->cb[2]) ;
  
    /* tortuosites anions et cations */
    tau_ani = tau_l*d_cl/(tortuosite_l(phi0)*do_cl) ;
    tau_cat = tau_ani/r_d ;
  
    /* sauvegarde */
    KD_W   += M_h2o*c_w*k_h ;
    KD_Cl  += c_cl*k_h ;
    KD_Na  += c_na*k_h ;
    
    KF_V   += tau_g*do_va*M_h2o/RT*p_vs ;
    KF_Cl  += tau_ani*do_cl ;
    KF_Na  += tau_cat*do_na ;

    LNA_W(i) = lna_w ;
    PHI(i)   = phi ;
  }

  KD_W  *= 0.5 ;
  KD_Cl *= 0.5 ;
  KD_Na *= 0.5 ;

  KF_V  *= 0.5 ;
  KF_Cl *= 0.5 ;
  KF_Na *= 0.5 ;

  return(0) ;
}

int ct27(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* v0 = Element_GetConstantTerm(&el) ;
  int    i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */

  phi0     = el.mat->pr[pm("porosite")] ;
   
  /* Contenus molaires */

  for(i=0;i<2;i++) {
    double c_cl_t = C_Cl_t(i) ; 
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ; 
    double p_v    = h_r*p_vs ;
    double c_cl,c_na,c_w ;
    double s_w,s_g ;
    double p_c,p_l ;
    double phi ;
    double lna_w ;
    double rho_v,rho_w ;

    /* concentrations */
    c_cl    = concentration_cl_l(c_cl_t) ;
    c_na    = c_cl ;
    c_w     = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;
  
    /* activite de l'eau */ 
    lna_w = activite(c_cl,c_na,c_w) ;

    /* pressions */
    p_l = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
    p_c = p_g - p_l ;

    /* masses volumiques */
    rho_v = M_h2o/RT*p_v ;
    rho_w = M_h2o*c_w ;

    /* saturations */
    s_w  = courbe(p_c,el.mat->cb[0]) ;
    s_g  = 1. - s_w ;
        
    /* porosite */
    phi0 = PHI0(i) ;
    phi  = porosite(phi0,s_w,c_cl_t) ;
    /* phi  = PHI(i) ; */

    /* contenus molaires */
    N_Cl(i)  = c_cl_t*phi0*s_w ;
    M_W(i)   = rho_w*phi*s_w + rho_v*phi*s_g ;

    /* sauvegarde */
    C_Cl_l(i) = c_cl ;
    P_l(i)    = p_l ;

    if(h_r <= 0 || phi < 0.) {
      printf("\n\
x       = %e\n\
h_r     = %e\n\
p_v     = %e\n\
c_cl_t  = %e\n\
s_w     = %e\n\
phi     = %e\n",x[i][0],h_r,p_v,c_cl_t,s_w,phi) ;
      return(-1) ;
    }
  }
  
  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
  
  return(0) ;
}

int mx27(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  double* v0 = Element_GetConstantTerm(&el) ;

  double dc_clsdc_cl_t[2],dp_lsdh_r[2],dp_lsdc_cl_t[2],dpsisdc_cl,c[2] ;
  double trd_w,trd_cl,trf_v,trf_cl,tre_cl ;
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  
  phi0     = el.mat->pr[pm("porosite")] ;

  /*
    CALCUL DE volume ET DE surf 
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = 2*M_PI*xm ; else surf = 1. ;

  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    double c_cl_t = C_Cl_t(i) ; 
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ; 
    double p_v    = h_r*p_vs ;
    double p_c,p_l ;
    double s_w,s_g ;
    double c_cl,c_na,c_w ;
    double phi ;
    double lna_w ;
    double rho_v,rho_w ;
    double ds_wsdp_c,ds_wsdh_r,ds_gsdh_r,ds_wsdc_cl_t,ds_gsdc_cl_t ;
    double drho_vsdh_r,drho_wsdc_cl_t ;
    double dphisdc_cl_t,dphisdh_r ;
    double dlna_wsdc_cl ;

    /* concentrations */
    c_cl    = concentration_cl_l(c_cl_t) ;
    c_na    = c_cl ;
    c_w     = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;
  
    /* activite de l'eau */ 
    lna_w = activite(c_cl,c_na,c_w) ;

    /* pressions */
    p_l = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
    p_c = p_g - p_l ;

    /*masses volumiques */
    rho_v = M_h2o/RT*p_v ;
    rho_w = M_h2o*c_w ;

    /* saturations */
    s_w  = courbe(p_c,el.mat->cb[0]) ;
    s_g  = 1. - s_w ;

    /* porosite */
    phi0 = PHI0(i) ;
    phi  = porosite(phi0,s_w,c_cl_t) ;
    /* phi  = PHI(i) ; */

    /* derivees */
    ds_wsdp_c = dcourbe(p_c,el.mat->cb[0]) ;
    {
      double dc_cl  = K_nacl*1.e-4 ;
      double c_cl1  = c_cl - 0.5*dc_cl,c_cl2 = c_cl + 0.5*dc_cl ;
      double c_na1  = c_cl1,c_na2 = c_cl2 ;
      double c_w1   = (1. - (c_cl1*v_cl + c_na1*v_na))/v_h2o ;
      double c_w2   = (1. - (c_cl2*v_cl + c_na2*v_na))/v_h2o ;
      double lna_w1 = activite(c_cl1,c_na1,c_w1) ;
      double lna_w2 = activite(c_cl2,c_na2,c_w2) ;
      dlna_wsdc_cl = (lna_w2 - lna_w1)/dc_cl ;
    }
    /* ... par rapport a h_r */	
    dp_lsdh_r[i] = RT/v_h2o/h_r ; 
    ds_wsdh_r    = -ds_wsdp_c*dp_lsdh_r[i] ;
    drho_vsdh_r  = M_h2o/RT*p_vs ;
    ds_gsdh_r    = -ds_wsdh_r ;
    {
      double ds_w = 1.e-4 ;
      double s_w1 = s_w - 0.5*ds_w,s_w2 = s_w + 0.5*ds_w ;
      double phi1 = porosite(phi0,s_w1,c_cl_t) ;
      double phi2 = porosite(phi0,s_w2,c_cl_t) ;
      dphisdh_r   = (phi2  - phi1 )/ds_w*ds_wsdh_r ;
      /* dphisdh_r   = 0. ; */
    }
    /* ... par rapport a c_cl_t */
    dc_clsdc_cl_t[i] = (c_cl_t <= K_nacl) ? 1. : 0 ;
    {
      double dc_cl_t = K_nacl*1.e-4 ;
      double c_cl_t1 = c_cl_t - 0.5*dc_cl_t,c_cl_t2 = c_cl_t + 0.5*dc_cl_t ;
      double c_cl1   = concentration_cl_l(c_cl_t1) ;
      double c_cl2   = concentration_cl_l(c_cl_t2) ;
      double phi1    = porosite(phi0,s_w,c_cl_t1) ;
      double phi2    = porosite(phi0,s_w,c_cl_t2) ;
      dc_clsdc_cl_t[i] = (c_cl2 - c_cl1)/dc_cl_t ;
      dphisdc_cl_t     = (phi2  - phi1 )/dc_cl_t ;
      /* dphisdc_cl_t     = 0. ; */
    }
    drho_wsdc_cl_t   = -M_h2o*(v_cl + v_na)/v_h2o*dc_clsdc_cl_t[i] ;
    dp_lsdc_cl_t[i]  = -RT/v_h2o*dlna_wsdc_cl*dc_clsdc_cl_t[i] ;
    ds_wsdc_cl_t     = -ds_wsdp_c*dp_lsdc_cl_t[i] ;
    ds_gsdc_cl_t     = -ds_wsdc_cl_t ;

    /*
      Conservation de Cl : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
    K(i*NEQ+E_Cl,i*NEQ+I_C_Cl) += volume[i]*phi0*(s_w + c_cl_t*ds_wsdc_cl_t) ;
    K(i*NEQ+E_Cl,i*NEQ+I_H_r)  += volume[i]*phi0*(c_cl_t*ds_wsdh_r) ;

    /*
      Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
    */
    K(i*NEQ+E_eau,i*NEQ+I_H_r)  += volume[i]*(phi*(drho_vsdh_r*s_g 
						   + rho_v*ds_gsdh_r
						   + rho_w*ds_wsdh_r)
					      + dphisdh_r*(rho_w*s_w 
							   + rho_v*s_g)) ;
    K(i*NEQ+E_eau,i*NEQ+I_C_Cl) += volume[i]*(dphisdc_cl_t*(rho_v*s_g
							    + rho_w*s_w)
					      + phi*drho_wsdc_cl_t*s_w
					      + phi*rho_w*ds_wsdc_cl_t
					      + phi*rho_v*ds_gsdc_cl_t) ;
  }

  /* termes d'ecoulement */
  {
    double tr  = dt*surf/dx ;
    double S_Dcz2  = KF_Cl*z_cl*z_cl + KF_Na*z_na*z_na ; 
    dpsisdc_cl = -(z_cl*KF_Cl + z_na*KF_Na)/S_Dcz2 ;
    
    trd_cl = tr*KD_Cl ;
    trf_cl = tr*KF_Cl ;
    tre_cl = tr*KF_Cl*z_cl ;

    trd_w  = tr*KD_W ;
    trf_v  = tr*KF_V ; 
  } 

  /*
    Conservation de Cl : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_cl*dp_lsdc_cl_t[i] + (trf_cl + tre_cl*dpsisdc_cl)*dc_clsdc_cl_t[i] ;
  }
  K(E_Cl,I_C_Cl)            += + c[0] ;
  K(E_Cl,I_C_Cl+NEQ)        += - c[1] ;
  K(E_Cl+NEQ,I_C_Cl)        += - c[0] ;
  K(E_Cl+NEQ,I_C_Cl+NEQ)    += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_cl*dp_lsdh_r[i] ;
  }
  K(E_Cl,I_H_r)             += + c[0] ;
  K(E_Cl,I_H_r+NEQ)         += - c[1] ;
  K(E_Cl+NEQ,I_H_r)         += - c[0] ;
  K(E_Cl+NEQ,I_H_r+NEQ)     += + c[1] ;
  
  /*
    Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_w*dp_lsdc_cl_t[i] ;
  }
  K(E_eau,I_C_Cl)           += + c[0] ;
  K(E_eau,I_C_Cl+NEQ)       += - c[1] ;
  K(E_eau+NEQ,I_C_Cl)       += - c[0] ;
  K(E_eau+NEQ,I_C_Cl+NEQ)   += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_w*dp_lsdh_r[i] + trf_v ;
  }
  K(E_eau,I_H_r)            += + c[0] ;
  K(E_eau,I_H_r+NEQ)        += - c[1] ;
  K(E_eau+NEQ,I_H_r)        += - c[0] ;
  K(E_eau+NEQ,I_H_r+NEQ)    += + c[1] ;
 
  return(0) ;

#undef K
}

void rs27(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define N_Cln(n)     (f_n[(n+2)])
#define M_Wn(n)      (f_n[(n)])

#define R(n,i)    (r[(n)*NEQ+(i)])

  double dx ,xm ;
  double volume[2],surf ;
  int    i ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0. ;

  if(el.dim < dim) return ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = 2*M_PI*xm ; else surf = 1. ;

  /*
    Conservation de Cl  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ;

  /*
    Conservation de H2O : (n_w1 - n_wn) + dt * div(w_V) = 0
  */
  R(0,E_eau) -= volume[0]*(M_W(0) - M_Wn(0)) + dt*surf*W_W ;
  R(1,E_eau) -= volume[1]*(M_W(1) - M_Wn(1)) - dt*surf*W_W ;

#undef N_Cln
#undef M_Wn

#undef R
}

int so27(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double* v0 = Element_GetConstantTerm(&el) ;
  int    i,j,nso ;
  double phi,c_w,lna_w;
  double s_w,s_g,p_l,p_c,c_cl,c_na,c_cl_t,h_r ;
  double w_wl,w_v,w_cl ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */

  phi0     = el.mat->pr[pm("porosite")] ;
 
  /* initialisation */
  nso = 12 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = 0. ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* concentrations */
  c_cl_t = param(u,h_s,el.nn,I_C_Cl) ;
  c_cl   = concentration_cl_l(c_cl_t) ;
  c_na   = c_cl ;
  c_w    = (1. - (c_cl*v_cl + c_na*v_na))/v_h2o ;

  /* activite de l'eau */ 
  lna_w = activite(c_cl,c_na,c_w) ;
  
  /* pression */
  h_r   = param(u,h_s,el.nn,I_H_r) ;
  p_l   = RT/v_h2o*(log(h_r) - lna_w) + p_g ;
  p_c   = p_g - p_l ;
  
  /* saturations */
  s_w   = courbe(p_c,el.mat->cb[0]) ;
  s_g   = 1. - s_w ;	 		

  /* porosite */
  phi0 = (s[0] < (x[0][0] + x[1][0])*0.5) ? PHI0(0) : PHI0(1) ;
  phi  = porosite(phi0,s_w,c_cl_t) ;

  {
    double dx = x[1][0] - x[0][0] ;
    /* Gradients */
    double grd_cl  = (C_Cl_l(1) - C_Cl_l(0))/dx ;
    double grd_na  = grd_cl ; 
    double grd_p_l = (P_l(1)  - P_l(0) )/dx ;
    double grd_h_r = (H_r(1)  - H_r(0) )/dx ; 
    double S_Dcz2  = KF_Cl*z_cl*z_cl + KF_Na*z_na*z_na ;  
    double grd_psi = -(z_cl*KF_Cl*grd_cl + z_na*KF_Na*grd_na)/S_Dcz2 ;

    /* Flux */
    w_wl   = - KD_W*grd_p_l ;
    w_v    = - KF_V*grd_h_r ;
    w_cl   = - KD_Cl*grd_p_l - KF_Cl*grd_cl - KF_Cl*z_cl*grd_psi ;
    strcpy(r[11].text,"Grad Psi") ; r[11].n = 1 ;
    r[11].v[0] = grd_psi ;
  }

  /* quantites exploitees */
  strcpy(r[0].text,"h_r") ; r[0].n = 1 ;
  r[0].v[0] = h_r ;
  strcpy(r[1].text,"Log(activite)") ; r[1].n = 1 ;
  r[1].v[0] = lna_w ;
  strcpy(r[2].text,"Saturation") ; r[2].n = 1 ;
  r[2].v[0] = s_w ;
  strcpy(r[3].text,"pression liquide") ; r[3].n = 1 ;
  r[3].v[0] = p_l ;
  strcpy(r[4].text,"flux eau liquide") ; r[4].n = 1 ;
  r[4].v[0] = w_wl ;
  strcpy(r[5].text,"flux vapeur") ; r[5].n = 1 ;
  r[5].v[0] = w_v ;
  strcpy(r[6].text,"Cl libre") ; r[6].n = 1 ;
  r[6].v[0] = c_cl ;
  strcpy(r[7].text,"NaCl solide") ; r[7].n = 1 ;
  r[7].v[0] = s_w*(phi0*c_cl_t - phi*c_cl) ;
  strcpy(r[8].text,"Cl total") ; r[8].n = 1 ;
  r[8].v[0] = phi0*s_w*c_cl_t ;
  strcpy(r[9].text,"porosite") ; r[9].n = 1 ;
  r[9].v[0] = phi ;
  strcpy(r[10].text,"Flux Cl") ; r[10].n = 1 ;
  r[10].v[0] = w_cl ;
  return(nso) ;
}

void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double dx = x[1][0] - x[0][0] ;
  /* Gradients */
  double grd_cl  = (C_Cl_l(1) - C_Cl_l(0))/dx ;
  double grd_p_l = (P_l(1)    - P_l(0)   )/dx ;
  double grd_h_r = (H_r(1)    - H_r(0)   )/dx ;
  double grd_na  = grd_cl ; 
  double S_Dcz2  = KF_Cl*z_cl*z_cl + KF_Na*z_na*z_na ;  
  double grd_psi = -(z_cl*KF_Cl*grd_cl + z_na*KF_Na*grd_na)/S_Dcz2 ;
 
  /* Flux */
  W_W    = - KD_W*grd_p_l  - KF_V*grd_h_r ;
  W_Cl   = - KD_Cl*grd_p_l - KF_Cl*grd_cl - KF_Cl*z_cl*grd_psi ;
}


double activite(double c_cl,double c_na,double c_w)
/* L'activite chimique de l'eau d'une solution de NaCl */
{
  double lna_w ;
  double m_cl,m_na,m_T ;
  double I,A,epsi ;

  double T = TEMPERATURE ;
  double T_0 = 273.15 ;
  double b0 = sqrt(M_h2o),S0 = pow(M_h2o,1.29) ; /* references */
  double b_na = 4.352/b0,b_cl = 1.827/b0 ; /* donnees intrinseques */
  double S_na = 26.448/S0,S_cl = 19.245/S0 ;

  /* molarites */
  if(c_cl < 0.) c_cl = 0. ;
  if(c_na < 0.) c_na = 0. ;
  
  epsi = 0.0007*(T - T_0)*(T - T_0) - 0.3918*(T - T_0) + 87.663 ;
  A = 1398779.816/pow(epsi*T,1.5)/b0 ;
  
  /* molalites*M_h2o (en moles/mole) */
  m_cl = c_cl/c_w ;
  m_na = c_na/c_w ;

  /* la force ionique */
  I = 0.5*(z_cl*z_cl*m_cl + z_na*z_na*m_na) ;
  
  if (I > 0.) {
    m_T =  m_cl + m_na ;

    lna_w = m_cl*lna_i(T,I,z_cl,b_cl,S_cl,A)
          + m_na*lna_i(T,I,z_na,b_na,S_na,A) ;

    /* selon Lin & Lee */
    /*
    lna_cl = lng_LinLee(T,I,z_cl,b_cl,S_cl,A) ;
    lna_na = lng_LinLee(T,I,z_na,b_na,S_na,A) ;
    */

    /* selon TQN */
    /*
    lna_cl = lng_TQN(T,I,z_cl,b_cl,S_cl,A,lna_w,m_T) ;
    lna_na = lng_TQN(T,I,z_na,b_na,S_na,A,lna_w,m_T) ;
    */

  } else {
    /*
    lna_cl = 0. ;
    lna_na = 0. ;
    */
    lna_w  = 0. ;
  }

  return(lna_w) ;
}

double tortuosite_l(double phi)
{
  double phi0 = 0.18 ;
  double dphi = (phi > phi0) ? phi - phi0 : 0. ;
  return(phi*(0.001 + 0.07*phi*phi + 1.8*dphi*dphi)) ;
}

double concentration_cl_l(double c_cl_t)
{
  return((c_cl_t < K_nacl) ? c_cl_t : K_nacl) ;
}

double porosite(double phi0,double s_w,double c_cl_t)
{
  /* r  = 1./(1. + s_w*(c_cl_t - K_nacl)*v_nacl) ; */
  double r    = (1. - s_w*c_cl_t*v_nacl)/(1. - s_w*K_nacl*v_nacl) ;
  double phi  = (c_cl_t < K_nacl) ? phi0 : phi0*r ;
  /* if(phi < 1.e-8*phi0) phi = 1.e-8*phi0 ; */
  if(phi < 0.) phi = 0. ;
  return(phi) ;
}

double lng_TQN(double T,double I,double z,double b,double S,double A,double lna_w,double m_t)
/* Le log du coefficient d'activite d'un ion (T.Q Nguyen) :
   lng_i = dGamma/dm_i = (dGamma/dm_i)_I - 0.5*z_i*z_i*(lna_w + m_t)/I 
   lna_w = - m_t - sum_i ( m_i*lng_i ) + Gamma */
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*2*log(1 + b*II)/b + S*pow(I,alpha)/(1+alpha)/T - 0.5*(lna_w + m_t)/I ;
  
  return(lng*z*z) ;
}

double lna_i(double T,double I,double z,double b,double S,double A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29,a1 = alpha/(1+alpha),II = sqrt(I) ;
  double lna ;
  
  lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/T ;
  
  return(-1 + lna*z*z) ;
}
