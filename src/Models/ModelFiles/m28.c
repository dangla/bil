#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "CommonModel.h"

#define MODELINDEX  28
#define TITLE "Cristallisation de sels (Jan. 2009)"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ 	(2)
#define NVI     (12)
#define NVE     (9)

#define E_eau	(0)
#define E_salt	(1)

#define I_H_r   (0)
#define I_C_s   (1)

#define H_r(n)      (u[(n)][I_H_r])
#define C_s(n)      (u[(n)][I_C_s])

#define M_W(n)      (f[(n)])
#define N_S(n)      (f[(n+2)])

#define W_W         (f[4])
#define W_S         (f[5])

#define P_l(n)      (f[(8+n)])
#define P_c(n)      (f[(10+n)])

#define KD_W        (va[(0)])
#define KF_V        (va[(1)])

#define KD_A        (va[(2)])
#define KD_C        (va[(3)])
#define KF_A        (va[(4)])
#define KF_C        (va[(5)])

#define KD_S        (va[(6)])
#define KF_S        (va[(7)])
#define KM_S        (va[(8)])

/* valences */
#define z_cl      (-1.)
#define z_na      (1.)
#define z_so4     (-2.)

/* coefficients stoechiometriques */
#define nu_h2o_nacl   (0)
#define nu_na_nacl    (1)
#define nu_cl_nacl    (1)
#define nu_h2o_na2so4 (10)
#define nu_na_na2so4  (2)
#define nu_so4_na2so4 (1)

/* volumes molaires partiels liquides (m3/mole) */
#define v_na      (1.87e-5)
#define v_cl      (2.52e-6)
#define v_so4     (-8.334e-6)
#define V_H2O     (1.80e-5)

/* volumes molaires solides (m3/mole) */
#define v_nacl    (24.5e-6)
#define v_na2so4  (220.e-6)

/* coefficients de diffusion moleculaire (m2/s) */
#define do_cl     (2.032e-9)
#define do_na     (1.334e-9)
#define do_so4    (1.065e-9)

#define do_va     (2.42e-5)

/* constante d'equilibre */
#define K_nacl    (6.e3)      /* Solubilite de NaCl (moles/m3) */
#define K_na2so4  (1.24166e3) /* Solubilite de Na2SO4.10H2O (moles/m3) */

/* tensions superficielles (N/m) */
#define Gamma_cl_nacl   (0.1)
#define Gamma_cl_na2so4 (0.1)

#define GAMMA_LG        (0.07)


/* constantes physiques */
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


/* Types de sel */
#define NaCl    (0)
#define Na2SO4  (1)
#define KNO3    (2)

/* Choix du sel */
#define SALT     Na2SO4

#if SALT == NaCl
#define Z_A       z_cl
#define Z_C       z_na
#define NU_A      nu_cl_nacl
#define NU_C      nu_na_nacl
#define NU_H2O    nu_h2o_nacl
#define V_A       v_cl
#define V_C       v_na
#define V_SALT    v_nacl
#define D_A       do_cl
#define D_C       do_na
#define K_SALT    K_nacl
#define GAMMA_CL  Gamma_cl_nacl
#elif SALT == Na2SO4
#define Z_A       z_so4
#define Z_C       z_na
#define NU_A      nu_so4_na2so4
#define NU_C      nu_na_na2so4
#define NU_H2O    nu_h2o_na2so4
#define V_A       v_so4
#define V_C       v_na
#define V_SALT    v_na2so4
#define D_A       do_so4
#define D_C       do_na
#define K_SALT    K_na2so4
#define GAMMA_CL  Gamma_cl_na2so4
#else
#error "Type de sel non prevu"
#endif

#define NU_AC     (NU_A + NU_C)
#define V_AC      (NU_A*V_A + NU_C*V_C)
#define V_ACH     (NU_H2O*V_H2O + V_AC)

/* Fonctions */
static int    pm(const char *s) ;
static double lna_i(double,double,double,double,double,double) ;
static double lng_TQN(double,double,double,double,double,double,double,double) ;
static double lng_LinLee(double,double,double,double,double,double) ;

static double tortuosite_l(double) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;

static double saturation_l(double,double,crbe_t) ;
static double saturation_c(double,double,crbe_t) ;

static double dsaturation_ll(double,double,crbe_t) ;
static double dsaturation_lc(double,double,crbe_t) ;
static double dsaturation_cl(double,double,crbe_t) ;
static double dsaturation_cc(double,double,crbe_t) ;

static double dsaturation(double,double,crbe_t,double (*)(double,double,crbe_t),char) ;

static double activite_w(double,double) ;
static double activite_s(double,double) ;
static double activite_w_ideal(double,double) ;
static double activite_s_ideal(double,double) ;

static double (*xactivite_w[])(double,double) = {activite_w,activite_w_ideal} ;
static double (*xactivite_s[])(double,double) = {activite_s,activite_s_ideal} ;

#define SATURATION_L(x,y)   saturation_l(x,y,el.mat->cb[0])
#define SATURATION_C(x,y)   saturation_c(x,y,el.mat->cb[0])

#define DSATURATION_LL(x,y) dsaturation_ll(x,y,el.mat->cb[0])
#define DSATURATION_LC(x,y) dsaturation_lc(x,y,el.mat->cb[0])
#define DSATURATION_CL(x,y) dsaturation_cl(x,y,el.mat->cb[0])
#define DSATURATION_CC(x,y) dsaturation_cc(x,y,el.mat->cb[0])


/*
#define DSATURATION_LL(x,y) (-dsaturation(x,y,el.mat->cb[0],saturation_l,'x'))
#define DSATURATION_LC(x,y) dsaturation(x,y,el.mat->cb[0],saturation_l,'y')
#define DSATURATION_CL(x,y) (-dsaturation(x,y,el.mat->cb[0],saturation_c,'x'))
#define DSATURATION_CC(x,y) dsaturation(x,y,el.mat->cb[0],saturation_c,'y')
*/

/*
#define ACTIVITE_W(a)   courbe(a,el.mat->cb[3])
#define ACTIVITE_S(a)   courbe(a,el.mat->cb[4])
*/
#define ACTIVITE_W(a)     xactivite_w[1](a,TEMPERATURE)
#define ACTIVITE_S(a)     xactivite_s[1](a,TEMPERATURE)

/* Parametres */
static double phi0,r_d,k_int ;
static double d_cl ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1);  
  else if(strcmp(s,"r_d") == 0) return (2);
  else if(strcmp(s,"k_int") == 0) return (3);
  else if(strcmp(s,"lna_w0") == 0) return (4);
  else if(strcmp(s,"lna_s0") == 0) return (5);
  else if(strcmp(s,"courbes") == 0) return (6);
  else return (-1) ;
}

int dm28(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    n_donnees = 7 ;

  if(dim > 1) arret("dm28 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;

  strcpy(mat->eqn[E_eau],"liq") ;
  strcpy(mat->eqn[E_salt],"sel") ;

  strcpy(mat->inc[I_H_r],"h_r") ;
  strcpy(mat->inc[I_C_s],"c_s") ;

  dmat(mat,ficd,pm,n_donnees) ;

#ifdef NOTDEFINED
  { /* on stocke les courbes lna_w et lna_s */
    int    n_points = 1000 ;
    int    n_courbes = 2 ;
    double c_s1 = 1.e-10*K_SALT,c_s2 = 10*K_SALT ;
    int    i ;

    for(i=0;i<n_courbes;i++) {
      int    j = mat->nc + i ;
      mat->cb[j].echelle = 'l' ;
      mat->cb[j].n = n_points ;
      mat->cb[j].a = mat->pr + mat->n ; /* memes abscisses pour les courbes */
      mat->cb[j].f = mat->cb[j].a + 2 + i*n_points ;
    }

    /* on met a jour le nb de proprietes */
    mat->n  += n_courbes*n_points + 2 ;

    for(i=0;i<n_points;i++) {
      int    n_i  = n_points - 1 ;
      double c_s  = c_s1*pow(c_s2/c_s1,((double) i)/n_i) ;
      int    j    = mat->nc ;
      mat->cb[j].f[i]   = xactivite_w[1](c_s,TEMPERATURE) ;
      mat->cb[j+1].f[i] = xactivite_s[1](c_s,TEMPERATURE) ;
    }
    mat->cb[mat->nc].a[0] = c_s1 ;
    mat->cb[mat->nc].a[1] = c_s2 ;

    mat->nc += n_courbes ;
  }
#endif

  {
    double c_s = K_SALT ;
    mat->pr[pm("lna_w0")] = ACTIVITE_W(c_s) ;
    mat->pr[pm("lna_s0")] = ACTIVITE_S(c_s) ;
  }

  return(mat->n) ;
}

int qm28(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 2 equations :\n\
\t 1. Conservation de la masse d\'eau  (h_r)\n\
\t 2. Conservation de la masse de sel (c_s)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.12  # Porosite\n") ;
  fprintf(ficd,"k_int = 1.e-20   # Permeabilite intrinseque (m2)\n") ;
  fprintf(ficd,"D_Cl = 6.25e-12  # Diffusion effective de Cl (m2/s)\n") ;
  fprintf(ficd,"r_d = 1.         # Rapport des tortuosites anions/cations\n") ;
  fprintf(ficd,"courbes = myfile # nom du fichier : p_c S_l k_rl tau_l\n") ;

  return(NEQ) ;
}

void tb28(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ; /* implicite */
  el->n_ve = NVE ; /* explicite */
}

void ch28(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;

}

void in28(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;

  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double c_s    = C_s(i) ;  
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ;  
    double p_v    = h_r*p_vs ;
    /* concentration en eau liquide */
    double c_w    = (1. - V_AC*c_s)/V_H2O ;
    /* activite de l'eau */
    double lna_w0 = el.mat->pr[pm("lna_w0")] ;
    double lna_w  = ACTIVITE_W(c_s) ;
    double dlna_w = lna_w - lna_w0 ;
    /* activite du sel */
    double lna_s0 = el.mat->pr[pm("lna_s0")] ;
    double lna_s  = ACTIVITE_S(c_s) ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
    double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
               + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
    /* masses volumiques */
    double rho_v = M_h2o/RT*p_v ;
    double rho_w = M_h2o*c_w ;
    /* saturations */
    double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
    double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;
    double s_g  = 1. - s_l - s_c ;

    /* contenus molaires en cations et anions */
    N_S(i)   = phi0*NU_AC*(s_l*c_s + s_c/V_SALT) ;
    /* masse d eau */
    M_W(i)   = phi0*(rho_w*s_l + rho_v*s_g) ;

    /* sauvegarde */
    P_l(i)    = p_l ;
    P_c(i)    = p_c ;
  }

  /* Coefficients de transfert */
  {
    ex_t ex28 ;
    ex28(x,u,f,va,el,dim,geom,0.) ; 
  }

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex28(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  int i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi0     = el.mat->pr[pm("porosite")] ;
  d_cl     = el.mat->pr[pm("D_Cl")] ;
  r_d      = el.mat->pr[pm("r_d")] ;
  k_int	   = el.mat->pr[pm("k_int")] ;

  for(i=0;i<NVE;i++) va[i] = 0. ;

  /* Contenus molaires */
  for(i=0;i<2;i++) {
    double c_s    = C_s(i) ;  
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ;
    /* concentrations */
    double c_w    = (1. - V_AC*c_s)/V_H2O ;
    /* activite de l'eau */
    double lna_w0 = el.mat->pr[pm("lna_w0")] ;
    double lna_w  = ACTIVITE_W(c_s) ;
    double dlna_w = lna_w - lna_w0 ;
    /* activite du sel */
    double lna_s0 = el.mat->pr[pm("lna_s0")] ;
    double lna_s  = ACTIVITE_S(c_s) ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
    double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
               + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
    /* saturations */
    double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
    double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;
    double s_g  = 1. - s_l - s_c ;
    /* permeabilite */
    double k_rl    = courbe(s_l,el.mat->cb[1]) ;
    double k_h     = k_int/mu_l*k_rl ;
    /* tortuosites gaz et liquide*/
    double phi = phi0 ;
    double a = 2.67,b = 4.67 ;
    double tau_g   = (s_g > 0.) ? pow(phi,a)*pow(s_g,b) : 0. ;
    double tau_l   = tortuosite_l(phi)*courbe(s_l,el.mat->cb[2]) ;
    /* tortuosites anions et cations */
    double tau_ani = tau_l*d_cl/(tortuosite_l(phi0)*do_cl) ;
    double tau_cat = tau_ani/r_d ;
  
    /* sauvegarde */
    KD_W   += M_h2o*c_w*k_h ;
    KD_A   += NU_A*c_s*k_h ;
    KD_C   += NU_C*c_s*k_h ;
    KD_S   += KD_A + KD_C ;

    KF_V   += tau_g*do_va*M_h2o/RT*p_vs ;
    KF_A   += tau_ani*D_A ;
    KF_C   += tau_cat*D_C ;
    KF_S   += NU_A*KF_A + NU_C*KF_C ;

    KM_S   += Z_A*KF_A  + Z_C*KF_C ;
  }

  /* on prend la moyenne */
  for(i=0;i<NVE;i++) va[i] *= 0.5 ;

  return(0) ;
}

int ct28(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */

  phi0     = el.mat->pr[pm("porosite")] ;
   
  /* Contenus molaires */

  for(i=0;i<2;i++) {
    double c_s    = C_s(i) ; 
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ; 
    double p_v    = h_r*p_vs ;
    /* concentration en eau liquide */
    double c_w    = (1. - V_AC*c_s)/V_H2O ;
    /* activite de l'eau */
    double lna_w0 = el.mat->pr[pm("lna_w0")] ;
    double lna_w  = ACTIVITE_W(c_s) ;
    double dlna_w = lna_w - lna_w0 ;
    /* activite du sel */
    double lna_s0 = el.mat->pr[pm("lna_s0")] ;
    double lna_s  = ACTIVITE_S(c_s) ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
    double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
               + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
    /* masses volumiques */
    double rho_v = M_h2o/RT*p_v ;
    double rho_w = M_h2o*c_w ;
    /* saturations */
    double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
    double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;
    double s_g  = 1. - s_l - s_c ;

    /* contenus molaires en cations et anions */
    N_S(i)   = phi0*NU_AC*(s_l*c_s + s_c/V_SALT) ;
    /* masse d eau */
    M_W(i)   = phi0*(rho_w*s_l + rho_v*s_g) ;

    /* sauvegarde */
    P_l(i)    = p_l ;
    P_c(i)    = p_c ;

    if(h_r <= 0 || c_s < 0. || c_w < 0.) {
      printf("\n") ;
      printf("x       = %e\n",x[i][0]) ;
      printf("h_r     = %e\n",h_r) ;
      printf("c_s     = %e\n",c_s) ;
      printf("c_w     = %e\n",c_w) ;
      printf("s_l     = %e\n",s_l) ;
      printf("s_c     = %e\n",s_c) ;
      return(-1) ;
    }
    /* printf("\nS_L = %e ; S_C = %e\n",s_l,s_c) ; */
    /* printf("\nc_s = %e\n",c_s) ; */
  }
  
  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
  
  return(0) ;
}

int mx28(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])

  double dx,xm ;
  double volume[2],surf ;
  int    i ;

  double Dp_lSDh_r[2] ;
  double Dp_lSDc_s[2] ;
  
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
    double c_s    = C_s(i) ; 
    double h_r    = H_r(i) ;
    double p_vs   = P_VS(TEMPERATURE) ; 
    double p_v    = h_r*p_vs ;
    /* concentration en eau liquide */
    double c_w    = (1. - V_AC*c_s)/V_H2O ;
    /* activite de l'eau */
    double lna_w0 = el.mat->pr[pm("lna_w0")] ;
    double lna_w  = ACTIVITE_W(c_s) ;
    double dlna_w = lna_w - lna_w0 ;
    /* activite du sel */
    double lna_s0 = el.mat->pr[pm("lna_s0")] ;
    double lna_s  = ACTIVITE_S(c_s) ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
    double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
               + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
    /* masses volumiques */
    double rho_v = M_h2o/RT*p_v ;
    double rho_w = M_h2o*c_w ;
    /* saturations */
    double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
    double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;
    double s_g  = 1. - s_l - s_c ;

    /* derivees */
    /* ... par rapport aux pressions */	
    double ds_lsdp_l = DSATURATION_LL(p_g - p_l,p_c - p_g) ;
    double ds_lsdp_c = DSATURATION_LC(p_g - p_l,p_c - p_g) ;
    double ds_csdp_l = DSATURATION_CL(p_g - p_l,p_c - p_g) ;
    double ds_csdp_c = DSATURATION_CC(p_g - p_l,p_c - p_g) ;
    double ds_gsdp_l = - ds_lsdp_l - ds_csdp_l ;
    double ds_gsdp_c = - ds_lsdp_c - ds_csdp_c ;
    /* ... par rapport a h_r */
    double dp_lsdh_r    = RT/V_H2O/h_r ;
    double dp_csdh_r    = V_ACH*dp_lsdh_r/V_SALT ;
    double ds_lsdh_r    = ds_lsdp_l*dp_lsdh_r + ds_lsdp_c*dp_csdh_r ;
    double ds_csdh_r    = ds_csdp_l*dp_lsdh_r + ds_csdp_c*dp_csdh_r ;
    double ds_gsdh_r    = ds_gsdp_l*dp_lsdh_r + ds_gsdp_c*dp_csdh_r ;
    double drho_vsdh_r  = M_h2o/RT*p_vs ;
    double dn_wsdh_r    = phi0*(drho_vsdh_r*s_g + rho_v*ds_gsdh_r	\
				+ rho_w*ds_lsdh_r) ;
    double dn_ssdh_r    = phi0*(ds_lsdh_r*c_s + ds_csdh_r/V_SALT) ;
    /* ... par rapport a c_s */
    double dc_s        = K_nacl*1.e-4 ;
    double c_s2        = c_s + dc_s ;
    double lna_w2      = ACTIVITE_W(c_s2) ;
    double lna_s2      = ACTIVITE_S(c_s2) ;
    double dlna_wsdc_s = (lna_w2 - lna_w)/dc_s ;
    double dlna_ssdc_s = (lna_s2 - lna_s)/dc_s ;
    double dp_lsdc_s   = -RT/V_H2O*dlna_wsdc_s ;
    double dp_csdc_s   = V_ACH*dp_lsdc_s/V_SALT				\
                       + RT*(NU_H2O*dlna_wsdc_s + dlna_ssdc_s)/V_SALT ;
    double ds_lsdc_s   = ds_lsdp_l*dp_lsdc_s + ds_lsdp_c*dp_csdc_s ;
    double ds_csdc_s   = ds_csdp_l*dp_lsdc_s + ds_csdp_c*dp_csdc_s ;
    double ds_gsdc_s   = ds_gsdp_l*dp_lsdc_s + ds_gsdp_c*dp_csdc_s ;
    double drho_wsdc_s = -M_h2o*V_AC/V_H2O ;
    double dn_wsdc_s   = phi0*(drho_wsdc_s*s_l + rho_w*ds_lsdc_s	\
			       + rho_v*ds_gsdc_s) ;
    double dn_ssdc_s   = phi0*(s_l + ds_lsdc_s*c_s + ds_csdc_s/V_SALT) ;
    /*
      Conservation du sel : (n_s1 - n_sn) + dt * div(w_s) = 0
    */
    K(i*NEQ+E_salt,i*NEQ+I_C_s)  += volume[i]*NU_AC*dn_ssdc_s ;
    K(i*NEQ+E_salt,i*NEQ+I_H_r)  += volume[i]*NU_AC*dn_ssdh_r ;

    /*
      Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
    */
    K(i*NEQ+E_eau,i*NEQ+I_H_r)  += volume[i]*dn_wsdh_r ;
    K(i*NEQ+E_eau,i*NEQ+I_C_s)  += volume[i]*dn_wsdc_s ;

    /* sauvegardes pour les termes de transport */
    Dp_lSDh_r[i] = dp_lsdh_r ;
    Dp_lSDc_s[i] = dp_lsdc_s ;
  }

  /* termes d'ecoulement */
  {
    double tr  = dt*surf/dx ;
    double dz2 = KF_A*Z_A*Z_A + KF_C*Z_C*Z_C ;
    double dpsisdc_s = -(Z_A*KF_A*NU_A + Z_C*KF_C*NU_C)/dz2 ;
    
    double trd_s = tr*KD_S ;
    double trf_s = tr*KF_S ;
    double tre_s = tr*KM_S ;

    double trd_w  = tr*KD_W ;
    double trf_v  = tr*KF_V ;

    double c[2] ;

  /*
    Conservation du sel : (n_s1 - n_sn) + dt * div(w_s) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_s*Dp_lSDc_s[i] + (trf_s + tre_s*dpsisdc_s) ;
  }
  K(E_salt,I_C_s)            += + c[0] ;
  K(E_salt,I_C_s+NEQ)        += - c[1] ;
  K(E_salt+NEQ,I_C_s)        += - c[0] ;
  K(E_salt+NEQ,I_C_s+NEQ)    += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_s*Dp_lSDh_r[i] ;
  }
  K(E_salt,I_H_r)             += + c[0] ;
  K(E_salt,I_H_r+NEQ)         += - c[1] ;
  K(E_salt+NEQ,I_H_r)         += - c[0] ;
  K(E_salt+NEQ,I_H_r+NEQ)     += + c[1] ;
  
  /*
    Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_w*Dp_lSDc_s[i] ;
  }
  K(E_eau,I_C_s)           += + c[0] ;
  K(E_eau,I_C_s+NEQ)       += - c[1] ;
  K(E_eau+NEQ,I_C_s)       += - c[0] ;
  K(E_eau+NEQ,I_C_s+NEQ)   += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_w*Dp_lSDh_r[i] + trf_v ;
  }
  K(E_eau,I_H_r)            += + c[0] ;
  K(E_eau,I_H_r+NEQ)        += - c[1] ;
  K(E_eau+NEQ,I_H_r)        += - c[0] ;
  K(E_eau+NEQ,I_H_r+NEQ)    += + c[1] ;
  }
 
  return(0) ;

#undef K
}

void rs28(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define M_Wn(n)      (f_n[(n)])
#define N_Sn(n)      (f_n[(n+2)])

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
    Conservation du sel  : (n_s1 - n_sn) + dt * div(w_s) = 0
  */
  R(0,E_salt) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S ;
  R(1,E_salt) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S ;

  /*
    Conservation de H2O : (n_w1 - n_wn) + dt * div(w_V) = 0
  */
  R(0,E_eau) -= volume[0]*(M_W(0) - M_Wn(0)) + dt*surf*W_W ;
  R(1,E_eau) -= volume[1]*(M_W(1) - M_Wn(1)) - dt*surf*W_W ;

#undef N_Sn
#undef M_Wn

#undef R
}

int so28(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;

  /* if(el.dim < dim) return(0) ; */
 
  /* initialisation */
  nso = 15 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = 0. ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /*
    Donnees
  */

  phi0     = el.mat->pr[pm("porosite")] ;

  {
    /* concentrations */
    double c_s    = param(u,h_s,el.nn,I_C_s) ;
    double h_r    = param(u,h_s,el.nn,I_H_r) ;
        /* activite de l'eau */
    double lna_w0 = el.mat->pr[pm("lna_w0")] ;
    double lna_w  = ACTIVITE_W(c_s) ;
    double dlna_w = lna_w - lna_w0 ;
    /* activite du sel */
    double lna_s0 = el.mat->pr[pm("lna_s0")] ;
    double lna_s  = ACTIVITE_S(c_s) ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
    double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
               + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
    /* saturations */
    double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
    double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;		
    
    double w_wl = 0,w_v = 0,w_ani = 0,w_cat = 0 ;

    if(el.dim == dim) {
      double dx = x[1][0] - x[0][0] ;
      /* Gradients */
      double grd_c_s   = (C_s(1) - C_s(0))/dx ;
      double grd_p_l   = (P_l(1) - P_l(0))/dx ;
      double grd_h_r   = (H_r(1) - H_r(0))/dx ;
      double grd_c_ani = NU_A*grd_c_s ;
      double grd_c_cat = NU_C*grd_c_s ;
      double dz2       = KF_A*Z_A*Z_A + KF_C*Z_C*Z_C ;  
      double grd_psi   = -(Z_A*KF_A*grd_c_ani + Z_C*KF_C*grd_c_cat)/dz2 ;
      
      /* Flux */
      w_wl   = - KD_W*grd_p_l ;
      w_v    = - KF_V*grd_h_r ;
      w_ani =  - KD_A*grd_p_l - KF_A*grd_c_ani - KF_A*Z_A*grd_psi ;
      w_cat =  - KD_C*grd_p_l - KF_C*grd_c_cat - KF_C*Z_C*grd_psi ;
    }
    
    /* quantites exploitees */
    i = 0 ;
    strcpy(r[i].text,"h_r") ; r[i].n = 1 ;
    r[i++].v[0] = h_r ;
    strcpy(r[i].text,"Concentration sel libre") ; r[i].n = 1 ;
    r[i++].v[0] = c_s ;
    strcpy(r[i].text,"Saturation_liquide") ; r[i].n = 1 ;
    r[i++].v[0] = s_l ;
    strcpy(r[i].text,"Saturation en sel solide") ; r[i].n = 1 ;
    r[i++].v[0] = s_c ;
    strcpy(r[i].text,"pression liquide") ; r[i].n = 1 ;
    r[i++].v[0] = p_l ;
    strcpy(r[i].text,"pression cristal") ; r[i].n = 1 ;
    r[i++].v[0] = p_c ;
    strcpy(r[i].text,"Log(a_w)") ; r[i].n = 1 ;
    r[i++].v[0] = lna_w ;
    strcpy(r[i].text,"Log(a_s)") ; r[i].n = 1 ;
    r[i++].v[0] = dlna_s ;
    strcpy(r[i].text,"Log(a_w) ideal") ; r[i].n = 1 ;
    r[i++].v[0] = activite_w_ideal(c_s,TEMPERATURE) ;
    strcpy(r[i].text,"Log(a_s) ideal") ; r[i].n = 1 ;
    {
    double c_s0 = K_SALT ;
    double lna_s_ideal  = activite_s_ideal(c_s,TEMPERATURE) ;
    double lna_s_ideal0 = activite_s_ideal(c_s0,TEMPERATURE) ;
    r[i++].v[0] = lna_s_ideal - lna_s_ideal0 ;
    }
    strcpy(r[i].text,"Sel total") ; r[i].n = 1 ;
    r[i++].v[0] = phi0*NU_AC*(s_l*c_s + s_c/V_SALT)  ;
    strcpy(r[i].text,"flux eau liquide") ; r[i].n = 1 ;
    r[i++].v[0] = w_wl ;
    strcpy(r[i].text,"flux vapeur") ; r[i].n = 1 ;
    r[i++].v[0] = w_v ;
    strcpy(r[i].text,"Flux anions") ; r[i].n = 1 ;
    r[i++].v[0] = w_ani ;
    strcpy(r[i].text,"Flux cations") ; r[i].n = 1 ;
    r[i++].v[0] = w_cat ;
  }
  return(nso) ;
}

void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double dx = x[1][0] - x[0][0] ;
  /* Gradients */
  double grd_c_s   = (C_s(1) - C_s(0))/dx ;
  double grd_p_l   = (P_l(1) - P_l(0))/dx ;
  double grd_h_r   = (H_r(1) - H_r(0))/dx ;
  double grd_c_ani = NU_A*grd_c_s ;
  double grd_c_cat = NU_C*grd_c_s ;
  double dz2       = KF_A*Z_A*Z_A + KF_C*Z_C*Z_C ;  
  double grd_psi   = -(Z_A*KF_A*grd_c_ani + Z_C*KF_C*grd_c_cat)/dz2 ;

  /*
  double w_ani =  - KD_A*grd_p_l - KF_A*grd_c_ani - KF_A*Z_A*grd_psi ;
  double w_cat =  - KD_C*grd_p_l - KF_C*grd_c_cat - KF_C*Z_C*grd_psi ;
  */
 
  /* Flux */
  W_W    = - KD_W*grd_p_l  - KF_V*grd_h_r ;
  W_S    = - KD_S*grd_p_l  - KF_S*grd_c_s - KM_S*grd_psi ;
}


double tortuosite_l(double phi)
{
  double phi0 = 0.18 ;
  double dphi = (phi > phi0) ? phi - phi0 : 0. ;
  return(phi*(0.001 + 0.07*phi*phi + 1.8*dphi*dphi)) ;
}

double saturation_l(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double s_l ;

  if(p_gl > p_cg*y) {
    s_l = courbe(p_gl,cb) ;
  } else {
    double p_cl = p_cg + p_gl ;
    s_l = courbe(p_cl*x,cb) ;
  }

  return(s_l) ;
}

double saturation_c(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double s_c ;

  if(p_gl > p_cg*y) {
    s_c = 1. - courbe(p_cg*y,cb) ;
  } else {
    double p_cl = p_cg + p_gl ;
    s_c = 1. - courbe(p_cl*x,cb) ;
  }

  return(s_c) ;
}


double dsaturation_ll(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_l ;

  if(p_gl > p_cg*y) {
    ds_l = - dcourbe(p_gl,cb) ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_l = - x*dcourbe(p_cl*x,cb) ;
  }

  return(ds_l) ;
}


double dsaturation_lc(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_l ;

  if(p_gl > p_cg*y) {
    ds_l = 0. ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_l = x*dcourbe(p_cl*x,cb) ;
  }

  return(ds_l) ;
}

double dsaturation_cc(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_c ;

  if(p_gl > p_cg*y) {
    ds_c = - y*dcourbe(p_cg*y,cb) ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_c = - x*dcourbe(p_cl*x,cb) ;
  }

  return(ds_c) ;
}

double dsaturation_cl(double p_gl,double p_cg,crbe_t cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_c ;

  if(p_gl > p_cg*y) {
    ds_c = 0. ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_c = x*dcourbe(p_cl*x,cb) ;
  }

  return(ds_c) ;
}

double dsaturation(double x,double y,crbe_t cb,double (*saturation)(double,double,crbe_t),char v)
{
  int    n_i = cb.n - 1 ;
  double a1 = cb.a[0],a2 = cb.a[1] ;
  double da ;
  double s1,s2,ds ;

  if(cb.echelle == 'n') {
    da = (a2 - a1)/n_i ;
  } else if(cb.echelle == 'l') {
    double loga1 = log10(a1),loga2 = log10(a2) ;
    double dloga = (loga2 - loga1)/n_i ;
    double a = (v == 'x') ? x : y ;
    double ada = a*pow(10.,dloga) ;
    da = ada - a ;
  } else arret("dsaturation_1 : option non prevue") ;

  if(v == 'x') {
    s1 = (*saturation)(x - da,y,cb) ;
    s2 = (*saturation)(x + da,y,cb) ;
  } else {
    s1 = (*saturation)(x,y - da,cb) ;
    s2 = (*saturation)(x,y + da,cb) ;
  }

  ds = (s2 - s1)*0.5/da ;

  return(ds) ;
}

double activite_w(double c_s,double Ta)
/* L'activite chimique de l'eau d'une solution saline */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  /* molalites (moles/mole) */
  double m_s   = (c_s > 0.) ? c_s/c_w : 0. ;
  double m_ani = NU_A*m_s ;
  double m_cat = NU_C*m_s ;
  /* force ionique */
  double I    = 0.5*(Z_A*Z_A*m_ani + Z_C*Z_C*m_cat) ;
  /* references */
  double T0   = 273.15 ;
  double b0   = sqrt(M_h2o),S0   = pow(M_h2o,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(Ta - T0)*(Ta - T0) - 0.3918*(Ta - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*Ta,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    break ;
  } 
  case(Na2SO4) : {
    b_cat = b_na_na2so4 ;
    b_ani = b_so4_na2so4 ;
    S_cat = S_na_na2so4 ;
    S_ani = S_so4_na2so4 ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }

  if (I > 0.) {
    double lna_w ;
    lna_w = m_ani*lna_i(Ta,I,Z_A,b_ani,S_ani,A)
          + m_cat*lna_i(Ta,I,Z_C,b_cat,S_cat,A) ;
    return(lna_w) ;
  } else {
    return(0.) ;
  }
}

double activite_s(double c_s,double Ta)
/* L'activite chimique du sel d'une solution saline */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  /* molalites (moles/mole) */
  double m_s   = (c_s/c_w > DBL_MIN) ? c_s/c_w : DBL_MIN ; 
  double m_ani = NU_A*m_s ;
  double m_cat = NU_C*m_s ;
  /* force ionique */
  double I     = 0.5*(Z_A*Z_A*m_ani + Z_C*Z_C*m_cat) ; 
  /* references */
  double T0   = 273.15 ;
  double b0   = sqrt(M_h2o),S0   = pow(M_h2o,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(Ta - T0)*(Ta - T0) - 0.3918*(Ta - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*Ta,1.5)/b0 ;
  
  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    break ;
  }
  case(Na2SO4) : {
    b_cat = b_na_na2so4 ;
    b_ani = b_so4_na2so4 ;
    S_cat = S_na_na2so4 ;
    S_ani = S_so4_na2so4 ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }

  {
    double m_t     = m_ani + m_cat ;
    double lna_w   = activite_w(c_s,Ta) ;
    double lng_ani = lng_TQN(Ta,I,Z_A,b_ani,S_ani,A,lna_w,m_t) ;
    double lng_cat = lng_TQN(Ta,I,Z_C,b_cat,S_cat,A,lna_w,m_t) ;
    double lna_s   = (NU_A*lng_ani + NU_C*lng_cat) + NU_AC*log(m_s) ;
    return(lna_s) ;
  }
}

double activite_w_ideal(double c_s,double Ta)
/* L'activite chimique de l'eau d'une solution ideale */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  double c_cat = NU_C*c_s ;
  double c_ani = NU_A*c_s ;
  double c_t   = c_w + c_cat + c_ani ;
  double lna_w = log(c_w/c_t) ;
  return(lna_w) ;
}

double activite_s_ideal(double c_s1,double Ta)
/* L'activite chimique du sel d'une solution ideale */
{
  double c_s     = (c_s1 > DBL_MIN/V_H2O) ? c_s1 : DBL_MIN/V_H2O ; 
  double c_w     = (1. - V_AC*c_s)/V_H2O ;
  double c_cat   = NU_C*c_s ;
  double c_ani   = NU_A*c_s ;
  double c_t     = c_w + c_cat + c_ani ;
  double lna_cat = log(c_cat/c_t) ;
  double lna_ani = log(c_ani/c_t) ;
  double lna_s   = NU_C*lna_cat + NU_A*lna_ani ;
  return(lna_s) ;
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
