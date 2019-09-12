#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  52
#define TITLE "Freezing and thawing of concrete with salt"
#define AUTHORS "Zeng"

#include "OldMethods.h"

/* Macros */
#define NEQ      (3)
#define NVI      (21)
#define NVE      (4)

#define E_h2o    (0)
#define E_salt   (1)
#define E_the    (2)

#define I_p_l    (0)
#define I_c_s    (1)
#define I_tem    (2)


#define P_l(n)   (u[(n)][I_p_l])
#define C_s(n)   (u[(n)][I_c_s])
#define TEM(n)   (u[(n)][I_tem])

#define M_H2O_l(n)  (f[(n)])
#define M_H2O_i(n)  (f[(n+2)])
#define M_SALT(n)   (f[(n+4)])
#define S(n)        (f[(n+6)])
#define S_H2O_l(n)  (f[(n+8)])
#define S_H2O_i(n)  (f[(n+10)])
#define S_SALT(n)   (f[(n+12)])
#define Rho(n)      (f[(n+14)])
#define Eps(n)      (f[(n+16)])
#define W_H2O       (f[(18)])
#define W_SALT      (f[(19)])
#define Q           (f[(20)])

#define P_ln(n)  (u_n[(n)][I_p_l])
#define C_sn(n)  (u_n[(n)][I_c_s])
#define TEM_n(n) (u_n[(n)][I_tem])

#define M_H2O_ln(n)  (f_n[(n)])
#define M_H2O_in(n)  (f_n[(n+2)])
#define M_SALTn(n)   (f_n[(n+4)])
#define S_n(n)       (f_n[(n+6)])
#define S_H2O_ln(n)  (f_n[(n+8)])
#define S_H2O_in(n)  (f_n[(n+10)])
#define S_SALTn(n)   (f_n[(n+12)])

#define KD_H2O      (va[(0)])
#define KD_salt     (va[(1)])
#define KF_salt     (va[(2)])
#define KTH         (va[(3)])


/* valences des ions */
#define z_na       (1.)
#define z_cl       (-1.)
#define z_ca       (2.)

/* volumes molaires liquides (m3/mole) */
#define V_h2o      (18.e-6)
#define V_na       (22.47e-6)
#define V_cl       (-0.35e-6)
#define V_ca       (-18.7e-6)

/* volumes molaires solides (m3/mole) */
#define V_ice      (19.63e-6)
#define V_nacl     (24.5e-6)
#define V_cacl2    (40.e-6)

/* Masses molaires (kg/m3) */
#define M_h2o      (18.e-3)
#define M_na       (23.e-3)
#define M_cl       (35.45e-3)
#define M_ca       (40.1e-3)
#define M_cacl2    (110.99e-3)
#define M_nacl     (58.45e-3)

/* coefficients de diffusion moleculaire (m2/s) */
#define d_ca       (7.92e-10)
#define d_na       (1.33e-9)
#define d_cl       (2.032e-9)
#define d_nacl     (1.e-9)
#define d_cacl2    (1.e-9)

/* viscosite */
#define mu_l       (1.79e-3)  /* Viscosite de l'eau (Pa.s) */

/* grandeurs de reference */
#define p_m        (1.e5)      /* Pression atmospherique (Pa)*/
#define T_m        (273.)      /* Temperature de fusion de la glace (K) */

/* Chaleurs specifiques (J/kg/K) */
#define C_l        (4180.)     /* Chaleur specifique du liquide */
#define C_i        (2000.)     /* Chaleur specifique de la glace */

/* entropie de fusion (J/mol/K) */
#define S_m        (23.54)

/* Conductivites thermiques (W/m/K) */
#define LAM_l      (0.6)
#define LAM_i      (2.2)

/* modules de compression (Pa) */
#define K_l        (1.8e9)
#define K_i        (7.8e9)

/* coefficients de dilatation thermique volumique (1/K) */
/* #define ALPHA_l    (-298.e-6) */
#define ALPHA_l(theta)    (-68.7e-6 + 13.877e-6*(theta))  /* coef secant D'apres Speedy (1987) */
#define ALPHA_i    (155.e-6)

/* constantes physiques */
#define R_g        (8.315)    /* Gaz parfait (J/mol/K) */


/* Type de sel */
#define NaCl   0
#define CaCl2  1

/* Choix du sel */
#define SALT   CaCl2

#if SALT == CaCl2
#define V_salt   V_cacl2
#define M_salt   M_cacl2
#define d_salt   d_cacl2
#elif SALT == NaCl
#define V_salt   V_nacl
#define M_salt   M_nacl
#define d_salt   d_nacl
#else
#define V_salt   V_nacl
#define M_salt   M_nacl
#define d_salt   d_nacl
#error "Type de sel non prevu"
#endif

/* Fonctions */
static int    pm(const char *s) ;
static int    flux(double**,double**,double**,double*,double*,double*,elem_t,int,geom_t) ;
static double activity(double,double) ;
static double tortuosity(double,double) ;
static double lna_i(double,double,double,double,double,double) ;
static double lng_LinLee(double,double,double,double,double,double) ;
static double lng_TQN(double,double,double,double,double,double,double,double) ;

/* Parametres */
/* static double phi,k_int,mu_l,R,lam_s,lam_l,lam_i,C_s,C_l,C_i,S_m,T_m,M_h2o,M_salt,V_h2o,V_salt,rho_h2o_i0,d0_salt,k_s,g_s,K_i,K_l,alpha_i,alpha_l,alpha_s,p_m,rho_h2o_l0 ; */
static double lam_l = LAM_l,lam_i = LAM_i ;
static double rho_h2o_i0 = M_h2o/V_ice,rho_h2o_l0 = M_h2o/V_h2o ;
/* static double alpha_i = ALPHA_i,alpha_l = ALPHA_l ; */
static double alpha_i = ALPHA_i ;
static double v_salt = (V_salt/M_salt - V_h2o/M_h2o) ;
static double phi,k_int,lam_s,C_s,k_s,g_s,alpha_s ;


int pm(const char *s)
{
  if(strcmp(s,"phi") == 0) return (0) ;
  else if(strcmp(s,"k_int") == 0) return (1) ;
  else if(strcmp(s,"C_s") == 0) return (2) ;
  else if(strcmp(s,"lam_s") == 0) return (3) ;
  else if(strcmp(s,"k_s") == 0) return (4) ;
  else if(strcmp(s,"g_s") == 0) return (5) ;
  else if(strcmp(s,"alpha_s") == 0) return (6) ;

  else if(strcmp(s,"rho_h2o_i") == 0) return (7) ;
  else if(strcmp(s,"mu_l") == 0) return (8) ;
  else if(strcmp(s,"T_m") == 0) return (9) ;
  else if(strcmp(s,"M_h2o") == 0) return (10) ;
  else if(strcmp(s,"M_salt") == 0) return (11) ;
  else if(strcmp(s,"V_salt") == 0) return (12) ;
  else if(strcmp(s,"V_h2o") == 0) return (13) ;
  else if(strcmp(s,"R") == 0) return (14) ;
  else if(strcmp(s,"lam_l") == 0) return (15) ;
  else if(strcmp(s,"lam_i") == 0) return (16) ;
  else if(strcmp(s,"C_l") == 0) return (17) ;
  else if(strcmp(s,"C_i") == 0) return (18) ;
  else if(strcmp(s,"S_m") == 0) return (19) ;
  else if(strcmp(s,"D_salt") == 0) return (20) ;
  else if(strcmp(s,"rho_h2o_l") == 0) return (21) ;
  else if(strcmp(s,"K_i") == 0) return (22) ;
  else if(strcmp(s,"K_l") == 0) return (23) ;
  else if(strcmp(s,"alpha_i") == 0) return (24) ;
  else if(strcmp(s,"alpha_l") == 0) return (25) ;
  else if(strcmp(s,"p_m") == 0) return (26) ;
  else if(strcmp(s,"courbes") == 0) return (27) ;
  else return (-1) ;
}


int dm52(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 28 ;

  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_h2o],"water") ;
  strcpy(mat->eqn[E_salt],"salt") ;
  strcpy(mat->eqn[E_the],"the") ;

  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_c_s],"c_s") ;
  strcpy(mat->inc[I_tem],"tem") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}


int qm52(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{ 
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 3 equations :\n\
\t 1. Conservation de la masse d\'eau  (p_l)\n\
\t 2. Conservation de la masse de sel (c_s)\n\
\t 2. Bilan d\'entropie (tem)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.2   # Porosite\n") ;
  fprintf(ficd,"k_int = 5.e-21   # Permeabilite intrinseque (m2)\n") ;
  fprintf(ficd,"C_s = 2e+06      # Chaleur volumique du solide\n") ;
  fprintf(ficd,"lam_s = 1.       # Conductivite thermique du solide\n") ;
  fprintf(ficd,"k_s = 3.18e10    # Module de compression du solide\n") ;
  fprintf(ficd,"g_s = 1.91e10    # Module cisaillement du solide\n") ;
  fprintf(ficd,"alpha_s = 54.e-6 # Dilatation thermique vol. du solide\n") ;
  fprintf(ficd,"courbes = myfile # nom du fichier : p_c S_l k_rl tau_l\n") ;

  return(NEQ) ;
}


void tb52(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch52(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i ;
  
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;

  for(i=0;i<el.nn*NEQ;i++) r[i] = -r[i] ;

  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_the)
  {
    R(0,E_the) /= TEM_n(0) ;
  }
    else if(!strcmp(cg.eqn,el.mat->eqn[E_the]))
  {
    R(0,E_the) /= TEM_n(0) ;
  }
#undef R
}

void in52(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;

  if(el.dim < dim) return ;

  /*Donnees */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  k_s     = el.mat->pr[pm("k_s")] ;
  g_s     = el.mat->pr[pm("g_s")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;
  
  /* masses of h2o, salt and entropies */
  for(i=0;i<2;i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;
    double c_s = C_s(i);
    double lna = activity(c_s,tem) ;
    /*pressure of ice crystallization */
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    /* capirally pressure */
    double p_c = p_i - p_l ;
    /*saturation degree  */
    double s_l = courbe(p_c,el.mat->cb[0]) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ; 
    double rho_h2o_l = rho_l - rho_salt ;

    M_H2O_l(i) = phi*s_l*rho_h2o_l ;
    M_H2O_i(i) = phi*s_i*rho_h2o_i ;
    M_SALT(i)  = phi*s_l*rho_salt ;
    S_H2O_l(i) = C_l*log(tem/T_m) - alpha_l/rho_h2o_i0*(p_l - p_m) ;
    S_H2O_i(i) = C_i*log(tem/T_m) - S_m/M_h2o - alpha_i/rho_h2o_i0*(p_i - p_m)  ;
    S_SALT(i)  = S_H2O_l(i) ;
    S(i)       = C_s*log(tem/T_m)
               + M_H2O_l(i)*S_H2O_l(i)
               + M_H2O_i(i)*S_H2O_i(i)
               + M_SALT(i)*S_SALT(i) ;
  }

  {
    ex_t ex52 ;
    ex52(x,u,f,va,el,dim,geom,0.) ;
  }
  /* flux */
  flux(x,u,u,f,f,va,el,dim,geom) ;
}

int ex52(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  int    i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  k_s     = el.mat->pr[pm("k_s")] ;
  g_s     = el.mat->pr[pm("g_s")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;

  for(i=0;i<NVE;i++) va[i] = 0. ;

  /* coefficient de transfert */
  for(i=0;i<2;i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;  
    double c_s = C_s(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;

    double s_l = courbe(p_c,el.mat->cb[0]) ;
    double s_i = 1. - s_l ;
    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ; 
    double rho_h2o_l = rho_l - rho_salt ;

    double k_l   = k_int/mu_l*courbe(p_c,el.mat->cb[1]) ;
    double tau_l = tortuosity(phi,s_l) ;
    double lam_h2o = s_l*lam_l + s_i*lam_i ;

    KD_H2O  += rho_h2o_l*k_l ;
    KD_salt += rho_salt*k_l ;
    KF_salt += M_salt*phi*s_l*tau_l*d_salt ;
    KTH     += lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
  }

  for(i=0;i<NVE;i++) va[i] *= 0.5 ;

  return(0) ;
}


int ct52(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  k_s     = el.mat->pr[pm("k_s")] ;
  g_s     = el.mat->pr[pm("g_s")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;

  /* masses of h2o, salt and entropies */
  for(i=0;i<2;i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;
    double c_s = C_s(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;

    double s_l = courbe(p_c,el.mat->cb[0]) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt);
    double rho_h2o_l = rho_l - rho_salt ;

    M_H2O_l(i) = phi*s_l*rho_h2o_l ;
    M_H2O_i(i) = phi*s_i*rho_h2o_i ;
    M_SALT(i)  = phi*s_l*rho_salt ;
    S_H2O_l(i) = C_l*log(tem/T_m) - alpha_l/rho_h2o_i0*(p_l - p_m) ;
    S_H2O_i(i) = C_i*log(tem/T_m) - S_m/M_h2o - alpha_i/rho_h2o_i0*(p_i - p_m) ;
    S_SALT(i)  = S_H2O_l(i) ;
    S(i)       = C_s*log(tem/T_m)
                + M_H2O_l(i)*S_H2O_l(i)
                + M_H2O_i(i)*S_H2O_i(i)
                + M_SALT(i)*S_SALT(i) ;
  }

  /* flux */
  flux(x,u,u_n,f,f_n,va,el,dim,geom) ;

  return(0) ;
} 


int mx52(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])

  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  k_s     = el.mat->pr[pm("k_s")] ;
  g_s     = el.mat->pr[pm("g_s")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++)
  {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    double p_l = P_l(i) ;
    double c_s = C_s(i) ;
    double tem = TEM(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;
    
    double s_l = courbe(p_c,el.mat->cb[0]) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. - v_salt*rho_salt + (p_l - p_m)/K_l - alpha_l*(tem - T_m)) ;
    double rho_h2o_l = rho_l - rho_salt ;

    /* pre-derivatives */
    double ds_lsdp_c = dcourbe(p_c,el.mat->cb[0]) ;
    /* derivatives with respect to ... */
    /* ... p_l */
    double dp_isdp_l = V_h2o/V_ice ;
    double dp_csdp_l = dp_isdp_l - 1. ;
    double drho_h2o_lsdp_l = rho_h2o_l0/K_l ;
    double drho_h2o_isdp_l = rho_h2o_i0/K_i*dp_isdp_l ;
    double dm_h2o_lsdp_l = phi*(ds_lsdp_c*dp_csdp_l*rho_h2o_l + s_l*drho_h2o_lsdp_l) ;
    double dm_h2o_isdp_l = phi*(-ds_lsdp_c*dp_csdp_l*rho_h2o_i + s_i*drho_h2o_isdp_l) ;
    double dm_saltsdp_l  = phi*ds_lsdp_c*dp_csdp_l*rho_salt ;
    double ds_h2o_lsdp_l = -alpha_l/rho_h2o_i0 ;
    double ds_h2o_isdp_l = -alpha_i/rho_h2o_i0*dp_isdp_l ;
    double ds_saltsdp_l  = ds_h2o_lsdp_l ;
    /* ... c_s */
    double dc_s = 1.e-5,c_s2 = c_s + dc_s ;
    double dlnasdc_s = (activity(c_s2,tem) - lna)/dc_s ;
    double dp_isdc_s = R_g*tem*dlnasdc_s/V_ice ;
    double dp_csdc_s = dp_isdc_s ;
    double drho_h2o_lsdc_s = -rho_h2o_l0*v_salt*M_salt - M_salt ;
    double drho_h2o_isdc_s = rho_h2o_i0/K_i*dp_isdc_s ;
    double dm_h2o_lsdc_s = phi*(ds_lsdp_c*dp_csdc_s*rho_h2o_l + s_l*drho_h2o_lsdc_s) ;
    double dm_h2o_isdc_s = phi*(-ds_lsdp_c*dp_csdc_s*rho_h2o_i + s_i*drho_h2o_isdc_s) ;
    double dm_saltsdc_s  = phi*(ds_lsdp_c*dp_csdc_s*rho_salt + s_l*M_salt) ;
    double ds_h2o_lsdc_s = 0. ;
    double ds_h2o_isdc_s = -alpha_i/rho_h2o_i0*dp_isdc_s ;
    double ds_saltsdc_s = ds_h2o_lsdc_s ;
    /* ... tem */
    double dtem = 1.,tem2 = tem + dtem ;
    double dlnasdtem = (activity(c_s,tem2) - lna)/dtem ;
    double dp_isdtem = (- S_m + R_g*lna + R_g*tem*dlnasdtem)/V_ice ;
    double dp_csdtem = dp_isdtem ;
    double drho_h2o_lsdtem = -rho_h2o_l0*alpha_l ;
    double drho_h2o_isdtem = rho_h2o_i0/K_i*dp_isdtem - rho_h2o_i0*alpha_i ;
    double dm_h2o_lsdtem = phi*(ds_lsdp_c*dp_csdtem*rho_h2o_l + s_l*drho_h2o_lsdtem) ;
    double dm_h2o_isdtem = phi*(-ds_lsdp_c*dp_csdtem*rho_h2o_i + s_i*drho_h2o_isdtem) ;
    double dm_saltsdtem  = phi*ds_lsdp_c*dp_csdtem*rho_salt ;
    double ds_h2o_lsdtem = C_l/tem ;
    double ds_h2o_isdtem = C_i/tem - alpha_i/rho_h2o_i0*dp_isdtem ;
    double ds_saltsdtem  = ds_h2o_lsdtem ;
    /*
      MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
    */
    K(E_h2o+i*NEQ,I_p_l+i*NEQ) += volume[i]*(dm_h2o_lsdp_l + dm_h2o_isdp_l) ;
    K(E_h2o+i*NEQ,I_c_s+i*NEQ) += volume[i]*(dm_h2o_lsdc_s + dm_h2o_isdc_s) ;
    K(E_h2o+i*NEQ,I_tem+i*NEQ) += volume[i]*(dm_h2o_lsdtem + dm_h2o_isdtem) ;
    /*
      MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
    */
    K(E_salt+i*NEQ,I_p_l+i*NEQ) += volume[i]*dm_saltsdp_l ;
    K(E_salt+i*NEQ,I_c_s+i*NEQ) += volume[i]*dm_saltsdc_s ;
    K(E_salt+i*NEQ,I_tem+i*NEQ) += volume[i]*dm_saltsdtem ;
    /*
      ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
    */
    K(E_the+i*NEQ,I_p_l+i*NEQ) += volume[i]*(dm_h2o_lsdp_l*S_H2O_l(i) + dm_h2o_isdp_l*S_H2O_i(i) + dm_saltsdp_l*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdp_l + M_H2O_i(i)*ds_h2o_isdp_l + M_SALT(i)*ds_saltsdp_l) ;
    K(E_the+i*NEQ,I_c_s+i*NEQ) += volume[i]*(dm_h2o_lsdc_s*S_H2O_l(i) + dm_h2o_isdc_s*S_H2O_i(i) + dm_saltsdc_s*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdc_s + M_H2O_i(i)*ds_h2o_isdc_s + M_SALT(i)*ds_saltsdc_s) ;
    K(E_the+i*NEQ,I_tem+i*NEQ) += volume[i]*(C_s/tem + dm_h2o_lsdtem*S_H2O_l(i) + dm_h2o_isdtem*S_H2O_i(i) + dm_saltsdtem*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdtem + M_H2O_i(i)*ds_h2o_isdtem + M_SALT(i)*ds_saltsdtem) ;
  }
  /*
    termes d'ecoulement
  */
  {
  double trd_h2o  = dt*surf/dx*KD_H2O ;
  double trf_h2o  = dt*surf/dx*(-KF_salt) ;
  double trd_salt = dt*surf/dx*KD_salt ;
  double trf_salt = dt*surf/dx*KF_salt ;
  double trth     = dt*surf/dx*KTH ;

  double tem      = (TEM_n(0)    + TEM_n(1)   )*0.5 ;
  double s_h2o_l  = (S_H2O_ln(0) + S_H2O_ln(1))*0.5 ;
  double s_salt   = (S_SALTn(0)  + S_SALTn(1 ))*0.5 ;
  double c ;
  /*
    MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
  */
  K(E_h2o,I_p_l)          += + trd_h2o ;
  K(E_h2o,I_p_l+NEQ)      += - trd_h2o ;
  K(E_h2o+NEQ,I_p_l)      += - trd_h2o ;
  K(E_h2o+NEQ,I_p_l+NEQ)  += + trd_h2o ;

  K(E_h2o,I_c_s)          += + trf_h2o ;
  K(E_h2o,I_c_s+NEQ)      += - trf_h2o ;
  K(E_h2o+NEQ,I_c_s)      += - trf_h2o ;
  K(E_h2o+NEQ,I_c_s+NEQ)  += + trf_h2o ;

  /*
    MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
  */
  K(E_salt,I_p_l)          += + trd_salt ;
  K(E_salt,I_p_l+NEQ)      += - trd_salt ;
  K(E_salt+NEQ,I_p_l)      += - trd_salt ;
  K(E_salt+NEQ,I_p_l+NEQ)  += + trd_salt ;

  K(E_salt,I_c_s)          += + trf_salt ;
  K(E_salt,I_c_s+NEQ)      += - trf_salt ;
  K(E_salt+NEQ,I_c_s)      += - trf_salt ;
  K(E_salt+NEQ,I_c_s+NEQ)  += + trf_salt ;

  /*
    ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  c = s_h2o_l*trd_h2o + s_salt*trd_salt ;
  K(E_the,I_p_l)           += + c ;
  K(E_the,I_p_l+NEQ)       += - c ;
  K(E_the+NEQ,I_p_l+NEQ)   += + c ;
  K(E_the+NEQ,I_p_l)       += - c ;

  c = s_h2o_l*trf_h2o + s_salt*trf_salt ;
  K(E_the,I_c_s)           += + c ;
  K(E_the,I_c_s+NEQ)       += - c ;
  K(E_the+NEQ,I_c_s+NEQ)   += + c ;
  K(E_the+NEQ,I_c_s)       += - c ;

  c = trth/tem ;
  K(E_the,I_tem)           += + c ;
  K(E_the,I_tem+NEQ)       += - c ;
  K(E_the+NEQ,I_tem+NEQ)   += + c ;
  K(E_the+NEQ,I_tem)       += - c ;
  }

  return(0) ;
#undef K
}

void rs52(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++)
  {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;

  /*
    MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
  */
  R(0,E_h2o) -= volume[0]*(M_H2O_l(0) + M_H2O_i(0) - M_H2O_ln(0) - M_H2O_in(0)) + dt*surf*W_H2O ;
  R(1,E_h2o) -= volume[1]*(M_H2O_l(1) + M_H2O_i(1) - M_H2O_ln(1) - M_H2O_in(1)) - dt*surf*W_H2O ;
  /*
    MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
  */
  R(0,E_salt) -= volume[0]*(M_SALT(0) - M_SALTn(0)) + dt*surf*W_SALT ;
  R(1,E_salt) -= volume[1]*(M_SALT(1) - M_SALTn(1)) - dt*surf*W_SALT ;
  /*
    ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  R(0,E_the) -= volume[0]*(S(0) - S_n(0)) + dt*surf*Q ;
  R(1,E_the) -= volume[1]*(S(1) - S_n(1)) - dt*surf*Q ;

#undef R
}

int so52(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  /* if(el.dim < dim) return(0) ; */

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  lam_s   = el.mat->pr[pm("lam_s")] ;
  C_s     = el.mat->pr[pm("C_s")] ;
  k_s     = el.mat->pr[pm("k_s")] ;
  g_s     = el.mat->pr[pm("g_s")] ;
  alpha_s = el.mat->pr[pm("alpha_s")] ;


  /* initialisation */
  nso = 16 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  {
    /* pressure */
    double p_l  =  param(u,h_s,el.nn,I_p_l) ;
    /* concentration */
    double c_s  =  param(u,h_s,el.nn,I_c_s) ;
    /* temperature */
    double tem  =  param(u,h_s,el.nn,I_tem) ;
    /* activity */
    double lna  =  activity(c_s,tem) ;
    /* ice pressure */
    double p_i1 = V_h2o*(p_l - p_m)/V_ice ;
    double p_i2 = S_m*(T_m - tem)/V_ice ;
    double p_i3 = R_g*tem*lna/V_ice ;
    double p_i  = p_m + p_i1 + p_i2 + p_i3 ;
    /* capillary pressure */
    double p_c  = p_i - p_l ;
    /* liquid saturation */
    double s_l  = courbe(p_c,el.mat->cb[0]) ;
    double s_i  = 1 - s_l ;
    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ;
    
    /* original Qiang !!!
    double K = k_s*(1 - phi*(1 + 3*k_s/(4*g_s))) ;
    double G = g_s*(1 - 5*phi*(4*g_s + 3*k_s)/(8*g_s + 9*k_s)) ;
    double b_i  = phi*s_i*(1 + 3*k_s/(4*g_s)) ;
    double b_l  = phi*s_l*(1 + 3*k_s/(4*g_s)) ;
    double N_ii = 3*phi*s_i/(4*g_s) ;
    double N_ll = 3*phi*s_l/(4*g_s) ; 
    double alpha_phi_l = alpha_s*(b_l - phi*s_l) ;
    double alpha_phi_i = alpha_s*(b_i - phi*s_i) ;
    double Epsi   = (b_l*p_l + b_i*p_i + alpha_s*(tem - T_m))/K ;
    double phi_l  = b_l*Epsi + p_l*N_ll - alpha_phi_l*(tem - T_m) ;
    double phi_i  = b_i*Epsi + p_i*N_ii - alpha_phi_i*(tem - T_m) ;
    double phi_p  = phi + phi_l + phi_i ;
    
    double Eps    = (b_l*(p_l - p_m) + b_i*(p_i - p_m) + alpha_s*K*(tem - T_m))/K ;
    double Eps1   = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double Eps2   = alpha_s*K*(tem - T_m)/K ;
    double Rho    = (phi_p<1) ? (phi_p*s_l*(p_l - p_m) + phi_p*s_i*(p_i - p_m))/(1-phi_p) : 0 ;
    double Eps_T  = alpha_s*K*(tem - T_m)/K ;
    double Eps_P  = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double v = 0.2 ;
    double Eps_xx = Eps*(1+v)/(1-v);
    double Rho_yy = 2*G*Eps_xx ;
    double Rho_zz = Rho_yy ;
    */
    
    /* coef poroelastiques */
    double K = k_s*(1 - phi*(1 + 0.75*k_s/g_s)) ;
    double G = g_s*(1 - 5*phi*(4*g_s + 3*k_s)/(8*g_s + 9*k_s)) ;
    double k_oedo = K + 4*G/3. ; /* Module oedometrique */
    double b_i  = phi*s_i*(1 + 0.75*k_s/g_s) ;
    double b_l  = phi*s_l*(1 + 0.75*k_s/g_s) ;
    double N_ii = 0.75*phi*s_i/g_s ;
    double N_ll = 0.75*phi*s_l/g_s ; 
    double alpha_phi_l = alpha_s*(b_l - phi*s_l) ;
    double alpha_phi_i = alpha_s*(b_i - phi*s_i) ;
    /* conditions Sig_xx = Eps_yy = Eps_zz = 0 */
    double Eps_T  = alpha_s*K*(tem - T_m)/k_oedo ;
    double Eps_P  = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/k_oedo ;
    double Eps_xx = Eps_P + Eps_T ;
    double Sig_yy = - 2*G*Eps_xx ;
    double Sig_zz = Sig_yy ;
    double phi_l  = b_l*Eps_xx + N_ll*(p_l - p_m) - alpha_phi_l*(tem - T_m) ;
    double phi_i  = b_i*Eps_xx + N_ii*(p_i - p_m) - alpha_phi_i*(tem - T_m) ;
    double phi_p  = phi + phi_l + phi_i ;
    
    i = 0 ;
    /* quantites exploitees */
    strcpy(r[i].text,"p_l") ; r[i].n = 1 ;
    r[i++].v[0] = p_l ;
    strcpy(r[i].text,"c_salt") ; r[i].n = 1 ;
    r[i++].v[0] = c_s ;
    strcpy(r[i].text,"Temperature") ; r[i].n = 1 ;
    r[i++].v[0] = tem ;
    strcpy(r[i].text,"Liquid-saturation") ; r[i].n = 1 ;
    r[i++].v[0] = s_l ;
    strcpy(r[i].text,"Flux_H2O") ; r[i].n = 1 ;
    r[i++].v[0] = W_H2O ;
    strcpy(r[i].text,"Flux_salt") ; r[i].n = 1 ;
    r[i++].v[0] = W_SALT ;
    strcpy(r[i].text,"p_i") ; r[i].n = 1 ;
    r[i++].v[0] = p_i ;
    strcpy(r[i].text,"Solution density") ; r[i].n = 1 ;
    r[i++].v[0] = rho_l ;
    strcpy(r[i].text,"Ice density") ; r[i].n = 1 ;
    r[i++].v[0] = rho_h2o_i ;
    strcpy(r[i].text,"strain") ; r[i].n = 1 ;
    r[i++].v[0] = Eps_xx ;
    strcpy(r[i].text,"stress") ; r[i].n = 1 ;
    r[i++].v[0] = Sig_yy ;
    strcpy(r[i].text,"strain_temperature") ; r[i].n = 1 ;
    r[i++].v[0] = Eps_T ;
    strcpy(r[i].text,"strain_pressure") ; r[i].n = 1 ;
    r[i++].v[0] = Eps_P ;
    strcpy(r[i].text,"porosity") ; r[i].n = 1 ;
    r[i++].v[0] = phi_p ;
    strcpy(r[i].text,"deformation of ice pores") ; r[i].n = 1 ;
    r[i++].v[0] = phi_i ;
    strcpy(r[i].text,"deformation of liquid pores") ; r[i].n = 1 ;
    r[i++].v[0] = phi_l ;
  }

  return (nso) ;
}


int flux(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom)
{

  if(el.dim < dim) return(0) ;

  {
    /* gradients */
    double dx      = x[1][0] - x[0][0] ;
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    double grd_c_s = (C_s(1) - C_s(0))/dx ;
    double grd_tem = (TEM(1) - TEM(0))/dx ;
    
    double tem     = (TEM_n(0)    + TEM_n(1)   )*0.5 ;
    double s_h2o_l = (S_H2O_ln(0) + S_H2O_ln(1))*0.5 ;
    double s_salt  = (S_SALTn(0)  + S_SALTn(1) )*0.5 ;
    
    /* flux */
    W_H2O  = - KD_H2O*grd_p_l  + KF_salt*grd_c_s ;
    W_SALT = - KD_salt*grd_p_l - KF_salt*grd_c_s ;
    Q      = - KTH/tem*grd_tem + s_h2o_l*W_H2O + s_salt*W_SALT ;
  }

  return(0) ;
} 

double activity(double c_s,double tem)
/* activity of water */
{
  double T_0  = T_m ;
  double b0   = sqrt(M_h2o),S0 = pow(M_h2o,1.29) ; /* references */
  /* NaCl (d'apres Lin et Lee) */
  double b_na_nacl = 4.352/b0,b_cl_nacl = 1.827/b0 ;
  double S_na_nacl = 26.448/S0,S_cl_nacl = 19.245/S0 ;
  /* CaCl2 (d'apres Lin et Lee) */
  double b_ca_cacl2 = 3.908/b0,b_cl_cacl2 = 2.085/b0 ;
  double S_ca_cacl2 = 18.321/S0,S_cl_cacl2 = 10.745/S0 ;

  double epsi = 0.0007*(tem - T_0)*(tem - T_0) - 0.3918*(tem - T_0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*tem,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;
  double c_ani ;
  double c_cat ;
  double z_ani ;
  double z_cat ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    c_ani = c_s ;
    c_cat = c_s ;
    z_ani = z_cl ;
    z_cat = z_na ;
    break ;
  }
  case(CaCl2) : {
    b_cat = b_ca_cacl2 ;
    b_ani = b_cl_cacl2 ;
    S_cat = S_ca_cacl2 ;
    S_ani = S_cl_cacl2 ;
    c_ani = 2*c_s ;
    c_cat = c_s ;
    z_ani = z_cl ;
    z_cat = z_ca ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }
  
  {
  /* concentrations */
  double c_h2o = (1 - c_s*V_salt)/V_h2o ;
  /* molalites*M_h2o */
  double m_ani  = c_ani/c_h2o ;
  double m_cat  = c_cat/c_h2o ;

  /* ion strenght */
  double I     =  0.5*(z_ani*z_ani*m_ani + z_cat*z_cat*m_cat);
  
  double II_ani   = lna_i(tem,I,z_ani,b_ani,S_ani,A) ;
  double II_cat   = lna_i(tem,I,z_cat,b_cat,S_cat,A) ;

  /* activity of water */
  double lna_h2o = m_ani*II_ani + m_cat*II_cat ;

  return(lna_h2o) ;

  /* linearised term */
  lna_h2o = - (m_ani + m_cat) ;

  return(lna_h2o) ;
  }
}

double tortuosity(double phi,double s_l)
{
  double tau_l_sat = 0.296e-3*exp(9.95*phi)/phi ;
  if(s_l > 0.) return(tau_l_sat/(s_l*(1 + 625*pow(1 - s_l,4)))) ;
  else return(0.) ;
}






double lng_LinLee(double T,double I,double z,double b,double S,double A)
/* Le log du coefficient d'activite d'un ion d'apres Lin & Lee */ 
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*(II/(1 + b*II) + 2*log(1 + b*II)/b) + S*pow(I,alpha)/T ;
  
  return(lng*z*z) ;
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
