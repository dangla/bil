/*

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  70
#define TITLE "Sulfuric acid attack (11/2011)"
#define AUTHORS "Yuan"

#include "OldMethods.h"

/* Macros */
#define NEQ     (3)
#define NVE     (36)
#define NVI     (21)
#define NVE_TR  (28)

#define E_S     (0)
#define E_q     (1)
#define E_Ca    (2)


#define I_H2SO4 (0)
#define I_psi   (1)
#define I_Ca_S  (2)


#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_H2SO4   LOG_RHO

#if (U_H2SO4 == LOG_RHO)
  #define C_H2SO4(n)   (exp(Ln10*u[(n)][I_H2SO4]))
  #define C_H2SO4n(n)  (exp(Ln10*u_n[(n)][I_H2SO4]))
#else
  #define C_H2SO4(n)   (u[(n)][I_H2SO4])
  #define C_H2SO4n(n)  (u_n[(n)][I_H2SO4])
#endif
#define ZN_Ca_S(n)  (u[(n)][I_Ca_S])
#define PSI(n)      (u[(n)][I_psi])


#define ZN_Ca_Sn(n) (u_n[(n)][I_Ca_S])
#define PSIn(n)     (u_n[(n)][I_psi])


#define N_S(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define N_K(n)     (f[(8+n)])
#define W_S        (f[10])
#define W_q        (f[11])
#define W_Ca       (f[12])
#define W_Si       (f[13])
#define W_K        (f[14])
#define N_CH(n)    (f[(15+n)])
#define N_CSH2(n)    (f[(17+n)])
#define N_Si_S(n)  (f[(19+n)])

#define N_Sn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_CHn(n)   (f_n[(15+n)])
#define N_CSH2n(n)   (f_n[(17+n)])
#define N_Si_Sn(n) (f_n[(19+n)])


#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_SO3      (va[(2)])
#define KF_H2SO4    (va[(3)])
#define KF_HSO4     (va[(4)])
#define KF_SO4      (va[(5)])
#define KF_Ca       (va[(6)])
#define KF_CaHSO4   (va[(7)])
#define KF_CaH3SiO4 (va[(8)])
#define KF_H3SiO4   (va[(9)])
#define KF_H4SiO4   (va[(10)])
#define KF_H2SiO4   (va[(11)])
#define KF_CaH2SiO4 (va[(12)])
#define KF_CaSO4aq  (va[(13)])
#define KF_CaOH     (va[(14)])
#define KF_K        (va[(15)])

#define Kpsi_OH       (va[(16)])
#define Kpsi_H        (va[(17)])
#define Kpsi_HSO4     (va[(18)])
#define Kpsi_SO4      (va[(19)])
#define Kpsi_Ca       (va[(20)])
#define Kpsi_CaHSO4   (va[(21)])
#define Kpsi_CaH3SiO4 (va[(22)])
#define Kpsi_H3SiO4   (va[(23)])
#define Kpsi_q        (va[(24)])
#define Kpsi_H2SiO4   (va[(25)])
#define Kpsi_CaOH     (va[(26)])
#define Kpsi_K        (va[(27)])

#define N_CH0(n)      (va[(28+n)])
#define N_CSH20(n)      (va[(30+n)])
#define N_Si_S0(n)    (va[(32+n)])
#define ZP_CH0(n)     (va[(34+n)])

/*
  Solution aqueuse
*/

/* les valences */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_hso4     (-1.)
#define z_so4      (-2.)
#define z_h3sio4   (-1.)
#define z_cahso4   (1.)
#define z_cah3sio4 (1.)
#define z_h2sio4   (-2.)
#define z_caoh     (1.)
#define z_k        (-1.)

/* volumes molaires partiels des ions (dm3/mole) */
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_so3      (32.81e-3)     /* unknown */
#define v_h2so4    (50.e-3)       /* unknown */
#define v_hso4     (24.21.e-3)    /* unknown */
#define v_so4      (22.9e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_caso4aq  (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_koh      (27.44e-3)    /* d'apres Antoine */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (1.22e-7)    /* 1.22e-7 (radius = 1.75e-10 m) */
#define d_h        (9.310e-7)    /* 4.76e-8 (radius = 4.5e-10 m) */
#define d_so3      (1.91e-7)      /* unknown */
#define d_h2so4    (1.5e-7)      
#define d_hso4     (1.385e-7)    /* 1.07e-7 (radius = 2e-10 m) */
#define d_so4      (1.065e-7)    /* 9.53e-8 (radius = 2.25e-10 m) */
#define d_ca       (7.92e-8)
#define d_cahso4   (1.07e-7)    /* (radius = 2e-10 m) */
#define d_sioh4    (xxx)
#define d_h4sio4   (1.07e-7)     /* */
#define d_h3sio4   (1.07e-7)   /* (radius = 2e-10 m) */
#define d_h2sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_cah2sio4 (1.07e-7)    /* a modifier */
#define d_cah3sio4 (1.07e-7)     /*(radius = 2e-10 m) */
#define d_caso4aq  (1.43e-7)    /* (radius = 1.5e-10 m) */
#define d_caoh     (1.07e-7)    /* (radius = 2e-10 m) */
#define d_k        (1.43e-7)   /* (radius = 1.5e-10 m) */
#define d_koh      (1.43e-7)   /* (radius = 1.5e-10 m) */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */
#define K_henry    (1.238)           /* cste de Henry du SO3 / RT */

/*#define K_h2so4    (1.7e-3)           SO3[0] + H2O = H2SO4 1.7e-3*/
#define K_hso4     (1.0e3)          /* H2SO4   = HSO4[-] + H[+] 2.5e-4*/
#define K_so4      (1.0e-2)           /* HSO4[-] = SO4[2-] + H[+] 4.69e-11*/

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O 4.68*/
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+] = H4SiO4 6.45e9*/
#define K_h3sio4   (1.55e-10)        /* H4SiO4    = H3SiO4[-] + H[+] 1.55e-10*/

#define K_cahso4   (1.276e+1)        /* Ca[2+] + HSO4[-]    = CaHSO4[+] 1.276e+1*/
#define K_caso4aq  (1.4e+3)          /* Ca[2+] + SO4[2-]    = CaSO4[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 3.98e+4*/
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] 1.58e+1*/
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-] = CaOH[+] */

/*
  Solides
  CH  = Portlandite
  CC  = Calcite
  CSH = Hydrated Calcium Silicates
  Sam = Amorphous Silica
*/

/* CH */
/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* CH = Ca[2+] + 2OH[-] */
/* volumes molaires solides (dm3/mole) */
#define V_CH       (0.e-3)      /* (33.e-3) */

/* CSH2 */
/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CSH2       (2.5e-5)         /* CSH2 = Ca[2+] + SO4[2-]  3.89e-9*/ 
/* volumes molaires solides (dm3/mole) */
#define V_CSH2       (0.e-3)      /* (75.e-3) */

/* C-S-H */
/* constantes d'equilibre (ref = 1 mole/L) */
/* CxSyHz = xCa[2+] + 2xOH[-] + yH4SiO4 + (z-x-2y)H2O */
/* Stoichiometric Coefficients 
#define X_CSH(q)    courbe(q,el.mat->cb[0])
#define DX_CSH(q)   dcourbe(q,el.mat->cb[0])
#define V_CSH(q)    courbe(q,el.mat->cb[2])
#define DV_CSH(q)   dcourbe(q,el.mat->cb[2])
#define K_Sam      (1.93642e-3)      
#define F_SAM(q)    courbe(q,el.mat->cb[3])
#define DF_SAM(q)   dcourbe(q,el.mat->cb[3])
#define N_SOLIDCSH  (el.mat->nc - 3)
#define F_CSH(q,i)  courbe(q,el.mat->cb[3+i])

#ifdef NOTDEFINED
#define V_Sam      (29.e-3)      
#define K_Tob1     (5.53e-29)      
#define X_Tob1     (2)
#define Y_Tob1     (2.4)
#define V_Tob1     (131.47e-3)     
#define K_Tob2     (6.42e-22)        
#define X_Tob2     (1.5)
#define Y_Tob2     (1.8)
#define V_Tob2     (98.6e-3) 
#define K_Jen      (2.39e-16)  
#define X_Jen      (1.5)
#define Y_Jen      (0.9)
#define V_Jen      (73.6e-3)
#endif*/

#define C_H2SO4_eq   (K_h2o*K_h2o*K_CSH2/(K_hso4*K_so4*K_CH))


/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,elem_t,double) ;
static double poly4(double,double,double,double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
static double quartic_equation(double,double,double,double,double,double *,double *);
static double cubic_equation(double,double,double,double);
static void   quadratic_equation(double,double,double,double *);

#define MIN(a,b)   ((a < b) ? a : b)
#define MAX(a,b)   ((a > b) ? a : b)
#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)
#define DNEGEXP(y) ((y < 0.) ? exp(y) : 0.)

/* Ion Activity Products */
#define IAP_CSH2(zc_h2so4,zn_ca_s)           (K_CSH2*NEGEXP(zn_ca_s)*MIN(zc_h2so4,1.))
#define IAP_CH(zc_h2so4,zn_ca_s)           (K_CH*NEGEXP(zn_ca_s)/MAX(zc_h2so4,1.))
/*#define IAP_SAM(zp_ch,zn_si_s)           (K_Sam*NEGEXP(zn_si_s)*F_SAM(zp_ch))
#define DIAP_SAMSDQ_CH(zp_ch,zn_si_s)    (K_Sam*NEGEXP(zn_si_s)*DF_SAM(zp_ch))
#define DIAP_SAMSDZN_Si_S(zp_ch,zn_si_s) (K_Sam*DNEGEXP(zn_si_s)*F_SAM(zp_ch))*/

/* Parametres */
static double phi0,c_h2so4_eq,t_ch,t_csh2 ;
static double n_ca_ref ;/*n_si_ref*/


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"C_H2SO4_eq") == 0) return (2) ;
  else if(strcmp(s,"T_CH") == 0)     return (3) ;
  else if(strcmp(s,"T_CSH2") == 0)     return (4) ;
  else return(-1) ;
}


int dm70(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 7 ;
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_S],   "carbone") ;
  strcpy(mat->eqn[E_Ca],  "calcium") ;
  strcpy(mat->eqn[E_q],   "charge") ;


#if (U_H2SO4 == LOG_RHO)
  strcpy(mat->inc[I_H2SO4],  "logc_h2so4") ;
#else
  strcpy(mat->inc[I_H2SO4],  "c_h2so4") ;
#endif
  strcpy(mat->inc[I_Ca_S], "z_ca") ;
  strcpy(mat->inc[I_psi],  "psi") ;


  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_csh2        = 0. ;
    double n_ca_ref    = 1. ;

    mat->pr[pm("N_CH")]  = n_ca_ref ;
    mat->pr[pm("T_CH")]  = t_ch ;
    mat->pr[pm("T_CSH2")]  = t_csh2 ;

    dmat(mat,ficd,pm,n_donnees) ;

    t_ch      = mat->pr[pm("T_CH")] ;
    t_csh2      = mat->pr[pm("T_CSH2")] ;

    if(t_csh2  == 0.) mat->pr[pm("T_CSH2")]  = t_ch ;

    mat->pr[pm("C_H2SO4_eq")] = C_H2SO4_eq ;
  }
  
  return(mat->n) ;
}


int qm70(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  printf("Le systeme est forme de 5 equations:\n") ;
#if (U_H2SO4 == LOG_RHO)
  printf("\t- la conservation de la masse de C      (logc_h2so4)\n") ;
#else
  printf("\t- la conservation de la masse de C      (c_h2so4)\n") ;
#endif
  printf("\t- la conservation de la charge          (psi)\n") ;
  printf("\t- la conservation de la masse de Ca     (zn_ca_s)\n") ;


  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n") ;

  printf("Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1       # Contenu en CH (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5      # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CSH2  = 1.e5      # Cinetique de dissolution de CSH2 (s)\n") ;

  return(NEQ) ;
}



void tb70(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch70(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in70(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ca_ref  = el.mat->pr[pm("N_CH")] ;
  c_h2so4_eq  = el.mat->pr[pm("C_H2SO4_eq")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarities */
    double c_h2so4      = (C_H2SO4(i) > 0.) ? C_H2SO4(i) : c_h2so4_eq ;
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    
    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    
    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    /*double c_h2so4    = K_h2so4*c_so3 ;*/
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ; 
    
    double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    
    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_csh2_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CSH2 */
    double n_ch       = (zc_h2so4 <= 1) ? n_ch_eq  : 0 ;
    double n_csh2       = (zc_h2so4 >  1) ? n_csh2_eq  : 0 ;
    double n_ca_s     = n_ch + n_csh2 ;

    /* porosity */
    double phi = phi0 ;

    /* molar contents */
    double n_hso4     = phi*c_hso4 ;
    double n_h2so4    = phi*c_h2so4 ;
    double n_so4      = phi*c_so4 ;
    double n_ca       = phi*c_ca ;
    double n_cahso4   = phi*c_cahso4 ;
    double n_caso4aq  = phi*c_caso4aq ;
    double n_caoh     = phi*c_caoh ;
    
    N_S(i)  = n_h2so4 + n_hso4 + n_so4 + n_cahso4 + n_csh2 + n_caso4aq ;
    N_Ca(i) = n_ca + n_cahso4 + n_caso4aq + n_caoh + n_ca_s ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hso4*c_hso4 + z_so4*c_so4 + z_ca*c_ca + z_cahso4*c_cahso4 \
    + z_caoh*c_caoh ;

#if (U_H2SO4 == LOG_RHO)
    u[i][I_H2SO4] = log(c_h2so4)/Ln10 ;
#else
    C_H2SO4(i)   = c_h2so4 ;
#endif
    ZN_Ca_S(i) = zn_ca_s ;

    N_CH(i)    = n_ch ;
    N_CSH2(i)    = n_csh2 ;

    N_CH0(i)   = n_ch ;
    N_CSH20(i)   = n_csh2 ;
    ZP_CH0(i)  = zp_ch ;
  }
  
  if(el.dim < dim) return ;

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex70(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  
  if(el.dim < dim) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ca_ref  = el.mat->pr[pm("N_CH")] ;
  c_h2so4_eq  = el.mat->pr[pm("C_H2SO4_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;
  t_csh2      = el.mat->pr[pm("T_CSH2")] ;
  
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarities */
    double c_h2so4      = C_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;

    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;

    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ;
    double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* kinetics */
    double n_chn      = N_CHn(i) ;
    double n_csh2n      = N_CSH2n(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2so4,-dt/t_ch) ;  /* if zc_h2so4 > 1 */
    double n_csh2_ci    = n_csh2n*pow(zc_h2so4,dt/t_csh2) ;   /* if zc_h2so4 < 1 */

    /* equilibriums */
    double n_ch_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    double n_csh2_eq    = n_ca_ref*MAX(zn_ca_s,0.) ;
    
    /* solid contents : CH, CSH2 */
    double n_ch       = (zc_h2so4 <= 1) ? n_ch_eq  : n_ch_ci ;
    double n_csh2       = (zc_h2so4 >  1) ? n_csh2_eq  : n_csh2_ci ;
    double n_ca_s     = n_ch + n_csh2 ;

    /* porosity */
    double n_ch0      = N_CH0(i) ;
    double n_csh20      = N_CSH20(i) ;
    double zp_ch0      = ZP_CH0(i) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CSH2*(n_csh20 - n_csh2);

    /* molar contents */
    double n_h2so4    = phi*c_h2so4 ;
    double n_hso4     = phi*c_hso4 ;
    double n_so4      = phi*c_so4 ;
    double n_ca       = phi*c_ca ;
    double n_cahso4   = phi*c_cahso4 ;
    double n_caso4aq  = phi*c_caso4aq ;
    double n_caoh     = phi*c_caoh ;

    N_S(i)  = n_h2so4 + n_hso4 + n_so4 + n_cahso4 + n_csh2 + n_caso4aq ;
    N_Ca(i) = n_ca + n_cahso4 + n_caso4aq + n_caoh + n_ca_s ;

    /* charge density */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hso4*c_hso4 + z_so4*c_so4 + z_ca*c_ca + z_cahso4*c_cahso4  + \
    z_caoh*c_caoh ;

    N_CH(i)    = n_ch ;
    N_CSH2(i)    = n_csh2 ;

    if(c_h2so4 < 0. || n_ca_s < 0.) {
      printf("x         = %e\n",x[i][0]) ;
      printf("c_h2so4     = %e\n",c_h2so4) ;
      printf("n_csh2      = %e\n",n_csh2) ;
      printf("n_ch      = %e\n",n_ch) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("zn_ca_s   = %e\n",zn_ca_s) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
    if(phi < 0.) {
      printf("phi = %e\n",phi) ;
      printf("x         = %e\n",x[i][0]) ;
      printf("c_h2so4     = %e\n",c_h2so4) ;
      printf("c_h2so4_eq     = %e\n",c_h2so4_eq) ;
      printf("n_csh2      = %e\n",n_csh2) ;
      printf("n_ch      = %e\n",n_ch) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("zn_ca_s   = %e\n",zn_ca_s) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
  }
  
  if(el.dim < dim) return(0) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;

  return(0) ;
}

int mx70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double c[2] ;

  double Dc_hSDc_h2so4[2]        ;
  double Dc_ohSDc_h2so4[2]       ;
  double Dc_h2so4SDc_h2so4[2]    ;
  double Dc_hso4SDc_h2so4[2]     ;
  double Dc_so4SDc_h2so4[2]      ;
  double Dc_caSDc_h2so4[2]       ;
  double Dc_cahso4SDc_h2so4[2]   ;
  double Dc_caso4aqSDc_h2so4[2]  ;  
  double Dc_caohSDc_h2so4[2]     ;
  

  double Dc_hSDzn_ca_s[2]        ;
  double Dc_ohSDzn_ca_s[2]       ;
  double Dc_hso4SDzn_ca_s[2]     ;
  double Dc_so4SDzn_ca_s[2]      ;
  double Dc_caSDzn_ca_s[2]       ;
  double Dc_cahso4SDzn_ca_s[2]   ;
  double Dc_caso4aqSDzn_ca_s[2]  ;  
  double Dc_caohSDzn_ca_s[2]     ;
  
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ca_ref  = el.mat->pr[pm("N_CH")] ;

  c_h2so4_eq  = el.mat->pr[pm("C_H2SO4_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;

  t_csh2      = el.mat->pr[pm("T_CSH2")] ;


  /*
    Initialisation 
  */
  for(i=0;i<4*NEQ*NEQ;i++) k[i] = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    CALCUL DE volume ET DE surf 
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])*0.5 ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)*0.5 ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = 2*M_PI*xm ; else surf = 1 ;
  /*
    termes d'accumulation
  */
  for(i=0;i<2;i++) {
    /* molarites */
    double c_h2so4      = C_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;

    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
 
    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ;
    double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* solid contents : CH, CSH2, CSH */
    double n_csh2       = N_CSH2(i) ;
    double n_ch       = N_CH(i) ;

    /* kinetics */
    double n_csh2n      = N_CSH2n(i) ;
    double n_chn      = N_CHn(i) ;
    double n_ch_ci    = n_chn*pow(zc_h2so4,-dt/t_ch) ;
    double n_csh2_ci    = n_csh2n*pow(zc_h2so4,dt/t_csh2) ;

    /* porosity */
    double n_ch0      = N_CH0(i) ;
    double n_csh20      = N_CSH20(i) ;
    double zp_ch0      = ZP_CH0(i) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CSH2*(n_csh20 - n_csh2);

    /* derivatives ... */
    /* ... with respect to c_h2so4 */
    double dc_h2so4             = 1.e-4*c_h2so4*((c_h2so4 > C_H2SO4n(i)) ? 1 : -1) ;
    double c_h2so42             = c_h2so4 + dc_h2so4 ;
    double c_oh2              = concentration_oh(c_h2so42,el,zn_ca_s) ;   
    double dc_ohsdc_h2so4       = (c_oh2 - c_oh)/dc_h2so4 ;
    double dc_hsdc_h2so4        = - c_h*dc_ohsdc_h2so4/c_oh ;
    double dc_h2so4sdc_h2so4    = 1 ;
    double dc_hso4sdc_h2so4     = c_hso4*(dc_h2so4sdc_h2so4/c_h2so4	\
                              + dc_ohsdc_h2so4/c_oh) ;
    double dc_so4sdc_h2so4      = c_so4*(dc_hso4sdc_h2so4/c_hso4 \
                              + dc_ohsdc_h2so4/c_oh) ;
			                     
    double dP_CSH2sdc_h2so4       = (zc_h2so4 < 1.) ? P_CSH2/c_h2so4 : 0. ;
    double dP_CHsdc_h2so4       = (zc_h2so4 < 1.) ? 0. : -P_CH/c_h2so4 ;
    
    double dc_casdc_h2so4       = (dP_CSH2sdc_h2so4 - c_ca*dc_so4sdc_h2so4)/c_so4 ;
    double dc_cahso4sdc_h2so4   = K_cahso4*(dc_casdc_h2so4*c_hso4 \
                              + c_ca*dc_hso4sdc_h2so4) ;


    double dc_caso4aqsdc_h2so4  = K_caso4aq*(dc_so4sdc_h2so4*c_ca \
                              + dc_casdc_h2so4*c_so4) ;
    double dc_caohsdc_h2so4     = K_caoh*(dc_casdc_h2so4*c_oh \
                              + dc_ohsdc_h2so4*c_ca) ;
		                    
    double dn_ch_cisdc_h2so4    = n_ch_ci*(-dt/t_ch)/c_h2so4 ;
    double dn_csh2_cisdc_h2so4    = n_csh2_ci*(dt/t_csh2)/c_h2so4 ;

    double dn_ch_eqsdc_h2so4    = 0 ;
    double dn_csh2_eqsdc_h2so4    = 0 ;


    double dn_chsdc_h2so4       = (zc_h2so4 <= 1) ? dn_ch_eqsdc_h2so4  : dn_ch_cisdc_h2so4 ;
    double dn_csh2sdc_h2so4       = (zc_h2so4 >  1) ? dn_csh2_eqsdc_h2so4  : dn_csh2_cisdc_h2so4 ;
    
    double dn_ca_ssdc_h2so4     = dn_chsdc_h2so4 + dn_csh2sdc_h2so4;

    double dphisdc_h2so4        = - V_CH*dn_chsdc_h2so4 - V_CSH2*dn_csh2sdc_h2so4 ;

    /*double dn_h2so4sdc_h2so4      = phi + dphisdc_h2so4*c_h2so4 ;*/
    double dn_h2so4sdc_h2so4    = phi*dc_h2so4sdc_h2so4 + dphisdc_h2so4*c_h2so4 ;
    double dn_hso4sdc_h2so4     = phi*dc_hso4sdc_h2so4 + dphisdc_h2so4*c_hso4 ;
    double dn_so4sdc_h2so4      = phi*dc_so4sdc_h2so4 + dphisdc_h2so4*c_so4 ;
    double dn_casdc_h2so4       = phi*dc_casdc_h2so4 + dphisdc_h2so4*c_ca ;
    double dn_cahso4sdc_h2so4   = phi*dc_cahso4sdc_h2so4 + dphisdc_h2so4*c_cahso4 ;
    double dn_caso4aqsdc_h2so4  = phi*dc_caso4aqsdc_h2so4 + dphisdc_h2so4*c_caso4aq ;
    double dn_caohsdc_h2so4     = phi*dc_caohsdc_h2so4 + dphisdc_h2so4*c_caoh ;			                 
    
    /* with respect to zn_ca_s */
    /* double dzn_ca_s             = ((zn_ca_s > 0.) ? 1 : -1)*1.e-2 ; */
    double dzn_ca_s            = 1.e-6*((zn_ca_s > ZN_Ca_Sn(i)) ? 1 : -1) ;
    double zn_ca_s2            = zn_ca_s + dzn_ca_s ;
    double c_oh3               = concentration_oh(c_h2so4,el,zn_ca_s2) ;
    double dc_ohsdzn_ca_s      = (c_oh3 - c_oh)/dzn_ca_s ;
    double dc_hsdzn_ca_s       = - c_h*dc_ohsdzn_ca_s/c_oh ;
    double dc_hso4sdzn_ca_s    = c_hso4*(dc_ohsdzn_ca_s/c_oh) ;
    double dc_so4sdzn_ca_s     = c_so4*(dc_hso4sdzn_ca_s/c_hso4 \
                               + dc_ohsdzn_ca_s/c_oh) ;

    double dP_CSH2sdzn_ca_s      = (zn_ca_s < 0) ? P_CSH2 : 0 ;
    double dP_CHsdzn_ca_s      = (zn_ca_s < 0) ? P_CH : 0 ;
    
    double dc_casdzn_ca_s      = (dP_CSH2sdzn_ca_s - c_ca*dc_so4sdzn_ca_s)/c_so4 ;
    double dc_cahso4sdzn_ca_s  = K_cahso4*(dc_casdzn_ca_s*c_hso4 \
                               + c_ca*dc_hso4sdzn_ca_s) ;
    

    double dc_caso4aqsdzn_ca_s = K_caso4aq*(dc_so4sdzn_ca_s*c_ca \
                               + dc_casdzn_ca_s*c_so4) ;
    double dc_caohsdzn_ca_s    = K_caoh*(dc_casdzn_ca_s*c_oh \
                               + dc_ohsdzn_ca_s*c_ca) ;					       

    double dn_ch_cisdzn_ca_s   = 0 ;
    double dn_csh2_cisdzn_ca_s   = 0 ;


    double dn_ch_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
    double dn_csh2_eqsdzn_ca_s   = (zn_ca_s > 0) ? n_ca_ref : 0 ;
     
    double dn_chsdzn_ca_s      = (zc_h2so4 <= 1) ? dn_ch_eqsdzn_ca_s  : dn_ch_cisdzn_ca_s ;
    double dn_csh2sdzn_ca_s      = (zc_h2so4 >  1) ? dn_csh2_eqsdzn_ca_s  : dn_csh2_cisdzn_ca_s ;
    
    double dn_ca_ssdzn_ca_s    = dn_chsdzn_ca_s + dn_csh2sdzn_ca_s ;

    double dphisdzn_ca_s       = - V_CH*dn_chsdzn_ca_s \
                               - V_CSH2*dn_csh2sdzn_ca_s  ;
              
    /*double dn_so3sdzn_ca_s     = dphisdzn_ca_s*c_so3 ;*/
    double dn_h2so4sdzn_ca_s   = dphisdzn_ca_s*c_h2so4 ;
    double dn_hso4sdzn_ca_s    = phi*dc_hso4sdzn_ca_s + dphisdzn_ca_s*c_hso4 ;
    double dn_so4sdzn_ca_s     = phi*dc_so4sdzn_ca_s + dphisdzn_ca_s*c_so4 ;
    double dn_casdzn_ca_s      = phi*dc_casdzn_ca_s + dphisdzn_ca_s*c_ca ;
    double dn_cahso4sdzn_ca_s  = phi*dc_cahso4sdzn_ca_s + dphisdzn_ca_s*c_cahso4 ;
    double dn_caso4aqsdzn_ca_s = phi*dc_caso4aqsdzn_ca_s + dphisdzn_ca_s*c_caso4aq ;
    double dn_caohsdzn_ca_s    = phi*dc_caohsdzn_ca_s + dphisdzn_ca_s*c_caoh ;
    
  

    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_S1 - n_Sn) + dt * div(w_S) = 0
    */
    K(E_S+j,I_H2SO4+j)   += volume[i]*(dn_h2so4sdc_h2so4 + dn_hso4sdc_h2so4 + dn_so4sdc_h2so4 + dn_cahso4sdc_h2so4 \
                        + dn_caso4aqsdc_h2so4) ;
    K(E_S+j,I_Ca_S+j)   += volume[i]*(dn_h2so4sdzn_ca_s + dn_hso4sdzn_ca_s + dn_so4sdzn_ca_s + dn_cahso4sdzn_ca_s \
                        + dn_caso4aqsdzn_ca_s + dn_csh2sdzn_ca_s) ;
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_H2SO4+j)  += volume[i]*(dn_casdc_h2so4 + dn_cahso4sdc_h2so4 + dn_ca_ssdc_h2so4 \
                        + dn_caso4aqsdc_h2so4 +dn_caohsdc_h2so4) ;
    K(E_Ca+j,I_Ca_S+j)  += volume[i]*(dn_casdzn_ca_s + dn_cahso4sdzn_ca_s + dn_ca_ssdzn_ca_s \
                        + dn_caso4aqsdzn_ca_s +dn_caohsdzn_ca_s ) ;
    
    /*
      Conservation de la charge  : div(w_q) = 0
    */

    /* sauvegardes pour les termes de transport */
    Dc_hSDc_h2so4[i]        = dc_hsdc_h2so4 ;
    Dc_ohSDc_h2so4[i]       = dc_ohsdc_h2so4 ;
    Dc_h2so4SDc_h2so4[i]    = dc_h2so4sdc_h2so4 ;
    Dc_hso4SDc_h2so4[i]     = dc_hso4sdc_h2so4 ;
    Dc_so4SDc_h2so4[i]      = dc_so4sdc_h2so4 ;
    Dc_caSDc_h2so4[i]       = dc_casdc_h2so4 ;
    Dc_cahso4SDc_h2so4[i]   = dc_cahso4sdc_h2so4 ;
    Dc_caso4aqSDc_h2so4[i]  = dc_caso4aqsdc_h2so4 ;
    Dc_caohSDc_h2so4[i]     = dc_caohsdc_h2so4 ;
    
    Dc_hSDzn_ca_s[i]        = dc_hsdzn_ca_s ;
    Dc_ohSDzn_ca_s[i]       = dc_ohsdzn_ca_s ;
    Dc_hso4SDzn_ca_s[i]     = dc_hso4sdzn_ca_s ;
    Dc_so4SDzn_ca_s[i]      = dc_so4sdzn_ca_s ;
    Dc_caSDzn_ca_s[i]       = dc_casdzn_ca_s ;      
    Dc_cahso4SDzn_ca_s[i]   = dc_cahso4sdzn_ca_s ;
    Dc_caso4aqSDzn_ca_s[i]  = dc_caso4aqsdzn_ca_s ;
    Dc_caohSDzn_ca_s[i]     = dc_caohsdzn_ca_s ;
  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trf_oh       = tr*KF_OH ;
  double trf_h        = tr*KF_H ;
  double trf_so3      = tr*KF_H2SO4 ;
  double trf_h2so4    = tr*KF_H2SO4 ;
  double trf_hso4     = tr*KF_HSO4 ;
  double trf_so4      = tr*KF_SO4 ;
  double trf_ca       = tr*KF_Ca ;
  double trf_cahso4   = tr*KF_CaHSO4 ;
  double trf_caso4aq  = tr*KF_CaSO4aq ;
  double trf_caoh     = tr*KF_CaOH ;
  
  double tre_hso4     = tr*Kpsi_HSO4 ;
  double tre_so4      = tr*Kpsi_SO4 ;
  double tre_ca       = tr*Kpsi_Ca ;
  double tre_cahso4   = tr*Kpsi_CaHSO4 ;
  double tre_caoh     = tr*Kpsi_CaOH ;

  double tre_q        = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_S1 - n_Sn) + dt * div(w_S) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h2so4*Dc_h2so4SDc_h2so4[i] + trf_hso4*Dc_hso4SDc_h2so4[i] + trf_so4*Dc_so4SDc_h2so4[i] + trf_cahso4*Dc_cahso4SDc_h2so4[i] \
           + trf_caso4aq*Dc_caso4aqSDc_h2so4[i] ;
  }
  K(E_S,I_H2SO4)          += + c[0] ;
  K(E_S,I_H2SO4+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_H2SO4)      += - c[0] ;
  K(E_S+NEQ,I_H2SO4+NEQ)  += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = trf_hso4*Dc_hso4SDzn_ca_s[i] + trf_so4*Dc_so4SDzn_ca_s[i] + trf_cahso4*Dc_cahso4SDzn_ca_s[i] + trf_caso4aq*Dc_caso4aqSDzn_ca_s[i] ;
  }
  K(E_S,I_Ca_S)           += + c[0] ;
  K(E_S,I_Ca_S+NEQ)       += - c[1] ;
  K(E_S+NEQ,I_Ca_S)       += - c[0] ;
  K(E_S+NEQ,I_Ca_S+NEQ)   += + c[1] ;
  

  for(i=0;i<2;i++){
    c[i] = tre_hso4 + tre_so4 + tre_cahso4 ;
  }
  K(E_S,I_psi)          += + c[0] ;
  K(E_S,I_psi+NEQ)      += - c[1] ;
  K(E_S+NEQ,I_psi)      += - c[0] ;
  K(E_S+NEQ,I_psi+NEQ)  += + c[1] ;

  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_h2so4[i] + trf_cahso4*Dc_cahso4SDc_h2so4[i]  \
    + trf_caoh*Dc_caohSDc_h2so4[i] + trf_caso4aq*Dc_caso4aqSDc_h2so4[i] ;
  }
  K(E_Ca,I_H2SO4)         += + c[0] ;
  K(E_Ca,I_H2SO4+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_H2SO4)     += - c[0] ;
  K(E_Ca+NEQ,I_H2SO4+NEQ) += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_ca_s[i] + trf_cahso4*Dc_cahso4SDzn_ca_s[i]  \
    + trf_caoh*Dc_caohSDzn_ca_s[i] + trf_caso4aq*Dc_caso4aqSDzn_ca_s[i] ;
  }
  K(E_Ca,I_Ca_S)         += + c[0] ;
  K(E_Ca,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Ca_S+NEQ) += + c[1] ;
  

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_cahso4  + tre_caoh;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDc_h2so4[i] + z_oh*trf_oh*Dc_ohSDc_h2so4[i] + z_hso4*trf_hso4*Dc_hso4SDc_h2so4[i] + z_so4*trf_so4*Dc_so4SDc_h2so4[i] \
           + z_ca*trf_ca*Dc_caSDc_h2so4[i] + z_cahso4*trf_cahso4*Dc_cahso4SDc_h2so4[i] \
           + z_caoh*trf_caoh*Dc_caohSDc_h2so4[i] ;
  }
  K(E_q,I_H2SO4)           += + c[0] ;
  K(E_q,I_H2SO4+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_H2SO4)       += - c[0] ;
  K(E_q+NEQ,I_H2SO4+NEQ)   += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_ca_s[i] + z_oh*trf_oh*Dc_ohSDzn_ca_s[i] + z_hso4*trf_hso4*Dc_hso4SDzn_ca_s[i] + z_so4*trf_so4*Dc_so4SDzn_ca_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_ca_s[i] + z_cahso4*trf_cahso4*Dc_cahso4SDzn_ca_s[i] \
           + z_caoh*trf_caoh*Dc_caohSDzn_ca_s[i] ;
  }
  K(E_q,I_Ca_S)           += + c[0] ;
  K(E_q,I_Ca_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Ca_S)       += - c[0] ;
  K(E_q+NEQ,I_Ca_S+NEQ)   += + c[1] ;
  

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  
  
 }

#if (U_H2SO4 == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_H2SO4)     *= Ln10*C_H2SO4(0) ;
    K(i,I_H2SO4+NEQ) *= Ln10*C_H2SO4(1) ;
  }
#endif
  

  return(0) ;

#undef K
}


void rs70(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ*2;i++) r[i] = zero ;

  if(el.dim < dim) return ;
  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    Conservation de C (carbone) : (n_S1 - n_Sn) + dt * div(w_S) = 0
  */
  R(0,E_S) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S ;
  R(1,E_S) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;

#undef R
}


int so70(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
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
  phi0     = el.mat->pr[pm("porosite")] ;
  c_h2so4_eq = el.mat->pr[pm("C_H2SO4_eq")] ;

  /* initialisation */
  nso = 18 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* output quantities */
  {
    /* molarities */
#if (U_H2SO4 == LOG_RHO)
    double c_h2so4      = exp(Ln10*param(u,h_s,el.nn,I_H2SO4)) ;
#else
    double c_h2so4      = param(u,h_s,el.nn,I_H2SO4) ;
#endif
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;
    double zn_ca_s    = param(u,h_s,el.nn,I_Ca_S) ;
    
    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
    
    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ;    
    double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* charge density */
    double c_q = z_h*c_h + z_oh*c_oh + z_hso4*c_hso4 + z_so4*c_so4 + z_ca*c_ca + z_cahso4*c_cahso4 \
    + z_caoh*c_caoh ;
    /* solid contents */
    int j = (el.dim < dim) ? 0 : ((s[0] < (x[0][0] + x[1][0])*0.5) ? 0 : 1) ;
    double n_ch       = N_CH(j) ;
    double n_csh2       = N_CSH2(j) ;
    
    /* porosity */
    double n_ch0      = N_CH0(j) ;
    double n_csh20      = N_CSH20(j) ;
    double zp_ch0     = ZP_CH0(j) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CSH2*(n_csh20 - n_csh2);

    double psi        = param(u,h_s,el.nn,I_psi) ;

    i = 0 ;
    strcpy(r[i].text,"c_h2so4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2so4 ;
    strcpy(r[i].text,"ph") ; r[i].n = 1 ;
    r[i++].v[0] = 14 + log(c_oh)/log(10.) ;
    strcpy(r[i].text,"c_h2so4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2so4 ;
    strcpy(r[i].text,"c_hso4") ; r[i].n = 1 ;
    r[i++].v[0] = c_hso4 ;
    strcpy(r[i].text,"c_so4") ; r[i].n = 1 ;
    r[i++].v[0] = c_so4 ;
    strcpy(r[i].text,"c_ca") ; r[i].n = 1 ;
    r[i++].v[0] = c_ca ;
    strcpy(r[i].text,"c_cahso4") ; r[i].n = 1 ;
    r[i++].v[0] = c_cahso4 ;
    strcpy(r[i].text,"n_ch") ; r[i].n = 1 ;
    r[i++].v[0] = n_ch ;
    strcpy(r[i].text,"n_csh2") ; r[i].n = 1 ;
    r[i++].v[0] = n_csh2 ;
    strcpy(r[i].text,"porosite") ; r[i].n = 1 ;
    r[i++].v[0] = phi ;
    strcpy(r[i].text,"c_oh") ; r[i].n = 1 ;
    r[i++].v[0] = c_oh ;
    strcpy(r[i].text,"potentiel_electrique") ; r[i].n = 1 ;
    r[i++].v[0] = psi ;
    strcpy(r[i].text,"charge") ; r[i].n = 1 ;
    r[i++].v[0] = c_q ;
    strcpy(r[i].text,"zn_ca_s") ; r[i].n = 1 ;
    r[i++].v[0] = zn_ca_s ;
    strcpy(r[i].text,"zc_h2so4") ; r[i].n = 1 ;
    r[i++].v[0] = zc_h2so4 ;
    strcpy(r[i].text,"c_caso4aq") ; r[i].n = 1 ;
    r[i++].v[0] = c_caso4aq ;
    strcpy(r[i].text,"c_caoh") ; r[i].n = 1 ;
    r[i++].v[0] = c_caoh ;
    strcpy(r[i].text,"pk_ch") ; r[i].n = 1 ;
    r[i++].v[0] = log10(zp_ch) ;

  }
  
  if(i != nso) arret("so53") ;

  return(nso) ;
}


void transfert(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom) 
/* Termes explicites (va)  */
{
  int    i ;
  /*
    Donnees
  */
  phi0     = el.mat->pr[pm("porosite")] ;
  c_h2so4_eq = el.mat->pr[pm("C_H2SO4_eq")] ;

  /* initialisation */
  for(i=0;i<NVE_TR;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* molarites */
    double c_h2so4      = C_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;

    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ; 
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;

    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    /* solid contents */
    double n_csh2       = N_CSH2(i) ;
    double n_ch       = N_CH(i) ;

    /* porosity */
    double n_csh20      = N_CSH20(i) ;
    double n_ch0      = N_CH0(i) ;
    double zp_ch0     = ZP_CH0(i) ;

    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CSH2*(n_csh20 - n_csh2) ;

    /* tortuosite liquide */
    double iff    = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_SO3        += d_so3*iff ;
    KF_H2SO4      += d_h2so4*iff ;
    KF_HSO4       += d_hso4*iff ;
    KF_SO4        += d_so4*iff ;

    KF_Ca         += d_ca*iff ;
    KF_CaHSO4     += d_cahso4*iff ;
    
    KF_CaSO4aq    += d_caso4aq*iff ;
    KF_CaOH       += d_caoh*iff ;


    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HSO4     += FsRT*KF_HSO4*z_hso4*c_hso4 ;
    Kpsi_SO4      += FsRT*KF_SO4*z_so4*c_so4 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHSO4   += FsRT*KF_CaHSO4*z_cahso4*c_cahso4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hso4*Kpsi_HSO4 + z_so4*Kpsi_SO4 + z_ca*Kpsi_Ca + z_cahso4*Kpsi_CaHSO4  \
                  + z_caoh*Kpsi_CaOH ;
  }
  
  /* moyenne */
  for(i=0;i<NVE_TR;i++) va[i] *= 0.5 ;
}


void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les flux (f) */
{
  double r_h[2],r_oh[2] ;
  double r_h2so4[2],r_hso4[2],r_so4[2],r_caso4aq[2],r_caoh[2] ;
  double r_ca[2],r_cahso4[2] ;

  int    i ;

  for(i=0;i<2;i++) {
    double c_h2so4      = C_H2SO4(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double zc_h2so4     = c_h2so4/c_h2so4_eq ;

    double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;    
    double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
    double zp_ch      = P_CH/K_CH ;
 
    double c_oh       = concentration_oh(c_h2so4,el,zn_ca_s) ;
    double c_h        = K_h2o/c_oh ;
    double c_hso4     = K_hso4*c_h2so4/c_h ;
    double c_so4      = K_so4*c_hso4/c_h ;
    double c_ca       = P_CSH2/c_so4 ;
    double c_cahso4   = K_cahso4*c_ca*c_hso4 ;
    double c_caso4aq  = K_caso4aq*c_ca*c_so4 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;

    r_oh[i]       = c_oh ;
    r_h[i]        = c_h ;
    r_h2so4[i]    = c_h2so4 ;
    r_hso4[i]     = c_hso4 ;
    r_so4[i]      = c_so4 ;
    r_ca[i]       = c_ca ;
    r_cahso4[i]   = c_cahso4 ;
    r_caoh[i]     = c_caoh ;
    r_caso4aq[i]  = c_caso4aq ;
  }

  /* Gradients */
  {
    double dx = x[1][0] - x[0][0] ;
    double grd_h        = (r_h[1]        - r_h[0]       )/dx ;
    double grd_oh       = (r_oh[1]       - r_oh[0]      )/dx ;
    double grd_h2so4    = (r_h2so4[1]    - r_h2so4[0]   )/dx ;
    double grd_hso4     = (r_hso4[1]     - r_hso4[0]    )/dx ;
    double grd_so4      = (r_so4[1]      - r_so4[0]     )/dx ;
    double grd_cahso4   = (r_cahso4[1]   - r_cahso4[0]  )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    double grd_caso4aq  = (r_caso4aq[1]  - r_caso4aq[0] )/dx ;
    double grd_caoh     = (r_caoh[1]     - r_caoh[0]    )/dx ;
    
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    
    /* Flux */
    double w_h2so4    = - KF_H2SO4*grd_h2so4   ;
    double w_hso4     = - KF_HSO4*grd_hso4          - Kpsi_HSO4*grd_psi ;
    double w_so4      = - KF_SO4*grd_so4            - Kpsi_SO4*grd_psi      ;
    double w_cahso4   = - KF_CaHSO4*grd_cahso4      - Kpsi_CaHSO4*grd_psi ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
    double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
    double w_caso4aq  = - KF_CaSO4aq*grd_caso4aq ;    
    
    double w_q        = - z_h*KF_H*grd_h		      \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hso4*KF_HSO4*grd_hso4             \
                        - z_so4*KF_SO4*grd_so4		      \
                        - z_ca*KF_Ca*grd_ca		      \
                        - z_cahso4*KF_CaHSO4*grd_cahso4	      \                   
                        - z_caoh*KF_CaOH*grd_caoh \
                        - Kpsi_q*grd_psi ;

    W_S     = w_h2so4 + w_hso4 + w_so4 + w_cahso4 + w_caso4aq ;
    W_Ca    = w_ca + w_cahso4 + w_caso4aq + w_caoh ;
    W_q     = w_q ;
  }
}


double concentration_oh(double c_h2so4,elem_t el,double zn_ca_s)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_h2so4, zn_si_s et zn_ca_s sont fixes */
  double c_h2so4_eq   = C_H2SO4_eq ;
  double zc_h2so4     = c_h2so4/c_h2so4_eq ;
  /* les produits sont donc aussi fixes */
  double P_CSH2       = IAP_CSH2(zc_h2so4,zn_ca_s) ;
  double P_CH       = IAP_CH(zc_h2so4,zn_ca_s) ;
  double zp_ch      = P_CH/K_CH ;

  /*
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                     : +1
     c_hso4     = K_hso4*c_h2so4/c_h             : -1
     c_so4      = K_so4*c_hso4/c_h               : -2
     c_ca       = P_CSH2/c_so4                     : +2
     c_cahso4   = K_cahso4*c_ca*c_hso4           : +1
     c_h4sio4   = P_Sam                          :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_caso4aq  = K_caso4aq*c_ca*c_so4 ;         :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;       : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;     :  0      
     c_caoh     = K_caoh*c_ca*c_oh ;             : +1       
  */
  double A_hso4     = K_hso4*c_h2so4 ;
  double A_so4      = K_so4*A_hso4 ;
  double A_ca       = P_CSH2/A_so4 ;
  double A_cahso4   = K_cahso4*A_ca*A_hso4 ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahso4*A_cahso4  + z_caoh*A_caoh ;
  double c = 0 ;
  double d = z_oh*K_h2o + z_hso4*A_hso4  ;
  double e = z_so4*A_so4 ;

  /* a > 0 ; b > 0 ; c > 0 ; d < 0 ; e < 0 */
  double c_h = poly4(a,b,c,d,e) ;
 
  return(K_h2o/c_h) ;
}

double poly4(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double iks[4],zed[2];
  double x0 = quartic_equation(a,b,c,d,e,iks,zed); /*quartic_equation(a,b,c,d,e,iks,zed)*/
  double x  = MAX(x0,0.) ;
  int    i = 0 ;
  /* 
     on neglige a et c
     solution de bx^3 + dx + e = 0 
     que l'on met sous la forme x^3 + px + q = 0
  double p = d/b,q = e/b ;
  double r = sqrt(q*q/4 + p*p*p/27) ;
  double uu = pow(-q/2. + r,1/3.), vv = pow(-q/2. - r,1/3.) ;
  double x = uu + vv ; *//* formule de Cardan */
  
  /* Newton */
  do {
    double x2 = x*x,x3 = x2*x,x4 = x2*x2 ;
    double f  = a*x4 + b*x3 + c*x2 + d*x + e ;
    double df = 4*a*x3 + 3*b*x2 + 2*c*x + d ;
    double dx = -f/df ;
    err = fabs(dx/x) ;
    x += dx ;
    if(i++ > 20) {
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      arret("poly_sansc : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}

double cubic_equation(double a,double b,double c,double d)
{ /* solution de ax^3 + bx^2+ cx + d = 0 */ 
  double x,y1,y2,y3,p,q,r,v,uu,vv,x1,x2,x3,k;  

  p = c/a-b*b/(3.0*a*a);
  q = d/a + 2.0*pow(b,3.0)/(27.0*pow(a,3.0)) - b*c/(3.0*a*a);
  v = q*q/4.0 + p*p*p/27.0 ;
  if (v>=0){
  
  r = sqrt(v) ;
  uu = pow(-q/2. + r,1/3.), vv = ((-q/2. - r)>0.) ? pow(-q/2. - r,1/3.) : -pow(q/2. + r,1/3.);
  
  y1 = uu + vv ;
  x1 = y1 - b/(3.0*a) ;
  x2 = 0. ;
  x3 = 0.;
  }  
  else {
  k = acos(-1.5*q*sqrt(-3.*p)/(p*p)); 
  y1 = 2.*sqrt(-3.*p)*cos(k/3.)/3.;
  y2 = -sqrt(-3.*p) * (cos(k/3.) + sqrt(3.)*sin(k/3.)) /3.;
  y3 = -sqrt(-3.*p) * (cos(k/3.) - sqrt(3.)*sin(k/3.)) /3.;
  x1 = y1 - b/(3.*a) ;
  x2 = y2 - b/(3.*a) ;
  x3 = y3 - b/(3.*a) ;
  }
  x = MAX(x1,x2) ;
  return(MAX(x,x3)) ;
}

void quadratic_equation(double a,double b,double c,double *zed)
{ /* solution de ax^2+ bx + c = 0 */ 
  double delta,temp1,temp2;  
  delta = b*b - 4.*a*c;
  if (delta>=0.){
  temp1 = -b/(2.*a);
  temp2 = sqrt(delta)/(2.*a);
  zed[0] = temp1 + temp2;
  zed[1] = temp1 - temp2; 
  } 
  else{
  zed[0] = 0.;
  zed[1] = 0.;
  }     
}

double quartic_equation(double a,double b,double c,double d,double e,double *iks,double *zed)
{ /* solution de ax^4 + bx^3+ cx^2 + dx + e = 0 */ 
  double y,a2,b2,c2,a4,b4,c4,a3,b3,c3,d3,a1,b1,c1,d1,e1,x1,x2,x3,x4,x5,x6,x;
  
  a1 = a/a;
  b1 = b/a;
  c1 = c/a;
  d1 = d/a;
  e1 = e/a; 
  if(b==0. && c==0. && d==0. && e==0.){
  iks[0]=0.;
  iks[1]=0.;
  iks[2]=0.;
  iks[3]=0.;
  } 
  else if(b==0. && d==0. && e==0.){
  quadratic_equation(a,0.,c,zed);
  iks[0] = zed[0];
  iks[1] = zed[1];
  iks[2] = 0.;
  iks[3] = 0.;
  }   
  else{
  a3 = 8.;
  b3 = -4.*c/a;
  c3 = 2.0*b*d/(a*a)- 8.*e/a;
  d3 = e*(4.*c/a-b*b/(a*a))/a - d*d/(a*a);
  y = cubic_equation(a3,b3,c3,d3);
  a2 = 1.;
  b2 = b1/2. - sqrt(8.*y + b1*b1 - 4.*c1)/2.;
  c2 = y - (b1*y - d1)/sqrt(8.*y + b1*b1 - 4.*c1);
  quadratic_equation(a2,b2,c2,zed);
  iks[0]=zed[0];
  iks[1]=zed[1];
  a4 = 1.;
  b4 = b1/2. + sqrt(8.*y + b1*b1 - 4.*c1)/2.;
  c4 = y + (b1*y - d1)/sqrt(8.*y + b1*b1 - 4.*c1);
  quadratic_equation(a4,b4,c4,zed);
  iks[2]=zed[0];
  iks[3]=zed[1];
  }
  x1 = (iks[0]>0. && iks[0]<1.)? iks[0] : 0. ;
  x2 = (iks[0]>0. && iks[0]<1.)? iks[1] : 0. ;
  x3 = (iks[0]>0. && iks[0]<1.)? iks[2] : 0. ;
  x4 = (iks[0]>0. && iks[0]<1.)? iks[3] : 0. ;
  x5 = MAX(x1,x2);
  x6 = MAX(x3,x4);
  x  = MAX(x5,x6);
  return(x);
  
}



