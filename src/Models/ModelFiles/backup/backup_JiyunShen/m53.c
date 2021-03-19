/*

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  53
#define TITLE "sc-Carbonation of saturated CBM (05/2009)"
#define AUTHORS "Shen"

#include "OldMethods.h"

/* Macros */
#define NEQ     (4)
#define NVE     (33)
#define NVI     (20)
#define NVE_TR  (25)

#define E_C     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)

#define I_CO2   (0)
#define I_psi   (1)
#define I_CHCC  (2)
#define I_CSH   (3)

#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_CO2   LOG_RHO

#if (U_CO2 == LOG_RHO)
  #define C_CO2(n)   (exp(Ln10*u[(n)][I_CO2]))
  #define C_CO2n(n)  (exp(Ln10*u_n[(n)][I_CO2]))
#else
  #define C_CO2(n)   (u[(n)][I_CO2])
  #define C_CO2n(n)  (u_n[(n)][I_CO2])
#endif
#define Z_CHCC(n)  (u[(n)][I_CHCC])
#define Z_CSH(n)   (u[(n)][I_CSH])
#define PSI(n)     (u[(n)][I_psi])

#define Z_CHCCn(n) (u_n[(n)][I_CHCC])
#define Z_CSHn(n)  (u_n[(n)][I_CSH])
#define PSIn(n)    (u_n[(n)][I_psi])

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define W_C        (f[8])
#define W_q        (f[9])
#define W_Ca       (f[10])
#define W_Si       (f[11])
#define N_CH(n)    (f[(12+n)])
#define N_CC(n)    (f[(14+n)])
#define N_CSH(n)   (f[(16+n)])
#define N_Sam(n)   (f[(18+n)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_CHn(n)   (f_n[(12+n)])
#define N_CCn(n)   (f_n[(14+n)])
#define N_CSHn(n)  (f_n[(16+n)])
#define N_Samn(n)  (f_n[(18+n)])


#define KF_OH         (va[(0)])
#define KF_H          (va[(1)])
#define KF_CO2        (va[(2)])
#define KF_H2CO3      (va[(3)])
#define KF_HCO3       (va[(4)])
#define KF_CO3        (va[(5)])
#define KF_Ca         (va[(6)])
#define KF_CaOH       (va[(7)])
#define KF_CaHCO3     (va[(8)])
#define KF_CaH2SiO4   (va[(9)])
#define KF_CaH3SiO4   (va[(10)])
#define KF_H2SiO4     (va[(11)])
#define KF_H3SiO4     (va[(12)])
#define KF_H4SiO4     (va[(13)])

#define Kpsi_OH       (va[(14)])
#define Kpsi_H        (va[(15)])
#define Kpsi_HCO3     (va[(16)])
#define Kpsi_CO3      (va[(17)])
#define Kpsi_Ca       (va[(18)])
#define Kpsi_CaOH     (va[(19)])
#define Kpsi_CaHCO3   (va[(20)])
#define Kpsi_CaH3SiO4 (va[(21)])
#define Kpsi_H2SiO4   (va[(22)])
#define Kpsi_H3SiO4   (va[(23)])
#define Kpsi_q        (va[(24)])

#define N_CH0(n)      (va[(25+n)])
#define N_CC0(n)      (va[(27+n)])
#define N_CSH0(n)     (va[(29+n)])
#define N_Sam0(n)     (va[(31+n)])

/*
  Solution aqueuse
*/

/* les valences */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_hco3     (-1.)
#define z_co3      (-2.)
#define z_h2sio4   (-2.)
#define z_h3sio4   (-1.)
#define z_caoh     (1.)
#define z_cahco3   (1.)
#define z_cah3sio4 (1.)

/* volumes molaires partiels des ions (dm3/mole) */
#define v_h        (-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       (-4.71e-3)    /* d'apres Lothenbach */
#define v_h2o      (18.e-3)
#define v_co2      (32.81e-3)
#define v_h2co3    (50.e-3)
#define v_hco3     (24.21.e-3)   /* d'apres Lothenbach */
#define v_co3      (-6.06e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.7e-3)    /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)
#define v_caoh     (26.20e-3)    /* a modifier */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (1.22e-7)    /* 1.22e-7 (radius = 1.75e-10 m) */
#define d_h        (4.76e-8)    /* 4.76e-8 (radius = 4.5e-10 m) */
#define d_co2      (1.43e-7)    /* 1.43e-7 (radius = 1.5e-10 m) */
#define d_h2co3    (7.2e-8)
#define d_hco3     (1.18e-7)    /* 1.07e-7 (radius = 2e-10 m) */
#define d_co3      (9.55e-8)    /* 9.53e-8 (radius = 2.25e-10 m) */
#define d_ca       (7.92e-8)
#define d_caoh     (1.07e-7)    /* (radius = 2e-10 m) */
#define d_cahco3   (1.07e-7)    /* (radius = 2e-10 m) */
#define d_h4sio4   (1.07e-7)    /*  */
#define d_h3sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_h2sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_cah2sio4 (1.07e-7)    /* (radius = 2e-10 m) */
#define d_cah3sio4 (1.07e-7)    /* (radius = 2e-10 m) */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */
#define K_henry    (1.238)           /* cste de Henry du CO2 / RT */

#define K_h2co3    (1.7e-3)          /* CO2[0] + H2O = H2CO3 */
#define K_hco3     (2.5e-4)          /* H2CO3   = HCO3[-] + H[+] */
#define K_co3      (4.69e-11)        /* HCO3[-] = CO3[2-] + H[+] */

#define K_h3sio4   (2.13675e13)      /* H2SiO4[2-] + H[+] = H3SiO4[-] */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-]  + H[+] = H4SiO4    */

#define K_caoh2    (1.)              /* Ca[2+] + 2OH[-]     = Ca(OH)2[0] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-]      = CaOH[+] */
#define K_cahco3   (1.276e+1)        /* Ca[2+] + HCO3[-]    = CaHCO3[+] */
#define K_caco3    (1.4e+3)          /* Ca[2+] + CO3[2-]    = CaCO3[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */

/*
  Solides
  CH  = Portlandite
  CC  = Calcite
  CSH = Hydrated Calcium Silicates
  Sam = Amorphous Silica
*/

/* volumes molaires solides (dm3/mole) */
#define v_ch       (33.e-3)      /* (33.e-3) */
#define v_cc       (37.e-3)      /* (37.e-3) */
#define v_csh      (V_TobI)
#define v_sam      (29.e-3)      /* d'apres Lothenbach */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* CH = Ca[2+] + 2OH[-] */
#define K_CC       (3.89e-9)         /* CC = Ca[2+] + CO3[2-]  */ 
#define K_Sam      (1.93642e-3)      /* S + 2H2O = H4SiO4 */
/* CxSyHz = xCa[2+] + 2xOH[-] + yH4SiO4 + (z-x-2y)H2O */
#define K_CSH      (K_TobI)
#define X_CSH      (X_TobI)
#define Y_CSH      (Y_TobI)

/* C-S-H */
/* Tobermorite I  (x,y,z) = (2,2.4,4.4) */
#define K_TobI     (5.53e-29)        /* Donnée GEMS */
#define X_TobI     (2)
#define Y_TobI     (2.4)
#define V_TobI     (118.e-3)         /* X_CSH*59.e-3 d'apres Lothenbach */
/* Tobermorite II (x,y,z) = (1.5,1.8,3.3) */
#define K_TobII    (6.42e-22)        /* Donnée GEMS */
#define X_TobII    (1.5)
#define Y_TobII    (1.8)
/* Jennite (x,y,z) = (1.5,0.9,2.4) */
#define K_Jen      (2.39e-16)        /* Donnée GEMS */
#define X_Jen      (1.5)
#define Y_Jen      (0.9)

#define C_CO2_eq   (K_h2o*K_h2o*K_CC/(K_h2co3*K_hco3*K_co3*K_CH))


/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,double,double) ;
static double poly4(double,double,double,double,double) ;
static double poly4_sansc(double,double,double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
#define MIN(a,b)  ((a < b) ? a : b)
#define MAX(a,b)  ((a > b) ? a : b)
#define NEGEXP(y) ((y < 0.) ? exp(y) : 1.)

#define Prod_CC(z_co2,z_chcc)  (K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc))
#define Prod_CH(z_co2,z_chcc)  (K_CH*MIN(z_co2,1.)*NEGEXP(z_chcc)/z_co2)
#define Prod_Sam(z_ch,z_csh)   (K_Sam*pow(MAX(z_ch,1.),-X_CSH/Y_CSH)*pow(NEGEXP(z_csh),1/Y_CSH))
#define Prod_CSH(P_CH,P_Sam)   (pow(P_CH,X_CSH)*pow(P_Sam,Y_CSH))

/* Parametres */
static double phi0,c_co2_eq,t_ch,t_csh,t_cc,t_sam,phi_min = -0.2 ;
static double n_ch_ref,n_csh_ref,n_sam_ref,n_cc_ref ;
static double P_CH_lim = pow(K_CSH/pow(K_Sam,Y_CSH),1./X_CSH) ;

static double k0_h2co3       = K_h2co3 ; 
/* CO2[0] + H2O              = H2CO3 */
static double k0_hco3        = K_hco3*K_h2co3/K_h2o ; 
/* CO2[0] + OH[-]            = HCO3[-] */
static double k0_co3         = K_co3*K_hco3*K_h2co3/(K_h2o*K_h2o) ; 
/* CO2[0] + 2OH[-]           = CO3[2-] + H2O */
/* Ca(OH)2[0]                = Ca[2+] + 2OH[-] */
static double k0_ca          = 1/(K_co3*K_hco3*K_h2co3*K_caco3) ;
/* CaCO3[0] + 2H[+]          = Ca[2+] + CO2[0] + H2O */
static double k0_caoh        = 0 ;
/* CaCO3[0] + H[+]           = CaOH[+] + CO2[0] */
static double k0_cahco3      = 0 ;
/* CaCO3[0] + H[+]           = CaHCO3[+] */
static double k0_h4sio4      = 0 ;
/* SiO2[0] + 2H2O            = H4SiO4[0] */
static double k0_h3sio4      = 0 ;
/* SiO2[0] + OH[-] + H2O     = H3SiO4[-] */
static double k0_h2sio4      = 0 ;
/* SiO2[0] + 2OH[-]          = H2SiO4[2-] */
static double k0_cah2sio4    = 0 ;
/* CaCO3[0] + SiO2[0] + H2O  = CaH2SiO4[0] + CO2[0] */
static double k0_cah3sio4    = 0 ;
/* CaCO3[0] + SiO2[0] + OH[-] = CaH3SiO4[-] + CO2[0] */

int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"N_CSH") == 0)    return (2) ;
  else if(strcmp(s,"N_CC") == 0)     return (3) ;
  else if(strcmp(s,"N_SAM") == 0)    return (4) ;
  else if(strcmp(s,"C_CO2_eq") == 0) return (5) ;
  else if(strcmp(s,"T_CH") == 0)     return (6) ;
  else if(strcmp(s,"T_CSH") == 0)    return (7) ;
  else if(strcmp(s,"T_CC") == 0)     return (8) ;
  else if(strcmp(s,"T_SAM") == 0)    return (9) ;
  else if(strcmp(s,"courbes") == 0)  return (10) ;
  else return(-1) ;
}


int dm53(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 11 ;
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_C],   "carbone") ;
  strcpy(mat->eqn[E_Ca],  "calcium") ;
  strcpy(mat->eqn[E_Si],  "Silicium") ;
  strcpy(mat->eqn[E_q],   "charge") ;

#if (U_CO2 == LOG_RHO)
  strcpy(mat->inc[I_CO2],  "logc_co2") ;
#else
  strcpy(mat->inc[I_CO2],  "c_co2") ;
#endif
  strcpy(mat->inc[I_CHCC], "z_chcc") ;
  strcpy(mat->inc[I_psi],  "psi") ;
  strcpy(mat->inc[I_CSH],  "z_csh") ;

  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_csh       = 600. ;
    double t_cc        = 0. ;
    double t_sam       = 0. ;
    double n_ch_ref    = 1. ;
    double n_csh_ref   = 1. ;
    double n_cc_ref    = 0. ;
    double n_sam_ref   = 0. ;

    mat->pr[pm("N_CH")]  = n_ch_ref ;
    mat->pr[pm("N_CSH")] = n_csh_ref ;
    mat->pr[pm("N_CC")]  = n_cc_ref ;
    mat->pr[pm("N_SAM")] = n_sam_ref ;
    mat->pr[pm("T_CH")]  = t_ch ;
    mat->pr[pm("T_CSH")] = t_csh ;
    mat->pr[pm("T_CC")]  = t_ch ;
    mat->pr[pm("T_SAM")] = t_sam ;

    dmat(mat,ficd,pm,n_donnees) ;

    n_ch_ref  = mat->pr[pm("N_CH")] ;
    n_csh_ref = mat->pr[pm("N_CSH")] ;
    n_cc_ref  = mat->pr[pm("N_CC")] ;
    n_sam_ref = mat->pr[pm("N_SAM")] ;

    t_ch      = mat->pr[pm("T_CH")] ;
    t_csh     = mat->pr[pm("T_CSH")] ;
    t_cc      = mat->pr[pm("T_CC")] ;
    t_sam     = mat->pr[pm("T_SAM")] ;

    if(n_cc_ref  == 0.) mat->pr[pm("N_CC")]  = n_ch_ref ;
    if(n_sam_ref == 0.) mat->pr[pm("N_SAM")] = Y_CSH*n_csh_ref ;
    if(t_cc  == 0.) mat->pr[pm("T_CC")]  = t_ch ;
    if(t_sam == 0.) mat->pr[pm("T_SAM")] = t_csh ;

    mat->pr[pm("C_CO2_eq")] = C_CO2_eq ;
  }
  
  return(mat->n) ;
}


int qm53(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  
  printf("Description:\n") ;

  printf("Carbonatation des m.c. satures en condition SC. Prise en compte d'un type de C-S-H par une approche discrete : Jennite, Tobermorite,...") ;
  
  printf("\n\n") ;
  printf("Le systeme est forme de 5 equations:\n") ;
#if (U_CO2 == LOG_RHO)
  printf("\t- la conservation de la masse de C      (logc_co2)\n") ;
#else
  printf("\t- la conservation de la masse de C      (c_co2)\n") ;
#endif
  printf("\t- la conservation de la charge          (psi)\n") ;
  printf("\t- la conservation de la masse de Ca     (z_chcc)\n") ;
  printf("\t- la conservation de la masse de Si     (z_csh)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n") ;

  printf("Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1    # Contenu en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"N_CSH = 2.4    # contenu en CSH (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5   # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"T_CSH = 1.e5   # Cinetique de dissolution des CSH (s)\n") ;

  return(NEQ) ;
}



void tb53(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}


void ch53(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in53(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;
  n_csh_ref = el.mat->pr[pm("N_CSH")] ;
  n_cc_ref  = el.mat->pr[pm("N_CC")] ;
  n_sam_ref = el.mat->pr[pm("N_SAM")] ;
  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double c_co2      = (C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;
    double z_co2      = c_co2/c_co2_eq ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ; 
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;
    
    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(z_chcc,0.) ;
    double n_cc_eq    = n_cc_ref*MAX(z_chcc,0.) ;
    double n_csh_eq   = n_csh_ref*MAX(z_csh,0.) ;
    double n_sam_eq   = n_sam_ref*MAX(z_csh,0.) ;
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch       = (z_co2 <= 1) ? n_ch_eq  : 0 ;
    double n_cc       = (z_co2 >  1) ? n_cc_eq  : 0 ;
    double n_csh      = (z_ch  >  1) ? n_csh_eq : 0 ;
    double n_sam      = (z_ch  <= 1) ? n_sam_eq : 0 ;

    /* solides */
    double n_si_s     = Y_CSH*n_csh + n_sam ;
    double n_ca_s     = n_ch + n_cc + X_CSH*n_csh ;

    /* porosite */
    double phi = phi0 ;

    /* contenus molaires */
    double n_co2      = phi*c_co2 ;
    double n_hco3     = phi*c_hco3 ;
    double n_h2co3    = phi*c_h2co3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_caoh     = phi*c_caoh ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h2sio4   = phi*c_h2sio4 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_cah2sio4 = phi*c_cah2sio4 ;

    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc ;
    N_Ca(i) = n_ca + n_caoh + n_cahco3 + n_cah2sio4 + n_cah3sio4 + n_ca_s ;
    N_Si(i) = n_h2sio4 + n_h3sio4 + n_h4sio4 + n_cah2sio4 + n_cah3sio4 + n_si_s ;

    /* densite de charge */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_caoh*c_caoh + z_cahco3*c_cahco3 + z_h2sio4*c_h2sio4 + z_h3sio4*c_h3sio4 + z_cah3sio4*c_cah3sio4 ;

#if (U_CO2 == LOG_RHO)
    u[i][I_CO2] = log(c_co2)/Ln10 ;
#else
    C_CO2(i)   = c_co2 ;
#endif
    Z_CSH(i)   = z_csh ;
    Z_CHCC(i)  = z_chcc ;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_CSH(i)   = n_csh ;
    N_Sam(i)   = n_sam ;

    N_CH0(i)   = n_ch ;
    N_CC0(i)   = n_cc ;
    N_CSH0(i)  = n_csh ;
    N_Sam0(i)  = n_sam ;
  }
  
  if(el.dim < dim) return ;

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex53(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  
  if(el.dim < dim) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct53(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;
  n_csh_ref = el.mat->pr[pm("N_CSH")] ;
  n_cc_ref  = el.mat->pr[pm("N_CC")] ;
  n_sam_ref = el.mat->pr[pm("N_SAM")] ;
  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;
  t_csh     = el.mat->pr[pm("T_CSH")] ;
  t_cc      = el.mat->pr[pm("T_CC")] ;
  t_sam     = el.mat->pr[pm("T_SAM")] ;
  
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double c_co2      = C_CO2(i) ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;

    /* cinetiques */
    double n_chn      = N_CHn(i) ;
    double n_ccn      = N_CCn(i) ;
    double n_cshn     = N_CSHn(i) ;
    double n_samn     = N_Samn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;  /* si z_co2 > 1 */
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;   /* si z_co2 < 1 */
    double n_csh_ci   = n_cshn*pow(z_ch,dt/t_csh) ;  /* si z_ch  < 1 */
    double n_sam_ci   = n_samn*pow(z_ch,-dt/t_sam) ; /* si z_ch  > 1 */
    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(z_chcc,0.) ;
    double n_cc_eq    = n_cc_ref*MAX(z_chcc,0.) ;
    double n_csh_eq   = n_csh_ref*MAX(z_csh,0.) ;
    double n_sam_eq   = n_sam_ref*MAX(z_csh,0.) ;
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch  = (z_co2 <= 1) ? n_ch_eq  : n_ch_ci ;
    double n_cc  = (z_co2 >  1) ? n_cc_eq  : n_cc_ci ;
    double n_csh = (z_ch  >  1) ? n_csh_eq : n_csh_ci ;
    double n_sam = (z_ch  <= 1) ? n_sam_eq : n_sam_ci ;

    /* Solides */
    double n_si_s     = Y_CSH*n_csh + n_sam ;
    double n_ca_s     = n_ch + n_cc + X_CSH*n_csh ;

    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_csh0     = N_CSH0(i) ;
    double n_sam0     = N_Sam0(i) ;
    double phi_th     = phi0 + v_ch*(n_ch0 - n_ch) + v_csh*(n_csh0 - n_csh) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
    double phi        = MAX(phi_th,phi_min) ;

    /* contenus molaires */
    double n_co2      = phi*c_co2 ;
    double n_h2co3    = phi*c_h2co3 ;
    double n_hco3     = phi*c_hco3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_caoh     = phi*c_caoh ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h2sio4   = phi*c_h2sio4 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_cah2sio4 = phi*c_cah2sio4 ;

    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc ;
    N_Ca(i) = n_ca + n_caoh + n_cahco3 + n_cah2sio4 + n_cah3sio4 + n_ca_s ;
    N_Si(i) = n_h2sio4 + n_h3sio4 + n_h4sio4 + n_cah2sio4 + n_cah3sio4 + n_si_s ;

    /* densite de charge */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_caoh*c_caoh + z_cahco3*c_cahco3 + z_h2sio4*c_h2sio4 + z_h3sio4*c_h3sio4 + z_cah3sio4*c_cah3sio4 ;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_CSH(i)   = n_csh ;
    N_Sam(i)   = n_sam ;

    if(c_oh < 0. || c_co2 < 0. || n_ca_s < 0. || n_si_s < 0. || phi < 0.) {
      printf("x         = %e\n",x[i][0]) ;
      printf("phi       = %e\n",phi) ;
      printf("c_co2     = %e\n",c_co2) ;
      printf("n_ch      = %e\n",n_ch) ;
      printf("n_cc      = %e\n",n_cc) ;
      printf("n_csh     = %e\n",n_csh) ;
      printf("n_sam     = %e\n",n_sam) ;
      printf("n_ca_s    = %e\n",n_ca_s) ;
      printf("n_si_s    = %e\n",n_si_s) ;
      printf("z_csh     = %e\n",z_csh) ;
      printf("z_chcc    = %e\n",z_chcc) ;
      printf("c_h2sio4  = %e\n",c_h2sio4) ;
      printf("c_h3sio4  = %e\n",c_h3sio4) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
  }
  
  if(el.dim < dim) return(0) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int mx53(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double c[2] ;

  double Dc_hSDc_co2[2]        ;
  double Dc_ohSDc_co2[2]       ;
  double Dc_h2co3SDc_co2[2]    ;
  double Dc_hco3SDc_co2[2]     ;
  double Dc_co3SDc_co2[2]      ;
  double Dc_caSDc_co2[2]       ;
  double Dc_caohSDc_co2[2]     ;
  double Dc_cahco3SDc_co2[2]   ;
  double Dc_h2sio4SDc_co2[2]   ;
  double Dc_h3sio4SDc_co2[2]   ;
  double Dc_h4sio4SDc_co2[2]   ;
  double Dc_cah2sio4SDc_co2[2] ;
  double Dc_cah3sio4SDc_co2[2] ;

  double Dc_hSDz_csh[2]        ;
  double Dc_ohSDz_csh[2]       ;
  double Dc_hco3SDz_csh[2]     ;
  double Dc_co3SDz_csh[2]      ;
  double Dc_caSDz_csh[2]       ;
  double Dc_caohSDz_csh[2]     ;
  double Dc_cahco3SDz_csh[2]   ;
  double Dc_h2sio4SDz_csh[2]   ;
  double Dc_h3sio4SDz_csh[2]   ;
  double Dc_h4sio4SDz_csh[2]   ;
  double Dc_cah2sio4SDz_csh[2] ;
  double Dc_cah3sio4SDz_csh[2] ;

  double Dc_hSDz_chcc[2]        ;
  double Dc_ohSDz_chcc[2]       ;
  double Dc_hco3SDz_chcc[2]     ;
  double Dc_co3SDz_chcc[2]      ;
  double Dc_caSDz_chcc[2]       ;
  double Dc_caohSDz_chcc[2]     ;
  double Dc_cahco3SDz_chcc[2]   ;
  double Dc_h2sio4SDz_chcc[2]   ;
  double Dc_h3sio4SDz_chcc[2]   ;
  double Dc_h4sio4SDz_chcc[2]   ;
  double Dc_cah2sio4SDz_chcc[2] ;
  double Dc_cah3sio4SDz_chcc[2] ;
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;
  n_csh_ref = el.mat->pr[pm("N_CSH")] ;
  n_cc_ref  = el.mat->pr[pm("N_CC")] ;
  n_sam_ref = el.mat->pr[pm("N_SAM")] ;
  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;
  t_csh     = el.mat->pr[pm("T_CSH")] ;
  t_cc      = el.mat->pr[pm("T_CC")] ;
  t_sam     = el.mat->pr[pm("T_SAM")] ;

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
    double c_co2      = C_CO2(i) ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;

    /* contenus en CH, CC, CSH, Sam */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_csh      = N_CSH(i) ;
    double n_sam      = N_Sam(i) ;

    /* cinetiques */
    double n_ccn      = N_CCn(i) ;
    double n_chn      = N_CHn(i) ;
    double n_cshn     = N_CSHn(i) ;
    double n_samn     = N_Samn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;
    double n_csh_ci   = n_cshn*pow(z_ch,dt/t_csh) ;
    double n_sam_ci   = n_samn*pow(z_ch,-dt/t_sam) ;

    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_csh0     = N_CSH0(i) ;
    double n_sam0     = N_Sam0(i) ;
    double phi_th     = phi0 + v_ch*(n_ch0 - n_ch) + v_csh*(n_csh0 - n_csh) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
    double phi        = MAX(phi_th,phi_min) ;

    /* derivees ... */
    /* ... par rapport a c_co2 */
    double dc_co2             = 1.e-4*c_co2*((c_co2 > C_CO2n(i)) ? 1 : -1) ;
    double c_co22             = c_co2 + dc_co2 ;
    double c_oh2              = concentration_oh(c_co22,z_csh,z_chcc) ;
    double dc_ohsdc_co2       = (c_oh2 - c_oh)/dc_co2 ;
    double dc_hsdc_co2        = - c_h*dc_ohsdc_co2/c_oh ;
    double dc_h2co3sdc_co2    = K_h2co3 ;
    double dc_hco3sdc_co2     = c_hco3*(dc_h2co3sdc_co2/c_h2co3	\
			      + dc_ohsdc_co2/c_oh) ;
    double dc_co3sdc_co2      = c_co3*(dc_hco3sdc_co2/c_hco3 \
			      + dc_ohsdc_co2/c_oh) ;

    double dP_CCsdc_co2       = (z_co2 < 1.) ? P_CC/c_co2 : 0. ;
    double dc_casdc_co2       = (dP_CCsdc_co2 - c_ca*dc_co3sdc_co2)/c_co3 ;
    double dc_caohsdc_co2     = K_caoh*(dc_casdc_co2*c_oh \
			      + c_ca*dc_ohsdc_co2) ;
    double dc_cahco3sdc_co2   = K_cahco3*(dc_casdc_co2*c_hco3 \
			      + c_ca*dc_hco3sdc_co2) ;

    double dz_chsdc_co2       = z_ch*(dP_CCsdc_co2/P_CC - 1./c_co2) ;
    double dP_Samsdc_co2      = (z_ch > 1.) ? P_Sam*(-X_CSH/Y_CSH)*dz_chsdc_co2/z_ch : 0 ;

    double dc_h4sio4sdc_co2   = dP_Samsdc_co2 ;
    double dc_h3sio4sdc_co2   = c_h3sio4*(dc_h4sio4sdc_co2/c_h4sio4 \
			      - dc_hsdc_co2/c_h) ;
    double dc_h2sio4sdc_co2   = c_h2sio4*(dc_h3sio4sdc_co2/c_h3sio4 \
			      - dc_hsdc_co2/c_h) ;
    double dc_cah3sio4sdc_co2 = K_cah3sio4*(dc_casdc_co2*c_h3sio4 \
			      + c_ca*dc_h3sio4sdc_co2) ;
    double dc_cah2sio4sdc_co2 = K_cah2sio4*(dc_casdc_co2*c_h2sio4 \
			      + c_ca*dc_h2sio4sdc_co2) ;

    double dn_ch_cisdc_co2   = n_ch_ci*(-dt/t_ch)/c_co2 ;
    double dn_cc_cisdc_co2   = n_cc_ci*(dt/t_cc)/c_co2 ;
    double dn_csh_cisdc_co2  = n_csh_ci*(dt/t_csh)*dz_chsdc_co2/z_ch ;
    double dn_sam_cisdc_co2  = n_sam_ci*(-dt/t_sam)*dz_chsdc_co2/z_ch ;

    double dn_ch_eqsdc_co2   = 0 ;
    double dn_cc_eqsdc_co2   = 0 ;
    double dn_csh_eqsdc_co2  = 0 ;
    double dn_sam_eqsdc_co2  = 0 ;

    double dn_chsdc_co2  = (z_co2 <= 1) ? dn_ch_eqsdc_co2  : dn_ch_cisdc_co2 ;
    double dn_ccsdc_co2  = (z_co2 >  1) ? dn_cc_eqsdc_co2  : dn_cc_cisdc_co2 ;
    double dn_cshsdc_co2 = (z_ch  >  1) ? dn_csh_eqsdc_co2 : dn_csh_cisdc_co2 ;
    double dn_samsdc_co2 = (z_ch  <= 1) ? dn_sam_eqsdc_co2 : dn_sam_cisdc_co2 ;

    double dn_ca_ssdc_co2 = dn_chsdc_co2 + dn_ccsdc_co2 + X_CSH*dn_cshsdc_co2 ;
    double dn_si_ssdc_co2 = Y_CSH*dn_cshsdc_co2 + dn_samsdc_co2 ;

    double dphi_thsdc_co2     = - v_ch*dn_chsdc_co2 - v_csh*dn_cshsdc_co2 - v_sam*dn_samsdc_co2 - v_cc*dn_ccsdc_co2 ;
    double dphisdc_co2        = (phi_th < phi_min) ? 0 : dphi_thsdc_co2 ;

    double dn_co2sdc_co2      = phi + dphisdc_co2*c_co2 ;
    double dn_h2co3sdc_co2    = phi*dc_h2co3sdc_co2 + dphisdc_co2*c_h2co3 ;
    double dn_hco3sdc_co2     = phi*dc_hco3sdc_co2 + dphisdc_co2*c_hco3 ;
    double dn_co3sdc_co2      = phi*dc_co3sdc_co2 + dphisdc_co2*c_co3 ;
    double dn_casdc_co2       = phi*dc_casdc_co2 + dphisdc_co2*c_ca ;
    double dn_caohsdc_co2     = phi*dc_caohsdc_co2 + dphisdc_co2*c_caoh ;
    double dn_cahco3sdc_co2   = phi*dc_cahco3sdc_co2 + dphisdc_co2*c_cahco3 ;
    double dn_cah3sio4sdc_co2 = phi*dc_cah3sio4sdc_co2 + dphisdc_co2*c_cah3sio4 ;
    double dn_cah2sio4sdc_co2 = phi*dc_cah2sio4sdc_co2 + dphisdc_co2*c_cah2sio4 ;
    double dn_h2sio4sdc_co2   = phi*dc_h2sio4sdc_co2 + dphisdc_co2*c_h2sio4 ;
    double dn_h3sio4sdc_co2   = phi*dc_h3sio4sdc_co2 + dphisdc_co2*c_h3sio4 ;
    double dn_h4sio4sdc_co2   = phi*dc_h4sio4sdc_co2 + dphisdc_co2*c_h4sio4 ;

    /* par rapport a z_csh */
    /* double dz_csh             = ((z_csh > 0.) ? 1 : -1)*1.e-2 ; */
    double dz_csh             = 1.e-6*((z_csh > Z_CSHn(i)) ? 1 : -1) ;
    double z_csh2             = z_csh + dz_csh ;
    double c_oh1              = concentration_oh(c_co2,z_csh2,z_chcc) ;
    double dc_ohsdz_csh       = (c_oh1 - c_oh)/dz_csh ;
    double dc_hsdz_csh        = - c_h*dc_ohsdz_csh/c_oh ;
    double dc_hco3sdz_csh     = c_hco3*(dc_ohsdz_csh/c_oh) ;
    double dc_co3sdz_csh      = c_co3*(dc_hco3sdz_csh/c_hco3 \
			      + dc_ohsdz_csh/c_oh) ;
    double dc_casdz_csh       = - c_ca*dc_co3sdz_csh/c_co3 ;
    double dc_caohsdz_csh     = K_caoh*(dc_casdz_csh*c_oh \
			      + c_ca*dc_ohsdz_csh) ;
    double dc_cahco3sdz_csh   = K_cahco3*(dc_casdz_csh*c_hco3 \
			      + c_ca*dc_hco3sdz_csh) ;

    double dP_Samsdz_csh      = (z_csh < 0) ? P_Sam/Y_CSH : 0 ;

    double dc_h4sio4sdz_csh   = dP_Samsdz_csh ;
    double dc_h3sio4sdz_csh   = c_h3sio4*(dc_h4sio4sdz_csh/c_h4sio4 \
			      - dc_hsdz_csh/c_h) ;
    double dc_h2sio4sdz_csh   = c_h2sio4*(dc_h3sio4sdz_csh/c_h3sio4 \
			      - dc_hsdz_csh/c_h) ;
    double dc_cah3sio4sdz_csh = K_cah3sio4*(dc_casdz_csh*c_h3sio4 \
			      + c_ca*dc_h3sio4sdz_csh) ;
    double dc_cah2sio4sdz_csh = K_cah2sio4*(dc_casdz_csh*c_h2sio4 \
			      + c_ca*dc_h2sio4sdz_csh) ;

    double dn_ch_cisdz_csh   = 0 ;
    double dn_cc_cisdz_csh   = 0 ;
    double dn_csh_cisdz_csh  = 0 ;
    double dn_sam_cisdz_csh  = 0 ;

    double dn_ch_eqsdz_csh   = 0 ;
    double dn_cc_eqsdz_csh   = 0 ;
    double dn_csh_eqsdz_csh  = (z_csh > 0) ? n_csh_ref : 0 ;
    double dn_sam_eqsdz_csh  = (z_csh > 0) ? n_sam_ref : 0 ;

    double dn_chsdz_csh  = (z_co2 <= 1) ? dn_ch_eqsdz_csh  : dn_ch_cisdz_csh ;
    double dn_ccsdz_csh  = (z_co2 >  1) ? dn_cc_eqsdz_csh  : dn_cc_cisdz_csh ;
    double dn_cshsdz_csh = (z_ch  >  1) ? dn_csh_eqsdz_csh : dn_csh_cisdz_csh ;
    double dn_samsdz_csh = (z_ch  <= 1) ? dn_sam_eqsdz_csh : dn_sam_cisdz_csh ;

    double dn_ca_ssdz_csh = dn_chsdz_csh + dn_ccsdz_csh + X_CSH*dn_cshsdz_csh ;
    double dn_si_ssdz_csh = Y_CSH*dn_cshsdz_csh + dn_samsdz_csh ;

    double dphi_thsdz_csh = - v_ch*dn_chsdz_csh - v_cc*dn_ccsdz_csh - v_csh*dn_cshsdz_csh - v_sam*dn_samsdz_csh ;
    double dphisdz_csh    = (phi_th < phi_min) ? 0 : dphi_thsdz_csh ;

    double dn_co2sdz_csh      = dphisdz_csh*c_co2 ;
    double dn_h2co3sdz_csh    = dphisdz_csh*c_h2co3 ;
    double dn_hco3sdz_csh     = phi*dc_hco3sdz_csh + dphisdz_csh*c_hco3 ;
    double dn_co3sdz_csh      = phi*dc_co3sdz_csh + dphisdz_csh*c_co3 ;
    double dn_casdz_csh       = phi*dc_casdz_csh + dphisdz_csh*c_ca ;
    double dn_caohsdz_csh     = phi*dc_caohsdz_csh + dphisdz_csh*c_caoh ;
    double dn_cahco3sdz_csh   = phi*dc_cahco3sdz_csh + dphisdz_csh*c_cahco3 ;
    double dn_cah3sio4sdz_csh = phi*dc_cah3sio4sdz_csh + dphisdz_csh*c_cah3sio4 ;
    double dn_cah2sio4sdz_csh = phi*dc_cah2sio4sdz_csh + dphisdz_csh*c_cah2sio4 ;
    double dn_h2sio4sdz_csh   = phi*dc_h2sio4sdz_csh + dphisdz_csh*c_h2sio4 ;
    double dn_h3sio4sdz_csh   = phi*dc_h3sio4sdz_csh + dphisdz_csh*c_h3sio4 ;
    double dn_h4sio4sdz_csh   = phi*dc_h4sio4sdz_csh + dphisdz_csh*c_h4sio4 ;

    /* par rapport a z_chcc */
    /* double dz_chcc             = ((z_chcc > 0.) ? 1 : -1)*1.e-2 ; */
    double dz_chcc             = 1.e-6*((z_chcc > Z_CHCCn(i)) ? 1 : -1) ;
    double z_chcc2             = z_chcc + dz_chcc ;
    double c_oh3               = concentration_oh(c_co2,z_csh,z_chcc2) ;
    double dc_ohsdz_chcc       = (c_oh3 - c_oh)/dz_chcc ;
    double dc_hsdz_chcc        = - c_h*dc_ohsdz_chcc/c_oh ;
    double dc_hco3sdz_chcc     = c_hco3*(dc_ohsdz_chcc/c_oh) ;
    double dc_co3sdz_chcc      = c_co3*(dc_hco3sdz_chcc/c_hco3 \
			       + dc_ohsdz_chcc/c_oh) ;

    double dP_CCsdz_chcc       = (z_chcc < 0) ? P_CC : 0 ;
    double dc_casdz_chcc       = (dP_CCsdz_chcc - c_ca*dc_co3sdz_chcc)/c_co3 ;
    double dc_caohsdz_chcc     = K_caoh*(dc_casdz_chcc*c_oh \
			       + c_ca*dc_ohsdz_chcc) ;
    double dc_cahco3sdz_chcc   = K_cahco3*(dc_casdz_chcc*c_hco3 \
			       + c_ca*dc_hco3sdz_chcc) ;

    double dz_chsdz_chcc  = z_ch*dP_CCsdz_chcc/P_CC ;
    double dP_Samsdz_chcc = (z_ch > 1) ? P_Sam*(-X_CSH/Y_CSH)*dP_CCsdz_chcc/P_CC : 0 ;

    double dc_h4sio4sdz_chcc   = dP_Samsdz_chcc ;
    double dc_h3sio4sdz_chcc   = c_h3sio4*(dc_h4sio4sdz_chcc/c_h4sio4 \
			       - dc_hsdz_chcc/c_h) ;
    double dc_h2sio4sdz_chcc   = c_h2sio4*(dc_h3sio4sdz_chcc/c_h3sio4 \
			       - dc_hsdz_chcc/c_h) ;
    double dc_cah3sio4sdz_chcc = K_cah3sio4*(dc_casdz_chcc*c_h3sio4 \
			       + c_ca*dc_h3sio4sdz_chcc) ;
    double dc_cah2sio4sdz_chcc = K_cah2sio4*(dc_casdz_chcc*c_h2sio4 \
			       + c_ca*dc_h2sio4sdz_chcc) ;

    double dn_ch_cisdz_chcc   = 0 ;
    double dn_cc_cisdz_chcc   = 0 ;
    double dn_csh_cisdz_chcc  = n_csh_ci*(dt/t_csh)*dz_chsdz_chcc/z_ch ;
    double dn_sam_cisdz_chcc  = n_sam_ci*(-dt/t_sam)*dz_chsdz_chcc/z_ch ;

    double dn_ch_eqsdz_chcc   = (z_chcc > 0) ? n_ch_ref : 0 ;
    double dn_cc_eqsdz_chcc   = (z_chcc > 0) ? n_cc_ref : 0 ;
    double dn_csh_eqsdz_chcc  = 0 ;
    double dn_sam_eqsdz_chcc  = 0 ;

    double dn_chsdz_chcc  = (z_co2 <= 1) ? dn_ch_eqsdz_chcc  : dn_ch_cisdz_chcc ;
    double dn_ccsdz_chcc  = (z_co2 >  1) ? dn_cc_eqsdz_chcc  : dn_cc_cisdz_chcc ;
    double dn_cshsdz_chcc = (z_ch  >  1) ? dn_csh_eqsdz_chcc : dn_csh_cisdz_chcc ;
    double dn_samsdz_chcc = (z_ch  <= 1) ? dn_sam_eqsdz_chcc : dn_sam_cisdz_chcc ;

    double dn_si_ssdz_chcc = Y_CSH*dn_cshsdz_chcc + dn_samsdz_chcc ;
    double dn_ca_ssdz_chcc = dn_chsdz_chcc + dn_ccsdz_chcc + X_CSH*dn_cshsdz_chcc ;

    double dphi_thsdz_chcc = - v_ch*dn_chsdz_chcc - v_cc*dn_ccsdz_chcc - v_csh*dn_cshsdz_chcc - v_sam*dn_samsdz_chcc ;
    double dphisdz_chcc    = (phi_th < phi_min) ? 0 : dphi_thsdz_chcc ;

    double dn_co2sdz_chcc      = dphisdz_chcc*c_co2 ;
    double dn_h2co3sdz_chcc    = dphisdz_chcc*c_h2co3 ;
    double dn_hco3sdz_chcc     = phi*dc_hco3sdz_chcc + dphisdz_chcc*c_hco3 ;
    double dn_co3sdz_chcc      = phi*dc_co3sdz_chcc + dphisdz_chcc*c_co3 ;
    double dn_casdz_chcc       = phi*dc_casdz_chcc + dphisdz_chcc*c_ca ;
    double dn_caohsdz_chcc     = phi*dc_caohsdz_chcc + dphisdz_chcc*c_caoh ;
    double dn_cahco3sdz_chcc   = phi*dc_cahco3sdz_chcc + dphisdz_chcc*c_cahco3 ;
    double dn_cah3sio4sdz_chcc = phi*dc_cah3sio4sdz_chcc + dphisdz_chcc*c_cah3sio4 ;
    double dn_cah2sio4sdz_chcc = phi*dc_cah2sio4sdz_chcc + dphisdz_chcc*c_cah2sio4 ;
    double dn_h2sio4sdz_chcc   = phi*dc_h2sio4sdz_chcc + dphisdz_chcc*c_h2sio4 ;
    double dn_h3sio4sdz_chcc   = phi*dc_h3sio4sdz_chcc + dphisdz_chcc*c_h3sio4 ;
    double dn_h4sio4sdz_chcc   = phi*dc_h4sio4sdz_chcc + dphisdz_chcc*c_h4sio4 ;


    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_CO2+j)   += volume[i]*(dn_co2sdc_co2 + dn_h2co3sdc_co2 + dn_hco3sdc_co2 + dn_co3sdc_co2 + dn_cahco3sdc_co2) ;
    K(E_C+j,I_CHCC+j)  += volume[i]*(dn_co2sdz_chcc + dn_h2co3sdz_chcc + dn_hco3sdz_chcc + dn_co3sdz_chcc + dn_cahco3sdz_chcc + dn_ccsdz_chcc) ;
    K(E_C+j,I_CSH+j)   += volume[i]*(dn_co2sdz_csh + dn_h2co3sdz_csh + dn_hco3sdz_csh + dn_co3sdz_csh + dn_cahco3sdz_csh) ;
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_CO2+j)  += volume[i]*(dn_casdc_co2 + dn_caohsdc_co2 + dn_cahco3sdc_co2 + dn_cah2sio4sdc_co2 + dn_cah3sio4sdc_co2 + dn_ca_ssdc_co2) ;
    K(E_Ca+j,I_CHCC+j) += volume[i]*(dn_casdz_chcc + dn_caohsdz_chcc + dn_cahco3sdz_chcc + dn_cah2sio4sdz_chcc + dn_cah3sio4sdz_chcc + dn_ca_ssdz_chcc) ;
    K(E_Ca+j,I_CSH+j)  += volume[i]*(dn_casdz_csh + dn_caohsdz_csh + dn_cahco3sdz_csh + dn_cah2sio4sdz_csh + dn_cah3sio4sdz_csh + dn_ca_ssdz_csh) ;
    /*
      Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
    */
    K(E_Si+j,I_CO2+j)  += volume[i]*(dn_h2sio4sdc_co2 + dn_h3sio4sdc_co2 + dn_h4sio4sdc_co2 + dn_cah2sio4sdc_co2 + dn_cah3sio4sdc_co2 + dn_si_ssdc_co2) ;
    K(E_Si+j,I_CHCC+j) += volume[i]*(dn_h2sio4sdz_chcc + dn_h3sio4sdz_chcc + dn_h4sio4sdz_chcc + dn_cah2sio4sdz_chcc + dn_cah3sio4sdz_chcc + dn_si_ssdz_chcc) ;
    K(E_Si+j,I_CSH+j)  += volume[i]*(dn_h2sio4sdz_csh + dn_h3sio4sdz_csh + dn_h4sio4sdz_csh + dn_cah2sio4sdz_csh + dn_cah3sio4sdz_csh + dn_si_ssdz_csh) ;
    /*
      Conservation de la charge  : div(w_q) = 0
    */


    /* sauvegardes pour les termes de transport */
    Dc_hSDc_co2[i]        = dc_hsdc_co2 ;
    Dc_ohSDc_co2[i]       = dc_ohsdc_co2 ;
    Dc_h2co3SDc_co2[i]    = dc_h2co3sdc_co2 ;
    Dc_hco3SDc_co2[i]     = dc_hco3sdc_co2 ;
    Dc_co3SDc_co2[i]      = dc_co3sdc_co2 ;
    Dc_caSDc_co2[i]       = dc_casdc_co2 ;
    Dc_caohSDc_co2[i]     = dc_caohsdc_co2 ;
    Dc_cahco3SDc_co2[i]   = dc_cahco3sdc_co2 ;
    Dc_h2sio4SDc_co2[i]   = dc_h2sio4sdc_co2 ;
    Dc_h3sio4SDc_co2[i]   = dc_h3sio4sdc_co2 ;
    Dc_h4sio4SDc_co2[i]   = dc_h4sio4sdc_co2 ;
    Dc_cah3sio4SDc_co2[i] = dc_cah3sio4sdc_co2 ;
    Dc_cah2sio4SDc_co2[i] = dc_cah2sio4sdc_co2 ;
    
    Dc_hSDz_csh[i]        = dc_hsdz_csh ;
    Dc_ohSDz_csh[i]       = dc_ohsdz_csh ;
    Dc_hco3SDz_csh[i]     = dc_hco3sdz_csh ;
    Dc_co3SDz_csh[i]      = dc_co3sdz_csh ;
    Dc_caSDz_csh[i]       = dc_casdz_csh ;
    Dc_caohSDz_csh[i]     = dc_caohsdz_csh ;
    Dc_cahco3SDz_csh[i]   = dc_cahco3sdz_csh ;
    Dc_h2sio4SDz_csh[i]   = dc_h2sio4sdz_csh ;
    Dc_h3sio4SDz_csh[i]   = dc_h3sio4sdz_csh ;
    Dc_h4sio4SDz_csh[i]   = dc_h4sio4sdz_csh ;
    Dc_cah3sio4SDz_csh[i] = dc_cah3sio4sdz_csh ;
    Dc_cah2sio4SDz_csh[i] = dc_cah2sio4sdz_csh ;
    
    Dc_hSDz_chcc[i]        = dc_hsdz_chcc ;
    Dc_ohSDz_chcc[i]       = dc_ohsdz_chcc ;
    Dc_hco3SDz_chcc[i]     = dc_hco3sdz_chcc ;
    Dc_co3SDz_chcc[i]      = dc_co3sdz_chcc ;
    Dc_caSDz_chcc[i]       = dc_casdz_chcc ;
    Dc_caohSDz_chcc[i]     = dc_caohsdz_chcc ;
    Dc_cahco3SDz_chcc[i]   = dc_cahco3sdz_chcc ;
    Dc_h2sio4SDz_chcc[i]   = dc_h2sio4sdz_chcc ;
    Dc_h3sio4SDz_chcc[i]   = dc_h3sio4sdz_chcc ;
    Dc_h4sio4SDz_chcc[i]   = dc_h4sio4sdz_chcc ;
    Dc_cah3sio4SDz_chcc[i] = dc_cah3sio4sdz_chcc ;
    Dc_cah2sio4SDz_chcc[i] = dc_cah2sio4sdz_chcc ;
  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trf_oh       = tr*KF_OH ;
  double trf_h        = tr*KF_H ;
  double trf_co2      = tr*KF_CO2 ;
  double trf_h2co3    = tr*KF_H2CO3 ;
  double trf_hco3     = tr*KF_HCO3 ;
  double trf_co3      = tr*KF_CO3 ;
  double trf_ca       = tr*KF_Ca ;
  double trf_caoh     = tr*KF_CaOH ;
  double trf_cahco3   = tr*KF_CaHCO3 ;
  double trf_cah2sio4 = tr*KF_CaH2SiO4 ;
  double trf_cah3sio4 = tr*KF_CaH3SiO4 ;
  double trf_h2sio4   = tr*KF_H2SiO4 ;
  double trf_h3sio4   = tr*KF_H3SiO4 ;
  double trf_h4sio4   = tr*KF_H4SiO4 ;

  double tre_hco3     = tr*Kpsi_HCO3 ;
  double tre_co3      = tr*Kpsi_CO3 ;
  double tre_ca       = tr*Kpsi_Ca ;
  double tre_caoh     = tr*Kpsi_CaOH ;
  double tre_cahco3   = tr*Kpsi_CaHCO3 ;
  double tre_cah3sio4 = tr*Kpsi_CaH3SiO4 ;
  double tre_h2sio4   = tr*Kpsi_H2SiO4 ;
  double tre_h3sio4   = tr*Kpsi_H3SiO4 ;

  double tre_q        = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_co2 + trf_h2co3*Dc_h2co3SDc_co2[i] + trf_hco3*Dc_hco3SDc_co2[i] + trf_co3*Dc_co3SDc_co2[i] + trf_cahco3*Dc_cahco3SDc_co2[i] ;
  }
  K(E_C,I_CO2)          += + c[0] ;
  K(E_C,I_CO2+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CO2)      += - c[0] ;
  K(E_C+NEQ,I_CO2+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDz_csh[i] + trf_co3*Dc_co3SDz_csh[i] + trf_cahco3*Dc_cahco3SDz_csh[i] ;
  }
  K(E_C,I_CSH)          += + c[0] ;
  K(E_C,I_CSH+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CSH)      += - c[0] ;
  K(E_C+NEQ,I_CSH+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDz_chcc[i] + trf_co3*Dc_co3SDz_chcc[i] + trf_cahco3*Dc_cahco3SDz_chcc[i] ;
  }
  K(E_C,I_CHCC)          += + c[0] ;
  K(E_C,I_CHCC+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CHCC)      += - c[0] ;
  K(E_C+NEQ,I_CHCC+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_cahco3 ;
  }
  K(E_C,I_psi)          += + c[0] ;
  K(E_C,I_psi+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_psi)      += - c[0] ;
  K(E_C+NEQ,I_psi+NEQ)  += + c[1] ;

  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_co2[i] + trf_caoh*Dc_caohSDc_co2[i] + trf_cahco3*Dc_cahco3SDc_co2[i] + trf_cah2sio4*Dc_cah2sio4SDc_co2[i] + trf_cah3sio4*Dc_cah3sio4SDc_co2[i] ;
  }
  K(E_Ca,I_CO2)         += + c[0] ;
  K(E_Ca,I_CO2+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CO2)     += - c[0] ;
  K(E_Ca+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDz_csh[i] + trf_caoh*Dc_caohSDz_csh[i] + trf_cahco3*Dc_cahco3SDz_csh[i] + trf_cah2sio4*Dc_cah2sio4SDz_csh[i] + trf_cah3sio4*Dc_cah3sio4SDz_csh[i] ;
  }
  K(E_Ca,I_CSH)         += + c[0] ;
  K(E_Ca,I_CSH+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CSH)     += - c[0] ;
  K(E_Ca+NEQ,I_CSH+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDz_chcc[i] + trf_caoh*Dc_caohSDz_chcc[i] + trf_cahco3*Dc_cahco3SDz_chcc[i] + trf_cah2sio4*Dc_cah2sio4SDz_chcc[i] + trf_cah3sio4*Dc_cah3sio4SDz_chcc[i] ;
  }
  K(E_Ca,I_CHCC)         += + c[0] ;
  K(E_Ca,I_CHCC+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CHCC)     += - c[0] ;
  K(E_Ca+NEQ,I_CHCC+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_caoh + tre_cahco3 + tre_cah3sio4 ;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h2sio4*Dc_h2sio4SDc_co2[i] + trf_h3sio4*Dc_h3sio4SDc_co2[i] + trf_h4sio4*Dc_h4sio4SDc_co2[i] + trf_cah2sio4*Dc_cah2sio4SDc_co2[i] + trf_cah3sio4*Dc_cah3sio4SDc_co2[i] ;
  }
  K(E_Si,I_CO2)         += + c[0] ;
  K(E_Si,I_CO2+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CO2)     += - c[0] ;
  K(E_Si+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h2sio4*Dc_h2sio4SDz_csh[i] + trf_h3sio4*Dc_h3sio4SDz_csh[i] + trf_h4sio4*Dc_h4sio4SDz_csh[i] + trf_cah2sio4*Dc_cah2sio4SDz_csh[i] + trf_cah3sio4*Dc_cah3sio4SDz_csh[i] ;
  }
  K(E_Si,I_CSH)         += + c[0] ;
  K(E_Si,I_CSH+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CSH)     += - c[0] ;
  K(E_Si+NEQ,I_CSH+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h2sio4*Dc_h2sio4SDz_chcc[i] + trf_h3sio4*Dc_h3sio4SDz_chcc[i] + trf_h4sio4*Dc_h4sio4SDz_chcc[i] + trf_cah2sio4*Dc_cah2sio4SDz_chcc[i] + trf_cah3sio4*Dc_cah3sio4SDz_chcc[i] ;
  }
  K(E_Si,I_CHCC)         += + c[0] ;
  K(E_Si,I_CHCC+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CHCC)     += - c[0] ;
  K(E_Si+NEQ,I_CHCC+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_h2sio4 + tre_h3sio4 + tre_cah3sio4 ;
  }
  K(E_Si,I_psi)          += + c[0] ;
  K(E_Si,I_psi+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_psi)      += - c[0] ;
  K(E_Si+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDc_co2[i] + z_oh*trf_oh*Dc_ohSDc_co2[i] + z_hco3*trf_hco3*Dc_hco3SDc_co2[i] + z_co3*trf_co3*Dc_co3SDc_co2[i] + z_ca*trf_ca*Dc_caSDc_co2[i] + z_caoh*trf_caoh*Dc_caohSDc_co2[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_co2[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_co2[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_co2[i] + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_co2[i] ;
  }
  K(E_q,I_CO2)           += + c[0] ;
  K(E_q,I_CO2+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CO2)       += - c[0] ;
  K(E_q+NEQ,I_CO2+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDz_csh[i] + z_oh*trf_oh*Dc_ohSDz_csh[i] + z_hco3*trf_hco3*Dc_hco3SDz_csh[i] + z_co3*trf_co3*Dc_co3SDz_csh[i] + z_ca*trf_ca*Dc_caSDz_csh[i] + z_caoh*trf_caoh*Dc_caohSDz_csh[i] + z_cahco3*trf_cahco3*Dc_cahco3SDz_csh[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDz_csh[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDz_csh[i] + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDz_csh[i] ;
  }
  K(E_q,I_CSH)           += + c[0] ;
  K(E_q,I_CSH+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CSH)       += - c[0] ;
  K(E_q+NEQ,I_CSH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDz_chcc[i] + z_oh*trf_oh*Dc_ohSDz_chcc[i] + z_hco3*trf_hco3*Dc_hco3SDz_chcc[i] + z_co3*trf_co3*Dc_co3SDz_chcc[i] + z_ca*trf_ca*Dc_caSDz_chcc[i] + z_caoh*trf_caoh*Dc_caohSDz_chcc[i] + z_cahco3*trf_cahco3*Dc_cahco3SDz_chcc[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDz_chcc[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDz_chcc[i] + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDz_chcc[i] ;
  }
  K(E_q,I_CHCC)           += + c[0] ;
  K(E_q,I_CHCC+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CHCC)       += - c[0] ;
  K(E_q+NEQ,I_CHCC+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  }


#if (U_CO2 == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_CO2)     *= Ln10*C_CO2(0) ;
    K(i,I_CO2+NEQ) *= Ln10*C_CO2(1) ;
  }
#endif
  

  return(0) ;

#undef K
}


void rs53(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;

#undef R
}


int so53(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
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
  c_co2_eq = el.mat->pr[pm("C_CO2_eq")] ;

  /* initialisation */
  nso = 26 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees */
  {
    /* molarites */
#if (U_CO2 == LOG_RHO)
    double c_co2      = exp(Ln10*param(u,h_s,el.nn,I_CO2)) ;
#else
    double c_co2      = param(u,h_s,el.nn,I_CO2) ;
#endif
    double z_co2      = c_co2/c_co2_eq ;
    double z_csh      = param(u,h_s,el.nn,I_CSH) ;
    double z_chcc     = param(u,h_s,el.nn,I_CHCC) ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;

    /* densite de charge */
    double c_q = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_caoh*c_caoh + z_cahco3*c_cahco3 + z_h2sio4*c_h2sio4 + z_h3sio4*c_h3sio4 + z_cah3sio4*c_cah3sio4 ;
  /* contenus solides */
    int j = (el.dim < dim) ? 0 : ((s[0] < (x[0][0] + x[1][0])*0.5) ? 0 : 1) ;
    double n_ch       = N_CH(j) ;
    double n_cc       = N_CC(j) ;
    double n_csh      = N_CSH(j)   ;
    double n_sam      = N_Sam(j)   ;
  /* porosite */
    double n_ch0      = N_CH0(j) ;
    double n_cc0      = N_CC0(j) ;
    double n_csh0     = N_CSH0(j) ;
    double n_sam0     = N_Sam0(j) ;
    double phi_th     = phi0 + v_ch*(n_ch0 - n_ch) + v_csh*(n_csh0 - n_csh) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
    double phi        = MAX(phi_th,phi_min) ;
    double psi        = param(u,h_s,el.nn,I_psi) ;

    i = 0 ;
    strcpy(r[i].text,"c_co2") ; r[i].n = 1 ;
    r[i++].v[0] = c_co2 ;
    strcpy(r[i].text,"ph") ; r[i].n = 1 ;
    r[i++].v[0] = 14 + log(c_oh)/log(10.) ;
    strcpy(r[i].text,"c_h2co3") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2co3 ;
    strcpy(r[i].text,"c_hco3") ; r[i].n = 1 ;
    r[i++].v[0] = c_hco3 ;
    strcpy(r[i].text,"c_co3") ; r[i].n = 1 ;
    r[i++].v[0] = c_co3 ;
    strcpy(r[i].text,"c_ca") ; r[i].n = 1 ;
    r[i++].v[0] = c_ca ;
    strcpy(r[i].text,"c_cahco3") ; r[i].n = 1 ;
    r[i++].v[0] = c_cahco3 ;
    strcpy(r[i].text,"c_cah3sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_cah3sio4 ;
    strcpy(r[i].text,"c_h3sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h3sio4 ;
    strcpy(r[i].text,"c_h4sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h4sio4 ;
    strcpy(r[i].text,"n_ch") ; r[i].n = 1 ;
    r[i++].v[0] = n_ch ;
    strcpy(r[i].text,"n_cc") ; r[i].n = 1 ;
    r[i++].v[0] = n_cc ;
    strcpy(r[i].text,"n_csh") ; r[i].n = 1 ;
    r[i++].v[0] = n_csh ;
    strcpy(r[i].text,"n_sam") ; r[i].n = 1 ;
    r[i++].v[0] = n_sam ;
    strcpy(r[i].text,"porosite") ; r[i].n = 1 ;
    r[i++].v[0] = phi ;
    strcpy(r[i].text,"c_oh") ; r[i].n = 1 ;
    r[i++].v[0] = c_oh ;
    strcpy(r[i].text,"potentiel_electrique") ; r[i].n = 1 ;
    r[i++].v[0] = psi ;
    strcpy(r[i].text,"charge") ; r[i].n = 1 ;
    r[i++].v[0] = c_q ;
    strcpy(r[i].text,"z_chcc") ; r[i].n = 1 ;
    r[i++].v[0] = z_chcc ;
    strcpy(r[i].text,"z_csh") ; r[i].n = 1 ;
    r[i++].v[0] = z_csh ;
    strcpy(r[i].text,"z_co2") ; r[i].n = 1 ;
    r[i++].v[0] = z_co2 ;
    strcpy(r[i].text,"P_CH/K_CH") ; r[i].n = 1 ;
    r[i++].v[0] = P_CH/K_CH ;
    strcpy(r[i].text,"P_Sam/K_Sam") ; r[i].n = 1 ;
    r[i++].v[0] = P_Sam/K_Sam ;
    strcpy(r[i].text,"c_h2sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2sio4 ;
    strcpy(r[i].text,"c_cah2sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_cah2sio4 ;
    strcpy(r[i].text,"c_caoh") ; r[i].n = 1 ;
    r[i++].v[0] = c_caoh ;
  }

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
  c_co2_eq = el.mat->pr[pm("C_CO2_eq")] ;

  /* initialisation */
  for(i=0;i<NVE_TR;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* molarites */
    double c_co2      = C_CO2(i) ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;

    /* solides */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_csh      = N_CSH(i) ;
    double n_sam      = N_Sam(i) ;

    /* porosite */
    double n_cc0      = N_CC0(i) ;
    double n_ch0      = N_CH0(i) ;
    double n_csh0     = N_CSH0(i) ;
    double n_sam0     = N_Sam0(i) ;
    double phi_th     = phi0 + v_ch*(n_ch0 - n_ch) + v_csh*(n_csh0 - n_csh) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
    double phi        = MAX(phi_th,phi_min) ;
    /* tortuosite liquide */
    double iff    = 2.9e-4*exp(9.95*phi) ;
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_CO2        += d_co2*iff ;
    KF_H2CO3      += d_h2co3*iff ;
    KF_HCO3       += d_hco3*iff ;
    KF_CO3        += d_co3*iff ;

    KF_Ca         += d_ca*iff ;
    KF_CaOH       += d_caoh*iff ;
    KF_CaHCO3     += d_cahco3*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaH3SiO4   += d_cah3sio4*iff ;

    KF_H2SiO4     += d_h2sio4*iff ;
    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    
    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HCO3     += FsRT*KF_HCO3*z_hco3*c_hco3 ;
    Kpsi_CO3      += FsRT*KF_CO3*z_co3*c_co3 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    Kpsi_CaHCO3   += FsRT*KF_CaHCO3*z_cahco3*c_cahco3 ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;

    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_caoh*Kpsi_CaOH + z_cahco3*Kpsi_CaHCO3 + z_h2sio4*Kpsi_H2SiO4 + z_h3sio4*Kpsi_H3SiO4 + z_cah3sio4*Kpsi_CaH3SiO4 ;
  }
  
  /* moyenne */
  for(i=0;i<NVE_TR;i++) va[i] *= 0.5 ;
}


void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les flux (f) */
{
  double r_h[2],r_oh[2] ;
  double r_co2[2],r_h2co3[2],r_hco3[2],r_co3[2] ;
  double r_ca[2],r_caoh[2],r_cahco3[2],r_cah2sio4[2],r_cah3sio4[2] ;
  double r_h2sio4[2],r_h3sio4[2],r_h4sio4[2] ;

  int    i ;

  for(i=0;i<2;i++) {
    double c_co2      = C_CO2(i) ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;

    double P_CC       = Prod_CC(z_co2,z_chcc) ;
    double P_CH       = Prod_CH(z_co2,z_chcc) ;
    double z_ch       = P_CH/P_CH_lim ;
    double P_Sam      = Prod_Sam(z_ch,z_csh) ;

    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_caoh     = K_caoh*c_ca*c_oh ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h4sio4   = P_Sam ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_h2sio4   = c_h3sio4/(K_h3sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4 ;

    r_co2[i]      = c_co2 ;
    r_oh[i]       = c_oh ;
    r_h[i]        = c_h ;
    r_h2co3[i]    = c_h2co3 ;
    r_hco3[i]     = c_hco3 ;
    r_co3[i]      = c_co3 ;
    r_ca[i]       = c_ca ;
    r_caoh[i]     = c_caoh ;
    r_cahco3[i]   = c_cahco3 ;
    r_h2sio4[i]   = c_h2sio4 ;
    r_h3sio4[i]   = c_h3sio4 ;
    r_h4sio4[i]   = c_h4sio4 ;
    r_cah2sio4[i] = c_cah2sio4 ;
    r_cah3sio4[i] = c_cah3sio4 ;
  }

  /* Gradients */
  {
    double dx = x[1][0] - x[0][0] ;
    double grd_h        = (r_h[1]        - r_h[0]       )/dx ;
    double grd_oh       = (r_oh[1]       - r_oh[0]      )/dx ;
    double grd_co2      = (r_co2[1]      - r_co2[0]     )/dx ;
    double grd_h2co3    = (r_h2co3[1]    - r_h2co3[0]   )/dx ;
    double grd_hco3     = (r_hco3[1]     - r_hco3[0]    )/dx ;
    double grd_co3      = (r_co3[1]      - r_co3[0]     )/dx ;
    double grd_caoh     = (r_caoh[1]     - r_caoh[0]    )/dx ;
    double grd_cahco3   = (r_cahco3[1]   - r_cahco3[0]  )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    double grd_cah2sio4 = (r_cah2sio4[1] - r_cah2sio4[0])/dx ;
    double grd_cah3sio4 = (r_cah3sio4[1] - r_cah3sio4[0])/dx ;
    double grd_h2sio4   = (r_h2sio4[1]   - r_h2sio4[0]  )/dx ;
    double grd_h3sio4   = (r_h3sio4[1]   - r_h3sio4[0]  )/dx ;
    double grd_h4sio4   = (r_h4sio4[1]   - r_h4sio4[0]  )/dx ;
    
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    
    /* Flux */
    double w_co2      = - KF_CO2*grd_co2       ; 
    double w_h2co3    = - KF_H2CO3*grd_h2co3   ;
    double w_hco3     = - KF_HCO3*grd_hco3          - Kpsi_HCO3*grd_psi     ;
    double w_co3      = - KF_CO3*grd_co3            - Kpsi_CO3*grd_psi      ;
    double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi     ;
    double w_cahco3   = - KF_CaHCO3*grd_cahco3      - Kpsi_CaHCO3*grd_psi   ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi       ;
    double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
    double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi   ;
    double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi   ;
    double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
    
    double w_q        = - z_h*KF_H*grd_h		      \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hco3*KF_HCO3*grd_hco3             \
                        - z_co3*KF_CO3*grd_co3		      \
                        - z_ca*KF_Ca*grd_ca		      \
                        - z_caoh*KF_CaOH*grd_caoh	      \
                        - z_cahco3*KF_CaHCO3*grd_cahco3	      \
                        - z_h2sio4*KF_H2SiO4*grd_h2sio4	      \
                        - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                        - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                        - Kpsi_q*grd_psi ;


    W_C     = w_co2 + w_h2co3 + w_hco3 + w_co3 + w_cahco3 ;
    W_Ca    = w_ca + w_caoh + w_cahco3 + w_cah2sio4 + w_cah3sio4 ;
    W_Si    = w_h2sio4 + w_h3sio4 + w_h4sio4 + w_cah2sio4 + w_cah3sio4 ;
    W_q     = w_q ;
  }
}


double concentration_oh(double c_co2,double z_csh,double z_chcc)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_co2, z_csh et z_chcc sont fixes */
  double c_co2_eq   = C_CO2_eq ;
  double z_co2      = c_co2/c_co2_eq ;
  /* les produits sont donc aussi fixes */
  double P_CC       = Prod_CC(z_co2,z_chcc) ;
  double P_CH       = Prod_CH(z_co2,z_chcc) ;
  double z_ch       = P_CH/P_CH_lim ;
  double P_Sam      = Prod_Sam(z_ch,z_csh) ;
  /* rappel des expressions c_i = A_i*(c_h)**n   : n = z_i
     c_h        = K_h2o/c_oh                     : +1
     c_h2co3    = K_h2co3*c_co2                  :  0
     c_hco3     = K_hco3*c_h2co3/c_h             : -1
     c_co3      = K_co3*c_hco3/c_h               : -2
     c_ca       = P_CC/c_co3                     : +2
     c_caoh     = K_caoh*c_ca*c_oh               : +1
     c_cahco3   = K_cahco3*c_ca*c_hco3           : +1
     c_h4sio4   = P_Sam                          :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_h2sio4   = c_h3sio4/(K_h3sio4*c_h)        : -2
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_cah2sio4 = K_cah2sio4*c_ca*c_h2sio4       : 0
  */

  double A_hco3     = K_hco3*K_h2co3*c_co2 ;
  double A_co3      = K_co3*A_hco3 ;
  double A_ca       = P_CC/A_co3 ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;
  double A_cahco3   = K_cahco3*A_ca*A_hco3 ;
  double A_h3sio4   = P_Sam/K_h4sio4 ;
  double A_h2sio4   = A_h3sio4/K_h3sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;

  double a = z_ca*A_ca ;
  double b = z_h + z_caoh*A_caoh + z_cahco3*A_cahco3 + z_cah3sio4*A_cah3sio4 ;
  /* double c = 0. ; *//* for the z_i = 0 */
  double d = z_oh*K_h2o + z_hco3*A_hco3 + z_h3sio4*A_h3sio4 ;
  double e = z_co3*A_co3 + z_h2sio4*A_h2sio4 ;

  double c_h = poly4_sansc(a,b,d,e) ;
 
  return(K_h2o/c_h) ;
}

double poly4(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double x0 = pow(-d/a,1./3) ;
  double x  = x0 ;
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
  } while(err > tol) ;
  return(x) ;
}

double poly4_sansc(double a,double b,double d,double e)
/* on resout ax^4 + bx^3 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double x0 = pow(-d/a,1./3) ;
  double x  = x0 ;
  int    i = 0 ;

  do {
    double x2 = x*x,x3 = x2*x ;
    double f  = x3 + (d*x + e)/(a*x + b) ;
    double df = 3*x2 + (d*b - e*a)/((a*x + b)*(a*x + b)) ;
    double dx = -f/df ;
    err = fabs(dx/x) ;
    x += dx ;
    if(i++ > 20) {
      printf("a = %e\n",a) ;
      printf("b = %e\n",b) ;
      printf("d = %e\n",d) ;
      printf("e = %e\n",e) ;
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      return(-fabs(x)) ;
    }
  } while(err > tol) ;
  return(x) ;
}
