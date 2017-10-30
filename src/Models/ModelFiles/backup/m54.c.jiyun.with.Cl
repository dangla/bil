/*

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../defs.h"
#include "../lib.h"

#define TITLE "Carbonatation du ciment sature en condition SC (mai 2010)"

/* Macros */
#define NEQ     (6)
#define NVE     (40)
#define NVI     (28)
#define NVE_TR  (30)

#define E_C     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)
#define E_Cl    (5)

#define I_CO2   (0)
#define I_psi   (1)
#define I_CHCC  (2)
#define I_CSH   (3)
#define I_K     (4)
#define I_Cl    (5)

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
#define X_K(n)     (u[(n)][I_K])
#define X_Cl(n)    (u[(n)][I_Cl])

#define Z_CHCCn(n) (u_n[(n)][I_CHCC])
#define Z_CSHn(n)  (u_n[(n)][I_CSH])
#define PSIn(n)    (u_n[(n)][I_psi])
#define X_Kn(n)    (u_n[(n)][I_K])
#define X_Cln(n)   (u_n[(n)][I_Cl])

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define N_K(n)     (f[(8+n)])
#define N_Cl(n)    (f[(10+n)])
#define W_C        (f[12])
#define W_q        (f[13])
#define W_Ca       (f[14])
#define W_Si       (f[15])
#define W_K        (f[16])
#define W_Cl       (f[17])
#define N_CH(n)    (f[(18+n)])
#define N_CC(n)    (f[(20+n)])
#define N_Jen(n)   (f[(22+n)])
#define N_Sam(n)   (f[(24+n)])
#define N_TobI(n)  (f[(26+n)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define N_CHn(n)   (f_n[(18+n)])
#define N_CCn(n)   (f_n[(20+n)])
#define N_Jenn(n)  (f_n[(22+n)])
#define N_Samn(n)  (f_n[(24+n)])
#define N_TobIn(n) (f_n[(26+n)])


#define KF_OH       (va[(0)])
#define KF_H        (va[(1)])
#define KF_CO2      (va[(2)])
#define KF_H2CO3    (va[(3)])
#define KF_HCO3     (va[(4)])
#define KF_CO3      (va[(5)])
#define KF_Ca       (va[(6)])
#define KF_CaHCO3   (va[(7)])
#define KF_CaH3SiO4 (va[(8)])
#define KF_H3SiO4   (va[(9)])
#define KF_H4SiO4   (va[(10)])
#define KF_H2SiO4   (va[(11)])
#define KF_CaH2SiO4 (va[(12)])
#define KF_CaCO3aq  (va[(13)])
#define KF_CaOH     (va[(14)])
#define KF_K        (va[(15)])
#define KF_Cl       (va[(16)])

#define Kpsi_OH       (va[(17)])
#define Kpsi_H        (va[(18)])
#define Kpsi_HCO3     (va[(19)])
#define Kpsi_CO3      (va[(20)])
#define Kpsi_Ca       (va[(21)])
#define Kpsi_CaHCO3   (va[(22)])
#define Kpsi_CaH3SiO4 (va[(23)])
#define Kpsi_H3SiO4   (va[(24)])
#define Kpsi_q        (va[(25)])
#define Kpsi_H2SiO4   (va[(26)])
#define Kpsi_CaOH     (va[(27)])
#define Kpsi_K        (va[(28)])
#define Kpsi_Cl       (va[(29)])

#define N_CH0(n)      (va[(30+n)])
#define N_CC0(n)      (va[(32+n)])
#define N_Jen0(n)     (va[(34+n)])
#define N_Sam0(n)     (va[(36+n)])
#define N_TobI0(n)    (va[(38+n)])

/*
  Solution aqueuse
*/

/* les valences */
#define z_ca       (2.)
#define z_h        (1.)
#define z_oh       (-1.)
#define z_hco3     (-1.)
#define z_co3      (-2.)
#define z_h3sio4   (-1.)
#define z_cahco3   (1.)
#define z_cah3sio4 (1.)
#define z_h2sio4   (-2.)
#define z_caoh     (1.)
#define z_k        (1.)
#define z_cl       (-1.)

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
#define v_caco3aq  (26.20e-3)    /* a modifier */
#define v_caoh     (26.20e-3)    /* a modifier */
#define v_k        (43.93e-3)    /* d'apres Antoine */
#define v_koh      (27.44e-3)    /* d'apres Antoine */
#define v_cl       (xxx)         /* d'apres Antoine */

/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s  */
#define d_oh       (5.273e-7)     /* chemistry handbook */
#define d_h        (9.311e-7)    /* chemistry handbook */
#define d_co2      (1.91e-7)    /* chemistry handbook */
#define d_h2co3    (7.2e-8)
#define d_hco3     (1.18e-7)    /* 1.07e-7 (radius = 2e-10 m) */
#define d_co3      (9.55e-8)    /* 9.53e-8 (radius = 2.25e-10 m) */
#define d_ca       (7.92e-8)
#define d_cahco3   (1.07e-7)    /* (radius = 2e-10 m) */
#define d_sioh4    (xxx)
#define d_h4sio4   (1.07e-7)    /*  */
#define d_h3sio4   (1.07e-7)    /*(radius = 2e-10 m) */
#define d_h2sio4   (1.07e-7)   /* (radius = 2e-10 m) */
#define d_cah2sio4 (1.07e-7)    /* a modifier */
#define d_cah3sio4 (1.07e-7)    /* (radius = 2e-10 m) */
#define d_caco3aq  (1.43e-7)    /* (radius = 1.5e-10 m) */
#define d_caoh     (1.07e-7)    /* (radius = 2e-10 m) */
#define d_k        (1.957e-7)   /* chemistry handbook */
#define d_koh      (1.96e-7)   
#define d_cl       (2.032e-7)   /* chemistry handbook */
#define d_na       (1.334e-7)   /* chemistry handbook */


/* constantes d'equilibre (ref = 1 mole/L) */
#define K_h2o      (1.e-14)          /* autoprotolyse de l'eau */
#define K_henry    (1.238)           /* cste de Henry du CO2 / RT */

#define K_h2co3    (1.7e-3)          /* CO2[0] + H2O = H2CO3 */
#define K_hco3     (2.5e-4)          /* H2CO3   = HCO3[-] + H[+] */
#define K_co3      (4.69e-11)        /* HCO3[-] = CO3[2-] + H[+] */

#define K_h2sio4   (4.68)            /* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define K_h4sio4   (6.45e9)          /* H3SiO4[-] + H[+] = H4SiO4 */
#define K_h3sio4   (1.55e-10)        /* H4SiO4    = H3SiO4[-] + H[+] */

#define K_cahco3   (1.276e+1)        /* Ca[2+] + HCO3[-]    = CaHCO3[+] */
#define K_caco3aq  (1.4e+3)          /* Ca[2+] + CO3[2-]    = CaCO3[0] */
#define K_cah2sio4 (3.98e+4)         /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define K_cah3sio4 (1.58e+1)         /* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define K_caoh     (1.66e+1)         /* Ca[2+] + OH[-] = CaOH[+] */

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
#define v_sam      (29.e-3)      /* d'apres Lothenbach */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* CH = Ca[2+] + 2OH[-] */
#define K_CC       (3.89e-9)         /* CC = Ca[2+] + CO3[2-]  */ 
#define K_Sam      (1.93642e-3)      /* S + 2H2O = H4SiO4 */
/* CxSyHz = xCa[2+] + 2xOH[-] + yH4SiO4 + (z-x-2y)H2O */
#define K_CSH      (K_Jen)
#define X_CSH      (X_Jen)
#define Y_CSH      (Y_Jen)
#define v_csh      (V_Jen)

/* C-S-H */
/* Tobermorite I  (x,y,z) = (2,2.4,4.4) */
#define K_TobI     (1.685e-12)        /* -28.256/2.4 */
#define X_TobI     (5./6.)
#define Y_TobI     (1.)
#define V_TobI     (54.78e-3)         /* X_CSH*59.e-3 d'apres Lothenbach */
/* Tobermorite II (x,y,z) = (1.5,1.8,3.3) */
#define K_TobII    (1.685e-12)        /* -21.192/1.8*/
#define X_TobII    (5./6.)
#define Y_TobII    (1.)
#define V_TobII    (54.78e-3)  
/* Jennite (x,y,z) = (5/3,1,8/3) */
#define K_Jen      (4.40e-18)        /* Donnée GEMS e(-15.621/0.9) */
#define X_Jen      (5./3.)
#define Y_Jen      (1.)
#define V_Jen      (81.78e-3)

#define C_CO2_eq   (K_h2o*K_h2o*K_CC/(K_h2co3*K_hco3*K_co3*K_CH))


/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,double,double,double,double) ;
static double poly4(double,double,double,double,double) ;
static double poly4_sansc(double,double,double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
static double product_Sam(double,double,double) ;
static double product_Jen(double,double,double) ;
static double product_TobI(double,double,double) ;
#define MIN(a,b)  ((a < b) ? a : b)
#define MAX(a,b)  ((a > b) ? a : b)
#define NEGEXP(y) ((y < 0.) ? exp(y) : 1.)

/* Parametres */
static double phi0,c_co2_eq,t_ch,t_cc ;
static double n_ch_ref,n_cc_ref ;


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"N_CSH") == 0)    return (2) ;
  else if(strcmp(s,"N_CC") == 0)     return (3) ;
  else if(strcmp(s,"N_SAM") == 0)    return (4) ;
  else if(strcmp(s,"C_CO2_eq") == 0) return (5) ;
  else if(strcmp(s,"T_CH") == 0)     return (6) ;
  else if(strcmp(s,"T_CC") == 0)     return (7) ;
  else if(strcmp(s,"courbes") == 0)  return (8) ;
  else return(-1) ;
}


int dm55(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 11 ;
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_C],   "carbone") ;
  strcpy(mat->eqn[E_Ca],  "calcium") ;
  strcpy(mat->eqn[E_Si],  "Silicium") ;
  strcpy(mat->eqn[E_q],   "charge") ;
  strcpy(mat->eqn[E_K],   "potassium") ;
  strcpy(mat->eqn[E_Cl],   "chloride") ;  

#if (U_CO2 == LOG_RHO)
  strcpy(mat->inc[I_CO2],  "logc_co2") ;
#else
  strcpy(mat->inc[I_CO2],  "c_co2") ;
#endif
  strcpy(mat->inc[I_CHCC], "z_chcc") ;
  strcpy(mat->inc[I_psi],  "psi") ;
  strcpy(mat->inc[I_CSH],  "z_csh") ;
  strcpy(mat->inc[I_K],    "c_k") ;
  strcpy(mat->inc[I_Cl],   "c_cl") ;

  {
    /* initialisation automatique */
    double t_ch        = 600. ;
    double t_cc        = 0. ;
    double n_ch_ref    = 1. ;
    double n_cc_ref    = 0. ;

    mat->pr[pm("N_CH")]  = n_ch_ref ;
    mat->pr[pm("N_CC")]  = n_cc_ref ;
    mat->pr[pm("T_CH")]  = t_ch ;
    mat->pr[pm("T_CC")]  = t_cc ;

    dmat(mat,ficd,pm,n_donnees) ;

    n_ch_ref  = mat->pr[pm("N_CH")] ;
    n_cc_ref  = mat->pr[pm("N_CC")] ;

    t_ch      = mat->pr[pm("T_CH")] ;
    t_cc      = mat->pr[pm("T_CC")] ;

    if(n_cc_ref  == 0.) mat->pr[pm("N_CC")]  = n_ch_ref ;
    if(t_cc  == 0.) mat->pr[pm("T_CC")]  = t_ch ;

    mat->pr[pm("C_CO2_eq")] = C_CO2_eq ;
  }
  
  return(mat->n) ;
}


int qm55(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  printf("Le systeme est forme de 6 equations:\n") ;
#if (U_CO2 == LOG_RHO)
  printf("\t- la conservation de la masse de C      (logc_co2)\n") ;
#else
  printf("\t- la conservation de la masse de C      (c_co2)\n") ;
#endif
  printf("\t- la conservation de la charge          (psi)\n") ;
  printf("\t- la conservation de la masse de Ca     (z_chcc)\n") ;
  printf("\t- la conservation de la masse de Si     (z_csh)\n") ;
  printf("\t- la conservation de la masse de K      (c_k)\n") ;
  printf("\t- la conservation de la masse de Cl     (c_cl)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n") ;

  printf("Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1    # Contenu en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"z_csh = 2.4    # contenu en CSH (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4    # contenu en K (moles/L)\n") ;  
  fprintf(ficd,"N_Cl  = 0.4    # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5   # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"courbes = solid   # Nom du fichier: log10(P_CH/K_CH) Xi\n") ;

  return(NEQ) ;
}



void tb55(elem_t *el,inte_t *fi,int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = (el->dim < dim) ? 0 : NVE ;
}


void ch55(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in55(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;
  n_cc_ref  = el.mat->pr[pm("N_CC")] ;
  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double c_co2      = (C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;
    double z_co2      = c_co2/c_co2_eq ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double c_k        = X_K(i);                 /* ajouter2 */
    double c_cl       = X_Cl(i); 
    
    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;
    double X_CH 	  = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
    double A_Jen 	  = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
    double A_TobI 	  = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
    double sumx 	  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1. ;
    double x_s        = NEGEXP(z_csh)/sumx;
    double x_jen      = A_Jen*pow(X_CH,X_Jen)*x_s ;
    double x_tobI	  = NEGEXP(z_csh) - x_s - x_jen ;

    double c_h4sio4   = x_s * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ; 
    
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */
    
    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(z_chcc,0.) ;
    double n_cc_eq    = n_cc_ref*MAX(z_chcc,0.) ;
    
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch       = (z_co2 <= 1) ? n_ch_eq  : 0 ;
    double n_cc       = (z_co2 >  1) ? n_cc_eq  : 0 ;
    double n_jen      = x_jen * MAX(z_csh,0) ;
    double n_sam      = x_s * MAX(z_csh,0) ;
    double n_tobI	  = x_tobI* MAX(z_csh,0) ;
      
    /* solides */
    double n_si_s     = Y_TobI*n_tobI + n_sam + Y_Jen*n_jen;
    double n_ca_s     = n_ch + n_cc + X_Jen*n_jen + X_TobI*n_tobI ;

    /* porosite */
    double phi = phi0 ;

    /* contenus molaires */
    double n_co2      = phi*c_co2 ;
    double n_hco3     = phi*c_hco3 ;
    double n_h2co3    = phi*c_h2co3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_h2sio4   = phi*c_h2sio4 ;          /* ajouter */
    double n_cah2sio4 = phi*c_cah2sio4 ;        /* ajouter */
    double n_caco3aq  = phi*c_caco3aq ;         /* ajouter */
    double n_caoh     = phi*c_caoh ;            /* ajouter */
    double n_k        = phi*c_k;                /* ajouter2 */
    double n_cl       = phi*c_cl;     
    
    
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc + n_caco3aq ;
    N_Ca(i) = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh + n_ca_s ;
    N_Si(i) = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4+ n_cah2sio4  + n_si_s ;
    N_K(i)  = n_k;                              /* ajouter2 */
    N_Cl(i) = n_cl;

    /* densite de charge */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 \
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k + z_cl*c_cl;

#if (U_CO2 == LOG_RHO)
    u[i][I_CO2] = log(c_co2)/Ln10 ;
#else
    C_CO2(i)   = c_co2 ;
#endif
    Z_CSH(i)   = z_csh ;
    Z_CHCC(i)  = z_chcc ;
    X_K(i)     = c_k;
    X_Cl(i)    = c_cl;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_Jen(i)   = n_jen ;
    N_Sam(i)   = n_sam ;
    N_TobI(i)  = n_tobI ;

    N_CH0(i)   = n_ch ;
    N_CC0(i)   = n_cc ;
    N_Jen0(i)  = n_jen ;
    N_Sam0(i)  = n_sam ;
    N_TobI0(i) = n_tobI ;
  }
  
  if(el.dim < dim) return ;

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex55(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  
  if(el.dim < dim) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct55(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  int    i ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;
  n_cc_ref  = el.mat->pr[pm("N_CC")] ;
  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;
  t_cc      = el.mat->pr[pm("T_CC")] ;
  
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double c_co2      = C_CO2(i) ;
    double z_csh      = Z_CSH(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;
    double c_k        = X_K(i);                 /* ajouter2 */
    double c_cl       = X_Cl(i); 

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;
    double X_CH 	  = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
    double A_Jen 	  = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
    double A_TobI 	  = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
    double sumx 	  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1. ;
    double x_s        = NEGEXP(z_csh)/sumx;
    double x_jen      = A_Jen*pow(X_CH,X_Jen)*x_s ;
    double x_tobI	  = NEGEXP(z_csh) - x_s - x_jen ;

    double c_h4sio4   = x_s * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl)  ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */

    /* cinetiques */
    double n_chn      = N_CHn(i) ;
    double n_ccn      = N_CCn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;  /* si z_co2 > 1 */
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;   /* si z_co2 < 1 */

    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(z_chcc,0.) ;
    double n_cc_eq    = n_cc_ref*MAX(z_chcc,0.) ;

    
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch   = (z_co2 <= 1) ? n_ch_eq  : n_ch_ci ;
    double n_cc   = (z_co2 >  1) ? n_cc_eq  : n_cc_ci ;
    double n_jen      = x_jen * MAX(z_csh,0) ;
    double n_sam      = x_s * MAX(z_csh,0) ;
    double n_tobI	  = x_tobI* MAX(z_csh,0) ;
      
    /* solides */
    double n_si_s     = Y_TobI*n_tobI + n_sam + Y_Jen*n_jen;
    double n_ca_s     = n_ch + n_cc + X_Jen*n_jen + X_TobI*n_tobI ;

    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_jen0     = N_Jen0(i) ;
    double n_tobI0	  = N_TobI0(i) ;    
    double n_sam0     = N_Sam0(i) ;
    double phi        = phi0 + v_ch*(n_ch0 - n_ch) + V_Jen*(n_jen0- n_jen) + V_TobI*(n_tobI0- n_tobI) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;

    /* contenus molaires */
    double n_co2      = phi*c_co2 ;
    double n_h2co3    = phi*c_h2co3 ;
    double n_hco3     = phi*c_hco3 ;
    double n_co3      = phi*c_co3 ;
    double n_ca       = phi*c_ca ;
    double n_cahco3   = phi*c_cahco3 ;
    double n_h3sio4   = phi*c_h3sio4 ;
    double n_h4sio4   = phi*c_h4sio4 ;
    double n_cah3sio4 = phi*c_cah3sio4 ;
    double n_h2sio4   = phi*c_h2sio4 ;          /* ajouter */
    double n_cah2sio4 = phi*c_cah2sio4 ;        /* ajouter */
    double n_caco3aq  = phi*c_caco3aq ;         /* ajouter */
    double n_caoh     = phi*c_caoh ;            /* ajouter */
    double n_k        = phi*c_k;                /* ajouter2 */
    double n_cl       = phi*c_cl;

    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc + n_caco3aq ;
    N_Ca(i) = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh + n_ca_s ;
    N_Si(i) = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4+ n_cah2sio4  + n_si_s ;
    N_K(i)  = n_k;
    N_Cl(i) = n_cl;

    /* densite de charge */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 + \
    z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k + z_cl*c_cl ;

    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_Jen(i)   = n_jen ;
    N_TobI(i)  = n_tobI ;
    N_Sam(i)   = n_sam ;

    if(c_co2 < 0. || n_ca_s < 0. || n_si_s < 0.) {
      printf("x         = %e\n",x[i][0]) ;
      printf("c_co2     = %e\n",c_co2) ;
      printf("n_cc      = %e\n",n_cc) ;
      printf("n_jen     = %e\n",n_jen) ;
      printf("n_ca_s    = %e\n",n_ca_s) ; 
      printf("n_si_s    = %e\n",n_si_s) ;
      printf("z_csh     = %e\n",z_csh) ;
      printf("z_chcc    = %e\n",z_chcc) ;
      printf("c_h3sio4  = %e\n",c_h3sio4) ;
      printf("c_oh      = %e\n",c_oh) ;
      return(-1) ;
    }
    if(phi < 0.) {
      printf("phi = %e\n",phi) ;
      return(-1) ;
    }
  }
  
  if(el.dim < dim) return(0) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;

  return(0) ;
}

int mx55(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
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
  double Dc_cahco3SDc_co2[2]   ;
  double Dc_h3sio4SDc_co2[2]   ;
  double Dc_h4sio4SDc_co2[2]   ;
  double Dc_cah3sio4SDc_co2[2] ;
  double Dc_h2sio4SDc_co2[2]   ;
  double Dc_cah2sio4SDc_co2[2] ;
  double Dc_caco3aqSDc_co2[2]  ;  
  double Dc_caohSDc_co2[2]     ;
  
  double Dc_hSDz_csh[2]        ;
  double Dc_ohSDz_csh[2]       ;
  double Dc_hco3SDz_csh[2]     ;
  double Dc_co3SDz_csh[2]      ;
  double Dc_caSDz_csh[2]       ;
  double Dc_cahco3SDz_csh[2]   ;
  double Dc_h3sio4SDz_csh[2]   ;
  double Dc_h4sio4SDz_csh[2]   ;
  double Dc_cah3sio4SDz_csh[2] ;
  double Dc_h2sio4SDz_csh[2]   ;
  double Dc_cah2sio4SDz_csh[2] ;
  double Dc_caco3aqSDz_csh[2]  ;
  double Dc_caohSDz_csh[2]     ;

  double Dc_hSDz_chcc[2]        ;
  double Dc_ohSDz_chcc[2]       ;
  double Dc_hco3SDz_chcc[2]     ;
  double Dc_co3SDz_chcc[2]      ;
  double Dc_caSDz_chcc[2]       ;
  double Dc_cahco3SDz_chcc[2]   ;
  double Dc_h3sio4SDz_chcc[2]   ;
  double Dc_h4sio4SDz_chcc[2]   ;
  double Dc_cah3sio4SDz_chcc[2] ;
  double Dc_h2sio4SDz_chcc[2]   ;
  double Dc_cah2sio4SDz_chcc[2] ;
  double Dc_caco3aqSDz_chcc[2]  ;  
  double Dc_caohSDz_chcc[2]     ;
  
  double Dc_hSDc_k[2]         ;
  double Dc_ohSDc_k[2]        ;
  double Dc_hco3SDc_k[2]      ;
  double Dc_co3SDc_k[2]       ;
  double Dc_caSDc_k[2]        ;
  double Dc_cahco3SDc_k[2]    ;
  double Dc_h3sio4SDc_k[2]    ;
  double Dc_cah3sio4SDc_k[2]  ;
  double Dc_h2sio4SDc_k[2]    ;
  double Dc_cah2sio4SDc_k[2]  ;
  double Dc_caco3aqSDc_k[2]   ;  
  double Dc_caohSDc_k[2]      ;
  
  double Dc_hSDc_cl[2]         ;
  double Dc_ohSDc_cl[2]        ;
  double Dc_hco3SDc_cl[2]      ;
  double Dc_co3SDc_cl[2]       ;
  double Dc_caSDc_cl[2]        ;
  double Dc_cahco3SDc_cl[2]    ;
  double Dc_h3sio4SDc_cl[2]    ;
  double Dc_cah3sio4SDc_cl[2]  ;
  double Dc_h2sio4SDc_cl[2]    ;
  double Dc_cah2sio4SDc_cl[2]  ;
  double Dc_caco3aqSDc_cl[2]   ;  
  double Dc_caohSDc_cl[2]      ;
  
  /*
    Donnees
  */
  phi0      = el.mat->pr[pm("porosite")] ;
  n_ch_ref  = el.mat->pr[pm("N_CH")] ;

  n_cc_ref  = el.mat->pr[pm("N_CC")] ;

  c_co2_eq  = el.mat->pr[pm("C_CO2_eq")] ;
  t_ch      = el.mat->pr[pm("T_CH")] ;

  t_cc      = el.mat->pr[pm("T_CC")] ;


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
    double z_csh      = Z_CSH(i);
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;
    double c_k        = X_K(i);                 /* ajouter2 */
    double c_cl       = X_Cl(i);  

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;
    double X_CH 	  = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
    double A_Jen 	  = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
    double A_TobI 	  = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
    double sumx 	  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1. ;
    double x_s        = NEGEXP(z_csh)/sumx;
    double x_jen      = A_Jen*pow(X_CH,X_Jen)*x_s ;
    double x_tobI	  = NEGEXP(z_csh) - x_s - x_jen ;

    double c_h4sio4   = x_s  * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */    

    /* contenus en CH, CC, CSH, Sam */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_tobI     = N_TobI(i) ;
    double n_jen	  = N_Jen(i) ;
    double n_sam      = N_Sam(i) ;

    /* cinetiques */
    double n_ccn      = N_CCn(i) ;
    double n_chn      = N_CHn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;


    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double n_tobI0    = N_TobI0(i) ;
    double n_jen0	  = N_Jen0(i) ;
    double n_sam0     = N_Sam0(i) ;
    double phi        = phi0 + v_ch*(n_ch0 - n_ch) + V_Jen*(n_jen0- n_jen) + V_TobI*(n_tobI0- n_tobI) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;


    /* derivees ... */
    /* ... par rapport a c_co2 */
    double dc_co2             = 1.e-4*c_co2*((c_co2 > C_CO2n(i)) ? 1 : -1) ;
    double c_co22             = c_co2 + dc_co2 ;
    double c_oh2              = concentration_oh(c_co22,z_csh,z_chcc,c_k,c_cl)  ;   
    double dc_ohsdc_co2       = (c_oh2 - c_oh)/dc_co2 ;
    double dc_hsdc_co2        = - c_h*dc_ohsdc_co2/c_oh ;
    double dc_h2co3sdc_co2    = K_h2co3 ;
    double dc_hco3sdc_co2     = c_hco3*(dc_h2co3sdc_co2/c_h2co3	\
			                    + dc_ohsdc_co2/c_oh) ;
    double dc_co3sdc_co2      = c_co3*(dc_hco3sdc_co2/c_hco3 \
			                    + dc_ohsdc_co2/c_oh) ;		                    
			                     
    double dP_CCsdc_co2       = (z_co2 < 1.) ? P_CC/c_co2 : 0. ;
    
    double dc_casdc_co2       = (dP_CCsdc_co2 - c_ca*dc_co3sdc_co2)/c_co3 ;
    double dc_cahco3sdc_co2   = K_cahco3*(dc_casdc_co2*c_hco3 \
			                    + c_ca*dc_hco3sdc_co2) ;
			                    
    double x_s2				  =  product_Sam(c_co22,z_chcc,z_csh) ;
    double x_jen2		      =  product_Jen(c_co22,z_chcc,z_csh) ;
    double x_tobI2		      =  NEGEXP(z_csh) - x_s2 - x_jen2 ;
    double dc_h4sio4sdc_co2   = (x_s2 - x_s )* K_Sam/dc_co2 ;
    double dc_h3sio4sdc_co2   = c_h3sio4*(dc_h4sio4sdc_co2/c_h4sio4 \
			                    - dc_hsdc_co2/c_h) ;
    double dc_cah3sio4sdc_co2 = K_cah3sio4*(dc_casdc_co2*c_h3sio4 \
			                    + c_ca*dc_h3sio4sdc_co2) ;
    double dc_h2sio4sdc_co2   = K_h2sio4*(dc_h3sio4sdc_co2*c_oh	\
			                    + dc_ohsdc_co2*c_h3sio4) ;
    double dc_cah2sio4sdc_co2 = K_cah2sio4*(dc_h2sio4sdc_co2*c_ca	\
			                    + dc_casdc_co2*c_h2sio4) ;
    double dc_caco3aqsdc_co2  = K_caco3aq*(dc_co3sdc_co2*c_ca \
			                    + dc_casdc_co2*c_co3) ;				                  
    double dc_caohsdc_co2     = K_caoh*(dc_casdc_co2*c_oh \
			                    + dc_ohsdc_co2*c_ca) ;			                  	
		                    
    double dn_ch_cisdc_co2   = n_ch_ci*(-dt/t_ch)/c_co2 ;
    double dn_cc_cisdc_co2   = n_cc_ci*(dt/t_cc)/c_co2 ;


    double dn_ch_eqsdc_co2   = 0 ;
    double dn_cc_eqsdc_co2   = 0 ;


    double dn_chsdc_co2   = (z_co2 <= 1) ? dn_ch_eqsdc_co2  : dn_ch_cisdc_co2 ;
    double dn_ccsdc_co2   = (z_co2 >  1) ? dn_cc_eqsdc_co2  : dn_cc_cisdc_co2 ;
    double dn_samsdc_co2  = (x_s2 - x_s) * MAX(z_csh,0) / dc_co2 ;
    double dn_jensdc_co2  = (x_jen2 - x_jen) * MAX(z_csh,0) / dc_co2 ;
    double dn_tobIsdc_co2 = (x_tobI2 - x_tobI) * MAX(z_csh,0) / dc_co2 ;

    double dn_ca_ssdc_co2 = dn_chsdc_co2 + dn_ccsdc_co2 + X_Jen*dn_jensdc_co2 + X_TobI*dn_tobIsdc_co2 ;
    double dn_si_ssdc_co2 = Y_Jen*dn_jensdc_co2 + Y_TobI*dn_tobIsdc_co2 + dn_samsdc_co2 ;

    double dphisdc_co2        = - v_ch*dn_chsdc_co2 - V_Jen*dn_jensdc_co2- V_TobI*dn_tobIsdc_co2 - v_sam*dn_samsdc_co2 - v_cc*dn_ccsdc_co2 ;

    double dn_co2sdc_co2      = phi + dphisdc_co2*c_co2 ;
    double dn_h2co3sdc_co2    = phi*dc_h2co3sdc_co2 + dphisdc_co2*c_h2co3 ;
    double dn_hco3sdc_co2     = phi*dc_hco3sdc_co2 + dphisdc_co2*c_hco3 ;
    double dn_co3sdc_co2      = phi*dc_co3sdc_co2 + dphisdc_co2*c_co3 ;
    double dn_casdc_co2       = phi*dc_casdc_co2 + dphisdc_co2*c_ca ;
    double dn_cahco3sdc_co2   = phi*dc_cahco3sdc_co2 + dphisdc_co2*c_cahco3 ;
    double dn_cah3sio4sdc_co2 = phi*dc_cah3sio4sdc_co2 + dphisdc_co2*c_cah3sio4 ;
    double dn_h3sio4sdc_co2   = phi*dc_h3sio4sdc_co2 + dphisdc_co2*c_h3sio4 ;
    double dn_h4sio4sdc_co2   = phi*dc_h4sio4sdc_co2 + dphisdc_co2*c_h4sio4 ;
    double dn_h2sio4sdc_co2   = phi*dc_h2sio4sdc_co2 + dphisdc_co2*c_h2sio4 ;
    double dn_cah2sio4sdc_co2 = phi*dc_cah2sio4sdc_co2 + dphisdc_co2*c_cah2sio4 ;
    double dn_caco3aqsdc_co2  = phi*dc_caco3aqsdc_co2 + dphisdc_co2*c_caco3aq ;
    double dn_caohsdc_co2     = phi*dc_caohsdc_co2 + dphisdc_co2*c_caoh ;
    double dn_ksdc_co2        = dphisdc_co2*c_k;
    double dn_clsdc_co2       = dphisdc_co2*c_cl;
			                    
    /* par rapport a z_csh */
    double dz_csh             = ((z_csh > Z_CSHn(i)) ? 1 : -1)*1.e-4 ; 
    double z_csh2             = z_csh + dz_csh ;
    
    double c_oh1              = concentration_oh(c_co2,z_csh2,z_chcc,c_k,c_cl)  ;
    double dc_ohsdz_csh       = (c_oh1 - c_oh)/dz_csh ;
    double dc_hsdz_csh        = - c_h*dc_ohsdz_csh/c_oh ;
    double dc_hco3sdz_csh     = c_hco3*(dc_ohsdz_csh/c_oh) ;
    double dc_co3sdz_csh      = c_co3*(dc_hco3sdz_csh/c_hco3 \
								+ dc_ohsdz_csh/c_oh) ;
    double dc_casdz_csh       = - c_ca*dc_co3sdz_csh/c_co3 ;
    double dc_cahco3sdz_csh   = K_cahco3*(dc_casdz_csh*c_hco3 \
								+ c_ca*dc_hco3sdz_csh) ;

    double x_s3         	  = product_Sam(c_co2,z_chcc,z_csh2) ;
    double x_jen3		      = product_Jen(c_co2,z_chcc,z_csh2) ;
    double x_tobI3		      = NEGEXP(z_csh2) - x_s3 - x_jen3 ;

    double dc_h4sio4sdz_csh   = (x_s3 - x_s)*K_Sam/dz_csh ;
    double dc_h3sio4sdz_csh   = c_h3sio4*(dc_h4sio4sdz_csh/c_h4sio4 \
								- dc_hsdz_csh/c_h) ;
    double dc_cah3sio4sdz_csh = K_cah3sio4*(dc_casdz_csh*c_h3sio4 \
								+ c_ca*dc_h3sio4sdz_csh) ;
    double dc_h2sio4sdz_csh   = K_h2sio4*(dc_h3sio4sdz_csh*c_oh	\
			                    + dc_ohsdz_csh*c_h3sio4) ;
    double dc_cah2sio4sdz_csh = K_cah2sio4*(dc_h2sio4sdz_csh*c_ca	\
			                    + dc_casdz_csh*c_h2sio4) ;
    double dc_caco3aqsdz_csh  = K_caco3aq*(dc_co3sdz_csh*c_ca \
			                    + dc_casdz_csh*c_co3) ;				                  
    double dc_caohsdz_csh     = K_caoh*(dc_casdz_csh*c_oh \
			                    + dc_ohsdz_csh*c_ca) ;				      

    double dn_samsdz_csh 	  = (x_s3 * MAX(z_csh2,0) - x_s * MAX(z_csh,0))/dz_csh ;
    double dn_jensdz_csh 	  = (x_jen3 * MAX(z_csh2,0)- x_jen* MAX(z_csh,0)) / dz_csh ;
    double dn_tobIsdz_csh	  = (x_tobI3 * MAX(z_csh2,0)- x_tobI* MAX(z_csh,0)) / dz_csh ;
    

    double dn_ca_ssdz_csh 	  = X_Jen*dn_jensdz_csh + X_TobI*dn_tobIsdz_csh;
    double dn_si_ssdz_csh 	  = Y_Jen*dn_jensdz_csh + Y_TobI*dn_tobIsdz_csh +dn_samsdz_csh ;

    double dphisdz_csh   	  = -V_Jen*dn_jensdz_csh -V_TobI*dn_tobIsdz_csh - v_sam*dn_samsdz_csh ;

    double dn_co2sdz_csh      = dphisdz_csh*c_co2 ;
    double dn_h2co3sdz_csh    = dphisdz_csh*c_h2co3 ;
    double dn_hco3sdz_csh     = phi*dc_hco3sdz_csh + dphisdz_csh*c_hco3 ;
    double dn_co3sdz_csh      = phi*dc_co3sdz_csh + dphisdz_csh*c_co3 ;
    double dn_casdz_csh       = phi*dc_casdz_csh + dphisdz_csh*c_ca ;
    double dn_cahco3sdz_csh   = phi*dc_cahco3sdz_csh + dphisdz_csh*c_cahco3 ;
    double dn_cah3sio4sdz_csh = phi*dc_cah3sio4sdz_csh + dphisdz_csh*c_cah3sio4 ;
    double dn_h3sio4sdz_csh   = phi*dc_h3sio4sdz_csh + dphisdz_csh*c_h3sio4 ;
    double dn_h4sio4sdz_csh   = phi*dc_h4sio4sdz_csh + dphisdz_csh*c_h4sio4 ;
    double dn_h2sio4sdz_csh   = phi*dc_h2sio4sdz_csh + dphisdz_csh*c_h2sio4 ;
    double dn_cah2sio4sdz_csh = phi*dc_cah2sio4sdz_csh + dphisdz_csh*c_cah2sio4 ;
    double dn_caco3aqsdz_csh  = phi*dc_caco3aqsdz_csh + dphisdz_csh*c_caco3aq ;
    double dn_caohsdz_csh     = phi*dc_caohsdz_csh + dphisdz_csh*c_caoh ;
    double dn_ksdz_csh        = dphisdz_csh*c_k;
    double dn_clsdz_csh       = dphisdz_csh*c_cl;
    
    /* par rapport a z_chcc */
    /* double dz_chcc             = ((z_chcc > 0.) ? 1 : -1)*1.e-2 ; */
    double dz_chcc             = 1.e-6*((z_chcc > Z_CHCCn(i)) ? 1 : -1) ;
    double z_chcc2             = z_chcc + dz_chcc ;
    double c_oh3               = concentration_oh(c_co2,z_csh,z_chcc2,c_k,c_cl) ;
    double dc_ohsdz_chcc       = (c_oh3 - c_oh)/dz_chcc ;
    double dc_hsdz_chcc        = - c_h*dc_ohsdz_chcc/c_oh ;
    double dc_hco3sdz_chcc     = c_hco3*(dc_ohsdz_chcc/c_oh) ;
    double dc_co3sdz_chcc      = c_co3*(dc_hco3sdz_chcc/c_hco3 \
								 + dc_ohsdz_chcc/c_oh) ;

    double dP_CCsdz_chcc       = (z_chcc < 0) ? P_CC : 0 ;
    double dc_casdz_chcc       = (dP_CCsdz_chcc - c_ca*dc_co3sdz_chcc)/c_co3 ;
    double dc_cahco3sdz_chcc   = K_cahco3*(dc_casdz_chcc*c_hco3 \
								 + c_ca*dc_hco3sdz_chcc) ;
    
    double x_s4         	   = product_Sam(c_co2,z_chcc2,z_csh) ;
    double x_jen4		       = product_Jen(c_co2,z_chcc2,z_csh) ;
    double x_tobI4		       = NEGEXP(z_csh) - x_s4 - x_jen4 ;

    double dc_h4sio4sdz_chcc   = (x_s4 - x_s)*K_Sam/dz_chcc ;
    double dc_h3sio4sdz_chcc   = c_h3sio4*(dc_h4sio4sdz_chcc/c_h4sio4 \
								 - dc_hsdz_chcc/c_h) ;
    double dc_cah3sio4sdz_chcc = K_cah3sio4*(dc_casdz_chcc*c_h3sio4 \
								 + c_ca*dc_h3sio4sdz_chcc) ;
    double dc_h2sio4sdz_chcc   = K_h2sio4*(dc_h3sio4sdz_chcc*c_oh	\
			                     + dc_ohsdz_chcc*c_h3sio4) ;
    double dc_cah2sio4sdz_chcc = K_cah2sio4*(dc_h2sio4sdz_chcc*c_ca	\
			                     + dc_casdz_chcc*c_h2sio4) ;
    double dc_caco3aqsdz_chcc  = K_caco3aq*(dc_co3sdz_chcc*c_ca \
			                     + dc_casdz_chcc*c_co3) ;				                  
    double dc_caohsdz_chcc     = K_caoh*(dc_casdz_chcc*c_oh \
			                     + dc_ohsdz_chcc*c_ca) ;					       

    double dn_ch_cisdz_chcc    = 0 ;
    double dn_cc_cisdz_chcc    = 0 ;


    double dn_ch_eqsdz_chcc    = (z_chcc > 0) ? n_ch_ref : 0 ;
    double dn_cc_eqsdz_chcc    = (z_chcc > 0) ? n_cc_ref : 0 ;

    double dn_chsdz_chcc       = (z_co2 <= 1) ? dn_ch_eqsdz_chcc  : dn_ch_cisdz_chcc ;
    double dn_ccsdz_chcc  	   = (z_co2 >  1) ? dn_cc_eqsdz_chcc  : dn_cc_cisdz_chcc ;    
    double dn_samsdz_chcc 	   = (x_s4 - x_s) * MAX(z_csh,0) / dz_chcc ;
    double dn_jensdz_chcc 	   = (x_jen4 - x_jen) * MAX(z_csh,0) / dz_chcc;
    double dn_tobIsdz_chcc 	   = (x_tobI4 - x_tobI) * MAX(z_csh,0) / dz_chcc;

    double dn_si_ssdz_chcc 	   = Y_Jen*dn_jensdz_chcc + Y_TobI*dn_tobIsdz_chcc + dn_samsdz_chcc ;
    double dn_ca_ssdz_chcc 	   = dn_chsdz_chcc + dn_ccsdz_chcc + X_Jen*dn_jensdz_chcc + X_TobI*dn_tobIsdz_chcc ;

    double dphisdz_chcc    	   = - v_ch*dn_chsdz_chcc - v_cc*dn_ccsdz_chcc - V_Jen*dn_jensdz_chcc - V_TobI*dn_tobIsdz_chcc - v_sam*dn_samsdz_chcc ;

    double dn_co2sdz_chcc      = dphisdz_chcc*c_co2 ;
    double dn_h2co3sdz_chcc    = dphisdz_chcc*c_h2co3 ;
    double dn_hco3sdz_chcc     = phi*dc_hco3sdz_chcc + dphisdz_chcc*c_hco3 ;
    double dn_co3sdz_chcc      = phi*dc_co3sdz_chcc + dphisdz_chcc*c_co3 ;
    double dn_casdz_chcc       = phi*dc_casdz_chcc + dphisdz_chcc*c_ca ;
    double dn_cahco3sdz_chcc   = phi*dc_cahco3sdz_chcc + dphisdz_chcc*c_cahco3 ;
    double dn_cah3sio4sdz_chcc = phi*dc_cah3sio4sdz_chcc + dphisdz_chcc*c_cah3sio4 ;
    double dn_h3sio4sdz_chcc   = phi*dc_h3sio4sdz_chcc + dphisdz_chcc*c_h3sio4 ;
    double dn_h4sio4sdz_chcc   = phi*dc_h4sio4sdz_chcc + dphisdz_chcc*c_h4sio4 ;
    double dn_h2sio4sdz_chcc   = phi*dc_h2sio4sdz_chcc + dphisdz_chcc*c_h2sio4 ;
    double dn_cah2sio4sdz_chcc = phi*dc_cah2sio4sdz_chcc + dphisdz_chcc*c_cah2sio4 ;
    double dn_caco3aqsdz_chcc  = phi*dc_caco3aqsdz_chcc + dphisdz_chcc*c_caco3aq ;
    double dn_caohsdz_chcc     = phi*dc_caohsdz_chcc + dphisdz_chcc*c_caoh ;
    double dn_ksdz_chcc        = dphisdz_chcc*c_k;
    double dn_clsdz_chcc       = dphisdz_chcc*c_cl;
    
    /* ... par rapport a c_k */
    double dc_k                = 0.4e-4 ;/*c_k*((c_k > X_Kn(i)) ? 1 : -1) ;*/
    double c_k2                = c_k + dc_k ;
    double c_oh4               = concentration_oh(c_co2,z_csh,z_chcc,c_k2,c_cl)  ;
    double dc_ohsdc_k          = (c_oh4 - c_oh)/dc_k ;
    double dc_hsdc_k           = - c_h*dc_ohsdc_k/c_oh ;
    double dc_hco3sdc_k        = c_hco3*(dc_ohsdc_k/c_oh) ;
    double dc_co3sdc_k         = c_co3*(dc_hco3sdc_k/c_hco3 + dc_ohsdc_k/c_oh) ; 
	                           
    double dc_casdc_k          = - c_ca*dc_co3sdc_k/c_co3 ;
    double dc_cahco3sdc_k      = K_cahco3*(dc_casdc_k*c_hco3 + c_ca*dc_hco3sdc_k) ;

    double dc_h3sio4sdc_k      = - c_h3sio4*dc_hsdc_k/c_h ;
    double dc_cah3sio4sdc_k    = K_cah3sio4*(dc_casdc_k*c_h3sio4 + c_ca*dc_h3sio4sdc_k) ;
    double dc_h2sio4sdc_k      = K_h2sio4*(dc_h3sio4sdc_k*c_oh + dc_ohsdc_k*c_h3sio4) ;
    double dc_cah2sio4sdc_k    = K_cah2sio4*(dc_h2sio4sdc_k*c_ca + dc_casdc_k*c_h2sio4) ;
    double dc_caco3aqsdc_k     = K_caco3aq*(dc_co3sdc_k*c_ca + dc_casdc_k*c_co3) ;				                  
    double dc_caohsdc_k        = K_caoh*(dc_casdc_k*c_oh + dc_ohsdc_k*c_ca) ;					       

    double dn_hco3sdc_k        = phi*dc_hco3sdc_k ;
    double dn_co3sdc_k         = phi*dc_co3sdc_k ;
    double dn_casdc_k          = phi*dc_casdc_k ;
    double dn_cahco3sdc_k      = phi*dc_cahco3sdc_k ;
    double dn_cah3sio4sdc_k    = phi*dc_cah3sio4sdc_k ;
    double dn_h3sio4sdc_k      = phi*dc_h3sio4sdc_k ;
    double dn_h2sio4sdc_k      = phi*dc_h2sio4sdc_k ;
    double dn_cah2sio4sdc_k    = phi*dc_cah2sio4sdc_k ;
    double dn_caco3aqsdc_k     = phi*dc_caco3aqsdc_k ;
    double dn_caohsdc_k        = phi*dc_caohsdc_k ;
    
    /* ... par rapport a c_cl */
    double dc_cl                = 0.4e-4 ;
    double c_cl2                = c_cl + dc_cl ;
    double c_oh5                = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl2) ;
    double dc_ohsdc_cl          = (c_oh5 - c_oh)/dc_cl ;
    double dc_hsdc_cl           = - c_h*dc_ohsdc_cl/c_oh ;
    double dc_hco3sdc_cl        = c_hco3*(dc_ohsdc_cl/c_oh) ;
    double dc_co3sdc_cl         = c_co3*(dc_hco3sdc_cl/c_hco3 + dc_ohsdc_cl/c_oh) ; 
	                           
    double dc_casdc_cl          = - c_ca*dc_co3sdc_cl/c_co3 ;
    double dc_cahco3sdc_cl      = K_cahco3*(dc_casdc_cl*c_hco3 + c_ca*dc_hco3sdc_cl) ;

    double dc_h3sio4sdc_cl      = - c_h3sio4*dc_hsdc_cl/c_h ;
    double dc_cah3sio4sdc_cl    = K_cah3sio4*(dc_casdc_cl*c_h3sio4 + c_ca*dc_h3sio4sdc_cl) ;
    double dc_h2sio4sdc_cl      = K_h2sio4*(dc_h3sio4sdc_cl*c_oh + dc_ohsdc_cl*c_h3sio4) ;
    double dc_cah2sio4sdc_cl    = K_cah2sio4*(dc_h2sio4sdc_cl*c_ca + dc_casdc_cl*c_h2sio4) ;
    double dc_caco3aqsdc_cl     = K_caco3aq*(dc_co3sdc_cl*c_ca + dc_casdc_cl*c_co3) ;				                  
    double dc_caohsdc_cl        = K_caoh*(dc_casdc_cl*c_oh + dc_ohsdc_cl*c_ca) ;					       

    double dn_hco3sdc_cl        = phi*dc_hco3sdc_cl ;
    double dn_co3sdc_cl         = phi*dc_co3sdc_cl ;
    double dn_casdc_cl          = phi*dc_casdc_cl ;
    double dn_cahco3sdc_cl      = phi*dc_cahco3sdc_cl ;
    double dn_cah3sio4sdc_cl    = phi*dc_cah3sio4sdc_cl ;
    double dn_h3sio4sdc_cl      = phi*dc_h3sio4sdc_cl ;
    double dn_h2sio4sdc_cl      = phi*dc_h2sio4sdc_cl ;
    double dn_cah2sio4sdc_cl    = phi*dc_cah2sio4sdc_cl ;
    double dn_caco3aqsdc_cl     = phi*dc_caco3aqsdc_cl ;
    double dn_caohsdc_cl        = phi*dc_caohsdc_cl ;

    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_CO2+j)   += volume[i]*(dn_co2sdc_co2 + dn_h2co3sdc_co2 + dn_hco3sdc_co2 + dn_co3sdc_co2 + dn_cahco3sdc_co2 \
                          + dn_caco3aqsdc_co2) ;
    K(E_C+j,I_CHCC+j)  += volume[i]*(dn_co2sdz_chcc + dn_h2co3sdz_chcc + dn_hco3sdz_chcc + dn_co3sdz_chcc + dn_cahco3sdz_chcc \
                          + dn_caco3aqsdz_chcc + dn_ccsdz_chcc) ;
    K(E_C+j,I_CSH+j)   += volume[i]*(dn_co2sdz_csh + dn_h2co3sdz_csh + dn_hco3sdz_csh + dn_co3sdz_csh + dn_cahco3sdz_csh \
                          + dn_caco3aqsdz_csh) ;
    K(E_C+j,I_K+j)     += volume[i]*(dn_hco3sdc_k + dn_co3sdc_k + dn_cahco3sdc_k + dn_caco3aqsdc_k) ; 
    K(E_C+j,I_Cl+j)    += volume[i]*(dn_hco3sdc_cl + dn_co3sdc_cl + dn_cahco3sdc_cl + dn_caco3aqsdc_cl) ;                          
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_CO2+j)  += volume[i]*(dn_casdc_co2 + dn_cahco3sdc_co2 + dn_cah3sio4sdc_co2 + dn_ca_ssdc_co2 \
                          + dn_cah2sio4sdc_co2 + dn_caco3aqsdc_co2 +dn_caohsdc_co2) ;
    K(E_Ca+j,I_CHCC+j) += volume[i]*(dn_casdz_chcc + dn_cahco3sdz_chcc + dn_cah3sio4sdz_chcc + dn_ca_ssdz_chcc \
                          + dn_cah2sio4sdz_chcc + dn_caco3aqsdz_chcc +dn_caohsdz_chcc ) ;
    K(E_Ca+j,I_CSH+j)  += volume[i]*(dn_casdz_csh + dn_cahco3sdz_csh + dn_cah3sio4sdz_csh + dn_ca_ssdz_csh \
                          + dn_cah2sio4sdz_csh + dn_caco3aqsdz_csh +dn_caohsdz_csh ) ;
    K(E_Ca+j,I_K+j)    += volume[i]*(dn_casdc_k + dn_cahco3sdc_k + dn_cah3sio4sdc_k \
                          + dn_cah2sio4sdc_k + dn_caco3aqsdc_k +dn_caohsdc_k ) ;
    K(E_Ca+j,I_Cl+j)   += volume[i]*(dn_casdc_cl + dn_cahco3sdc_cl + dn_cah3sio4sdc_cl \
                          + dn_cah2sio4sdc_cl + dn_caco3aqsdc_cl +dn_caohsdc_cl ) ;
    /*
      Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
    */
    K(E_Si+j,I_CO2+j)  += volume[i]*(dn_h3sio4sdc_co2 + dn_h4sio4sdc_co2 + dn_cah3sio4sdc_co2 + dn_si_ssdc_co2\
                          + dn_h2sio4sdc_co2 + dn_cah2sio4sdc_co2) ;
    K(E_Si+j,I_CHCC+j) += volume[i]*(dn_h3sio4sdz_chcc + dn_h4sio4sdz_chcc + dn_cah3sio4sdz_chcc + dn_si_ssdz_chcc \
                          + dn_h2sio4sdz_chcc + dn_cah2sio4sdz_chcc) ;
    K(E_Si+j,I_CSH+j)  += volume[i]*(dn_h3sio4sdz_csh + dn_h4sio4sdz_csh + dn_cah3sio4sdz_csh + dn_si_ssdz_csh \
                          + dn_h2sio4sdz_csh + dn_cah2sio4sdz_csh) ;
    K(E_Si+j,I_K+j)    += volume[i]*(dn_h3sio4sdc_k + dn_cah3sio4sdc_k + dn_h2sio4sdc_k + dn_cah2sio4sdc_k) ;
    K(E_Si+j,I_Cl+j)   += volume[i]*(dn_h3sio4sdc_cl + dn_cah3sio4sdc_cl + dn_h2sio4sdc_cl + dn_cah2sio4sdc_cl) ;    
    /*
      Conservation de la charge  : div(w_q) = 0
    */
    /*
      Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(E_K+j,I_CO2+j)   += volume[i]*(dn_ksdc_co2) ; 
    K(E_K+j,I_CHCC+j)  += volume[i]*(dn_ksdz_chcc) ; 
    K(E_K+j,I_CSH+j)   += volume[i]*(dn_ksdz_csh) ;    
    K(E_K+j,I_K+j)     += volume[i]*(phi*1.) ;
    /*      Conservation de Cl (chloride)  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
    K(E_Cl+j,I_CO2+j)   += volume[i]*(dn_clsdc_co2) ; 
    K(E_Cl+j,I_CHCC+j)  += volume[i]*(dn_clsdz_chcc) ; 
    K(E_Cl+j,I_CSH+j)   += volume[i]*(dn_clsdz_csh) ;    
    K(E_Cl+j,I_Cl+j)    += volume[i]*(phi*1.) ;

    /* sauvegardes pour les termes de transport */
    Dc_hSDc_co2[i]        = dc_hsdc_co2 ;
    Dc_ohSDc_co2[i]       = dc_ohsdc_co2 ;
    Dc_h2co3SDc_co2[i]    = dc_h2co3sdc_co2 ;
    Dc_hco3SDc_co2[i]     = dc_hco3sdc_co2 ;
    Dc_co3SDc_co2[i]      = dc_co3sdc_co2 ;
    Dc_caSDc_co2[i]       = dc_casdc_co2 ;
    Dc_cahco3SDc_co2[i]   = dc_cahco3sdc_co2 ;
    Dc_h3sio4SDc_co2[i]   = dc_h3sio4sdc_co2 ;
    Dc_h4sio4SDc_co2[i]   = dc_h4sio4sdc_co2 ;
    Dc_cah3sio4SDc_co2[i] = dc_cah3sio4sdc_co2 ;
    Dc_h2sio4SDc_co2[i]   = dc_h2sio4sdc_co2 ;
    Dc_cah2sio4SDc_co2[i] = dc_cah2sio4sdc_co2 ;
    Dc_caco3aqSDc_co2[i]  = dc_caco3aqsdc_co2 ;
    Dc_caohSDc_co2[i]     = dc_caohsdc_co2 ;
    
    Dc_hSDz_csh[i]        = dc_hsdz_csh ;
    Dc_ohSDz_csh[i]       = dc_ohsdz_csh ;
    Dc_hco3SDz_csh[i]     = dc_hco3sdz_csh ;
    Dc_co3SDz_csh[i]      = dc_co3sdz_csh ;   
    Dc_caSDz_csh[i]       = dc_casdz_csh ;
    Dc_cahco3SDz_csh[i]   = dc_cahco3sdz_csh ;
    Dc_h3sio4SDz_csh[i]   = dc_h3sio4sdz_csh ;
    Dc_h4sio4SDz_csh[i]   = dc_h4sio4sdz_csh ;
    Dc_cah3sio4SDz_csh[i] = dc_cah3sio4sdz_csh ;
    Dc_h2sio4SDz_csh[i]   = dc_h2sio4sdz_csh ;
    Dc_cah2sio4SDz_csh[i] = dc_cah2sio4sdz_csh ;
    Dc_caco3aqSDz_csh[i]  = dc_caco3aqsdz_csh ;
    Dc_caohSDz_csh[i]     = dc_caohsdz_csh ;  

    Dc_hSDz_chcc[i]        = dc_hsdz_chcc ;
    Dc_ohSDz_chcc[i]       = dc_ohsdz_chcc ;
    Dc_hco3SDz_chcc[i]     = dc_hco3sdz_chcc ;
    Dc_co3SDz_chcc[i]      = dc_co3sdz_chcc ;
    Dc_caSDz_chcc[i]       = dc_casdz_chcc ;      
    Dc_cahco3SDz_chcc[i]   = dc_cahco3sdz_chcc ;
    Dc_h3sio4SDz_chcc[i]   = dc_h3sio4sdz_chcc ;
    Dc_h4sio4SDz_chcc[i]   = dc_h4sio4sdz_chcc ;
    Dc_cah3sio4SDz_chcc[i] = dc_cah3sio4sdz_chcc ;
    Dc_h2sio4SDz_chcc[i]   = dc_h2sio4sdz_chcc ;
    Dc_cah2sio4SDz_chcc[i] = dc_cah2sio4sdz_chcc ;
    Dc_caco3aqSDz_chcc[i]  = dc_caco3aqsdz_chcc ;
    Dc_caohSDz_chcc[i]     = dc_caohsdz_chcc ;
    
    Dc_hSDc_k[i]        = dc_hsdc_k ;
    Dc_ohSDc_k[i]       = dc_ohsdc_k ;
    Dc_hco3SDc_k[i]     = dc_hco3sdc_k ;
    Dc_co3SDc_k[i]      = dc_co3sdc_k ;
    Dc_caSDc_k[i]       = dc_casdc_k ;      
    Dc_cahco3SDc_k[i]   = dc_cahco3sdc_k ;
    Dc_h3sio4SDc_k[i]   = dc_h3sio4sdc_k ;
    Dc_cah3sio4SDc_k[i] = dc_cah3sio4sdc_k ;
    Dc_h2sio4SDc_k[i]   = dc_h2sio4sdc_k ;
    Dc_cah2sio4SDc_k[i] = dc_cah2sio4sdc_k ;
    Dc_caco3aqSDc_k[i]  = dc_caco3aqsdc_k ;
    Dc_caohSDc_k[i]     = dc_caohsdc_k ;
    
    Dc_hSDc_cl[i]        = dc_hsdc_cl ;
    Dc_ohSDc_cl[i]       = dc_ohsdc_cl ;
    Dc_hco3SDc_cl[i]     = dc_hco3sdc_cl ;
    Dc_co3SDc_cl[i]      = dc_co3sdc_cl ;
    Dc_caSDc_cl[i]       = dc_casdc_cl ;      
    Dc_cahco3SDc_cl[i]   = dc_cahco3sdc_cl ;
    Dc_h3sio4SDc_cl[i]   = dc_h3sio4sdc_cl ;
    Dc_cah3sio4SDc_cl[i] = dc_cah3sio4sdc_cl ;
    Dc_h2sio4SDc_cl[i]   = dc_h2sio4sdc_cl ;
    Dc_cah2sio4SDc_cl[i] = dc_cah2sio4sdc_cl ;
    Dc_caco3aqSDc_cl[i]  = dc_caco3aqsdc_cl ;
    Dc_caohSDc_cl[i]     = dc_caohsdc_cl ;     
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
  double trf_cahco3   = tr*KF_CaHCO3 ;
  double trf_cah3sio4 = tr*KF_CaH3SiO4 ;
  double trf_h3sio4   = tr*KF_H3SiO4 ;
  double trf_h4sio4   = tr*KF_H4SiO4 ;
  double trf_h2sio4   = tr*KF_H2SiO4 ;
  double trf_cah2sio4 = tr*KF_CaH2SiO4 ;
  double trf_caco3aq  = tr*KF_CaCO3aq ;
  double trf_caoh     = tr*KF_CaOH ;
  double trf_k        = tr*KF_K ;                   /* ajouter2 */
  double trf_cl       = tr*KF_Cl ;  
  
  double tre_hco3     = tr*Kpsi_HCO3 ;
  double tre_co3      = tr*Kpsi_CO3 ;
  double tre_ca       = tr*Kpsi_Ca ;
  double tre_cahco3   = tr*Kpsi_CaHCO3 ;
  double tre_cah3sio4 = tr*Kpsi_CaH3SiO4 ;
  double tre_h3sio4   = tr*Kpsi_H3SiO4 ;
  double tre_h2sio4   = tr*Kpsi_H2SiO4 ;
  double tre_caoh     = tr*Kpsi_CaOH ;
  double tre_k        = tr*Kpsi_K ;                  /* ajouter2 */ 
  double tre_cl       = tr*Kpsi_Cl ;    

  double tre_q        = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_co2 + trf_h2co3*Dc_h2co3SDc_co2[i] + trf_hco3*Dc_hco3SDc_co2[i] + trf_co3*Dc_co3SDc_co2[i] + trf_cahco3*Dc_cahco3SDc_co2[i] \
           + trf_caco3aq*Dc_caco3aqSDc_co2[i] ;
  }
  K(E_C,I_CO2)          += + c[0] ;
  K(E_C,I_CO2+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CO2)      += - c[0] ;
  K(E_C+NEQ,I_CO2+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDz_csh[i] + trf_co3*Dc_co3SDz_csh[i] + trf_cahco3*Dc_cahco3SDz_csh[i] + trf_caco3aq*Dc_caco3aqSDz_csh[i] ;
  }
  K(E_C,I_CSH)          += + c[0] ;
  K(E_C,I_CSH+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CSH)      += - c[0] ;
  K(E_C+NEQ,I_CSH+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDz_chcc[i] + trf_co3*Dc_co3SDz_chcc[i] + trf_cahco3*Dc_cahco3SDz_chcc[i] + trf_caco3aq*Dc_caco3aqSDz_chcc[i] ;
  }
  K(E_C,I_CHCC)          += + c[0] ;
  K(E_C,I_CHCC+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CHCC)      += - c[0] ;
  K(E_C+NEQ,I_CHCC+NEQ)  += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDc_k[i] + trf_co3*Dc_co3SDc_k[i] + trf_cahco3*Dc_cahco3SDc_k[i] + trf_caco3aq*Dc_caco3aqSDc_k[i] ;
  }
  K(E_C,I_K)          += + c[0] ;
  K(E_C,I_K+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_K)      += - c[0] ;
  K(E_C+NEQ,I_K+NEQ)  += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDc_cl[i] + trf_co3*Dc_co3SDc_cl[i] + trf_cahco3*Dc_cahco3SDc_cl[i] + trf_caco3aq*Dc_caco3aqSDc_cl[i] ;
  }
  K(E_C,I_Cl)          += + c[0] ;
  K(E_C,I_Cl+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_Cl)      += - c[0] ;
  K(E_C+NEQ,I_Cl+NEQ)  += + c[1] ;

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
    c[i] = trf_ca*Dc_caSDc_co2[i] + trf_cahco3*Dc_cahco3SDc_co2[i] + trf_cah3sio4*Dc_cah3sio4SDc_co2[i] + trf_cah2sio4*Dc_cah2sio4SDc_co2[i] \
    + trf_caoh*Dc_caohSDc_co2[i] + trf_caco3aq*Dc_caco3aqSDc_co2[i] ;
  }
  K(E_Ca,I_CO2)         += + c[0] ;
  K(E_Ca,I_CO2+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CO2)     += - c[0] ;
  K(E_Ca+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDz_csh[i] + trf_cahco3*Dc_cahco3SDz_csh[i] + trf_cah3sio4*Dc_cah3sio4SDz_csh[i] + trf_cah2sio4*Dc_cah2sio4SDz_csh[i] \
    + trf_caoh*Dc_caohSDz_csh[i] + trf_caco3aq*Dc_caco3aqSDz_csh[i] ;
  }
  K(E_Ca,I_CSH)         += + c[0] ;
  K(E_Ca,I_CSH+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CSH)     += - c[0] ;
  K(E_Ca+NEQ,I_CSH+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDz_chcc[i] + trf_cahco3*Dc_cahco3SDz_chcc[i] + trf_cah3sio4*Dc_cah3sio4SDz_chcc[i] + trf_cah2sio4*Dc_cah2sio4SDz_chcc[i] \
    + trf_caoh*Dc_caohSDz_chcc[i] + trf_caco3aq*Dc_caco3aqSDz_chcc[i] ;
  }
  K(E_Ca,I_CHCC)         += + c[0] ;
  K(E_Ca,I_CHCC+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CHCC)     += - c[0] ;
  K(E_Ca+NEQ,I_CHCC+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_k[i] + trf_cahco3*Dc_cahco3SDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i] + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_caoh*Dc_caohSDc_k[i] + trf_caco3aq*Dc_caco3aqSDc_k[i] ;
  }
  K(E_Ca,I_K)         += + c[0] ;
  K(E_Ca,I_K+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_K)     += - c[0] ;
  K(E_Ca+NEQ,I_K+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDc_cl[i] + trf_cahco3*Dc_cahco3SDc_cl[i] + trf_cah3sio4*Dc_cah3sio4SDc_cl[i] + trf_cah2sio4*Dc_cah2sio4SDc_cl[i] \
    + trf_caoh*Dc_caohSDc_cl[i] + trf_caco3aq*Dc_caco3aqSDc_cl[i] ;
  }
  K(E_Ca,I_Cl)         += + c[0] ;
  K(E_Ca,I_Cl+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Cl)     += - c[0] ;
  K(E_Ca+NEQ,I_Cl+NEQ) += + c[1] ;   

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_cahco3 + tre_cah3sio4 + tre_caoh;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_co2[i] + trf_h4sio4*Dc_h4sio4SDc_co2[i] + trf_cah3sio4*Dc_cah3sio4SDc_co2[i] + trf_cah2sio4*Dc_cah2sio4SDc_co2[i] \
    + trf_h2sio4*Dc_h2sio4SDc_co2[i] ;
  }
  K(E_Si,I_CO2)         += + c[0] ;
  K(E_Si,I_CO2+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CO2)     += - c[0] ;
  K(E_Si+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDz_csh[i] + trf_h4sio4*Dc_h4sio4SDz_csh[i] + trf_cah3sio4*Dc_cah3sio4SDz_csh[i]  + trf_cah2sio4*Dc_cah2sio4SDz_csh[i] \
    + trf_h2sio4*Dc_h2sio4SDz_csh[i];
  }
  K(E_Si,I_CSH)         += + c[0] ;
  K(E_Si,I_CSH+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CSH)     += - c[0] ;
  K(E_Si+NEQ,I_CSH+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDz_chcc[i] + trf_h4sio4*Dc_h4sio4SDz_chcc[i] + trf_cah3sio4*Dc_cah3sio4SDz_chcc[i]  + trf_cah2sio4*Dc_cah2sio4SDz_chcc[i] \
    + trf_h2sio4*Dc_h2sio4SDz_chcc[i];
  }
  K(E_Si,I_CHCC)         += + c[0] ;
  K(E_Si,I_CHCC+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CHCC)     += - c[0] ;
  K(E_Si+NEQ,I_CHCC+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_k[i] + trf_cah3sio4*Dc_cah3sio4SDc_k[i]  + trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_Si,I_K)         += + c[0] ;
  K(E_Si,I_K+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_K)     += - c[0] ;
  K(E_Si+NEQ,I_K+NEQ) += + c[1] ; 
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDc_cl[i] + trf_cah3sio4*Dc_cah3sio4SDc_cl[i]  + trf_cah2sio4*Dc_cah2sio4SDc_cl[i] \
    + trf_h2sio4*Dc_h2sio4SDc_cl[i];
  }
  K(E_Si,I_Cl)         += + c[0] ;
  K(E_Si,I_Cl+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_Cl)     += - c[0] ;
  K(E_Si+NEQ,I_Cl+NEQ) += + c[1] ;   

  for(i=0;i<2;i++){
    c[i] = tre_h3sio4 + tre_cah3sio4 + tre_h2sio4 ;
  }
  K(E_Si,I_psi)          += + c[0] ;
  K(E_Si,I_psi+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_psi)      += - c[0] ;
  K(E_Si+NEQ,I_psi+NEQ)  += + c[1] ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDc_co2[i] + z_oh*trf_oh*Dc_ohSDc_co2[i] + z_hco3*trf_hco3*Dc_hco3SDc_co2[i] + z_co3*trf_co3*Dc_co3SDc_co2[i] \
           + z_ca*trf_ca*Dc_caSDc_co2[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_co2[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_co2[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_co2[i] + z_caoh*trf_caoh*Dc_caohSDc_co2[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_co2[i];
  }
  K(E_q,I_CO2)           += + c[0] ;
  K(E_q,I_CO2+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CO2)       += - c[0] ;
  K(E_q+NEQ,I_CO2+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDz_csh[i] + z_oh*trf_oh*Dc_ohSDz_csh[i] + z_hco3*trf_hco3*Dc_hco3SDz_csh[i] + z_co3*trf_co3*Dc_co3SDz_csh[i] \
           + z_ca*trf_ca*Dc_caSDz_csh[i] + z_cahco3*trf_cahco3*Dc_cahco3SDz_csh[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDz_csh[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDz_csh[i] + z_caoh*trf_caoh*Dc_caohSDz_csh[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDz_csh[i];
  }
  K(E_q,I_CSH)           += + c[0] ;
  K(E_q,I_CSH+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CSH)       += - c[0] ;
  K(E_q+NEQ,I_CSH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDz_chcc[i] + z_oh*trf_oh*Dc_ohSDz_chcc[i] + z_hco3*trf_hco3*Dc_hco3SDz_chcc[i] + z_co3*trf_co3*Dc_co3SDz_chcc[i] \
           + z_ca*trf_ca*Dc_caSDz_chcc[i] + z_cahco3*trf_cahco3*Dc_cahco3SDz_chcc[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDz_chcc[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDz_chcc[i] + z_caoh*trf_caoh*Dc_caohSDz_chcc[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDz_chcc[i];
  }
  K(E_q,I_CHCC)           += + c[0] ;
  K(E_q,I_CHCC+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CHCC)       += - c[0] ;
  K(E_q+NEQ,I_CHCC+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_k*trf_k + z_h*trf_h*Dc_hSDc_k[i] + z_oh*trf_oh*Dc_ohSDc_k[i] + z_hco3*trf_hco3*Dc_hco3SDc_k[i] + z_co3*trf_co3*Dc_co3SDc_k[i] \
           + z_ca*trf_ca*Dc_caSDc_k[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_k[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_k[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_k[i] + z_caoh*trf_caoh*Dc_caohSDc_k[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_k[i];
  }
  K(E_q,I_K)           += + c[0] ;
  K(E_q,I_K+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_K)       += - c[0] ;
  K(E_q+NEQ,I_K+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_cl*trf_cl + z_h*trf_h*Dc_hSDc_cl[i] + z_oh*trf_oh*Dc_ohSDc_cl[i] + z_hco3*trf_hco3*Dc_hco3SDc_cl[i] + z_co3*trf_co3*Dc_co3SDc_cl[i] \
           + z_ca*trf_ca*Dc_caSDc_cl[i] + z_cahco3*trf_cahco3*Dc_cahco3SDc_cl[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDc_cl[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDc_cl[i] + z_caoh*trf_caoh*Dc_caohSDc_cl[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDc_cl[i];
  }
  K(E_q,I_Cl)           += + c[0] ;
  K(E_q,I_Cl+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Cl)       += - c[0] ;
  K(E_q+NEQ,I_Cl+NEQ)   += + c[1] ;     

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  
    /*
  Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  for(i=0;i<2;i++){
    c[i] = tre_k ;
  }
  K(E_K,I_psi)          += + c[0] ;
  K(E_K,I_psi+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_psi)      += - c[0] ;
  K(E_K+NEQ,I_psi+NEQ)  += + c[1] ;
    
   
  for(i=0;i<2;i++){
    c[i] = trf_k ;
  }
  K(E_K,I_K)          += + c[0] ;
  K(E_K,I_K+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_K)      += - c[0] ;
  K(E_K+NEQ,I_K+NEQ)  += + c[1] ; 
  
  /*
  Conservation de Cl (chloride)  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  for(i=0;i<2;i++){
    c[i] = tre_cl ;
  }
  K(E_Cl,I_psi)          += + c[0] ;
  K(E_Cl,I_psi+NEQ)      += - c[1] ;
  K(E_Cl+NEQ,I_psi)      += - c[0] ;
  K(E_Cl+NEQ,I_psi+NEQ)  += + c[1] ;
    
   
  for(i=0;i<2;i++){
    c[i] = trf_cl ;
  }
  K(E_Cl,I_Cl)          += + c[0] ;
  K(E_Cl,I_Cl+NEQ)      += - c[1] ;
  K(E_Cl+NEQ,I_Cl)      += - c[0] ;
  K(E_Cl+NEQ,I_Cl+NEQ)  += + c[1] ;   
  
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


void rs55(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
      Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ; 
  /*
  Conservation de Cl (chloride)  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ;  
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;

#undef R
}


int so55(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
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
  nso = 30 ;
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
    double z_chcc     = param(u,h_s,el.nn,I_CHCC) ;
    double c_k        = param(u,h_s,el.nn,I_K) ;
    double z_csh      = param(u,h_s,el.nn,I_CSH);
    double c_cl       = param(u,h_s,el.nn,I_Cl) ;
    
    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;           
    double pk_ch      = (MIN(z_chcc,0) - log(MAX(z_co2,1.)))/log(10); /*log10(P_CH/K_CH);  ajouter3 */
    double x_s        = product_Sam(c_co2,z_chcc,z_csh);
    
    double c_h4sio4   = product_Sam(c_co2,z_chcc,z_csh) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl)  ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;    
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */

    /* densite de charge */
    double c_q = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 \
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k ;
  /* contenus solides */
    int j = (el.dim < dim) ? 0 : ((s[0] < (x[0][0] + x[1][0])*0.5) ? 0 : 1) ;
    double n_ch       = N_CH(j) ;
    double n_cc       = N_CC(j) ;
    double n_jen      = N_Jen(j)   ;
    double n_sam      = N_Sam(j)   ;
    double n_tobI     = N_TobI(j)   ;
  /* porosite */
    double n_ch0      = N_CH0(j) ;
    double n_cc0      = N_CC0(j) ;
    double n_jen0     = N_Jen0(j)   ;
    double n_sam0     = N_Sam0(j)   ;
    double n_tobI0    = N_TobI0(j)   ;
    double phi        = phi0 + v_ch*(n_ch0 - n_ch) + V_Jen*(n_jen0- n_jen) + V_TobI*(n_tobI0- n_tobI) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
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
    strcpy(r[i].text,"n_jen") ; r[i].n = 1 ;
    r[i++].v[0] = n_jen ;
    strcpy(r[i].text,"n_tobI") ; r[i].n = 1 ;
    r[i++].v[0] = n_tobI ;
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
    strcpy(r[i].text,"z_co2") ; r[i].n = 1 ;
    r[i++].v[0] = z_co2 ;
    strcpy(r[i].text,"c_caco3aq") ; r[i].n = 1 ;
    r[i++].v[0] = c_caco3aq ;
    strcpy(r[i].text,"c_caoh") ; r[i].n = 1 ;
    r[i++].v[0] = c_caoh ;
    strcpy(r[i].text,"c_h2sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2sio4 ;   
    strcpy(r[i].text,"c_cah2sio4") ; r[i].n = 1 ;
    r[i++].v[0] = c_cah2sio4 ;
    strcpy(r[i].text,"c_k") ; r[i].n = 1 ;
    r[i++].v[0] = c_k ; 
    strcpy(r[i].text,"pk_ch") ; r[i].n = 1 ;
    r[i++].v[0] = pk_ch ;
    strcpy(r[i].text,"x_s") ; r[i].n = 1 ;
    r[i++].v[0] = x_s ;  
    strcpy(r[i].text,"z_csh") ; r[i].n = 1 ;
    r[i++].v[0] = z_csh ; 
    strcpy(r[i].text,"c_cl") ; r[i].n = 1 ;
    r[i++].v[0] = c_cl ;   
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
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;
    double c_k        = X_K(i);
    double z_csh      = Z_CSH(i);
    double c_cl       = X_Cl(i);    

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ; 

    double c_h4sio4   = product_Sam(c_co2,z_chcc,z_csh) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */

    /* solides */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_jen      = N_Jen(i) ;
    double n_tobI     = N_TobI(i) ;
    double n_sam      = N_Sam(i) ;

    /* porosite */
    double n_cc0      = N_CC0(i) ;
    double n_ch0      = N_CH0(i) ;
    double n_jen0     = N_Jen0(i) ;
    double n_tobI0    = N_TobI0(i) ;
    double n_sam0     = N_Sam0(i) ;
    double phi        = phi0 + v_ch*(n_ch0 - n_ch) + V_Jen*(n_jen0- n_jen) + V_TobI*(n_tobI0- n_tobI) + v_sam*(n_sam0 - n_sam) + v_cc*(n_cc0 - n_cc) ;
    /* tortuosite liquide */
    double iff    = phi*(0.0004 + 0.03*phi*phi + 1.7*(phi-0.17)*(phi-0.17)*(((phi-0.17)>0)? 1:0));
    
    KF_OH         += d_oh*iff ;
    KF_H          += d_h*iff ;
    KF_CO2        += d_co2*iff ;
    KF_H2CO3      += d_h2co3*iff ;
    KF_HCO3       += d_hco3*iff ;
    KF_CO3        += d_co3*iff ;

    KF_Ca         += d_ca*iff ;
    KF_CaHCO3     += d_cahco3*iff ;
    KF_CaH3SiO4   += d_cah3sio4*iff ;

    KF_H3SiO4     += d_h3sio4*iff ;
    KF_H4SiO4     += d_h4sio4*iff ;
    KF_H2SiO4     += d_h2sio4*iff ;
    KF_CaH2SiO4   += d_cah2sio4*iff ;
    KF_CaCO3aq    += d_caco3aq*iff ;
    KF_CaOH       += d_caoh*iff ;
    
    KF_K          += d_k*iff;
    KF_Cl         += d_cl*iff;

    
    Kpsi_H        += FsRT*KF_H*z_h*c_h ;
    Kpsi_OH       += FsRT*KF_OH*z_oh*c_oh ;
    Kpsi_HCO3     += FsRT*KF_HCO3*z_hco3*c_hco3 ;
    Kpsi_CO3      += FsRT*KF_CO3*z_co3*c_co3 ;

    Kpsi_Ca       += FsRT*KF_Ca*z_ca*c_ca ;
    Kpsi_CaHCO3   += FsRT*KF_CaHCO3*z_cahco3*c_cahco3 ;
    Kpsi_CaH3SiO4 += FsRT*KF_CaH3SiO4*z_cah3sio4*c_cah3sio4 ;
    Kpsi_H3SiO4   += FsRT*KF_H3SiO4*z_h3sio4*c_h3sio4 ;
    Kpsi_H2SiO4   += FsRT*KF_H2SiO4*z_h2sio4*c_h2sio4 ;
    Kpsi_CaOH     += FsRT*KF_CaOH*z_caoh*c_caoh ;
    
    Kpsi_K        += FsRT*KF_K*z_k*c_k ;
    Kpsi_Cl       += FsRT*KF_Cl*z_cl*c_cl ;

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_cahco3*Kpsi_CaHCO3 + z_h3sio4*Kpsi_H3SiO4 \
                     + z_cah3sio4*Kpsi_CaH3SiO4 + z_caoh*Kpsi_CaOH + z_h2sio4*Kpsi_H2SiO4 + z_k*Kpsi_K + z_cl*Kpsi_Cl;
  }
  
  /* moyenne */
  for(i=0;i<NVE_TR;i++) va[i] *= 0.5 ;
}


void flux(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Les flux (f) */
{
  double r_h[2],r_oh[2] ;
  double r_co2[2],r_h2co3[2],r_hco3[2],r_co3[2],r_caco3aq[2],r_caoh[2] ;
  double r_ca[2],r_cahco3[2],r_h3sio4[2],r_h4sio4[2],r_cah3sio4[2],r_cah2sio4[2],r_h2sio4[2] ;
  double r_k[2],r_cl[2] ;

  int    i ;

  for(i=0;i<2;i++) {
    double c_co2      = C_CO2(i) ;
    double z_chcc     = Z_CHCC(i) ;
    double z_co2      = c_co2/c_co2_eq ;
    double c_k        = X_K(i);                 /* ajouter2 */
    double z_csh      = Z_CSH(i);
    double c_cl       = X_Cl(i);  

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;    
 
    double c_h4sio4   = product_Sam(c_co2,z_chcc,z_csh) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,z_csh,z_chcc,c_k,c_cl) ;
    double c_h        = K_h2o/c_oh ;
    double c_h2co3    = K_h2co3*c_co2 ;
    double c_hco3     = K_hco3*c_h2co3/c_h ;
    double c_co3      = K_co3*c_hco3/c_h ;
    double c_ca       = P_CC/c_co3 ;
    double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
    double c_h3sio4   = c_h4sio4/(K_h4sio4*c_h) ;
    double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;               /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    double c_caoh     = K_caoh*c_ca*c_oh ;                    /* ajouter */

    r_co2[i]      = c_co2 ;
    r_oh[i]       = c_oh ;
    r_h[i]        = c_h ;
    r_h2co3[i]    = c_h2co3 ;
    r_hco3[i]     = c_hco3 ;
    r_co3[i]      = c_co3 ;
    r_ca[i]       = c_ca ;
    r_cahco3[i]   = c_cahco3 ;
    r_h3sio4[i]   = c_h3sio4 ;
    r_h4sio4[i]   = c_h4sio4 ;
    r_cah3sio4[i] = c_cah3sio4 ;
    r_h2sio4[i]   = c_h2sio4 ;
    r_cah2sio4[i] = c_cah2sio4 ;
    r_caoh[i]     = c_caoh ;
    r_caco3aq[i]  = c_caco3aq ;
    r_k[i]        = c_k ;
    r_cl[i]       = c_cl ;
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
    double grd_cahco3   = (r_cahco3[1]   - r_cahco3[0]  )/dx ;
    double grd_ca       = (r_ca[1]       - r_ca[0]      )/dx ;
    double grd_cah3sio4 = (r_cah3sio4[1] - r_cah3sio4[0])/dx ;
    double grd_h3sio4   = (r_h3sio4[1]   - r_h3sio4[0]  )/dx ;
    double grd_h4sio4   = (r_h4sio4[1]   - r_h4sio4[0]  )/dx ;
    double grd_h2sio4   = (r_h2sio4[1]   - r_h2sio4[0]  )/dx ;
    double grd_cah2sio4 = (r_cah2sio4[1] - r_cah2sio4[0])/dx ;
    double grd_caco3aq  = (r_caco3aq[1]  - r_caco3aq[0] )/dx ;
    double grd_caoh     = (r_caoh[1]     - r_caoh[0]    )/dx ;
    double grd_k        = (r_k[1]        - r_k[0]       )/dx ;
    double grd_cl       = (r_cl[1]       - r_cl[0]       )/dx ;
    
    double grd_psi      = (PSI(1)        - PSI(0)       )/dx ;
    
    /* Flux */
    double w_co2      = - KF_CO2*grd_co2       ; 
    double w_h2co3    = - KF_H2CO3*grd_h2co3   ;
    double w_hco3     = - KF_HCO3*grd_hco3          - Kpsi_HCO3*grd_psi     ;
    double w_co3      = - KF_CO3*grd_co3            - Kpsi_CO3*grd_psi      ;
    double w_cahco3   = - KF_CaHCO3*grd_cahco3      - Kpsi_CaHCO3*grd_psi   ;
    double w_ca       = - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi       ;
    double w_cah3sio4 = - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    double w_h3sio4   = - KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi   ;
    double w_h2sio4   = - KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi   ;
    double w_caoh     = - KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi     ;
    double w_h4sio4   = - KF_H4SiO4*grd_h4sio4 ;
    double w_cah2sio4 = - KF_CaH2SiO4*grd_cah2sio4 ;
    double w_caco3aq  = - KF_CaCO3aq*grd_caco3aq ;    
    double w_k        = - KF_K*grd_k                - Kpsi_K*grd_psi        ; 
    double w_cl       = - KF_Cl*grd_cl              - Kpsi_Cl*grd_psi        ; 
    
    double w_q        = - z_h*KF_H*grd_h		      \
                        - z_oh*KF_OH*grd_oh		      \
                        - z_hco3*KF_HCO3*grd_hco3             \
                        - z_co3*KF_CO3*grd_co3		      \
                        - z_ca*KF_Ca*grd_ca		      \
                        - z_cahco3*KF_CaHCO3*grd_cahco3	      \
                        - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                        - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                        - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                        - z_caoh*KF_CaOH*grd_caoh \
                        - z_k*KF_K*grd_k \
                        - z_cl*KF_Cl*grd_cl \
                        - Kpsi_q*grd_psi ;


    W_C     = w_co2 + w_h2co3 + w_hco3 + w_co3 + w_cahco3 + w_caco3aq ;
    W_Ca    = w_ca + w_cahco3 + w_cah3sio4 + w_caco3aq + w_caoh + w_cah2sio4 ;
    W_Si    = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
    W_q     = w_q ;
    W_K     = w_k ;
    W_Cl    = w_cl ;
  }
}


double concentration_oh(double c_co2,double z_csh,double z_chcc,double c_k,double c_cl)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_co2, z_csh et z_chcc sont fixes */
  double c_co2_eq   = C_CO2_eq ;
  double z_co2      = c_co2/c_co2_eq ;
  /* les produits sont donc aussi fixes */
  double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(z_chcc) ;
  
  double c_h4sio4   = product_Sam(c_co2,z_chcc,z_csh) * K_Sam ;

  /* double P_Sam      = K_Sam*pow(MAX(z_ch,1.),-X_CSH/Y_CSH)*pow(NEGEXP(z_csh),1/Y_CSH) ;
  double P_CSH      = pow(P_CH,X_CSH)*pow(P_Sam,Y_CSH) ;
  rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                     : +1
     c_h2co3    = K_h2co3*c_co2                  :  0
     c_hco3     = K_hco3*c_h2co3/c_h             : -1
     c_co3      = K_co3*c_hco3/c_h               : -2
     c_ca       = P_CC/c_co3                     : +2
     c_cahco3   = K_cahco3*c_ca*c_hco3           : +1
     c_h4sio4   = P_Sam                          :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_caco3aq  = K_caco3aq*c_ca*c_co3 ;         :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;        : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;     :  0      
     c_caoh     = K_caoh*c_ca*c_oh ;             : +1       
  */

  double A_hco3     = K_hco3*K_h2co3*c_co2 ;
  double A_co3      = K_co3*A_hco3 ;
  double A_ca       = P_CC/A_co3 ;
  double A_cahco3   = K_cahco3*A_ca*A_hco3 ;
  double A_h3sio4   = c_h4sio4/K_h4sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahco3*A_cahco3 + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k + z_cl*c_cl ; /* 0.04 ;*/
  double d = z_oh*K_h2o + z_hco3*A_hco3 + z_h3sio4*A_h3sio4 ;
  double e = z_co3*A_co3 + z_h2sio4*A_h2sio4 ;

  double c_h = poly4(a,b,c,d,e) ;
 
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
      printf("x0 = %e\n",x0) ;
      printf("x  = %e\n",x) ;
      arret("poly_sansc : non convergence") ;
    }
  } while(err > tol) ;
  return(x) ;
}



double product_Sam(double c_co2,double z_chcc,double z_csh)
{
  double z_co2 =  c_co2/c_co2_eq ;
  double X_CH = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
  double A_Jen = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
  double A_TobI = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
  double sumx  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1 ;
  return(NEGEXP(z_csh)/sumx) ;
}

double product_Jen(double c_co2,double z_chcc,double z_csh)
{
  double z_co2 =  c_co2/c_co2_eq ;
  double X_CH = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
  double A_Jen = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
  double A_TobI = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
  double sumx  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1 ;
  return(A_Jen*pow(X_CH,X_Jen)*NEGEXP(z_csh)/sumx) ;
}

double product_TobI(double c_co2,double z_chcc,double z_csh)
{
  double z_co2 =  c_co2/c_co2_eq ;
  double X_CH = NEGEXP(z_chcc)/MAX(z_co2,1.) ;
  double A_Jen = pow(K_CH,X_Jen)*pow(K_Sam,Y_Jen)/K_Jen ;
  double A_TobI = pow(K_CH,X_TobI)*pow(K_Sam,Y_TobI)/K_TobI ;
  double sumx  = A_Jen*pow(X_CH,X_Jen) + A_TobI*pow(X_CH,X_TobI) + 1 ;
  return(A_TobI*pow(X_CH,X_TobI)*NEGEXP(z_csh)/sumx) ;
}




