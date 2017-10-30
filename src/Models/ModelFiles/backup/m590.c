/*
courbes for csh
introduce total mass
introduce scco2
introduce water in scCO2
with real co2 transport properties
with real co2 concentration
transport parameters as a function of (P,T) 
able to calculate with alkalis and cloride 
diffusion coeffition as a function of T
fix a bug for z_ca_s
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  590
#define TITLE   "Carbonatation du ciment sature en condition SC (mai 2010)"
#define AUTHORS "Shen"

#include "OldMethods.h"


/* Macros */
#define NEQ     (7)
#define NVE     (58)
#define NVI     (27)
#define NVE_TR  (50)

#define E_C     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)
#define E_Cl    (5)
#define E_mass  (6)

#define I_CO2   (0)
#define I_psi   (1)
#define I_Ca_S  (2)
#define I_Si_S   (3)
#define I_K     (4)
#define I_Cl    (5)
#define I_P_l   (6)

#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_CO2   LOG_RHO
#define U_PL    RHO

#if (U_CO2 == LOG_RHO)
  #define P_CO2(n)   (exp(Ln10*u[(n)][I_CO2]))
  #define P_CO2n(n)  (exp(Ln10*u_n[(n)][I_CO2]))
#else
  #define P_CO2(n)   (u[(n)][I_CO2])
  #define P_CO2n(n)  (u_n[(n)][I_CO2])
#endif
#if (U_PL == LOG_RHO)
  #define P_l(n)   (exp(Ln10*u[(n)][I_P_l]))
  #define P_ln(n)  (exp(Ln10*u_n[(n)][I_P_l]))
#else
  #define P_l(n)   (u[(n)][I_P_l])
  #define P_ln(n)  (u_n[(n)][I_P_l])
#endif
#define ZN_Ca_S(n)  (u[(n)][I_Ca_S])
#define ZN_Si_S(n)  (u[(n)][I_Si_S])
#define PSI(n)     (u[(n)][I_psi])
#define X_K(n)     (u[(n)][I_K])
#define X_Cl(n)    (u[(n)][I_Cl])

#define ZN_Ca_Sn(n) (u_n[(n)][I_Ca_S])
#define ZN_Si_Sn(n) (u_n[(n)][I_Si_S])
#define PSIn(n)    (u_n[(n)][I_psi])
#define X_Kn(n)    (u_n[(n)][I_K])
#define X_Cln(n)   (u_n[(n)][I_Cl])

#define N_C(n)     (f[(n)])
#define N_q(n)     (f[(2+n)])
#define N_Ca(n)    (f[(4+n)])
#define N_Si(n)    (f[(6+n)])
#define N_K(n)     (f[(8+n)])
#define N_Cl(n)    (f[(10+n)])
#define M(n)       (f[(12+n)])
#define W_C        (f[14])
#define W_q        (f[15])
#define W_Ca       (f[16])
#define W_Si       (f[17])
#define W_K        (f[18])
#define W_Cl       (f[19])
#define W_M        (f[20])
#define N_CH(n)    (f[(21+n)])
#define N_CC(n)    (f[(23+n)])
#define N_Si_S(n)  (f[(25+n)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define M_n(n)     (f_n[(12+n)])
#define N_CHn(n)   (f_n[(21+n)])
#define N_CCn(n)   (f_n[(23+n)])
#define N_Si_Sn(n) (f_n[(25+n)])


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

#define KD_OH       (va[(30)])
#define KD_H        (va[(31)])
#define KD_CO2      (va[(32)])
#define KD_H2CO3    (va[(33)])
#define KD_HCO3     (va[(34)])
#define KD_CO3      (va[(35)])
#define KD_Ca       (va[(36)])
#define KD_CaHCO3   (va[(37)])
#define KD_CaH3SiO4 (va[(38)])
#define KD_H3SiO4   (va[(39)])
#define KD_H4SiO4   (va[(40)])
#define KD_H2SiO4   (va[(41)])
#define KD_CaH2SiO4 (va[(42)])
#define KD_CaCO3aq  (va[(43)])
#define KD_CaOH     (va[(44)])
#define KD_K        (va[(45)])
#define KD_Cl       (va[(46)])
#define KD_M        (va[(47)])

#define KD_GCO2       (va[(48)])    /*gas transport*/
#define KD_GH2O       (va[(49)])

#define N_CH0(n)      (va[(50+n)])
#define N_CC0(n)      (va[(52+n)])
#define ZP_CH0(n)     (va[(54+n)])
#define N_Si_S0(n)    (va[(56+n)])


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
#define v_co2      (32.81e-3)    /* d'apres Lothenbach */
#define v_h2co3    (50.e-3)
#define v_hco3     (24.21e-3)    /* d'apres Lothenbach */
#define v_co3      (-6.06e-3)    /* d'apres Lothenbach */
#define v_ca       (-18.44e-3)   /* d'apres Lothenbach */
#define v_sioh4    (xxx)
#define v_h3sio4   (4.53e-3)     /* d'apres Lothenbach */
#define v_h2sio4   (34.13e-3)    /* d'apres Lothenbach */
#define v_cah2sio4 (15.69e-3)    /* d'apres Lothenbach */
#define v_cah3sio4 (-6.74e-3)    /* d'apres Lothenbach */
#define v_caco3aq  (-15.65e-3)   /* d'apres Lothenbach */
#define v_caoh     (5.76e-3)     /* d'apres Lothenbach */
#define v_k        (9.06e-3)     /* d'apres Inorgqnic species in geologic fluids*/
#define v_koh      (3.511e-3)    /* a modifier*/
#define v_cl       (17.34e-3)    /* d'apres kulik 2001*/
#define v_h4sio4   (16.06e-3)    /* d'apres Lothenbach */
#define v_cahco3   (13.33e-3)    /* d'apres Lothenbach */
#define v_na       (-1.21e-3)    /* d'apres Lothenbach */
#define v_naoh     (3.511e-3)    /* d'apres Lothenbach */
#define v_naco3    (-0.42e-3)    /* d'apres Lothenbach */
#define v_nahco3   (32.32e-3)    /* d'apres Lothenbach */


/* coefficients de diffusion moleculaire (dm2/s) */
/* Stokes-Einstein ~ kT/(6*pi*mu*r)  kT/(6*pi*mu) = 2.1451e-19 m3/s T=298K  */
#define T			temperature
#define T0			(298.0)
#define d_oh       (5.273e-7*T/T0)     /* chemistry handbook */
#define d_h        (9.311e-7*T/T0)    /* chemistry handbook (9.311e-7)*/
#define d_co2      (1.91e-7*T/T0)    /* chemistry handbook (1.91e-7)*/
#define d_h2co3    (7.2e-8*T/T0)
#define d_hco3     (1.185e-7*T/T0)    /* chemistry handbook */
#define d_co3      (9.55e-8*T/T0)    /* 9.53e-8 (radius = 2.25e-10 m) */
#define d_ca       (7.92e-8*T/T0)    /*7.92e-8*/
#define d_cahco3   (1.07e-7*T/T0)    /* (radius = 2e-10 m) 1.07e-7*/
#define d_sioh4    (xxx)
#define d_h4sio4   (1.07e-7*T/T0)    /* 1.07e-7 a modifier*/
#define d_h3sio4   (1.07e-7*T/T0)    /*(radius = 2e-10 m) a modifier*/
#define d_h2sio4   (1.07e-7*T/T0)   /* (radius = 2e-10 m) a modifier*/
#define d_cah2sio4 (1.07e-7*T/T0)    /* a modifier */
#define d_cah3sio4 (1.07e-7*T/T0)    /* (radius = 2e-10 m) a modifier*/
#define d_caco3aq  (1.43e-7*T/T0)    /* (radius = 1.5e-10 m) */
#define d_caoh     (1.07e-7*T/T0)    /* (radius = 2e-10 m) 1.07e-7*/
#define d_k        (1.957e-7*T/T0)   /* chemistry handbook (1.957e-7)*/
#define d_koh      (1.96e-7*T/T0)   
#define d_cl       (2.032e-7*T/T0)   /* chemistry handbook (2.032e-7)*/
#define d_na       (1.334e-7*T/T0)   /* chemistry handbook */

/* viscosite (Pa.s) */
#define T_mu		((1.3272*(293.-T)-0.001053*(T-293.)*(T-293.))/(T-168.))
#define mu_l       	(1.002e-3*pow(10.,T_mu))      /*0.547e-3 unit pa.s t=50°c,www.engineeringtoolbox.com, equation from chemistry handbook*/

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
#define V_CH       (33.e-3)      /* (33.e-3) */
#define V_CC       (37.e-3)      /* (37.e-3) */

/* constantes d'equilibre (ref = 1 mole/L) */
#define K_CH       (6.456e-6)        /* logk = -5.2 CH = Ca[2+] + 2OH[-] */
#define K_CC       (3.89e-9)         /* logk = -8.36 CC = Ca[2+] + CO3[2-]  */ 

/* C-S-H */
/* constantes d'equilibre (ref = 1 mole/L) */
/* CxSyHz = xCa[2+] + 2xOH[-] + yH4SiO4 + (z-x-2y)H2O */
/* Stoichiometric Coefficients */
#define X_CSH(q)    courbe(q,el.mat->cb[3])
#define DX_CSH(q)   dcourbe(q,el.mat->cb[3])
/* Water content */
#define Z_CSH(q)    courbe(q,el.mat->cb[4])
#define DZ_CSH(q)   dcourbe(q,el.mat->cb[4])
/* Molar Volume */
#define V_CSH(q)    courbe(q,el.mat->cb[5])
#define DV_CSH(q)   dcourbe(q,el.mat->cb[5])
/* Amorphous Silica (ref = 1 mole/L) */
#define K_Sam      (1.93642e-3)      /* S + 2H2O = H4SiO4 */
/* Amorphous Silica Solid fraction */
#define F_SAM(q)    courbe(q,el.mat->cb[6])
#define DF_SAM(q)   dcourbe(q,el.mat->cb[6])
/* Other Solid Fractions */
#define N_SOLIDCSH  (el.mat->nc - 6)
#define F_CSH(q,i)  courbe(q,el.mat->cb[6+i])


/* Masses molaires (unite arbitraire = M_H) */
#define M_Ca    		 (40.)
#define M_H2CO3  		 (62.)
#define M_HCO3   		 (61.)
#define M_CO3    		 (60.)
#define M_OH     		 (17.)
#define M_H      		 (1.)
#define M_H2O    		 (18.)
#define M_Na	 		 (23.)
#define M_NaOH	  	 	 (40.)
#define M_NaHCO3  		 (84.)
#define M_NaCO3	  		 (83.)
#define M_CO2     	 	 (44.)
#define M_CaOH2   	 	 (74.)
#define M_CaCO3   		 (100.)
#define M_CaO 			 (56.)
#define M_SiO2 			 (60.)
#define M_Jen       	 (X_Jen*M_CaO+ Y_Jen*M_SiO2 + Z_Jen*M_H2O)  
#define M_TobII    		 (X_TobII*M_CaO + Y_TobII*M_SiO2 + Z_TobII*M_H2O)  
#define M_TobI     		 (X_TobI*M_CaO + Y_TobI*M_SiO2 + Z_TobI*M_H2O)  
#define M_Sam     		 (X_Sam*M_CaO + Y_Sam*M_SiO2 + Z_Sam*M_H2O)   
#define M_K      	 	 (39.)  
#define M_KOH    		 (56.)
#define M_CaOH    		 (57.)
#define M_CaHCO3  		 (101.)
#define M_CaCO3aq 		 (100.)
#define M_CaOH2aq  		 (74.)
#define M_CaH2SiO4 		 (134.)
#define M_CaH3SiO4  	 (135.)
#define M_H3SiO4  		 (95.)
#define M_H2SiO4  		 (94.)
#define M_H4SiO4  		 (96.)
#define M_Cl	  		 (35.5)
#define M_H2Og	  		 (18.)

#define C_CO2_eq   (K_h2o*K_h2o*K_CC/(K_h2co3*K_hco3*K_co3*K_CH))
#define P_CO2_eq   (C_CO2_eq/dproductco2_co2(1.e4,1.e7,temperature))/*(K_h2o*K_h2o*K_CC/(K_h2co3*K_hco3*K_co3*K_CH)*1.e8)*/

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (mC/mole = Pa.dm3/V/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (Pa.dm3/mole) */
#define FsRT      (3.961e1)   /* F/RT (1/V) */


/* Fonctions */
static int    pm(char *s) ;
static double concentration_oh(double,double,double,double,double,elem_t) ;
static double poly4(double,double,double,double,double) ;
static double poly4_sansc(double,double,double,double) ;
static void   transfert(double**,double**,double*,double*,elem_t,int,geom_t) ;
static void   flux(double**,double**,double*,double*,elem_t,int,geom_t) ;
static double product_Sam(double,double,elem_t) ;
static double productx_co2(double,double,double) ;
static double product_co2(double,double,double) ;
static double producty_h2o(double,double,double) ;
static double product_gco2(double,double) ;
static double dproduct_gco2(double,double) ;
static double product_gh2o(double,double,double) ;
static double dproductgh2o_h2o(double,double,double) ;
static double dproductgh2o_co2(double,double,double) ;
static double dproductco2_co2(double,double,double) ;
static double dproductco2_h2o(double,double,double) ;
static double quartic_equation(double,double,double,double,double,double *,double *);
static double cubic_equation(double,double,double,double);
static double product_muco2(double,double) ;
static void   quadratic_equation(double,double,double,double *);
#define MIN(a,b)  ((a < b) ? a : b)
#define MAX(a,b)  ((a > b) ? a : b)
#define NEGEXP(y) ((y < 0.) ? exp(y) : 1.)

/* Parametres */
static double phi0,c_co2_eq,t_ch,t_cc,k_intl,k_intg,temperature ;
static double n_ch_ref,n_cc_ref ;
static double p_l0 = 0.;/*-4.e7 ;*/


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
  else if(strcmp(s,"k_intl") == 0)    return (9) ;
  else if(strcmp(s,"k_intg") == 0)    return (10) ;
  else if(strcmp(s,"temperature") == 0)    return (11) ;
  else return(-1) ;
}


int dm590(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 12 ;            /* ? */
  
  mat->neq      = NEQ ;

  strcpy(mat->eqn[E_C],     "carbone") ;
  strcpy(mat->eqn[E_Ca],    "calcium") ;
  strcpy(mat->eqn[E_Si],   "Silicium") ;
  strcpy(mat->eqn[E_q],      "charge") ;
  strcpy(mat->eqn[E_K],   "potassium") ;
  strcpy(mat->eqn[E_Cl],   "chloride") ;
  strcpy(mat->eqn[E_mass],    "masse") ;  

#if (U_CO2 == LOG_RHO)
  strcpy(mat->inc[I_CO2],  "logp_co2") ;
#else
  strcpy(mat->inc[I_CO2],  "p_co2") ;
#endif
#if (U_PL == LOG_RHO)
  strcpy(mat->inc[I_P_l],  "logp_l") ;
#else
  strcpy(mat->inc[I_P_l],  "p_l") ; 
#endif  
  strcpy(mat->inc[I_Ca_S], "z_ca") ;
  strcpy(mat->inc[I_psi],  "psi") ;
  strcpy(mat->inc[I_Si_S],  "z_si") ;
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


int qm590(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n") ;
  printf("Le systeme est forme de 6 equations:\n") ;
#if (U_CO2 == LOG_RHO)
  printf("\t- la conservation de la masse de C      (logp_co2)\n") ;
#else
  printf("\t- la conservation de la masse de C      (p_co2)\n") ;
#endif
#if (U_PL == LOG_RHO)
  printf("\t- la conservation de la masse totale      (logp_l)\n") ;
#else
  printf("\t- la conservation de la masse totale      (p_l)\n") ;
#endif
  printf("\t- la conservation de la charge          (psi)\n") ;
  printf("\t- la conservation de la masse de Ca     (zn_ca_s)\n") ;
  printf("\t- la conservation de la masse de Si     (zn_si_s)\n") ;
  printf("\t- la conservation de la masse de K      (c_k)\n") ;
  printf("\t- la conservation de la masse de Cl     (c_cl)\n") ;
  printf("\t- la conservation de la masse totale    (p_l)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n") ;

  printf("Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1    # Contenu en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"zn_si_s = 2.4    # contenu en CSH (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4    # contenu en K (moles/L)\n") ;  
  fprintf(ficd,"N_Cl  = 0.4    # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5   # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"courbes = solid   # Nom du fichier: log10(P_CH/K_CH) Xi\n") ;

  return(NEQ) ;
}



void tb590(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = (el->dim < dim) ? 0 : NVE ;
}


void ch590(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in590(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
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
  k_intl     = el.mat->pr[pm("k_intl")] ;
  k_intg     = el.mat->pr[pm("k_intg")] ;
  temperature 	= el.mat->pr[pm("temperature")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double p_co2	  = (P_CO2(i) > 0.) ? P_CO2(i) : P_CO2_eq ;
    double p_l     	  = P_l(i) ;
    double gc_co2     = product_gco2(p_co2,temperature) ;/*(co2 concentration in gas)*/
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ;
    double c_co2      = product_co2(p_co2,p_l,temperature) ;/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s    = ZN_Si_S(i) ;
    double zn_ca_s    = ZN_Ca_S(i) ;
    double c_k        = X_K(i);                 
    double c_cl       = X_Cl(i); 
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0 ;
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ;
    
    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;
    double zp_ch      = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ;
    double x_s        = product_Sam(c_co2,zn_ca_s,el);

    double c_h4sio4   = NEGEXP(zn_si_s)*x_s * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el) ;
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
    double c_h2o      = (1. - (c_co2*v_co2 + c_h2co3*v_h2co3 + c_hco3*v_hco3 + c_co3*v_co3 + c_h*v_h + c_oh*v_oh + c_ca*v_ca \
						+ c_caoh*v_caoh + v_cahco3*c_cahco3 + c_caco3aq*v_caco3aq \
						+ c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 + c_cah2sio4*v_cah2sio4  \
						+ c_cah3sio4*v_cah3sio4+ c_k*v_k + c_cl*v_cl))/v_h2o ;
						/* + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + \
						 * x_koh*v_koh + c_caoh2aq*v_caoh2aq \*/

    /* masse volumique liquide */
    double rho_l   = M_CO2*c_co2 + M_H*c_h + M_OH*c_oh + M_H2O*c_h2o + M_H2CO3*c_h2co3 + M_HCO3*c_hco3 + M_CO3*c_co3 + M_Ca*c_ca \
			+ M_CaOH*c_caoh + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 \
			+ M_H2SiO4*c_h2sio4 + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 + M_Cl*c_cl + M_K*c_k ;
			/* + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_KOH*x_koh + M_CaOH2aq*x_caoh2aq */
    
    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(zn_ca_s,0.) ;
    double n_cc_eq    = n_cc_ref*MAX(zn_ca_s,0.) ;
    
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch       = (z_co2 <= 1) ? n_ch_eq  : 0 ;
    double n_cc       = (z_co2 >  1) ? n_cc_eq  : 0 ;
    double x_csh      = X_CSH(zp_ch) ;
    double z_csh      = Z_CSH(zp_ch);
      
    /* solides */
    double n_si_s     = MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;

    /* porosite */
    double phi = phi0 ;

    /* contenus molaires */
    double n_co2      = phi*s_l*c_co2 + phi*s_g*gc_co2 ;
    double n_hco3     = s_l*phi*c_hco3 ;
    double n_h2co3    = s_l*phi*c_h2co3 ;
    double n_co3      = s_l*phi*c_co3 ;
    double n_ca       = s_l*phi*c_ca ;
    double n_cahco3   = s_l*phi*c_cahco3 ;
    double n_h3sio4   = s_l*phi*c_h3sio4 ;
    double n_h4sio4   = s_l*phi*c_h4sio4 ;
    double n_cah3sio4 = s_l*phi*c_cah3sio4 ;
    double n_h2sio4   = s_l*phi*c_h2sio4 ;          /* ajouter */
    double n_cah2sio4 = s_l*phi*c_cah2sio4 ;        /* ajouter */
    double n_caco3aq  = s_l*phi*c_caco3aq ;         /* ajouter */
    double n_caoh     = s_l*phi*c_caoh ;            /* ajouter */
    double n_k        = s_l*phi*c_k;                /* ajouter2 */
    double n_cl       = s_l*phi*c_cl;     
    
    
    N_C(i)  = n_co2 + n_h2co3 + n_hco3 + n_co3 + n_cahco3 + n_cc + n_caco3aq ;
    N_Ca(i) = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh + n_ca_s ;
    N_Si(i) = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4+ n_cah2sio4  + n_si_s ;
    N_K(i)  = n_k;                              
    N_Cl(i) = n_cl;
    
    /* masse totale */
    /*double m_csh = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;*/
    M(i)    = M_H2Og*phi*s_g*gc_h2o + M_CO2*phi*s_g*gc_co2 + s_l*phi*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s;

    /* densite de charge */
    N_q(i)  = z_h*c_h + z_oh*c_oh + z_hco3*c_hco3 + z_co3*c_co3 + z_ca*c_ca + z_cahco3*c_cahco3 + z_h3sio4*c_h3sio4 \
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k + z_cl*c_cl;

#if (U_CO2 == LOG_RHO)
    u[i][I_CO2] = log(p_co2)/Ln10 ;
#else
    P_CO2(i)   = p_co2 ;
#endif
#if (U_PL == LOG_RHO)
    u[i][I_P_l] = log(p_l)/Ln10 ;
#else
    P_l(i)   = p_l ;
#endif
    ZN_Si_S(i)   = zn_si_s ;
    ZN_Ca_S(i)  = zn_ca_s ;
    X_K(i)     = c_k;
    X_Cl(i)    = c_cl;
    
    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;
    N_Si_S(i)  = n_si_s ;

    N_CH0(i)   = n_ch ;
    N_CC0(i)   = n_cc ;
    ZP_CH0(i)  = zp_ch ;
    N_Si_S0(i) = n_si_s ;
  }
  
  if(el.dim < dim) return ;

  /* Coefficient de transfert */
  transfert(x,u,f,va,el,dim,geom) ;

  /* Flux */
  flux(x,u,f,va,el,dim,geom) ;
}


int ex590(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Thermes explicites (va)  */
{
  
  if(el.dim < dim) return(0) ;
  
  /*
    Coefficients de transfert
  */
  transfert(x,u,f,va,el,dim,geom) ;

  return(0) ;
}


int ct590(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
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
  temperature 	= el.mat->pr[pm("temperature")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    /* molarites */
    double p_co2	  = (P_CO2(i) > 0.) ? P_CO2(i) : P_CO2_eq ;
    double p_l     	  = P_l(i) ;    
    double gc_co2     = product_gco2(p_co2,temperature) ;/*(co2 concentration in gas)*/
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ;
    double c_co2      = product_co2(p_co2,p_l,temperature) ;/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s      = ZN_Si_S(i) ;
    double zn_ca_s     = ZN_Ca_S(i) ;
    double c_k        = X_K(i);                 
    double c_cl       = X_Cl(i); 
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0 ; ;
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ;

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;
    double zp_ch      = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ;
    double x_s        = product_Sam(c_co2,zn_ca_s,el);

    double c_h4sio4   = NEGEXP(zn_si_s)*x_s * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el)  ;
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

    double c_h2o      = (1. - (c_co2*v_co2 + c_h2co3*v_h2co3 + c_hco3*v_hco3 + c_co3*v_co3 + c_h*v_h + c_oh*v_oh + c_ca*v_ca \
						+ c_caoh*v_caoh + v_cahco3*c_cahco3 + c_caco3aq*v_caco3aq \
						+ c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 + c_cah2sio4*v_cah2sio4  \
						+ c_cah3sio4*v_cah3sio4+ c_k*v_k + c_cl*v_cl))/v_h2o ;
						/* + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + \
						 * x_koh*v_koh + c_caoh2aq*v_caoh2aq \*/

    /* masse volumique liquide */
    double rho_l   = M_CO2*c_co2 + M_H*c_h + M_OH*c_oh + M_H2O*c_h2o + M_H2CO3*c_h2co3 + M_HCO3*c_hco3 + M_CO3*c_co3 + M_Ca*c_ca \
			+ M_CaOH*c_caoh + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 \
			+ M_H2SiO4*c_h2sio4 + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 + M_Cl*c_cl + M_K*c_k ;
			/* + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_KOH*x_koh + M_CaOH2aq*x_caoh2aq */
    
    /* cinetiques */
    double n_chn      = N_CHn(i) ;
    double n_ccn      = N_CCn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;  /* si z_co2 > 1, s_l*n_chn*pow(z_co2,-dt/t_ch)*/
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;   /* si z_co2 < 1 */

    /* equilibres */
    double n_ch_eq    = n_ch_ref*MAX(zn_ca_s,0.) - n_cc_ci ;
    double n_cc_eq    = n_cc_ref*MAX(zn_ca_s,0.) - n_ch_ci ;

    
    /* contenus en CH, CC, CSH, S(am) */
    double n_ch   = (z_co2 <= 1) ? n_ch_eq  : n_ch_ci ;
    double n_cc   = (z_co2 >  1) ? n_cc_eq  : n_cc_ci ;
    double x_csh      = X_CSH(zp_ch) ;
    double v_csh      = V_CSH(zp_ch) ;
    double z_csh      = Z_CSH(zp_ch) ;
      
    /* solides */
    double n_si_s     = MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch_ref*MAX(zn_ca_s,0.) + x_csh*n_si_s ;

    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double zp_ch0     = ZP_CH0(i) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double n_si_s0	  = N_Si_S0(i) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) + v_csh0*n_si_s0- v_csh*n_si_s;   
    /* masse totale */
    /*double m_csh = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;*/
    M(i)    = M_H2Og*phi*s_g*gc_h2o + M_CO2*phi*s_g*gc_co2 + s_l*phi*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s;
    
    /* contenus molaires */
    double n_co2      = phi*s_l*c_co2 + phi*s_g*gc_co2 ;
    double n_hco3     = s_l*phi*c_hco3 ;
    double n_h2co3    = s_l*phi*c_h2co3 ;
    double n_co3      = s_l*phi*c_co3 ;
    double n_ca       = s_l*phi*c_ca ;
    double n_cahco3   = s_l*phi*c_cahco3 ;
    double n_h3sio4   = s_l*phi*c_h3sio4 ;
    double n_h4sio4   = s_l*phi*c_h4sio4 ;
    double n_cah3sio4 = s_l*phi*c_cah3sio4 ;
    double n_h2sio4   = s_l*phi*c_h2sio4 ;          /* ajouter */
    double n_cah2sio4 = s_l*phi*c_cah2sio4 ;        /* ajouter */
    double n_caco3aq  = s_l*phi*c_caco3aq ;         /* ajouter */
    double n_caoh     = s_l*phi*c_caoh ;            /* ajouter */
    double n_k        = s_l*phi*c_k;                /* ajouter2 */
    double n_cl       = s_l*phi*c_cl;    

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
    N_Si_S(i)  = n_si_s ;

    if(c_co2 < 0. || n_ca_s < 0. || n_si_s < 0.) {
      printf("x         = %e\n",x[i][0]) ;
      printf("c_co2     = %e\n",c_co2) ;
      printf("n_cc      = %e\n",n_cc) ;
      printf("n_ca_s    = %e\n",n_ca_s) ; 
      printf("n_si_s    = %e\n",n_si_s) ;
      printf("zn_si_s     = %e\n",zn_si_s) ;
      printf("zn_ca_s    = %e\n",zn_ca_s) ;
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

int mx590(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  double dx,xm ;
  double volume[2],surf ;
  int    i,j ;
  double c[2] ;

  
  double Dc_hSDp_co2[2]        ;
  double Dc_ohSDp_co2[2]       ;
  double Dc_h2co3SDp_co2[2]    ;
  double Dc_hco3SDp_co2[2]     ;
  double Dc_co3SDp_co2[2]      ;
  double Dc_caSDp_co2[2]       ;
  double Dc_cahco3SDp_co2[2]   ;
  double Dc_h3sio4SDp_co2[2]   ;
  double Dc_h4sio4SDp_co2[2]   ;
  double Dc_cah3sio4SDp_co2[2] ;
  double Dc_h2sio4SDp_co2[2]   ;
  double Dc_cah2sio4SDp_co2[2] ;
  double Dc_caco3aqSDp_co2[2]  ;  
  double Dc_caohSDp_co2[2]     ;
  double Dc_co2SDp_co2[2]     ;
  
  double Dc_hSDp_l[2]        ;
  double Dc_ohSDp_l[2]       ;
  double Dc_h2co3SDp_l[2]    ;
  double Dc_hco3SDp_l[2]     ;
  double Dc_co3SDp_l[2]      ;
  double Dc_caSDp_l[2]       ;
  double Dc_cahco3SDp_l[2]   ;
  double Dc_h3sio4SDp_l[2]   ;
  double Dc_h4sio4SDp_l[2]   ;
  double Dc_cah3sio4SDp_l[2] ;
  double Dc_h2sio4SDp_l[2]   ;
  double Dc_cah2sio4SDp_l[2] ;
  double Dc_caco3aqSDp_l[2]  ;  
  double Dc_caohSDp_l[2]     ;
  double Dc_co2SDp_l[2]     ;
  
  double Dc_hSDzn_si_s[2]        ;
  double Dc_ohSDzn_si_s[2]       ;
  double Dc_hco3SDzn_si_s[2]     ;
  double Dc_co3SDzn_si_s[2]      ;
  double Dc_caSDzn_si_s[2]       ;
  double Dc_cahco3SDzn_si_s[2]   ;
  double Dc_h3sio4SDzn_si_s[2]   ;
  double Dc_h4sio4SDzn_si_s[2]   ;
  double Dc_cah3sio4SDzn_si_s[2] ;
  double Dc_h2sio4SDzn_si_s[2]   ;
  double Dc_cah2sio4SDzn_si_s[2] ;
  double Dc_caco3aqSDzn_si_s[2]  ;
  double Dc_caohSDzn_si_s[2]     ;

  double Dc_hSDzn_ca_s[2]        ;
  double Dc_ohSDzn_ca_s[2]       ;
  double Dc_hco3SDzn_ca_s[2]     ;
  double Dc_co3SDzn_ca_s[2]      ;
  double Dc_caSDzn_ca_s[2]       ;
  double Dc_cahco3SDzn_ca_s[2]   ;
  double Dc_h3sio4SDzn_ca_s[2]   ;
  double Dc_h4sio4SDzn_ca_s[2]   ;
  double Dc_cah3sio4SDzn_ca_s[2] ;
  double Dc_h2sio4SDzn_ca_s[2]   ;
  double Dc_cah2sio4SDzn_ca_s[2] ;
  double Dc_caco3aqSDzn_ca_s[2]  ;  
  double Dc_caohSDzn_ca_s[2]     ;
  
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
  temperature 	= el.mat->pr[pm("temperature")] ;

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
    double p_co2	  = (P_CO2(i) > 0.) ? P_CO2(i) : P_CO2_eq ;
    double p_l     	  = P_l(i) ;
    double gc_co2     = product_gco2(p_co2,temperature);/*(co2 concentration in gas)*/
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ;    
    double c_co2      = product_co2(p_co2,p_l,temperature) ;/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s      = ZN_Si_S(i) ;
    double zn_ca_s     = ZN_Ca_S(i) ;
    double c_k        = X_K(i);                 
    double c_cl       = X_Cl(i);    
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0 ; 
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ; 

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;
    double zp_ch      = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ;
    double x_s        = product_Sam(c_co2,zn_ca_s,el);

    double c_h4sio4   = NEGEXP(zn_si_s)*x_s * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el) ;
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
    double c_h2o      = (1. - (c_co2*v_co2 + c_h2co3*v_h2co3 + c_hco3*v_hco3 + c_co3*v_co3 + c_h*v_h + c_oh*v_oh + c_ca*v_ca \
						+ c_caoh*v_caoh + v_cahco3*c_cahco3 + c_caco3aq*v_caco3aq \
						+ c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 + c_cah2sio4*v_cah2sio4  \
						+ c_cah3sio4*v_cah3sio4+ c_k*v_k + c_cl*v_cl))/v_h2o ;
						/* + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + \
						 * x_koh*v_koh + c_caoh2aq*v_caoh2aq \*/
	
	double x_csh 			  = X_CSH(zp_ch) ;
    double v_csh      		  = V_CSH(zp_ch) ;
    double z_csh 			  = Z_CSH(zp_ch);

    /* masse volumique liquide */
    double rho_l   = M_CO2*c_co2 + M_H*c_h + M_OH*c_oh + M_H2O*c_h2o + M_H2CO3*c_h2co3 + M_HCO3*c_hco3 + M_CO3*c_co3 + M_Ca*c_ca \
			+ M_CaOH*c_caoh + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 \
			+ M_H2SiO4*c_h2sio4 + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 + M_Cl*c_cl + M_K*c_k ;
			/* + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_KOH*x_koh + M_CaOH2aq*x_caoh2aq */    

    /* contenus en CH, CC, CSH, Sam */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;

    /* cinetiques */
    double n_ccn      = N_CCn(i) ;
    double n_chn      = N_CHn(i) ;
    double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;/*n_chn*(1.-s_l*(1.-pow(z_co2,-dt/t_ch))) */
    double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;
    
    /* solides */
    double n_si_s     = MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch_ref*MAX(zn_ca_s,0.) + x_csh*n_si_s ;


    /* porosite */
    double n_ch0      = N_CH0(i) ;
    double n_cc0      = N_CC0(i) ;
    double zp_ch0     = ZP_CH0(i) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double n_si_s0	  = N_Si_S0(i) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) + v_csh0*n_si_s0- v_csh*n_si_s;   
    /* masse totale */

    M(i)    = M_H2Og*phi*s_g*gc_h2o + M_CO2*phi*s_g*gc_co2 + s_l*phi*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s;

    /* derivees ... */
    /* ... par rapport a p_co2 */
    double dc_co2sdp_co2      = dproductco2_co2(p_co2,p_l,temperature);
    double dgc_h2osdp_co2   = dproductgh2o_co2(p_co2,p_l,temperature);
    double dc_co2             = 1.e-4*c_co2; /*((c_co2 > C_CO2n(i)) ? 1 : -1) */
    double c_co22             = c_co2 + dc_co2 ;
    double z_co22			  = c_co22/c_co2_eq ;
    double c_oh2              = concentration_oh(c_co22,zn_si_s,zn_ca_s,c_k,c_cl,el) ;   
    double dc_ohsdp_co2       = dc_co2sdp_co2*(c_oh2 - c_oh)/dc_co2 ;
    double dc_hsdp_co2        = -c_h*dc_ohsdp_co2/c_oh ;
    double dc_h2co3sdp_co2    = dc_co2sdp_co2*K_h2co3 ;
    double dc_hco3sdp_co2     = c_hco3*(dc_h2co3sdp_co2/c_h2co3 + dc_ohsdp_co2/c_oh) ;
    double dc_co3sdp_co2      = c_co3*(dc_hco3sdp_co2/c_hco3 + dc_ohsdp_co2/c_oh) ;		                    
			                     
    double dP_CCsdp_co2       = (z_co2 < 1.) ? dc_co2sdp_co2*P_CC/c_co2 : 0. ;
    
    double dc_casdp_co2       = (dP_CCsdp_co2 - c_ca*dc_co3sdp_co2)/c_co3 ;
    double dc_cahco3sdp_co2   = K_cahco3*(dc_casdp_co2*c_hco3 + c_ca*dc_hco3sdp_co2) ;
			                    
    double zp_ch2			  = NEGEXP(zn_ca_s)/MAX(z_co22,1.) ;
    double x_s2				  = product_Sam(c_co22,zn_ca_s,el) ;
    double x_csh2 			  = X_CSH(zp_ch2);
    double v_csh2 			  = V_CSH(zp_ch2) ;
    double z_csh2 			  = Z_CSH(zp_ch2) ;
    
    double dc_h4sio4sdp_co2   = NEGEXP(zn_si_s)*dc_co2sdp_co2*(x_s2 - x_s )* K_Sam/dc_co2 ;
    double dc_h3sio4sdp_co2   = c_h3sio4*(dc_h4sio4sdp_co2/c_h4sio4 - dc_hsdp_co2/c_h) ;
    double dc_cah3sio4sdp_co2 = K_cah3sio4*(dc_casdp_co2*c_h3sio4 + c_ca*dc_h3sio4sdp_co2) ;
    double dc_h2sio4sdp_co2   = K_h2sio4*(dc_h3sio4sdp_co2*c_oh + dc_ohsdp_co2*c_h3sio4) ;
    double dc_cah2sio4sdp_co2 = K_cah2sio4*(dc_h2sio4sdp_co2*c_ca	+ dc_casdp_co2*c_h2sio4) ;
    double dc_caco3aqsdp_co2  = K_caco3aq*(dc_co3sdp_co2*c_ca + dc_casdp_co2*c_co3) ;				                  
    double dc_caohsdp_co2     = K_caoh*(dc_casdp_co2*c_oh + dc_ohsdp_co2*c_ca) ;			                  	
		                    
    double dn_ch_cisdp_co2   = dc_co2sdp_co2*n_ch_ci*(-dt/t_ch)/c_co2 ;
    double dn_cc_cisdp_co2   = dc_co2sdp_co2*n_cc_ci*(dt/t_cc)/c_co2 ;


    double dn_ch_eqsdp_co2   = -dn_cc_cisdp_co2 ;
    double dn_cc_eqsdp_co2   = -dn_ch_cisdp_co2 ;


    double dn_chsdp_co2   = (z_co2 <= 1) ? dn_ch_eqsdp_co2  : dn_ch_cisdp_co2 ;
    double dn_ccsdp_co2   = (z_co2 >  1) ? dn_cc_eqsdp_co2  : dn_cc_cisdp_co2 ;
    double dx_cshsdp_co2  = dc_co2sdp_co2*(x_csh2- x_csh) / dc_co2 ;
    double dv_cshsdp_co2  = dc_co2sdp_co2*(v_csh2- v_csh) / dc_co2 ;
    double dz_cshsdp_co2  = dc_co2sdp_co2*(z_csh2- z_csh) / dc_co2 ;

    double dn_ca_ssdp_co2 = dx_cshsdp_co2 * MAX(zn_si_s,0.) ;
    double dn_si_ssdp_co2 = 0. ;

    double dphisdp_co2        = - V_CH*dn_chsdp_co2 - dv_cshsdp_co2*n_si_s - V_CC*dn_ccsdp_co2 ;
    double dc_h2osdp_co2      = -(dc_co2sdp_co2*v_co2 + v_h2co3*dc_h2co3sdp_co2 + v_hco3*dc_hco3sdp_co2 + v_co3*dc_co3sdp_co2 + v_h*dc_hsdp_co2 \
								+ v_oh*dc_ohsdp_co2 + v_ca*dc_casdp_co2 + v_caoh*dc_caohsdp_co2 + v_cahco3*dc_cahco3sdp_co2 + v_caco3aq*dc_caco3aqsdp_co2 \
								+ v_h3sio4*dc_h3sio4sdp_co2 + v_h4sio4*dc_h4sio4sdp_co2 + v_h2sio4*dc_h2sio4sdp_co2 + v_cah2sio4*dc_cah2sio4sdp_co2 \
								+ v_cah3sio4*dc_cah3sio4sdp_co2)/v_h2o ;
	double drho_lsdp_co2	  = M_H2O*dc_h2osdp_co2 + dc_co2sdp_co2*M_CO2 + M_H2CO3*dc_h2co3sdp_co2 + M_HCO3*dc_hco3sdp_co2 + M_CO3*dc_co3sdp_co2 + M_H*dc_hsdp_co2 \
								+ M_OH*dc_ohsdp_co2 + M_Ca*dc_casdp_co2 + M_CaOH*dc_caohsdp_co2 + M_CaHCO3*dc_cahco3sdp_co2 + M_CaCO3*dc_caco3aqsdp_co2 \
								+ M_H3SiO4*dc_h3sio4sdp_co2 + M_H4SiO4*dc_h4sio4sdp_co2 + M_H2SiO4*dc_h2sio4sdp_co2 + M_CaH2SiO4*dc_cah2sio4sdp_co2 \
								+ M_CaH3SiO4*dc_cah3sio4sdp_co2 ;   
	double dgc_co2sdp_co2     = dproduct_gco2(p_co2,temperature) ;
	double ds_lsdp_co2        = dcourbe(p_c,el.mat->cb[0]) ;
	double ds_gsdp_co2        = 0.- dcourbe(p_c,el.mat->cb[0]) ;	
								
    /*double dn_co2sdp_co2      = phi + dphisdp_co2*c_co2 ;
    double dn_h2co3sdp_co2    = phi*dc_h2co3sdp_co2 + dphisdp_co2*c_h2co3 ;
    double dn_hco3sdp_co2     = phi*dc_hco3sdp_co2 + dphisdp_co2*c_hco3 ;
    double dn_co3sdp_co2      = phi*dc_co3sdp_co2 + dphisdp_co2*c_co3 ;
    double dn_casdp_co2       = phi*dc_casdp_co2 + dphisdp_co2*c_ca ;
    double dn_cahco3sdp_co2   = phi*dc_cahco3sdp_co2 + dphisdp_co2*c_cahco3 ;
    double dn_cah3sio4sdp_co2 = phi*dc_cah3sio4sdp_co2 + dphisdp_co2*c_cah3sio4 ;
    double dn_h3sio4sdp_co2   = phi*dc_h3sio4sdp_co2 + dphisdp_co2*c_h3sio4 ;
    double dn_h4sio4sdp_co2   = phi*dc_h4sio4sdp_co2 + dphisdp_co2*c_h4sio4 ;
    double dn_h2sio4sdp_co2   = phi*dc_h2sio4sdp_co2 + dphisdp_co2*c_h2sio4 ;
    double dn_cah2sio4sdp_co2 = phi*dc_cah2sio4sdp_co2 + dphisdp_co2*c_cah2sio4 ;
    double dn_caco3aqsdp_co2  = phi*dc_caco3aqsdp_co2 + dphisdp_co2*c_caco3aq ;
    double dn_caohsdp_co2     = phi*dc_caohsdp_co2 + dphisdp_co2*c_caoh ;
    double dn_ksdp_co2        = dphisdp_co2*c_k;
    double dn_clsdp_co2       = dphisdp_co2*c_cl;	*/
     
	/* ... par rapport a p_co2 */

	double dn_co2sdp_co2      = dphisdp_co2*c_co2*s_l + phi*ds_lsdp_co2*c_co2 + phi*s_l*dc_co2sdp_co2 \
							  + dphisdp_co2*gc_co2*s_g + phi*ds_gsdp_co2*gc_co2 + phi*s_g*dgc_co2sdp_co2 ;
	double dn_hco3sdp_co2     = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_hco3 + phi*s_l*dc_hco3sdp_co2 ;
    double dn_h2co3sdp_co2    = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_h2co3 + phi*s_l*dc_h2co3sdp_co2 ;
    double dn_co3sdp_co2      = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_co3 + phi*s_l*dc_co3sdp_co2  ;
    double dn_casdp_co2       = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_ca + phi*s_l*dc_casdp_co2  ;
    double dn_cahco3sdp_co2   = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_cahco3 + phi*s_l*dc_cahco3sdp_co2 ;
    double dn_h3sio4sdp_co2   = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_h3sio4 + phi*s_l*dc_h3sio4sdp_co2 ;
    double dn_h4sio4sdp_co2   = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_h4sio4 + phi*s_l*dc_h4sio4sdp_co2 ;
    double dn_cah3sio4sdp_co2 = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_cah3sio4 + phi*s_l*dc_cah3sio4sdp_co2 ;
    double dn_h2sio4sdp_co2   = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_h2sio4 + phi*s_l*dc_h2sio4sdp_co2 ;       
    double dn_cah2sio4sdp_co2 = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_cah2sio4 + phi*s_l*dc_cah2sio4sdp_co2 ;    
    double dn_caco3aqsdp_co2  = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_caco3aq + phi*s_l*dc_caco3aqsdp_co2 ;        
    double dn_caohsdp_co2     = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_caoh + phi*s_l*dc_caohsdp_co2 ;       
    double dn_ksdp_co2        = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_k ;             
    double dn_clsdp_co2       = (dphisdp_co2*s_l + phi*ds_lsdp_co2)*c_cl ; 
    
    	/* ... par rapport a p_l */
    double dc_co2sdp_l      = dproductco2_h2o(p_co2,p_l,temperature);
    double dgc_h2osdp_l 	= dproductgh2o_h2o(p_co2,p_l,temperature);
    double dc_ohsdp_l       = dc_co2sdp_l*(c_oh2 - c_oh)/dc_co2 ;
    double dc_hsdp_l        = -c_h*dc_ohsdp_l/c_oh ;
    double dc_h2co3sdp_l    = dc_co2sdp_l*K_h2co3 ;
    double dc_hco3sdp_l     = c_hco3*(dc_h2co3sdp_l/c_h2co3 + dc_ohsdp_l/c_oh) ;
    double dc_co3sdp_l      = c_co3*(dc_hco3sdp_l/c_hco3 + dc_ohsdp_l/c_oh) ;		  

    double dP_CCsdp_l       = (z_co2 < 1.) ? dc_co2sdp_l*P_CC/c_co2 : 0. ;
    
    double dc_casdp_l       = (dP_CCsdp_l - c_ca*dc_co3sdp_l)/c_co3 ;
    double dc_cahco3sdp_l   = K_cahco3*(dc_casdp_l*c_hco3 + c_ca*dc_hco3sdp_l) ;

    double dc_h4sio4sdp_l   = NEGEXP(zn_si_s)*dc_co2sdp_l*(x_s2 - x_s )* K_Sam/dc_co2 ;
    double dc_h3sio4sdp_l   = c_h3sio4*(dc_h4sio4sdp_l/c_h4sio4 - dc_hsdp_l/c_h) ;
    double dc_cah3sio4sdp_l = K_cah3sio4*(dc_casdp_l*c_h3sio4 + c_ca*dc_h3sio4sdp_l) ;
    double dc_h2sio4sdp_l   = K_h2sio4*(dc_h3sio4sdp_l*c_oh + dc_ohsdp_l*c_h3sio4) ;
    double dc_cah2sio4sdp_l = K_cah2sio4*(dc_h2sio4sdp_l*c_ca	+ dc_casdp_l*c_h2sio4) ;
    double dc_caco3aqsdp_l  = K_caco3aq*(dc_co3sdp_l*c_ca + dc_casdp_l*c_co3) ;				                  
    double dc_caohsdp_l     = K_caoh*(dc_casdp_l*c_oh + dc_ohsdp_l*c_ca) ;			                  	
		                    
    double dn_ch_cisdp_l   = dc_co2sdp_l*n_ch_ci*(-dt/t_ch)/c_co2 ;
    double dn_cc_cisdp_l   = dc_co2sdp_l*n_cc_ci*(dt/t_cc)/c_co2 ;


    double dn_ch_eqsdp_l   = -dn_cc_cisdp_l ;
    double dn_cc_eqsdp_l   = -dn_ch_cisdp_l ;
    
    double dn_chsdp_l   = (z_co2 <= 1) ? dn_ch_eqsdp_l  : dn_ch_cisdp_l ;
    double dn_ccsdp_l   = (z_co2 >  1) ? dn_cc_eqsdp_l  : dn_cc_cisdp_l ;
    double dx_cshsdp_l  = dc_co2sdp_l*(x_csh2- x_csh) / dc_co2 ;
    double dv_cshsdp_l  = dc_co2sdp_l*(v_csh2- v_csh) / dc_co2 ;
    double dz_cshsdp_l  = dc_co2sdp_l*(z_csh2- z_csh) / dc_co2 ;

    double dn_ca_ssdp_l = dx_cshsdp_l * MAX(zn_si_s,0.) ;
    double dn_si_ssdp_l = 0. ;

    double dphisdp_l        = - V_CH*dn_chsdp_l - dv_cshsdp_l*n_si_s - V_CC*dn_ccsdp_l ;
    double dc_h2osdp_l      = -(dc_co2sdp_l*v_co2 + v_h2co3*dc_h2co3sdp_l + v_hco3*dc_hco3sdp_l + v_co3*dc_co3sdp_l + v_h*dc_hsdp_l \
								+ v_oh*dc_ohsdp_l + v_ca*dc_casdp_l + v_caoh*dc_caohsdp_l + v_cahco3*dc_cahco3sdp_l + v_caco3aq*dc_caco3aqsdp_l \
								+ v_h3sio4*dc_h3sio4sdp_l + v_h4sio4*dc_h4sio4sdp_l + v_h2sio4*dc_h2sio4sdp_l + v_cah2sio4*dc_cah2sio4sdp_l \
								+ v_cah3sio4*dc_cah3sio4sdp_l)/v_h2o ;
	double drho_lsdp_l	  = M_H2O*dc_h2osdp_l + dc_co2sdp_l*M_CO2 + M_H2CO3*dc_h2co3sdp_l + M_HCO3*dc_hco3sdp_l + M_CO3*dc_co3sdp_l + M_H*dc_hsdp_l \
								+ M_OH*dc_ohsdp_l + M_Ca*dc_casdp_l + M_CaOH*dc_caohsdp_l + M_CaHCO3*dc_cahco3sdp_l + M_CaCO3*dc_caco3aqsdp_l \
								+ M_H3SiO4*dc_h3sio4sdp_l + M_H4SiO4*dc_h4sio4sdp_l + M_H2SiO4*dc_h2sio4sdp_l + M_CaH2SiO4*dc_cah2sio4sdp_l \
								+ M_CaH3SiO4*dc_cah3sio4sdp_l ;  			                    
    	
    double ds_lsdp_l        = 0.- dcourbe(p_c,el.mat->cb[0]) ;
	double ds_gsdp_l        = dcourbe(p_c,el.mat->cb[0]) ;		

	double dn_co2sdp_l      = dphisdp_l*c_co2*s_l + phi*ds_lsdp_l*c_co2 + phi*s_l*dc_co2sdp_l \
							  + dphisdp_l*gc_co2*s_g + phi*ds_gsdp_l*gc_co2  ;
	double dn_hco3sdp_l     = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_hco3 + phi*s_l*dc_hco3sdp_l ;
    double dn_h2co3sdp_l    = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_h2co3 + phi*s_l*dc_h2co3sdp_l ;
    double dn_co3sdp_l      = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_co3 + phi*s_l*dc_co3sdp_l  ;
    double dn_casdp_l       = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_ca + phi*s_l*dc_casdp_l  ;
    double dn_cahco3sdp_l   = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_cahco3 + phi*s_l*dc_cahco3sdp_l ;
    double dn_h3sio4sdp_l   = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_h3sio4 + phi*s_l*dc_h3sio4sdp_l ;
    double dn_h4sio4sdp_l   = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_h4sio4 + phi*s_l*dc_h4sio4sdp_l ;
    double dn_cah3sio4sdp_l = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_cah3sio4 + phi*s_l*dc_cah3sio4sdp_l ;
    double dn_h2sio4sdp_l   = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_h2sio4 + phi*s_l*dc_h2sio4sdp_l ;       
    double dn_cah2sio4sdp_l = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_cah2sio4 + phi*s_l*dc_cah2sio4sdp_l ;    
    double dn_caco3aqsdp_l  = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_caco3aq + phi*s_l*dc_caco3aqsdp_l ;        
    double dn_caohsdp_l     = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_caoh + phi*s_l*dc_caohsdp_l ;       
    double dn_ksdp_l        = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_k ;             
    double dn_clsdp_l       = (dphisdp_l*s_l + phi*ds_lsdp_l)*c_cl ; 
    		                    
    /* par rapport a zn_si_s */
    double dzn_si_s             = ((zn_si_s > ZN_Si_Sn(i)) ? 1 : -1)*1.e-4 ; 
    double zn_si_s2             = zn_si_s + dzn_si_s ;
    
    double c_oh1              = concentration_oh(c_co2,zn_si_s2,zn_ca_s,c_k,c_cl,el)  ;
    double dc_ohsdzn_si_s       = (c_oh1 - c_oh)/dzn_si_s ;
    double dc_hsdzn_si_s        = - c_h*dc_ohsdzn_si_s/c_oh ;
    double dc_hco3sdzn_si_s     = c_hco3*(dc_ohsdzn_si_s/c_oh) ;
    double dc_co3sdzn_si_s      = c_co3*(dc_hco3sdzn_si_s/c_hco3 \
								+ dc_ohsdzn_si_s/c_oh) ;
    double dc_casdzn_si_s       = - c_ca*dc_co3sdzn_si_s/c_co3 ;
    double dc_cahco3sdzn_si_s   = K_cahco3*(dc_casdzn_si_s*c_hco3 \
								+ c_ca*dc_hco3sdzn_si_s) ;

    double x_s3         	  = product_Sam(c_co2,zn_ca_s,el) ;

    double dc_h4sio4sdzn_si_s   = (NEGEXP(zn_si_s2)*x_s3 *K_Sam- NEGEXP(zn_si_s)*x_s *K_Sam)/dzn_si_s ;
    double dc_h3sio4sdzn_si_s   = c_h3sio4*(dc_h4sio4sdzn_si_s/c_h4sio4 \
								- dc_hsdzn_si_s/c_h) ;
    double dc_cah3sio4sdzn_si_s = K_cah3sio4*(dc_casdzn_si_s*c_h3sio4 \
								+ c_ca*dc_h3sio4sdzn_si_s) ;
    double dc_h2sio4sdzn_si_s   = K_h2sio4*(dc_h3sio4sdzn_si_s*c_oh	\
			                    + dc_ohsdzn_si_s*c_h3sio4) ;
    double dc_cah2sio4sdzn_si_s = K_cah2sio4*(dc_h2sio4sdzn_si_s*c_ca	\
			                    + dc_casdzn_si_s*c_h2sio4) ;
    double dc_caco3aqsdzn_si_s  = K_caco3aq*(dc_co3sdzn_si_s*c_ca \
			                    + dc_casdzn_si_s*c_co3) ;				                  
    double dc_caohsdzn_si_s     = K_caoh*(dc_casdzn_si_s*c_oh \
			                    + dc_ohsdzn_si_s*c_ca) ;				      
    
	double dn_si_ssdzn_si_s 	  = (zn_si_s > 0) ? 1. : 0. ;
    double dn_ca_ssdzn_si_s 	  = x_csh*dn_si_ssdzn_si_s;
    
    double dphisdzn_si_s   	  = -v_csh*dn_si_ssdzn_si_s ;

    double dn_co2sdzn_si_s      = dphisdzn_si_s*(c_co2*s_l+gc_co2*s_g) ;
    double dn_h2co3sdzn_si_s    = s_l*dphisdzn_si_s*c_h2co3 ;
    double dn_hco3sdzn_si_s     = s_l*phi*dc_hco3sdzn_si_s + s_l*dphisdzn_si_s*c_hco3 ;
    double dn_co3sdzn_si_s      = s_l*phi*dc_co3sdzn_si_s + s_l*dphisdzn_si_s*c_co3 ;
    double dn_casdzn_si_s       = s_l*phi*dc_casdzn_si_s + s_l*dphisdzn_si_s*c_ca ;
    double dn_cahco3sdzn_si_s   = s_l*phi*dc_cahco3sdzn_si_s + s_l*dphisdzn_si_s*c_cahco3 ;
    double dn_cah3sio4sdzn_si_s = s_l*phi*dc_cah3sio4sdzn_si_s + s_l*dphisdzn_si_s*c_cah3sio4 ;
    double dn_h3sio4sdzn_si_s   = s_l*phi*dc_h3sio4sdzn_si_s + s_l*dphisdzn_si_s*c_h3sio4 ;
    double dn_h4sio4sdzn_si_s   = s_l*phi*dc_h4sio4sdzn_si_s + s_l*dphisdzn_si_s*c_h4sio4 ;
    double dn_h2sio4sdzn_si_s   = s_l*phi*dc_h2sio4sdzn_si_s + s_l*dphisdzn_si_s*c_h2sio4 ;
    double dn_cah2sio4sdzn_si_s = s_l*phi*dc_cah2sio4sdzn_si_s + s_l*dphisdzn_si_s*c_cah2sio4 ;
    double dn_caco3aqsdzn_si_s  = s_l*phi*dc_caco3aqsdzn_si_s + s_l*dphisdzn_si_s*c_caco3aq ;
    double dn_caohsdzn_si_s     = s_l*phi*dc_caohsdzn_si_s + s_l*dphisdzn_si_s*c_caoh ;
    double dn_ksdzn_si_s        = s_l*dphisdzn_si_s*c_k;
    double dn_clsdzn_si_s       = s_l*dphisdzn_si_s*c_cl;
    double dc_h2osdzn_si_s      = -(v_hco3*dc_hco3sdzn_si_s + v_co3*dc_co3sdzn_si_s + v_h*dc_hsdzn_si_s \
								+ v_oh*dc_ohsdzn_si_s + v_ca*dc_casdzn_si_s + v_caoh*dc_caohsdzn_si_s + v_cahco3*dc_cahco3sdzn_si_s + v_caco3aq*dc_caco3aqsdzn_si_s \
								+ v_h3sio4*dc_h3sio4sdzn_si_s + v_h4sio4*dc_h4sio4sdzn_si_s + v_h2sio4*dc_h2sio4sdzn_si_s + v_cah2sio4*dc_cah2sio4sdzn_si_s \
								+ v_cah3sio4*dc_cah3sio4sdzn_si_s)/v_h2o ;
	double drho_lsdzn_si_s	  = M_H2O*dc_h2osdzn_si_s + M_HCO3*dc_hco3sdzn_si_s + M_CO3*dc_co3sdzn_si_s + M_H*dc_hsdzn_si_s \
								+ M_OH*dc_ohsdzn_si_s + M_Ca*dc_casdzn_si_s + M_CaOH*dc_caohsdzn_si_s + M_CaHCO3*dc_cahco3sdzn_si_s + M_CaCO3*dc_caco3aqsdzn_si_s \
								+ M_H3SiO4*dc_h3sio4sdzn_si_s + M_H4SiO4*dc_h4sio4sdzn_si_s + M_H2SiO4*dc_h2sio4sdzn_si_s + M_CaH2SiO4*dc_cah2sio4sdzn_si_s \
								+ M_CaH3SiO4*dc_cah3sio4sdzn_si_s ; 
    
    /* par rapport a zn_ca_s */
    /* double dzn_ca_s             = ((zn_ca_s > 0.) ? 1 : -1)*1.e-2 ; */
    double dzn_ca_s             = 1.e-6*((zn_ca_s > ZN_Ca_Sn(i)) ? 1 : -1) ;
    double zn_ca_s2             = zn_ca_s + dzn_ca_s ;
    double c_oh3               = concentration_oh(c_co2,zn_si_s,zn_ca_s2,c_k,c_cl,el) ;
    double dc_ohsdzn_ca_s       = (c_oh3 - c_oh)/dzn_ca_s ;
    double dc_hsdzn_ca_s        = - c_h*dc_ohsdzn_ca_s/c_oh ;
    double dc_hco3sdzn_ca_s     = c_hco3*(dc_ohsdzn_ca_s/c_oh) ;
    double dc_co3sdzn_ca_s      = c_co3*(dc_hco3sdzn_ca_s/c_hco3 \
								 + dc_ohsdzn_ca_s/c_oh) ;

    double dP_CCsdzn_ca_s       = (zn_ca_s < 0) ? P_CC : 0 ;
    double dc_casdzn_ca_s       = (dP_CCsdzn_ca_s - c_ca*dc_co3sdzn_ca_s)/c_co3 ;
    double dc_cahco3sdzn_ca_s   = K_cahco3*(dc_casdzn_ca_s*c_hco3 \
								 + c_ca*dc_hco3sdzn_ca_s) ;
    
    double zp_ch4			   = NEGEXP(zn_ca_s2)/MAX(z_co2,1.) ;
    double x_s4         	   = product_Sam(c_co2,zn_ca_s2,el) ;
    double x_csh4 			   = X_CSH(zp_ch4) ;
    double v_csh4 			   = V_CSH(zp_ch4) ;
    double z_csh4 			   = Z_CSH(zp_ch4);

    double dc_h4sio4sdzn_ca_s   = NEGEXP(zn_si_s)*(x_s4 - x_s)*K_Sam/dzn_ca_s ;
    double dc_h3sio4sdzn_ca_s   = c_h3sio4*(dc_h4sio4sdzn_ca_s/c_h4sio4 \
								 - dc_hsdzn_ca_s/c_h) ;
    double dc_cah3sio4sdzn_ca_s = K_cah3sio4*(dc_casdzn_ca_s*c_h3sio4 \
								 + c_ca*dc_h3sio4sdzn_ca_s) ;
    double dc_h2sio4sdzn_ca_s   = K_h2sio4*(dc_h3sio4sdzn_ca_s*c_oh	\
			                     + dc_ohsdzn_ca_s*c_h3sio4) ;
    double dc_cah2sio4sdzn_ca_s = K_cah2sio4*(dc_h2sio4sdzn_ca_s*c_ca	\
			                     + dc_casdzn_ca_s*c_h2sio4) ;
    double dc_caco3aqsdzn_ca_s  = K_caco3aq*(dc_co3sdzn_ca_s*c_ca \
			                     + dc_casdzn_ca_s*c_co3) ;				                  
    double dc_caohsdzn_ca_s     = K_caoh*(dc_casdzn_ca_s*c_oh \
			                     + dc_ohsdzn_ca_s*c_ca) ;					       

    double dn_ch_cisdzn_ca_s    = 0. ;
    double dn_cc_cisdzn_ca_s    = 0. ;


    double dn_ch_eqsdzn_ca_s    = (zn_ca_s > 0.) ? n_ch_ref : 0. ;
    double dn_cc_eqsdzn_ca_s    = (zn_ca_s > 0.) ? n_cc_ref : 0. ;

    double dn_chsdzn_ca_s       = (z_co2 <= 1.) ? dn_ch_eqsdzn_ca_s  : dn_ch_cisdzn_ca_s ;
    double dn_ccsdzn_ca_s  	   = (z_co2 >  1.) ? dn_cc_eqsdzn_ca_s  : dn_cc_cisdzn_ca_s ;    
    double dx_cshsdzn_ca_s	   = (x_csh4 - x_csh) / dzn_ca_s;
    double dz_cshsdzn_ca_s	   = (z_csh4 - z_csh) / dzn_ca_s;
    double dv_cshsdzn_ca_s	   = (v_csh4 - v_csh) / dzn_ca_s;

    double dn_si_ssdzn_ca_s 	   = 0. ;
    double dn_ca_ssdzn_ca_s 	   = dn_chsdzn_ca_s + dn_ccsdzn_ca_s + dx_cshsdzn_ca_s*MAX(zn_si_s,0.) ;

    double dphisdzn_ca_s    	   = - V_CH*dn_chsdzn_ca_s - V_CC*dn_ccsdzn_ca_s - dv_cshsdzn_ca_s*n_si_s;

    double dn_co2sdzn_ca_s      = dphisdzn_ca_s*(c_co2*s_l+gc_co2*s_g) ;
    double dn_h2co3sdzn_ca_s    = s_l*dphisdzn_ca_s*c_h2co3 ;
    double dn_hco3sdzn_ca_s     = s_l*phi*dc_hco3sdzn_ca_s + s_l*dphisdzn_ca_s*c_hco3 ;
    double dn_co3sdzn_ca_s      = s_l*phi*dc_co3sdzn_ca_s + s_l*dphisdzn_ca_s*c_co3 ;
    double dn_casdzn_ca_s       = s_l*phi*dc_casdzn_ca_s + s_l*dphisdzn_ca_s*c_ca ;
    double dn_cahco3sdzn_ca_s   = s_l*phi*dc_cahco3sdzn_ca_s + s_l*dphisdzn_ca_s*c_cahco3 ;
    double dn_cah3sio4sdzn_ca_s = s_l*phi*dc_cah3sio4sdzn_ca_s + s_l*dphisdzn_ca_s*c_cah3sio4 ;
    double dn_h3sio4sdzn_ca_s   = s_l*phi*dc_h3sio4sdzn_ca_s + s_l*dphisdzn_ca_s*c_h3sio4 ;
    double dn_h4sio4sdzn_ca_s   = s_l*phi*dc_h4sio4sdzn_ca_s + s_l*dphisdzn_ca_s*c_h4sio4 ;
    double dn_h2sio4sdzn_ca_s   = s_l*phi*dc_h2sio4sdzn_ca_s + s_l*dphisdzn_ca_s*c_h2sio4 ;
    double dn_cah2sio4sdzn_ca_s = s_l*phi*dc_cah2sio4sdzn_ca_s + s_l*dphisdzn_ca_s*c_cah2sio4 ;
    double dn_caco3aqsdzn_ca_s  = s_l*phi*dc_caco3aqsdzn_ca_s + s_l*dphisdzn_ca_s*c_caco3aq ;
    double dn_caohsdzn_ca_s     = s_l*phi*dc_caohsdzn_ca_s + s_l*dphisdzn_ca_s*c_caoh ;
    double dn_ksdzn_ca_s        = s_l*dphisdzn_ca_s*c_k;
    double dn_clsdzn_ca_s       = s_l*dphisdzn_ca_s*c_cl;
    double dc_h2osdzn_ca_s      = -(v_hco3*dc_hco3sdzn_ca_s + v_co3*dc_co3sdzn_ca_s + v_h*dc_hsdzn_ca_s \
								+ v_oh*dc_ohsdzn_ca_s + v_ca*dc_casdzn_ca_s + v_caoh*dc_caohsdzn_ca_s + v_cahco3*dc_cahco3sdzn_ca_s + v_caco3aq*dc_caco3aqsdzn_ca_s \
								+ v_h3sio4*dc_h3sio4sdzn_ca_s + v_h4sio4*dc_h4sio4sdzn_ca_s + v_h2sio4*dc_h2sio4sdzn_ca_s + v_cah2sio4*dc_cah2sio4sdzn_ca_s \
								+ v_cah3sio4*dc_cah3sio4sdzn_ca_s)/v_h2o ;
	double drho_lsdzn_ca_s	   = M_H2O*dc_h2osdzn_ca_s + M_HCO3*dc_hco3sdzn_ca_s + M_CO3*dc_co3sdzn_ca_s + M_H*dc_hsdzn_ca_s \
								+ M_OH*dc_ohsdzn_ca_s + M_Ca*dc_casdzn_ca_s + M_CaOH*dc_caohsdzn_ca_s + M_CaHCO3*dc_cahco3sdzn_ca_s + M_CaCO3*dc_caco3aqsdzn_ca_s \
								+ M_H3SiO4*dc_h3sio4sdzn_ca_s + M_H4SiO4*dc_h4sio4sdzn_ca_s + M_H2SiO4*dc_h2sio4sdzn_ca_s + M_CaH2SiO4*dc_cah2sio4sdzn_ca_s \
								+ M_CaH3SiO4*dc_cah3sio4sdzn_ca_s ; 
    
    /* ... par rapport a c_k */
    double dc_k                = 0.4e-4 ;/*c_k*((c_k > X_Kn(i)) ? 1 : -1) ;*/
    double c_k2                = c_k + dc_k ;
    double c_oh4               = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k2,c_cl,el)  ;
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

    double dn_hco3sdc_k        = s_l*phi*dc_hco3sdc_k ;
    double dn_co3sdc_k         = s_l*phi*dc_co3sdc_k ;
    double dn_casdc_k          = s_l*phi*dc_casdc_k ;
    double dn_cahco3sdc_k      = s_l*phi*dc_cahco3sdc_k ;
    double dn_cah3sio4sdc_k    = s_l*phi*dc_cah3sio4sdc_k ;
    double dn_h3sio4sdc_k      = s_l*phi*dc_h3sio4sdc_k ;
    double dn_h2sio4sdc_k      = s_l*phi*dc_h2sio4sdc_k ;
    double dn_cah2sio4sdc_k    = s_l*phi*dc_cah2sio4sdc_k ;
    double dn_caco3aqsdc_k     = s_l*phi*dc_caco3aqsdc_k ;
    double dn_caohsdc_k        = s_l*phi*dc_caohsdc_k ;
    double dc_h2osdc_k         = -(v_hco3*dc_hco3sdc_k + v_co3*dc_co3sdc_k + v_h*dc_hsdc_k \
								+ v_oh*dc_ohsdc_k + v_ca*dc_casdc_k + v_caoh*dc_caohsdc_k + v_cahco3*dc_cahco3sdc_k + v_caco3aq*dc_caco3aqsdc_k \
								+ v_h3sio4*dc_h3sio4sdc_k + v_h2sio4*dc_h2sio4sdc_k + v_cah2sio4*dc_cah2sio4sdc_k \
								+ v_cah3sio4*dc_cah3sio4sdc_k + v_k)/v_h2o ;
	double drho_lsdc_k	       = M_H2O*dc_h2osdc_k + M_HCO3*dc_hco3sdc_k + M_CO3*dc_co3sdc_k + M_H*dc_hsdc_k \
								+ M_OH*dc_ohsdc_k + M_Ca*dc_casdc_k + M_CaOH*dc_caohsdc_k + M_CaHCO3*dc_cahco3sdc_k + M_CaCO3*dc_caco3aqsdc_k \
								+ M_H3SiO4*dc_h3sio4sdc_k + M_H2SiO4*dc_h2sio4sdc_k + M_CaH2SiO4*dc_cah2sio4sdc_k \
								+ M_CaH3SiO4*dc_cah3sio4sdc_k + M_K; 
    
    /* ... par rapport a c_cl */
    double dc_cl                = 0.4e-4 ;
    double c_cl2                = c_cl + dc_cl ;
    double c_oh5                = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl2,el) ;
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

    double dn_hco3sdc_cl        = s_l*phi*dc_hco3sdc_cl ;
    double dn_co3sdc_cl         = s_l*phi*dc_co3sdc_cl ;
    double dn_casdc_cl          = s_l*phi*dc_casdc_cl ;
    double dn_cahco3sdc_cl      = s_l*phi*dc_cahco3sdc_cl ;
    double dn_cah3sio4sdc_cl    = s_l*phi*dc_cah3sio4sdc_cl ;
    double dn_h3sio4sdc_cl      = s_l*phi*dc_h3sio4sdc_cl ;
    double dn_h2sio4sdc_cl      = s_l*phi*dc_h2sio4sdc_cl ;
    double dn_cah2sio4sdc_cl    = s_l*phi*dc_cah2sio4sdc_cl ;
    double dn_caco3aqsdc_cl     = s_l*phi*dc_caco3aqsdc_cl ;
    double dn_caohsdc_cl        = s_l*phi*dc_caohsdc_cl ;
    double dc_h2osdc_cl         = -(v_hco3*dc_hco3sdc_cl + v_co3*dc_co3sdc_cl + v_h*dc_hsdc_cl \
								+ v_oh*dc_ohsdc_cl + v_ca*dc_casdc_cl + v_caoh*dc_caohsdc_cl + v_cahco3*dc_cahco3sdc_cl + v_caco3aq*dc_caco3aqsdc_cl \
								+ v_h3sio4*dc_h3sio4sdc_cl + v_h2sio4*dc_h2sio4sdc_cl + v_cah2sio4*dc_cah2sio4sdc_cl \
								+ v_cah3sio4*dc_cah3sio4sdc_cl + v_cl)/v_h2o ;
	double drho_lsdc_cl	        = M_H2O*dc_h2osdc_cl + M_HCO3*dc_hco3sdc_cl + M_CO3*dc_co3sdc_cl + M_H*dc_hsdc_cl \
								+ M_OH*dc_ohsdc_cl + M_Ca*dc_casdc_cl + M_CaOH*dc_caohsdc_cl + M_CaHCO3*dc_cahco3sdc_cl + M_CaCO3*dc_caco3aqsdc_cl \
								+ M_H3SiO4*dc_h3sio4sdc_cl + M_H2SiO4*dc_h2sio4sdc_cl + M_CaH2SiO4*dc_cah2sio4sdc_cl \
								+ M_CaH3SiO4*dc_cah3sio4sdc_cl + M_Cl ;     

    j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_CO2+j)   += volume[i]*(dn_co2sdp_co2 + dn_h2co3sdp_co2 + dn_hco3sdp_co2 + dn_co3sdp_co2 + dn_cahco3sdp_co2 \
                          + dn_caco3aqsdp_co2 + dn_ccsdp_co2 ) ;
    K(E_C+j,I_P_l+j)     += volume[i]*(dn_co2sdp_l + dn_h2co3sdp_l + dn_hco3sdp_l + dn_co3sdp_l + dn_cahco3sdp_l \
                          + dn_caco3aqsdp_l + dn_ccsdp_l) ;                 
    K(E_C+j,I_Ca_S+j)  += volume[i]*(dn_co2sdzn_ca_s + dn_h2co3sdzn_ca_s + dn_hco3sdzn_ca_s + dn_co3sdzn_ca_s + dn_cahco3sdzn_ca_s \
                          + dn_caco3aqsdzn_ca_s + dn_ccsdzn_ca_s) ;
    K(E_C+j,I_Si_S+j)   += volume[i]*(dn_co2sdzn_si_s + dn_h2co3sdzn_si_s + dn_hco3sdzn_si_s + dn_co3sdzn_si_s + dn_cahco3sdzn_si_s \
                          + dn_caco3aqsdzn_si_s) ;
    K(E_C+j,I_K+j)     += volume[i]*(dn_hco3sdc_k + dn_co3sdc_k + dn_cahco3sdc_k + dn_caco3aqsdc_k) ; 
    K(E_C+j,I_Cl+j)    += volume[i]*(dn_hco3sdc_cl + dn_co3sdc_cl + dn_cahco3sdc_cl + dn_caco3aqsdc_cl) ;                          
    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_CO2+j)  += volume[i]*(dn_casdp_co2 + dn_cahco3sdp_co2 + dn_cah3sio4sdp_co2  \
                          + dn_cah2sio4sdp_co2 + dn_caco3aqsdp_co2 + dn_caohsdp_co2 + dn_ca_ssdp_co2) ;
    K(E_Ca+j,I_P_l+j)  += volume[i]*(dn_casdp_l + dn_cahco3sdp_l + dn_cah3sio4sdp_l  \
                          + dn_cah2sio4sdp_l + dn_caco3aqsdp_l + dn_caohsdp_l + dn_ca_ssdp_l) ;
    K(E_Ca+j,I_Ca_S+j) += volume[i]*(dn_casdzn_ca_s + dn_cahco3sdzn_ca_s + dn_cah3sio4sdzn_ca_s + dn_ca_ssdzn_ca_s \
                          + dn_cah2sio4sdzn_ca_s + dn_caco3aqsdzn_ca_s +dn_caohsdzn_ca_s ) ;
    K(E_Ca+j,I_Si_S+j)  += volume[i]*(dn_casdzn_si_s + dn_cahco3sdzn_si_s + dn_cah3sio4sdzn_si_s + dn_ca_ssdzn_si_s \
                          + dn_cah2sio4sdzn_si_s + dn_caco3aqsdzn_si_s +dn_caohsdzn_si_s ) ;
    K(E_Ca+j,I_K+j)    += volume[i]*(dn_casdc_k + dn_cahco3sdc_k + dn_cah3sio4sdc_k \
                          + dn_cah2sio4sdc_k + dn_caco3aqsdc_k +dn_caohsdc_k ) ;
    K(E_Ca+j,I_Cl+j)   += volume[i]*(dn_casdc_cl + dn_cahco3sdc_cl + dn_cah3sio4sdc_cl \
                          + dn_cah2sio4sdc_cl + dn_caco3aqsdc_cl +dn_caohsdc_cl ) ;
    /*
      Conservation de Si (silicium) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
    */
    K(E_Si+j,I_CO2+j)  += volume[i]*(dn_h3sio4sdp_co2 + dn_h4sio4sdp_co2 + dn_cah3sio4sdp_co2 \
                          + dn_h2sio4sdp_co2 + dn_cah2sio4sdp_co2+ dn_si_ssdp_co2) ;
    K(E_Si+j,I_P_l+j)  += volume[i]*(dn_h3sio4sdp_l + dn_h4sio4sdp_l + dn_cah3sio4sdp_l \
                          + dn_h2sio4sdp_l + dn_cah2sio4sdp_l + dn_si_ssdp_l) ;
    K(E_Si+j,I_Ca_S+j) += volume[i]*(dn_h3sio4sdzn_ca_s + dn_h4sio4sdzn_ca_s + dn_cah3sio4sdzn_ca_s + dn_si_ssdzn_ca_s \
                          + dn_h2sio4sdzn_ca_s + dn_cah2sio4sdzn_ca_s) ;
    K(E_Si+j,I_Si_S+j)  += volume[i]*(dn_h3sio4sdzn_si_s + dn_h4sio4sdzn_si_s + dn_cah3sio4sdzn_si_s + dn_si_ssdzn_si_s \
                          + dn_h2sio4sdzn_si_s + dn_cah2sio4sdzn_si_s) ;
    K(E_Si+j,I_K+j)    += volume[i]*(dn_h3sio4sdc_k + dn_cah3sio4sdc_k + dn_h2sio4sdc_k + dn_cah2sio4sdc_k) ;
    K(E_Si+j,I_Cl+j)   += volume[i]*(dn_h3sio4sdc_cl + dn_cah3sio4sdc_cl + dn_h2sio4sdc_cl + dn_cah2sio4sdc_cl) ;    
    /*
      Conservation de la charge  : div(w_q) = 0
    */
    /*
      Conservation de K (potassium)  : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(E_K+j,I_CO2+j)   += volume[i]*(dn_ksdp_co2) ; 
    K(E_K+j,I_P_l+j)   += volume[i]*(dn_ksdp_l) ; 
    K(E_K+j,I_Ca_S+j)  += volume[i]*(dn_ksdzn_ca_s) ; 
    K(E_K+j,I_Si_S+j)   += volume[i]*(dn_ksdzn_si_s) ;    
    K(E_K+j,I_K+j)     += volume[i]*(phi*1.) ;
    
    /*      Conservation de Cl (chloride)  : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
    K(E_Cl+j,I_CO2+j)   += volume[i]*(dn_clsdp_co2) ; 
    K(E_Cl+j,I_P_l+j)   += volume[i]*(dn_clsdp_l) ; 
    K(E_Cl+j,I_Ca_S+j)  += volume[i]*(dn_clsdzn_ca_s) ; 
    K(E_Cl+j,I_Si_S+j)   += volume[i]*(dn_clsdzn_si_s) ;    
    K(E_Cl+j,I_Cl+j)    += volume[i]*(phi*1.) ;
    
    /*      Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
    */
    K(E_mass+j,I_CO2+j)   += volume[i]*(s_l*rho_l*dphisdp_co2 + s_l*phi*drho_lsdp_co2 + phi*rho_l*ds_lsdp_co2 + M_CO2*dphisdp_co2*s_g*gc_co2 \
							+ M_CO2*phi*ds_gsdp_co2*gc_co2 + M_CO2*phi*s_g*dgc_co2sdp_co2+ M_CaOH2*dn_chsdp_co2 + M_CaCO3*dn_ccsdp_co2\
							+ M_H2Og*dphisdp_co2*s_g*gc_h2o + M_H2Og*phi*ds_gsdp_co2*gc_h2o + M_H2Og*phi*s_g*dgc_h2osdp_co2\
							+ n_si_s*(M_CaO*dx_cshsdp_co2 + M_H2O*dz_cshsdp_co2)) ;
	K(E_mass+j,I_P_l+j)   += volume[i]*(s_l*rho_l*dphisdp_l + s_l*phi*drho_lsdp_l + phi*rho_l*ds_lsdp_l + M_CO2*dphisdp_l*s_g*gc_co2 \
							+ M_CO2*phi*ds_gsdp_l*gc_co2 + M_CaOH2*dn_chsdp_l + M_CaCO3*dn_ccsdp_l\
							+ M_H2Og*dphisdp_l*s_g*gc_h2o + M_H2Og*phi*ds_gsdp_l*gc_h2o+ M_H2Og*phi*s_g*dgc_h2osdp_l \
							+ n_si_s*(M_CaO*dx_cshsdp_l + M_H2O*dz_cshsdp_l)) ;		    
    K(E_mass+j,I_Ca_S+j)  += volume[i]*(s_l*rho_l*dphisdzn_ca_s + s_l*phi*drho_lsdzn_ca_s + M_CaOH2*dn_chsdzn_ca_s + M_CaCO3*dn_ccsdzn_ca_s\
							+ n_si_s*(M_CaO*dx_cshsdzn_ca_s + M_H2O*dz_cshsdzn_ca_s)); 
    K(E_mass+j,I_Si_S+j)   += volume[i]*(s_l*rho_l*dphisdzn_si_s + s_l*phi*drho_lsdzn_si_s \
							+ dn_si_ssdzn_si_s*(M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)); 
    K(E_mass+j,I_K+j)     += volume[i]*(s_l*phi*drho_lsdc_k) ;
    K(E_mass+j,I_Cl+j)    += volume[i]*(s_l*phi*drho_lsdc_cl) ;
							 
    /* sauvegardes pour les termes de transport */
    Dc_hSDp_co2[i]        = dc_hsdp_co2 ;
    Dc_ohSDp_co2[i]       = dc_ohsdp_co2 ;
    Dc_h2co3SDp_co2[i]    = dc_h2co3sdp_co2 ;
    Dc_hco3SDp_co2[i]     = dc_hco3sdp_co2 ;
    Dc_co3SDp_co2[i]      = dc_co3sdp_co2 ;
    Dc_caSDp_co2[i]       = dc_casdp_co2 ;
    Dc_cahco3SDp_co2[i]   = dc_cahco3sdp_co2 ;
    Dc_h3sio4SDp_co2[i]   = dc_h3sio4sdp_co2 ;
    Dc_h4sio4SDp_co2[i]   = dc_h4sio4sdp_co2 ;
    Dc_cah3sio4SDp_co2[i] = dc_cah3sio4sdp_co2 ;
    Dc_h2sio4SDp_co2[i]   = dc_h2sio4sdp_co2 ;
    Dc_cah2sio4SDp_co2[i] = dc_cah2sio4sdp_co2 ;
    Dc_caco3aqSDp_co2[i]  = dc_caco3aqsdp_co2 ;
    Dc_caohSDp_co2[i]     = dc_caohsdp_co2 ;
    Dc_co2SDp_co2[i]      = dc_co2sdp_co2 ;

    Dc_hSDp_l[i]        = dc_hsdp_l ;
    Dc_ohSDp_l[i]       = dc_ohsdp_l ;
    Dc_h2co3SDp_l[i]    = dc_h2co3sdp_l ;
    Dc_hco3SDp_l[i]     = dc_hco3sdp_l ;
    Dc_co3SDp_l[i]      = dc_co3sdp_l ;
    Dc_caSDp_l[i]       = dc_casdp_l ;
    Dc_cahco3SDp_l[i]   = dc_cahco3sdp_l ;
    Dc_h3sio4SDp_l[i]   = dc_h3sio4sdp_l ;
    Dc_h4sio4SDp_l[i]   = dc_h4sio4sdp_l ;
    Dc_cah3sio4SDp_l[i] = dc_cah3sio4sdp_l ;
    Dc_h2sio4SDp_l[i]   = dc_h2sio4sdp_l ;
    Dc_cah2sio4SDp_l[i] = dc_cah2sio4sdp_l ;
    Dc_caco3aqSDp_l[i]  = dc_caco3aqsdp_l ;
    Dc_caohSDp_l[i]     = dc_caohsdp_l ;
    Dc_co2SDp_l[i]      = dc_co2sdp_l ;    
    
    Dc_hSDzn_si_s[i]        = dc_hsdzn_si_s ;
    Dc_ohSDzn_si_s[i]       = dc_ohsdzn_si_s ;
    Dc_hco3SDzn_si_s[i]     = dc_hco3sdzn_si_s ;
    Dc_co3SDzn_si_s[i]      = dc_co3sdzn_si_s ;   
    Dc_caSDzn_si_s[i]       = dc_casdzn_si_s ;
    Dc_cahco3SDzn_si_s[i]   = dc_cahco3sdzn_si_s ;
    Dc_h3sio4SDzn_si_s[i]   = dc_h3sio4sdzn_si_s ;
    Dc_h4sio4SDzn_si_s[i]   = dc_h4sio4sdzn_si_s ;
    Dc_cah3sio4SDzn_si_s[i] = dc_cah3sio4sdzn_si_s ;
    Dc_h2sio4SDzn_si_s[i]   = dc_h2sio4sdzn_si_s ;
    Dc_cah2sio4SDzn_si_s[i] = dc_cah2sio4sdzn_si_s ;
    Dc_caco3aqSDzn_si_s[i]  = dc_caco3aqsdzn_si_s ;
    Dc_caohSDzn_si_s[i]     = dc_caohsdzn_si_s ;  

    Dc_hSDzn_ca_s[i]        = dc_hsdzn_ca_s ;
    Dc_ohSDzn_ca_s[i]       = dc_ohsdzn_ca_s ;
    Dc_hco3SDzn_ca_s[i]     = dc_hco3sdzn_ca_s ;
    Dc_co3SDzn_ca_s[i]      = dc_co3sdzn_ca_s ;
    Dc_caSDzn_ca_s[i]       = dc_casdzn_ca_s ;      
    Dc_cahco3SDzn_ca_s[i]   = dc_cahco3sdzn_ca_s ;
    Dc_h3sio4SDzn_ca_s[i]   = dc_h3sio4sdzn_ca_s ;
    Dc_h4sio4SDzn_ca_s[i]   = dc_h4sio4sdzn_ca_s ;
    Dc_cah3sio4SDzn_ca_s[i] = dc_cah3sio4sdzn_ca_s ;
    Dc_h2sio4SDzn_ca_s[i]   = dc_h2sio4sdzn_ca_s ;
    Dc_cah2sio4SDzn_ca_s[i] = dc_cah2sio4sdzn_ca_s ;
    Dc_caco3aqSDzn_ca_s[i]  = dc_caco3aqsdzn_ca_s ;
    Dc_caohSDzn_ca_s[i]     = dc_caohsdzn_ca_s ;
    
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
  
  double trg_co2        = tr*KD_GCO2 ;
  double trg_h2o        = tr*KD_GH2O ;
  
  double trd_h2co3 	 	= tr*KD_H2CO3 ;
  double trd_hco3  	 	= tr*KD_HCO3 ;
  double trd_co3  	  	= tr*KD_CO3 ;
  double trd_oh    	 	= tr*KD_OH ;
  double trd_h     	 	= tr*KD_H ;
  double trd_ca   	  	= tr*KD_Ca ;
  double trd_m     	 	= tr*KD_M ;
  double trd_k     		= tr*KD_K ;
  double trd_cl    		= tr*KD_Cl ;
  double trd_caoh  	 	= tr*KD_CaOH ;
  double trd_cahco3	 	= tr*KD_CaHCO3 ;
  double trd_caco3aq	= tr*KD_CaCO3aq ;
  double trd_h3sio4 	= tr*KD_H3SiO4 ;
  double trd_h2sio4 	= tr*KD_H2SiO4 ;
  double trd_h4sio4 	= tr*KD_H4SiO4 ;
  double trd_cah2sio4	= tr*KD_CaH2SiO4 ;
  double trd_cah3sio4	= tr*KD_CaH3SiO4 ;
  double trd_co2        = tr*KD_CO2 ;

  
  /*
  double trd_na    = tr*KD_Na ;
  double trd_naoh  = tr*KD_NaOH ;
  double trd_nahco3= tr*KD_NaHCO3 ;
  double trd_naco3 = tr*KD_NaCO3 ;
  double trd_koh   = tr*KD_KOH ;
  double trd_caoh2aq= tr*KD_CaOH2aq ;
  */
  
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
  
  double tre_h		  = tr*Kpsi_H ;
  double tre_oh		  = tr*Kpsi_OH ;	
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
    c[i] = trd_co2 + trd_h2co3 + trd_hco3 + trd_co3 + trd_cahco3 + trd_caco3aq + trf_co2*Dc_co2SDp_l[i] + trf_h2co3*Dc_h2co3SDp_l[i] + trf_hco3*Dc_hco3SDp_l[i] + trf_co3*Dc_co3SDp_l[i] \
           + trf_cahco3*Dc_cahco3SDp_l[i] + trf_caco3aq*Dc_caco3aqSDp_l[i];
  }
  K(E_C,I_P_l)          += + c[0] ;
  K(E_C,I_P_l+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_P_l)      += - c[0] ;
  K(E_C+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trg_co2 + trf_co2*Dc_co2SDp_co2[i] + trf_h2co3*Dc_h2co3SDp_co2[i] + trf_hco3*Dc_hco3SDp_co2[i] + trf_co3*Dc_co3SDp_co2[i] \
           + trf_cahco3*Dc_cahco3SDp_co2[i] + trf_caco3aq*Dc_caco3aqSDp_co2[i];
  }
  K(E_C,I_CO2)          += + c[0] ;
  K(E_C,I_CO2+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CO2)      += - c[0] ;
  K(E_C+NEQ,I_CO2+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDzn_si_s[i] + trf_co3*Dc_co3SDzn_si_s[i] + trf_cahco3*Dc_cahco3SDzn_si_s[i] + trf_caco3aq*Dc_caco3aqSDzn_si_s[i] ;
  }
  K(E_C,I_Si_S)          += + c[0] ;
  K(E_C,I_Si_S+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_Si_S)      += - c[0] ;
  K(E_C+NEQ,I_Si_S+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3*Dc_hco3SDzn_ca_s[i] + trf_co3*Dc_co3SDzn_ca_s[i] + trf_cahco3*Dc_cahco3SDzn_ca_s[i] + trf_caco3aq*Dc_caco3aqSDzn_ca_s[i] ;
  }
  K(E_C,I_Ca_S)          += + c[0] ;
  K(E_C,I_Ca_S+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_Ca_S)      += - c[0] ;
  K(E_C+NEQ,I_Ca_S+NEQ)  += + c[1] ;
  
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
    c[i] = trd_ca + trd_caoh + trd_cahco3 + trd_caco3aq + trd_cah2sio4 + trd_cah3sio4 +trf_ca*Dc_caSDp_l[i] + trf_cahco3*Dc_cahco3SDp_l[i] + trf_cah3sio4*Dc_cah3sio4SDp_l[i] + trf_cah2sio4*Dc_cah2sio4SDp_l[i] \
    + trf_caoh*Dc_caohSDp_l[i] + trf_caco3aq*Dc_caco3aqSDp_l[i];
  }
  K(E_Ca,I_P_l)          += + c[0] ;
  K(E_Ca,I_P_l+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_P_l)      += - c[0] ;
  K(E_Ca+NEQ,I_P_l+NEQ)  += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDp_co2[i] + trf_cahco3*Dc_cahco3SDp_co2[i] + trf_cah3sio4*Dc_cah3sio4SDp_co2[i] + trf_cah2sio4*Dc_cah2sio4SDp_co2[i] \
    + trf_caoh*Dc_caohSDp_co2[i] + trf_caco3aq*Dc_caco3aqSDp_co2[i] ;
  }
  K(E_Ca,I_CO2)         += + c[0] ;
  K(E_Ca,I_CO2+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_CO2)     += - c[0] ;
  K(E_Ca+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_si_s[i] + trf_cahco3*Dc_cahco3SDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_caoh*Dc_caohSDzn_si_s[i] + trf_caco3aq*Dc_caco3aqSDzn_si_s[i] ;
  }
  K(E_Ca,I_Si_S)         += + c[0] ;
  K(E_Ca,I_Si_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Si_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dc_caSDzn_ca_s[i] + trf_cahco3*Dc_cahco3SDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_caoh*Dc_caohSDzn_ca_s[i] + trf_caco3aq*Dc_caco3aqSDzn_ca_s[i] ;
  }
  K(E_Ca,I_Ca_S)         += + c[0] ;
  K(E_Ca,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Ca+NEQ,I_Ca_S+NEQ) += + c[1] ;
  
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
    c[i] = trd_h3sio4 + trd_h2sio4 + trd_h4sio4 + trd_cah2sio4 + trd_cah3sio4 +trf_h3sio4*Dc_h3sio4SDp_l[i] + trf_h4sio4*Dc_h4sio4SDp_l[i] + trf_cah3sio4*Dc_cah3sio4SDp_l[i] + trf_cah2sio4*Dc_cah2sio4SDp_l[i] \
    + trf_h2sio4*Dc_h2sio4SDp_l[i];
  }
  K(E_Si,I_P_l)          += + c[0] ;
  K(E_Si,I_P_l+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_P_l)      += - c[0] ;
  K(E_Si+NEQ,I_P_l+NEQ)  += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDp_co2[i] + trf_h4sio4*Dc_h4sio4SDp_co2[i] + trf_cah3sio4*Dc_cah3sio4SDp_co2[i] + trf_cah2sio4*Dc_cah2sio4SDp_co2[i] \
    + trf_h2sio4*Dc_h2sio4SDp_co2[i] ;
  }
  K(E_Si,I_CO2)         += + c[0] ;
  K(E_Si,I_CO2+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_CO2)     += - c[0] ;
  K(E_Si+NEQ,I_CO2+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_si_s[i] + trf_h4sio4*Dc_h4sio4SDzn_si_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_Si,I_Si_S)         += + c[0] ;
  K(E_Si,I_Si_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_Si_S)     += - c[0] ;
  K(E_Si+NEQ,I_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] + trf_h4sio4*Dc_h4sio4SDzn_ca_s[i] + trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i]  + trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_Si,I_Ca_S)         += + c[0] ;
  K(E_Si,I_Ca_S+NEQ)     += - c[1] ;
  K(E_Si+NEQ,I_Ca_S)     += - c[0] ;
  K(E_Si+NEQ,I_Ca_S+NEQ) += + c[1] ;
  
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
    c[i] = z_h*trf_h*Dc_hSDp_l[i] + z_oh*trf_oh*Dc_ohSDp_l[i] + z_hco3*trf_hco3*Dc_hco3SDp_l[i] + z_co3*trf_co3*Dc_co3SDp_l[i] \
           + z_ca*trf_ca*Dc_caSDp_l[i] + z_cahco3*trf_cahco3*Dc_cahco3SDp_l[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDp_l[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDp_l[i] + z_caoh*trf_caoh*Dc_caohSDp_l[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDp_l[i];
  }
  K(E_q,I_P_l)           += + c[0] ;
  K(E_q,I_P_l+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_P_l)       += - c[0] ;
  K(E_q+NEQ,I_P_l+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDp_co2[i] + z_oh*trf_oh*Dc_ohSDp_co2[i] + z_hco3*trf_hco3*Dc_hco3SDp_co2[i] + z_co3*trf_co3*Dc_co3SDp_co2[i] \
           + z_ca*trf_ca*Dc_caSDp_co2[i] + z_cahco3*trf_cahco3*Dc_cahco3SDp_co2[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDp_co2[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDp_co2[i] + z_caoh*trf_caoh*Dc_caohSDp_co2[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDp_co2[i];
  }
  K(E_q,I_CO2)           += + c[0] ;
  K(E_q,I_CO2+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_CO2)       += - c[0] ;
  K(E_q+NEQ,I_CO2+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_si_s[i] + z_oh*trf_oh*Dc_ohSDzn_si_s[i] + z_hco3*trf_hco3*Dc_hco3SDzn_si_s[i] + z_co3*trf_co3*Dc_co3SDzn_si_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_si_s[i] + z_cahco3*trf_cahco3*Dc_cahco3SDzn_si_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_si_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_si_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_si_s[i];
  }
  K(E_q,I_Si_S)           += + c[0] ;
  K(E_q,I_Si_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Si_S)       += - c[0] ;
  K(E_q+NEQ,I_Si_S+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dc_hSDzn_ca_s[i] + z_oh*trf_oh*Dc_ohSDzn_ca_s[i] + z_hco3*trf_hco3*Dc_hco3SDzn_ca_s[i] + z_co3*trf_co3*Dc_co3SDzn_ca_s[i] \
           + z_ca*trf_ca*Dc_caSDzn_ca_s[i] + z_cahco3*trf_cahco3*Dc_cahco3SDzn_ca_s[i] + z_h3sio4*trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] \
           + z_cah3sio4*trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + z_caoh*trf_caoh*Dc_caohSDzn_ca_s[i] + z_h2sio4*trf_h2sio4*Dc_h2sio4SDzn_ca_s[i];
  }
  K(E_q,I_Ca_S)           += + c[0] ;
  K(E_q,I_Ca_S+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Ca_S)       += - c[0] ;
  K(E_q+NEQ,I_Ca_S+NEQ)   += + c[1] ;
  
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
    c[i] = trd_k ;
  }
  K(E_K,I_P_l)          += + c[0] ;
  K(E_K,I_P_l+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_P_l)      += - c[0] ;
  K(E_K+NEQ,I_P_l+NEQ)  += + c[1] ;
  
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
    c[i] = trd_cl ;
  }
  K(E_Cl,I_P_l)          += + c[0] ;
  K(E_Cl,I_P_l+NEQ)      += - c[1] ;
  K(E_Cl+NEQ,I_P_l)      += - c[0] ;
  K(E_Cl+NEQ,I_P_l+NEQ)  += + c[1] ;
  
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
  
  /*
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  K(E_mass,I_P_l)          += + trd_m ;
  K(E_mass,I_P_l+NEQ)      += - trd_m ;
  K(E_mass+NEQ,I_P_l)      += - trd_m ;
  K(E_mass+NEQ,I_P_l+NEQ)  += + trd_m ;
  
  for(i=0;i<2;i++){
    c[i] = M_H2Og*trg_h2o + M_CO2*trg_co2 ; /* + M_CO2*trf_co2*Dc_co2SDp_co2[i] ;*/
  }
  K(E_mass,I_CO2)         += + c[0] ;
  K(E_mass,I_CO2+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_CO2)     += - c[0] ;
  K(E_mass+NEQ,I_CO2+NEQ) += + c[1] ;

    
/*  for(i=0;i<2;i++){
    c[i] = M_H*tre_h + M_OH*tre_oh + M_HCO3*tre_hco3 + M_CO3*tre_co3 + M_Ca*tre_ca + M_CaHCO3*tre_cahco3 + M_CaH3SiO4*tre_cah3sio4 \
			+ M_H3SiO4*tre_h3sio4  + M_H2SiO4*tre_h2sio4 + M_CaOH*tre_caoh + M_K*tre_k + M_Cl*tre_cl ;
  }
  K(E_mass,I_psi)          += + c[0] ;
  K(E_mass,I_psi+NEQ)      += - c[1] ;
  K(E_mass+NEQ,I_psi)      += - c[0] ;
  K(E_mass+NEQ,I_psi+NEQ)  += + c[1] ;  
  
  for(i=0;i<2;i++){
    c[i] = M_H*trf_h*Dc_hSDp_co2[i] + M_OH*trf_oh*Dc_ohSDp_co2[i] + M_H2O*trg_h2o + M_CO2*trg_co2 + M_CO2*trf_co2*Dc_co2SDp_co2[i] + M_H2CO3*trf_h2co3*Dc_h2co3SDp_co2[i] + M_HCO3*trf_hco3*Dc_hco3SDp_co2[i] \
    + M_CO3*trf_co3*Dc_co3SDp_co2[i] + M_CaHCO3*trf_cahco3*Dc_cahco3SDp_co2[i] + M_Ca*trf_ca*Dc_caSDp_co2[i] \
    + M_CaH3SiO4*trf_cah3sio4*Dc_cah3sio4SDp_co2[i] + M_H3SiO4*trf_h3sio4*Dc_h3sio4SDp_co2[i] + M_H2SiO4*trf_h2sio4*Dc_h2sio4SDp_co2[i]\
    + M_CaOH*trf_caoh*Dc_caohSDp_co2[i] + M_H4SiO4*trf_h4sio4*Dc_h4sio4SDp_co2[i] + M_CaH2SiO4*trf_cah2sio4*Dc_cah2sio4SDp_co2[i] \
    + M_CaCO3aq*trf_caco3aq*Dc_caco3aqSDp_co2[i] ;
  }
  K(E_mass,I_CO2)         += + c[0] ;
  K(E_mass,I_CO2+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_CO2)     += - c[0] ;
  K(E_mass+NEQ,I_CO2+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = M_H*trf_h*Dc_hSDzn_ca_s[i] + M_OH*trf_oh*Dc_ohSDzn_ca_s[i] + M_HCO3*trf_hco3*Dc_hco3SDzn_ca_s[i] \
    + M_CO3*trf_co3*Dc_co3SDzn_ca_s[i] + M_CaHCO3*trf_cahco3*Dc_cahco3SDzn_ca_s[i] + M_Ca*trf_ca*Dc_caSDzn_ca_s[i] \
    + M_CaH3SiO4*trf_cah3sio4*Dc_cah3sio4SDzn_ca_s[i] + M_H3SiO4*trf_h3sio4*Dc_h3sio4SDzn_ca_s[i] + M_H2SiO4*trf_h2sio4*Dc_h2sio4SDzn_ca_s[i]\
    + M_CaOH*trf_caoh*Dc_caohSDzn_ca_s[i] + M_H4SiO4*trf_h4sio4*Dc_h4sio4SDzn_ca_s[i] + M_CaH2SiO4*trf_cah2sio4*Dc_cah2sio4SDzn_ca_s[i] \
    + M_CaCO3aq*trf_caco3aq*Dc_caco3aqSDzn_ca_s[i] ;
  }
  K(E_mass,I_Ca_S)         += + c[0] ;
  K(E_mass,I_Ca_S+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_Ca_S)     += - c[0] ;
  K(E_mass+NEQ,I_Ca_S+NEQ) += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = M_H*trf_h*Dc_hSDzn_si_s[i] + M_OH*trf_oh*Dc_ohSDzn_si_s[i] + M_HCO3*trf_hco3*Dc_hco3SDzn_si_s[i] \
    + M_CO3*trf_co3*Dc_co3SDzn_si_s[i] + M_CaHCO3*trf_cahco3*Dc_cahco3SDzn_si_s[i] + M_Ca*trf_ca*Dc_caSDzn_si_s[i] \
    + M_CaH3SiO4*trf_cah3sio4*Dc_cah3sio4SDzn_si_s[i] + M_H3SiO4*trf_h3sio4*Dc_h3sio4SDzn_si_s[i] + M_H2SiO4*trf_h2sio4*Dc_h2sio4SDzn_si_s[i]\
    + M_CaOH*trf_caoh*Dc_caohSDzn_si_s[i] + M_H4SiO4*trf_h4sio4*Dc_h4sio4SDzn_si_s[i] + M_CaH2SiO4*trf_cah2sio4*Dc_cah2sio4SDzn_si_s[i] \
    + M_CaCO3aq*trf_caco3aq*Dc_caco3aqSDzn_si_s[i] ;
  }
  K(E_mass,I_Si_S)         += + c[0] ;
  K(E_mass,I_Si_S+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_Si_S)     += - c[0] ;
  K(E_mass+NEQ,I_Si_S+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = M_K*trf_k + M_H*trf_h*Dc_hSDc_k[i] + M_OH*trf_oh*Dc_ohSDc_k[i] + M_HCO3*trf_hco3*Dc_hco3SDc_k[i] \
    + M_CO3*trf_co3*Dc_co3SDc_k[i] + M_CaHCO3*trf_cahco3*Dc_cahco3SDc_k[i] + M_Ca*trf_ca*Dc_caSDc_k[i] \
    + M_CaH3SiO4*trf_cah3sio4*Dc_cah3sio4SDc_k[i] + M_H3SiO4*trf_h3sio4*Dc_h3sio4SDc_k[i] + M_H2SiO4*trf_h2sio4*Dc_h2sio4SDc_k[i]\
    + M_CaOH*trf_caoh*Dc_caohSDc_k[i] + M_CaH2SiO4*trf_cah2sio4*Dc_cah2sio4SDc_k[i] \
    + M_CaCO3aq*trf_caco3aq*Dc_caco3aqSDc_k[i] ;
  }
  K(E_mass,I_K)         += + c[0] ;
  K(E_mass,I_K+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_K)     += - c[0] ;
  K(E_mass+NEQ,I_K+NEQ) += + c[1] ; 

for(i=0;i<2;i++){
    c[i] = M_Cl*trf_cl + M_H*trf_h*Dc_hSDc_cl[i] + M_OH*trf_oh*Dc_ohSDc_cl[i] + M_HCO3*trf_hco3*Dc_hco3SDc_cl[i] \
    + M_CO3*trf_co3*Dc_co3SDc_cl[i] + M_CaHCO3*trf_cahco3*Dc_cahco3SDc_cl[i] + M_Ca*trf_ca*Dc_caSDc_cl[i] \
    + M_CaH3SiO4*trf_cah3sio4*Dc_cah3sio4SDc_cl[i] + M_H3SiO4*trf_h3sio4*Dc_h3sio4SDc_cl[i] + M_H2SiO4*trf_h2sio4*Dc_h2sio4SDc_cl[i]\
    + M_CaOH*trf_caoh*Dc_caohSDc_cl[i] + M_CaH2SiO4*trf_cah2sio4*Dc_cah2sio4SDc_cl[i] \
    + M_CaCO3aq*trf_caco3aq*Dc_caco3aqSDc_cl[i] ;
  }
  K(E_mass,I_Cl)         += + c[0] ;
  K(E_mass,I_Cl+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_Cl)     += - c[0] ;
  K(E_mass+NEQ,I_Cl+NEQ) += + c[1] ; 
  
  for(i=0;i<2;i++){
    c[i] = M_H*trd_h + M_OH*trd_oh + M_CO2*trd_co2 + M_H2CO3*trd_h2co3+ M_HCO3*trd_hco3 \
    + M_CO3*trd_co3 + M_CaHCO3*trd_cahco3 + M_Ca*trd_ca \
    + M_CaH3SiO4*trd_cah3sio4 + M_H3SiO4*trd_h3sio4 + M_H2SiO4*trd_h2sio4\
    + M_CaOH*trd_caoh + M_H4SiO4*trd_h4sio4 + M_CaH2SiO4*trd_cah2sio4 \
    + M_CaCO3aq*trd_caco3aq + M_K*trd_k + M_Cl*trd_cl;
  }
  K(E_mass,I_P_l)         += + c[0] ;
  K(E_mass,I_P_l+NEQ)     += - c[1] ;
  K(E_mass+NEQ,I_P_l)     += - c[0] ;
  K(E_mass+NEQ,I_P_l+NEQ) += + c[1] ;*/
  

 /*
 						- M_CO2*KF_CO2*grd_co2 					- M_H2CO3*KF_H2CO3*grd_h2co3 \
						- M_HCO3*KF_HCO3*grd_hco3				- M_CO3*KF_CO3*grd_co3 \
						- M_CaHCO3*KF_CaHCO3*grd_cahco3 		- M_Ca*KF_Ca*grd_ca\
						- M_CaH3SiO4*KF_CaH3SiO4*grd_cah3sio4	- M_H3SiO4*KF_H3SiO4*grd_h3sio4\
						- M_H2SiO4*KF_H2SiO4*grd_h2sio4\
						- M_CaOH*KF_CaOH*grd_caoh				- M_H4SiO4*KF_H4SiO4*grd_h4sio4\
						- M_CaH2SiO4*KF_CaH2SiO4*grd_cah2sio4	- M_CaCO3aq*KF_CaCO3aq*grd_caco3aq\
						- M_K*KF_K*grd_k 						- M_Cl*KF_Cl*grd_cl; */ 
 }

#if (U_CO2 == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_CO2)     *= Ln10*P_CO2(0) ;
    K(i,I_CO2+NEQ) *= Ln10*P_CO2(1) ;
  }
#endif
#if (U_PL == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_P_l)     *= Ln10*P_l(0) ;
    K(i,I_P_l+NEQ) *= Ln10*P_l(1) ;
  }
#endif

  

  return(0) ;

#undef K
}


void rs590(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  R(0,E_mass) -= volume[0]*(M(0) - M_n(0)) + dt*surf*W_M ;
  R(1,E_mass) -= volume[1]*(M(1) - M_n(1)) - dt*surf*W_M ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;

#undef R
}


int so590(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
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
  phi0     		= el.mat->pr[pm("porosite")] ;
  c_co2_eq 		= el.mat->pr[pm("C_CO2_eq")] ;
  temperature 	= el.mat->pr[pm("temperature")] ;

  /* initialisation */
  nso = 46 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees */
  {
    /* molarites */
#if (U_CO2 == LOG_RHO)
    double p_co2      = exp(Ln10*param(u,h_s,el.nn,I_CO2)) ;
#else
    double p_co2      = param(u,h_s,el.nn,I_CO2) ;
#endif
#if (U_PL == LOG_RHO)
    double p_l      = exp(Ln10*param(u,h_s,el.nn,I_P_l)) ;
#else
    double p_l      = param(u,h_s,el.nn,I_P_l) ;
#endif
    double gc_co2     = product_gco2(p_co2,temperature) ;/*(co2 concentration in gas)*/
    double y_h2o 	  = producty_h2o(p_co2,p_l,temperature) ;
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ;
    double c_co2      = product_co2(p_co2,p_l,temperature);/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_ca_s     = param(u,h_s,el.nn,I_Ca_S) ;
    double c_k        = param(u,h_s,el.nn,I_K) ;
    double zn_si_s      = param(u,h_s,el.nn,I_Si_S);
    double c_cl       = param(u,h_s,el.nn,I_Cl) ;
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0 ; 
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ;     
    
    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;           
    double pk_ch      = (MIN(zn_ca_s,0) - log(MAX(z_co2,1.)))/log(10); /*log10(P_CH/K_CH);  ajouter3 */
    double zp_ch      = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ;
    double x_s        = product_Sam(c_co2,zn_ca_s,el);
    
    double c_h4sio4   = NEGEXP(zn_si_s)*product_Sam(c_co2,zn_ca_s,el) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el)  ;
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
    + z_cah3sio4*c_cah3sio4 + z_h2sio4*c_h2sio4 + z_caoh*c_caoh + z_k*c_k + z_cl*c_cl ;
  /* contenus solides */
    int j = (el.dim < dim) ? 0 : ((s[0] < (x[0][0] + x[1][0])*0.5) ? 0 : 1) ;
    double n_ch       = N_CH(j) ;
    double n_cc       = N_CC(j) ;
    double x_csh      = X_CSH(zp_ch) ;
    double v_csh      = V_CSH(zp_ch);
    double z_csh      = Z_CSH(zp_ch);
        /* solides */
    double n_si_s     = MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  /* porosite */
    double n_ch0      = N_CH0(j) ;
    double n_cc0      = N_CC0(j) ;

    double zp_ch0     = ZP_CH0(j) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double n_si_s0	  = N_Si_S0(j) ;
    double phi        = phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) + v_csh0*n_si_s0- v_csh*n_si_s;   
    double psi        = param(u,h_s,el.nn,I_psi) ;
    double c_h2o      = (1. - (c_co2*v_co2 + c_h2co3*v_h2co3 + c_hco3*v_hco3 + c_co3*v_co3 + c_h*v_h + c_oh*v_oh + c_ca*v_ca \
						+ c_caoh*v_caoh + v_cahco3*c_cahco3 + c_caco3aq*v_caco3aq \
						+ c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 + c_cah2sio4*v_cah2sio4  \
						+ c_cah3sio4*v_cah3sio4+ c_k*v_k + c_cl*v_cl))/v_h2o ;
						/* + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + \
						 * x_koh*v_koh + c_caoh2aq*v_caoh2aq \*/

    /* masse volumique liquide */
    double rho_l   = M_CO2*c_co2 + M_H*c_h + M_OH*c_oh + M_H2O*c_h2o + M_H2CO3*c_h2co3 + M_HCO3*c_hco3 + M_CO3*c_co3 + M_Ca*c_ca \
			+ M_CaOH*c_caoh + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 \
			+ M_H2SiO4*c_h2sio4 + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 + M_Cl*c_cl + M_K*c_k ;
			/* + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_K*x_k + M_KOH*x_koh + M_CaOH2aq*x_caoh2aq */ 
   
   /* masse totale */
    double M    = M_H2Og*phi*s_g*gc_h2o + M_CO2*phi*s_g*gc_co2 + s_l*phi*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
   
   double mu_co2   = product_muco2(p_co2,temperature) ; 

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
    strcpy(r[i].text,"zn_si_s") ; r[i].n = 1 ;
    r[i++].v[0] = zn_si_s ; 
    strcpy(r[i].text,"c_cl") ; r[i].n = 1 ;
    r[i++].v[0] = c_cl ;  
    strcpy(r[i].text,"M") ; r[i].n = 1 ;
    r[i++].v[0] = M ; 
    strcpy(r[i].text,"p_l") ; r[i].n = 1 ;
    r[i++].v[0] = p_l ; 
    strcpy(r[i].text,"c_h2o") ; r[i].n = 1 ;
    r[i++].v[0] = c_h2o ; 
    strcpy(r[i].text,"W_M") ; r[i].n = 1 ;
    r[i++].v[0] = W_M ;
    strcpy(r[i].text,"rho_l") ; r[i].n = 1 ;
    r[i++].v[0] = rho_l ;
    strcpy(r[i].text,"M_L") ; r[i].n = 1 ;
    r[i++].v[0] = s_l*phi*rho_l ;
    strcpy(r[i].text,"M_S") ; r[i].n = 1 ;
    r[i++].v[0] = M_CaOH2*n_ch + M_CaCO3*n_cc + (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
    strcpy(r[i].text,"p_co2") ; r[i].n = 1 ;
    r[i++].v[0] = p_co2 ;
    strcpy(r[i].text,"s_l") ; r[i].n = 1 ;
    r[i++].v[0] = s_l ;
    strcpy(r[i].text,"s_g") ; r[i].n = 1 ;
    r[i++].v[0] = s_g ;
    strcpy(r[i].text,"gc_co2") ; r[i].n = 1 ;
    r[i++].v[0] = gc_co2 ;
    strcpy(r[i].text,"x_csh") ; r[i].n = 1 ;
    r[i++].v[0] = x_csh ;
    strcpy(r[i].text,"z_csh") ; r[i].n = 1 ;
    r[i++].v[0] = z_csh ;
    strcpy(r[i].text,"v_csh") ; r[i].n = 1 ;
    r[i++].v[0] = v_csh ;
    strcpy(r[i].text,"temperature") ; r[i].n = 1 ;
    r[i++].v[0] = temperature ;
    strcpy(r[i].text,"mu_co2") ; r[i].n = 1 ;
    r[i++].v[0] = mu_co2 ;
    strcpy(r[i].text,"gc_h2o") ; r[i].n = 1 ;
    r[i++].v[0] = gc_h2o ;
    strcpy(r[i].text,"y_h2o") ; r[i].n = 1 ;
    r[i++].v[0] = y_h2o ;
    strcpy(r[i].text,"mu_l") ; r[i].n = 1 ;
    r[i++].v[0] = mu_l ;
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
  k_intl    = el.mat->pr[pm("k_intl")] ;
  k_intg    = el.mat->pr[pm("k_intg")] ;
  temperature 	= el.mat->pr[pm("temperature")] ;


  /* initialisation */
  for(i=0;i<NVE_TR;i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i=0;i<2;i++) {
    /* molarites */
    double p_co2	  = (P_CO2(i) > 0.) ? P_CO2(i) : P_CO2_eq ;
    double p_l     	  = P_l(i) ;
    double gc_co2     = product_gco2(p_co2,temperature) ;/*(co2 concentration in gas)*/
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ;    
    double c_co2      = product_co2(p_co2,p_l,temperature) ;/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s      = ZN_Si_S(i) ;
    double zn_ca_s     = ZN_Ca_S(i) ;
    double c_k        = X_K(i);                 
    double c_cl       = X_Cl(i); 
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0;
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ; 

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ; 
    double zp_ch      = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ;

    double c_h4sio4   = NEGEXP(zn_si_s)*product_Sam(c_co2,zn_ca_s,el) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el) ;
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
    double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;                /* ajouter */
    double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;            /* ajouter */
    
    double c_h2o      = (1. - (c_co2*v_co2 + c_h2co3*v_h2co3 + c_hco3*v_hco3 + c_co3*v_co3 + c_h*v_h + c_oh*v_oh + c_ca*v_ca \
						+ c_caoh*v_caoh + v_cahco3*c_cahco3 + c_caco3aq*v_caco3aq \
						+ c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 + c_cah2sio4*v_cah2sio4  \
						+ c_cah3sio4*v_cah3sio4+ c_k*v_k + c_cl*v_cl))/v_h2o ;
						/* + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 + x_k*v_k + \
						 * x_koh*v_koh + c_caoh2aq*v_caoh2aq \*/

    /* masse volumique liquide */
    double rho_l   = M_CO2*c_co2 + M_H*c_h + M_OH*c_oh + M_H2O*c_h2o + M_H2CO3*c_h2co3 + M_HCO3*c_hco3 + M_CO3*c_co3 + M_Ca*c_ca \
			+ M_CaOH*c_caoh + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 \
			+ M_H2SiO4*c_h2sio4 + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 + M_Cl*c_cl + M_K*c_k ;
			/* + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 + M_KOH*x_koh + M_CaOH2aq*x_caoh2aq */
			
    /* solides */
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double x_csh      = X_CSH(zp_ch) ;
    double v_csh      = V_CSH(zp_ch) ;
    double z_csh      = Z_CSH(zp_ch);
    
    /* solides */
    double n_si_s     = MAX(zn_si_s,0.) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;

    /* porosite */
    double n_cc0      = N_CC0(i) ;
    double n_ch0      = N_CH0(i) ;
    double zp_ch0     = ZP_CH0(i) ;
    double v_csh0     = V_CSH(zp_ch0) ;
    double n_si_s0	  = N_Si_S0(i) ;
    double phi        = MIN((phi0 + V_CH*(n_ch0 - n_ch) + V_CC*(n_cc0 - n_cc) + v_csh0*n_si_s0- v_csh*n_si_s),0.8);   
       
    /* permeabilite */
    double k_l  = (k_intl/mu_l)*courbe(p_c,el.mat->cb[1])*pow(phi/phi0,3.)*pow(((1.-phi0)/(1.-phi)),2.) ;
    /* permeabilite de gas */
    double mu_co2   = product_muco2(p_co2,temperature) ;
    double k_g  = (k_intg/mu_co2)*courbe(p_c,el.mat->cb[2])*pow(phi/phi0,3.)*pow(((1.-phi0)/(1.-phi)),2.) ;
    
    /* tortuosite gaz */
    /*double tau  = pow(phi,1/3)*pow(s_g,7/3) ;*/
    double tau  = pow(phi,1.74)*pow(s_g,3.20) ;
    /* tortuosite liquide */
    double iff    = 0.00029*exp(9.95*phi)/(1.+625.*pow((1-s_l),4.)) ;/*0.00029*exp(9.95*phi)/(1.+625.*pow((1-s_l),4.)) ;*/
    
    KD_Ca     	+= c_ca*k_l ;
    KD_H2CO3  	+= c_h2co3*k_l ;
    KD_HCO3   	+= c_hco3*k_l ;
    KD_CO3    	+= c_co3*k_l ;
    KD_OH     	+= c_oh*k_l ;
    KD_H      	+= c_h*k_l ;
    KD_K      	+= c_k*k_l ;
    KD_Cl       += c_cl*k_l ;
    KD_CaOH   	+= c_caoh*k_l ;
    KD_CaHCO3 	+= c_cahco3*k_l ;
    KD_CaCO3aq	+= c_caco3aq*k_l ;
    KD_H3SiO4 	+= c_h3sio4*k_l ;
    KD_H2SiO4 	+= c_h2sio4*k_l ;
    KD_H4SiO4 	+= c_h4sio4*k_l ;
    KD_CaH3SiO4 += c_cah3sio4*k_l ;
    KD_CaH2SiO4 += c_cah2sio4*k_l ;
    KD_M      	+= rho_l*k_l ;
    KD_CO2      += c_co2*k_l ;
    
    KD_GCO2     += gc_co2*k_g ;
    KD_GH2O     += gc_h2o*k_g ;
    /*
    KD_Na     	+= c_na*k_l ;
    KD_NaOH   	+= c_naoh*k_l ;
    KD_NaHCO3 	+= c_nahco3*k_l ;
    KD_NaCO3  	+= c_naco3*k_l ;
    KD_CaOH2aq	+= c_caoh2aq*k_l ;
    KD_KOH    	+= c_koh*k_l ; 
    */
    
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
  /*
    Donnees
  */
  temperature 	= el.mat->pr[pm("temperature")] ;  

  for(i=0;i<2;i++) {
    double p_co2	  = (P_CO2(i) > 0.) ? P_CO2(i) : P_CO2_eq ;
    double p_l     	  = P_l(i) ;
    double gc_co2     = product_gco2(p_co2,temperature) ;/*(co2 concentration in gas)*/
    double gc_h2o     = product_gh2o(p_co2,p_l,temperature) ; 
    double c_co2      = product_co2(p_co2,p_l,temperature) ;/*(co2 concentration in water)(C_CO2(i) > 0.) ? C_CO2(i) : c_co2_eq ;*/
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s      = ZN_Si_S(i) ;
    double zn_ca_s     = ZN_Ca_S(i) ;
    double c_k        = X_K(i);                 
    double c_cl       = X_Cl(i); 
    double p_g        = p_co2 ; 
    double p_c        = p_g - p_l - p_l0;
    double s_l        = courbe(p_c,el.mat->cb[0]) ;
    double s_g     	  = 1. - s_l ; 

    double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;    
 
    double c_h4sio4   = NEGEXP(zn_si_s)*product_Sam(c_co2,zn_ca_s,el) * K_Sam ;
    double c_oh       = concentration_oh(c_co2,zn_si_s,zn_ca_s,c_k,c_cl,el) ;
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
    double grd_p_g      = (P_CO2(1)      - P_CO2(0)         )/dx ;
    double grd_p_l      = (P_l(1)      - P_l(0)         )/dx ;
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
    double w_co2      = - KD_CO2*grd_p_l 		- KF_CO2*grd_co2            -KD_GCO2*grd_p_g ; 
    double w_h2co3    = - KD_H2CO3*grd_p_l 		- KF_H2CO3*grd_h2co3   ;
    double w_hco3     = - KD_HCO3*grd_p_l 		- KF_HCO3*grd_hco3          - Kpsi_HCO3*grd_psi     ;
    double w_co3      = - KD_CO3*grd_p_l 		- KF_CO3*grd_co3            - Kpsi_CO3*grd_psi      ;
    double w_cahco3   = - KD_CaHCO3*grd_p_l 	- KF_CaHCO3*grd_cahco3      - Kpsi_CaHCO3*grd_psi   ;
    double w_ca       = - KD_Ca*grd_p_l 		- KF_Ca*grd_ca              - Kpsi_Ca*grd_psi       ;
    double w_cah3sio4 = - KD_CaH3SiO4*grd_p_l 	- KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
    double w_h3sio4   = - KD_H3SiO4*grd_p_l 	- KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi   ;
    double w_h2sio4   = - KD_H2SiO4*grd_p_l 	- KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi   ;
    double w_caoh     = - KD_CaOH*grd_p_l 		- KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi     ;
    double w_h4sio4   = - KD_H4SiO4*grd_p_l 	- KF_H4SiO4*grd_h4sio4 ;
    double w_cah2sio4 = - KD_CaH2SiO4*grd_p_l 	- KF_CaH2SiO4*grd_cah2sio4 ;
    double w_caco3aq  = - KD_CaCO3aq*grd_p_l 	- KF_CaCO3aq*grd_caco3aq ;    
    double w_k        = - KD_K*grd_p_l 			- KF_K*grd_k                - Kpsi_K*grd_psi        ; 
    double w_cl       = - KD_Cl*grd_p_l 		- KF_Cl*grd_cl              - Kpsi_Cl*grd_psi        ; 
    double w_m     	  = - KD_M*grd_p_l 					- M_CO2*KD_GCO2*grd_p_g	- M_H2Og*KD_GH2O*grd_p_g ;
						/*- M_H*Kpsi_H*grd_psi			- M_OH*Kpsi_OH*grd_psi	\
						- M_HCO3*Kpsi_HCO3*grd_psi		- M_CO3*Kpsi_CO3*grd_psi	\
						- M_Ca*Kpsi_Ca*grd_psi			- M_CaHCO3*Kpsi_CaHCO3*grd_psi	\
						- M_H3SiO4*Kpsi_H3SiO4*grd_psi	- M_CaH3SiO4*Kpsi_CaH3SiO4*grd_psi	\
						- M_H2SiO4*Kpsi_H2SiO4*grd_psi	- M_CaOH*Kpsi_CaOH*grd_psi	\
						- M_Cl*Kpsi_Cl*grd_psi			- M_K*Kpsi_K*grd_psi \
						- M_H*KF_H*grd_h						- M_OH*KF_OH*grd_oh \
						- M_CO2*KF_CO2*grd_co2 					- M_H2CO3*KF_H2CO3*grd_h2co3 \
						- M_HCO3*KF_HCO3*grd_hco3				- M_CO3*KF_CO3*grd_co3 \
						- M_CaHCO3*KF_CaHCO3*grd_cahco3 		- M_Ca*KF_Ca*grd_ca\
						- M_CaH3SiO4*KF_CaH3SiO4*grd_cah3sio4	- M_H3SiO4*KF_H3SiO4*grd_h3sio4\
						- M_H2SiO4*KF_H2SiO4*grd_h2sio4\
						- M_CaOH*KF_CaOH*grd_caoh				- M_H4SiO4*KF_H4SiO4*grd_h4sio4\
						- M_CaH2SiO4*KF_CaH2SiO4*grd_cah2sio4	- M_CaCO3aq*KF_CaCO3aq*grd_caco3aq\
						\
						- M_K*KF_K*grd_k 						- M_Cl*KF_Cl*grd_cl;*/
						
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
    W_M     = w_m ;
  }
}


double concentration_oh(double c_co2,double zn_si_s,double zn_ca_s,double c_k,double c_cl,elem_t el)
/* on resout l'electroneutralie : SUM(z_i c_i) = 0
   racine de ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_co2, zn_si_s et zn_ca_s sont fixes */
  double c_co2_eq   = C_CO2_eq ;
  double z_co2      = c_co2/c_co2_eq ;
  /* les produits sont donc aussi fixes */
  double P_CC       = K_CC*MIN(z_co2,1.)*NEGEXP(zn_ca_s) ;
  
  double c_h4sio4   = NEGEXP(zn_si_s)*product_Sam(c_co2,zn_ca_s,el) * K_Sam ;
  /* double P_Sam      = K_Sam*pow(MAX(z_ch,1.),-X_CSH/Y_CSH)*pow(NEGEXP(zn_si_s),1/Y_CSH) ;
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
  double c =  z_k*c_k + z_cl*c_cl; /* z_k*c_k + z_cl*c_cl ;*/
  double d = z_oh*K_h2o + z_hco3*A_hco3 + z_h3sio4*A_h3sio4 ;
  double e = z_co3*A_co3 + z_h2sio4*A_h2sio4 ;

  double c_h = poly4(a,b,c,d,e) ;
      if(c_h < 0) {
      printf("c_h  = %e\n",c_h) ;
      printf("a  = %e\n",a) ;
      printf("b  = %e\n",b) ;
      printf("c  = %e\n",c) ;
      printf("d  = %e\n",d) ;
      printf("c_co2 = %e\n",c_co2) ;
      arret("concentration_oh : c_h<0") ;
    }
 
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

double poly4_sansc(double a,double b,double d,double e)
/* on resout ax^4 + bx^3 + dx + e = 0 */
{
  double err,tol = 1e-8 ;
  double x0 = pow(-d/a,1./3) ;
  double x  = MAX(x0,0.) ;
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

double product_muco2(double p_co2,double T)  
/* 	to produce co2 viscosity,after A.fenghour1997
	unit: Pa.s,t in K,P in Pa,T between 200-1500K
	pressure tp to 300MPa for T<1000K, 30MPa for T>1000K
	
  */
{
  double p_gas		= MAX(p_co2,1.e5) ; 		/*when p_co2 is smaller than 1e5Pa,mu_co2=product_muco2(1e5,T)*/
  double gc_co2 	= product_gco2(p_gas,T) ;    /*mol/L */
  double rho    	= 44.* gc_co2 ;				/*kg/m3*/
  double d11  		= 0.4071119e-2 ;
  double d21		= 0.7198037e-4 ;
  double d64 		= 0.2411697e-16 ;
  double d81		= 0.2971072e-22 ;
  double d82		= -0.1627888e-22 ;
  double T_red		= T/251.196 ;
  double logT_red	= log(T_red) ;
  double a0			= 0.235156 ;
  double a1			= -0.491266 ;
  double a2			= 5.211155e-2 ;
  double a3			= 5.347906e-2 ;
  double a4			= -1.537102e-2 ;
  double A_mu		= a0 + a1*logT_red + a2*pow(logT_red,2.) + a3*pow(logT_red,3.)  + a4*pow(logT_red,4.) ;
  double g_mu   	= exp(A_mu) ;
  double mu_zero	=1.00697*pow(T,0.5)/g_mu ; /*zero-density viscosity, uPa.s*/
  double mu_excess 	= d11*rho+ d21*pow(rho,2.) + d64*pow(rho,6.)/pow(T_red,3.) + d81*pow(rho,8.) + d82*pow(rho,8.)/T_red;/* in uPa.s*/
  double mu 		= mu_zero + mu_excess; /* uPa.s*/
  double mu_co2		= mu*1.e-6 ;           /* Pa.s*/
  return(mu_co2) ;
}

double product_Sam(double c_co2,double zn_ca_s,elem_t el)
{
  double z_co2 =  c_co2/c_co2_eq ;
  double zp_ch = NEGEXP(zn_ca_s)/MAX(z_co2,1.) ; 
  double x_sam = courbe(zp_ch,el.mat->cb[6]) ;
  return(x_sam) ;
}

double productx_co2(double pascal_co2,double p_l,double T) 
/*co2 mole fraction in h2o rich phase,after Nicolas Spycher 2003*/
{
  double p_co2     = pascal_co2/1.e5;   /*unit of pressure:bar*/
  double p_h2ol	   = MAX(p_l/1.e5,1.) ;
  double R         = 83.14472 ;   /*cm3*bar*K-1*mol-1*/
  double A_CO2     = 7.54e7 - 4.13e4*T ;/*bar*cm⁶*K0.5*mol^-2*/
  double B_CO2     = 27.80 ; /*cm³/mol*/
  double B_H2O     = 18.18 ; /*cm³/mol*/
  double A_H2OCO2  = 7.89e7 ;   /*bar*cm⁶*K0.5*mol^-2*/
  double av        = 1.0 ;
  double bv        = -(R*T/p_co2) ;
  double cv        = -(R*T*B_CO2/p_co2 - A_CO2/(p_co2*pow(T,0.5)) + B_CO2*B_CO2) ;
  double dv		   = -(A_CO2*B_CO2/(p_co2*pow(T,0.5)));
  double V         = cubic_equation(av,bv,cv,dv) ;           /*cm³/mol*/
  double Vav_gco2  = 32.6 ; /*cm³/mol*/
  double Vav_h2o   = 18.1 ;/*cm³/mol*/
  double Vav_lco2  = 32.6 ; /*cm³/mol*/
  double ak_h2o    = -2.209 ;
  double bk_h2o    = 3.097e-2 ;
  double ck_h2o    = -1.098e-4 ;
  double dk_h2o    = 2.048e-7 ;
  double ak_gco2   = 1.189 ;
  double bk_gco2   = 1.304e-2 ;
  double ck_gco2   = -5.446e-5 ;
  double dk_gco2   = 0. ;
  double ak_lco2   = 1.169 ;
  double bk_lco2   = 1.368e-2 ;
  double ck_lco2   = -5.380e-5 ;
  double dk_lco2   = 0. ;
  double P0        = 1. ;
  double P         = p_co2 ;
  double p_tot     = p_co2 ;
  
  double logK0_gco2   = ak_gco2 + bk_gco2*(T-273.0) + ck_gco2*(T-273.0)*(T-273.0) + dk_gco2*(T-273.0)*(T-273.0)*(T-273.0) ;
  double logK0_lco2   = ak_lco2 + bk_lco2*(T-273.0) + ck_lco2*(T-273.0)*(T-273.0) + dk_lco2*(T-273.0)*(T-273.0)*(T-273.0) ;
  double logK0_h2o    = ak_h2o + bk_h2o*(T-273.0) + ck_h2o*(T-273.0)*(T-273.0) + dk_h2o*(T-273.0)*(T-273.0)*(T-273.0) ;
  double K0_gco2   = pow(10,logK0_gco2) ;
  double K0_h2o    = pow(10,logK0_h2o) ;
  double lco2      = log(V/(V-B_CO2)) + B_CO2/(V-B_CO2) - 2*(A_CO2)*log((V + B_CO2)/V)/(R*pow(T,1.5)*B_CO2) \
			         + (A_CO2*B_CO2/(R*pow(T,1.5)*B_CO2*B_CO2))*(log((V + B_CO2)/V) - (B_CO2/(V+B_CO2))) - log(p_co2*V/(R*T)) ;
  double lh2o      = log(V/(V-B_CO2)) + B_H2O/(V-B_CO2) - 2*(A_H2OCO2)*log((V + B_CO2)/V)/(R*pow(T,1.5)*B_CO2) \
			         + (A_CO2*B_H2O/(R*pow(T,1.5)*B_CO2*B_CO2))*(log((V + B_CO2)/V) - (B_CO2/(V+B_CO2))) - log(p_co2*V/(R*T)) ;
  double fi_co2    = exp(lco2);
  double fi_h2o    = exp(lh2o);
  double A         = K0_h2o*exp((p_h2ol-P0)*Vav_h2o/(R*T))/(p_tot*fi_h2o) ;
  double B         = p_tot*fi_co2*exp(-(p_h2ol-P0)*Vav_gco2/(R*T))/(55.508*K0_gco2) ;
  double y_h2o     = MIN((1.-B)/(1./A-B),0.1) ;
  double x_co2     = B*(1.-y_h2o) ;
  return(x_co2) ;
}

double producty_h2o(double pascal_co2,double p_l,double T) 
/*co2 mole fraction in h2o rich phase,after Nicolas Spycher 2003*/
{
  double p_co2     = pascal_co2/1.e5;   /*unit of pressure:bar*/
  double p_h2ol	   = MAX(p_l/1.e5,1.) ;
  double R         = 83.14472 ;   /*cm3*bar*K-1*mol-1*/
  double A_CO2     = 7.54e7 - 4.13e4*T ;/*bar*cm⁶*K0.5*mol^-2*/
  double B_CO2     = 27.80 ; /*cm³/mol*/
  double B_H2O     = 18.18 ; /*cm³/mol*/
  double A_H2OCO2  = 7.89e7 ;   /*bar*cm⁶*K0.5*mol^-2*/
  double av        = 1.0 ;
  double bv        = -(R*T/p_co2) ;
  double cv        = -(R*T*B_CO2/p_co2 - A_CO2/(p_co2*pow(T,0.5)) + B_CO2*B_CO2) ;
  double dv		   = -(A_CO2*B_CO2/(p_co2*pow(T,0.5)));
  double V         = cubic_equation(av,bv,cv,dv) ;           /*cm³/mol*/
  double Vav_gco2  = 32.6 ; /*cm³/mol*/
  double Vav_h2o   = 18.1 ;/*cm³/mol*/
  double Vav_lco2  = 32.6 ; /*cm³/mol*/
  double ak_h2o    = -2.209 ;
  double bk_h2o    = 3.097e-2 ;
  double ck_h2o    = -1.098e-4 ;
  double dk_h2o    = 2.048e-7 ;
  double ak_gco2   = 1.189 ;
  double bk_gco2   = 1.304e-2 ;
  double ck_gco2   = -5.446e-5 ;
  double dk_gco2   = 0. ;
  double ak_lco2   = 1.169 ;
  double bk_lco2   = 1.368e-2 ;
  double ck_lco2   = -5.380e-5 ;
  double dk_lco2   = 0. ;
  double P0        = 1. ;
  double P         = p_co2 ;
  double p_tot     = p_co2 ;
  
  double logK0_gco2   = ak_gco2 + bk_gco2*(T-273.0) + ck_gco2*(T-273.0)*(T-273.0) + dk_gco2*(T-273.0)*(T-273.0)*(T-273.0) ;
  double logK0_lco2   = ak_lco2 + bk_lco2*(T-273.0) + ck_lco2*(T-273.0)*(T-273.0) + dk_lco2*(T-273.0)*(T-273.0)*(T-273.0) ;
  double logK0_h2o    = ak_h2o + bk_h2o*(T-273.0) + ck_h2o*(T-273.0)*(T-273.0) + dk_h2o*(T-273.0)*(T-273.0)*(T-273.0) ;
  double K0_gco2   = pow(10,logK0_gco2) ;
  double K0_h2o    = pow(10,logK0_h2o) ;
  double lco2      = log(V/(V-B_CO2)) + B_CO2/(V-B_CO2) - 2*(A_CO2)*log((V + B_CO2)/V)/(R*pow(T,1.5)*B_CO2) \
			         + (A_CO2*B_CO2/(R*pow(T,1.5)*B_CO2*B_CO2))*(log((V + B_CO2)/V) - (B_CO2/(V+B_CO2))) - log(p_co2*V/(R*T)) ;
  double lh2o      = log(V/(V-B_CO2)) + B_H2O/(V-B_CO2) - 2*(A_H2OCO2)*log((V + B_CO2)/V)/(R*pow(T,1.5)*B_CO2) \
			         + (A_CO2*B_H2O/(R*pow(T,1.5)*B_CO2*B_CO2))*(log((V + B_CO2)/V) - (B_CO2/(V+B_CO2))) - log(p_co2*V/(R*T)) ;
  double fi_co2    = exp(lco2);
  double fi_h2o    = exp(lh2o);
  double A         = K0_h2o*exp((p_h2ol-P0)*Vav_h2o/(R*T))/(p_tot*fi_h2o) ;
  double B         = p_tot*fi_co2*exp(-(p_h2ol-P0)*Vav_gco2/(R*T))/(55.508*K0_gco2) ;
  double y_h2o     = MIN((1.-B)/(1./A-B),0.1) ;
  double x_co2     = B*(1.-y_h2o) ;
  
  return(y_h2o) ;
}

double product_co2(double p_co2,double p_l,double T)    /*co2 concentration in water*/
{
  double c_co2 ;
  double x_co2 	= productx_co2(p_co2,p_l,T) ;
  double p0		= 1.e5 ;
  double c_co20 = 55.508*productx_co2((p0+53000),p_l,T) ;/*fit for p_co2=1atm*/
  double p1		= 4.e5 ;
  double x_co21 = productx_co2(p1,p_l,T) ;
  if (p_co2 < p0) {c_co2 = c_co20*p_co2/p0 ; }
  else if (p_co2 < p1) {c_co2 = c_co20 + (p_co2-p0)*(55.508*x_co21-c_co20)/(p1-p0) ; }
  else {c_co2 = 55.508*x_co2 ; }
  return(c_co2) ;
}

double dproductco2_co2(double p_co2,double p_l,double T)
{
  double dp_co2 = 1.e-4*p_co2 ;
  double p_co22 = p_co2 + dp_co2 ;
  double c_co22 = product_co2(p_co22,p_l,T) ;
  double c_co2  = product_co2(p_co2,p_l,T) ;
  double dc_co2 = (c_co22-c_co2)/dp_co2 ;
  return(dc_co2) ;
}

double dproductco2_h2o(double p_co2,double p_l,double T)
{
  double dp_l   = 1.e-4*p_l ;
  double p_l2 	= p_l + dp_l ;
  double c_co22 = product_co2(p_co2,p_l2,T) ;
  double c_co2  = product_co2(p_co2,p_l,T) ;
  double dc_co2 = (c_co22-c_co2)/dp_l ;
  return(dc_co2) ;
}

double product_gco2(double p_co2,double T) /*redilich-kwong equation(Spycher2003)*/
{
  double p         = p_co2/1.e5 ; /*unit of p: bar*/
  double p0        = 1. ;
  double R         = 83.14472 ;   /*cm3*bar*K-1*mol-1*/
  double A_CO2     = 7.54e7 - 4.13e4*T ;/*bar*cm⁶*K0.5*mol^-2*/
  double B_CO2     = 27.80 ; /*cm³/mol*/
  double a         = 1.0 ;
  double b         = -(R*T/p) ;
  double c         = -(R*T*B_CO2/p - A_CO2/(p*pow(T,0.5)) + B_CO2*B_CO2) ;
  double d         = -(A_CO2*B_CO2/(p*pow(T,0.5)));
  double a0         = 1.0 ;
  double b0         = -(R*T/p0) ;
  double c0         = -(R*T*B_CO2/p0 - A_CO2/(p0*pow(T,0.5)) + B_CO2*B_CO2) ;
  double d0         = -(A_CO2*B_CO2/(p0*pow(T,0.5)));
  double V         = cubic_equation(a,b,c,d) ;          /*cm³/mol*/
  double V0        = cubic_equation(a0,b0,c0,d0) ;
  double gc_co2    = (p_co2 > 1.e5) ? 1.e3/V : p*1.e3/V0  ;                             /*mol/dm³*/
  return(gc_co2) ;
}

double dproduct_gco2(double p_co2,double T)
{
  double p_co20 = 1.e5 ;
  double gc_co20 = product_gco2(p_co20,T) ;
  double dgc_co20 = gc_co20/p_co20 ;
  double gc_co2 = product_gco2(p_co2,T) ;
  double dp_co2 = 1.e-4*p_co2 ;
  double p_co22 = p_co2 + dp_co2 ;
  double gc_co22 = product_gco2(p_co22,T) ;
  double dgc_co2 = (p_co2 > 1.e5) ? (gc_co22-gc_co2)/dp_co2 : dgc_co20 ; 
  
  return(dgc_co2) ;
}

double product_gh2o(double p_co2,double p_l,double T)    /*co2 concentration in water*/
{
  double y_h2o 	= producty_h2o(p_co2,p_l,T) ;
  double g_co2  = product_gco2(p_co2,T) ;
  double gc_h2o 	= y_h2o*g_co2 ;
  return(gc_h2o) ;
}

double dproductgh2o_co2(double p_co2,double p_l,double T)   
{
  double gc_h2o 	= product_gh2o(p_co2,p_l,T) ;
  double dp_co2 = 1.e-4*p_co2 ;
  double p_co22 = p_co2 + dp_co2 ;
  double gc_h2o2 = product_gh2o(p_co22,p_l,T) ;
  double dgc_h2o = (gc_h2o2-gc_h2o)/dp_co2 ;   
  return(dgc_h2o) ;
}

double dproductgh2o_h2o(double p_co2,double p_l,double T)    
{
  double gc_h2o 	= product_gh2o(p_co2,p_l,T) ;
  double dp_l = 1.e-4*p_l ;
  double p_l2 = p_l + dp_l ;
  double gc_h2o2 = product_gh2o(p_co2,p_l2,T) ;
  double dgc_h2o = (gc_h2o2-gc_h2o)/dp_l ;   
  return(dgc_h2o) ;
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



