/* General features of the model:
 * Curves for CSH:
 *   - C/S ratio
 *   - H/S ratio
 *   - Molar Volume
 * Alkalis (sodium et potassium) : Na+, K+, NaOH, KOH, NaHCO3, NaCO3-
 * Dissolution kinetics for CH based on 
 * spherical crystal coated by a calcite layer.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* The Finite Volume Method */
#include "FVM.h"

#define TITLE   "Carbonation Of CBM with kinetics on C-S-H (2013)"
#define AUTHORS "Morandeau-Thiery-Dangla"

#include "PredefinedMethods.h"

/* Macros */

#define NEQ    	(8)
#define NVE    	(58)
#define NVI     (32)
#define NV0     (2)

#define E_C       (0)
#define E_q       (1)
#define E_mass    (2)
#define E_Ca      (3)
#define E_el      (4)
#define E_Na      (5)
#define E_K       (6)
#define E_Si      (7)

#define I_C_CO2   (0)
#define I_C_OH    (4)
#define I_P_L     (2)
#define I_N_CC    (3)
#define I_PSI     (1)
#define I_C_Na    (5)
#define I_C_K     (6)
#define I_N_Si_S  (7)

#define NOLOG_U   1
#define LOG_U     2
#define Ln10      2.302585093
#define U_CO2     LOG_U

#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])

#if (U_CO2 == LOG_U)
#define C_CO2(n)	(exp(Ln10*UNKNOWN(n,I_C_CO2)))
#define C_CO2n(n)	(exp(Ln10*UNKNOWNn(n,I_C_CO2)))
#else
#define C_CO2(n)	(UNKNOWN(n,I_C_CO2))
#define C_CO2n(n)	(UNKNOWNn(n,I_C_CO2))
#endif

#define N_Si_S(n)   (UNKNOWN(n,I_N_Si_S))
#define C_OH(n)     (UNKNOWN(n,I_C_OH))
#define N_CC(n)     (UNKNOWN(n,I_N_CC))
#define P_L(n)      (UNKNOWN(n,I_P_L))
#define PSI(n)      (UNKNOWN(n,I_PSI))
#define C_Na(n)     (UNKNOWN(n,I_C_Na))
#define C_K(n)      (UNKNOWN(n,I_C_K))

#define N_Si_Sn(n)   (UNKNOWNn(n,I_N_Si_S))
#define C_OHn(n)     (UNKNOWNn(n,I_C_OH))
#define N_CCn(n)     (UNKNOWNn(n,I_N_CC))
#define P_Ln(n)      (UNKNOWNn(n,I_P_L))
#define PSIn(n)      (UNKNOWNn(n,I_PSI))
#define C_Nan(n)     (UNKNOWNn(n,I_C_Na))
#define C_Kn(n)      (UNKNOWNn(n,I_C_K))


#define N_C(n)        (f[(n)])
#define N_q(n)        (f[(2+n)])
#define Mass(n)       (f[(4+n)])
#define N_Ca(n)       (f[(6+n)])
#define N_Na(n)       (f[(10+n)])
#define N_K(n)        (f[(12+n)])
#define N_Si(n)       (f[(16+n)])
#define W_C           (f[20])
#define W_q           (f[21])
#define W_m           (f[22])
#define W_Ca          (f[23])
#define W_Na          (f[25])
#define W_K           (f[26])
#define W_Si          (f[27])
#define N_CH(n)       (f[(28+n)])
#define S_CH_EQ(n)    (f[(30+n)])

#define N_Cn(n)       (f_n[(n)])
#define N_qn(n)       (f_n[(2+n)])
#define Mass_n(n)     (f_n[(4+n)])
#define N_Can(n)      (f_n[(6+n)])
#define N_Nan(n)      (f_n[(10+n)])
#define N_Kn(n)       (f_n[(12+n)])
#define N_Sin(n)      (f_n[(16+n)])
#define N_CHn(n)      (f_n[(28+n)])
#define S_CH_EQn(n)   (f_n[(30+n)])

#define KD_Ca           (va[(0)])
#define KD_OH           (va[(1)])
#define KD_H            (va[(2)])
#define KD_H2CO3        (va[(3)])
#define KD_HCO3         (va[(4)])
#define KD_CO3          (va[(5)])
#define KD_Na           (va[(6)])
#define KD_NaOH         (va[(7)])
#define KD_NaHCO3       (va[(8)])
#define KD_NaCO3        (va[(9)])
#define KD_m            (va[(10)])
#define KD_K            (va[(11)])
#define KD_KOH          (va[(12)])
#define KD_CaOH         (va[(13)])
#define KD_CaHCO3       (va[(14)])
#define KD_CaCO3aq      (va[(15)])
#define KD_CaOH2aq      (va[(16)])
#define KD_H3SiO4       (va[(17)])
#define KD_H2SiO4       (va[(18)])
#define KD_H4SiO4       (va[(19)])
#define KD_CaH2SiO4     (va[(20)])
#define KD_CaH3SiO4     (va[(21)])

#define KF_CO2     	    (va[(22)])
#define KF_Ca      		  (va[(23)])
#define KF_OH      		  (va[(24)])
#define KF_H       	  	(va[(25)])
#define KF_H2CO3   	    (va[(26)])
#define KF_HCO3    	    (va[(27)])
#define KF_CO3     	    (va[(28)])
#define KF_Na      		  (va[(29)])
#define KF_NaOH    	    (va[(30)])
#define KF_NaHCO3     	(va[(31)])
#define KF_NaCO3      	(va[(32)])
#define KF_K       	  	(va[(33)])
#define KF_KOH     	    (va[(34)])
#define KF_CaOH    	    (va[(35)])
#define KF_CaHCO3  	    (va[(36)])
#define KF_CaCO3aq 	    (va[(37)])
#define KF_CaOH2aq 	    (va[(38)])
#define KF_H3SiO4  	    (va[(39)])
#define KF_H2SiO4  	    (va[(40)])
#define KF_H4SiO4  	    (va[(41)])
#define KF_CaH2SiO4 	  (va[(42)])
#define KF_CaH3SiO4 	  (va[(43)])

#define Kpsi_Ca     		(va[(44)])
#define Kpsi_OH     		(va[(45)])
#define Kpsi_H      		(va[(46)])
#define Kpsi_HCO3   		(va[(47)])
#define Kpsi_CO3    		(va[(48)])
#define Kpsi_Na     		(va[(49)])
#define Kpsi_NaCO3  		(va[(50)])
#define Kpsi_q      		(va[(51)])
#define Kpsi_K      		(va[(52)])
#define Kpsi_CaOH   		(va[(53)])
#define Kpsi_CaHCO3 		(va[(54)])
#define Kpsi_H3SiO4 		(va[(55)])
#define Kpsi_H2SiO4 		(va[(56)])
#define Kpsi_CaH3SiO4	  (va[(57)])

#define V_S0(n)         (v0[(0+n)])


/* les valences */
#define z_ca    	 (2.)
#define z_h     	 (1.)
#define z_oh    	 (-1.)
#define z_hco3  	 (-1.)
#define z_co3   	 (-2.)
#define z_na    	 (1.)
#define z_naco3 	 (-1.)
#define z_k      	 (1.)
#define z_caoh   	 (1.)
#define z_cahco3 	 (1.)
#define z_h3sio4 	 (-1.)
#define z_h2sio4 	 (-2.)
#define z_cah3sio4 (1.)


/* Stoechimetrie des C-S-H */
#define x_jen 		  (1.5)
#define y_jen 		  (0.9)
#define z_jen 		  (0.9)
#define zz_jen 		  (1.83)
#define x_tobII 		(1.5)
#define y_tobII 		(1.8)
#define z_tobII 		(1.2)
#define zz_tobII		(2.7)
#define x_tobI 		  (2)
#define y_tobI 		  (2.4)
#define z_tobI 		  (2)
#define zz_tobI 		(3.6)
#define x_sil			  (0.)
#define y_sil			  (1.)
#define z_sil			  (2.)


/* volumes molaires partiels des ions (dm3/mole) */
#define v_h       		(-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       		(23.50e-3)    /* (23.50e-3)  d'apres TQN */
#define v_h2o      		(18.e-3)
#define v_h2co3    	  (50.e-3)
#define v_hco3     		(50.e-3)
#define v_co3      		(-2.3e-3)     /* (-2.3e-3) d'apres 3w */
#define v_ca       		(-18.7e-3)    /* (-18.7e-3 d'apres */
#define v_na       		(22.4e-3)
#define v_naoh     		(22.4e-3)
#define v_nahco3  	  (22.4e-3) 
#define v_naco3    	  (22.4e-3) 
#define v_k        		(43.93e-3)
#define v_koh      		(27.44e-3)
#define v_caoh     		(26.20e-3)  /*a modifier*/
#define v_cahco3   	  (26.20e-3)  /*a modifier*/
#define v_caco3aq  	  (26.20e-3)  /*a modifier*/
#define v_caoh2aq  	  (26.20e-3)  /*a modifier*/
#define v_cah2sio4 	  (26.20e-3)  /*a modifier*/
#define v_cah3sio4 	  (26.20e-3)  /*a modifier*/
#define v_h3sio4   	  (50.e-3)    /*a modifier*/
#define v_h2sio4   	  (50.e-3)    /*a modifier*/
#define v_h4sio4   	  (50.e-3)    /*a modifier*/

/* volumes molaires solides (dm3/mole) */
#define V_CH	    	  (33.e-3)
#define V_CC	      	(37.e-3)

/* Masses molaires (unite arbitraire = M_H) */
#define M_Ca    		 (40.1)
#define M_H2CO3  	   (62.)
#define M_HCO3   		 (61.)
#define M_CO3    		 (60.)
#define M_OH     		 (17.)
#define M_H      		 (1.)
#define M_H2O    		 (18.)
#define M_Na	 		   (23.)
#define M_NaOH	  	 (40.)
#define M_NaHCO3  	 (84.)
#define M_NaCO3	  	 (83.)
#define M_CO2     	 (44.)
#define M_CaOH2   	 (74.)
#define M_CaCO3   	 (100.)
#define M_CaCO3CSH 	 (100.)
#define M_CaO 		   (56.)
#define M_SiO2 		   (60.)
#define M_K      	 	 (39.)  
#define M_KOH    		 (56.)
#define M_CaOH    	 (57.)
#define M_CaHCO3  	 (101.)
#define M_CaCO3aq 	 (100.)
#define M_CaOH2aq  	 (74.)
#define M_CaH2SiO4 	 (134.)
#define M_CaH3SiO4   (135.)
#define M_H3SiO4  	 (95)
#define M_H2SiO4  	 (94)
#define M_H4SiO4  	 (96)

/* coefficients de diffusion moleculaire (dm2/s) */
#define d_oh       		(5.273e-7)    /* (5.273e-7) d'apres TQN */
#define d_h       		(9.310e-7)
#define d_ca       		(7.92e-8)
#define d_h2co3    	  (7.2e-8)
#define d_hco3     		(11.8e-8)
#define d_co3      		(9.55e-8)
#define d_na       		(1.33e-7)
#define d_naoh     		(1.33e-7)
#define d_nahco3  	  (1.33e-7)
#define d_naco3    	  (1.33e-7)
#define d_k        		(1.957e-7)
#define d_koh      		(1.957e-7)
#define d_caoh     		(7.92e-8) /*a modifier */
#define d_cahco3  	  (7.92e-8) /*a modifier */
#define d_caco3aq 	  (7.92e-8) /*a modifier */
#define d_caoh2aq  	  (7.92e-8) /*a modifier */
#define d_cah2sio4   	(7.92e-8) /*a modifier */
#define d_cah3sio4   	(7.92e-8) /*a modifier */
#define d_h3sio4   	  (11.8e-8) /*a modifier */
#define d_h2sio4   	  (11.8e-8) /*a modifier */
#define d_h4sio4     	(11.8e-8) /*a modifier */

#define d_co2      		(1.6e-3)

/* constantes d'equilibre (ref = 1 mole/L) */
#define k_e       		(1.e-14)                /* autoprotolyse de l'eau */
#define k_h        		(1.)                    /* cste de Henry */
#define k_co3      		(4.570881896148751e3)   /* Equilibre de HCO3  <-> CO3 */
#define k_ca       		(3.890451449942805e-9)  /* Equilibre de CaCO3 */
#define k_1        		(2.187761623949552e-8)  /* Equilibre de H2CO3 <-> HCO3 */
#define k_2        		(6.456542290346550e-6)  /* Equilibre de Ca(OH)2 */
#define k_naoh     		(6.60693448e-15)        /* Equilibre de Na <-> NaOH */
#define k_nahco3     	(0.5623413252)          /* Equilibre de Na <-> NaHCO3  */
#define k_naco3    	  (8.72971368e-10)        /* Equilibre de Na <-> NaCO3 */
#define k_caoh     		(1.65958691e-13)        /* Equilibre de Ca <-> CaOH */
#define k_cahco3   	  (12.76438809)           /* Equilibre de Ca <-> CaHCO3 */
#define k_caco3      	(7.852356346e-8)        /* Equilibre de Ca <-> CaCO3aqueux */
#define k_koh      		(3.4673685e-15)         /* Equilibre de K + H20 <-> KOH + H */
#define k_caoh2    	  (1.)                    /* Equilibre de Ca(OH)2 en solution */
#define k_h2sio4     	(2.13796209e+13)        /* Equilibre de H2SiO4-- + H+ <-> H3SiO4- */
#define k_h4sio4     	(1.380384265e+23)       /* Equilibre de H2SiO4 -- + 2H+ <-> H4SiO4 */
#define k_jen      		(1.412537545e-12)       /* Solubilité des CSH-jennite */
#define k_tobII    		(2.238721139e-14)       /* Solubilité des CSH-tobermoriteII */
#define k_tobI     		(6.30957344e-19)        /* Solubilité des CSH-tobermoriteII */
#define k_sil      		(1.936421964e-3)        /* Solubilité de la silice amorphe */
#define k_cah2sio4 	  (39810.71)              /* Ca <-> CaH2SiO4*/
#define k_cah3sio4 	  (15.84)                 /* Ca <-> CaH3SiO4*/

/* viscosite (Pa.s) */
#define mu_l       	(1.e-3)

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (C/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (J/mole) */

/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp   	   	 (0.5)

/*Loi de Davies pour la prise en compte de l'activite ionique*/
#define A_DAVIES  (0.5)
#define B_DAVIES  (0.24)


/* Material Properties */
#define SaturationDegree(p)              (Curve_ComputeValue(Element_GetCurve(el),p))
#define RelativePermeabilityToLiquid(p)  (Curve_ComputeValue(Element_GetCurve(el) + 1,p))


/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(s)      (Curve_ComputeValue(Element_GetCurve(el) + 2,s))
#define dCalciumSiliconRatioInCSH(s)     (Curve_ComputeDerivative(Element_GetCurve(el) + 2,s))
#define WaterSiliconRatioInCSH(s)        (Curve_ComputeValue(Element_GetCurve(el) + 3,s))
#define MolarVolumeOfCSH(s)              (Curve_ComputeValue(Element_GetCurve(el) + 4,s))


/* CH Properties */
#define C_CO2_eq                         (k_ca*k_1/(k_h*k_co3*k_2))
#define SaturationDegreeOfCH(c_co2)      (C_CO2_eq/c_co2)
#define SaturationDegreeOfCHAtEq         ComputeS_CH_EQ


/* S-H Properties */
/* Saturation degree of dissolved S-H */
#define S_SHeq(s)                        (Curve_ComputeValue(Element_GetCurve(el) + 5,s))
#define LNS_SHEQ(s,s_ch)                 (log(S_SHeq(s)) + CalciumSiliconRatioInCSH(s)*(log(s) - log(s_ch)))
#define SaturationDegreeOfSH(s,s_ch)     (S_SHeq(s)*pow(s/s_ch,CalciumSiliconRatioInCSH(s)))
/* #define SaturationDegreeOfSH(s,s_ch)     exp(LNS_SHEQ(s,s_ch)) */
#define IonActivityProductOfSH(s,s_ch)   (k_sil*SaturationDegreeOfSH(s,s_ch))


/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(char *s) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    concentration1(double x_co2, double x_na_tot, double x_k_tot, double *pointeur_x_oh, double *pointeur_x_h3sio4, double *pointeur_x_na, double *pointeur_x_k) ;

static double  dn1_caoh2sdt(double,double) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static double  ComputeS_CH_EQ(Element_t*,double,double,double) ;
static double  ComputeS_CH_EQ1(Element_t*,double,double,double) ;


/* Parametres */
static double phii,k_int,a_1,a_2,c_2,n_ch0,x_na0,x_k0,tau ;
static double p_g = 0. ;

#define NbOfComponents    (52)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_C_H          (9)
#define I_C_H2O        (10)

#define I_C_HCO3       (8)
#define I_C_H2CO3      (11)
#define I_C_CO3        (12)

#define I_C_Ca         (13)
#define I_C_CaOH       (14)
#define I_C_CaHCO3     (15)
#define I_C_CaCO3aq    (16)
#define I_C_CaOH2aq    (17)

#define I_C_H2SiO4     (18)
#define I_C_H3SiO4     (19)
#define I_C_H4SiO4     (20)

#define I_C_CaH2SiO4   (21)
#define I_C_CaH3SiO4   (22)

#define I_C_NaOH       (23)
#define I_C_NaHCO3     (24)
#define I_C_NaCO3      (25)

#define I_C_KOH        (26)

#define I_S_CH         (27)
#define I_RHO_L        (28)

#define I_N_C          (36)
#define I_N_Ca         (37)
#define I_N_Si         (38)
#define I_N_K          (39)
#define I_N_Na         (40)
#define I_Mass         (41)
#define I_N_Q          (42)

#define I_N_CH         (43)
#define I_V_S          (45)

#define I_N_CHn        (46)
#define I_V_S0         (48)

#define I_Phi          (49)

#define I_S_CH_EQ      (50)
#define I_S_CH_EQn     (51)

#define NbOfComponentFluxes    (7)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_C           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_Na          (3)
#define I_W_K           (4)
#define I_W_m           (5)
#define I_W_q           (6)


int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)     return (0) ;
  else if(strcmp(s,"k_int") == 0)   return (1) ;
  else if(strcmp(s,"N_CaOH2") == 0) return (2) ;
  else if(strcmp(s,"N_Si") == 0)    return (4) ;
  else if(strcmp(s,"X_K") == 0)     return (5) ;
  else if(strcmp(s,"X_Na") == 0)    return (6) ;
  else if(strcmp(s,"A_1") == 0)     return (7) ;
  else if(strcmp(s,"A_2") == 0)     return (8) ;
  else if(strcmp(s,"C_2") == 0)     return (9) ;
  else if(strcmp(s,"R_CaOH2") == 0) return (10) ;
  else if(strcmp(s,"D") == 0) 	    return (11) ;
  else if(strcmp(s,"Tau") == 0) 	  return (12) ;
  else return(-1) ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_C,"carbone") ;
  Model_CopyNameOfEquation(model,E_q,"charge") ;
  Model_CopyNameOfEquation(model,E_mass,"masse") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_el,"E_el") ;
  Model_CopyNameOfEquation(model,E_Na,"sodium") ;
  Model_CopyNameOfEquation(model,E_K,"potassium") ;
  Model_CopyNameOfEquation(model,E_Si,"silice") ;
  
  
#if (U_CO2 == LOG_U)
  Model_CopyNameOfUnknown(model,I_C_CO2,"logc_co2") ;
#else
  Model_CopyNameOfUnknown(model,I_C_CO2,"c_co2") ;
#endif
  Model_CopyNameOfUnknown(model,I_N_Si_S,"n_si_s") ;
  Model_CopyNameOfUnknown(model,I_C_OH,"c_oh") ;
  Model_CopyNameOfUnknown(model,I_P_L,"p_l") ;
  Model_CopyNameOfUnknown(model,I_N_CC,"c_caco3") ;
  Model_CopyNameOfUnknown(model,I_PSI,"psi") ;
  Model_CopyNameOfUnknown(model,I_C_Na,"c_na") ;
  Model_CopyNameOfUnknown(model,I_C_K,"c_k") ;
  
  return(0) ;
}

int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 13 ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Initialisation automatique */
  {
    double h   = 5.6e-6 ;  /* (moles/dm2/s) */
    double R_0 = Material_GetProperty(mat)[pm("R_CaOH2")] ;
    double D   = Material_GetProperty(mat)[pm("D")] ;
    n_ch0 = Material_GetProperty(mat)[pm("N_CaOH2")] ; /* contenu molaire initial en CaOH2 */
    a_2 = 3*h/R_0*n_ch0*V_CH ; /* (dm/mole/s) these MT p 227 */
    c_2 = h*R_0/D ;     /* (1/dm) these MT p 228 */
  }
  
  Material_GetProperty(mat)[pm("A_2")] = a_2 ;
  Material_GetProperty(mat)[pm("C_2")] = c_2 ;
  
  return(n_donnees) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 7 equations:\n\
\t- la conservation de la masse de C      (c_co2)\n\
\t- la conservation de la charge          (psi)\n\
\t- la conservation de la masse totale    (p_l)\n\
\t- la conservation de la masse de Ca     (c_caco3)\n\
\t- 1 equation de cinetique               (c_hco3)\n\
\t- Electroneutralite                     (c_oh)\n\
\t- la conservation de la masse de Na     (c_na)\n\
\t- la conservation de la masse de K      (c_k)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n\
\t pression : Pa !\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"k_int = 5e-19     # Permeabilite intrinseque (dm2)\n") ;
  fprintf(ficd,"N_CaOH2 = 6.1     # Contenu initial en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"N_Si = 2.4       # contenu initial en CSH (moles/L)\n") ;
  fprintf(ficd,"T_csh = 1.8e-4    # k_h/T_csh = temps caracteristique de carbo des CSH (1/s)\n") ;
  fprintf(ficd,"A_1 = 150         # Coef de la cinetique 1 (dm/mole/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Coef de la cinetique 2 (dm/mole/s)\n") ;
  fprintf(ficd,"C_2 = 0.14e6      # Coef de la cinetique 2 (1/dm)\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;  

  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = (Element_IsSubmanifold(el)) ? 0 : NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double *f = Element_GetImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  x_na0    = GetProperty("X_Na") ;
  x_k0     = GetProperty("X_K") ;
  tau      = GetProperty("Tau") ;
  
  {
    double x_na_tot 	= x_na0 ;
    double x_k_tot  	= x_k0 ;
    double x_co2   	  = C_CO2_eq ;
    double x_na,x_k,x_h3sio4,x_oh ;
    
    concentration1(x_co2,x_na_tot,x_k_tot,&x_oh,&x_h3sio4,&x_na,&x_k) ;

    for(i = 0 ; i < nn ; i++) {
      double s_ch       = SaturationDegreeOfCH(x_co2) ;
      double s_cheq     = s_ch ;
      double n_cc       = N_CC(i) ;
      double n_si_s     = N_Si_S(i) ;
      double v_csh      = MolarVolumeOfCSH(s_cheq) ;
      double v_s0       = V_CH*n_ch0 + V_CC*n_cc + v_csh*n_si_s ;
      
#if (U_CO2 == LOG_U)
      UNKNOWN(i,I_C_CO2) = log(x_co2)/Ln10 ;
#else
      C_CO2(i)    = x_co2 ;
#endif
      C_Na(i)     = x_na ;
      C_K(i)      = x_k ;
      C_OH(i)     = x_oh ;
      
      /* Solid contents */
      V_S0(i)    = v_s0 ;
      N_CH(i)    = n_ch0 ;
      /* Liquid components */
      S_CH_EQ(i) = s_cheq ;
    }
  }
  

  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x       = ComputeComponents(el,u,f,0,i) ;
    
    /* Back up */
    N_C(i)  = x[I_N_C] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Na(i) = x[I_N_Na] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ; 
    Mass(i) = x[I_Mass] ;
    N_q(i)  = x[I_N_Q] ;
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el,u,f) ;

  /* Flux */
  ComputeFluxes(el,u) ;
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double *f = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetPreviousNodalUnknown(el,i) ;
  }
  
  /*
    Coefficients de transfert
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetCurrentNodalUnknown(el,i) ;
  }
  
  
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  tau      = GetProperty("Tau") ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x      = ComputeComponents(el,u,f_n,dt,i) ;
    
    /* Back up */
    N_C(i)  = x[I_N_C] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Na(i) = x[I_N_Na] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ; 
    Mass(i) = x[I_Mass] ;
    N_q(i)  = x[I_N_Q] ;

    /* contenus solides */
    N_CH(i) = x[I_N_CH] ;
    /* Liquid components */
    S_CH_EQ(i) = x[I_S_CH_EQ] ;

    {
      double x_co2      = x[I_C_CO2] ;
      double x_oh    	  = x[I_C_OH] ;
      double x_na    	  = x[I_C_Na] ;
      double x_k     	  = x[I_C_K] ;
      double x_ca    	  = x[I_C_Ca] ;
      double x_h2o      = x[I_C_H2O] ;
      double n_si_s     = x[I_N_Si_S] ;
      double s_ch       = x[I_S_CH] ;
      double s_cheq     = x[I_S_CH_EQ] ;
      double s_sh       = SaturationDegreeOfSH(s_cheq,s_ch) ;
      if(x_co2 < 0 || x_oh <= 0 || x_h2o <= 0 || x_na < 0 || x_k < 0 || x_ca < 0 || n_si_s < 0.) {
        double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        double n_cc       = x[I_N_CC] ;
        double x_naoh    	= x[I_C_NaOH] ;
        double x_nahco3  	= x[I_C_NaHCO3] ;
        double x_naco3 	  = x[I_C_NaCO3] ;
        printf("\n") ;
        printf("en x     = %e\n",x0) ;
        printf("x_co2    = %e\n",x_co2) ;
        printf("x_oh     = %e\n",x_oh) ;
        printf("x_h2o    = %e\n",x_h2o) ;
        printf("n_cc     = %e\n",n_cc) ;
        printf("x_na     = %e\n",x_na) ;
        printf("x_k      = %e\n",x_k) ;
        printf("x_naoh   = %e\n",x_naoh) ;
        printf("x_nahco3 = %e\n",x_nahco3) ;
        printf("x_naco3  = %e\n",x_naco3) ;
        return(1) ;
      }
      if(s_sh > 1 || s_ch < 0) {
        double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        printf("\n") ;
        printf("en x     = %e\n",x0) ;
        printf("s_sh     = %e\n",s_sh) ;
        arret("ComputeImplicitTerms") ;
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;
  /* Flux */
  ComputeFluxes(el,u) ;

  return(0) ;
}



int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  double *u[Element_MaxNbOfNodes] ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double c[4*NEQ*NEQ] ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data 
  */
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  tau      = GetProperty("Tau") ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }

#if (U_CO2 == LOG_U)
  for(i=0;i<2*NEQ;i++){
    K(i,I_C_CO2)     *= Ln10*C_CO2(0) ;
    K(i,I_C_CO2+NEQ) *= Ln10*C_CO2(1) ;
  }
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double *surface = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = surface[1] ;
  }
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;
  /*
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  R(0,E_mass) -= volume[0]*(Mass(0) - Mass_n(0)) + dt*surf*W_m ;
  R(1,E_mass) -= volume[1]*(Mass(1) - Mass_n(1)) - dt*surf*W_m ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Electroneutralite : q = 0
  */
  R(0,E_el) -= volume[0]*N_q(0) ;
  R(1,E_el) -= volume[1]*N_q(1) ;
  /*
    Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*W_Na ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*W_Na ;
  /*
    Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ;

  /*
    Conservation de Si (silice) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    nso = 36 ;
  int    i ;
  double *u[MAX_NOEUDS] ;
#define UNKNOWN_s(i)   UNKNOWN(j,i)


  if(Element_IsSubmanifold(el)) return(0) ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
   
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  tau      = GetProperty("Tau") ;

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double xm = (x0 + x1)*0.5 ;
    int j = (Element_IsSubmanifold(el)) ? 0 : ((s[0] < xm) ? 0 : 1) ;
    /* pression */
    double p_l     =  UNKNOWN_s(I_P_L) ;
    /* saturation */
    double p_c     = p_g - p_l ;
    double s_l     = SaturationDegree(p_c) ;
    /* molarites */
    double *x = ComputeComponents(el,u,f,0,j) ;
    
    double x_co2      = x[I_C_CO2] ;
    double x_oh    	  = x[I_C_OH] ;
    double x_hco3  	  = x[I_C_HCO3] ;
    double x_na    	  = x[I_C_Na] ;
    double x_k     	  = x[I_C_K] ;
    double x_co3   	  = x[I_C_CO3] ;
    double x_h     	  = x[I_C_H] ;
    double x_ca    	  = x[I_C_Ca] ;
    double x_naoh    	= x[I_C_NaOH] ;
    double x_nahco3  	= x[I_C_NaHCO3] ;
    double x_naco3 	  = x[I_C_NaCO3] ;
    double x_koh   	  = x[I_C_KOH] ;
    double x_caoh  	  = x[I_C_CaOH] ;
    double x_cahco3  	= x[I_C_CaHCO3] ;
    double x_caco3aq 	= x[I_C_CaCO3aq] ;
    double x_caoh2aq 	= x[I_C_CaOH2aq] ;
    double x_h4sio4   = x[I_C_H4SiO4] ;
    double x_h3sio4   = x[I_C_H3SiO4] ;
    double x_h2sio4   = x[I_C_H2SiO4] ;
    double x_cah2sio4 = x[I_C_CaH2SiO4];
    double x_cah3sio4 = x[I_C_CaH3SiO4] ;
    
    /* Force Ionique */
    double I = 0.5*(z_ca*z_ca*x_ca + z_h2sio4*z_h2sio4*x_h2sio4 + z_co3*z_co3*x_co3   + z_cahco3*z_cahco3*x_cahco3 + z_caoh*z_caoh*x_caoh + z_k*z_k*x_k + z_na*z_na*x_na + z_h*z_h*x_h + z_h3sio4*z_h3sio4*x_h3sio4 + z_naco3*z_naco3*x_naco3 + z_hco3*z_hco3*x_hco3 + z_oh*z_oh*x_oh ) ;

    /* densite de charge */
    double x_q = x[I_N_Q] ;

    /* contenus solides */
    double n_si_s     = UNKNOWN_s(I_N_Si_S) ;
    double n_cc 	    = UNKNOWN_s(I_N_CC) ;
    double n_ch       = (N_CH(0) + N_CH(1))*0.5 ;
    double s_cheq     = x[I_S_CH_EQ] ;
    double x_csh      = CalciumSiliconRatioInCSH(s_cheq) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  
    double CsurS      = (x_csh + n_ch/n_si_s) ;
    

    /* porosite */
    double phi        = x[I_Phi] ;

    double ph = 14 + log(x_oh)/log(10.) ;
    double n_Na = 0.5*(N_Na(0) + N_Na(1)) ;
    double n_Ca = 0.5*(N_Ca(0) + N_Ca(1)) ;
    double n_Si = 0.5*(N_Si(0) + N_Si(1)) ; 


    /* Transferts */
    double dx        = x1 - x0 ;
    double grd_psi   = (PSI(1) - PSI(0))/dx ;

    /* quantites exploitees */
    i = 0 ;
    Result_Store(r + i++,&x_co2,"x_co2",1) ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&n_si_s,"n_Si_s",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&n_ch,"n_CH",1) ;
    Result_Store(r + i++,&x_ca,"x_ca",1) ;
    Result_Store(r + i++,&x_co3,"x_co3",1) ;
    Result_Store(r + i++,&x_hco3,"x_hco3",1) ;
    Result_Store(r + i++,&n_cc,"n_CC",1) ;
    Result_Store(r + i++,&x_h,"x_h",1) ;
    Result_Store(r + i++,&x_oh,"x_oh",1) ;
    Result_Store(r + i++,&s_l,"saturation",1) ;
    Result_Store(r + i++,&grd_psi,"grad_psi",1) ;
    Result_Store(r + i++,&x_q,"charge",1) ;
    Result_Store(r + i++,&x_na,"x_na",1) ;
    Result_Store(r + i++,&x_naoh,"x_naoh",1) ;
    Result_Store(r + i++,&x_nahco3,"x_nahco3",1) ;
    Result_Store(r + i++,&x_naco3,"x_naco3",1) ;
    Result_Store(r + i++,&x_k,"x_k",1) ;
    Result_Store(r + i++,&x_koh,"x_koh",1) ;
    Result_Store(r + i++,&x_caoh,"x_caoh",1) ;
    Result_Store(r + i++,&x_cahco3,"x_cahco3",1) ;
    Result_Store(r + i++,&x_caco3aq,"x_caco3aq",1) ;
    Result_Store(r + i++,&x_caoh2aq,"x_caoh2aq",1) ;
    Result_Store(r + i++,&p_l,"p_l",1) ;
    Result_Store(r + i++,&x_h3sio4,"x_h3sio4",1) ;
    Result_Store(r + i++,&n_Na,"n_Na",1) ;
    Result_Store(r + i++,&n_Ca,"n_Ca",1) ;
    Result_Store(r + i++,&n_Si,"n_Si",1) ;
    Result_Store(r + i++,&n_ca_s,"n_Ca_s",1) ;
    Result_Store(r + i++,&x_cah2sio4,"x_cah2sio4",1) ;
    Result_Store(r + i++,&x_cah3sio4,"x_cah3sio4",1) ;
    Result_Store(r + i++,&CsurS,"CsurS",1) ;
    Result_Store(r + i++,&x_h4sio4,"x_h4sio4",1) ;
    Result_Store(r + i++,&x_h2sio4,"x_h2sio4",1) ;
    Result_Store(r + i++,&I,"I",1) ;
  }
  
  
  if(i != nso) arret("so50") ;
  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int    i ; 
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  tau      = GetProperty("Tau") ;

  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;
    /* pressions */
    double p_l     = P_L(i) ;
    double p_c     = p_g - p_l ;
    /* saturation */
    double s_l     = SaturationDegree(p_c) ;
    
    double x_h    	  = x[I_C_H] ;
    double x_oh    	  = x[I_C_OH] ;
    double x_hco3  	  = x[I_C_HCO3] ;
    double x_na    	  = x[I_C_Na] ;
    double x_k     	  = x[I_C_K] ;
    double x_h2co3 	  = x[I_C_H2CO3] ;
    double x_co3   	  = x[I_C_CO3] ;
    double x_ca    	  = x[I_C_Ca] ;
    double x_naoh    	= x[I_C_NaOH] ;
    double x_nahco3  	= x[I_C_NaHCO3] ;
    double x_naco3 	  = x[I_C_NaCO3] ;
    double x_koh   	  = x[I_C_KOH] ;
    double x_caoh  	  = x[I_C_CaOH] ;
    double x_cahco3  	= x[I_C_CaHCO3] ;
    double x_caco3aq 	= x[I_C_CaCO3aq] ;
    double x_caoh2aq 	= x[I_C_CaOH2aq] ;
    double x_h4sio4   = x[I_C_H4SiO4] ;
    double x_h3sio4   = x[I_C_H3SiO4] ;
    double x_h2sio4   = x[I_C_H2SiO4] ;
    double x_cah2sio4 = x[I_C_CaH2SiO4];
    double x_cah3sio4 = x[I_C_CaH3SiO4] ;
    
    /* masse volumique liquide */
    double rho_l      = x[I_RHO_L] ;

    /* Porosity */
    double phi        = x[I_Phi] ;


    /* permeabilite */
    double k_l  = (k_int/mu_l)*RelativePermeabilityToLiquid(p_c)*pow(phi/phii,3.)*pow(((1-phii)/(1-phi)),2.) ;
    /* tortuosite gaz */
    /*double tau  = pow(phi,1/3)*pow(s_g,7/3) ;*/
    /* double tau  = pow(phi,1.74)*pow(s_g,3.20) ; */
    /* tortuosite liquide */
    double iff    = 0.00029*exp(9.95*phi)/(1+625*pow((1-s_l),4)) ;

    /* Humidité relative */
    double hr = exp(-p_c*M_H2O/(RT*rho_l)) ;
 
    KD_Ca     	+= x_ca*k_l ;
    KD_H2CO3  	+= x_h2co3*k_l ;
    KD_HCO3   	+= x_hco3*k_l ;
    KD_CO3    	+= x_co3*k_l ;
    KD_OH     	+= x_oh*k_l ;
    KD_H      	+= x_h*k_l ;
    KD_Na     	+= x_na*k_l ;
    KD_NaOH   	+= x_naoh*k_l ;
    KD_NaHCO3 	+= x_nahco3*k_l ;
    KD_NaCO3  	+= x_naco3*k_l ;
    KD_K      	+= x_k*k_l ;
    KD_KOH    	+= x_koh*k_l ;
    KD_CaOH   	+= x_caoh*k_l ;
    KD_CaHCO3 	+= x_cahco3*k_l ;
    KD_CaCO3aq	+= x_caco3aq*k_l ;
    KD_CaOH2aq	+= x_caoh2aq*k_l ;
    KD_H3SiO4 	+= x_h3sio4*k_l ;
    KD_H2SiO4 	+= x_h2sio4*k_l ;
    KD_H4SiO4 	+= x_h4sio4*k_l ;
    KD_CaH3SiO4 += x_cah3sio4*k_l ;
    KD_CaH2SiO4 += x_cah2sio4*k_l ;
    KD_m      	+= rho_l*k_l ;
    
    /* KF_CO2    += phi_g*tau*d_co2 ; */
    KF_CO2    	+= (1.6e-3)*pow(phi,1.8)*pow(1-hr,2.2) ;
    KF_Ca     	+= d_ca*iff ;
    KF_OH     	+= d_oh*iff ;
    KF_H      	+= d_h*iff ;
    KF_H2CO3  	+= d_h2co3*iff ;
    KF_HCO3   	+= d_hco3*iff ;
    KF_CO3    	+= d_co3*iff ;
    KF_Na     	+= d_na*iff ;
    KF_NaOH   	+= d_naoh*iff ;
    KF_NaHCO3 	+= d_nahco3*iff ;
    KF_NaCO3  	+= d_naco3*iff ;
    KF_K      	+= d_k*iff ;
    KF_KOH    	+= d_koh*iff ;
    KF_CaOH   	+= d_caoh*iff ;
    KF_CaHCO3 	+= d_cahco3*iff ;
    KF_CaCO3aq	+= d_caco3aq*iff ;
    KF_CaOH2aq	+= d_caoh2aq*iff ;
    KF_H3SiO4 	+= d_h3sio4*iff ;
    KF_H2SiO4 	+= d_h2sio4*iff ;
    KF_H4SiO4 	+= d_h4sio4*iff ;
    KF_CaH3SiO4 += d_cah3sio4*iff ;
    KF_CaH2SiO4 += d_cah2sio4*iff ;
    
    Kpsi_Ca   	 += FARADAY/RT*KF_Ca*z_ca*x_ca ;
    Kpsi_HCO3 	 += FARADAY/RT*KF_HCO3*z_hco3*x_hco3 ;
    Kpsi_CO3  	 += FARADAY/RT*KF_CO3*z_co3*x_co3 ;
    Kpsi_OH   	 += FARADAY/RT*KF_OH*z_oh*x_oh ;
    Kpsi_H    	 += FARADAY/RT*KF_H*z_h*x_h ;
    Kpsi_Na    	 += FARADAY/RT*KF_Na*z_na*x_na ;
    Kpsi_NaCO3 	 += FARADAY/RT*KF_NaCO3*z_naco3*x_naco3 ;
    Kpsi_K     	 += FARADAY/RT*KF_K*z_k*x_k ;
    Kpsi_CaOH 	 += FARADAY/RT*KF_CaOH*z_caoh*x_caoh ;
    Kpsi_CaHCO3	 +=FARADAY/RT*KF_CaHCO3*z_cahco3*x_cahco3 ;
    Kpsi_H3SiO4	 +=FARADAY/RT*KF_H3SiO4*z_h3sio4*x_h3sio4 ;
    Kpsi_H2SiO4	 +=FARADAY/RT*KF_H2SiO4*z_h2sio4*x_h2sio4 ;
    Kpsi_CaH3SiO4+=FARADAY/RT*KF_CaH3SiO4*z_cah3sio4*x_cah3sio4 ;
    Kpsi_q    	 += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_na*Kpsi_Na + z_naco3*Kpsi_NaCO3 + z_k*Kpsi_K + z_caoh*Kpsi_CaOH + z_cahco3*Kpsi_CaHCO3 + z_h3sio4*Kpsi_H3SiO4 + z_h2sio4*Kpsi_H2SiO4 + z_cah3sio4*Kpsi_CaH3SiO4 ;
  }
  
  /* moyenne */
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}


void ComputeFluxes(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *grd = dComponents ;

  /* Gradients (electric potential included) */
  {
    double *x1 = ComputeComponents(el,u,f,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,f,0.,0) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] -= x0[i] ;
  }
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] /= dx ;
  }
  
  /* Fluxes */
  {
    double *w = Fluxes(el,grd) ;

    W_C     = w[I_W_C] ;
    W_Ca    = w[I_W_Ca] ;
    W_Na    = w[I_W_Na] ;
    W_Si    = w[I_W_Si] ;
    W_q     = w[I_W_q] ;
    W_K     = w[I_W_K] ;
    W_m     = w[I_W_m] ;
  }
    
}


double* Fluxes(Element_t *el,double *grd)
{
  double *va = Element_GetExplicitTerm(el) ;
  double *w  = ComponentFluxes ;

  /* Gradients */
  double grd_p_l      = grd[I_P_L] ;
  double grd_co2      = grd[I_C_CO2] ;
  double grd_oh       = grd[I_C_OH] ;
  double grd_na       = grd[I_C_Na] ;
  double grd_k        = grd[I_C_K] ;
  double grd_ca       = grd[I_C_Ca] ;
  double grd_h        = grd[I_C_H] ;
  double grd_h2co3    = grd[I_C_H2CO3] ;
  double grd_hco3     = grd[I_C_HCO3] ;
  double grd_co3      = grd[I_C_CO3] ;
  double grd_naoh     = grd[I_C_NaOH] ;
  double grd_nahco3   = grd[I_C_NaHCO3] ;
  double grd_naco3    = grd[I_C_NaCO3] ;
  double grd_koh      = grd[I_C_KOH] ;
  double grd_caoh     = grd[I_C_CaOH] ;
  double grd_cahco3   = grd[I_C_CaHCO3] ;
  double grd_caco3aq  = grd[I_C_CaCO3aq] ;
  double grd_caoh2aq  = grd[I_C_CaOH2aq] ;
  double grd_h2sio4   = grd[I_C_H2SiO4] ;
  double grd_h3sio4   = grd[I_C_H3SiO4] ;
  double grd_h4sio4   = grd[I_C_H4SiO4] ;
  double grd_cah2sio4 = grd[I_C_CaH2SiO4] ;
  double grd_cah3sio4 = grd[I_C_CaH3SiO4] ;
  double grd_psi      = grd[I_PSI] ;

    /* Flux */
  double w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca       - Kpsi_Ca*grd_psi    ;
  double w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3   - Kpsi_HCO3*grd_psi  ;
  double w_h3sio4= - KD_H3SiO4*grd_p_l- KF_H3SiO4*grd_h3sio4 - Kpsi_H3SiO4*grd_psi ;
  double w_h2sio4= - KD_H2SiO4*grd_p_l- KF_H2SiO4*grd_h2sio4 - Kpsi_H2SiO4*grd_psi ;
  double w_h4sio4= - KD_H4SiO4*grd_p_l- KF_H4SiO4*grd_h4sio4  ;
  double w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3     - Kpsi_CO3*grd_psi   ; 
  double w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3                      ;
  double w_na    = - KD_Na*grd_p_l    - KF_Na*grd_na       - Kpsi_Na*grd_psi    ;
  double w_naoh  = - KD_NaOH*grd_p_l  - KF_NaOH*grd_naoh                        ;
  double w_nahco3= - KD_NaHCO3*grd_p_l- KF_NaHCO3*grd_nahco3                    ;
  double w_naco3 = - KD_NaCO3*grd_p_l - KF_NaCO3*grd_naco3 - Kpsi_NaCO3*grd_psi ;
  double w_k     = - KD_K*grd_p_l     - KF_K*grd_k         - Kpsi_K*grd_psi     ;
  double w_koh   = - KD_KOH*grd_p_l   - KF_KOH*grd_koh                          ;
  double w_caoh  = - KD_CaOH*grd_p_l  - KF_CaOH*grd_caoh   - Kpsi_CaOH*grd_psi  ;
  double w_cahco3= - KD_CaHCO3*grd_p_l- KF_CaHCO3*grd_cahco3-Kpsi_CaHCO3*grd_psi;
  double w_caco3aq=- KD_CaCO3aq*grd_p_l-KF_CaCO3aq*grd_caco3aq                  ;
  double w_caoh2aq=- KD_CaOH2aq*grd_p_l-KF_CaOH2aq*grd_caoh2aq                  ;
  double w_cah3sio4= - KD_CaH3SiO4*grd_p_l- KF_CaH3SiO4*grd_cah3sio4 - Kpsi_CaH3SiO4*grd_psi  ;
  double w_cah2sio4= - KD_CaH2SiO4*grd_p_l- KF_CaH2SiO4*grd_cah2sio4            ;
  double w_co2   =                    - KF_CO2*grd_co2                          ;
  double w_m     = - KD_m*grd_p_l + M_CO2*w_co2 ;
  double w_q     =                    - z_h*KF_H*grd_h             \
                                        - z_oh*KF_OH*grd_oh          \
                                        - z_hco3*KF_HCO3*grd_hco3    \
                                        - z_co3*KF_CO3*grd_co3       \
                                        - z_ca*KF_Ca*grd_ca          \
                                        - z_na*KF_Na*grd_na          \
                                        - z_naco3*KF_NaCO3*grd_naco3 \
                                        - z_k*KF_K*grd_k             \
                                        - z_caoh*KF_CaOH*grd_caoh    \
                                        - z_cahco3*KF_CaHCO3*grd_cahco3 \
                                        - z_h3sio4*KF_H3SiO4*grd_h3sio4 \
                                        - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                                        - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                                                    - Kpsi_q*grd_psi ;


  w[I_W_C ] = w_co2 + w_h2co3 + w_hco3 + w_co3 + w_nahco3 + w_naco3 + w_cahco3 + w_caco3aq ;
  w[I_W_Ca] = w_ca + w_caoh + w_cahco3 + w_caco3aq + w_caoh2aq + w_cah2sio4 + w_cah3sio4 ;
  w[I_W_Na] = w_na + w_naoh + w_nahco3 + w_naco3 ;
  w[I_W_m ] = w_m ;
  w[I_W_Si] = w_h3sio4 + w_h4sio4 + w_h2sio4 + w_cah2sio4 + w_cah3sio4 ;
  w[I_W_q ] = w_q ;
  w[I_W_K ] = w_k + w_koh ;
    
  return(w) ;
}




int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[Element_MaxNbOfNodes] ;
  double *u_n[Element_MaxNbOfNodes] ;
  double x1    = Element_GetNodeCoordinate(el,1)[0] ;
  double x0    = Element_GetNodeCoordinate(el,0)[0] ;
  double dij   = x1 - x0 ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
    u_n[i] = Element_GetPreviousNodalUnknown(el,i) ;
  }
  
  /* Initialization */
  for(i = 0 ; i < nn*nn*NEQ*NEQ ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < 2 ; i++) {
    int j = 1 - i ;
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    /* Transfer terms from node i to node j */
    double *cij = c + (i*nn + j)*NEQ*NEQ ;
    /* Liquid and gas components */
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    dxi[I_C_CO2 ] = x[I_C_CO2]*1.e-2 ;
    dxi[I_C_OH  ] = x[I_C_OH ]*1.e-2 ;
    dxi[I_C_Na  ] = 1.e-4 ;
    dxi[I_C_K   ] = 1.e-4 ;
    dxi[I_N_CC  ] = x[I_N_CC  ]*1.e-2 ;
    dxi[I_N_Si_S] = x[I_N_Si_S]*1.e-2 ;
    dxi[I_P_L   ] = P_Ln(i)*1.e-6 ;
    dxi[I_PSI   ] = 1. ;
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      cii[E_C*NEQ    + k] = dx[I_N_C] ;
      cii[E_Ca*NEQ   + k] = dx[I_N_Ca] ;
      cii[E_Na*NEQ   + k] = dx[I_N_Na] ;
      cii[E_Si*NEQ   + k] = dx[I_N_Si] ;
      cii[E_K*NEQ    + k] = dx[I_N_K] ;
      cii[E_mass*NEQ + k] = dx[I_Mass] ;
      cii[E_el*NEQ   + k] = dx[I_N_Q] ;
      
      cij[E_C*NEQ    + k] = - dtdij*dw[I_W_C] ;
      cij[E_Ca*NEQ   + k] = - dtdij*dw[I_W_Ca] ;
      cij[E_Na*NEQ   + k] = - dtdij*dw[I_W_Na] ;
      cij[E_Si*NEQ   + k] = - dtdij*dw[I_W_Si] ;
      cij[E_K*NEQ    + k] = - dtdij*dw[I_W_K] ;
      cij[E_mass*NEQ + k] = - dtdij*dw[I_W_m] ;
      cij[E_q*NEQ    + k] = - dtdij*dw[I_W_q] ;
    }
  }

  return(dec) ;
}


void concentration1(double x_co2, double x_na_tot, double x_k_tot, double *pointeur_x_oh, double *pointeur_x_h3sio4, double *pointeur_x_na, double *pointeur_x_k)
{
/*
Resolution du systeme a quatre equations, quatre inconnues
avec f(x_oh,x_h3sio4,x_na,x_k)=0 electroneutralitée
     g(x_oh,x_h3sio4)=0 phases du CSH
     h_na(x_oh,x_na)=0 conservation Na
     k_k(x_oh,x_k)=0 conservation K
*/
  double x_oh_nouveau, x_h3sio4_nouveau, x_na_nouveau, x_k_nouveau ;
  double x_na = x_na_tot ;
  double x_k   = x_k_tot ;
  double x_oh = x_na + x_k +0.02 ;
  double x_co2_eq = x_co2 ;

  double x_h3sio4 = x_oh*pow(k_jen,1/y_jen)/pow(k_2,x_jen/y_jen) ;
  
  /*
  double x_h4sio4   = IonActivityProductOfSH(s_ch) ;
  double x_h3sio4   = (x_h4sio4*x_oh)/(k_e*k_h4sio4/k_h2sio4) ;
  */


  /* fonction f(x_oh,x_h3sio4,x_na,x_k) = A*x_oh^4+B*x_oh^3+ C*x_oh^2 + D*x_oh^1 + E*x_oh^0 
  + F*x_na*x_oh^2 + Gx_k*x_oh^2
  + K*x_na*x_oh^4
  + I*x_h3sio4*xoh^2 + H*x_h3sio4*x_oh^3 + J*x_h3sio4
 */
  double A = z_co3*x_co2_eq*k_co3*k_h/k_1 ;
  double B = z_oh + z_hco3*x_co2_eq*k_h/k_1 ;
  double C = 0 ;
  double D = z_h*k_e + z_caoh*k_caoh*k_2/k_e + z_cahco3*k_cahco3*k_h*k_2*x_co2_eq/k_1 ;
  double E = z_ca*k_2 ;
  double F = z_na ;
  double G = z_k ;
  double I   = z_h3sio4 ;
  double J  = z_cah3sio4*k_cah3sio4*k_2 ;
  double H  = z_h2sio4/(k_e*k_h2sio4) ;
  double K  = z_naco3*k_naco3*k_h*x_co2_eq/(k_e*k_1) ;

  double x_oh2= x_oh*x_oh ;
  double x_oh3= x_oh2*x_oh ;
  double x_oh4= x_oh2*x_oh2 ;
  double f ;
  double dfsdx_oh ;
  double dfsdx_h3sio4 ;
  double dfsdx_na ;
  double dfsdx_k ;

 
  /* fonction g(x_oh,x_h3sio4) = jen*x_oh^-y_jen*x_h3sio4^y_jen 
                               + tobII*x_oh^-y_tobII*x_h3sio4^y_tobII
                               + tobI*x_oh^-y_tobI*x_h3sio4^y_tobI
                               + sil*x_oh^-y_sil*x_h3sio4^y_sil - 1 */
  double jen = pow(k_2,x_jen)/k_jen ;
  double tobII = pow(k_2,x_tobII)/k_tobII  ;
  double tobI = pow(k_2,x_tobI)/k_tobI  ;
  double sil = k_h4sio4*k_e/(k_sil*k_h2sio4) ;
  double g  ;
  double dgsdx_oh ;
  double dgsdx_h3sio4  ;

  /* fonction h_na(x_oh,x_na) =  A_Na*x_na + B_Na*x_na*x_oh + C_Na*x_na*x_oh^2 -x_na_tot */
  double A_Na = 1 ;
  double B_Na = k_naoh/k_e + k_nahco3*k_h*x_co2_eq/k_1 ;
  double C_Na = k_naco3*k_h*x_co2_eq/(k_1*k_e) ;
  double h_na ;
  double dh_nasdx_oh ;
  double dh_nasdx_na ;

  /* fonction h_k(x_oh,x_k) =  A_K*x_k + B_K*x_k*x_oh - x_k_tot */
  double A_K = 1 ;
  double B_K = k_koh/k_e  ;
  double h_k ;
  double dh_ksdx_oh ;
  double dh_ksdx_k ;

  double a11,a12,a13,a14,a21,a22,a31,a33,a41,a44 ;
  
  double det ;
  
  double tolerance = 1e-3 ;
  double delta_oh  = 1 ;
  double delta_h3sio4 = 1 ;
  double delta_na  = 1 ;
  double delta_k  = 1 ;

  double f_jen,f_tobI,f_tobII,f_sil ;
 

  do
  {
  x_oh2= x_oh*x_oh ;
  x_oh3= x_oh2*x_oh ;
  x_oh4= x_oh2*x_oh2 ;

  f = A*x_oh4 + B*x_oh3+ C*x_oh2 + D*x_oh + E + F*x_na*x_oh2 + G*x_k*x_oh2 + K*x_na*x_oh4 + I*x_h3sio4*x_oh2 + H*x_h3sio4*x_oh3 + J*x_h3sio4 ;
  dfsdx_oh     	= 4*A*x_oh3 + 3*B*x_oh2 + 2*C*x_oh + D + 2*(F*x_na + G*x_k)*x_oh + 4*K*x_na*x_oh3 + 2*I*x_h3sio4*x_oh + 3*H*x_h3sio4*x_oh2 ;
  dfsdx_na     	= (F*x_oh2 + K*x_oh4) ;
  dfsdx_k     	= G*x_oh2 ;
  dfsdx_h3sio4 = I*x_oh2 + H*x_oh3 + J ;
   
    double x_hco3 = k_h*x_co2_eq*x_oh/k_1 ;
    double x_co3   = k_co3*x_oh*x_hco3 ;
    double x_h     = k_e/x_oh ;
    double x_ca    = k_ca/x_co3 ;
    double x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    double x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    double x_cahco3  = k_cahco3*x_hco3*x_ca ;
    double x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    double x_h2sio4  = x_oh*x_h3sio4/(k_e*k_h2sio4) ;
    double x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca ;

    double f_bis = (z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 + z_na*x_na + z_naco3*x_naco3 + z_k*x_k + z_caoh*x_caoh + z_cahco3*x_cahco3 + z_cah3sio4*x_cah3sio4 + z_h2sio4*x_h2sio4 + z_h3sio4*x_h3sio4)*x_oh2 ;


  g = jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen) + tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII) + tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI) + sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil) - 1 ;
  dgsdx_oh     = - y_jen*jen*pow(x_oh,-y_jen-1)*pow(x_h3sio4,y_jen) - y_tobII*tobII*pow(x_oh,-y_tobII-1)*pow(x_h3sio4,y_tobII) -y_tobI*tobI*pow(x_oh,-y_tobI-1)*pow(x_h3sio4,y_tobI) - y_sil*sil*pow(x_oh,-y_sil-1)*pow(x_h3sio4,y_sil) ;
  dgsdx_h3sio4 = y_jen*jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen-1) + y_tobII*tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII-1) + y_tobI*tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI-1) + y_sil*sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil-1) ;

  f_jen = jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen) ;
  f_tobII = tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII) ;
  f_tobI = tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI) ;
  f_sil = sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil) ;

   h_na = A_Na*x_na + B_Na*x_na*x_oh + C_Na*x_na*x_oh2 - x_na_tot ;
  dh_nasdx_na = A_Na + B_Na*x_oh + C_Na*x_oh2 ;
  dh_nasdx_oh = B_Na*x_na + 2*C_Na*x_na*x_oh ;

   h_k = A_K*x_k + B_K*x_k*x_oh - x_k_tot ;
  dh_ksdx_k = A_K + B_K*x_oh ;
  dh_ksdx_oh = B_K*x_k ;

  a11 = dfsdx_oh ;
  a12 = dfsdx_h3sio4 ;
  a13 = dfsdx_na ;
  a14 = dfsdx_k ;

  a21 = dgsdx_oh ;
  a22 = dgsdx_h3sio4 ;

  a31 = dh_nasdx_oh ;
  a33 = dh_nasdx_na ;

  a41 = dh_ksdx_oh ;
  a44 = dh_ksdx_k ;
  
  det =  a11*a22*a33*a44-a12*a21*a33*a44-a13*a22*a31*a44-a14*a22*a33*a41 ;




  x_oh_nouveau     = x_oh - (-a13*a22*a44*h_na-a14*a22*a33*h_k-a12*a33*a44*g+a22*a33*a44*f)/det ;
  x_h3sio4_nouveau = x_h3sio4 - (a13*a21*a44*h_na+a14*a21*a33*h_k+((a11*a33-a13*a31)*a44-a14*a33*a41)*g-a21*a33*a44*f)/det ;
  x_na_nouveau = x_na - (((a11*a22-a12*a21)*a44-a14*a22*a41)*h_na+a14*a22*a31*h_k+a12*a31*a44*g-a22*a31*a44*f)/det ;
  x_k_nouveau = x_k - (a13*a22*a41*h_na+((a11*a22-a12*a21)*a33-a13*a22*a31)*h_k+a12*a33*a41*g-a22*a33*a41*f)/det ; 

  delta_oh     = fabs((x_oh_nouveau     - x_oh)/x_oh) ;
  delta_h3sio4 = fabs((x_h3sio4_nouveau - x_h3sio4)/x_h3sio4) ;
  delta_na     = fabs((x_na_nouveau     - x_na)/x_na) ;
  delta_k     = fabs((x_k_nouveau     - x_k)/x_k) ;
  
  x_oh = x_oh_nouveau ;
  x_h3sio4 = x_h3sio4_nouveau ;
  x_na = x_na_nouveau ;
  x_k = x_k_nouveau ;


  
  } while(delta_oh > tolerance && delta_h3sio4 > tolerance && delta_na > tolerance && delta_k > tolerance) ;
  
printf("\n electro= %e \n", f) ;
printf("CSH= %e \n", g) ;
printf("concer Na= %e \n", h_na) ;
printf("concer K= %e \n", h_k) ;
printf("\n x_oh= %e \n", x_oh) ;
printf("x_h3sio4= %e \n", x_h3sio4) ;
printf("x_na= %e \n", x_na) ;
printf("x_k= %e \n", x_k) ;
printf("f_jen= %e \n", f_jen) ;
printf("f_tobII= %e \n", f_tobII) ;
printf("f_tobI= %e \n", f_tobI) ;
printf("f_sil= %e \n", f_sil) ;
printf("Tolerances \n") ;
printf("delta_oh= %e \n", delta_oh) ;
printf("delta_h3sio4= %e \n", delta_h3sio4) ;
printf("delta_na= %e \n", delta_na) ;
printf("delta_k= %e \n", delta_k) ;


  *pointeur_x_oh = x_oh ;
  *pointeur_x_h3sio4 = x_h3sio4 ;
  *pointeur_x_na = x_na ;
  *pointeur_x_k = x_k ;
 
}

double dn1_caoh2sdt(double av,double c)
{
#define AV         ((av < 1.) ? av : 1.)
  double rp = (av < 1.) ? pow(1 - av,1./3) : 0. ;
  double rc = pow(1 - AV + V_CC/V_CH*AV,1./3) ;
  /*
  double alpha2 = -1./3*av*av - 2./3*av + 1 ;
  double alpha  = -5.29478*av*av*av*av + 8.6069*av*av*av - 4.2444*av*av + 0.9325*av ;
  */

  return((rc > 0.) ? rp*rp/(1 + c*rp*(1 - rp/rc)) : 0.) ;
  /* return(alpha2/(1 + c_2*alpha)) ; */
#undef AV
}



double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *v0 = Element_GetConstantTerm(el) ;
  double *x = Components ;
  
  /* Primary Variables */
  x[I_C_CO2 ] = C_CO2(n) ;
  x[I_C_OH  ] = C_OH(n) ;
  x[I_C_Na  ] = C_Na(n) ;
  x[I_C_K   ] = C_K(n) ;
  x[I_N_CC  ] = N_CC(n) ;
  x[I_N_Si_S] = N_Si_S(n) ;
  x[I_P_L   ] = P_L(n) ;
  x[I_PSI   ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn]  = N_CHn(n) ;
  x[I_S_CH_EQn]  = S_CH_EQn(n) ;
  x[I_V_S0 ]  = V_S0(n) ;
  
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  dx[I_C_CO2 ] = x[I_C_CO2] ;
  dx[I_C_OH  ] = x[I_C_OH ] ;
  dx[I_C_Na  ] = x[I_C_Na ] ;
  dx[I_C_K   ] = x[I_C_K  ] ;
  dx[I_N_CC  ] = x[I_N_CC   ] ;
  dx[I_N_Si_S] = x[I_N_Si_S ] ;
  dx[I_P_L   ] = x[I_P_L  ] ;
  dx[I_PSI   ] = x[I_PSI  ] ;

  dx[I_N_CHn] = x[I_N_CHn] ;
  dx[I_S_CH_EQn]  = x[I_S_CH_EQn] ;
  dx[I_V_S0 ] = x[I_V_S0] ;
  
  dx[i] += dxi ;
  
  ComputeSecondaryComponents(el,dt,dx) ;
  
  for(j = 0 ; j < NbOfComponents ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}


void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)                   
{
  double x_co2      = x[I_C_CO2 ] ;
  double x_oh       = x[I_C_OH  ] ;
  double x_na       = x[I_C_Na  ] ;
  double x_k        = x[I_C_K   ] ;
  double n_cc       = x[I_N_CC  ] ;
  double n_si_s     = x[I_N_Si_S] ;
  double p_l        = x[I_P_L   ] ;
  
  /* Liquid components */
  double x_h        = k_e/x_oh ;
  
  double x_h2co3    = k_h*x_co2 ;
  double x_hco3     = x_oh*x_h2co3/k_1 ;
  double x_co3      = k_co3*x_oh*x_hco3 ;
  
  double x_ca       = k_ca/x_co3 ;
  double x_caoh     = k_caoh/k_e*x_oh*x_ca ;
  double x_cahco3   = k_cahco3*x_hco3*x_ca ;
  double x_caco3aq  = k_caco3*k_ca/(k_e*k_co3);
  double x_caoh2aq  = k_caoh2*x_ca*x_oh*x_oh ;
  
  double s_ch       = SaturationDegreeOfCH(x_co2) ;
  double s_cheqn    = x[I_S_CH_EQn] ;
  double s_cheq     = SaturationDegreeOfCHAtEq(el,dt,s_ch,s_cheqn) ;
  double x_h4sio4   = IonActivityProductOfSH(s_cheq,s_ch) ;
  double x_h3sio4   = (x_h4sio4*x_oh)/(k_e*k_h4sio4/k_h2sio4) ;
  double x_h2sio4   = x_oh*x_h3sio4/(k_e*k_h2sio4) ;
  double x_cah2sio4 = k_cah2sio4*x_h2sio4*x_ca ;
  double x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca ;
  
  double x_naoh     = k_naoh/k_e*x_na*x_oh ;
  double x_nahco3   = k_nahco3*x_na*x_hco3 ;
  double x_naco3    = k_naco3/k_e*x_na*x_oh*x_hco3 ;
  double x_koh      = k_koh/k_e*x_k*x_oh ;
    
  double x_h2o      = (1 - (x_h*v_h + x_oh*v_oh \
         + x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 \
         + x_ca*v_ca + x_caoh*v_caoh + x_caoh2aq*v_caoh2aq \
         + x_cahco3*v_cahco3 + x_caco3aq*v_caco3aq \
         + x_h3sio4*v_h3sio4 + x_h4sio4*v_h4sio4 + x_h2sio4*v_h2sio4 \
         + x_cah2sio4*v_cah2sio4 + x_cah3sio4*v_cah3sio4 \
         + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 \
         + x_k*v_k + x_koh*v_koh))/v_h2o ;
  
  double x_q        = z_h*x_h + z_oh*x_oh \
         + z_hco3*x_hco3 + z_co3*x_co3 \
         + z_ca*x_ca + z_caoh*x_caoh + z_cahco3*x_cahco3 \
         + z_h3sio4*x_h3sio4 + z_h2sio4*x_h2sio4 \
         + z_cah3sio4*x_cah3sio4 \
         + z_na*x_na + z_naco3*x_naco3 \
         + z_k*x_k ;
       
    
  /* Solid contents */
  /* ... as components: CH, CC, CSH */
  double n_chn      = x[I_N_CHn] ;
  double x_csh      = CalciumSiliconRatioInCSH(s_cheq) ;
  double z_csh      = WaterSiliconRatioInCSH(s_cheq) ;
  /* double av = n_cc/n_ch0 ; */
  double av = 1 - n_chn/n_ch0 ;
  double dn1sdt = a_2*dn1_caoh2sdt(av,c_2) ;
  /*cinétique de dissolution de la portlandite en loi log*/
  double dn_chsdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ;    
  /*cinétique de dissolution de la portlandite en loi puissance*/ 
  /* double dn_chsdt = - dn1sdt*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),Xp)) ; */
  double n_ch    = MAX(n_chn + dt*dn_chsdt , 0.) ;
  /* ... as elements: C, Ca, Si */
  double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  double n_c_s      = n_cc ;
  /* ... as mass */
  double m_csh      = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
  double m_s        = M_CaOH2*n_ch + M_CaCO3*n_cc + m_csh ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(s_cheq) ;
  double v_s        = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;
  
  /* Porosity */
  double v_s0     = x[I_V_S0] ;
  double phi      = phii + v_s0 - v_s ;
  
  /* Saturation */
  double p_c      = p_g - p_l ;
  double s_l      = SaturationDegree(p_c) ;
  double s_g      = 1 - s_l ;
  
  /* Liquid contents */
  double phi_l    = phi*s_l ;
  /* ... as elements: C, Ca, Si */
  double x_c_l  = x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq ;
  double x_ca_l = x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq + x_cah2sio4 + x_cah3sio4 ;
  double x_na_l = x_na + x_naoh + x_nahco3 + x_naco3 ;
  double x_k_l  = x_k + x_koh ;
  double x_si_l = x_h2sio4 + x_h3sio4 + x_h4sio4 + x_cah2sio4 + x_cah3sio4 ;
  double n_c_l  = phi_l*x_c_l ;
  double n_ca_l = phi_l*x_ca_l ;
  double n_na_l = phi_l*x_na_l ;
  double n_k_l  = phi_l*x_k_l ;
  double n_si_l = phi_l*x_si_l ;
  /* ... as mass */
  double rho_l  = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o \
       + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 \
       + M_Ca*x_ca + M_CaOH*x_caoh + M_CaHCO3*x_cahco3 \
       + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq \
       + M_H3SiO4*x_h3sio4 + M_H4SiO4*x_h4sio4 + M_H2SiO4*x_h2sio4 \
       + M_CaH2SiO4*x_cah2sio4 + M_CaH3SiO4*x_cah3sio4 \
       + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 \
       + M_K*x_k + M_KOH*x_koh ;
  double m_l    = phi_l*rho_l ;
       
  /* Gas contents */
  double phi_g  = phi*s_g ;
  /* ... as elements */
  double n_c_g  = phi_g*x_co2 ;
  /* ... as mass */
  double rho_g  = M_CO2*x_co2 ;
  double m_g    = phi_g*rho_g ;

  /* Back up */
  
  /* Liquid components */
  x[I_C_H       ] = x_h ;
  x[I_C_OH      ] = x_oh ;
  x[I_C_H2O     ] = x_h2o ;
  
  x[I_C_CO2     ] = x_co2 ;
  x[I_C_HCO3    ] = x_hco3 ;
  x[I_C_H2CO3   ] = x_h2co3 ;
  x[I_C_CO3     ] = x_co3 ;
  
  x[I_C_Ca      ] = x_ca ;
  x[I_C_CaOH    ] = x_caoh ;
  x[I_C_CaHCO3  ] = x_cahco3 ;
  x[I_C_CaCO3aq ] = x_caco3aq ;
  x[I_C_CaOH2aq ] = x_caoh2aq ;
  
  x[I_C_H4SiO4  ] = x_h4sio4 ;
  x[I_C_H3SiO4  ] = x_h3sio4 ;
  x[I_C_H2SiO4  ] = x_h2sio4 ;
  x[I_C_CaH2SiO4] = x_cah2sio4 ;
  x[I_C_CaH3SiO4] = x_cah3sio4 ;
  
  x[I_C_Na      ] = x_na ;
  x[I_C_NaOH    ] = x_naoh ;
  x[I_C_NaHCO3  ] = x_nahco3 ;
  x[I_C_NaCO3   ] = x_naco3 ;
  
  x[I_C_K       ] = x_k ;
  x[I_C_KOH     ] = x_koh ;
  
  x[I_S_CH      ] = s_ch ;
  
  x[I_S_CH_EQ   ] = s_cheq ;
  
  x[I_RHO_L     ] = rho_l ;
  
  /* Solid components */
  x[I_N_CH    ] = n_ch ;
  x[I_V_S     ] = v_s ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  /* Element contents */
  x[I_N_C ]  = n_c_l  + n_c_s  + n_c_g ;
  x[I_N_Ca]  = n_ca_l + n_ca_s ;
  x[I_N_Na]  = n_na_l ; 
  x[I_N_K ]  = n_k_l  ;
  x[I_N_Si]  = n_si_l + n_si_s ;
  /* Total mass */
  x[I_Mass]  = m_g + m_l + m_s ;
  /* Charge density */
  x[I_N_Q]   = x_q ;
    
  return ;
}




double ComputeS_CH_EQ(Element_t *el,double dt,double s_ch,double s_cheq_n)
{
  double err,tol = 1.e-4 ;
  double x_n   = CalciumSiliconRatioInCSH(s_cheq_n) ;
  double x_ch  = CalciumSiliconRatioInCSH(s_ch) ;
  double lns_cheq_n  = log(s_cheq_n) ;
  double lns_ch  = log(s_ch) ;
  double lns ;
  double s ;
  int    i = 0 ;
  
  if(dt == 0 || s_ch == s_cheq_n) {
    return(s_cheq_n) ;
  } else {
    /* Assume a linear x-lns relationship */
    double sum = tau*(x_n - x_ch) +  dt*(lns_cheq_n - lns_ch) ;
    double prd = tau*(x_n - x_ch) *  (lns_cheq_n - lns_ch) ;
    lns = lns_ch - prd/sum ;
    s = exp(lns) ;
  }
  
  do {
    double x  = CalciumSiliconRatioInCSH(s) ;
    double dx = dCalciumSiliconRatioInCSH(s)*s ;
    double f  = tau*(x  - x_n) + dt*(lns  - lns_ch) ;
    double df = tau*dx + dt ; 
    double dlns = (df != 0) ? - f/df : 0 ;
    lns += dlns ;
    s  = exp(lns) ;
    /*
    if(lns > MAX(lns_ch,lns_cheq_n)) lns = MAX(lns_ch,lns_cheq_n) ;
    if(lns < MIN(lns_ch,lns_cheq_n)) lns = MIN(lns_ch,lns_cheq_n) ;
    */
    err = fabs(dlns) ;
    if(i++ > 50) {
      printf("\n") ;
      printf("s_cheq_n   = %e\n",s_cheq_n) ;
      printf("s_ch = %e\n",s_ch) ;
      printf("x_n  = %e\n",x_n) ;
      printf("x_ch = %e\n",x_ch) ;
      printf("x    = %e\n",x) ;
      printf("s    = %e\n",s) ;
      printf("err  = %e\n",err) ;
      printf("f    = %e\n",f) ;
      printf("df   = %e\n",df) ;
      printf("dlns = %e\n",dlns) ;
      arret("ComputeS_CH_EQ : no convergence") ;
    }
  } while(err > tol) ;
  
  return(s) ;
}




double ComputeS_CH_EQ1(Element_t *el,double dt,double s_ch,double s_cheq_n)
{
  double err,tol = 1.e-6 ;
  double x_n     = CalciumSiliconRatioInCSH(s_cheq_n) ;
  double lns_ch  = log(s_ch) ;
  double s1    = s_ch ;
  double s2    = s_cheq_n ;
  int    i = 0 ;
  
  if(s_ch == s_cheq_n) {
    return(s_ch) ;
  }
  
  do {
    double x1   = CalciumSiliconRatioInCSH(s1) ;
    double x2   = CalciumSiliconRatioInCSH(s2) ;
    double lns1 = log(s1) ;
    double lns2 = log(s2) ;
    double f1   = tau*(x1 - x_n) + dt*(lns1 - lns_ch) ;
    double f2   = tau*(x2 - x_n) + dt*(lns2 - lns_ch) ;
    double lns  = (f1 == f2) ? lns1 : (lns2*f1 - lns1*f2)/(f1 - f2) ;
    double s    = exp(lns) ;
    double x    = CalciumSiliconRatioInCSH(s) ;
    double f    = tau*(x  - x_n) + dt*(lns  - lns_ch) ;
    
    s1 = s ;
    
    err = fabs(f) ;
    if(i++ > 20) {
      printf("\n") ;
      printf("s_cheq_n   = %e\n",s_cheq_n) ;
      printf("s_ch = %e\n",s_ch) ;
      printf("x_n  = %e\n",x_n) ;
      printf("x    = %e\n",x) ;
      printf("s    = %e\n",s) ;
      printf("err  = %e\n",err) ;
      printf("f    = %e\n",f) ;
      arret("ComputeS_CH_EQ1 (2) : no convergence") ;
    }
  } while(err > tol) ;
  
  return(s1) ;
}
