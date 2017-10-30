/* General features of the model:
 * Curves for CSH:
 *   - C/S ratio
 *   - H/S ratio
 *   - Molar Volume
 * Alkalis (as sodium and potassium compounds):
 * (Na[+], K[+], NaOH[0], KOH[0], NaHCO3[0], NaCO3[-])
 * Dissolution kinetics for CH based on spherical crystal coated 
 * by a calcite layer.
 * Dissolution and continuous decalcification of CSH
 * Split CC in Vaterite and Calcite
 * Precipitation/Dissolution of CC Vaterite
 * Precipitation kinetics of CC Calcite
 * Use of Zeta unknowns for Calcium and Silicon
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

/* The Finite Volume Method */
#include "FVM.h"

#define TITLE   "Carbonation Of CBM: CC as vaterite (equ.) and calcite (kin.) (2015)"
#define AUTHORS "Morandeau-Thiery-Dangla"

#include "PredefinedMethods.h"

/* Macros */

#define NEQ    	  (7)
#define NVE    	  (58)
#define NVI       (25)
#define NV0       (2)

#define E_C       (0)
#define E_q       (1)
#define E_mass    (2)
#define E_Ca      (3)
#define E_Na      (5)
#define E_K       (6)
#define E_Si      (4)

#define U_C_CO2   (0)
#define U_P_L     (2)
#define U_ZN_Ca_S (3)
#define U_PSI     (1)
#define U_C_Na    (5)
#define U_C_K     (6)
#define U_ZN_Si_S (4)

#define NOLOG_U   1
#define LOG_U     2
#define Ln10      Math_Ln10
#define U_CO2     LOG_U

#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])



#if (U_CO2 == LOG_U)
#define LogC_CO2(n)	  (UNKNOWN(n,U_C_CO2))
#define LogC_CO2n(n)	(UNKNOWNn(n,U_C_CO2))
#define C_CO2(n)	    (pow(10,UNKNOWN(n,U_C_CO2)))
#define C_CO2n(n)	    (pow(10,UNKNOWNn(n,U_C_CO2)))
#else
#define C_CO2(n)	    (UNKNOWN(n,U_C_CO2))
#define C_CO2n(n)	    (UNKNOWNn(n,U_C_CO2))
#define LogC_CO2(n)	  (log10(UNKNOWN(n,U_C_CO2)))
#define LogC_CO2n(n)	(log10(UNKNOWNn(n,U_C_CO2)))
#endif

#define ZN_Si_S(n)    (UNKNOWN(n,U_ZN_Si_S))
#define ZN_Ca_S(n)    (UNKNOWN(n,U_ZN_Ca_S))
#define P_L(n)        (UNKNOWN(n,U_P_L))
#define PSI(n)        (UNKNOWN(n,U_PSI))
#define C_Na(n)       (UNKNOWN(n,U_C_Na))
#define C_K(n)        (UNKNOWN(n,U_C_K))

#define ZN_Si_Sn(n)   (UNKNOWNn(n,U_ZN_Si_S))
#define ZN_Ca_Sn(n)   (UNKNOWNn(n,U_ZN_Ca_S))
#define P_Ln(n)       (UNKNOWNn(n,U_P_L))
#define PSIn(n)       (UNKNOWNn(n,U_PSI))
#define C_Nan(n)      (UNKNOWNn(n,U_C_Na))
#define C_Kn(n)       (UNKNOWNn(n,U_C_K))


#define N_C(n)        (f[(n)])
#define N_q(n)        (f[(2+n)])
#define Mass(n)       (f[(4+n)])
#define N_Ca(n)       (f[(6+n)])
#define N_Na(n)       (f[(8+n)])
#define N_K(n)        (f[(10+n)])
#define N_Si(n)       (f[(12+n)])
#define W_C           (f[14])
#define W_q           (f[15])
#define W_m           (f[16])
#define W_Ca          (f[17])
#define W_Na          (f[18])
#define W_K           (f[19])
#define W_Si          (f[20])
#define N_CH(n)       (f[(21+n)])
#define N_Calcite(n)  (f[(23+n)])

#define N_Cn(n)       (f_n[(n)])
#define N_qn(n)       (f_n[(2+n)])
#define Mass_n(n)     (f_n[(4+n)])
#define N_Can(n)      (f_n[(6+n)])
#define N_Nan(n)      (f_n[(8+n)])
#define N_Kn(n)       (f_n[(10+n)])
#define N_Sin(n)      (f_n[(12+n)])
#define N_CHn(n)      (f_n[(21+n)])
#define N_Calciten(n) (f_n[(23+n)])

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

#define KF_CO2          (va[(22)])
#define KF_Ca           (va[(23)])
#define KF_OH           (va[(24)])
#define KF_H            (va[(25)])
#define KF_H2CO3        (va[(26)])
#define KF_HCO3         (va[(27)])
#define KF_CO3          (va[(28)])
#define KF_Na           (va[(29)])
#define KF_NaOH         (va[(30)])
#define KF_NaHCO3     	(va[(31)])
#define KF_NaCO3      	(va[(32)])
#define KF_K            (va[(33)])
#define KF_KOH          (va[(34)])
#define KF_CaOH         (va[(35)])
#define KF_CaHCO3       (va[(36)])
#define KF_CaCO3aq      (va[(37)])
#define KF_CaOH2aq      (va[(38)])
#define KF_H3SiO4       (va[(39)])
#define KF_H2SiO4       (va[(40)])
#define KF_H4SiO4       (va[(41)])
#define KF_CaH2SiO4     (va[(42)])
#define KF_CaH3SiO4     (va[(43)])

#define Kpsi_Ca         (va[(44)])
#define Kpsi_OH         (va[(45)])
#define Kpsi_H          (va[(46)])
#define Kpsi_HCO3       (va[(47)])
#define Kpsi_CO3        (va[(48)])
#define Kpsi_Na         (va[(49)])
#define Kpsi_NaCO3      (va[(50)])
#define Kpsi_q          (va[(51)])
#define Kpsi_K          (va[(52)])
#define Kpsi_CaOH       (va[(53)])
#define Kpsi_CaHCO3     (va[(54)])
#define Kpsi_H3SiO4     (va[(55)])
#define Kpsi_H2SiO4     (va[(56)])
#define Kpsi_CaH3SiO4   (va[(57)])

#define V_S0(n)         (v0[(0+n)])


#define TEMPERATURE  (298)


/* Charge of ions */
#include "ElectricChargeOfIonInWater.h"
#define z_h           ElectricChargeOfIonInWater(H)
#define z_oh          ElectricChargeOfIonInWater(OH)
#define z_hco3        ElectricChargeOfIonInWater(HCO3)
#define z_co3         ElectricChargeOfIonInWater(CO3)
#define z_ca          ElectricChargeOfIonInWater(Ca)
#define z_caoh        ElectricChargeOfIonInWater(CaOH)
#define z_cahco3      ElectricChargeOfIonInWater(CaHCO3)
#define z_h3sio4      ElectricChargeOfIonInWater(H3SiO4)
#define z_h2sio4      ElectricChargeOfIonInWater(H2SiO4)
#define z_cahco3      ElectricChargeOfIonInWater(CaHCO3)
#define z_cah3sio4    ElectricChargeOfIonInWater(CaH3SiO4)
#define z_na          ElectricChargeOfIonInWater(Na)
#define z_naco3       ElectricChargeOfIonInWater(NaCO3)
#define z_k           ElectricChargeOfIonInWater(K)


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


/* Molar Masses */
#include "MolarMassOfMolecule.h"
#define gr             (1.e3)
#define M_Ca           (MolarMassOfMolecule(Ca)*gr)
#define M_H2CO3        (MolarMassOfMolecule(H2CO3)*gr)
#define M_HCO3         (MolarMassOfMolecule(HCO3)*gr)
#define M_CO3          (MolarMassOfMolecule(CO3)*gr)
#define M_OH           (MolarMassOfMolecule(OH)*gr)
#define M_H            (MolarMassOfMolecule(H)*gr)
#define M_H2O          (MolarMassOfMolecule(H2O)*gr)
#define M_Na           (MolarMassOfMolecule(Na)*gr)
#define M_NaOH         (MolarMassOfMolecule(NaOH)*gr)
#define M_NaHCO3       (MolarMassOfMolecule(NaHCO3)*gr)
#define M_NaCO3        (MolarMassOfMolecule(NaCO3)*gr)
#define M_CO2          (MolarMassOfMolecule(CO2)*gr)
#define M_CaOH2        (MolarMassOfMolecule(CaO2H2)*gr)
#define M_CaCO3        (MolarMassOfMolecule(CaCO3)*gr)
#define M_CaO          (MolarMassOfMolecule(CaO)*gr)
#define M_SiO2         (MolarMassOfMolecule(SiO2)*gr)
#define M_K            (MolarMassOfMolecule(K)*gr)
#define M_KOH          (MolarMassOfMolecule(KOH)*gr)
#define M_CaOH         (MolarMassOfMolecule(CaOH)*gr)
#define M_CaHCO3       (MolarMassOfMolecule(CaHCO3)*gr)
#define M_CaCO3aq      (MolarMassOfMolecule(CaCO3)*gr)
#define M_CaOH2aq      (MolarMassOfMolecule(CaO2H2)*gr)
#define M_CaH2SiO4     (MolarMassOfMolecule(CaH2SiO4)*gr)
#define M_CaH3SiO4     (MolarMassOfMolecule(CaH3SiO4)*gr)
#define M_H3SiO4       (MolarMassOfMolecule(H3SiO4)*gr)
#define M_H2SiO4       (MolarMassOfMolecule(H2SiO4)*gr)
#define M_H4SiO4       (MolarMassOfMolecule(H4SiO4)*gr)


/* Equilibrium constant of homogeneous aqueous reactions */
#define k_e       		(1.e-14)                /* H2O = OH[-] + H[+] */

#define k_h        		(1.)                    /* CO2(g) + H2O = H2CO3 */
#define k_co3      		(4.570881896148751e3)   /* CO3[2-] + H2O = HCO3[-] + OH[-] */
#define k_1        		(2.187761623949552e-8)  /* HCO3[-] + H2O = H2CO3 + OH[-] */

#define k_naoh     		(6.60693448e-15)        /* Na[+] + OH[-] = NaOH */
#define k_nahco3     	(0.5623413252)          /* Na[+] + HCO3[-] = NaHCO3  */
#define k_naco3    	  (8.72971368e-10)        /* Na[+] + HCO3[-] = NaCO3[-] + H[+] */

#define k_caoh     		(1.65958691e-13)        /* Ca[2+] + H2O = CaOH[+] + H[+] */
#define k_cahco3   	  (12.76438809)           /* Ca[2+] + HCO3[-] = CaHCO3[+] */
#define k_caco3      	(7.852356346e-8)        /* Ca[2+] + HCO3[-] = CaCO3[0] + H[+] */
#define k_caoh2    	  (1.)                    /* Ca(OH)2[0] = Ca[2+] + 2OH[-] */

#define k_koh      		(3.4673685e-15)         /* K[+] + H2O = KOH[-] + H[+] */

#define k_h2sio4     	(2.13796209e+13)        /* H2SiO4[2-] + H[+] = H3SiO4[-] */
#define k_h4sio4     	(1.380384265e+23)       /* H2SiO4[2-] + 2H[+] = H4SiO4 */

#define k_cah2sio4 	  (39810.71)              /* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define k_cah3sio4 	  (15.84)                 /* Ca[2+] + H3SiO4[-] = CaH3SiO4[+] */


/* Threshold of CO2 concentration */
#define C_CO2_eq                         (k_ca*k_1/(k_h*k_co3*k_2))


/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp   	   	 (0.5)


/*Loi de Davies pour la prise en compte de l'activite ionique*/
#define A_DAVIES  (0.5)
#define B_DAVIES  (0.24)


/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CC   = Calcium Carbonate (Calcite)
  CSH  = Calcium Silicate Hydrate
  SH   = Amorphous Silica Gel
*/

/* Material Properties */
#define SaturationDegree(p)              (Curve_ComputeValue(Element_GetCurve(el),p))
#define RelativePermeabilityToLiquid(p)  (Curve_ComputeValue(Element_GetCurve(el) + 1,p))


/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(q)      (Curve_ComputeValue(Element_GetCurve(el) + 2,q))
#define WaterSiliconRatioInCSH(q)        (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define MolarVolumeOfCSH(q)              (Curve_ComputeValue(Element_GetCurve(el) + 4,q))


/* S-H Properties */
/* Equilibrium constant */
#define k_sil      		(1.936421964e-3)        /* S-H = H4SiO4 + (t - 2)H2O */
/* Saturation degree of dissolved S-H */
#define S_SHeq(q)                        (Curve_ComputeValue(Element_GetCurve(el) + 5,q))
#define SaturationDegreeOfSH(s_ch,zn_si_s)    (NEGEXP(zn_si_s)*S_SHeq(s_ch))
/* Ion Activity Product of dissolved S-H */
#define IonActivityProductOfSH(s_sh)     (k_sil*(s_sh))


/* CH Properties */
/* Equilibrium constant */
#define k_2        		(6.456542290346550e-6)  /* CH = Ca[2+] + 2OH[-] */
/* Molar volume of CH solid (dm3/mole) */
#define V_CH	    	  (33.e-3)
/* Saturation Degree of Dissolved CH */
#define SaturationDegreeOfCH(z_co2,zn_ca_s)   (NEGEXP(zn_ca_s)/MAX((z_co2),1.))



/* CC Properties */
/* Equilibrium constant */
#define k_ca       		(12.218e-9)  /* CC = Ca[2+] + CO3[2-] */
/* Molar volume of CC (dm3/mole) */
#define V_CC	      	(37.e-3)
/* Saturation Degree of Dissolved CC */
/*
#define SaturationDegreeOfCC(z_co2,zn_ca_s)   (NEGEXP(zn_ca_s)*MIN((z_co2),1.))
*/
#define SaturationDegreeOfCC(z_co2,s_ch)      ((s_ch)*(z_co2))
/* Ion Activity Product of Dissolved CC */
#define IonActivityProductOfCC(s_cc)        (k_ca*(s_cc))
/* Calcite Properties */
/* Equilibrium constant */
#define k_calcite      (3.3113e-9)  /* CC = Ca[2+] + CO3[2-] */
#define r_calcite      (3.69)       /* ratio of equilibrium constants K_Vaterite/K_Calcite */
#define SaturationDegreeOfCalcite(s_cc)      (r_calcite*(s_cc))


/* Element contents in solid phases  */
#define n_ca_ref                           (n_ch0)
#define n_si_ref                           (n_csh0)
#define CalciumContentInCHAndCC(zn_ca_s)   (n_ca_ref*MAX(zn_ca_s,0.))
#define SiliconContentInCSH(zn_si_s)       (n_si_ref*MAX(zn_si_s,0.))


/* Concentration of OH computed from electroneutrality */
#define ConcentrationOfOHInLiquid(A,B,C,D,E)  concentration_oh(A,B,C,D,E)


#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)


/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static double  dn1_caoh2sdt(double,double) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;


static double  PermeabilityCoefficient(Element_t*,double) ;
static void    ComputePhysicoChemicalProperties(double) ;

static double concentration_oh(double,double,double,double,double) ;
static double concentrations_oh_na_k(double,double,double,double*,double*) ;
static double poly4(double,double,double,double,double) ;


/* Internal parameters */
static double phii ;
static double k_int,frac,phi_r ;
static double a_2,c_2 ;
static double tau_ch ;
static double n_ch0,n_csh0,x_na0,x_k0 ;
static double n_refcalcite, tau_calcite ;
static double c_co2_eq ;
static double p_g = 0. ;

static double d_h,d_oh ;
static double d_co2,d_h2co3,d_hco3,d_co3 ;
static double d_ca,d_caoh,d_caoh2aq ;
static double d_cahco3,d_caco3aq ;
static double d_h4sio4,d_h3sio4,d_h2sio4 ;
static double d_cah2sio4,d_cah3sio4 ;
static double d_k,d_koh ;
static double d_na,d_naoh,d_nahco3,d_naco3 ;

static double mu_l ;
static double FRT ;


#include "WaterViscosity.h"
#include "DiffusionCoefficientOfMoleculeInWater.h"
#include "PhysicalConstant.h"
#define dm2            (1.e2)


void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (dm2/s) */
  d_oh         = DiffusionCoefficientOfMoleculeInWater(OH,TK)*dm2 ;
  d_h          = DiffusionCoefficientOfMoleculeInWater(H,TK)*dm2 ;
  d_hco3       = DiffusionCoefficientOfMoleculeInWater(HCO3,TK)*dm2 ;
  d_h2co3      = DiffusionCoefficientOfMoleculeInWater(H2CO3,TK)*dm2 ;
  d_co3        = DiffusionCoefficientOfMoleculeInWater(CO3,TK)*dm2 ;
  
  d_ca         = DiffusionCoefficientOfMoleculeInWater(Ca,TK)*dm2 ;
  d_caoh       = DiffusionCoefficientOfMoleculeInWater(CaOH,TK)*dm2 ;
  d_cahco3     = DiffusionCoefficientOfMoleculeInWater(CaHCO3,TK)*dm2 ;
  d_caco3aq    = DiffusionCoefficientOfMoleculeInWater(CaCO3,TK)*dm2 ;
  d_caoh2aq  	 = 7.92e-8 ;
  
  d_h4sio4     = DiffusionCoefficientOfMoleculeInWater(H4SiO4,TK)*dm2 ;
  d_h3sio4     = DiffusionCoefficientOfMoleculeInWater(H3SiO4,TK)*dm2 ;
  d_h2sio4     = DiffusionCoefficientOfMoleculeInWater(H2SiO4,TK)*dm2 ;
  
  d_cah2sio4   = DiffusionCoefficientOfMoleculeInWater(CaH2SiO4,TK)*dm2 ;
  d_cah3sio4   = DiffusionCoefficientOfMoleculeInWater(CaH3SiO4,TK)*dm2 ;
  
  d_na         = DiffusionCoefficientOfMoleculeInWater(Na,TK)*dm2*100 ;  /* (1.33e-7) */
  d_naoh       = 1.33e-7 ;  /* (1.33e-7) */
  d_nahco3     = 1.33e-7 ;  /* (1.33e-7) */
  d_naco3      = 1.33e-7 ;  /* (1.33e-7) */
  d_k          = DiffusionCoefficientOfMoleculeInWater(K,TK)*dm2*100 ; /* (1.957e-7) */
  d_koh        = DiffusionCoefficientOfMoleculeInWater(KOH,TK)*dm2*100 ; /* (1.957e-7) */

  /* Diffusion Coefficient Of Molecules In Air (dm2/s) */
  d_co2      	 = 1.6e-3 ;
  
  /* Viscosity (Pa.s) */
  mu_l       = WaterViscosity(TK) ;
  
  /* Physical constants */
  {
    double RT      = PhysicalConstant(PerfectGasConstant)*TK*1.e3 ;
    double Faraday = PhysicalConstant(Faraday)*1.e3 ;
    FRT     = Faraday/RT ;
  }
}

#define NbOfComponents    (49)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;


#define I_C_OH         (7)
#define I_C_H          (8)
#define I_C_H2O        (9)

#define I_C_CO2        (10)
#define I_C_HCO3       (11)
#define I_C_H2CO3      (12)
#define I_C_CO3        (13)

#define I_C_Ca         (14)
#define I_C_CaOH       (15)
#define I_C_CaHCO3     (16)
#define I_C_CaCO3aq    (17)
#define I_C_CaOH2aq    (18)

#define I_C_H2SiO4     (19)
#define I_C_H3SiO4     (20)
#define I_C_H4SiO4     (21)

#define I_C_CaH2SiO4   (22)
#define I_C_CaH3SiO4   (23)

#define I_C_Na         (24)
#define I_C_NaOH       (25)
#define I_C_NaHCO3     (26)
#define I_C_NaCO3      (27)

#define I_C_K          (28)
#define I_C_KOH        (29)

#define I_S_CH         (30)

#define I_P_L          (31)
#define I_RHO_L        (32)

#define I_N_C          (33)
#define I_N_Ca         (34)
#define I_N_Si         (35)
#define I_N_K          (36)
#define I_N_Na         (37)
#define I_Mass         (38)
#define I_N_Q          (39)
#define I_N_CC         (40)
#define I_N_Si_S       (41)

#define I_N_CH         (42)
#define I_V_S          (43)

#define I_N_CHn        (44)
#define I_V_S0         (45)

#define I_Phi          (46)

#define I_PSI          (47)

#define I_N_Calciten   (48)

#define NbOfComponentFluxes    (7)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_C           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_Na          (3)
#define I_W_K           (4)
#define I_W_m           (5)
#define I_W_q           (6)


int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0)     return (0) ;
  else if(strcmp(s,"k_int") == 0)   return (1) ;
  else if(strcmp(s,"N_CaOH2") == 0) return (2) ;
  else if(strcmp(s,"C_CO2_eq") == 0) return (3) ;
  else if(strcmp(s,"N_Si") == 0)    return (4) ;
  else if(strcmp(s,"X_K") == 0)     return (5) ;
  else if(strcmp(s,"X_Na") == 0)    return (6) ;
  else if(strcmp(s,"A_1") == 0)     return (7) ;
  else if(strcmp(s,"A_2") == 0)     return (8) ;
  else if(strcmp(s,"C_2") == 0)     return (9) ;
  else if(strcmp(s,"R_CaOH2") == 0) return (10) ;
  else if(strcmp(s,"D") == 0) 	    return (11) ;
  else if(strcmp(s,"Tau") == 0)     return (12) ;
  else if(strcmp(s,"frac") == 0)    return (13) ;
  else if(strcmp(s,"phi_r") == 0)   return (14) ;
  else if(strcmp(s,"tau_calcite") == 0)   return (15) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_csh0   = GetProperty("N_Si") ;
  x_na0    = GetProperty("X_Na") ;
  x_k0     = GetProperty("X_K") ;
  frac     = GetProperty("frac") ;
  phi_r    = GetProperty("phi_r") ;
  c_co2_eq = GetProperty("C_CO2_eq") ;
  tau_ch   = GetProperty("Tau") ;
  tau_calcite = GetProperty("tau_calcite") ;
  n_refcalcite = n_ch0 ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_C   ,"carbone") ;
  Model_CopyNameOfEquation(model,E_q   ,"charge") ;
  Model_CopyNameOfEquation(model,E_mass,"masse") ;
  Model_CopyNameOfEquation(model,E_Ca  ,"calcium") ;
  Model_CopyNameOfEquation(model,E_Na  ,"sodium") ;
  Model_CopyNameOfEquation(model,E_K   ,"potassium") ;
  Model_CopyNameOfEquation(model,E_Si  ,"silicon") ;
  
  
#if (U_CO2 == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_CO2 ,"logc_co2") ;
#else
  Model_CopyNameOfUnknown(model,U_C_CO2 ,"c_co2") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,U_P_L    ,"p_l") ;
  Model_CopyNameOfUnknown(model,U_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,U_PSI    ,"psi") ;
  Model_CopyNameOfUnknown(model,U_C_Na   ,"c_na") ;
  Model_CopyNameOfUnknown(model,U_C_K    ,"c_k") ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 16 ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Initialisation automatique */
  {
    double h   = 5.6e-6 ;  /* (mol/dm2/s) these MT p 223 */
    double R_0 = Material_GetProperty(mat)[pm("R_CaOH2")] ; /* (dm) */
    double D   = Material_GetProperty(mat)[pm("D")] ; /* (mol/dm/s) */
    
    if(R_0 == 0.) R_0 = 40e-5 ;
    
    if(D == 0.) D = 7e-15 ;
    
    /* contenu molaire de reference en CaOH2 */
    n_ch0 = Material_GetProperty(mat)[pm("N_CaOH2")] ;
    
    {
      double t_ch = Material_GetProperty(mat)[pm("Tau")] ; /* (s) */
      /* double t_ch = R_0/(3*h*V_CH) ; */ /* (s) approx 721.5 s */
      /* double t_ch = R_0*R_0/(3*V_CH*D) ; */ /* (s) approx 2.3e8 s */
      
      if(t_ch == 0) {
        t_ch = R_0/(3*h*V_CH) ;     /* (s) approx 721.5 s */
        /* t_ch = R_0*R_0/(3*V_CH*D) ; */ /* (s) approx 2.3e8 s */
        Material_GetProperty(mat)[pm("Tau")] = t_ch ;
      }
      
      a_2 = n_ch0/t_ch ;  /* (mol/dm3/s) these MT p 227 */
    }
    
    c_2 = h*R_0/D ;     /* (no dim) approx 3.2e5 these MT p 228 */
    
    /*
    a_1 = Material_GetProperty(mat)[pm("A_1")] ;
    if(a_1 == 0.) a_1 = 6000. ;
    */
  }
  
  /* Material_GetProperty(mat)[pm("A_1")] = a_1 ; */
  Material_GetProperty(mat)[pm("A_2")] = a_2 ;
  Material_GetProperty(mat)[pm("C_2")] = c_2 ;
  
    
  /* contenu molaire de reference en CSH */
  n_csh0 = Material_GetProperty(mat)[pm("N_Si")] ;
  if(n_csh0 == 0) n_csh0 = 1. ;
  Material_GetProperty(mat)[pm("N_Si")] = n_csh0 ;
  
  {
    c_co2_eq = C_CO2_eq ;
    Material_GetProperty(mat)[pm("C_CO2_eq")] = c_co2_eq ;
  }
  
  return(n_donnees) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The set of 7 equations is:\n") ;
  printf("\t- Mass balance of C      (carbone)\n") ;
  printf("\t- Mass balance of Ca     (calcium)\n") ;
  printf("\t- Mass balance of Si     (silicon)\n") ;
  printf("\t- Mass balance of Na     (sodium)\n") ;
  printf("\t- Mass balance of K      (potassium)\n") ;
  printf("\t- Total mass balance     (mass)\n") ;
  printf("\t- Charge balance         (charge)\n") ;
  
  printf("\n") ;
  printf("The 7 primary unknowns are:\n") ;
  printf("\t- Liquid pressure                  (p_l)\n") ;
  printf("\t- Electric potential               (psi) \n") ;
  printf("\t- Carbon dioxide gas concentration (c_co2)\n") ;
  printf("\t- Potassium concentration          (c_k)\n") ;
  printf("\t- Sodium concentration             (c_na)\n") ;
  printf("\t- Zeta unknown for calcium         (z_ca)\n") ;
  printf("\t   \t z_ca is defined as:\n") ;
  printf("\t   \t z_ca = n_ch/n0 + log(s_ch)  for c_co2 < c_co2_eq\n") ;
  printf("\t   \t z_ca = n_cc/n0 + log(s_cc)  for c_co2 > c_co2_eq\n") ;
  printf("\t- Zeta unknown for silicon         (z_si)\n") ;
  printf("\t   \t z_si is defined as:\n") ;
  printf("\t   \t z_si = n_si/n0 + log(s_sh/s_sh_eq)\n") ;
  
  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length    : dm !\n") ;
  printf("\t time      : s !\n") ;
  printf("\t pressure  : Pa !\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;


  fprintf(ficd,"porosite = 0.38   # Porosity\n") ;
  fprintf(ficd,"k_int = 1.4e-17   # Intrinsic permeability (dm2)\n") ;
  fprintf(ficd,"N_CaOH2 = 3.9     # Initial content in Ca(OH)2 (mol/L)\n") ;
  fprintf(ficd,"R_CaOH2 = 40.e-5  # Portlandite crystal radius \n") ;
  fprintf(ficd,"N_Si = 2.4        # Initial content in CSH (mol/L)\n") ;
  fprintf(ficd,"T_csh = 1.8e-4    # k_h/T_csh = Characteristic time of CSH carbonation (1/s)\n") ;
  fprintf(ficd,"X_Na = 0.019      # Total content in Na (mol/L)\n") ;
  fprintf(ficd,"X_K  = 0.012      # Total content in K  (mol/L)\n") ;
  fprintf(ficd,"A_1 = 6000        # Kinetic coef 1 (dm/mol/s)\n") ;
  fprintf(ficd,"D = 7.e-15        # Diffusion coef in CC (dm/mol/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Kinetic coef 2 (dm/mol/s)\n") ;
  fprintf(ficd,"frac = 0.8        # Fractionnal length of pore bodies\n") ;
  fprintf(ficd,"phi_r = 0.7       # Porosity for which permeability vanishes\n") ;
  fprintf(ficd,"Curves = my_file  # File name: p_c S_l k_rl C/S H/S V_csh\n") ;  

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
  double* f  = Element_GetImplicitTerm(el) ;
  double* v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int i ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  {
    double x_na_tot 	= x_na0 ;
    double x_k_tot  	= x_k0 ;

    for(i = 0 ; i < nn ; i++) {
      double x_na       = x_na_tot ;
      double x_k        = x_k_tot ;
      double x_co2      = C_CO2(i) ;
      double z_co2      = x_co2/c_co2_eq ;
      double zn_ca_s    = ZN_Ca_S(i) ;
      double zn_si_s    = ZN_Si_S(i) ;
      double s_ch       = SaturationDegreeOfCH(z_co2,zn_ca_s) ;
      double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
      double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
      double n_ch       = (z_co2 > 1) ? 0       : n_ch_cc ;
      double n_cc       = (z_co2 > 1) ? n_ch_cc : 0 ;
      double n_si_s     = SiliconContentInCSH(zn_si_s) ;
      double v_csh      = MolarVolumeOfCSH(s_ch) ;
      double v_s0       = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;
    
      concentrations_oh_na_k(x_co2,s_ch,s_sh,&x_na,&x_k) ;

      C_Na(i)     = x_na ;
      C_K(i)      = x_k ;
      
      /* Solid contents */
      V_S0(i)    = v_s0 ;
      N_CH(i)    = n_ch ;
      N_Calcite(i)    = 0 ;
      
      if(z_co2 > 1) {
        printf("z_co2    = %e\n",z_co2) ;
        arret("ComputeInitialState (MTD3)") ;
      }
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
  double*  f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /*
    Coefficients de transfert
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int i ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
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

    {
      double x_co2      = x[I_C_CO2] ;
      double x_oh    	  = x[I_C_OH] ;
      double x_na    	  = x[I_C_Na] ;
      double x_k     	  = x[I_C_K] ;
      double x_ca    	  = x[I_C_Ca] ;
      double x_h2o      = x[I_C_H2O] ;
      double n_si_s     = x[I_N_Si_S] ;
      
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
        printf("x_ca     = %e\n",x_ca) ;
        printf("n_si_s   = %e\n",n_si_s) ;
        printf("x_naoh   = %e\n",x_naoh) ;
        printf("x_nahco3 = %e\n",x_nahco3) ;
        printf("x_naco3  = %e\n",x_naco3) ;
        return(1) ;
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
  FVM_t *fvm = FVM_GetInstance(el) ;
  double c[4*NEQ*NEQ] ;
  int    i ;
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }

#if (U_CO2 == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_CO2)     *= Ln10*C_CO2(0) ;
      K(i,U_C_CO2+NEQ) *= Ln10*C_CO2(1) ;
    }
  }
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *f   = Element_GetCurrentImplicitTerm(el) ;
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
    Conservation of elements (C,Ca,Si,Na,K), total mass and charge
  */
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
  
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;
  
  R(0,E_mass) -= volume[0]*(Mass(0) - Mass_n(0)) + dt*surf*W_m ;
  R(1,E_mass) -= volume[1]*(Mass(1) - Mass_n(1)) - dt*surf*W_m ;
  
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*W_Na ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*W_Na ;
  
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ;

  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    nso = 48 ;
  int    i ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double *x = ComputeComponents(el,u,f,0,j) ;

    double p_l        = x[I_P_L] ;
    double p_c        = p_g - p_l ;
    /* saturation */
    double s_l        = SaturationDegree(p_c) ;
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
    double zn_si_s    = x[U_ZN_Si_S] ;
    double n_si_s     = x[I_N_Si_S] ;
    double n_cc 	    = x[I_N_CC] ;
    double n_ch       = x[I_N_CH] ;
    double s_ch       = x[I_S_CH] ;
    double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
    
    double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double v_solide_csh   = v_csh*n_si_s ;
    double v_solide_ch    = V_CH*n_ch ;
    double v_solide_cc    = V_CC*n_cc ;

    double n_chn      = x[I_N_CHn] ;
    double av         = 1 - n_chn/n_ch0 ;
    double dn1sdt     = a_2*dn1_caoh2sdt(av,c_2) ;
    double dn_chsdt   = dn1sdt*log(s_ch) ;
    double coeff_dnCH = log(s_ch) ;
  
    double CsurS      = (x_csh + n_ch/n_si_s) ;
    

    /* porosite */
    double phi        = x[I_Phi] ;

    double ph = 14 + log10(x_oh) ;
    double n_Na = 0.5*(N_Na(0) + N_Na(1)) ;
    double n_Ca = 0.5*(N_Ca(0) + N_Ca(1)) ;
    double n_Si = 0.5*(N_Si(0) + N_Si(1)) ; 


    /* Transferts */
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx        = x1 - x0 ;
    double grd_psi   = (PSI(1) - PSI(0))/dx ;
    
    double coeff_permeability = PermeabilityCoefficient(el,phi) ;
	  double k_l  = (k_int/mu_l)*RelativePermeabilityToLiquid(p_c)*coeff_permeability ;


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
    
    /* Added by A. Morandeau */
    Result_Store(r + i++,&x_csh,"x_csh",1) ;
    Result_Store(r + i++,&n_si_s,"n_si_s",1) ;
    Result_Store(r + i++,&s_ch,"s_ch",1) ;
    Result_Store(r + i++,&s_sh,"s_sh",1) ;
    Result_Store(r + i++,&k_l,"k_l",1) ;
    Result_Store(r + i++,&coeff_permeability,"verma-pruess",1) ;
    Result_Store(r + i++,&dn_chsdt,"dn_chsdt",1) ;
    Result_Store(r + i++,&dn1sdt,"dn1sdt",1) ;
    Result_Store(r + i++,&coeff_dnCH,"coeff_dnCH",1) ;
    Result_Store(r + i++,&v_solide_csh,"v_csh",1) ;
    Result_Store(r + i++,&v_solide_ch,"v_ch",1) ;
    Result_Store(r + i++,&v_solide_cc,"v_cc",1) ;
  }
  
  
  if(i != nso) arret("ComputeOutputs") ;
  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int    i ; 

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
	
    /* Permeability */
    double coeff_permeability = PermeabilityCoefficient(el,phi) ; /* Change here */
    double k_l  = (k_int/mu_l)*RelativePermeabilityToLiquid(p_c)*coeff_permeability ;
    
    
    /* tortuosite gaz */
    double s_g = 1 - s_l ;
    double phi_g = phi*s_g ;
   /* double tau  = pow(phi,1/3)*pow(s_g,7/3) ; */
    double tau  = pow(phi,1.74)*pow(s_g,3.20) ;

    
    /* tortuosite liquide */
    double phi_cap = phi/2  ;
    double phi_c   = 0.17 ; /*Percolation capilar porosity*/
    
    /*Diffusivity according to Oh and Jang, CCR203*/
		double n = 2.7 ; 		/* OPC n = 2.7  --------------  Fly ash n = 4.5 */
    double ds_norm = 5e-5 ;	/* OPC ds_norm = 1e-4 --------  Fly ash ds_norm = 5e-5 */
		double m_phi = 0.5*( pow(ds_norm,1/n) + phi_cap/(1-phi_c)*(1 - pow(ds_norm,1/n)) - phi_c/(1-phi_c)) ;
    double iff =  pow( m_phi + sqrt( m_phi*m_phi +  pow(ds_norm,1/n)*phi_c/(1-phi_c)),n)*pow(s_l,4.5) ;
    
    /*Diffusivity : ITZ for mortars and concrete */

    /*Diffusivity according to Bazant et Najjar */
		/*double iff    = 0.00029*exp(9.95*phi)/(1+625*pow((1-s_l),4)) ;*/

    /* Humidité relative */
    /* double hr = exp(-p_c*M_H2O/(RT*rho_l)) ; */
 
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
    
    KF_CO2      += phi_g*tau*d_co2 ;
    /* KF_CO2    	+= (1.6e-3)*pow(phi,1.8)*pow(1-hr,2.2) ; */
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
    
    Kpsi_Ca   	 += FRT*KF_Ca*z_ca*x_ca ;
    Kpsi_HCO3 	 += FRT*KF_HCO3*z_hco3*x_hco3 ;
    Kpsi_CO3  	 += FRT*KF_CO3*z_co3*x_co3 ;
    Kpsi_OH   	 += FRT*KF_OH*z_oh*x_oh ;
    Kpsi_H    	 += FRT*KF_H*z_h*x_h ;
    Kpsi_Na    	 += FRT*KF_Na*z_na*x_na ;
    Kpsi_NaCO3 	 += FRT*KF_NaCO3*z_naco3*x_naco3 ;
    Kpsi_K     	 += FRT*KF_K*z_k*x_k ;
    Kpsi_CaOH 	 += FRT*KF_CaOH*z_caoh*x_caoh ;
    Kpsi_CaHCO3	 += FRT*KF_CaHCO3*z_cahco3*x_cahco3 ;
    Kpsi_H3SiO4	 += FRT*KF_H3SiO4*z_h3sio4*x_h3sio4 ;
    Kpsi_H2SiO4	 += FRT*KF_H2SiO4*z_h2sio4*x_h2sio4 ;
    Kpsi_CaH3SiO4+= FRT*KF_CaH3SiO4*z_cah3sio4*x_cah3sio4 ;
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
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  double x1    = Element_GetNodeCoordinate(el,1)[0] ;
  double x0    = Element_GetNodeCoordinate(el,0)[0] ;
  double dij   = fabs(x1 - x0) ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
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
    
    /*
    dxi[U_C_CO2  ] = x[U_C_CO2]*1.e-2 ;
    dxi[U_C_Na   ] = 1.e-4 ;
    dxi[U_C_K    ] = 1.e-4 ;
    dxi[U_ZN_Ca_S] = x[U_ZN_Ca_S]*((x[U_ZN_Ca_S] > ZN_Ca_Sn(i)) ? 1 : -1)*1.e-4 ;
    dxi[U_ZN_Si_S] = x[U_ZN_Si_S]*((x[U_ZN_Si_S] > ZN_Si_Sn(i)) ? 1 : -1)*1.e-4 ;
    dxi[U_P_L    ] = P_Ln(i)*1.e-6 ;
    dxi[U_PSI    ] = 1. ;
    */

    dxi[U_C_CO2  ] =  1.e-4*ObVal_GetValue(obval + U_C_CO2) ;
    dxi[U_C_Na   ] =  1.e-3*ObVal_GetValue(obval + U_C_Na) ;
    dxi[U_C_K    ] =  1.e-3*ObVal_GetValue(obval + U_C_K) ;
    dxi[U_ZN_Ca_S] =  1.e-4*ObVal_GetValue(obval + U_ZN_Ca_S) ;
    dxi[U_ZN_Si_S] =  1.e-4*ObVal_GetValue(obval + U_ZN_Si_S) ;
    dxi[U_P_L    ] =  1.e-4*ObVal_GetValue(obval + U_P_L) ;
    dxi[U_PSI    ] =  1.e+0*ObVal_GetValue(obval + U_PSI) ;
    
    dxi[U_ZN_Si_S] *= ((x[U_ZN_Si_S] > ZN_Si_Sn(i)) ? 1 : -1) ; 
    dxi[U_ZN_Ca_S] *= ((x[U_ZN_Ca_S] > ZN_Ca_Sn(i)) ? 1 : -1) ;
    
    #if (U_CO2 == LOG_U)
    dxi[U_C_CO2  ] *=  C_CO2n(i) ;
    #endif
    
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



double concentrations_oh_na_k(double x_co2,double s_ch,double s_sh,double *px_na,double *px_k)
{
/* Solve a set of 3 equations:
 * 1. Electroneutralilty
 * 2. Mass balance of Na
 * 3. Mass balance of K
 * Unknowns: x_oh, x_na, x_k.
 * On input, px_na and px_k point to total contents of Na and K
 */
  double x_na_tot = *px_na ;
  double x_k_tot  = *px_k ;
  
  /* Initialization */
  double x_na = x_na_tot ;
  double x_k  = x_k_tot ;
  double x_oh0 = x_na + x_k ;
  double x_oh = x_oh0 ;
  
  /* x_na_tot =  x_na * (A_Na + B_Na*x_oh + C_Na*x_oh*x_oh) */
  double A_Na = 1 ;
  double B_Na = k_naoh/k_e + k_nahco3*k_h*x_co2/k_1 ;
  double C_Na = k_naco3*k_h*x_co2/(k_1*k_e) ;

  /* x_k_tot =  x_k * (A_K + B_K*x_oh) */
  double A_K = 1 ;
  double B_K = k_koh/k_e  ;
  
  double err,tol = 1.e-8 ;
  int i = 0 ;
  
  do {
    double dx_oh = - x_oh ;
    
    x_oh = ConcentrationOfOHInLiquid(x_co2,s_ch,s_sh,x_na,x_k) ;
    
    x_na = x_na_tot/(A_Na + B_Na*x_oh + C_Na*x_oh*x_oh) ;
    
    x_k  = x_k_tot/(A_K + B_K*x_oh) ;
    
    dx_oh += x_oh ;
    
    err = fabs(dx_oh/x_oh) ;
    
    if(i++ > 20) {
      printf("x_na_tot = %e\n",x_na_tot) ;
      printf("x_na     = %e\n",x_na) ;
      printf("x_k_tot  = %e\n",x_k_tot) ;
      printf("x_oh0    = %e\n",x_oh0) ;;
      printf("x_oh     = %e\n",x_oh) ;
      arret("concentrations_oh_na_k : non convergence") ;
    }

  } while(err > tol || x_oh < 0) ;
  
  /*
  printf("\n") ;
  printf("x_oh = %e \n", x_oh) ;
  printf("x_na = %e \n", x_na) ;
  printf("x_k  = %e \n", x_k) ;
  */
  
  *px_na = x_na ;
  *px_k  = x_k ;
  return(x_oh) ;
}


double dn1_caoh2sdt(double av0,double c)
{
  double av = ((av0 > 0) ? ((av0 < 1) ? av0 : 1) : 0) ; /* av = 1 - n_ch/n_ch0 */
  double rp = (av < 1) ? pow(1 - av,1./3) : 0 ; /* rp = Rp/R0 */
  double rc = pow(1 - av + V_CC/V_CH*av,1./3) ; /* rc = Rc/R0 */
  double width = rc - rp ;
  double dn1dt = (rc > 0.) ? rp*rp/(1 + c*width*rp/rc) : 0 ;
  
  return(dn1dt) ;
}


double PermeabilityCoefficient(Element_t *el,double phi)
{
  double coeff_permeability ;
  
  {
    /* permeability Kozeny-Carman */
    double kozeny_carman  = pow(phi/phii,3.)*pow(((1 - phii)/(1 - phi)),2.) ;

    /* permeability Verma Pruess 1988 */
    /* frac = fractionnal length of pore bodies (0.8) */
    /* phi_r = fraction of initial porosity (phi/phi0) at which permeability is 0 (0.9) */

	  double S_s =  (phii - phi)/phii    ; /* saturation en solide */
	  double w = 1 + (1/frac)/(1/phi_r - 1) ;
    double t = (1 - S_s - phi_r)/(1 - phi_r) ;
	  double verma_pruess = (t > 0) ? t*t*(1 - frac + (frac/(w*w)))/(1 - frac + frac*(pow(t/(t + w - 1),2.))) : 0 ;
	
    /* permeability coefficient */
    coeff_permeability = verma_pruess ;
  }
  
  return(coeff_permeability) ;
}



double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *v0 = Element_GetConstantTerm(el) ;
  double *x = Components ;
  
  /* Primary Variables */
  x[U_C_CO2  ] = C_CO2(n) ;
  x[U_C_Na   ] = C_Na(n) ;
  x[U_C_K    ] = C_K(n) ;
  x[U_ZN_Ca_S] = ZN_Ca_S(n) ;
  x[U_ZN_Si_S] = ZN_Si_S(n) ;
  x[U_P_L    ] = P_L(n) ;
  x[U_PSI    ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn]  = N_CHn(n) ;
  x[I_V_S0 ]  = V_S0(n) ;
  x[I_N_Calciten] = N_Calciten(n) ;
  
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NEQ ; j++) {
    dx[j] = x[j] ;
  }

  /* Needed variables to compute secondary components */
  dx[I_N_CHn] = x[I_N_CHn] ;
  dx[I_V_S0 ] = x[I_V_S0] ;
  dx[I_N_Calciten] = x[I_N_Calciten] ;
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryComponents(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfComponents ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}


void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)
{
  double x_co2      = x[U_C_CO2  ] ;
  double x_na       = x[U_C_Na   ] ;
  double x_k        = x[U_C_K    ] ;
  double zn_ca_s    = x[U_ZN_Ca_S] ;
  double zn_si_s    = x[U_ZN_Si_S] ;
  double p_l        = x[U_P_L    ] ;
  
  
  /* Liquid components */
  double z_co2      = x_co2/c_co2_eq ;
  double s_ch       = SaturationDegreeOfCH(z_co2,zn_ca_s) ;
  double s_cc       = SaturationDegreeOfCC(z_co2,s_ch) ;
  double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
  
  double x_oh       = ConcentrationOfOHInLiquid(x_co2,s_ch,s_sh,x_na,x_k) ;
  double x_h        = k_e/x_oh ;
  
  double x_h2co3    = k_h*x_co2 ;
  double x_hco3     = x_oh*x_h2co3/k_1 ;
  double x_co3      = k_co3*x_oh*x_hco3 ;
  
  double Q_CC       = IonActivityProductOfCC(s_cc) ;
  double x_ca       = Q_CC/x_co3 ;
  double x_caoh     = k_caoh/k_e*x_oh*x_ca ;
  double x_cahco3   = k_cahco3*x_hco3*x_ca ;
  double x_caco3aq  = k_caco3*Q_CC/(k_e*k_co3);
  double x_caoh2aq  = k_caoh2*x_ca*x_oh*x_oh ;
  
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  double x_h4sio4   = Q_SH ;
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
  double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
  double n_chn      = x[I_N_CHn] ;
  double av         = 1 - n_chn/n_ch0 ;
  double dn1sdt     = a_2*dn1_caoh2sdt(av,c_2) ;
  double dn_chsdt   = dn1sdt*log(s_ch) ; /* Kinetics */
  double n_ch_ki    = MAX(n_chn + dt*dn_chsdt , 0.) ;
  double n_ch       = (z_co2 > 1) ? n_ch_ki : n_ch_cc ;
  double n_vaterite = (z_co2 > 1) ? n_ch_cc - n_ch_ki : 0 ;
  
  double n_calciten = x[I_N_Calciten] ;
  double s_calcite  = SaturationDegreeOfCalcite(s_cc) ;
  double dn_calcitesdt = n_refcalcite/tau_calcite*log(s_calcite) ;
  double n_calcite  = MAX(n_calciten + dt*dn_calcitesdt , 0.) ;
  
  double n_cc       = n_vaterite + n_calcite ;
  
  /* ... as elements: C, Ca, Si */
  double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double n_si_s     = SiliconContentInCSH(zn_si_s) ;
  double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  double n_c_s      = n_cc ;
  /* ... as mass */
  double z_csh      = WaterSiliconRatioInCSH(s_ch) ;
  double m_csh      = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
  double m_s        = M_CaOH2*n_ch + M_CaCO3*n_cc + m_csh ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(s_ch) ;
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
  
  x[I_RHO_L     ] = rho_l ;
  x[I_P_L       ] = p_l ;
  
  /* Solid components */
  x[I_N_CH    ] = n_ch ;
  x[I_V_S     ] = v_s ;
  x[I_N_Si_S  ] = n_si_s ;
  x[I_N_CC    ] = n_cc ;
  
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
  
  /* Electric potential */
  x[I_PSI]   = x[U_PSI] ;
    
  return ;
}



double concentration_oh(double x_co2,double s_ch,double s_sh,double x_na,double x_k)
/* Solve electroneutrality, SUM(z_i c_i) = 0, for x_h
   as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* The primary variables are considered as constant */
  double z_co2      = x_co2/c_co2_eq ;
  /* Ion activity products are constant as well */
  double s_cc       = SaturationDegreeOfCC(z_co2,s_ch) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  double Q_CC       = IonActivityProductOfCC(s_cc) ;
  /* Other constant concentrations */
  double x_h4sio4   = Q_SH ;
  double x_h2co3    = k_h*x_co2 ;
  /* Electroneutrality is written as  Sum z_i*x_i = 0, with
   * ------------------------------------------------------------------
   * x_i        =  A_i*(x_h)**n                              : n    z_i
   * ------------------------------------------------------------------
   * x_h        = x_h                                        : +1   +1
   * x_oh       = k_e/x_h                                    : -1   -1
   * 
   * x_h2co3    = k_h*x_co2                                  :  0    0
   * x_hco3     = k_e*x_h2co3/k_1/x_h                        : -1   -1
   * x_co3      = k_e*k_co3*x_hco3/x_h                       : -2   -2
   * 
   * x_ca       = Q_CC/x_co3                                 : +2   +2
   * x_caoh     = k_caoh*x_ca/x_h                            : +1   +1
   * x_cahco3   = k_cahco3*x_hco3*x_ca                       : +1   +1
   * x_caco3aq  = k_caco3*Q_CC/(k_e*k_co3)                   :  0    0
   * x_caoh2aq  = k_caoh2*Q_CH                               :  0    0
   * 
   * x_h4sio4   = Q_SH                                       :  0    0
   * x_h3sio4   = x_h4sio4/(k_h4sio4/k_h2sio4)/x_h           : -1   -1
   * x_h2sio4   = x_h3sio4/(k_h2sio4)/x_h                    : -2   -2
   * 
   * x_cah2sio4 = k_cah2sio4*x_h2sio4*x_ca                   :  0    0
   * x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca                   : +1   +1
   * 
   * x_na       = x_na                                       :  0   +1
   * x_naoh     = k_naoh*x_na/x_h                            : -1    0
   * x_nahco3   = k_nahco3*x_na*x_hco3                       : -1    0
   * x_naco3    = k_naco3*x_na*x_hco3/x_h                    : -2   -1
   * 
   * x_k        = x_k                                        :  0   +1
   * x_koh      = k_koh*x_k/x_h                              : -1    0
   *
   * We compute the A_i for i having non zero z_i
   */
  double A_hco3     = x_h2co3/k_1*k_e ;
  double A_co3      = k_e*k_co3*A_hco3 ;
  double A_ca       = Q_CC/A_co3 ;
  double A_caoh     = k_caoh*A_ca ;
  double A_cahco3   = k_cahco3*A_hco3*A_ca ;
  double A_h3sio4   = x_h4sio4/(k_h4sio4/k_h2sio4) ;
  double A_h2sio4   = A_h3sio4/(k_h2sio4) ;
  double A_cah3sio4 = k_cah3sio4*A_ca*A_h3sio4 ;
  
  double A_naco3    = k_naco3*x_na*A_hco3 ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahco3*A_cahco3 + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_na*x_na + z_k*x_k ;
  double d = z_oh*k_e + z_hco3*A_hco3 + z_h3sio4*A_h3sio4 ;
  double e = z_co3*A_co3 + z_h2sio4*A_h2sio4 + z_naco3*A_naco3 ;

  double x_h = poly4(a,b,c,d,e) ;
  
  if(x_h < 0) {
      printf("x_h   = %e\n",x_h) ;
      printf("a     = %e\n",a) ;
      printf("b     = %e\n",b) ;
      printf("c     = %e\n",c) ;
      printf("d     = %e\n",d) ;
      printf("e     = %e\n",e) ;
      printf("x_co2 = %e\n",x_co2) ;
      printf("s_ch  = %e\n",s_ch) ;
      printf("s_sh  = %e\n",s_sh) ;
      printf("x_na  = %e\n",x_na) ;
      printf("x_k   = %e\n",x_k) ;
      arret("concentration_oh : x_h<0") ;
  }
 
  return(k_e/x_h) ;
}


double poly4(double a,double b,double c,double d,double e)
/* Solve ax^4 + bx^3 + cx^2 + dx + e = 0 for x in the range [0:1] */
{
  double tol = 1e-4 ;
  double y[5] ;
  double x ;
  
  y[0] = a ;
  y[1] = b ;
  y[2] = c ;
  y[3] = d ;
  y[4] = e ;
  
  {
    int n = Math_ComputePolynomialEquationRoots(y,4) ;
    int i ;
    for(i = 0 ; i < n ; i++) {
      if((x = y[i]) < 1) break ;
    }
  }
  
  y[0] = a ;
  y[1] = b ;
  y[2] = c ;
  y[3] = d ;
  y[4] = e ;
  
  x = Math_PolishPolynomialEquationRoot(y,4,x,tol*x,20) ;
  
  return(x) ;
}
