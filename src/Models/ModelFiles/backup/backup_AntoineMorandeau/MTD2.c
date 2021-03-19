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
#include "CommonModel.h"

/* The Finite Volume Method */
#include "FVM.h"

#define TITLE   "Carbonation Of CBM with kinetics on C-S-H (2014)"
#define AUTHORS "Morandeau-Thiery-Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ    	(7)
#define NVE    	(58)
#define NVI     (32)
#define NV0     (2)

#define E_C       (0)
#define E_q       (1)
#define E_mass    (2)
#define E_Ca      (3)
#define E_Na      (5)
#define E_K       (6)
#define E_Si      (4)

#define U_C_CO2   (0)
#define U_P_L     (2)
#define U_N_CC    (3)
#define U_PSI     (1)
#define U_C_Na    (5)
#define U_C_K     (6)
#define U_N_Si_S  (4)

#define NOLOG_U   1
#define LOG_U     2
#define Ln10      2.302585093
#define U_CO2     LOG_U

#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


#if (U_CO2 == LOG_U)
#define LogC_CO2(n)	  (UNKNOWN(n,U_C_CO2))
#define LogC_CO2n(n)	(UNKNOWNn(n,U_C_CO2))
#define C_CO2(n)	    (exp(Ln10*UNKNOWN(n,U_C_CO2)))
#define C_CO2n(n)	    (exp(Ln10*UNKNOWNn(n,U_C_CO2)))
#else
#define C_CO2(n)	    (UNKNOWN(n,U_C_CO2))
#define C_CO2n(n)	    (UNKNOWNn(n,U_C_CO2))
#define LogC_CO2(n)	  (log10(UNKNOWN(n,U_C_CO2)))
#define LogC_CO2n(n)	(log10(UNKNOWNn(n,U_C_CO2)))
#endif


#define N_Si_S(n)   (UNKNOWN(n,U_N_Si_S))
#define N_CC(n)     (UNKNOWN(n,U_N_CC))
#define P_L(n)      (UNKNOWN(n,U_P_L))
#define PSI(n)      (UNKNOWN(n,U_PSI))
#define C_Na(n)     (UNKNOWN(n,U_C_Na))
#define C_K(n)      (UNKNOWN(n,U_C_K))

#define N_Si_Sn(n)   (UNKNOWNn(n,U_N_Si_S))
#define N_CCn(n)     (UNKNOWNn(n,U_N_CC))
#define P_Ln(n)      (UNKNOWNn(n,U_P_L))
#define PSIn(n)      (UNKNOWNn(n,U_PSI))
#define C_Nan(n)     (UNKNOWNn(n,U_C_Na))
#define C_Kn(n)      (UNKNOWNn(n,U_C_K))


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


/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp   	   	 (0.5)

/*Loi de Davies pour la prise en compte de l'activite ionique*/
#define A_DAVIES  (0.5)
#define B_DAVIES  (0.24)


/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CC   = Calcium Carbonate (Calcite)
  CSH  = Calcium Silicates Hydrate
  SH   = Amorphous Silica Gel
*/

/* Material Properties */
#define SaturationDegree(p)              (Curve_ComputeValue(Element_GetCurve(el),p))
#define RelativePermeabilityToLiquid(p)  (Curve_ComputeValue(Element_GetCurve(el) + 1,p))


/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(q)      (Curve_ComputeValue(Element_GetCurve(el) + 2,q))
#define dCalciumSiliconRatioInCSH(s)     (Curve_ComputeDerivative(Element_GetCurve(el) + 2,s))
#define WaterSiliconRatioInCSH(q)        (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define MolarVolumeOfCSH(q)              (Curve_ComputeValue(Element_GetCurve(el) + 4,q))
#define InverseOfCalciumSiliconRatioInCSH(x)      (Curve_ComputeValue(Element_GetCurve(el) + 6,x))


/* CH Properties */
/* Equilibrium constant */
#define k_2        		(6.456542290346550e-6)  /* CH = Ca[2+] + 2OH[-] */
/* Molar volume of CH solid (dm3/mole) */
#define V_CH	    	  (33.e-3)
#define C_CO2_eq                         (k_ca*k_1/(k_h*k_co3*k_2))
#define SaturationDegreeOfCH(c_co2)      (C_CO2_eq/c_co2)
#define SaturationDegreeOfCHAtEq(dt,s_ch,s_n,r)         ComputeS_CH_EQ(el,dt,s_ch,s_n,r)


/* S-H Properties */
/* Equilibrium constant */
#define k_sil      		(1.936421964e-3)        /* S-H = H4SiO4 + (t - 2)H2O */
/* Saturation degree of dissolved S-H */
#define S_SHeq(q)                        (Curve_ComputeValue(Element_GetCurve(el) + 5,q))
#define SaturationDegreeOfSHAtEq(q)      S_SHeq(q)
#define SaturationDegreeOfSH(dt,s,r)     ComputeS_SH(el,dt,s,r)
/* Ion Activity Product of dissolved S-H */
#define IonActivityProductOfSH(s_sh)     (k_sil*s_sh)


/* CC Properties */
/* Equilibrium constant */
#define k_ca       		(3.890451449942805e-9)  /* CC = Ca[2+] + CO3[2-] */
/* Molar volume of CC (dm3/mole) */
#define V_CC	      	(37.e-3)


/* Concentration of OH computed from electroneutrality */
#define ConcentrationOfOHInLiquid(A,B,C,D)  concentration_oh(A,B,C,D,el)


/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static double  dn1_caoh2sdt(double,double) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;


static double  PermeabilityCoefficient(Element_t*,double) ;
static void    ComputePhysicoChemicalProperties(double) ;

static double concentration_oh(double,double,double,double,Element_t*) ;
static double concentrations_oh_na_k(double,double*,double*,Element_t*) ;
static double poly4(double,double,double,double,double) ;

static double  ComputeS_CH_EQ(Element_t*,double,double,double,double) ;
static double  ComputeS_CH_EQ0(Element_t*,double,double,double) ;
static double  ComputeS_CH_EQ1(Element_t*,double,double,double,double) ;
static double  ComputeS_CH_EQ2(Element_t*,double,double,double,double) ;
static double  ComputeS_CH_EQ3(Element_t*,double,double,double,double) ;
static double  ComputeS_SH(Element_t*,double,double,double) ;


/* Parametres */
static double phii ;
static double k_int,frac,phi_r ;
static double a_2,c_2 ;
static double tau_ch ;
static double n_ch0,n_csh0,x_na0,x_k0 ;
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

static double tau_ca,tau_si ;


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

#define NbOfComponents    (54)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_C_CO2        (0)
#define I_P_L          (2)
#define I_N_CC         (3)
#define I_PSI          (1)
#define I_C_Na         (5)
#define I_C_K          (6)
#define I_N_Si_S       (4)

#define I_C_OH         (7)
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
#define I_N_Si_Sn      (52)

#define I_S_SH         (53)

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
  else if(strcmp(s,"tau_ca") == 0)   return (15) ;
  else if(strcmp(s,"tau_si") == 0)   return (16) ;
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
  tau_ca   = GetProperty("tau_ca") ;
  tau_si   = GetProperty("tau_si") ;
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
  Model_CopyNameOfUnknown(model,U_N_Si_S,"n_si_s") ;
  Model_CopyNameOfUnknown(model,U_P_L   ,"p_l") ;
  Model_CopyNameOfUnknown(model,U_N_CC  ,"c_caco3") ;
  Model_CopyNameOfUnknown(model,U_PSI   ,"psi") ;
  Model_CopyNameOfUnknown(model,U_C_Na  ,"c_na") ;
  Model_CopyNameOfUnknown(model,U_C_K   ,"c_k") ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 18 ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Initialisation automatique */
  {
    double h   = 5.6e-6 ;  /* (mol/dm2/s) these MT p 223 */
    double R_0 = Material_GetProperty(mat)[pm("R_CaOH2")] ; /* (dm) */
    double D   = Material_GetProperty(mat)[pm("D")] ; /* (mol/dm/s) */
    
    if(R_0 == 0.) R_0 = 40e-5 ;
    
    if(D == 0.) D = 7e-15 ;
    
    n_ch0 = Material_GetProperty(mat)[pm("N_CaOH2")] ; /* contenu molaire initial en CaOH2 */
    
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
    
    c_2 = h*R_0/D ;     /* (no dim) these MT p 228 */
    
  }
  
  Material_GetProperty(mat)[pm("A_2")] = a_2 ;
  Material_GetProperty(mat)[pm("C_2")] = c_2 ;
  
  /* Create the inverse of the C/S curve */
  {
    Curves_t* curves  = Material_GetCurves(mat) ;
    Curve_t*  curve   = Material_GetCurve(mat) ;
    Curve_t*  csratio = curve + 2 ;
    int saturation_ch = Curves_CreateInverse(curves,csratio,'n') ;
    
    if(saturation_ch != 6) {
      arret("ReadMatProp") ;
    }
    
    {
      int NbOfCurves = Material_GetNbOfCurves(mat) ;
      
      if(NbOfCurves > Material_MaxNbOfCurves) {
        arret("ReadMatProp") ;
      }
    }
    
    /*{
      char* filename = Curve_PrintInFile(curve + saturation_ch) ;
      Message_Direct("Curve printed in %s\n",filename) ;
    }*/
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
  printf("\t- Solid calcium carbonate content  (c_caco3)\n") ;
  printf("\t- Solid silicon content            (n_si_s)\n") ;
  
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
  double *f = Element_GetImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int i ;
  
  
  
  /*
    Donnees
  */
  GetProperties(el) ;
  
  {
    double x_na_tot 	= x_na0 ;
    double x_k_tot  	= x_k0 ;
    double x_co2   	  = C_CO2_eq ;
    double x_na = x_na_tot ;
    double x_k  = x_k_tot ;
    
    concentrations_oh_na_k(x_co2,&x_na,&x_k,el) ;


    for(i = 0 ; i < nn ; i++) {
      double s_ch       = SaturationDegreeOfCH(x_co2) ;
      double s_cheq     = s_ch ;
      double n_cc       = N_CC(i) ;
      double n_si_s     = N_Si_S(i) ;
      double v_csh      = MolarVolumeOfCSH(s_cheq) ;
      double v_s0       = V_CH*n_ch0 + V_CC*n_cc + v_csh*n_si_s ;
      
#if (U_CO2 == LOG_U)
      LogC_CO2(i) = log10(x_co2) ;
#else
      C_CO2(i)    = x_co2 ;
#endif
      C_Na(i)     = x_na ;
      C_K(i)      = x_k ;
      
      /* Solid contents */
      V_S0(i)    = v_s0 ;
      N_CH(i)    = n_ch0 ;
      S_CH_EQ(i) = s_ch ;
    }
  }
  

  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x       = ComputeComponents(el,u,u,f,0,i) ;
    
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
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
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
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int i ;
  
  
  /*
    Donnees
  */
  GetProperties(el) ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x      = ComputeComponents(el,u,u_n,f_n,dt,i) ;
    
    /* Back up */
    N_C(i)  = x[I_N_C] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Na(i) = x[I_N_Na] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ; 
    Mass(i) = x[I_Mass] ;
    N_q(i)  = x[I_N_Q] ;

    /* Liquid components */
    S_CH_EQ(i) = x[I_S_CH_EQ] ;
    
    /* Solid components */
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
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
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
  for(i=0;i<2*NEQ;i++){
    K(i,U_C_CO2)     *= Ln10*C_CO2(0) ;
    K(i,U_C_CO2+NEQ) *= Ln10*C_CO2(1) ;
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
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    nso = 48 ;
  int    i ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;


  if(Element_IsSubmanifold(el)) return(0) ;
  
   
  /*
    Donnees
  */
  GetProperties(el) ;
  

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double *x = ComputeComponents(el,u,u,f,0,j) ;
    
    /* pression */
    double p_l        =  x[I_P_L] ;
    /* saturation */
    double p_c        = p_g - p_l ;
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
    double n_si_s     = x[I_N_Si_S] ;
    double n_cc 	    = x[I_N_CC] ;
    double n_ch       = x[I_N_CH] ;
    double s_ch       = x[I_S_CH] ;
    double s_cheq     = x[I_S_CH_EQ] ;
    double x_csh      = CalciumSiliconRatioInCSH(s_cheq) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
    
    double s_sh       = x[I_S_SH] ;
  
    double v_csh      = MolarVolumeOfCSH(s_cheq) ;
    double v_solide_csh   = v_csh*n_si_s ;
    double v_solide_ch    = V_CH*n_ch ;
    double v_solide_cc    = V_CC*n_cc ;

    double n_chn      = x[I_N_CHn] ;
    /* double av = n_cc/n_ch0 ; */
    double av = 1 - n_chn/n_ch0 ;
    double dn1sdt = a_2*dn1_caoh2sdt(av,c_2) ;
    /* cinetique de dissolution de la portlandite en loi log(s_ch) */
    /* double dn_chsdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ; */
    double dn_chsdt = dn1sdt*log(s_ch) ;
    /* double coeff_dnCH = log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ; */
    double coeff_dnCH = log(s_ch) ;
  
    double CsurS      = (x_csh + n_ch/n_si_s) ;
    

    /* porosite */
    double phi        = x[I_Phi] ;

    double ph = 14 + log(x_oh)/log(10.) ;
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
    double *x = ComputeComponents(el,u,u,f,0,i) ;
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
    double *x1 = ComputeComponents(el,u,u,f,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,u,f,0.,0) ;
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
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
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
    double *x         = ComputeComponents(el,u,u_n,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    
    dxi[U_C_CO2 ] = x[U_C_CO2]*1.e-2 ;
    dxi[U_C_Na  ] = 1.e-4 ;
    dxi[U_C_K   ] = 1.e-4 ;
    dxi[U_N_CC  ] = x[U_N_CC  ]*1.e-2 ;
    dxi[U_N_Si_S] = -x[U_N_Si_S]*1.e-2 ;
    dxi[U_P_L   ] = P_Ln(i)*1.e-6 ;
    dxi[U_PSI   ] = 1. ;
    
    /*
    dxi[U_C_CO2 ] =  ObVal_GetValue(obval + U_C_CO2) ;
#if (U_CO2 == LOG_U)
    dxi[U_C_CO2 ] *=  x[U_C_CO2] ;
#endif
    dxi[U_C_Na  ] =  ObVal_GetValue(obval + U_C_Na) ;
    dxi[U_C_K   ] =  ObVal_GetValue(obval + U_C_K) ;
    dxi[U_N_CC  ] =  ObVal_GetValue(obval + U_N_CC) ;
    dxi[U_N_Si_S] = -ObVal_GetValue(obval + U_N_Si_S) ;
    dxi[U_P_L   ] =  ObVal_GetValue(obval + U_P_L) ;
    dxi[U_PSI   ] =  ObVal_GetValue(obval + U_PSI) ;
    */
    
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



double concentrations_oh_na_k(double x_co2,double *px_na,double *px_k,Element_t *el)
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
  
  double s_ch       = SaturationDegreeOfCH(x_co2) ;
  double s_cheq     = s_ch ;
  double s_sheq     = SaturationDegreeOfSHAtEq(s_cheq) ;
  double s_sh       = s_sheq ;
  double x_h4sio4   = IonActivityProductOfSH(s_sh) ;
  
  /* Initialization */
  double x_na = x_na_tot ;
  double x_k  = x_k_tot ;
  double x_oh0 = x_na + x_k ;
  double x_oh = x_oh0 ;
  
  /* x_na_tot =  A_Na*x_na + B_Na*x_na*x_oh + C_Na*x_na*x_oh^2 */
  double A_Na = 1 ;
  double B_Na = k_naoh/k_e + k_nahco3*k_h*x_co2/k_1 ;
  double C_Na = k_naco3*k_h*x_co2/(k_1*k_e) ;

  /* x_k_tot =  A_K*x_k + B_K*x_k*x_oh */
  double A_K = 1 ;
  double B_K = k_koh/k_e  ;
  
  double err,tol = 1.e-8 ;
  int i = 0 ;
  
  
  
  do {
    double dx_oh = - x_oh ;
    
    x_oh = concentration_oh(x_co2,x_na,x_k,x_h4sio4,el) ;
    
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
  
  /*
    Data
  */
  
  phii     = GetProperty("porosite") ;
  frac     = GetProperty("frac") ;
  phi_r    = GetProperty("phi_r") ;
  
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



double* ComputeComponents(Element_t *el,double **u,double **u_n,double *f_n,double dt,int n)
{
  double *v0 = Element_GetConstantTerm(el) ;
  double *x = Components ;
  
  /* Primary Variables */
  x[U_C_CO2 ] = C_CO2(n) ;
  x[U_C_Na  ] = C_Na(n) ;
  x[U_C_K   ] = C_K(n) ;
  x[U_N_CC  ] = N_CC(n) ;
  x[U_N_Si_S] = N_Si_S(n) ;
  x[U_P_L   ] = P_L(n) ;
  x[U_PSI   ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn   ]  = N_CHn(n) ;
  x[I_S_CH_EQn]  = S_CH_EQn(n) ;
  x[I_N_Si_Sn ]  = N_Si_Sn(n) ;
  x[I_V_S0    ]  = V_S0(n) ;
  
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  dx[U_C_CO2 ] = x[U_C_CO2] ;
  dx[U_C_Na  ] = x[U_C_Na ] ;
  dx[U_C_K   ] = x[U_C_K  ] ;
  dx[U_N_CC  ] = x[U_N_CC   ] ;
  dx[U_N_Si_S] = x[U_N_Si_S ] ;
  dx[U_P_L   ] = x[U_P_L  ] ;
  dx[U_PSI   ] = x[U_PSI  ] ;

  dx[I_N_CHn   ] = x[I_N_CHn] ;
  dx[I_S_CH_EQn] = x[I_S_CH_EQn] ;
  dx[I_N_Si_Sn ] = x[I_N_Si_Sn] ;
  dx[I_V_S0    ] = x[I_V_S0] ;
  
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
  double x_co2      = x[U_C_CO2 ] ;
  double x_na       = x[U_C_Na  ] ;
  double x_k        = x[U_C_K   ] ;
  double n_cc       = x[U_N_CC  ] ;
  double n_si_s     = x[U_N_Si_S] ;
  double p_l        = x[U_P_L   ] ;
  
  /* Liquid components */
  double s_ch       = SaturationDegreeOfCH(x_co2) ;
  double s_cheqn    = x[I_S_CH_EQn] ;
  double n_si_sn    = x[I_N_Si_Sn] ;
  double r_si_s     = n_si_s/n_si_sn ;
  double s_cheq     = SaturationDegreeOfCHAtEq(dt,s_ch,s_cheqn,r_si_s) ;
  
  double s_sh       = SaturationDegreeOfSH(dt,s_cheq,r_si_s) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  
  double x_h4sio4   = Q_SH ;
  
  double x_oh       = ConcentrationOfOHInLiquid(x_co2,x_na,x_k,x_h4sio4) ;
  double x_h        = k_e/x_oh ;
  
  double x_h2co3    = k_h*x_co2 ;
  double x_hco3     = x_oh*x_h2co3/k_1 ;
  double x_co3      = k_co3*x_oh*x_hco3 ;
  
  double x_ca       = k_ca/x_co3 ;
  double x_caoh     = k_caoh/k_e*x_oh*x_ca ;
  double x_cahco3   = k_cahco3*x_hco3*x_ca ;
  double x_caco3aq  = k_caco3*k_ca/(k_e*k_co3);
  double x_caoh2aq  = k_caoh2*x_ca*x_oh*x_oh ;
  
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
  /* Dissolution kinetics of portlandite: A*log(s_ch) */
  /* double dn_chsdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ; */
  double dn_chsdt = dn1sdt*log(s_ch) ;
  /* ... (power law) */ 
  /* double dn_chsdt = - dn1sdt*pow(fabs(1 - s_ch),Xp) ; */
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
  x[I_S_SH      ] = s_sh ;
  
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



double concentration_oh(double x_co2,double x_na,double x_k,double x_h4sio4,Element_t *el)
/* Solve for electroneutrality : SUM(z_i c_i) = 0
   as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* The primary variables are considered as constant */
  /* Ion activity products are constant as well */
  /* Other constant concentrations */
  double x_h2co3    = k_h*x_co2 ;
  /* Electroneutrality is written as  Sum z_i*x_i = 0, with
   * ------------------------------------------------------------------
   * x_i        =  A_i*(x_h)**n                              : n    z_i
   * ------------------------------------------------------------------
   * x_h                                                     : +1   +1
   * x_oh       = k_e/x_h ;                                  : -1   -1
   * 
   * x_h2co3    = k_h*x_co2 ;                                :  0    0
   * x_hco3     = k_e*x_h2co3/k_1/x_h ;                      : -1   -1
   * x_co3      = k_e*k_co3*x_hco3/x_h ;                     : -2   -2
   * 
   * x_ca       = k_ca/x_co3 ;                               : +2   +2
   * x_caoh     = k_caoh*x_ca/x_h ;                          : +1   +1
   * x_cahco3   = k_cahco3*x_hco3*x_ca ;                     : +1   +1
   * x_caco3aq  = k_caco3*k_ca/(k_e*k_co3);                  :  0    0
   * x_caoh2aq  = k_caoh2*k_e*k_e*x_ca/(x_h*x_h) ;           :  0    0
   * 
   * x_h4sio4   = IonActivityProductOfSH(s_sh) ;             :  0    0
   * x_h3sio4   = x_h4sio4/(k_h4sio4/k_h2sio4)/x_h ;         : -1   -1
   * x_h2sio4   = x_h3sio4/(k_h2sio4)/x_h ;                  : -2   -2
   * 
   * x_cah2sio4 = k_cah2sio4*x_h2sio4*x_ca ;                 :  0    0
   * x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca ;                 : +1   +1
   * 
   * x_na                                                    :  0   +1
   * x_naoh     = k_naoh*x_na/x_h ;                          : -1    0
   * x_nahco3   = k_nahco3*x_na*x_hco3 ;                     : -1    0
   * x_naco3    = k_naco3*x_na*x_hco3/x_h ;                  : -2   -1
   * 
   * x_k                                                     :  0   +1
   * x_koh      = k_koh*x_k/x_h ;                            : -1    0
   *
   * We compute the A_i for i having non zero z_i
   */
  double A_hco3     = x_h2co3/k_1*k_e ;
  double A_co3      = k_e*k_co3*A_hco3 ;
  double A_ca       = k_ca/A_co3 ;
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
      printf("x_h  = %e\n",x_h) ;
      printf("a  = %e\n",a) ;
      printf("b  = %e\n",b) ;
      printf("c  = %e\n",c) ;
      printf("d  = %e\n",d) ;
      printf("x_co2 = %e\n",x_co2) ;
      arret("concentration_oh : x_h<0") ;
  }
 
  return(k_e/x_h) ;
}


double poly4(double a,double b,double c,double d,double e)
/* on resout ax^4 + bx^3 + cx^2 + dx + e = 0 */
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
  
  x = Math_PolishPolynomialEquationRoot(y,4,x,tol*x,10) ;
  
  return(x) ;
}





double ComputeS_CH_EQ(Element_t *el,double dt,double s_ch,double s_cheqn,double r_si_s)
/* Solve for s the two following equations: 
 * (L) tau * [ ln(x) - ln(x0) ]  = dt * [ ln(s0) - ln(s) ]
 * (C) x = chi(s) 
 * for the two given values s0 = s_ch and x0 = xn/r_si_s.
 * (L) is a line of slope dt/tau in the plane [ln(s),ln(x)], 
 * passing through the point L0 = (ln(s0),ln(x0)).
 * The solution is the intersection of (L) with (C) in [ln(s),ln(x)]. 
 * It is noticed that the solution is invariant by any transformation 
 * of (s0,x0) defined by tau * ln(x0) + dt * ln(s0) = constant.
 * Return s.
 */
{
  double tau = tau_ca ;
  double xn = CalciumSiliconRatioInCSH(s_cheqn) ;
  double s0  = s_ch ;
  double x0  = xn/r_si_s ;
  
  
  /* (L) is a vertical line */
  if(tau == 0) {
    return(s0) ;
  }
  
  
  /* (L) is a horizontal line */
  if((dt == 0)) {
    double x_max = CalciumSiliconRatioInCSH(1) ;
    
    if(x0 > x_max) {
      Message_Warning("ComputeS_CH_EQ(1): x > x_max, theoretically the solution is inf") ;
      return(1) ;
    } else {
      double s_ini = s_cheqn ;
      double s = ComputeS_CH_EQ0(el,0,s_ini,x0) ;
      
      return(s) ;
    }
  }


  /* General case
   * We look for another point L1 = (ln(s1),ln(x1)) of (L) so that 
   * the signs of (x0 - chi(s0)) and (x1 - chi(s1)) are opposite. This is 
   * either (s1,x1) with s1 < 1  and  x1 = chi(s0) 
   * or     (s1,x1) with s1 = 1  and  x1 > chi(s0).
   * So we compute the coordinates (s1,x1) of L1 as follows.
   */
  {
    /* Point C0 = (s0,xs0) on the curve associated with s0 */
    double xs0  = CalciumSiliconRatioInCSH(s0) ;
    /* Point L1 = (s1,x1) with x1 = xs0 */
    double lns1 = log(s0) + tau/dt * log(x0/xs0) ;
    double s1   = exp(lns1) ;
    /* Point L1 = (s1,x1) with s1 = 1 */
    double lnx1 = log(x0) + dt/tau * log(s0) ;
    double x1   = exp(lnx1) ;
    

    /* Case L1 = (s1,x1) with s1 = 1 and x1 > xs0 */
    if(lns1 > 0) {
        
      s1   = 1 ;
      lns1 = 0 ;
        
    /* Case L1 = (s1,x1) with s1 < 1 and x1 = xs0 */
    } else {
        
      x1 = xs0 ;
      lnx1 = log(xs0) ;
      
      /* Check that lns1 is not too small and s1 > 0 */
      {
        double s1min = InverseOfCalciumSiliconRatioInCSH(0) ;
        double lns1min = log(s1min) ;
      
        if(lns1 < lns1min) {
          s1   = s1min ;
          lns1 = lns1min ;
          lnx1 = log(x0) + dt/tau * log(s0/s1) ;
          x1   = exp(lnx1) ;
        }
      }
    }
    
    
    /* Check if the solution is already found */
    {
      double xs1  = CalciumSiliconRatioInCSH(s1) ;
      double tol = 1.e-10 ;
      double err ;
      double s ;
      
      if(fabs(x0 -xs0) < fabs(x1 - xs1)) {
        err = fabs(x0 - xs0) ;
        s = s0 ;
      } else {
        err = fabs(x1 - xs1) ;
        s = s1 ;
      }
      
      if(err < tol) return(s) ;
    }
    
    
    /* Case A: L0 is above (C) ie x0/xs0 > 1 (and s1 > s0) */
    if(x0 > xs0) {
      /* If s0 >= 1 the solution should be 
       * ln(s) = ln(s0) + tau/dt * ln(x0/xs0) no matter if s > 1 
       * But to avoid a solution = inf we take s0 */
      if(s0 >= 1) {
        /*
        Message_Warning("ComputeS_CH_EQ(2): s0 = %e",s0) ;
        double lns = log(s0) + tau/dt * log(x0/xs0) ;
        double s   = exp(lns) ;
        return(s) ;
        */
        return(s0) ;
      }

      /* If x1 >= x_max we get the solution */
      {
        double x_max = CalciumSiliconRatioInCSH(1) ;
      
        /* It doesn't matter if s > 1 eventually
         * But to avoid a solution = inf we take 1 */
        if(x1 >= x_max) {
          /*
          double lns = log(s0) + tau/dt * log(x1/x_max) ;
          double s   = exp(lns) ;
          */
      
          return(1) ;
        }
      }
      
      /* Check if L1 is below (C).  */
      {
        /* Point C1 = (s1,xs1) on the curve associated with s1 */
        double xs1  = CalciumSiliconRatioInCSH(s1) ;
        
        if(x1 > xs1) {
          arret("ComputeS_CH_EQ(3): no opposite signs") ;
        }
      }
      
      
    /* Case B: L0 is below (C) ie x0/xs0 < 1 (and s1 < s0).
     * We find a point L1 above (C)  */
    } else {
      /* L1 = (s1,x1) with x1 = xs0 and s1 < s0.
       * If s1 >= 1 this is the solution, no matter if s1 > 1 */
      if(lns1 >= 0) {
        return(s1) ;
      }
      
      /* Check if L1 is above (C).  */
      {
        /* Point C1 = (s1,xs1) on the curve associated with s1 */
        double xs1  = CalciumSiliconRatioInCSH(s1) ;
        
        if(x1 < xs1) {
          arret("ComputeS_CH_EQ(4): no opposite signs") ;
        }
      }
      
    }
    
    /* Test of signs */
    {
      double xs1 = CalciumSiliconRatioInCSH(s1) ;
      
      if(((x0 - xs0)*(x1 - xs1)) > 0) {
        Message_Warning("ComputeS_CH_EQ(5): no opposite signs") ;
        
        if(fabs(x0 - xs0) < fabs(x1 - xs1)) {
          double s = ComputeS_CH_EQ0(el,dt/tau,s0,x0) ;
          
          return(s) ;
        } else {
          double s = ComputeS_CH_EQ0(el,dt/tau,s1,x1) ;
          
          return(s) ;
        }
      
      }
    }
    
    {
      double s = ComputeS_CH_EQ2(el,s0,x0,s1,x1) ;
    
      return(s) ;
    }
  }
  
}





double ComputeS_SH(Element_t *el,double dt,double s_cheq,double r_si_s)
/* Solve for s_sh, the following equation:
 * tau * ln(r_si_s) = dt * [ ln(s_sh) - ln(s_sheq) ]
 * with given r_si_s and s_sheq.
 */
{
  double tau = tau_si ;
  double s_sheq = SaturationDegreeOfSHAtEq(s_cheq) ;
  
  if(tau == 0 || r_si_s == 1) {
    return(s_sheq) ;
  }
  
  if(dt == 0) {
    Message_Warning("ComputeS_SH: abnormal situation") ;
    
    return(s_sheq) ;
  }
  
  {
    double lns_sheq = log(s_sheq) ;
    double lnr_si_s = log(r_si_s) ;
    double lns_sh = lns_sheq + tau/dt*lnr_si_s ;
    double s_sh = exp(lns_sh) ;
    
    /*
    if(lnr_si_s > 0 && lns_sh > 0) {
      Message_Warning("ComputeS_SH: lnr_si_s = %e lns_sh = %e > 0",lnr_si_s,lns_sh) ;
    }
    */
    
    return(s_sh) ;
  }
  
}



double ComputeS_CH_EQ0(Element_t *el,double b,double s0,double x0)
/* Solve for s, the two following equations: 
 * (L) ln(x) - ln(x0)  = b * [ ln(s0) - ln(s) ]
 * (C) x = chi(s) 
 * (L) is a line in the plane  [ln(s),ln(x)].
 * (C) is a curve
 * L0 = (ln(s0),ln(x0)) is a given arbitrary point of (L).
 * The solution is the intersection of (L) with (C) in [ln(s),ln(x)].
 * Newton method.
 * Return s.
 */
{
  
  {
    double lns0 = log(s0) ;
    double lnx0 = log(x0) ;
    double err,tol = 1.e-4 ;
    double s = s0 ;
    double lns = log(s0) ;
    int    i = 0 ;
    
    do {
      double x    = CalciumSiliconRatioInCSH(s) ;
      double lnx  = log(x) ;
      double dlnx = dCalciumSiliconRatioInCSH(s)*s/x ;
      double f  = lnx  - lnx0 + b*(lns  - lns0) ;
      double df = dlnx + b ; 
      double dlns = (df != 0) ? - f/df : 0 ;
      
      lns += dlns ;
      s  = exp(lns) ;
      
      err = fabs(dlns) ;
      
      if(i++ > 50) {
        double y0  = CalciumSiliconRatioInCSH(s0) ;
        
        printf("\n") ;
        printf("s0 = %e\n",s0) ;
        printf("x0  = %e\n",x0) ;
        printf("chi(s0) = %e\n",y0) ;
        printf("x    = %e\n",x) ;
        printf("s    = %e\n",s) ;
        printf("err  = %e\n",err) ;
        printf("f    = %e\n",f) ;
        printf("df   = %e\n",df) ;
        printf("dlns = %e\n",dlns) ;
        arret("ComputeS_CH_EQ0: no convergence") ;
      }
    } while(err > tol) ;
    
    return(s) ;
  }
  
}



double ComputeS_CH_EQ1(Element_t *el,double s0,double x0,double s1,double x1)
/* Solve for s, the two following equations: 
 * (C) x = chi(s) 
 * (L) ln(x) - ln(x0) = b * [ ln(s) - ln(s0) ]
 * with b given by ln(x1) - ln(x0) = b * [ ln(s1) - ln(s0) ].
 * (L) is the line joining the 2 points L0 and L1 of coordinates
 * (s0,x0) and (s1,x1) in the plane [ln(s),ln(x)].
 * (C) is a assumed as given by an monotoneous increasing function chi(s).
 * The solution is the intersection of (L) with (C) in [ln(s),ln(x)]. 
 * We assume that the solution lies between s0 and s1.
 * Return s.
 */
{
  double b = log(x1/x0)/log(s1/s0) ;

  
  {
    double s0ini = s0 ;
    double x0ini = x0 ;
    double s1ini = s1 ;
    double x1ini = x1 ;
    double err,tol = 1.e-8 ;
    double s ;
    int    i = 0 ;
    
    /* Li is a point (si,xi) of the line (L) in the plane [ln(s),ln(x)]
     * Ci is a point (si,xsi) of the curve (C) in the plane [ln(s),ln(x)]
     */
    do {
      
      
      /* The two points C0 and C1 on the curve associated with s0 and s1 */
      double xs0    = CalciumSiliconRatioInCSH(s0) ;
      double xs1    = CalciumSiliconRatioInCSH(s1) ;
      /* The slope of the line joining C0 and C1 */
      double a01    = log(xs1/xs0)/log(s1/s0) ;
      /* L2 = intersection of the lines of slope a01 and (L) */
      double lns2 = log(s0) + log(xs0/x0) / (b - a01) ;
      double lnx2 = log(x0) + b * (lns2 - log(s0)) ;
      double s2   = exp(lns2) ;
      double x2   = exp(lnx2) ;
      /* C2 = point on the curve at s2 */
      double xs2  = CalciumSiliconRatioInCSH(s2) ;
      /* The slope of the tangent at C2 */
      double t2   = dCalciumSiliconRatioInCSH(s2)*s2/xs2 ;
      /* L3 = intersection of the tangent t2 and the line (L) */
      double lns3 = log(s2) + log(xs2/x2) / (b - t2) ;
      double lnx3 = log(x2) + b * (lns3 - log(s2)) ;
      double s3   = exp(lns3) ;
      double x3   = exp(lnx3) ;
      /* C3 = point on the curve at s3 */
      double xs3  = CalciumSiliconRatioInCSH(s3) ;
      
      
      if(fabs(x0 - xs0) < fabs(x1 - xs1)) {
        err = fabs(x0 - xs0) ;
        s = s0 ;
      } else {
        err = fabs(x1 - xs1) ;
        s = s1 ;
      }
      
      if(err < tol) return(s) ;
      
      
      /* Test the signs */
      if(((x2 - xs2)*(x3 - xs3)) > 0) {
        Message_Warning("ComputeS_CH_EQ1: no opposite signs") ;
        if(fabs(x2 - xs2) < fabs(x3 - xs3)) {
          return(ComputeS_CH_EQ0(el,-b,s2,x2)) ;
        } else {
          return(ComputeS_CH_EQ0(el,-b,s3,x3)) ;
        }
      }
      
      /* We actualize L0 as L2 or L3 */
      if(((x0 - xs0)*(x2 - xs2)) >= 0) {
        s0 = s2 ;
        x0 = x2 ;
      } else if(((x0 - xs0)*(x3 - xs3)) >= 0)  {
        s0 = s3 ;
        x0 = x3 ;
      } else {
        arret("ComputeS_CH_EQ1: something's wrong") ;
      }
      
      /* We actualize L1 as l2 or L3 */
      if(((x1 - xs1)*(x3 - xs3)) >= 0) {
        s1 = s3 ;
        x1 = x3 ;
      } else if(((x1 - xs1)*(x2 - xs2)) >= 0)  {
        s1 = s2 ;
        x1 = x2 ;
      } else {
        arret("ComputeS_CH_EQ1: something's wrong") ;
      }
      
      
    } while(i++ < 50) ;
    
    if(1) {
      printf("\n") ;
      printf("s0ini = %e\n",s0ini) ;
      printf("x0ini  = %e\n",x0ini) ;
      printf("s1ini = %e\n",s1ini) ;
      printf("x1ini  = %e\n",x1ini) ;
      printf("s0 = %e\n",s0) ;
      printf("x0  = %e\n",x0) ;
      printf("s1 = %e\n",s1) ;
      printf("x1  = %e\n",x1) ;
      printf("err  = %e\n",err) ;
      arret("ComputeS_CH_EQ1: no convergence") ;
    }
    
    return(s) ;
  }
  
}



double ComputeS_CH_EQ2(Element_t *el,double s0,double x0,double s1,double x1)
/* Solve for s, the two following equations: 
 * (C) x = chi(s) 
 * (L) ln(x) - ln(x0) = b * [ ln(s) - ln(s0) ]
 * with b given by ln(x1) - ln(x0) = b * [ ln(s1) - ln(s0) ].
 * (L) is the line joining the 2 points L0 and L1 of coordinates
 * (ln(s0),ln(x0)) and (ln(s1),ln(x1)) in the plane [ln(s),ln(x)].
 * (C) is assumed as given by a monotoneous increasing function chi(s).
 * The solution is the intersection of (L) with (C) in [ln(s),ln(x)]. 
 * We assume that the solution lies between s0 and s1.
 * Dichotomy method.
 * Return s.
 */
{
  
  {
    double s0ini = s0 ;
    double x0ini = x0 ;
    double s1ini = s1 ;
    double x1ini = x1 ;
    double err,tol = 1.e-8 ;
    double s ;
    int    i = 0 ;
    
    /* Li is a point (si,xi) of the line (L) in the plane [ln(s),ln(x)]
     * Ci is a point (si,xsi) of the curve (C) in the plane [ln(s),ln(x)]
     */
    do {
      /* The two points C0 and C1 on the curve associated with s0 and s1 */
      double xs0    = CalciumSiliconRatioInCSH(s0) ;
      double xs1    = CalciumSiliconRatioInCSH(s1) ;
      /* L2 = Point on (L) between L0 and L1 obtained by dichotomy */
      double a = 0.5 ;
      double lns2 = a*log(s1) + (1 - a)*log(s0) ;
      double lnx2 = a*log(x1) + (1 - a)*log(x0) ;
      double s2   = exp(lns2) ;
      double x2   = exp(lnx2) ;
      /* C2 = point on the curve at s2 */
      double xs2  = CalciumSiliconRatioInCSH(s2) ;
      
      if(fabs(x0 - xs0) < fabs(x1 - xs1)) {
        err = fabs(x0 - xs0) ;
        s = s0 ;
      } else {
        err = fabs(x1 - xs1) ;
        s = s1 ;
      }
      
      if(err < tol) return(s) ;
      
      /* We actualize L0 or L1 as L2 */
      if(((x0 - xs0)*(x2 - xs2)) >= 0) {
        s0 = s2 ;
        x0 = x2 ;
      } else  {
        s1 = s2 ;
        x1 = x2 ;
      }
      
      
    } while(i++ < 50) ;
    
    if(1) {
      printf("\n") ;
      printf("s0ini = %e\n",s0ini) ;
      printf("x0ini  = %e\n",x0ini) ;
      printf("s1ini = %e\n",s1ini) ;
      printf("x1ini  = %e\n",x1ini) ;
      printf("s0 = %e\n",s0) ;
      printf("x0  = %e\n",x0) ;
      printf("s1 = %e\n",s1) ;
      printf("x1  = %e\n",x1) ;
      printf("err  = %e\n",err) ;
      arret("ComputeS_CH_EQ2: no convergence") ;
    }
    
    return(s) ;
  }
  
}






double ComputeS_CH_EQ3(Element_t *el,double dt,double s_ch,double s_cheqn,double r_si_s)
/* Solve for s the two following equations: 
 * (L) tau * [ ln(x) - ln(x0) ]  = dt * [ ln(s0) - ln(s) ]
 * (C) x = chi(s) 
 * for the two given values s0 = s_ch and x0 = xn/r_si_s.
 * (L) is a line of slope dt/tau in the plane [ln(s),ln(x)], 
 * passing through the point L0 = (ln(s0),ln(x0)).
 * The solution is the intersection of (L) with (C) in [ln(s),ln(x)]. 
 * It is noticed that the solution is invariant by any transformation 
 * of (s0,x0) defined by tau * ln(x0) + dt * ln(s0) = constant.
 * Return s.
 */
{
  double tau = tau_ca ;
  double xn = CalciumSiliconRatioInCSH(s_cheqn) ;
  double s0  = s_ch ;
  double x0  = xn/r_si_s ;
  double s0ini = s0 ;
  double x0ini = x0 ;
  
  
  /* (L) is a vertical line */
  if(tau == 0) {
    return(s0) ;
  }
  
  
  /* (L) is a horizontal line */
  if((dt == 0)) {
    double x_max = CalciumSiliconRatioInCSH(1) ;
    
    if(x0 > x_max) {
      Message_Warning("ComputeS_CH_EQ(1): x > x_max, theoretically the solution is inf") ;
      return(1) ;
    } else {
      double s = InverseOfCalciumSiliconRatioInCSH(x0) ;
      
      return(s) ;
    }
  }


  /* General case
   * We look for another point L1 = (ln(s1),ln(x1)) of (L)
   * so that (s1,x1) with s1 = 1.
   * So we compute the coordinates (s1,x1) of L1 as follows.
   */
  {
    /* Point L1 = (s1,x1) with s1 = 1 */
    double lnx1 = log(x0) + dt/tau * log(s0) ;
    double x1   = exp(lnx1) ;
    double s1   = 1 ;
    
    
    /* Check if the solution is already found */
    {
      double xs0  = CalciumSiliconRatioInCSH(s0) ;
      double xs1  = CalciumSiliconRatioInCSH(s1) ;
      double tol = 1.e-10 ;
      double err ;
      double s ;
      
      if(fabs(x0 -xs0) < fabs(x1 - xs1)) {
        err = fabs(x0 - xs0) ;
        s = s0 ;
      } else {
        err = fabs(x1 - xs1) ;
        s = s1 ;
      }
      
      if(err < tol) return(s) ;
    }
    
    
    /* If x1 >= x_max we get the solution */
    {
      double x_max = CalciumSiliconRatioInCSH(1) ;
      
      /* It doesn't matter if s > 1 eventually
        * But to avoid a solution = inf we take 1 */
      if(x1 >= x_max) {
        /*
        double lns = log(s0) + tau/dt * log(x1/x_max) ;
        double s   = exp(lns) ;
        */
      
        return(1) ;
      }
    }
    
    s0 = s1 ;
    x0 = x1 ;
  }
  
  {
    double a = dt/tau ;
    double err ;
    double tol = 1.e-8 ;
    int j = 0 ;
    
    do {
      double s1 = InverseOfCalciumSiliconRatioInCSH(x0) ;
      double x1 = CalciumSiliconRatioInCSH(s1) ;
      double s2 = s0 ;
      double x2 = CalciumSiliconRatioInCSH(s2) ;
      
      err = fabs(x1 - x2) ;
      
      if(fabs(s1 - s2) > 0) {
        double b = (log(x2) - log(x1))/(log(s2) - log(s1)) ;
        double lns = log(s0) + log(x0/x2)/(b + a) ;
        double lnx = log(x0) + log(x2/x0)*a/(b + a) ;
        
        s0 = exp(lns) ;
        x0 = exp(lnx) ;
      } else {
        double s = ComputeS_CH_EQ0(el,a,s2,x2) ;
          
        return(s) ;
      }
      
      if(j++ > 50) {
        printf("\n") ;
        printf("s0ini = %e\n",s0ini) ;
        printf("x0ini = %e\n",x0ini) ;
        printf("s0 = %e\n",s0) ;
        printf("x0 = %e\n",x0) ;
        printf("s1 = %e\n",s1) ;
        printf("x1 = %e\n",x1) ;
        printf("s2 = %e\n",s2) ;
        printf("x2 = %e\n",x2) ;
        printf("err  = %e\n",err) ;
        arret("ComputeS_CH_EQ3") ;
      }
      
    } while(err > tol) ;
    
    return(s0) ;
  }
  
}
