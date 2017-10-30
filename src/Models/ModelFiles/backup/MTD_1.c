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

#define TITLE   "Carbonation Of CBM (2013)"
#define AUTHORS "Morandeau-Thiery-Dangla"

#include "PredefinedMethods.h"

/* Macros */

#define NEQ    	(8)
#define NVE    	(58)
#define NVI     (30)
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

#define N_Cn(n)       (f_n[(n)])
#define N_qn(n)       (f_n[(2+n)])
#define Mass_n(n)     (f_n[(4+n)])
#define N_Can(n)      (f_n[(6+n)])
#define N_Nan(n)      (f_n[(10+n)])
#define N_Kn(n)       (f_n[(12+n)])
#define N_Sin(n)      (f_n[(16+n)])
#define N_CHn(n)      (f_n[(28+n)])

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
#include "ElectricChargeOfIonsInWater.h"
#define z_h           ElectricChargeOfIonsInWater(H)
#define z_oh          ElectricChargeOfIonsInWater(OH)
#define z_hco3        ElectricChargeOfIonsInWater(HCO3)
#define z_co3         ElectricChargeOfIonsInWater(CO3)
#define z_ca          ElectricChargeOfIonsInWater(Ca)
#define z_caoh        ElectricChargeOfIonsInWater(CaOH)
#define z_cahco3      ElectricChargeOfIonsInWater(CaHCO3)
#define z_h3sio4      ElectricChargeOfIonsInWater(H3SiO4)
#define z_h2sio4      ElectricChargeOfIonsInWater(H2SiO4)
#define z_cahco3      ElectricChargeOfIonsInWater(CaHCO3)
#define z_cah3sio4    ElectricChargeOfIonsInWater(CaH3SiO4)
#define z_na          ElectricChargeOfIonsInWater(Na)
#define z_naco3       ElectricChargeOfIonsInWater(NaCO3)
#define z_k           ElectricChargeOfIonsInWater(K)


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


/* Molar Masses */
#include "MolarMassOfMolecules.h"
#define gr             (1.e3)
#define M_Ca           (MolarMassOfMolecules(Ca)*gr)
#define M_H2CO3        (MolarMassOfMolecules(H2CO3)*gr)
#define M_HCO3         (MolarMassOfMolecules(HCO3)*gr)
#define M_CO3          (MolarMassOfMolecules(CO3)*gr)
#define M_OH           (MolarMassOfMolecules(OH)*gr)
#define M_H            (MolarMassOfMolecules(H)*gr)
#define M_H2O          (MolarMassOfMolecules(H2O)*gr)
#define M_Na           (MolarMassOfMolecules(Na)*gr)
#define M_NaOH         (MolarMassOfMolecules(NaOH)*gr)
#define M_NaHCO3       (MolarMassOfMolecules(NaHCO3)*gr)
#define M_NaCO3        (MolarMassOfMolecules(NaCO3)*gr)
#define M_CO2          (MolarMassOfMolecules(CO2)*gr)
#define M_CaOH2        (MolarMassOfMolecules(CaOH2)*gr)
#define M_CaCO3        (MolarMassOfMolecules(CaCO3)*gr)
#define M_CaO          (MolarMassOfMolecules(CaO)*gr)
#define M_SiO2         (MolarMassOfMolecules(SiO2)*gr)
#define M_K            (MolarMassOfMolecules(K)*gr)
#define M_KOH          (MolarMassOfMolecules(KOH)*gr)
#define M_CaOH         (MolarMassOfMolecules(CaOH)*gr)
#define M_CaHCO3       (MolarMassOfMolecules(CaHCO3)*gr)
#define M_CaCO3aq      (MolarMassOfMolecules(CaCO3)*gr)
#define M_CaOH2aq      (MolarMassOfMolecules(CaOH2)*gr)
#define M_CaH2SiO4     (MolarMassOfMolecules(CaH2SiO4)*gr)
#define M_CaH3SiO4     (MolarMassOfMolecules(CaH3SiO4)*gr)
#define M_H3SiO4       (MolarMassOfMolecules(H3SiO4)*gr)
#define M_H2SiO4       (MolarMassOfMolecules(H2SiO4)*gr)
#define M_H4SiO4       (MolarMassOfMolecules(H4SiO4)*gr)

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
#define k_sil      		(1.936421964e-3)        /* Solubilité de la silice amorphe */
#define k_cah2sio4 	  (39810.71)              /* Ca <-> CaH2SiO4*/
#define k_cah3sio4 	  (15.84)                 /* Ca <-> CaH3SiO4*/


/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp   	   	 (0.5)

/*Loi de Davies pour la prise en compte de l'activite ionique*/
#define A_DAVIES  (0.5)
#define B_DAVIES  (0.24)


/* Material Properties */
#define SaturationDegree(p)              (Curve_ComputeValue(Element_GetCurve(el),p))
#define RelativePermeabilityToLiquid(p)  (Curve_ComputeValue(Element_GetCurve(el) + 1,p))


/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(q)      (Curve_ComputeValue(Element_GetCurve(el) + 2,q))
#define WaterSiliconRatioInCSH(q)        (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define MolarVolumeOfCSH(q)              (Curve_ComputeValue(Element_GetCurve(el) + 4,q))


/* CH Properties */
#define C_CO2_eq                         (k_ca*k_1/(k_h*k_co3*k_2))
#define SaturationDegreeOfCH(c_co2)      (C_CO2_eq/c_co2)


/* S-H Properties */
/* Saturation degree of dissolved S-H */
#define S_SHeq(q)                        (Curve_ComputeValue(Element_GetCurve(el) + 5,q))
#define SaturationDegreeOfSH(q)          S_SHeq(q)
/* Ion Activity Product of dissolved S-H */
#define IonActivityProductOfSH(q)        (k_sil*SaturationDegreeOfSH(q))


/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(char *s) ;

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

static double concentration_oh(double,double,double,Element_t*) ;
static double concentrations_oh_na_k(double,double*,double*,Element_t*) ;
static double poly4(double,double,double,double,double) ;


/* Parametres */
static double phii,k_int,a_1,a_2,c_2,n_ch0,x_na0,x_k0,frac,phi_r ;
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
static double Faraday ;
static double RT ;


#include "WaterViscosity.h"
#include "DiffusionCoefficientOfMoleculesInWater.h"
#include "EquilibriumConstantOfHomogeneousReactionsInWater.h"
#include "DissociationConstantOfCementHydrationProducts.h"
#include "PhysicalConstants.h"
#define dm2            (1.e2)


void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (dm2/s) */
  d_oh         = DiffusionCoefficientOfMoleculesInWater(OH,TK)*dm2 ;
  d_h          = DiffusionCoefficientOfMoleculesInWater(H,TK)*dm2 ;
  d_hco3       = DiffusionCoefficientOfMoleculesInWater(HCO3,TK)*dm2 ;
  d_h2co3      = DiffusionCoefficientOfMoleculesInWater(H2CO3,TK)*dm2 ;
  d_co3        = DiffusionCoefficientOfMoleculesInWater(CO3,TK)*dm2 ;
  
  d_ca         = DiffusionCoefficientOfMoleculesInWater(Ca,TK)*dm2 ;
  d_caoh       = DiffusionCoefficientOfMoleculesInWater(CaOH,TK)*dm2 ;
  d_cahco3     = DiffusionCoefficientOfMoleculesInWater(CaHCO3,TK)*dm2 ;
  d_caco3aq    = DiffusionCoefficientOfMoleculesInWater(CaCO3,TK)*dm2 ;
  d_caoh2aq  	 = 7.92e-8 ;
  
  d_h4sio4     = DiffusionCoefficientOfMoleculesInWater(H4SiO4,TK)*dm2 ;
  d_h3sio4     = DiffusionCoefficientOfMoleculesInWater(H3SiO4,TK)*dm2 ;
  d_h2sio4     = DiffusionCoefficientOfMoleculesInWater(H2SiO4,TK)*dm2 ;
  
  d_cah2sio4   = DiffusionCoefficientOfMoleculesInWater(CaH2SiO4,TK)*dm2 ;
  d_cah3sio4   = DiffusionCoefficientOfMoleculesInWater(CaH3SiO4,TK)*dm2 ;
  
  d_na         = DiffusionCoefficientOfMoleculesInWater(Na,TK)*dm2*100 ;  /* (1.33e-7) */
  d_naoh       = 1.33e-5 ;  /* (1.33e-7) */
  d_nahco3     = 1.33e-5 ;  /* (1.33e-7) */
  d_naco3      = 1.33e-5 ;  /* (1.33e-7) */
  d_k          = DiffusionCoefficientOfMoleculesInWater(K,TK)*dm2*100 ; /* (1.957e-7) */
  d_koh        = DiffusionCoefficientOfMoleculesInWater(KOH,TK)*dm2*100 ; /* (1.957e-7) */

  /* Diffusion Coefficient Of Molecules In Air (dm2/s) */
  d_co2      	 = 1.6e-3 ;
  
  /* Viscosity (Pa.s) */
  mu_l       = WaterViscosity(TK) ;
  
  /* Physical constants */
  {
    RT      = PhysicalConstants(PerfectGasConstant)*TK*1.e3 ;
    Faraday = PhysicalConstants(Faraday)*1.e3 ;
  }
}

#define NbOfComponents    (50)
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
  else if(strcmp(s,"Tau") == 0)     return (12) ;
  else if(strcmp(s,"frac") == 0)    return (13) ;
  else if(strcmp(s,"phi_r") == 0)   return (14) ;
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
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}

int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 15 ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Initialisation automatique */
  {
    double h   = 5.6e-6 ;  /* (mol/dm2/s) these MT p 223 */
    double R_0 = Material_GetProperty(mat)[pm("R_CaOH2")] ;
    double D   = Material_GetProperty(mat)[pm("D")] ; /* (mol/dm/s) */
    double t_ch = R_0/(3*h*V_CH) ; /* (s) approx 7.215 s */
    n_ch0 = Material_GetProperty(mat)[pm("N_CaOH2")] ; /* contenu molaire initial en CaOH2 */
    /* a_2 = 3*h/R_0*n_ch0*V_CH ; */ /* (mol/dm3/s) these MT p 227 */
    a_2 = n_ch0/t_ch ;  /* (mol/dm3/s) these MT p 227 */
    c_2 = h*R_0/D ;     /* (no dim) these MT p 228 */
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
  double *u[MAX_NOEUDS] ;
  int i ;
  
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
  x_na0    = GetProperty("X_Na") ;
  x_k0     = GetProperty("X_K") ;
  
  {
    double x_na_tot 	= x_na0 ;
    double x_k_tot  	= x_k0 ;
    double x_co2   	  = C_CO2_eq ;
    double x_na = x_na_tot ;
    double x_k  = x_k_tot ;
    double x_oh = concentrations_oh_na_k(x_co2,&x_na,&x_k,el) ;
    
    /*
    {
      double x_h3sio4 ;
      concentration1(x_co2,x_na_tot,x_k_tot,&x_oh,&x_h3sio4,&x_na,&x_k) ;
    }
    */

    for(i = 0 ; i < 2 ; i++) {
      double s_ch       = SaturationDegreeOfCH(x_co2) ;
      double n_cc       = N_CC(i) ;
      double n_si_s     = N_Si_S(i) ;
      double v_csh      = MolarVolumeOfCSH(s_ch) ;
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
    }
  }
  

  for(i = 0 ; i < 2 ; i++) {
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
  int    nso = 48 ;
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
  frac     = GetProperty("frac") ;
  phi_r    = GetProperty("phi_r") ;

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
    double s_ch       = x[I_S_CH] ;
    double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
    
    double s_sh       = SaturationDegreeOfSH(s_ch) ;
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
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
  frac     = GetProperty("frac") ;
  phi_r    = GetProperty("phi_r") ;

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
    
    Kpsi_Ca   	 += Faraday/RT*KF_Ca*z_ca*x_ca ;
    Kpsi_HCO3 	 += Faraday/RT*KF_HCO3*z_hco3*x_hco3 ;
    Kpsi_CO3  	 += Faraday/RT*KF_CO3*z_co3*x_co3 ;
    Kpsi_OH   	 += Faraday/RT*KF_OH*z_oh*x_oh ;
    Kpsi_H    	 += Faraday/RT*KF_H*z_h*x_h ;
    Kpsi_Na    	 += Faraday/RT*KF_Na*z_na*x_na ;
    Kpsi_NaCO3 	 += Faraday/RT*KF_NaCO3*z_naco3*x_naco3 ;
    Kpsi_K     	 += Faraday/RT*KF_K*z_k*x_k ;
    Kpsi_CaOH 	 += Faraday/RT*KF_CaOH*z_caoh*x_caoh ;
    Kpsi_CaHCO3	 += Faraday/RT*KF_CaHCO3*z_cahco3*x_cahco3 ;
    Kpsi_H3SiO4	 += Faraday/RT*KF_H3SiO4*z_h3sio4*x_h3sio4 ;
    Kpsi_H2SiO4	 += Faraday/RT*KF_H2SiO4*z_h2sio4*x_h2sio4 ;
    Kpsi_CaH3SiO4+= Faraday/RT*KF_CaH3SiO4*z_cah3sio4*x_cah3sio4 ;
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
    
    x_oh = concentration_oh(x_co2,x_na,x_k,el) ;
    
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
    /* double frac  = 0.8 ; fractionnal length of pore bodies */
    /* double phi_r = 0.9 ;  phi_r = phi/phi0 original porosity for which permeability = 0 */

	  double S_s =  (phii - phi)/phii    ; /*saturation en solide */
	  double w = 1 + (1/frac)/(1/phi_r - 1) ;
    double t = (1 - S_s - phi_r)/(1 - phi_r) ;
	  double verma_pruess = (t > 0) ? t*t*(1 - frac + (frac/(w*w)))/(1 - frac + frac*(pow(t/(t + w - 1),2.))) : 0 ;
	
    /* permeability coefficient */
    coeff_permeability = verma_pruess ; /* Change here */
  }
  
  return(coeff_permeability) ;
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
  double x_h4sio4   = IonActivityProductOfSH(s_ch) ;
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
  double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double z_csh      = WaterSiliconRatioInCSH(s_ch) ;
  /* double av = n_cc/n_ch0 ; */
  double av = 1 - n_chn/n_ch0 ;
  double dn1sdt = a_2*dn1_caoh2sdt(av,c_2) ;
  /* Dissolution kinetics of portlandite: A log(s_ch) */
  /* double dn_chsdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ; */
  double dn_chsdt = dn1sdt*log(s_ch) ;
  /* ... (power law) */ 
  /* double dn_chsdt = - dn1sdt*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),Xp)) ; */
  double n_ch    = MAX(n_chn + dt*dn_chsdt , 0.) ;
  /* ... as elements: C, Ca, Si */
  double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  double n_c_s      = n_cc ;
  /* ... as mass */
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



double concentration_oh(double x_co2,double x_na,double x_k,Element_t *el)
/* Solve for electroneutrality : SUM(z_i c_i) = 0
   as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* The primary variables are considered as constant */
  /* Ion activity products are constant as well */
  double s_ch       = SaturationDegreeOfCH(x_co2) ;
  double Q_SH       = IonActivityProductOfSH(s_ch) ;
  /* Other constant concentrations */
  double x_h4sio4   = Q_SH ;
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
   * x_h4sio4   = IonActivityProductOfSH(s_ch) ;             :  0    0
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
  
  x = Math_PolishPolynomialEquationRoot(y,4,x,tol*x,20) ;
  
  return(x) ;
}
