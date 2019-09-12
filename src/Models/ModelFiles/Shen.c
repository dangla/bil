/* General features of the model:
 * Curves for CSH:
 *   - C/S ratio
 *   - H/S ratio
 *   - Molar Volume
 * CO2-H2O gas mixture as a supercritical gas
 * Temperature depencies for:
 *   - Diffusion coefficients
 *   - Chemical equilibrium constants
 *   - Transport parameters
 *   - Liquid and Gas viscosities
 * Alkalis and Chloride
 * Dissolution/Precipitation of CH, CSH and CC
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the numerical method */
#include "FVM.h"

#define TITLE   "Carbonation of CBM with scCO2 (2012)"
#define AUTHORS "Shen"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ     (7)
#define NVE     (51)
#define NVI     (28)
#define NV0     (2)

/* Equation Indexes */
#define E_C     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_Si    (3)
#define E_K     (4)
#define E_Cl    (5)
#define E_mass  (6)

/* Primary Unknown Indexes */
#define U_P_G      (0)
#define U_PSI      (1)
#define U_ZN_Ca_S  (2)
#define U_ZN_Si_S  (3)
#define U_C_K      (4)
#define U_C_Cl     (5)
#define U_P_L      (6)

#define NOLOG_U  1
#define LOG_U    2
#define Ln10     (2.302585093)
#define U_PG     LOG_U
#define U_PL     NOLOG_U

/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])

#if (U_PG == LOG_U)
  #define P_G(n)     (exp(Ln10*UNKNOWN(n,U_P_G)))
  #define P_Gn(n)    (exp(Ln10*UNKNOWNn(n,U_P_G)))
#else
  #define P_G(n)     (UNKNOWN(n,U_P_G))
  #define P_Gn(n)    (UNKNOWNn(n,U_P_G))
#endif
#if (U_PL == LOG_U)
  #define P_L(n)     (exp(Ln10*UNKNOWN(n,U_P_L)))
  #define P_Ln(n)    (exp(Ln10*UNKNOWNn(n,U_P_L)))
#else
  #define P_L(n)     (UNKNOWN(n,U_P_L))
  #define P_Ln(n)    (UNKNOWNn(n,U_P_L))
#endif
#define ZN_Ca_S(n)  (UNKNOWN(n,U_ZN_Ca_S))
#define ZN_Si_S(n)  (UNKNOWN(n,U_ZN_Si_S))
#define PSI(n)      (UNKNOWN(n,U_PSI))
#define C_K(n)      (UNKNOWN(n,U_C_K))
#define C_Cl(n)     (UNKNOWN(n,U_C_Cl))

#define ZN_Ca_Sn(n) (UNKNOWNn(n,U_ZN_Ca_S))
#define ZN_Si_Sn(n) (UNKNOWNn(n,U_ZN_Si_S))
#define PSIn(n)     (UNKNOWNn(n,U_PSI))
#define C_Kn(n)     (UNKNOWNn(n,U_C_K))
#define C_Cln(n)    (UNKNOWNn(n,U_C_Cl))

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
#define W_K_D      (f[(25)])
#define W_K_psi    (f[(26)])
#define W_K_dif    (f[(27)])

#define N_Cn(n)    (f_n[(n)])
#define N_qn(n)    (f_n[(2+n)])
#define N_Can(n)   (f_n[(4+n)])
#define N_Sin(n)   (f_n[(6+n)])
#define N_Kn(n)    (f_n[(8+n)])
#define N_Cln(n)   (f_n[(10+n)])
#define M_n(n)     (f_n[(12+n)])
#define N_CHn(n)   (f_n[(21+n)])
#define N_CCn(n)   (f_n[(23+n)])



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

#define KD_GCO2     (va[(48)])
#define KD_GH2O     (va[(49)])
#define KD_GM       (va[(50)])

#define V_S0(n)     (v0[(0+n)])


/*
  Aqueous solution
*/

/* Charge of ions */
#include "ElectricChargeOfIonInWater.h"
#define z_ca          ElectricChargeOfIonInWater(Ca)
#define z_h           ElectricChargeOfIonInWater(H)
#define z_oh          ElectricChargeOfIonInWater(OH)
#define z_hco3        ElectricChargeOfIonInWater(HCO3)
#define z_co3         ElectricChargeOfIonInWater(CO3)
#define z_h3sio4      ElectricChargeOfIonInWater(H3SiO4)
#define z_cahco3      ElectricChargeOfIonInWater(CaHCO3)
#define z_cah3sio4    ElectricChargeOfIonInWater(CaH3SiO4)
#define z_h2sio4      ElectricChargeOfIonInWater(H2SiO4)
#define z_caoh        ElectricChargeOfIonInWater(CaOH)
#define z_k           ElectricChargeOfIonInWater(K)
#define z_cl          ElectricChargeOfIonInWater(Cl)

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
#define M_Cl           (MolarMassOfMolecule(Cl)*gr)

/* Partial Molar Volume of Ions (dm3/mole) */
#include "PartialMolarVolumeOfMoleculeInWater.h"
#define dm3            (1.e3)
#define v_h            (PartialMolarVolumeOfMoleculeInWater(H)*dm3)
#define v_oh           (PartialMolarVolumeOfMoleculeInWater(OH)*dm3)
#define v_h2o          (PartialMolarVolumeOfMoleculeInWater(H2O)*dm3)
#define v_co2          (PartialMolarVolumeOfMoleculeInWater(CO2)*dm3)
#define v_h2co3        (PartialMolarVolumeOfMoleculeInWater(H2CO3)*dm3)
#define v_hco3         (PartialMolarVolumeOfMoleculeInWater(HCO3)*dm3)
#define v_co3          (PartialMolarVolumeOfMoleculeInWater(CO3)*dm3)
#define v_ca           (PartialMolarVolumeOfMoleculeInWater(Ca)*dm3)
#define v_h3sio4       (PartialMolarVolumeOfMoleculeInWater(H3SiO4)*dm3)
#define v_h2sio4       (PartialMolarVolumeOfMoleculeInWater(H2SiO4)*dm3)
#define v_cah2sio4     (PartialMolarVolumeOfMoleculeInWater(CaH2SiO4)*dm3)
#define v_cah3sio4     (PartialMolarVolumeOfMoleculeInWater(CaH3SiO4)*dm3)
#define v_caco3aq      (PartialMolarVolumeOfMoleculeInWater(CaCO3)*dm3)
#define v_caoh         (PartialMolarVolumeOfMoleculeInWater(CaOH)*dm3)
#define v_k            (PartialMolarVolumeOfMoleculeInWater(K)*dm3)
#define v_koh          (PartialMolarVolumeOfMoleculeInWater(KOH)*dm3)
#define v_cl           (PartialMolarVolumeOfMoleculeInWater(Cl)*dm3)
#define v_h4sio4       (PartialMolarVolumeOfMoleculeInWater(H4SiO4)*dm3)
#define v_cahco3       (PartialMolarVolumeOfMoleculeInWater(CaHCO3)*dm3)
#define v_na           (PartialMolarVolumeOfMoleculeInWater(Na)*dm3)
#define v_naoh         (PartialMolarVolumeOfMoleculeInWater(NaOH)*dm3)
#define v_naco3        (PartialMolarVolumeOfMoleculeInWater(NaCO3)*dm3)
#define v_nahco3       (PartialMolarVolumeOfMoleculeInWater(NaHCO3)*dm3)




/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CC   = Calcium Carbonate (Calcite)
  CSH  = Calcium Silicate Hydrate
  SH   = Amorphous Silica Gel
*/

/* Material Properties */
#define SaturationDegree(x)                Curve_ComputeValue(Element_GetCurve(el),x)
#define RelativePermeabilityToLiquid(x)    Curve_ComputeValue(Element_GetCurve(el) + 1,x)
#define RelativePermeabilityToGas(x)       Curve_ComputeValue(Element_GetCurve(el) + 2,x)
#define TortuosityToGas(phi,sg)            (pow(phi,1.74)*pow(sg,3.20))
#define TortuosityToLiquid1(phi,sl)        (0.00029*exp(9.95*phi)/(1+625*pow((1-sl),4)))
#define TortuosityToLiquid(phi,sl)         (MIN(phi/0.25,1)*TortuosityToLiquid1(phi,sl))
#define PermeabilityChange(phi,phi0)       (pow(phi/phi0,3)*pow(((1-phi0)/(1-phi)),2))


/* C-S-H Properties */
#define CalciumSiliconRatioInCSH(x)        Curve_ComputeValue(Element_GetCurve(el) + 3,x)
#define WaterSiliconRatioInCSH(x)          Curve_ComputeValue(Element_GetCurve(el) + 4,x)
#define MolarVolumeOfCSH(x)                Curve_ComputeValue(Element_GetCurve(el) + 5,x)


/* S-H Properties */
/* Saturation Degree of Dissolved S-H */
#define S_SHeq(x)                             Curve_ComputeValue(Element_GetCurve(el) + 6,x)
#define SaturationDegreeOfSH(s_ch,zn_si_s)    (NEGEXP(zn_si_s)*S_SHeq(s_ch))
/* Ion Activity Product of Dissolved S-H */
#define IonActivityProductOfSH(s_sh)          (K_SH*(s_sh))


/* CH Properties */
/* Molar Volume of CH (dm3/mole) */
#define V_CH       (33.e-3)      /* (33.e-3) */
/* Threshold of CO2 concentration */
#define C_CO2_eq                              (K_h2o*K_h2o*K_CC/(K_hco3*K_co3*K_CH))
/* Saturation Degree of Dissolved CH */
#define SaturationDegreeOfCH(z_co2,zn_ca_s)   (NEGEXP(zn_ca_s)/MAX(z_co2,1.))
/* Ion Activity Product of Dissolved CH */
#define IonActivityProductOfCH(s_ch)          (K_CH*(s_ch))


/* CC Properties */
/* Molar Volume of CC (dm3/mole) */
#define V_CC       (37.e-3)      /* (37.e-3) */
/* Saturation Degree of Dissolved CC */
#define SaturationDegreeOfCC(z_co2,zn_ca_s)   (NEGEXP(zn_ca_s)*MIN(z_co2,1.))
/* Ion Activity Product of Dissolved CC */
#define IonActivityProductOfCC(s_cc)          (K_CC*(s_cc))


/* Element contents in solid phases  */
#define CalciumContentInCHAndCC(zn_ca_s)   (n_ca_ref*MAX(zn_ca_s,0.))
#define SiliconContentInCSH(zn_si_s)       (n_si_ref*MAX(zn_si_s,0.))


/* Pure CO2 Gas Properties */
#define MolarDensityOfCO2(P,T)   Curve_ComputeValue(Element_GetCurve(el) + 7,P) /* product_gco2(P,T) */
#define ViscosityOfCO2(P,T)      Curve_ComputeValue(Element_GetCurve(el) + 8,P) /* product_muco2(P,T) */


/* CO2-H2O Gas Mixture Properties */
#define ConcentrationOfCO2InGas               MolarDensityOfCO2 /* assumed as pure co2 gas */
#define ConcentrationOfH2OInGas(PG,PL,T)      (MoleFractionOfH2OInGas(PG,PL,T)*MolarDensityOfCO2(PG,T))
#define MoleFractionOfH2OInGas(PG,PL,T)       molefractionH2O(el,PG,PL,T)
#define FugacityCoefficientOfCO2InGas(PG,T)   fugacityCO2(el,PG,T)
#define FugacityCoefficientOfH2OInGas(PG,T)   fugacityH2O(el,PG,T)
    
    
/* Solubility of CO2 */
#define ConcentrationOfCO2InLiquid(PG,PL,T)   solubilityCO2b(el,PG,PL,T)
#define EquilibriumLiquidGasOfCO2             equilibriumCO2
#define EquilibriumLiquidGasOfH2O             equilibriumH2O


/* Concentration of OH computed from electroneutrality */
#define ConcentrationOfOHInLiquid(A,B,C,D,E)  concentration_oh(A,B,C,D,E,el)


#define NEGEXP(y)  ((y < 0.) ? exp(y) : 1.)



/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])



/* Internal Functions */
static int    pm(const char *s) ;
static double concentration_oh(double,double,double,double,double,Element_t*) ;
static double poly4(double,double,double,double,double) ;
static void   ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void   ComputeFluxes(elem_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;
static int    TangentCoefficients(Element_t*,double,double*) ;
static double fugacityCO2(Element_t*,double,double) ;
static double fugacityH2O(Element_t*,double,double) ;
static double equilibriumH2O(double,double) ;
static double equilibriumCO2(double,double) ;
static double (solubilityCO2)(Element_t*,double,double,double) ;
static double (solubilityCO2b)(Element_t*,double,double,double) ;
static double (molefractionH2O)(Element_t*,double,double,double) ;
static void   ComputePhysicoChemicalProperties(double) ;



/* Internal Parameters */
static double temperature ;
static double phi0 ;
static double k_intl,k_intg ;
static double t_ch,t_cc ;
static double c_co2_eq ;
static double n_ca_ref,n_si_ref ;

static double K_h2o ;
static double K_hco3,K_co3 ;
static double K_h2sio4,K_h3sio4 ;
static double K_cahco3,K_caco3aq,K_caoh ;
static double K_cah2sio4,K_cah3sio4 ;

static double K_CH,K_CC,K_SH ;

static double d_h,d_oh ;
static double d_co2,d_hco3,d_co3 ;
static double d_ca,d_caoh ;
static double d_cahco3,d_caco3aq ;
static double d_h4sio4,d_h3sio4,d_h2sio4 ;
static double d_cah2sio4,d_cah3sio4 ;
static double d_k,d_koh ;
static double d_cl ;
static double d_na ;

static double mu_l ;
static double FsRT ;


#include "WaterViscosity.h"
#include "DiffusionCoefficientOfMoleculeInWater.h"
#include "EquilibriumConstantOfHomogeneousReactionInWater.h"
#include "DissociationConstantOfCementHydrationProduct.h"
#include "DissociationConstantOfCalciumCarbonate.h"
#include "PhysicalConstant.h"
#define dm2            (1.e2)

void ComputePhysicoChemicalProperties(double TK)
{
  K_h2o      = EquilibriumConstantOfHomogeneousReactionInWater(H2O__H_OH,TK) ;
  K_hco3     = EquilibriumConstantOfHomogeneousReactionInWater(CO2_H2O__HCO3_H,TK) ;
  K_co3      = EquilibriumConstantOfHomogeneousReactionInWater(HCO3__H_CO3,TK) ;
  K_h3sio4   = EquilibriumConstantOfHomogeneousReactionInWater(H4SiO4__H3SiO4_H,TK) ;
  K_h2sio4   = EquilibriumConstantOfHomogeneousReactionInWater(H3SiO4_OH__H2SiO4_H2O,TK) ;
  K_cahco3   = EquilibriumConstantOfHomogeneousReactionInWater(Ca_HCO3__CaHCO3,TK) ;
  K_caco3aq  = EquilibriumConstantOfHomogeneousReactionInWater(Ca_CO3__CaCO3,TK) ;
  K_cah2sio4 = EquilibriumConstantOfHomogeneousReactionInWater(Ca_H2SiO4__CaH2SiO4,TK) ;
  K_cah3sio4 = EquilibriumConstantOfHomogeneousReactionInWater(Ca_H3SiO4__CaH3SiO4,TK) ;
  K_caoh     = EquilibriumConstantOfHomogeneousReactionInWater(Ca_OH__CaOH,TK) ;

  K_CH = DissociationConstantOfCementHydrationProduct(CH__Ca_2OH,TK) ;
  K_CC = DissociationConstantOfCalciumCarbonate(Calcite__Ca_CO3,TK) ;
  K_SH = DissociationConstantOfCementHydrationProduct(S_2H2O__H4SiO4,TK) ;
  
  /* Diffusion Coefficient Of Molecules In Water (dm2/s) */
  d_oh         = DiffusionCoefficientOfMoleculeInWater(OH,TK)*dm2 ;
  d_h          = DiffusionCoefficientOfMoleculeInWater(H,TK)*dm2 ;
  d_co2        = DiffusionCoefficientOfMoleculeInWater(CO2,TK)*dm2 ;
  d_hco3       = DiffusionCoefficientOfMoleculeInWater(HCO3,TK)*dm2 ;
  d_co3        = DiffusionCoefficientOfMoleculeInWater(CO3,TK)*dm2 ;
  d_ca         = DiffusionCoefficientOfMoleculeInWater(Ca,TK)*dm2 ;

  d_caoh       = DiffusionCoefficientOfMoleculeInWater(CaOH,TK)*dm2 ;
  d_cahco3     = DiffusionCoefficientOfMoleculeInWater(CaHCO3,TK)*dm2 ;
  d_h4sio4     = DiffusionCoefficientOfMoleculeInWater(H4SiO4,TK)*dm2 ;
  d_h3sio4     = DiffusionCoefficientOfMoleculeInWater(H3SiO4,TK)*dm2 ;
  d_h2sio4     = DiffusionCoefficientOfMoleculeInWater(H2SiO4,TK)*dm2 ;
  d_cah2sio4   = DiffusionCoefficientOfMoleculeInWater(CaH2SiO4,TK)*dm2 ;
  d_cah3sio4   = DiffusionCoefficientOfMoleculeInWater(CaH3SiO4,TK)*dm2 ;

  d_caco3aq    = DiffusionCoefficientOfMoleculeInWater(CaCO3,TK)*dm2 ;
  d_k          = DiffusionCoefficientOfMoleculeInWater(K,TK)*dm2 ;
  d_koh        = DiffusionCoefficientOfMoleculeInWater(KOH,TK)*dm2 ;
  d_cl         = DiffusionCoefficientOfMoleculeInWater(Cl,TK)*dm2 ;
  d_na         = DiffusionCoefficientOfMoleculeInWater(Na,TK)*dm2 ;
  
  /* Viscosity (Pa.s) */
  mu_l       = WaterViscosity(TK) ;
  
  /* Physical constants */
  {
    double RT  = PhysicalConstant(PerfectGasConstant)*TK ;
    FsRT       = PhysicalConstant(Faraday)/RT ;
  }
}


#define NbOfComponents    (52)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;


#define I_C_H          (7)
#define I_C_OH         (8)
#define I_C_H2O        (9)

#define I_C_CO2        (10)
#define I_C_HCO3       (11)
#define I_C_H2CO3      (xxx)
#define I_C_CO3        (13)

#define I_C_Ca         (14)
#define I_C_CaOH       (15)
#define I_C_CaHCO3     (16)
#define I_C_CaCO3aq    (17)
#define I_C_CaOH2aq    (xxx)

#define I_C_H2SiO4     (19)
#define I_C_H3SiO4     (20)
#define I_C_H4SiO4     (21)

#define I_C_CaH2SiO4   (22)
#define I_C_CaH3SiO4   (23)

#define I_C_K          (24)
#define I_C_Cl         (25)

#define I_S_CH         (26)
#define I_Q_SH         (27)
#define I_Q_CC         (28)

#define I_C_CO2_G      (29)
#define I_C_H2O_G      (30)

#define I_RHO_L        (31)
#define I_RHO_G        (32)
#define I_P_G          (33)
#define I_P_L          (34)

#define I_N_CH         (35)
#define I_N_CC         (36)
#define I_V_S          (37)

#define I_N_C          (38)
#define I_N_Ca         (39)
#define I_N_Si         (40)
#define I_N_K          (41)
#define I_N_Cl         (42)
#define I_Mass         (43)
#define I_N_Q          (44)

#define I_Phi          (45)

#define I_PSI          (46)

#define I_ZN_Ca_S      (47)
#define I_ZN_Si_S      (48)

#define I_N_CHn        (49)
#define I_N_CCn        (50)
#define I_V_S0         (51)

#define NbOfComponentFluxes    (7)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_C           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_K           (3)
#define I_W_Cl          (4)
#define I_W_Mass        (5)
#define I_W_Q           (6)

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0)      return (0) ;
  else if(strcmp(s,"N_CH") == 0)     return (1) ;
  else if(strcmp(s,"N_Si") == 0)     return (2) ;
  else if(strcmp(s,"T_CH") == 0)     return (3) ;
  else if(strcmp(s,"T_CC") == 0)     return (4) ;
  else if(strcmp(s,"k_intl") == 0)    return (5) ;
  else if(strcmp(s,"k_intg") == 0)    return (6) ;
  else if(strcmp(s,"temperature") == 0)    return (7) ;
  else if(strcmp(s,"C_CO2_eq") == 0) return (8) ;
  else return(-1) ;
}

int SetModelProp(Model_t *model)
/** Set the model properties 
 *  Return 0 */
{
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_C,"carbone") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_Si,"silicium") ;
  Model_CopyNameOfEquation(model,E_q,"charge") ;
  Model_CopyNameOfEquation(model,E_K,"potassium") ;
  Model_CopyNameOfEquation(model,E_Cl,"chlorine") ;
  Model_CopyNameOfEquation(model,E_mass,"mass") ;
  
  /** Names of the main (nodal) unknowns */
#if (U_PG == LOG_U)
  Model_CopyNameOfUnknown(model,U_P_G,"logp_co2") ;
#else
  Model_CopyNameOfUnknown(model,U_P_G,"p_co2") ;
#endif
#if (U_PL == LOG_U)
  Model_CopyNameOfUnknown(model,U_P_L,"logp_l") ;
#else
  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,U_ZN_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,U_PSI,"psi") ;
  Model_CopyNameOfUnknown(model,U_C_K,"c_k") ;
  Model_CopyNameOfUnknown(model,U_C_Cl,"c_cl") ;
  
  return(0) ;
}



int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 9 ;

  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("N_CH")] = 1 ;
    Material_GetProperty(mat)[pm("N_Si")] = 1 ;
    Material_GetProperty(mat)[pm("N_CC")] = 0 ;
    Material_GetProperty(mat)[pm("T_CH")] = 600 ;
    Material_GetProperty(mat)[pm("T_CC")] = 0 ;
  
    Material_ScanProperties(mat,datafile,pm) ;

    t_ch      = Material_GetProperty(mat)[pm("T_CH")] ;
    t_cc      = Material_GetProperty(mat)[pm("T_CC")] ;

    if(t_cc  == 0.) Material_GetProperty(mat)[pm("T_CC")]  = t_ch ;
    
    temperature  = Material_GetProperty(mat)[pm("temperature")] ;
    ComputePhysicoChemicalProperties(temperature) ;
    Material_GetProperty(mat)[pm("C_CO2_eq")] = C_CO2_eq ;
  }
  
  return(NbOfProp) ;
}




int PrintModelProp(Model_t *model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;

  printf("Description:\n") ;
  printf("We take into account the following properties:\n ") ;
  printf("\t - alcalis\n") ;
  printf("\t - general approach for CSH\n") ;
  printf("\t - water release\n") ;
  
  printf("\n") ;
  printf("The set of 7 equations is:\n") ;
  printf("\t - Mass balance of C      (carbone)\n") ;
  printf("\t - Mass balance of Ca     (calcium)\n") ;
  printf("\t - Mass balance of Si     (silicium)\n") ;
  printf("\t - Mass balance of K      (potassium)\n") ;
  printf("\t - Mass balance of Cl     (chlorine)\n") ;
  printf("\t - Total mass balance     (mass)\n") ;
  printf("\t - Charge balance         (charge)\n") ;

  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length : dm !\n") ;
  printf("\t time   : s !\n") ;
  
  printf("\n") ;
  printf("The 7 primary unknowns are:\n") ;
  printf("\t - Gas pressure                 (p_co2)\n") ;
  printf("\t - Liquid pressure              (p_l)\n") ;
  printf("\t - Electric potential           (psi)\n") ;
  printf("\t - Zeta unknown for calcium     (zn_ca_s)\n") ;
  printf("\t   \t zn_ca_s is defined as:\n") ;
  printf("\t   \t zn_ca_s = n_ch/n0 + log(s_ch)  for c_co2 < c_co2_eq\n") ;
  printf("\t   \t zn_ca_s = n_cc/n0 + log(s_cc)  for c_co2 > c_co2_eq\n") ;
  printf("\t - Zeta unknown for silicon     (zn_si_s)\n") ;
  printf("\t   \t zn_si_s is defined as:\n") ;
  printf("\t   \t zn_si_s = n_si/n0 + log(s_sh/s_sh_eq)\n") ;
  printf("\t - Potassium concentration      (c_k)\n") ;
  printf("\t - Chlorine concentration       (c_cl)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"N_CH  = 6.1    # Contenu en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"N_Si  = 2.4    # contenu en CSH (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4    # contenu en K (moles/L)\n") ;  
  fprintf(ficd,"N_Cl  = 0.4    # contenu en Cl (moles/L)\n") ;
  fprintf(ficd,"T_CH  = 1.e5   # Cinetique de dissolution de CH (s)\n") ;
  fprintf(ficd,"Curves = csh-properties # Nom du fichier: S_CH X Z V\n") ;

  return(NEQ) ;
}




int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
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
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double *f  = Element_GetImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    i ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature = GetProperty("temperature") ;
  
  
  /* Contents */
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x         = ComputeComponents(el,u,f,0,i) ;
    
    /* Liquid contents */
    double c_co2      = x[I_C_CO2] ;
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    
    double s_ch       = x[I_S_CH] ;
    

    /* Solid contents */
    /* ... as components: CH, CC, CSH */
    double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
    double n_ch       = (z_co2 <= 1) ? n_ch_cc  : 0 ;
    double n_cc       = (z_co2 >  1) ? n_ch_cc  : 0 ;
    /* ... as elements: Ca, Si, C */
    double n_si_s     = SiliconContentInCSH(zn_si_s) ;
    /* ... as volume */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double v_s        = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;

    
    N_CH(i)    = n_ch ;
    N_CC(i)    = n_cc ;

    V_S0(i)    = v_s ;
    
    ComputeComponents(el,u,f,0,i) ;
    
    /* Back up */
    N_C(i)  = x[I_N_C] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ;                              
    N_Cl(i) = x[I_N_Cl] ;
    M(i)    = x[I_Mass] ;
    N_q(i)  = x[I_N_Q] ;
    
  }
  
  /* We skip if the element is a submanifold */
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
  int i ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature = GetProperty("temperature") ;
  
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    
    /* Back up */
    N_C(i)  = x[I_N_C] ;
    N_Ca(i) = x[I_N_Ca] ;
    N_Si(i) = x[I_N_Si] ;
    N_K(i)  = x[I_N_K] ;                              
    N_Cl(i) = x[I_N_Cl] ;
    M(i)    = x[I_Mass] ;
    N_q(i)  = x[I_N_Q] ;

    N_CH(i) = x[I_N_CH] ;
    N_CC(i) = x[I_N_CC] ;

    {
      double c_co2      = x[I_C_CO2] ;
      double p_g	      = x[I_P_G] ;
      if(c_co2 < 0 || p_g < 0) {
        double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        double c_oh       = x[I_C_OH] ;
        double c_h3sio4   = x[I_C_H3SiO4] ;
        double zn_si_s    = x[I_ZN_Si_S] ;
        double zn_ca_s    = x[I_ZN_Ca_S] ;
        double s_ch       = x[I_S_CH] ;
        double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
        double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
        double n_si_s     = SiliconContentInCSH(zn_si_s) ;
        double n_ca_s     = n_ch_cc + x_csh*n_si_s ;
        double n_cc       = x[I_N_CC] ;
        printf("x         = %e\n",x0) ;
        printf("c_co2     = %e\n",c_co2) ;
        printf("n_cc      = %e\n",n_cc) ;
        printf("n_ca_s    = %e\n",n_ca_s) ; 
        printf("n_si_s    = %e\n",n_si_s) ;
        printf("zn_si_s   = %e\n",zn_si_s) ;
        printf("zn_ca_s   = %e\n",zn_ca_s) ;
        printf("c_h3sio4  = %e\n",c_h3sio4) ;
        printf("c_oh      = %e\n",c_oh) ;
        printf("p_g       = %e\n",p_g) ;
        return(-1) ;
      }
    }
    
    {
      double phi        = x[I_Phi] ;
      if(phi < 0) {
        printf("phi = %e\n",phi) ;
        return(-1) ;
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
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }

#if (U_PG == LOG_U)
  for(i=0;i<2*NEQ;i++){
    K(i,U_P_G)     *= Ln10*P_G(0) ;
    K(i,U_P_G+NEQ) *= Ln10*P_G(1) ;
  }
#endif
#if (U_PL == LOG_U)
  for(i=0;i<2*NEQ;i++){
    K(i,U_P_L)     *= Ln10*P_L(0) ;
    K(i,U_P_L+NEQ) *= Ln10*P_L(1) ;
  }
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
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
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nso = 60 ;
  int    i ;

  /* if(el.dim < dim) return(0) ; */

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }
  
  /*
    Input Data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature = GetProperty("temperature") ;
  
  /* Quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Liquid and gas components */
    double *x         = ComputeComponents(el,u,f,0.,j) ;
    
    double p_l     	  = x[I_P_L] ;
    double p_g	      = x[I_P_G] ;
    double gc_co2     = x[I_C_CO2_G] ;
    double gc_h2o     = x[I_C_H2O_G] ;
    double c_co2      = x[I_C_CO2] ;
    double z_co2      = c_co2/c_co2_eq ;
    double zn_si_s    = x[I_ZN_Si_S] ;
    double zn_ca_s    = x[I_ZN_Ca_S] ;
    double c_k        = x[I_C_K] ;                 
    double c_cl       = x[I_C_Cl] ;
    double p_c        = p_g - p_l ;
    double s_l        = SaturationDegree(p_c) ;
    double s_g     	  = 1 - s_l ;
    
    double s_ch       = x[I_S_CH] ;

    double c_h4sio4   = x[I_C_H4SiO4] ;
    double c_oh       = x[I_C_OH] ;
    double c_hco3     = x[I_C_HCO3] ;
    double c_co3      = x[I_C_CO3] ;
    double c_ca       = x[I_C_Ca] ;
    double c_cahco3   = x[I_C_CaHCO3] ; 
    
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_caco3aq  = x[I_C_CaCO3aq] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_cah2sio4 = x[I_C_CaH2SiO4] ;
    double c_caoh     = x[I_C_CaOH] ;

    double rho_l      = x[I_RHO_L] ;

    double rho_g      = x[I_RHO_G] ;
    
    double y_h2o 	    = MoleFractionOfH2OInGas(p_g,p_l,temperature) ;
    double y_h2os 	  = MoleFractionOfH2OInGas(p_g,p_g,temperature) ;
    double gc_h2os    = ConcentrationOfH2OInGas(p_g,p_g,temperature) ;
    double RH		      = MIN(ConcentrationOfH2OInGas(p_l,p_l,temperature)/gc_h2os,1.) ;

    /* Charge density */
    double c_q = x[I_N_Q] ;
    
    /* Solid contents */
    /* ... as components: CH, CC, CSH */
    double n_ch       = N_CH(j) ;
    double n_cc       = N_CC(j) ;
    double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
    double z_csh      = WaterSiliconRatioInCSH(s_ch) ;
    /* ... as elements: Ca, Si */
    double n_si_s     = SiliconContentInCSH(zn_si_s) ;
    /* ... as mass */
    double m_csh      = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
    double m_s        = M_CaOH2*n_ch + M_CaCO3*n_cc + m_csh ;
    /* ... as volume */
    double v_csh      = MolarVolumeOfCSH(s_ch) ;
    double v_s        = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;

    /* Porosity */
    double v_s0       = V_S0(j) ;
    double phi        = phi0 + v_s0 - v_s ;
    double phimax 	  = MAX(phi,0.) ; 
    
    
    double psi        = x[I_PSI] ;
    
    double c_h2o      = x[I_C_H2O] ;
   
    /* Total mass */
    double M    = s_g*phi*rho_g + s_l*phi*rho_l + m_s ;
   
    double mu_co2   = ViscosityOfCO2(p_g,temperature) ; 
    double ph = 14 + log(c_oh)/log(10.) ;

    i = 0 ;
    Result_Store(r + i++,&c_co2,"c_co2",1) ;
    Result_Store(r + i++,&ph,"ph",1) ;
    {double k_hco3 = K_hco3 ; Result_Store(r + i++,&k_hco3,"k_hco3",1) ; }
    Result_Store(r + i++,&c_hco3,"c_hco3",1) ;
    Result_Store(r + i++,&c_co3,"c_co3",1) ;
    Result_Store(r + i++,&c_ca,"c_ca",1) ;
    Result_Store(r + i++,&c_cahco3,"c_cahco3",1) ;
    Result_Store(r + i++,&c_cah3sio4,"c_cah3sio4",1) ;
    Result_Store(r + i++,&c_h3sio4,"c_h3sio4",1) ;
    Result_Store(r + i++,&c_h4sio4,"c_h4sio4",1) ;
    Result_Store(r + i++,&n_ch,"n_ch",1) ;
    Result_Store(r + i++,&n_cc,"n_cc",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&c_oh,"c_oh",1) ;
    Result_Store(r + i++,&psi,"potentiel_electrique",1) ;
    Result_Store(r + i++,&c_q,"charge",1) ;
    Result_Store(r + i++,&zn_ca_s,"zn_ca_s",1) ;
    Result_Store(r + i++,&z_co2,"z_co2",1) ;
    Result_Store(r + i++,&c_caco3aq,"c_caco3aq",1) ;
    Result_Store(r + i++,&c_caoh,"c_caoh",1) ;
    Result_Store(r + i++,&c_h2sio4,"c_h2sio4",1) ;
    Result_Store(r + i++,&c_cah2sio4,"c_cah2sio4",1) ;
    Result_Store(r + i++,&c_k,"c_k",1) ;
    {double pk_ch =log10(s_ch) ; Result_Store(r + i++,&pk_ch,"pk_ch",1) ; }
    {double s_cheq = S_SHeq(s_ch) ; Result_Store(r + i++,&s_cheq,"x_s",1) ; }
    Result_Store(r + i++,&zn_si_s,"zn_si_s",1) ;
    Result_Store(r + i++,&c_cl,"c_cl",1) ;
    Result_Store(r + i++,&M,"M",1) ;
    Result_Store(r + i++,&p_l,"p_l",1) ;
    Result_Store(r + i++,&c_h2o,"c_h2o",1) ;
    Result_Store(r + i++,&W_M,"W_M",1) ;
    Result_Store(r + i++,&rho_l,"rho_l",1) ;
    {double m_l = s_l*phi*rho_l ; Result_Store(r + i++,&m_l,"M_L",1) ; }
    Result_Store(r + i++,&m_s,"M_S",1) ;
    Result_Store(r + i++,&p_g,"p_co2",1) ;
    Result_Store(r + i++,&s_l,"s_l",1) ;
    Result_Store(r + i++,&s_g,"s_g",1) ;
    Result_Store(r + i++,&gc_co2,"gc_co2",1) ;
    Result_Store(r + i++,&x_csh,"x_csh",1) ;
    Result_Store(r + i++,&z_csh,"z_csh",1) ;
    Result_Store(r + i++,&v_csh,"v_csh",1) ;
    Result_Store(r + i++,&temperature,"temperature",1) ;
    Result_Store(r + i++,&mu_co2,"mu_co2",1) ;
    Result_Store(r + i++,&gc_h2o,"gc_h2o",1) ;
    Result_Store(r + i++,&y_h2o,"y_h2o",1) ;
    {double a = WaterViscosity(298) ; Result_Store(r + i++,&a,"mu_l25",1) ; }
    {double a = d_oh ; Result_Store(r + i++,&a,"d_oh",1) ; }
    {double eta = NEGEXP(zn_si_s) ; Result_Store(r + i++,&eta,"eta",1) ; }
    {double a = K_h2sio4 ; Result_Store(r + i++,&a,"K_h2sio4",1) ; }
    Result_Store(r + i++,&phimax,"phimax",1) ;
    Result_Store(r + i++,&RH,"RH",1) ;
    {double a = K_CH ; Result_Store(r + i++,&a,"kch",1) ; }
    Result_Store(r + i++,&y_h2os,"yh2os",1) ;
    Result_Store(r + i++,&W_Ca,"W_Ca",1) ;
    Result_Store(r + i++,&W_Si,"W_Si",1) ;
    Result_Store(r + i++,&W_C,"W_C",1) ;
    Result_Store(r + i++,&W_K,"W_K",1) ;
    Result_Store(r + i++,&W_K_D,"W_K_D ",1) ;
    Result_Store(r + i++,&W_K_psi,"W_K_psi ",1) ;
    Result_Store(r + i++,&W_K_dif,"W_K_dif ",1) ;
  }

  return(nso) ;
}




void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int    i ;
  
  /*
    Input data
  */
  k_intl    = GetProperty("k_intl") ;
  k_intg    = GetProperty("k_intg") ;
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature = GetProperty("temperature") ;


  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* Liquid and gas components */
    double *x         = ComputeComponents(el,u,f,0.,i) ;
    
    double p_l     	  = x[I_P_L] ;
    double p_g	      = x[I_P_G] ;
    double gc_co2     = x[I_C_CO2_G] ;
    double gc_h2o     = x[I_C_H2O_G] ;
    double c_co2      = x[I_C_CO2] ;
    double c_k        = x[I_C_K] ;                 
    double c_cl       = x[I_C_Cl] ;
    double p_c        = p_g - p_l ;
    double s_l        = SaturationDegree(p_c) ;

    double c_h4sio4   = x[I_C_H4SiO4] ;
    double c_oh       = x[I_C_OH] ;
    double c_h        = x[I_C_H] ;
    double c_hco3     = x[I_C_HCO3] ;
    double c_co3      = x[I_C_CO3] ;
    double c_ca       = x[I_C_Ca] ;
    double c_cahco3   = x[I_C_CaHCO3] ; 
    
    double c_h3sio4   = x[I_C_H3SiO4] ;
    double c_cah3sio4 = x[I_C_CaH3SiO4] ;
    double c_caco3aq  = x[I_C_CaCO3aq] ;
    double c_h2sio4   = x[I_C_H2SiO4] ;
    double c_cah2sio4 = x[I_C_CaH2SiO4] ;
    double c_caoh     = x[I_C_CaOH] ;

    double rho_l      = x[I_RHO_L] ;
    double rho_g      = x[I_RHO_G] ;

    /* Porosity */
    double phi        = x[I_Phi] ;
       
    /* Evolution of permeabilities */
    double perm       = PermeabilityChange(phi,phi0) ;
    /* Permeability to liquid */
    double k_l  = (k_intl/mu_l)*RelativePermeabilityToLiquid(p_c)*perm ;
    /* Permeability to gas */
    double mu_co2   = ViscosityOfCO2(p_g,temperature) ;
    double k_g  = (k_intg/mu_co2)*RelativePermeabilityToGas(p_c)*perm ;
    
    /* Tortuosities */
    /* double tau  = TortuosityToGas(phi,s_g) ; */
    double iff  = TortuosityToLiquid(phi,s_l) ;
    
    KD_Ca     	+= c_ca*k_l ;
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
    KD_CO2      += c_co2*k_l ;
    
    KD_M      	+= rho_l*k_l ;
    
    KD_GCO2     += gc_co2*k_g ;
    KD_GH2O     += gc_h2o*k_g ;
    
    KD_GM       += rho_g*k_g ;
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

    Kpsi_q        += z_h*Kpsi_H + z_oh*Kpsi_OH \
                   + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 \
                   + z_ca*Kpsi_Ca + z_caoh*Kpsi_CaOH \
                   + z_cahco3*Kpsi_CaHCO3 \
                   + z_h3sio4*Kpsi_H3SiO4 + z_h2sio4*Kpsi_H2SiO4 \
                   + z_cah3sio4*Kpsi_CaH3SiO4 \
                   + z_k*Kpsi_K + z_cl*Kpsi_Cl ;
  }
  
  /* Averaging */
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}



double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *v0 = Element_GetConstantTerm(el) ;
  double *x = Components ;
  
  /* Primary Variables */
  x[U_P_L    ] = P_L(n) ;
  x[U_P_G    ] = P_G(n) ;
  x[U_ZN_Si_S] = ZN_Si_S(n) ;
  x[U_ZN_Ca_S] = ZN_Ca_S(n) ;
  x[U_C_Cl   ] = C_Cl(n) ;
  x[U_C_K    ] = C_K(n) ;
  x[U_PSI    ] = PSI(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn]  = N_CHn(n) ;
  x[I_N_CCn]  = N_CCn(n) ;
  x[I_V_S0 ]  = V_S0(n) ;
  
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  dx[U_P_L    ]  = x[U_P_L    ] ;
  dx[U_P_G    ]  = x[U_P_G    ] ;
  dx[U_ZN_Si_S]  = x[U_ZN_Si_S] ;
  dx[U_ZN_Ca_S]  = x[U_ZN_Ca_S] ;
  dx[U_C_Cl   ]  = x[U_C_Cl   ] ;
  dx[U_C_K    ]  = x[U_C_K    ] ;
  dx[U_PSI    ]  = x[U_PSI    ] ;
  
  dx[I_N_CHn]  = x[I_N_CHn] ;
  dx[I_N_CCn]  = x[I_N_CCn] ;
  dx[I_V_S0 ]  = x[I_V_S0 ] ;
  
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
  double p_l        = x[U_P_L] ;
  double p_g        = x[U_P_G] ; 
  double zn_si_s    = x[U_ZN_Si_S] ;
  double zn_ca_s    = x[U_ZN_Ca_S] ;
  double c_k        = x[U_C_K] ;                 
  double c_cl       = x[U_C_Cl] ;
  
  /* Gas components */
  double gc_co2     = ConcentrationOfCO2InGas(p_g,temperature) ;
  double gc_h2o     = ConcentrationOfH2OInGas(p_g,p_l,temperature) ;
  
  /* Liquid components */
  double c_co2      = ConcentrationOfCO2InLiquid(p_g,p_l,temperature) ;
  double z_co2      = c_co2/c_co2_eq ;

  double s_cc       = SaturationDegreeOfCC(z_co2,zn_ca_s) ;
  double Q_CC       = IonActivityProductOfCC(s_cc) ;
  double s_ch       = SaturationDegreeOfCH(z_co2,zn_ca_s) ;
  double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;

  double c_h4sio4   = Q_SH ;
  double c_oh       = ConcentrationOfOHInLiquid(c_co2,zn_si_s,zn_ca_s,c_k,c_cl) ;
  double c_h        = K_h2o/c_oh ;
  /*double c_h2co3    = K_h2co3*c_co2 ;*/
  double c_hco3     = K_hco3*c_co2/c_h ;
  double c_co3      = K_co3*c_hco3/c_h ;
  double c_ca       = Q_CC/c_co3 ;
  double c_cahco3   = K_cahco3*c_ca*c_hco3 ;
  double c_h3sio4   = K_h3sio4*c_h4sio4/c_h ;
  double c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4 ;
  double c_caco3aq  = K_caco3aq*c_ca*c_co3 ;
  double c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;
  double c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;
  double c_caoh     = K_caoh*c_ca*c_oh ;

  double c_h2o      = (1 - (c_co2*v_co2 + c_hco3*v_hco3 + c_co3*v_co3 \
                    + c_h*v_h + c_oh*v_oh \
                    + c_ca*v_ca + c_caoh*v_caoh \
                    + c_cahco3*v_cahco3 + c_caco3aq*v_caco3aq \
                    + c_h3sio4*v_h3sio4 + c_h4sio4*v_h4sio4 + c_h2sio4*v_h2sio4 \
                    + c_cah2sio4*v_cah2sio4 + c_cah3sio4*v_cah3sio4 \
                    + c_k*v_k + c_cl*v_cl))/v_h2o ;
                 
  double c_q        = z_h*c_h + z_oh*c_oh \
                    + z_hco3*c_hco3 + z_co3*c_co3 \
                    + z_ca*c_ca + z_caoh*c_caoh + z_cahco3*c_cahco3 \
                    + z_h2sio4*c_h2sio4 + z_h3sio4*c_h3sio4 \
                    + z_cah3sio4*c_cah3sio4 \
                    + z_k*c_k + z_cl*c_cl ;
    
  /* Solid contents */
  /* ... as components: CH, CC, CSH */
  double n_chn      = x[I_N_CHn] ;
  double n_ccn      = x[I_N_CCn] ;
  double n_ch_ci    = n_chn*pow(z_co2,-dt/t_ch) ;
  double n_cc_ci    = n_ccn*pow(z_co2,dt/t_cc) ;
  double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
  double n_ch       = (z_co2 <= 1) ? n_ch_cc - n_cc_ci : n_ch_ci ;
  double n_cc       = (z_co2 >  1) ? n_ch_cc - n_ch_ci : n_cc_ci ;
  double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double z_csh      = WaterSiliconRatioInCSH(s_ch) ;
  /* ... as elements: C, Ca, Si */
  double n_si_s     = SiliconContentInCSH(zn_si_s) ;
  double n_ca_s     = n_ch_cc + x_csh*n_si_s ;
  double n_c_s      = n_cc ;
  /* ... as mass */
  double m_csh      = (M_CaO*x_csh + M_SiO2 + M_H2O*z_csh)*n_si_s ;
  double m_s        = M_CaOH2*n_ch + M_CaCO3*n_cc + m_csh ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(s_ch) ;
  double v_s        = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;

  /* Porosity */
  double v_s0       = x[I_V_S0] ;
  double phi        = phi0 + v_s0 - v_s ;
  
  /* Saturation */
  double p_c        = p_g - p_l ;
  double s_l        = SaturationDegree(p_c) ;
  double s_g     	  = 1 - s_l ;
  
  /* Liquid contents */
  double phi_l      = s_l*phi ;
  /* ... as components */
  double n_co2      = phi_l*c_co2 ;
  double n_hco3     = phi_l*c_hco3 ;
  double n_co3      = phi_l*c_co3 ;
  double n_ca       = phi_l*c_ca ;
  double n_cahco3   = phi_l*c_cahco3 ;
  double n_h3sio4   = phi_l*c_h3sio4 ;
  double n_h4sio4   = phi_l*c_h4sio4 ;
  double n_cah3sio4 = phi_l*c_cah3sio4 ;
  double n_h2sio4   = phi_l*c_h2sio4 ;
  double n_cah2sio4 = phi_l*c_cah2sio4 ;
  double n_caco3aq  = phi_l*c_caco3aq ;
  double n_caoh     = phi_l*c_caoh ;
  double n_k        = phi_l*c_k ;
  double n_cl       = phi_l*c_cl ;
  /* ... as elements: C, Ca, Si */
  double n_c_l      = n_co2 + n_hco3 + n_co3 + n_cahco3 + n_caco3aq ;
  double n_ca_l     = n_ca + n_cahco3 + n_cah3sio4 + n_cah2sio4 + n_caco3aq + n_caoh ;
  double n_si_l     = n_h3sio4 + n_h4sio4 + n_cah3sio4 + n_h2sio4+ n_cah2sio4 ;
  /* ... as mass */
  double rho_l      = M_H2O*c_h2o + M_H*c_h + M_OH*c_oh \
                    + M_CO2*c_co2 + M_HCO3*c_hco3 + M_CO3*c_co3 \
                    + M_Ca*c_ca + M_CaOH*c_caoh \
                    + M_CaHCO3*c_cahco3 + M_CaCO3aq*c_caco3aq \
                    + M_H3SiO4*c_h3sio4 + M_H4SiO4*c_h4sio4 + M_H2SiO4*c_h2sio4 \
                    + M_CaH2SiO4*c_cah2sio4 + M_CaH3SiO4*c_cah3sio4 \
                    + M_Cl*c_cl + M_K*c_k ;
  double m_l        = phi_l*rho_l ;
    
  /* Gas contents */
  double phi_g      = s_g*phi ;
  /* ... as components */
  double n_co2_g    = phi_g*gc_co2 ;
  /* ... as elements */
  double n_c_g      = n_co2_g ;
  /* ... as mass */
  double rho_g      = M_CO2*gc_co2 + M_H2O*gc_h2o ;
  double m_g        = phi_g*rho_g ;

  /* Back up */
  
  
  /* Liquid components */
  x[I_C_H       ] = c_h ;
  x[I_C_OH      ] = c_oh ;
  x[I_C_H2O     ] = c_h2o ;

  x[I_C_CO2     ] = c_co2 ;
  x[I_C_HCO3    ] = c_hco3 ;
  x[I_C_CO3     ] = c_co3 ;

  x[I_C_Ca      ] = c_ca ;
  x[I_C_CaOH    ] = c_caoh ;

  x[I_C_CaHCO3  ] = c_cahco3 ;
  x[I_C_CaCO3aq ] = c_caco3aq ;

  x[I_C_H2SiO4  ] = c_h2sio4 ;
  x[I_C_H3SiO4  ] = c_h3sio4 ;
  x[I_C_H4SiO4  ] = c_h4sio4 ;

  x[I_C_CaH2SiO4] = c_cah2sio4 ;
  x[I_C_CaH3SiO4] = c_cah3sio4 ;
  
  x[I_C_K       ] = c_k ;
  
  x[I_C_Cl      ] = c_cl ;

  x[I_Q_SH      ] = Q_SH ;
  x[I_Q_CC      ] = Q_CC ;
  x[I_S_CH      ] = s_ch ;

  x[I_RHO_L     ] = rho_l ;
  
  x[I_P_L       ] = p_l ;
  x[I_P_G       ] = p_g ;
  
  /* Gas components */
  x[I_C_CO2_G ] = gc_co2 ;
  x[I_C_H2O_G ] = gc_h2o ;

  x[I_RHO_G   ] = rho_g ;
  
  /* Solid components */
  x[I_N_CH    ] = n_ch ;
  x[I_N_CC    ] = n_cc ;
  x[I_V_S     ] = v_s ;
  
  x[I_ZN_Si_S ] = zn_si_s ;
  x[I_ZN_Ca_S ] = zn_ca_s ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  /* Electric potential */
  x[I_PSI     ] = x[U_PSI] ;
  
  /* Element contents */
  x[I_N_C     ] = n_c_l  + n_c_g  + n_c_s ;
  x[I_N_Ca    ] = n_ca_l + n_ca_s ;
  x[I_N_Si    ] = n_si_l + n_si_s ;
  x[I_N_K     ] = n_k ;                              
  x[I_N_Cl    ] = n_cl ;
  x[I_Mass    ] = m_l + m_g + m_s ;
  x[I_N_Q     ] = c_q ;

  return ;
}



void ComputeFluxes(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double *grd = dComponents ;
  
  /*
    Input data
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature = GetProperty("temperature") ;

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
    W_Si    = w[I_W_Si] ;
    W_q     = w[I_W_Q] ;
    W_K     = w[I_W_K] ;
    W_Cl    = w[I_W_Cl] ;
    W_M     = w[I_W_Mass] ;
    W_K_D   = -KD_K*grd[I_P_L] ;
    W_K_psi  = -Kpsi_K*grd[I_PSI] ;
    W_K_dif  = - KF_K*grd[I_C_K] ;
  }
    
}



double* Fluxes(Element_t *el,double *grd)
{
  double *va = Element_GetExplicitTerm(el) ;
  double *w  = ComponentFluxes ;
  
  /* Gradients */
  double grd_p_g      = grd[I_P_G] ;
  double grd_p_l      = grd[I_P_L] ;
  double grd_h        = grd[I_C_H ] ;
  double grd_oh       = grd[I_C_OH] ;
  double grd_co2      = grd[I_C_CO2] ;
  double grd_hco3     = grd[I_C_HCO3] ;
  double grd_co3      = grd[I_C_CO3 ] ;
  double grd_cahco3   = grd[I_C_CaHCO3] ;
  double grd_ca       = grd[I_C_Ca] ;
  double grd_cah3sio4 = grd[I_C_CaH3SiO4] ;
  double grd_h3sio4   = grd[I_C_H3SiO4] ;
  double grd_h4sio4   = grd[I_C_H4SiO4] ;
  double grd_h2sio4   = grd[I_C_H2SiO4] ;
  double grd_cah2sio4 = grd[I_C_CaH2SiO4] ;
  double grd_caco3aq  = grd[I_C_CaCO3aq] ;
  double grd_caoh     = grd[I_C_CaOH] ;
  double grd_k        = grd[I_C_K] ;
  double grd_cl       = grd[I_C_Cl] ;
  double grd_psi      = grd[I_PSI] ;
    
    /* Fluxes in the gas phase */
  double wg_co2     = - KD_GCO2*grd_p_g ;
  double wg_m     	= - KD_GM*grd_p_g ;
  
    /* Fluxes in the liquid phase */
  double w_co2      = - KD_CO2*grd_p_l 		  - KF_CO2*grd_co2 ; 
  double w_hco3     = - KD_HCO3*grd_p_l 		- KF_HCO3*grd_hco3          - Kpsi_HCO3*grd_psi ;
  double w_co3      = - KD_CO3*grd_p_l 		  - KF_CO3*grd_co3            - Kpsi_CO3*grd_psi ;
  
  double w_ca       = - KD_Ca*grd_p_l 		  - KF_Ca*grd_ca              - Kpsi_Ca*grd_psi ;
  double w_caoh     = - KD_CaOH*grd_p_l 		- KF_CaOH*grd_caoh          - Kpsi_CaOH*grd_psi ;
  
  double w_cahco3   = - KD_CaHCO3*grd_p_l 	- KF_CaHCO3*grd_cahco3      - Kpsi_CaHCO3*grd_psi ;
  double w_caco3aq  = - KD_CaCO3aq*grd_p_l 	- KF_CaCO3aq*grd_caco3aq ;    
  
  double w_h2sio4   = - KD_H2SiO4*grd_p_l 	- KF_H2SiO4*grd_h2sio4      - Kpsi_H2SiO4*grd_psi ;
  double w_h3sio4   = - KD_H3SiO4*grd_p_l 	- KF_H3SiO4*grd_h3sio4      - Kpsi_H3SiO4*grd_psi ;
  double w_h4sio4   = - KD_H4SiO4*grd_p_l 	- KF_H4SiO4*grd_h4sio4 ;
  
  double w_cah2sio4 = - KD_CaH2SiO4*grd_p_l - KF_CaH2SiO4*grd_cah2sio4 ;
  double w_cah3sio4 = - KD_CaH3SiO4*grd_p_l - KF_CaH3SiO4*grd_cah3sio4  - Kpsi_CaH3SiO4*grd_psi ;
  
  double w_k        = - KD_K*grd_p_l			  - KF_K*grd_k                - Kpsi_K*grd_psi ; 
  double w_cl       = - KD_Cl*grd_p_l 		  - KF_Cl*grd_cl              - Kpsi_Cl*grd_psi ; 
    
  /* Changed compared to what Shen J. did */
  double w_m     	  = - KD_M*grd_p_l ;
  /* This is the law originally introduced by Shen */
  /*
  double w_m     	  = - KD_M*grd_p_l 				- M_CO2*KD_GCO2*grd_p_g	- M_H2O*KD_GH2O*grd_p_g \
						- M_H*Kpsi_H*grd_psi			- M_OH*Kpsi_OH*grd_psi	\
						- M_HCO3*Kpsi_HCO3*grd_psi		- M_CO3*Kpsi_CO3*grd_psi	\
						- M_Ca*Kpsi_Ca*grd_psi			- M_CaHCO3*Kpsi_CaHCO3*grd_psi	\
						- M_H3SiO4*Kpsi_H3SiO4*grd_psi	- M_CaH3SiO4*Kpsi_CaH3SiO4*grd_psi	\
						- M_H2SiO4*Kpsi_H2SiO4*grd_psi	- M_CaOH*Kpsi_CaOH*grd_psi	\
						- M_Cl*Kpsi_Cl*grd_psi			- M_K*Kpsi_K*grd_psi \
						- M_H*KF_H*grd_h						- M_OH*KF_OH*grd_oh \
						- M_CO2*KF_CO2*grd_co2 					 \
						- M_HCO3*KF_HCO3*grd_hco3				- M_CO3*KF_CO3*grd_co3 \
						- M_CaHCO3*KF_CaHCO3*grd_cahco3 		- M_Ca*KF_Ca*grd_ca\
						- M_CaH3SiO4*KF_CaH3SiO4*grd_cah3sio4	- M_H3SiO4*KF_H3SiO4*grd_h3sio4\
						- M_H2SiO4*KF_H2SiO4*grd_h2sio4\
						- M_CaOH*KF_CaOH*grd_caoh				- M_H4SiO4*KF_H4SiO4*grd_h4sio4\
						- M_CaH2SiO4*KF_CaH2SiO4*grd_cah2sio4	- M_CaCO3aq*KF_CaCO3aq*grd_caco3aq\
						- M_K*KF_K*grd_k 						- M_Cl*KF_Cl*grd_cl\
						+ M_H2O*(Kpsi_H*grd_psi + KF_H*grd_h)*v_h/v_h2o\
						+ M_H2O*(Kpsi_OH*grd_psi + KF_OH*grd_oh)*v_oh/v_h2o\
						+ M_H2O*(Kpsi_CaHCO3*grd_psi + KF_CaHCO3*grd_cahco3)*v_cahco3/v_h2o\
						+ M_H2O*(Kpsi_CO3*grd_psi + KF_CO3*grd_co3)*v_co3/v_h2o\
						+ M_H2O*(Kpsi_Ca*grd_psi + KF_Ca*grd_ca)*v_ca/v_h2o\
						+ M_H2O*(Kpsi_CaH3SiO4*grd_psi + KF_CaH3SiO4*grd_cah3sio4)*v_cah3sio4/v_h2o\
						+ M_H2O*(Kpsi_H3SiO4*grd_psi + KF_H3SiO4*grd_h3sio4)*v_h3sio4/v_h2o\
						+ M_H2O*(Kpsi_H2SiO4*grd_psi + KF_H2SiO4*grd_h2sio4)*v_h2sio4/v_h2o\
						+ M_H2O*(Kpsi_CaOH*grd_psi + KF_CaOH*grd_caoh)*v_caoh/v_h2o\
						+ M_H2O*(KF_CaH2SiO4*grd_cah2sio4)*v_cah2sio4/v_h2o\
						+ M_H2O*(KF_H4SiO4*grd_h4sio4)*v_h4sio4/v_h2o\
						+ M_H2O*(KF_CaCO3aq*grd_caco3aq)*v_caco3aq/v_h2o\
						+ M_H2O*(Kpsi_K*grd_psi + KF_K*grd_k)*v_k/v_h2o\
						+ M_H2O*(Kpsi_Cl*grd_psi + KF_Cl*grd_cl)*v_cl/v_h2o\
						+ M_H2O*(Kpsi_HCO3*grd_psi + KF_HCO3*grd_hco3)*v_hco3/v_h2o\
						+ M_H2O*KF_CO2*grd_co2*v_co2/v_h2o ;
  */
						
  double w_q        = - z_h*KF_H*grd_h		                  \
                      - z_oh*KF_OH*grd_oh		                \
                      - z_hco3*KF_HCO3*grd_hco3             \
                      - z_co3*KF_CO3*grd_co3		            \
                      - z_ca*KF_Ca*grd_ca		                \
                      - z_cahco3*KF_CaHCO3*grd_cahco3	      \
                      - z_h3sio4*KF_H3SiO4*grd_h3sio4	      \
                      - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                      - z_h2sio4*KF_H2SiO4*grd_h2sio4       \
                      - z_caoh*KF_CaOH*grd_caoh             \
                      - z_k*KF_K*grd_k                      \
                      - z_cl*KF_Cl*grd_cl                   \
                      - Kpsi_q*grd_psi ;

  /* Back up */

  w[I_W_C    ] = wg_co2 + w_co2 + w_hco3 + w_co3 + w_cahco3 + w_caco3aq ;
  w[I_W_Ca   ] = w_ca + w_cahco3 + w_cah3sio4 + w_caco3aq + w_caoh + w_cah2sio4 ;
  w[I_W_Si   ] = w_h3sio4 + w_h4sio4 + w_cah3sio4 + w_cah2sio4 + w_h2sio4 ;
  w[I_W_Q    ] = w_q ;
  w[I_W_K    ] = w_k ;
  w[I_W_Cl   ] = w_cl ;
  w[I_W_Mass ] = wg_m + w_m ;
    
  return(w) ;
}



int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  double x1    = Element_GetNodeCoordinate(el,1)[0] ;
  double x0    = Element_GetNodeCoordinate(el,0)[0] ;
  double dij   = x1 - x0 ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  

  /* Initialization */
  for(i = 0 ; i < nn*nn*NEQ*NEQ ; i++) c[i] = 0. ;


  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Donnees
  */
  phi0      = GetProperty("porosite") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_Si") ;
  c_co2_eq  = GetProperty("C_CO2_eq") ;
  t_ch      = GetProperty("T_CH") ;
  t_cc      = GetProperty("T_CC") ;
  temperature 	= GetProperty("temperature") ;
  
  
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
    
    dxi[U_P_G    ] = P_Gn(i)*1.e-6 ;
    dxi[U_P_L    ] = P_Ln(i)*1.e-6 ;
    dxi[U_ZN_Si_S] = ((x[U_ZN_Si_S] > ZN_Si_Sn(i)) ? 1 : -1)*1.e-4 ; 
    dxi[U_ZN_Ca_S] = ((x[U_ZN_Ca_S] > ZN_Ca_Sn(i)) ? 1 : -1)*1.e-6 ;
    dxi[U_C_K    ] = 0.4e-4 ;
    dxi[U_C_Cl   ] = 0.4e-4 ;
    dxi[U_PSI    ] = 1. ;
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      cii[E_C*NEQ    + k] = dx[I_N_C] ;
      cii[E_Ca*NEQ   + k] = dx[I_N_Ca] ;
      cii[E_Si*NEQ   + k] = dx[I_N_Si] ;
      cii[E_K*NEQ    + k] = dx[I_N_K] ;
      cii[E_Cl*NEQ   + k] = dx[I_N_Cl] ;
      cii[E_mass*NEQ + k] = dx[I_Mass] ;
      
      cij[E_C*NEQ    + k] = - dtdij*dw[I_W_C] ;
      cij[E_Ca*NEQ   + k] = - dtdij*dw[I_W_Ca] ;
      cij[E_Si*NEQ   + k] = - dtdij*dw[I_W_Si] ;
      cij[E_K*NEQ    + k] = - dtdij*dw[I_W_K] ;
      cij[E_Cl*NEQ   + k] = - dtdij*dw[I_W_Cl] ;
      cij[E_mass*NEQ + k] = - dtdij*dw[I_W_Mass] ;
      cij[E_q*NEQ    + k] = - dtdij*dw[I_W_Q] ;
    }
  }

  return(dec) ;
}





double concentration_oh(double c_co2,double zn_si_s,double zn_ca_s,double c_k,double c_cl,Element_t *el)
/* Solve for electroneutrality : SUM(z_i c_i) = 0
   as root of ax^4 + bx^3 + cx^2 + dx + e = 0 */
{
  /* c_co2, zn_si_s and zn_ca_s are considered as constant */
  double z_co2      = c_co2/c_co2_eq ;
  /* Ion activity products are constant as well */
  double s_cc       = SaturationDegreeOfCC(z_co2,zn_ca_s) ;
  double Q_CC       = IonActivityProductOfCC(s_cc) ;
  double s_ch       = SaturationDegreeOfCH(z_co2,zn_ca_s) ;
  double s_sh       = SaturationDegreeOfSH(s_ch,zn_si_s) ;
  double Q_SH       = IonActivityProductOfSH(s_sh) ;
  
  double c_h4sio4   = Q_SH ;
  /* rappel des expressions c_i = A_i*(c_h)**n   : n
     c_h        = K_h2o/c_oh                     : +1
     c_h2co3    = K_h2co3*c_co2                  :  0
     c_hco3     = K_hco3*c_h2co3/c_h             : -1
     c_co3      = K_co3*c_hco3/c_h               : -2
     c_ca       = Q_CC/c_co3                     : +2
     c_cahco3   = K_cahco3*c_ca*c_hco3           : +1
     c_h4sio4   = Q_SH                           :  0
     c_h3sio4   = c_h4sio4/(K_h4sio4*c_h)        : -1
     c_cah3sio4 = K_cah3sio4*c_ca*c_h3sio4       : +1
     c_caco3aq  = K_caco3aq*c_ca*c_co3 ;         :  0      
     c_h2sio4   = K_h2sio4*c_h3sio4*c_oh ;       : -2       
     c_cah2sio4 = K_cah2sio4*c_h2sio4*c_ca ;     :  0      
     c_caoh     = K_caoh*c_ca*c_oh ;             : +1       
  */

  double A_hco3     = K_hco3*c_co2 ;
  double A_co3      = K_co3*A_hco3 ;
  double A_ca       = Q_CC/A_co3 ;
  double A_cahco3   = K_cahco3*A_ca*A_hco3 ;
  double A_h3sio4   = c_h4sio4*K_h3sio4 ;
  double A_cah3sio4 = K_cah3sio4*A_ca*A_h3sio4 ;
  double A_h2sio4   = K_h2sio4*A_h3sio4*K_h2o ;
  double A_caoh     = K_caoh*A_ca*K_h2o ;

  double a = z_ca*A_ca ;
  double b = z_h + z_cahco3*A_cahco3 + z_cah3sio4*A_cah3sio4 + z_caoh*A_caoh ;
  double c = z_k*c_k + z_cl*c_cl ;
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



/* Models for scCO2 and scCO2-H2O mixture */

double (solubilityCO2)(Element_t *el,double p_g,double p_l,double TK)
/* Solubility of CO2 in liquid phase, after Nicolas Spycher 2003 */
{
  double pb     = p_g/1.e5 ;   /* unit of bar */
  double f_co2  = FugacityCoefficientOfCO2InGas(p_g,TK) ;
  double f_h2o  = FugacityCoefficientOfH2OInGas(p_g,TK) ;
  double k_h2o  = EquilibriumLiquidGasOfH2O(p_l,TK) ;
  double k_co2  = EquilibriumLiquidGasOfCO2(p_l,TK) ;
  double A      = k_h2o/(pb*f_h2o) ;
  double B      = (pb*f_co2)/k_co2 ;
  double a_h2o  = (55.508 - B)/(55.508 - A*B) ;
  double y_h2o  = MIN(A*a_h2o,1) ;
  double c_co2  = B*(1 - y_h2o) ;
  return(c_co2) ;
}

double (solubilityCO2b)(Element_t *el,double p_g,double p_l,double TK)
/* CO2 concentration in liquid */
{
  double p0		= 1.e5 ;
  double p1		= 4.e5 ;
  double c_co2 ;
  
  if (p_g < p0) {
    double c_co20 = solubilityCO2(el,(p0+53000),p_l,TK) ; /*fit for p_g=1atm*/
    c_co2  = c_co20*p_g/p0 ;
  } else if (p_g < p1) {
    double c_co20 = solubilityCO2(el,(p0+53000),p_l,TK) ;
    double c_co21 = solubilityCO2(el,p1,p_l,TK) ;
    c_co2  = c_co20 + (c_co21 - c_co20)*(p_g - p0)/(p1 - p0) ;
  } else {
    c_co2 	= solubilityCO2(el,p_g,p_l,TK) ;
  }
  
  return(c_co2) ;
}


double (molefractionH2O)(Element_t *el,double p_g,double p_l,double TK)
/* Mole fraction of H2O in gas phase, after Nicolas Spycher 2003 */
{
  double pb     = p_g/1.e5 ;   /* unit of bar */
  double f_co2  = FugacityCoefficientOfCO2InGas(p_g,TK) ;
  double f_h2o  = FugacityCoefficientOfH2OInGas(p_g,TK) ;
  double k_h2o  = EquilibriumLiquidGasOfH2O(p_l,TK) ;
  double k_co2  = EquilibriumLiquidGasOfCO2(p_l,TK) ;
  double A      = k_h2o/(pb*f_h2o) ;
  double B      = (pb*f_co2)/k_co2 ;
  double a_h2o  = (55.508 - B)/(55.508 - A*B) ;
  double y_h2o  = MIN(A*a_h2o,1) ;
  return(y_h2o) ;
}

double equilibriumCO2(double pl,double TK)
{
  double T0        = 273. ;
  double TC        = TK - T0 ;
  double TC2       = TC*TC ;
  double R         = 83.14472 ;   /*cm3*bar*K-1*mol-1*/
  double RTK       = R*TK ;
  double p0        = 1. ;
  double pb        = MAX(pl/1.e5,p0) ;
  double a         = 1.189 ;
  double b         = 1.304e-2 ;
  double c         = -5.446e-5 ;
  double v0_co2    = 32.6 ; /*cm³/mol*/
  double logk0     = a + b*TC + c*TC2 ;
  double k0        = pow(10,logk0) ;
  double k_co2     = k0*exp((pb - p0)*v0_co2/RTK) ;
  return(k_co2) ;
}

double equilibriumH2O(double pl,double TK)
{
  double T0        = 273. ;
  double TC        = TK - T0 ;
  double TC2       = TC*TC ;
  double TC3       = TC*TC2 ;
  double R         = 83.14472 ;   /*cm3*bar*K-1*mol-1*/
  double RTK       = R*TK ;
  double p0        = 1. ;
  double pb        = MAX(pl/1.e5,p0) ;
  double a         = -2.209 ;
  double b         = 3.097e-2 ;
  double c         = -1.098e-4 ;
  double d         = 2.048e-7 ;
  double v0_h2o    = 18.1 ; /* cm3/mol */
  double logk0     = a + b*TC + c*TC2 + d*TC3 ;
  double k0        = pow(10,logk0) ;
  double k_h2o     = k0*exp((pb - p0)*v0_h2o/RTK) ;
  return(k_h2o) ;
}

double fugacityCO2(Element_t *el,double pg,double TK)
{
  double pb        = pg/1.e5 ;    /* unit in bar */
  double R         = 83.14472 ;   /* cm3*bar/mol/K */
  double RTK       = R*TK ;       /* cm3*bar/mol */
  double T05       = sqrt(TK) ;
  double A_CO2     = 7.54e7 - 4.13e4*TK ; /* bar*cm6*K0.5/mol^2 */
  double B_CO2     = 27.80 ;              /* cm3/mol */
  double C_CO2     = A_CO2/(RTK*T05*B_CO2) ;
  double V         = 1.e3/MolarDensityOfCO2(pg,TK) ; /* cm3/mol */
  double zg        = pb*V/RTK ; /* Compression factor */
  double lco2      = - log(zg*(1 - B_CO2/V)) + B_CO2/(V - B_CO2) \
			               - C_CO2*(log(1 + B_CO2/V) + B_CO2/(V + B_CO2)) ;
  return(exp(lco2)) ;
}

double fugacityH2O(Element_t *el,double pg,double TK)
{
  double pb        = pg/1.e5 ;    /* unit in bar */
  double R         = 83.14472 ;   /* cm3*bar/mol/K */
  double RTK       = R*TK ;       /* cm3*bar/mol */
  double T05       = sqrt(TK) ;
  double A_CO2     = 7.54e7 - 4.13e4*TK ; /* bar*cm6*K0.5/mol^2 */
  double B_CO2     = 27.80 ;              /* cm3/mol */
  double B_H2O     = 18.18 ;              /* cm3/mol */
  double A_H2OCO2  = 7.89e7 ;             /* bar*cm6*K0.5*mol^-2 */
  double C_H2O     = A_CO2*B_H2O/(RTK*T05*B_CO2*B_CO2) ;
  double V         = 1.e3/MolarDensityOfCO2(pg,TK) ; /* cm3/mol */
  double zg        = pb*V/RTK ; /* Compression factor */
  double lh2o      = - log(zg*(1 - B_CO2/V)) + B_H2O/(V - B_CO2) \
                     - 2*A_H2OCO2/(RTK*T05*B_CO2)*log(1 + B_CO2/V) \
			               - C_H2O*(- log(1 + B_CO2/V) + B_CO2/(V + B_CO2)) ;
  return(exp(lh2o)) ;
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
  
  Math_PolishPolynomialEquationRoot(y,4,&x,tol*x,20) ;
  
  return(x) ;
}
