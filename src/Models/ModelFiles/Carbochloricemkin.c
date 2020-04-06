/* General features of the model:
 * Curves for CSH:
 *   - C/S ratio
 *   - H/S ratio
 *   - Molar Volume
 * Alkalis (as sodium and potassium compounds)
 * Dissolution kinetics for CH based on spherical crystal 
 * coated by a calcite layer.
 * Dissolution and continuous decalcification of CSH
 * Precipitation/Dissolution of CC (possibly with kinetics)
 * Chloride ingress
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* The Finite Volume Method */
#include "FVM.h"

/* Cement chemistry */
#include "HardenedCementChemistry.h"
#include "CementSolutionDiffusion.h"

#define TITLE   "Carbonation of CBM (2019)"
#define AUTHORS "Dangla and many others"

#include "PredefinedMethods.h"




/* Indices of equations/unknowns */
enum {
  E_Carbon    ,
  E_charge    ,
  E_mass      ,
  E_Calcium   ,
  E_Silicon   ,
  E_Sodium    ,
  E_Potassium ,
  /* Uncomment/comment the next two lines to consider/suppress electroneutrality */
  E_eneutral  ,
  #define E_eneutral E_eneutral 
  /* Uncomment/Comment the next two lines to consider/suppress chlorine */
  E_Chlorine  ,
  #define E_Chlorine E_Chlorine
  /* Uncomment/Comment the next two lines to consider/suppress air */
  //E_Air    ,
  //#define E_air E_air
  E_Last
} ;


/* Nb of equations */
#define NbOfEquations      E_Last
#define NEQ                NbOfEquations






/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Generic names of nodal unknowns */
#define U_Carbon(n)     (UNKNOWN(n,E_Carbon))
#define Un_Carbon(n)    (UNKNOWNn(n,E_Carbon))

#define U_charge(n)     (UNKNOWN(n,E_charge))
#define Un_charge(n)    (UNKNOWNn(n,E_charge))

#define U_mass(n)       (UNKNOWN(n,E_mass))
#define Un_mass(n)      (UNKNOWNn(n,E_mass))

#define U_Calcium(n)    (UNKNOWN(n,E_Calcium))
#define Un_Calcium(n)   (UNKNOWNn(n,E_Calcium))

#define U_Silicon(n)    (UNKNOWN(n,E_Silicon))
#define Un_Silicon(n)   (UNKNOWNn(n,E_Silicon))

#define U_Sodium(n)     (UNKNOWN(n,E_Sodium))
#define Un_Sodium(n)    (UNKNOWNn(n,E_Sodium))

#define U_Potassium(n)  (UNKNOWN(n,E_Potassium))
#define Un_Potassium(n) (UNKNOWNn(n,E_Potassium))

#define U_eneutral(n)   (UNKNOWN(n,E_eneutral))
#define Un_eneutral(n)  (UNKNOWNn(n,E_eneutral))

#define U_Chlorine(n)   (UNKNOWN(n,E_Chlorine))
#define Un_Chlorine(n)  (UNKNOWNn(n,E_Chlorine))

#define U_Air(n)        (UNKNOWN(n,E_Air))
#define Un_Air(n)       (UNKNOWNn(n,E_Air))




/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */

/* Carbon: unknown either C or logC */
//#define U_C_CO2     U_Carbon
#define U_LogC_CO2  U_Carbon

/* Mass: the liquid pressure */
#define U_P_L       U_mass

/* Calcium:
 * - U_ZN_Ca_S: dissolution kinetics of CH and Cc at equilibrium 
 * - U_LogS_CH: dissolution kinetics of CH and precipitation kinetics of Cc */
//#define U_ZN_Ca_S   U_Calcium
#define U_LogS_CH     U_Calcium

/* Silicon: */
#define U_ZN_Si_S   U_Silicon

/* charge: */
#define U_PSI       U_charge

/* Sodium: unknown either C or logC */
//#define U_C_Na      U_Sodium
#define U_LogC_Na   U_Sodium

/* Potassium: unknown either C or logC */
//#define U_C_K       U_Potassium
#define U_LogC_K    U_Potassium

/* Electroneutrality: unknown either C_OH, logC_OH or Z_OH = C_H - C_OH */
//#define U_C_OH      U_eneutral
#define U_LogC_OH   U_eneutral
//#define U_Z_OH      U_eneutral

/* Chlorine: unknown either C or logC */
#define U_LogC_Cl   U_Chlorine
//#define U_C_Cl      U_Chlorine

/* Air: the gas pressure */
#define U_P_G       U_Air





/* Names of nodal unknowns */
#if defined (U_LogC_CO2) && !defined (U_C_CO2)
  #define LogC_CO2(n)   U_Carbon(n)
  #define LogC_CO2n(n)  Un_Carbon(n)
  #define C_CO2(n)      (pow(10,LogC_CO2(n)))
  #define C_CO2n(n)     (pow(10,LogC_CO2n(n)))
#elif defined (U_C_CO2) && !defined (U_LogC_CO2)
  #define C_CO2(n)      U_Carbon(n)
  #define C_CO2n(n)     Un_Carbon(n)
  #define LogC_CO2(n)   (log10(C_CO2(n)))
  #define LogC_CO2n(n)  (log10(C_CO2n(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif


#define P_L(n)        U_mass(n)
#define P_Ln(n)       Un_mass(n)

#define PSI(n)        U_charge(n)
#define PSIn(n)       Un_charge(n)

#if defined (U_LogC_Na) && !defined (U_C_Na)
  #define LogC_Na(n)    U_Sodium(n)
  #define LogC_Nan(n)   Un_Sodium(n)
  #define C_Na(n)       (pow(10,LogC_Na(n)))
  #define C_Nan(n)      (pow(10,LogC_Nan(n)))
#elif defined (U_C_Na) && !defined (U_LogC_Na)
  #define C_Na(n)       U_Sodium(n)
  #define C_Nan(n)      Un_Sodium(n)
  #define LogC_Na(n)    (log10(C_Na(n)))
  #define LogC_Nan(n)   (log10(C_Nan(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif

#if defined (U_LogC_K) && !defined (U_C_K)
  #define LogC_K(n)     U_Potassium(n)
  #define LogC_Kn(n)    Un_Potassium(n)
  #define C_K(n)        (pow(10,LogC_K(n)))
  #define C_Kn(n)       (pow(10,LogC_Kn(n)))
#elif defined (U_C_K) && !defined (U_LogC_K)
  #define C_K(n)        U_Potassium(n)
  #define C_Kn(n)       Un_Potassium(n)
  #define LogC_K(n)     (log10(C_K(n)))
  #define LogC_Kn(n)    (log10(C_Kn(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif

#ifdef E_eneutral
  #if defined (U_LogC_OH) && !defined (U_C_OH) && !defined (U_Z_OH)
    #define LogC_OH(n)   U_eneutral(n)
    #define LogC_OHn(n)  Un_eneutral(n)
    #define C_OH(n)      (pow(10,LogC_OH(n)))
    #define C_OHn(n)     (pow(10,LogC_OHn(n)))
  #elif defined (U_C_OH) && !defined (U_LogC_OH) && !defined (U_Z_OH)
    #define C_OH(n)      U_eneutral(n)
    #define C_OHn(n)     Un_eneutral(n)
    #define LogC_OH(n)   (log10(C_OH(n)))
    #define LogC_OHn(n)  (log10(C_OHn(n)))
  #elif defined (U_Z_OH) && !defined (U_LogC_OH) && !defined (U_C_OH)
    #define Z_OH(n)      U_eneutral(n)
    #define Z_OHn(n)     Un_eneutral(n)
    #define Z_OHDefinition(c_oh,c_h)   ((c_h) - (c_oh))
    #define C_OHDefinition(z_oh)       (0.5 * (sqrt((z_oh)*(z_oh) + 4*K_w) - (z_oh)))
    #define C_OH(n)      C_OHDefinition(Z_OH(n))
    #define C_OHn(n)     C_OHDefinition(Z_OHn(n))
    #define LogC_OH(n)   (log10(C_OH(n)))
    #define LogC_OHn(n)  (log10(C_OHn(n)))
    #define dLogC_OHdZ_OH(z_oh)    (-1/(Ln10 * sqrt((z_oh)*(z_oh) + 4*K_w)))
  #else
    #error "Ambiguous or undefined unknown"
  #endif
#endif



#ifdef E_Chlorine
  #if defined (U_LogC_Cl) && !defined (U_C_Cl)
    #define LogC_Cl(n)     U_Chlorine(n)
    #define LogC_Cln(n)    Un_Chlorine(n)
    #define C_Cl(n)        (pow(10,LogC_Cl(n)))
    #define C_Cln(n)       (pow(10,LogC_Cln(n)))
  #elif defined (U_C_Cl) && !defined (U_LogC_Cl)
    #define C_Cl(n)        U_Chlorine(n)
    #define C_Cln(n)       Un_Chlorine(n)
    #define LogC_Cl(n)     (log10(C_Cl(n)))
    #define LogC_Cln(n)    (log10(C_Cln(n)))
  #else
    #error "Ambiguous or undefined unknown"
  #endif
#endif



#ifdef E_Air
  #define P_G(n)        U_Air(n)
  #define P_Gn(n)       Un_Air(n)
#endif



/* Nb of nodes (el must be used below) */
#define NbOfNodes               Element_GetNbOfNodes(el)


/* Nb of terms */
#define NbOfExplicitTerms       ((12 + CementSolutionDiffusion_NbOfConcentrations)*NbOfNodes)
#define NbOfImplicitTerms       (9*NbOfNodes*NbOfNodes + 3*NbOfNodes)
#define NbOfConstantTerms       (2)


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NbOfNodes + (j))])

#define NW_C          (f)
#define NW_Cn         (f_n)
#define N_C(i)        MassAndFlux(NW_C,i,i)
#define N_Cn(i)       MassAndFlux(NW_Cn,i,i)
#define W_C(i,j)      MassAndFlux(NW_C,i,j)

#define NW_q          (f   + NbOfNodes*NbOfNodes)
#define NW_qn         (f_n + NbOfNodes*NbOfNodes)
#define N_q(i)        MassAndFlux(NW_q,i,i)
#define N_qn(i)       MassAndFlux(NW_qn,i,i)
#define W_q(i,j)      MassAndFlux(NW_q,i,j)

#define MW_tot        (f   + 2*NbOfNodes*NbOfNodes)
#define MW_totn       (f_n + 2*NbOfNodes*NbOfNodes)
#define M_tot(i)      MassAndFlux(MW_tot,i,i)
#define M_totn(i)     MassAndFlux(MW_totn,i,i)
#define W_tot(i,j)    MassAndFlux(MW_tot,i,j)

#define NW_Ca         (f   + 3*NbOfNodes*NbOfNodes)
#define NW_Can        (f_n + 3*NbOfNodes*NbOfNodes)
#define N_Ca(i)       MassAndFlux(NW_Ca,i,i)
#define N_Can(i)      MassAndFlux(NW_Can,i,i)
#define W_Ca(i,j)     MassAndFlux(NW_Ca,i,j)

#define NW_Na         (f   + 4*NbOfNodes*NbOfNodes)
#define NW_Nan        (f_n + 4*NbOfNodes*NbOfNodes)
#define N_Na(i)       MassAndFlux(NW_Na,i,i)
#define N_Nan(i)      MassAndFlux(NW_Nan,i,i)
#define W_Na(i,j)     MassAndFlux(NW_Na,i,j)

#define NW_K          (f   + 5*NbOfNodes*NbOfNodes)
#define NW_Kn         (f_n + 5*NbOfNodes*NbOfNodes)
#define N_K(i)        MassAndFlux(NW_K,i,i)
#define N_Kn(i)       MassAndFlux(NW_Kn,i,i)
#define W_K(i,j)      MassAndFlux(NW_K,i,j)

#define NW_Si         (f   + 6*NbOfNodes*NbOfNodes)
#define NW_Sin        (f_n + 6*NbOfNodes*NbOfNodes)
#define N_Si(i)       MassAndFlux(NW_Si,i,i)
#define N_Sin(i)      MassAndFlux(NW_Sin,i,i)
#define W_Si(i,j)     MassAndFlux(NW_Si,i,j)

#define NW_Cl         (f   + 7*NbOfNodes*NbOfNodes)
#define NW_Cln        (f_n + 7*NbOfNodes*NbOfNodes)
#define N_Cl(i)       MassAndFlux(NW_Cl,i,i)
#define N_Cln(i)      MassAndFlux(NW_Cln,i,i)
#define W_Cl(i,j)     MassAndFlux(NW_Cl,i,j)

#define MW_Air         (f   + 8*NbOfNodes*NbOfNodes)
#define MW_Airn        (f_n + 8*NbOfNodes*NbOfNodes)
#define M_Air(i)       MassAndFlux(MW_Air,i,i)
#define M_Airn(i)      MassAndFlux(MW_Airn,i,i)
#define W_Air(i,j)     MassAndFlux(MW_Air,i,j)

#define N_CH(i)       (f   + 9*NbOfNodes*NbOfNodes)[i]
#define N_CHn(i)      (f_n + 9*NbOfNodes*NbOfNodes)[i]

#define N_CC(i)       (f   + 9*NbOfNodes*NbOfNodes + NbOfNodes)[i]
#define N_CCn(i)      (f_n + 9*NbOfNodes*NbOfNodes + NbOfNodes)[i]

#ifndef E_eneutral
  #define C_OH(i)       (f   + 9*NbOfNodes*NbOfNodes + 2*NbOfNodes)[i]
  #define C_OHn(i)      (f_n + 9*NbOfNodes*NbOfNodes + 2*NbOfNodes)[i]
  
  #define LogC_OH(n)    (log10(C_OH(n)))
  #define LogC_OHn(n)   (log10(C_OHn(n)))
#endif




/* Names used for explicit terms */
#define TransferCoefficient(va,n)  ((va) + (n)*NbOfNodes)

#define KF_CO2          TransferCoefficient(va,0)

#define KD_L            TransferCoefficient(va,1)

#define KC_C_L          TransferCoefficient(va,2)
#define KC_Ca_L         TransferCoefficient(va,3)
#define KC_Na_L         TransferCoefficient(va,4)
#define KC_K_L          TransferCoefficient(va,5)
#define KC_Si_L         TransferCoefficient(va,6)
#define KC_Cl_L         TransferCoefficient(va,7)

#define KF_H2O          TransferCoefficient(va,8)

#define TORTUOSITY      TransferCoefficient(va,9)

#define KC_CO2_G        TransferCoefficient(va,10)
#define KC_H2O_G        TransferCoefficient(va,11)

#define KD_G            TransferCoefficient(va,12)

#define CONCENTRATION(i) \
        (TransferCoefficient(va,12) + (i)*CementSolutionDiffusion_NbOfConcentrations)



/* Names used for constant terms */
#define V_S0(n)         (v0[(0+n)])





/* Math constants */
#define Ln10      Math_Ln10




/* Units
 * ----- */
#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define dm    (0.1*InternationalSystemOfUnits_OneMeter)
#define cm    (0.01*InternationalSystemOfUnits_OneMeter)
#define dm2   (dm*dm)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)
#define MPa   (1.e6*InternationalSystemOfUnits_OnePascal)
#define GPa   (1.e3*MPa)
#define mol   InternationalSystemOfUnits_OneMole
#define sec   InternationalSystemOfUnits_OneSecond
#define kg    InternationalSystemOfUnits_OneKilogram
#define gr    (0.001*kg)


#define TEMPERATURE  (298)




#include "MolarMassOfMolecule.h"


/* Water property
 * -------------- */
 /* Molar mass */
#define M_H2O          MolarMassOfMolecule(H2O)
/* Molar volume of liquid water */
#define V_H2O          (18 * cm3)
/* Mass density */
#define MassDensityOfWaterVapor(p_v)   (M_H2O*(p_v)/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(p_l)          (exp(V_H2O/RT*((p_l) - p_l0)))
#define VaporPressure(p_l)             (p_v0*RelativeHumidity(p_l))
//#define LiquidPressure(hr)             (p_l0 + RT/V_H2O*(log(hr)))



/* Liquid properties
 * ----------------- */
//#undef  HardenedCementChemistry_GetLiquidMassDensity
//#define HardenedCementChemistry_GetLiquidMassDensity(hcc)    (rho_l0)



/* Dry air properties
 * ------------------ */
/* Molar mass */
#define M_AIR          (28.8 * gr)
/* Mass density */
#define MassDensityOfDryAir(p_a)       (M_AIR*(p_a)/RT)


/* CO2 gas properties
 * ------------------ */
#define M_CO2          MolarMassOfMolecule(CO2)
/* Partial pressure of CO2 */
#define PartialPressureOfCO2(rho_co2)   ((rho_co2)*RT)
/* Henry's laww constant for the solubility of CO2 gas */
#define k_h           (0.9983046)                /* CO2(g) = CO2(aq) (T = 293K)*/




/* Material Properties
 * ------------------- */
#define SATURATION_CURVE                 (saturationcurve)
#if defined (E_Air)
  #define SaturationDegree(p)              (saturationdegree(p,p_c3,SATURATION_CURVE))
#else
  #define SaturationDegree(p)              (Curve_ComputeValue(SATURATION_CURVE,p))
#endif
#define RELATIVEPERMLIQ_CURVE            (relativepermliqcurve)
#define RelativePermeabilityToLiquid(s)  (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,s))
#ifdef E_Air
  #define RELATIVEPERMGAS_CURVE            (relativepermgascurve)
  #define RelativePermeabilityToGas(s)     (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,s))
#else
  #define RelativePermeabilityToGas(s)     (1)
#endif
//#define TortuosityToLiquid               TortuosityToLiquid_Xie
#define TortuosityToLiquid               TortuosityToLiquid_OhJang
#define PermeabilityCoefficient          PermeabilityCoefficient_VermaPruess




/* Calcium Silicate Hydrate Properties (C-S-H)
 * ------------------------------------------- */
#define M_CaO          MolarMassOfMolecule(CaO)
#define M_SiO2         MolarMassOfMolecule(SiO2)
#define MolarMassOfCSH(x,z)     (M_CaO*(x) + M_SiO2 + M_H2O*(z))
#define MOLARVOLUMEOFCSH_CURVE           (molarvolumeofcshcurve)
#define MolarVolumeOfCSH(x_ch)           (Curve_ComputeValue(MOLARVOLUMEOFCSH_CURVE,x_ch))
//#define MolarVolumeOfCSH(s_ch)           (Curve_ComputeValue(Element_GetCurve(el) + 4,s_ch))
#define V_CSH         (78 * cm3)
#define V_SH          (43 * cm3)
//#define MolarVolumeOfCSH(x)    ((x)/1.7*V_CSH + (1 - (x)/1.7)*V_SH)
/* Below is how to manage dissolution/precipitation */
/* Definition of U_Silicon = ZN_Si_S in the next lines:
 * ZN_Si_S = N/N0 + log(S_CSH)
 * with N = silicon content in CSH
 * and S_CSH = saturation index of CSH */
/* Log of saturation index */
#define Log10SaturationIndexOfCSH(zn_si_s)  MIN(zn_si_s,0.)
#define SiliconContentInCSH(zn_si_s)        (n_si_ref*MAX(zn_si_s,0.))
#define CSHSolidContent(zn_si_s)            SiliconContentInCSH(zn_si_s)



/* Calcium Hydroxide (Portlandite) Properties (CH)
 * ----------------------------------------------- */
#define M_CaOH2        MolarMassOfMolecule(CaO2H2)
/* Molar volume of CH solid */
#define V_CH           (33 * cm3)
/* Below is how to manage dissolution/precipitation kinetics */
#if defined (U_ZN_Ca_S)
  /* Definition of U_Calcium = ZN_Ca_S in the next lines: 
   * ZN_Ca_S = N/N0 + log(S) 
   * with N = calcium content in CH and CC  (>= 0)
   * and  S = saturation index of CcH ie max(log(S_CH),log(S_Cc)) (<= 0) */
  /* Log of saturation index of CcH = CaO-CO2-H2O  */
  #define Log10SaturationIndexOfCcH(zn_ca_s)   MIN(zn_ca_s,0.)
  /* Calcium solid content in CcH ie in CH and Cc */
  #define CalciumContentInCcH(zn_ca_s)        (n_ca_ref*MAX(zn_ca_s,0.))
  /* Initial solid content of CH */
  #define InitialCHSolidContent(zn_ca_s,s_ch,s_cc) \
          (((s_cc) > (s_ch)) ? 0 : CalciumContentInCcH(zn_ca_s))
  /* Current solid content of CH */
  #define CHSolidContent(zn_ca_s,n_chn,n_ccn,s_ch,s_cc,dt) \
          (((s_cc) > (s_ch)) ? CHSolidContent_kin1(n_chn,s_ch,dt) : \
          CalciumContentInCcH(zn_ca_s))
#elif defined (U_LogS_CH)
  /* Definition of U_Calcium = LogS_CH in the next lines:
   * LogS_CH = log(S_CH)
   * with S_CH = saturation index of CH */
  /* Initial solid content of CH */
  #define InitialCHSolidContent(logs_ch,s_ch,s_cc)     (n_ca_ref)
  /* Log of saturation index of CH */
  #define Log10SaturationIndexOfCH(logs_ch)            (logs_ch)
  /* Current solid content of CH */
  #define CHSolidContent(logs_ch,n_chn,n_ccn,s_ch,s_cc,dt) \
          CHSolidContent_kin1(n_chn,s_ch,dt)
#endif



/* Calcium Carbonate (Calcite) Properties (CC)
 * ------------------------------------------- */
#define M_CaCO3        MolarMassOfMolecule(CaCO3)
/* Molar volume of CC */
#define V_CC           (37 * cm3)
/* Below is how to manage dissolution/precipitation kinetics */
#if defined (U_ZN_Ca_S)
  /* See above the definition of ZN_Ca_S */
  /* Initial solid content of Cc */
  #define InitialCCSolidContent(zn_ca_s,s_ch,s_cc) \
          (((s_cc) > (s_ch)) ? CalciumContentInCcH(zn_ca_s) : 0)
  /* Current solid content of Cc */
  #define CCSolidContent(zn_ca_s,n_chn,n_ccn,s_ch,s_cc,dt) \
          (CalciumContentInCcH(zn_ca_s) - CHSolidContent(zn_ca_s,n_chn,n_ccn,s_ch,s_cc,dt))
#elif defined (U_LogS_CH)
  /* See above the definition of LogS_CH */
  /* Initial solid content of Cc */
  #define InitialCCSolidContent(logs_ch,s_ch,s_cc)   (0)
  /* Current solid content */
  #define CCSolidContent_kin(n,s,dt)        MAX((n + dt*rate_cc*(s - 1)),0.)
  #define CCSolidContent(logs_ch,n_chn,n_ccn,s_ch,s_cc,dt) \
          CCSolidContent_kin(n_ccn,s_cc,dt)
#endif






/* Chloride properties
 * ------------------- */
#define M_Cl        MolarMassOfMolecule(Cl)
#define M_OH        MolarMassOfMolecule(OH)
/* Chloride adsorption curve */
#define AlphaCoef_CURVE     adsorbedchloridecurve_a
#define BetaCoef_CURVE      adsorbedchloridecurve_b
#define AlphaCoef(x) \
        (Curve_ComputeValue(AlphaCoef_CURVE,x))
#define BetaCoef(x) \
        (Curve_ComputeValue(BetaCoef_CURVE,x))
#ifdef E_Chlorine
#define AdsorbedChloridePerUnitMoleOfCSH(c_cl,x) \
        (AlphaCoef(x) * (c_cl) / (1. + BetaCoef(x) * (c_cl)))
#else
#define AdsorbedChloridePerUnitMoleOfCSH(c_cl,x) \
        (0)
#endif



/* Sodium adsorption curve 
 * ------------------------- */
#define RNa(x) \
        (Curve_ComputeValue(Element_GetCurve(el) + 4,x))
#define AdsorbedSodiumPerUnitMoleOfCSH(c_na,x) \
        ((c_na < 0.3) ? (RNa(x) * (c_na)) : (RNa(x) * 0.3))



/* Potassium adsorption curve 
 * ------------------------- */
#define RK(x) \
        (Curve_ComputeValue(Element_GetCurve(el) + 5,x))
#define AdsorbedPotassiumPerUnitMoleOfCSH(c_k,x) \
        ((c_k < 0.3) ? (RK(x) * (c_k)) : (RK(x) * 0.3)) 



/* Element contents in solid phases  */
#define n_ca_ref                           (n_ch0)
#define n_si_ref                           (n_csh0)



/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(const char* s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static int     ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;

static double  dn1_caoh2sdt(double,double) ;
static double  CHSolidContent_kin1(double,double,double) ;

static int     ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,double**,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;

static int     TangentCoefficients(Element_t*,double,double*) ;


static void    ComputePhysicoChemicalProperties(double) ;

static int     concentrations_oh_na_k(double,double,double,double,double,double) ;

static double  PermeabilityCoefficient_KozenyCarman(Element_t*,double) ;
static double  PermeabilityCoefficient_VermaPruess(Element_t*,double) ;
static double  TortuosityToLiquid_OhJang(double,double) ;
static double  TortuosityToLiquid_BazantNajjar(double,double) ;
static double  TortuosityToLiquid_Xie(double,double) ;
static double  TortuosityToGas(double,double) ;

static double  saturationdegree(double,double,Curve_t*) ;


/* Internal parameters */
static double phi0 ;
static double phi_min ;
static double kl_int ;
static double kg_int ;
static double frac,phi_r ;
static Curve_t* saturationcurve ;
static Curve_t* relativepermliqcurve ;
static Curve_t* relativepermgascurve ;
static Curve_t* molarvolumeofcshcurve ;
#ifdef E_Chlorine
static Curve_t* adsorbedchloridecurve_a ;
static Curve_t* adsorbedchloridecurve_b ;
#endif
static double a_2,c_2 ;
static double rate_cc ;
static double n_ch0,n_csh0,c_na0,c_k0 ;

static double p_g0 ;
static double p_l0 ;
static double p_v0 ;

static double d_co2 ;
static double d_vap ;

static double mu_l ;
static double mu_g ;

static double p_c3 ;

static double RT ;

static double K_w ;

static double rho_l0 ;

static CementSolutionDiffusion_t* csd = NULL ;
static HardenedCementChemistry_t* hcc = NULL ;



#include "PhysicalConstant.h"
#include "AtmosphericPressure.h"
#include "WaterViscosity.h"
#include "AirViscosity.h"
#include "WaterVaporPressure.h"
#include "DiffusionCoefficientOfMoleculeInAir.h"
#include "EquilibriumConstantOfHomogeneousReactionInWater.h"


void ComputePhysicoChemicalProperties(double TK)
{

  /* Diffusion Coefficient Of Molecules In Air (dm2/s) */
  d_co2   = DiffusionCoefficientOfMoleculeInAir(CO2,TK) ;
  d_vap   = DiffusionCoefficientOfMoleculeInAir(H2O,TK) ;
  
  /* Viscosities */
  mu_l    = WaterViscosity(TK) ;
  mu_g    = AirViscosity(TK) ;
  
  /* Water vapor pressure */
  p_v0    = WaterVaporPressure(TK) ;
  
  /* Reference pressures */
  p_l0    = 0 ; //AtmosphericPressure ;
  p_g0    = 0 ; //AtmosphericPressure ;
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
  
  /* Chemical constants */
  K_w = EquilibriumConstantOfHomogeneousReactionInWater(H2O__H_OH,TK) ;
  
  /* Liquid mass density */
  rho_l0 = 1 * kg/dm3 ;
}



enum {
I_P_L  = NEQ   ,

I_N_C          ,
I_N_Ca         ,
I_N_Si         ,
I_N_K          ,
I_N_Na         ,
I_Mass         ,
I_N_Q          ,
I_N_Cl         ,

I_N_Si_S       ,
I_N_Ca_S       ,

I_N_CH         ,
I_N_CC         ,
I_N_CSH        ,

I_V_S          ,
I_V_S0         ,

I_Phi          ,

I_V_CSH        ,

I_C_OH         ,

I_RHO_H2O_g    ,

I_M_Air        ,

I_P_G          ,
I_C_CO2        ,

I_S_L          ,
I_Last
} ;


#define NbOfVariables    (I_Last)
static double Variables[Element_MaxNbOfNodes][2*NbOfVariables] ;
//static double Variables_n[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;
#define Variables_n(x)   ((x) + NbOfVariables)


enum {
I_W_C           ,
I_W_Ca          ,
I_W_Si          ,
I_W_Na          ,
I_W_K           ,
I_W_tot         ,
I_W_q           ,
I_W_Cl          ,
I_W_Air         ,
I_W_Last
} ;


#define NbOfVariableFluxes    (I_W_Last)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;



int pm(const char* s)
{
       if(strcmp(s,"porosity") == 0)      return (0) ;
  else if(strcmp(s,"k_int") == 0)         return (1) ;
  else if(strcmp(s,"kl_int") == 0)        return (1) ;
  else if(strcmp(s,"N_CH") == 0)          return (2) ;
  else if(strcmp(s,"N_CSH") == 0)         return (4) ;
  else if(strcmp(s,"C_K") == 0)           return (5) ;
  else if(strcmp(s,"C_Na") == 0)          return (6) ;
  else if(strcmp(s,"A_2") == 0)           return (8) ;
  else if(strcmp(s,"C_2") == 0)           return (9) ;
  else if(strcmp(s,"Radius_CH") == 0)     return (10) ;
  else if(strcmp(s,"D") == 0)             return (11) ;
  else if(strcmp(s,"Tau") == 0)           return (12) ;
  else if(strcmp(s,"frac") == 0)          return (13) ;
  else if(strcmp(s,"phi_r") == 0)         return (14) ;
  else if(strcmp(s,"porosity_min") == 0)  return (15) ;
  else if(strcmp(s,"Rate_Calcite") == 0)  return (16) ;
  else if(strcmp(s,"kg_int") == 0)        return (17) ;
  else if(strcmp(s,"p_c3") == 0)          return (18) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  phi0     = GetProperty("porosity") ;
  kl_int   = GetProperty("kl_int") ;
  kg_int   = GetProperty("kg_int") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CH") ;
  n_csh0   = GetProperty("N_CSH") ;
  c_na0    = GetProperty("C_Na") ;
  c_k0     = GetProperty("C_K") ;
  frac     = GetProperty("frac") ;
  phi_r    = GetProperty("phi_r") ;
  phi_min  = GetProperty("porosity_min") ;
  rate_cc  = GetProperty("Rate_Calcite") ;
  p_c3     = GetProperty("p_c3") ;
  
  saturationcurve         = Element_FindCurve(el,"s_l") ;
  relativepermliqcurve    = Element_FindCurve(el,"kl_r") ;
#ifdef E_Air
  relativepermgascurve    = Element_FindCurve(el,"kg_r") ;
#endif
  molarvolumeofcshcurve   = Element_FindCurve(el,"v_csh") ;
#ifdef E_Chlorine
  adsorbedchloridecurve_a = Element_FindCurve(el,"alpha") ;
  adsorbedchloridecurve_b = Element_FindCurve(el,"beta") ;
#endif
}


int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Carbon   ,"carbon") ;
  Model_CopyNameOfEquation(model,E_charge   ,"charge") ;
  Model_CopyNameOfEquation(model,E_mass     ,"mass") ;
  Model_CopyNameOfEquation(model,E_Calcium  ,"calcium") ;
  Model_CopyNameOfEquation(model,E_Sodium   ,"sodium") ;
  Model_CopyNameOfEquation(model,E_Potassium,"potassium") ;
  Model_CopyNameOfEquation(model,E_Silicon  ,"silicon") ;
#ifdef E_eneutral
  Model_CopyNameOfEquation(model,E_eneutral ,"electroneutrality") ;
#endif
#ifdef E_Chlorine
  Model_CopyNameOfEquation(model,E_Chlorine ,"chlorine") ;
#endif
#ifdef E_Air
  Model_CopyNameOfEquation(model,E_Air      ,"air") ;
#endif
  
  
#ifdef U_LogC_CO2
  Model_CopyNameOfUnknown(model,E_Carbon ,"logc_co2") ;
#else
  Model_CopyNameOfUnknown(model,E_Carbon ,"c_co2") ;
#endif

  Model_CopyNameOfUnknown(model,E_Silicon,"z_si") ;
  Model_CopyNameOfUnknown(model,E_mass    ,"p_l") ;

#if defined (U_ZN_Ca_S)
  Model_CopyNameOfUnknown(model,E_Calcium,"z_ca") ;
#elif defined (U_LogS_CH)
  Model_CopyNameOfUnknown(model,E_Calcium,"logs_ch") ;
#endif

  Model_CopyNameOfUnknown(model,E_charge    ,"psi") ;
  
#ifdef U_LogC_Na
  Model_CopyNameOfUnknown(model,E_Sodium   ,"logc_na") ;
#else
  Model_CopyNameOfUnknown(model,E_Sodium   ,"c_na") ;
#endif

#ifdef U_LogC_K
  Model_CopyNameOfUnknown(model,E_Potassium    ,"logc_k") ;
#else
  Model_CopyNameOfUnknown(model,E_Potassium    ,"c_k") ;
#endif

#ifdef E_eneutral
  #if defined (U_LogC_OH)
    Model_CopyNameOfUnknown(model,E_eneutral, "logc_oh") ;
  #elif defined (U_C_OH)
    Model_CopyNameOfUnknown(model,E_eneutral, "c_oh") ;
  #else
    Model_CopyNameOfUnknown(model,E_eneutral, "z_oh") ;
  #endif
#endif

#ifdef E_Chlorine
  #ifdef U_LogC_Cl
    Model_CopyNameOfUnknown(model,E_Chlorine, "logc_cl") ;
  #else
    Model_CopyNameOfUnknown(model,E_Chlorine, "c_cl") ;
  #endif
#endif

#ifdef E_Air
  Model_CopyNameOfUnknown(model,E_Air, "p_g") ;
#endif
  
  //Model_GetNbOfVariables(model) = NbOfVariables ;
  //Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 19 ;
  
  InternationalSystemOfUnits_UseAsLength("decimeter") ;
  InternationalSystemOfUnits_UseAsMass("hectogram") ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Default initialization */
  {
    double h   = 5.6e-6 * (mol/dm2/sec) ;  /* (mol/dm2/s) these MT p 223 */
    double R_0 = Material_GetProperty(mat)[pm("Radius_CH")] ; /* (dm) */
    double D   = Material_GetProperty(mat)[pm("D")] ; /* (mol/dm/s) */
    
    if(R_0 == 0.) R_0 = 40.e-5 * dm ;
    
    if(D == 0.) D = 7e-15 * (mol/dm/sec) ;
    
    /* Initial (or reference) CH molar content */
    n_ch0 = Material_GetProperty(mat)[pm("N_CH")] ;
    
    {
      double t_ch = Material_GetProperty(mat)[pm("Tau")] ; /* (s) */
      
      if(t_ch == 0) {
        t_ch = R_0/(3*h*V_CH) ;     /* (s) approx 721.5 s */
        /* t_ch = R_0*R_0/(3*V_CH*D) ; */ /* (s) approx 2.3e8 s */
        Material_GetProperty(mat)[pm("Tau")] = t_ch ;
      }
      
      a_2 = n_ch0/t_ch ;  /* (mol/dm3/s) M. Thiery, PhD thesis, p 227 */
    }
    
    c_2 = h*R_0/D ;     /* (no dim) M. Thiery, PhD thesis p 228 */
  
    Material_GetProperty(mat)[pm("A_2")] = a_2 ;
    Material_GetProperty(mat)[pm("C_2")] = c_2 ;
  }
  
  {
    /* Initial (or reference) CSH molar content */
    n_csh0 = Material_GetProperty(mat)[pm("N_CSH")] ;
    if(n_csh0 == 0) n_csh0 = 1. ;
    
    Material_GetProperty(mat)[pm("N_CSH")] = n_csh0 ;
  }
  
  {
    frac = Material_GetProperty(mat)[pm("frac")] ;
    if(frac == 0) frac = 0.8 ;
    
    Material_GetProperty(mat)[pm("frac")] = frac ;
  }
  
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;

  {
    if(!csd) csd = CementSolutionDiffusion_Create() ;
    if(!hcc) hcc = HardenedCementChemistry_Create() ;
    
    HardenedCementChemistry_SetRoomTemperature(hcc,TEMPERATURE) ;
    
    CementSolutionDiffusion_SetRoomTemperature(csd,TEMPERATURE) ;
  }
  
  
  {
      Curves_t* curves = Material_GetCurves(mat) ;
      int i ;

      if((i = Curves_FindCurveIndex(curves,"s_l")) < 0) {
        arret("ReadMatProp: no s_l - p_c curve") ;
      }

      if((i = Curves_FindCurveIndex(curves,"kl_r")) < 0) {
        arret("ReadMatProp: no kl_r - p_c curve") ;
      }

#ifdef E_Air
      if((i = Curves_FindCurveIndex(curves,"kg_r")) < 0) {
        arret("ReadMatProp: no kg_r - p_c curve") ;
      }
#endif

      if((i = Curves_FindCurveIndex(curves,"v_csh")) < 0) {
        arret("ReadMatProp: no v_csh - x_csh curve") ;
      }

      if((i = Curves_FindCurveIndex(curves,"X_CSH")) >= 0) {
        Curve_t* curve = Curves_GetCurve(curves) + i ;
      
        HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(hcc) = curve ;
      }

      if((i = Curves_FindCurveIndex(curves,"Z_CSH")) >= 0) {
        Curve_t* curve = Curves_GetCurve(curves) + i ;
      
        HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(hcc) = curve ;
      }

      if((i = Curves_FindCurveIndex(curves,"S_SH")) >= 0) {
        Curve_t* curve = Curves_GetCurve(curves) + i ;
      
        HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(hcc) = curve ;
      }
  }
  
  return(NbOfProp) ;
}


int PrintModelChar(Model_t* model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The set of 7 equations is:\n") ;
  printf("\t- Mass balance of C      (carbon)\n") ;
  printf("\t- Mass balance of Ca     (calcium)\n") ;
  printf("\t- Mass balance of Si     (silicon)\n") ;
  printf("\t- Mass balance of Na     (sodium)\n") ;
  printf("\t- Mass balance of K      (potassium)\n") ;
#ifdef E_Chlorine
  printf("\t- Mass balance of Cl     (chlorine)\n") ;
#endif
  printf("\t- Total mass balance     (mass)\n") ;
  printf("\t- Charge balance         (charge)\n") ;
#ifdef E_eneutral
  printf("\t- Electroneutrality      (electroneutrality)\n") ;
#endif
#ifdef E_Air
  printf("\t- Mass blance of air     (air)\n") ;
#endif
  
  printf("\n") ;
  printf("The 7-10 primary unknowns are:\n") ;
  printf("\t- Liquid pressure                  (p_l)\n") ;
#ifdef E_Air
  printf("\t- Gas pressure                     (p_g)\n") ;
#endif
  printf("\t- Electric potential               (psi) \n") ;
  printf("\t- Carbon dioxide gas concentration (c_co2 or logc_co2)\n") ;
  printf("\t- Potassium concentration          (c_k or logc_k)\n") ;
  printf("\t- Sodium concentration             (c_na or logc_na)\n") ;
#if defined (U_ZN_Ca_S)
  printf("\t- Zeta unknown for calcium         (z_ca)\n") ;
  printf("\t   \t z_ca is defined as:\n") ;
  printf("\t   \t z_ca = n_ch/n0 + log(s_ch)  for c_co2 < c_co2_eq\n") ;
  printf("\t   \t z_ca = n_cc/n0 + log(s_cc)  for c_co2 > c_co2_eq\n") ;
#elif defined (U_LogS_CH)
  printf("\t- Log10 of saturation index of CH  (logs_ch)\n") ;
#endif
  printf("\t- Zeta unknown for silicon         (z_si)\n") ;
  printf("\t   \t z_si is defined as:\n") ;
  printf("\t   \t z_si = n_si/n0 + log(s_sh/s_sh_eq)\n") ;
#ifdef E_Chlorine
  printf("\t- Chloride ion concentration       (c_cl or logc_cl)\n") ;
#endif
#ifdef E_eneutral
  printf("\t- Hydroxide ion concentration     (c_oh or logc_oh or z_oh)\n") ;
#endif
  
  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length    : dm !\n") ;
  printf("\t time      : s !\n") ;
  printf("\t mass      : hg !\n") ;
  printf("\t pressure  : Pa !\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;


  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;
  fprintf(ficd,"kl_int = 1.4e-17   # Intrinsic permeability (dm2)\n") ;
  fprintf(ficd,"N_CH = 3.9        # Initial content in Ca(OH)2 (mol/L)\n") ;
  fprintf(ficd,"Radius_CH = 40.e-5  # Portlandite crystal radius \n") ;
  fprintf(ficd,"N_CSH = 2.4        # Initial content in CSH (mol/L)\n") ;
  fprintf(ficd,"C_Na = 0.019      # Total content in Na (mol/L)\n") ;
  fprintf(ficd,"C_K  = 0.012      # Total content in K  (mol/L)\n") ;
  fprintf(ficd,"D = 7.e-15        # Diffusion coef in CC (dm/mol/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Kinetic coef 2 (dm/mol/s)\n") ;
  fprintf(ficd,"frac = 0.8        # Fractionnal length of pore bodies\n") ;
  fprintf(ficd,"phi_r = 0.7       # Porosity for which permeability vanishes\n") ;
  fprintf(ficd,"Curves = my_file  # File name: p_c S_l kl_r\n") ;  
  fprintf(ficd,"Curves = my_file  # File name: si_ch C/S H/S V_csh\n") ;  
  fprintf(ficd,"Curves = my_file  # File name: s_l kl_r kg_r\n") ;  

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NbOfImplicitTerms ;
  Element_GetNbOfExplicitTerms(el) = (Element_IsSubmanifold(el)) ? 0 : NbOfExplicitTerms ;
  Element_GetNbOfConstantTerms(el) = NbOfConstantTerms ;
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
/* Initialise les variables du systeme (f,va) */ 
{
  double* f  = Element_GetImplicitTerm(el) ;
  double* v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  {
    double c_na_tot = c_na0 ;
    double c_k_tot  = c_k0 ;
    int i ;

    for(i = 0 ; i < nn ; i++) {
      double c_na       = C_Na(i) ;
      double c_k        = C_K(i) ;
      double c_co2      = C_CO2(i) ;
      double u_calcium  = U_Calcium(i) ;
      double u_silicon  = U_Silicon(i) ;
      #ifdef E_Chlorine
      double c_cl       = C_Cl(i) ;
      #else
      double c_cl       = 1.e-99 ;
      #endif
      
      if(c_na_tot > 0 && c_k_tot > 0) {
        c_na   = c_na_tot ;
        c_k    = c_k_tot ;

        /* Compute the concentrations of alkalis Na and K */
        concentrations_oh_na_k(c_co2,u_calcium,u_silicon,c_cl,c_na_tot,c_k_tot) ;
  
        c_na = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Na) ;
        c_k  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,K) ;


#ifdef U_LogC_Na
        LogC_Na(i)  = log10(c_na) ;
#else
        C_Na(i)     = c_na ;
#endif
#ifdef U_LogC_K
        LogC_K(i)   = log10(c_k) ;
#else
        C_K(i)      = c_k ;
#endif
    
      /* Solve cement chemistry */
      } else {
        double c_co2aq    = k_h*c_co2 ;
        double logc_co2aq = log10(c_co2aq) ;
        double logc_na    = log10(c_na) ;
        double logc_k     = log10(c_k) ;
        double logc_cl    = log10(c_cl) ;
        double logc_oh    = -7 ;
        double psi        = 0 ;

        #if defined (U_ZN_Ca_S)
        {
          double si_ch_cc = Log10SaturationIndexOfCcH(u_calcium) ;
          
          HardenedCementChemistry_SetInput(hcc,SI_CH_CC,si_ch_cc) ;
        }
        #elif defined (U_LogS_CH)
        {
          double si_ch = Log10SaturationIndexOfCH(u_calcium) ;
          
          HardenedCementChemistry_SetInput(hcc,SI_CH,si_ch) ;
        }
        #endif
        
        {
          double si_csh = Log10SaturationIndexOfCSH(u_silicon) ;
          
          HardenedCementChemistry_SetInput(hcc,SI_CSH,si_csh) ;
        }
        
        HardenedCementChemistry_SetInput(hcc,LogC_CO2,logc_co2aq) ;
        HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
        HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
        HardenedCementChemistry_SetInput(hcc,LogC_OH,logc_oh) ;
        HardenedCementChemistry_GetElectricPotential(hcc) = psi ;
        HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = c_cl ;
        HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,Cl) = logc_cl ;
  
        HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_CO2_H2O) ;

       {
         int k = HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      
         if(k < 0) return(1) ;
       }
      }
      
      /* pH */
      {
        double c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
        double c_h  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,H) ;

        #ifdef E_eneutral
          #if defined (U_LogC_OH)
            LogC_OH(i) = log10(c_oh) ;
          #elif defined (U_C_OH)
            C_OH(i)    = c_oh ;
          #elif defined (U_Z_OH)
            Z_OH(i)    = Z_OHDefinition(c_oh,c_h) ;
          #endif
        #else
          C_OH(i)    = c_oh ;
        #endif
      }
      
      /* Solid contents */
      {
        double s_ch       = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
        double s_cc       = HardenedCementChemistry_GetSaturationIndexOf(hcc,CC) ;
        double n_ch       = InitialCHSolidContent(u_calcium,s_ch,s_cc) ;
        double n_cc       = InitialCCSolidContent(u_calcium,s_ch,s_cc) ;
        double n_csh      = CSHSolidContent(u_silicon) ;
        double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
        double v_csh      = MolarVolumeOfCSH(x_csh) ;
        double v_s0       = V_CH*n_ch + V_CC*n_cc + v_csh*n_csh ;
        
        V_S0(i)    = v_s0 ;
        N_CH(i)    = n_ch ;
        N_CC(i)    = n_cc ;
      }
    }
  }
  

  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x   = ComputeVariables(el,u,u,f,0,0,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
      
      if(!x) return(1) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
    
      /* Back up */
      N_C(i)   = x[I_N_C] ;
      N_Ca(i)  = x[I_N_Ca] ;
      N_Na(i)  = x[I_N_Na] ;
      N_Si(i)  = x[I_N_Si] ;
      N_K(i)   = x[I_N_K] ; 
      M_tot(i) = x[I_Mass] ;
      N_q(i)   = x[I_N_Q] ;
      N_Cl(i)  = x[I_N_Cl] ;
      M_Air(i) = x[I_M_Air] ;

      /* Solid contents */
      N_CH(i) = x[I_N_CH] ;
      N_CC(i) = x[I_N_CC] ;
      
      /* pH */
      #ifndef E_eneutral
        C_OH(i)    = x[I_C_OH] ;
      #endif
    }
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  {
    int i = ComputeTransferCoefficients(el,u,f) ;
    
    if(i) return(1) ;
  }


  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,u,i,j) ;
        
        W_C(i,j)     = w[I_W_C] ;
        W_Ca(i,j)    = w[I_W_Ca] ;
        W_Na(i,j)    = w[I_W_Na] ;
        W_Si(i,j)    = w[I_W_Si] ;
        W_q(i,j)     = w[I_W_q] ;
        W_K(i,j)     = w[I_W_K] ;
        W_tot(i,j)   = w[I_W_tot] ;
        W_Cl(i,j)    = w[I_W_Cl] ;
        W_Air(i,j)   = w[I_W_Air] ;
        
        W_C(j,i)     = - w[I_W_C] ;
        W_Ca(j,i)    = - w[I_W_Ca] ;
        W_Na(j,i)    = - w[I_W_Na] ;
        W_Si(j,i)    = - w[I_W_Si] ;
        W_q(j,i)     = - w[I_W_q] ;
        W_K(j,i)     = - w[I_W_K] ;
        W_tot(j,i)   = - w[I_W_tot] ;
        W_Cl(j,i)    = - w[I_W_Cl] ;
        W_Air(j,i)   = - w[I_W_Air] ;
      }
    }
  }
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
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
  {
    int i = ComputeTransferCoefficients(el,u,f) ;
    
    if(i) return(1) ;
  }

  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Molar contents */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x   = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
      
      if(!x) return(1) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
    
      /* Back up */
      N_C(i)   = x[I_N_C] ;
      N_Ca(i)  = x[I_N_Ca] ;
      N_Na(i)  = x[I_N_Na] ;
      N_Si(i)  = x[I_N_Si] ;
      N_K(i)   = x[I_N_K] ; 
      M_tot(i) = x[I_Mass] ;
      N_q(i)   = x[I_N_Q] ;
      N_Cl(i)  = x[I_N_Cl] ;
      M_Air(i) = x[I_M_Air] ;

      /* Solid contents */
      N_CH(i)  = x[I_N_CH] ;
      N_CC(i)  = x[I_N_CC] ;
      
      /* pH */
#ifndef E_eneutral
      C_OH(i)  = x[I_C_OH] ;
#endif

      {
        double c_co2      = x[I_C_CO2] ;
        
        double c_h2o = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,H2O) ;
        double c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
        double c_ca= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Ca) ;
        double c_na= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Na) ;
        double c_k = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,K) ;
        
        double n_csh     = x[I_N_CSH] ;
        double n_ch      = x[I_N_CH] ;
      
        if(c_co2 < 0 || c_oh <= 0 || c_h2o <= 0 || c_na < 0 || c_k < 0 || c_ca < 0 || n_csh < 0. || n_ch < 0.) {
          double x0 = Element_GetNodeCoordinate(el,i)[0] ;
          double n_cc     = x[I_N_CC] ;
          double c_naoh   = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaOH) ;
          double c_nahco3 = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaHCO3) ;
          double c_naco3  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaCO3) ;
          double c_cl     = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) ;
          printf("\n") ;
          printf("en x     = %e\n",x0) ;
          printf("c_co2    = %e\n",c_co2) ;
          printf("c_oh     = %e\n",c_oh) ;
          printf("c_h2o    = %e\n",c_h2o) ;
          printf("n_cc     = %e\n",n_cc) ;
          printf("c_na     = %e\n",c_na) ;
          printf("c_k      = %e\n",c_k) ;
          printf("c_ca     = %e\n",c_ca) ;
          printf("n_csh    = %e\n",n_csh) ;
          printf("c_naoh   = %e\n",c_naoh) ;
          printf("c_nahco3 = %e\n",c_nahco3) ;
          printf("c_naco3  = %e\n",c_naco3) ;
          printf("c_cl     = %e\n",c_cl) ;
          return(1) ;
        }
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;
  

  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,u,i,j) ;
        
        W_C(i,j)     = w[I_W_C] ;
        W_Ca(i,j)    = w[I_W_Ca] ;
        W_Na(i,j)    = w[I_W_Na] ;
        W_Si(i,j)    = w[I_W_Si] ;
        W_q(i,j)     = w[I_W_q] ;
        W_K(i,j)     = w[I_W_K] ;
        W_tot(i,j)   = w[I_W_tot] ;
        W_Cl(i,j)    = w[I_W_Cl] ;
        W_Air(i,j)   = w[I_W_Air] ;
        
        W_C(j,i)     = - w[I_W_C] ;
        W_Ca(j,i)    = - w[I_W_Ca] ;
        W_Na(j,i)    = - w[I_W_Na] ;
        W_Si(j,i)    = - w[I_W_Si] ;
        W_q(j,i)     = - w[I_W_q] ;
        W_K(j,i)     = - w[I_W_K] ;
        W_tot(j,i)   = - w[I_W_tot] ;
        W_Cl(j,i)    = - w[I_W_Cl] ;
        W_Air(j,i)   = - w[I_W_Air] ;
      }
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
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
  
  if(TangentCoefficients(el,dt,c) < 0) return(1) ;
  
  {
    double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }


/* On output TangentCoefficients has computed the derivatives wrt
 * LogC_CO2, LogC_Na, LogC_K, LogC_OH, LogC_Cl
 * (see ComputeVariables and ComputeVariablesDerivatives). */
 
#ifdef U_C_CO2
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_Carbon)     /= Ln10*C_CO2(0) ;
      K(i,E_Carbon+NEQ) /= Ln10*C_CO2(1) ;
    }
  }
#endif

#ifdef U_C_Na
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_Sodium)     /= Ln10*C_Na(0) ;
      K(i,E_Sodium+NEQ) /= Ln10*C_Na(1) ;
    }
  }
#endif

#ifdef U_C_K
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_Potassium)     /= Ln10*C_K(0) ;
      K(i,E_Potassium+NEQ) /= Ln10*C_K(1) ;
    }
  }
#endif
  
#ifdef E_eneutral
  #if defined (U_C_OH)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_eneutral)     /= Ln10*C_OH(0) ;
      K(i,E_eneutral+NEQ) /= Ln10*C_OH(1) ;
    }
  }
  #elif defined (U_Z_OH)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_eneutral)     *= dLogC_OHdZ_OH(Z_OH(0)) ;
      K(i,E_eneutral+NEQ) *= dLogC_OHdZ_OH(Z_OH(1)) ;
    }
  }
  #endif
#endif
  
#ifdef E_Chlorine
  #ifdef U_C_Cl
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_Chlorine)     /= Ln10*C_Cl(0) ;
      K(i,E_Chlorine+NEQ) /= Ln10*C_Cl(1) ;
    }
  }
  #endif
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  
  /*
    Initialization
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Conservation of element C: (N_C - N_Cn) + dt * div(W_C) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_C,NW_Cn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Carbon) -= r1[i] ;
    }
  }
  
  /*
    Conservation of charge: div(W_q) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = 0 ;
        } else {
          //g[i*nn + j] = dt * W_q(i,j) ;
          g[i*nn + j] = W_q(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_charge) -= r1[i] ;
      }
    }
  }
  
  /*
    Conservation of total mass: (M_tot - M_totn) + dt * div(W_tot) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,MW_tot,MW_totn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_mass) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element Ca: (N_Ca - N_Can) + dt * div(W_Ca) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Ca,NW_Can,dt) ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Calcium) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element Na: (N_Na - N_Nan) + dt * div(W_Na) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Na,NW_Nan,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Sodium) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element K: (N_K - N_Kn) + dt * div(W_K) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_K,NW_Kn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Potassium) -= r1[i] ;
    }
  }

  /*
    Conservation of element Si: (N_Si - N_Sin) + dt * div(W_Si) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Si,NW_Sin,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Silicon) -= r1[i] ;
    }
  }

#ifdef E_Chlorine
  /*
    Conservation of element Cl: (N_Cl - N_Cln) + dt * div(W_Cl) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Cl,NW_Cln,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Chlorine) -= r1[i] ;
    }
  }
#endif
  
  
#ifdef E_eneutral
  /*
    Electroneutrality
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = N_q(i) ;
        } else {
          g[i*nn + j] = 0 ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_eneutral) -= r1[i] ;
      }
    }
  }
#endif


#ifdef E_Air
  /*
    Conservation of dry air mass: (M_Air - M_Airn) + dt * div(W_Air) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,MW_Air,MW_Airn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Air) -= r1[i] ;
    }
  }
#endif

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nso = 68 ;
  int    i ;

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
    double* x = ComputeVariables(el,u,u,f,t,0,j) ;
      
    if(!x) return(0) ;
    
    /* Macros */
#define ptC(CPD)   &(HardenedCementChemistry_GetAqueousConcentrationOf(hcc,CPD))
#define ptEC(CPD)  &(HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,CPD))
#define ptS(CPD)   &(HardenedCementChemistry_GetSaturationIndexOf(hcc,CPD))
#define ptPSI      &(HardenedCementChemistry_GetElectricPotential(hcc))
#define ptX_CSH    &(HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc))
#define ptZ_CSH    &(HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc))



    /* Outputs */
    i = 0 ;
    
    /* Liquid pressure */
    Result_Store(r + i++,x + I_P_L,"p_l",1) ;
    
    /* Liquid saturation degree */
    Result_Store(r + i++,x + I_S_L,"saturation",1) ;
    
    Result_Store(r + i++,x + I_Phi,"porosity",1) ;
    
    /* Concentration in gas phase */
    Result_Store(r + i++,x + I_C_CO2,"c_co2",1) ;
    
    /* Element concentrations in liquid phase */
    Result_Store(r + i++,ptEC(Ca ),"c_ca_l",1) ;
    Result_Store(r + i++,ptEC(Si ),"c_si_l",1) ;
    Result_Store(r + i++,ptEC(Na ),"c_na_l",1) ;
    Result_Store(r + i++,ptEC(K  ),"c_k_l" ,1) ;
    Result_Store(r + i++,ptEC(C  ),"c_c_l" ,1) ;
    Result_Store(r + i++,ptEC(Cl ),"c_cl_l",1) ;
    
    /* Portlandite */
    Result_Store(r + i++,x + I_N_CH,"n_CH",1) ;
    Result_Store(r + i++,ptS(CH),"s_ch",1) ;
    
    /* C-S-H */
    Result_Store(r + i++,x + I_N_CSH,"n_CSH",1) ;
    Result_Store(r + i++,ptX_CSH,"x_csh",1) ;
    Result_Store(r + i++,ptS(SH),"s_sh",1) ;
    
    /* Calcite */
    Result_Store(r + i++,x + I_N_CC,"n_CC",1) ;
    Result_Store(r + i++,ptS(CC),"s_cc",1) ;
    
    
    /* Ion concentrations in liquid phase */
    Result_Store(r + i++,ptC(H ),"c_h",1) ;
    Result_Store(r + i++,ptC(OH),"c_oh",1) ;
    {
      double c_h       = *(ptC(H )) ;
      double ph        = - log10(c_h) ;
      
      Result_Store(r + i++,&ph,"ph",1) ;
    }
    
    Result_Store(r + i++,ptC(Ca  ),"c_ca",1) ;
    Result_Store(r + i++,ptC(CaOH),"c_caoh",1) ;
    
    Result_Store(r + i++,ptC(H2SiO4),"c_h2sio4",1) ;
    Result_Store(r + i++,ptC(H3SiO4),"c_h3sio4",1) ;
    Result_Store(r + i++,ptC(H4SiO4),"c_h4sio4",1) ;
    
    Result_Store(r + i++,ptC(Na  ),"c_na",1) ;
    Result_Store(r + i++,ptC(NaOH),"c_naoh",1) ;
    
    Result_Store(r + i++,ptC(K  ),"c_k",1) ;
    Result_Store(r + i++,ptC(KOH),"c_koh",1) ;
    
    Result_Store(r + i++,ptC(CO3 ),"c_co3",1) ;
    Result_Store(r + i++,ptC(HCO3),"c_hco3",1) ;
    
    Result_Store(r + i++,ptC(CaH2SiO4),"c_cah2sio4",1) ;
    Result_Store(r + i++,ptC(CaH3SiO4),"c_cah3sio4",1) ;
    
    Result_Store(r + i++,ptC(CaHCO3),"c_cahco3",1) ;
    Result_Store(r + i++,ptC(CaCO3),"c_caco3aq",1) ;
    Result_Store(r + i++,ptC(CaO2H2),"c_caoh2aq",1) ;
    
    Result_Store(r + i++,ptC(NaHCO3),"c_nahco3",1) ;
    Result_Store(r + i++,ptC(NaCO3),"c_naco3",1) ;
    
    Result_Store(r + i++,ptC(Cl),"c_cl",1) ;
    
    /* Total element contents */
    Result_Store(r + i++,x + I_N_Ca,"n_Ca",1) ;
    Result_Store(r + i++,x + I_N_Si,"n_Si",1) ;
    Result_Store(r + i++,x + I_N_Na,"n_Na",1) ;
    Result_Store(r + i++,x + I_N_K ,"n_K" ,1) ;
    Result_Store(r + i++,x + I_N_C ,"n_C" ,1) ;
    Result_Store(r + i++,x + I_N_Cl,"n_Cl" ,1) ;
    
    /* Total mass content */
    Result_Store(r + i++,x + I_Mass,"total mass",1) ;
    
    /* Mass flows */
    Result_Store(r + i++,&(W_tot(0,1)),"total mass flow",1) ;
    Result_Store(r + i++,&(W_C(0,1)),"carbon mass flow",1) ;
    Result_Store(r + i++,&(W_Ca(0,1)),"calcium mass flow",1) ;
    Result_Store(r + i++,&(W_Si(0,1)),"silicon mass flow",1) ;
    Result_Store(r + i++,&(W_Na(0,1)),"sodium mass flow",1) ;
    Result_Store(r + i++,&(W_K(0,1)),"potassium mass flow",1) ;
    Result_Store(r + i++,&(W_Cl(0,1)),"chlorine mass flow",1) ;
    
    
    /* Miscellaneous */
    {
      double CS = x[I_N_Ca_S]/x[I_N_Si_S] ;
      
      Result_Store(r + i++,&CS,"Ca/Si ratio",1) ;
    }
    
    {
      double psi = PSI(j) ;
      
      Result_Store(r + i++,&psi,"Electric potential",1) ;
    }
    
    Result_Store(r + i++,x + I_N_Q,"charge",1) ;
    
    {
      double I = HardenedCementChemistry_GetIonicStrength(hcc) ;
      
      Result_Store(r + i++,&I,"I",1) ;
    }
    
    /* Molar volumes */
    {
      double v_solide_csh   = x[I_V_CSH] * x[I_N_CSH] ;
      
      Result_Store(r + i++,&v_solide_csh,"v_csh",1) ;
    }
    {
      double v_solide_ch    = V_CH * x[I_N_CH] ;
      
      Result_Store(r + i++,&v_solide_ch,"v_ch",1) ;
    }
    {
      double v_solide_cc    = V_CC * x[I_N_CC] ;
      
      Result_Store(r + i++,&v_solide_cc,"v_cc",1) ;
    }
    
    /* Gas pressures */
    {
      double p_l        = x[I_P_L] ;
      double p_g        = x[I_P_G] ;
      double p_v        = VaporPressure(p_l) ;
      double h_r        = RelativeHumidity(p_l) ;
      double p_co2      = x[I_C_CO2] * RT ;
      double p_atm      = 101325. ;
      double c_co2      = p_co2 / p_atm * 1.e6 ;
      double p_air      = p_g - p_v - p_co2 ;
      
      Result_Store(r + i++,&p_air,"air pressure",1) ;
      Result_Store(r + i++,&h_r,"humidity",1) ;
      Result_Store(r + i++,&c_co2,"CO2 ppm",1) ;
      Result_Store(r + i++,&p_g,"gas pressure",1) ;
    }
      

    /* Additional outputs */
    {
      double s_l        = x[I_S_L] ;
      double phi        = x[I_Phi] ;
      double coeff_permeability = PermeabilityCoefficient(el,phi) ;
      double k_l  = (kl_int/mu_l)*RelativePermeabilityToLiquid(s_l)*coeff_permeability ;

      Result_Store(r + i++,&k_l,"permeability to liquid",1) ;
      Result_Store(r + i++,&coeff_permeability,"permeability coef",1) ;
    }
    
    /* Adsorbed chloride */
    {
      double n_csh  = x[I_N_CSH] ;
      double c_cl   = (ptC(Cl))[0] ;
      double x_csh  = (ptX_CSH)[0] ;
      double n_cl_s = n_csh * AdsorbedChloridePerUnitMoleOfCSH(c_cl,x_csh) ;
      
      Result_Store(r + i++,&n_cl_s,"adsorbed chloride",1) ;
    }
    
    /* Liquid mass density */
    {
      double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
      
      Result_Store(r + i++,&rho_l,"liquid mass density",1) ;
    }
  }
  
  
  if(i != nso) arret("ComputeOutputs") ;
  return(nso) ;
}


int ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Termes explicites (va)  */
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    i ; 

  /* initialization */
  for(i = 0 ; i < NbOfExplicitTerms ; i++) va[i] = 0. ;
  
  /*
    Transfer coefficients
  */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;
    
    if(!x) return(1) ;


    /* Transport in liquid phase */
    {
      /* saturation */
      double s_l    = x[I_S_L] ;

      /* Porosity */
      double phi    = x[I_Phi] ;

      /* Permeability */
      double coeff_permeability = PermeabilityCoefficient(el,phi) ;
      double k_l  = (kl_int/mu_l)*RelativePermeabilityToLiquid(s_l)*coeff_permeability ;
  
      double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
      double c_c_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,C) ;
      double c_ca_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Ca) ;
      double c_na_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Na) ;
      double c_k_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
      double c_si_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Si) ;
      double c_cl_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Cl) ;
    
      KD_L[i]    = rho_l * k_l ;

      KC_C_L[i]  = c_c_l  / rho_l ;
      KC_Ca_L[i] = c_ca_l / rho_l ;
      KC_Na_L[i] = c_na_l / rho_l ;
      KC_K_L[i]  = c_k_l  / rho_l ;
      KC_Si_L[i] = c_si_l / rho_l ;
      KC_Cl_L[i] = c_cl_l / rho_l ;
    }
    
    
    /* Transport in gas phase (diffusion coef) */
    {
      /* saturation */
      double s_l     = x[I_S_L] ;
      double s_g     = 1 - s_l ;
      /* Porosities */
      double phi     = x[I_Phi] ;
      double phi_g   = phi * s_g ;
      /* tortuosity gas */
      double taugas  = TortuosityToGas(phi,s_l) ;
      
      KF_CO2[i]      = phi_g * taugas * d_co2 ;
      KF_H2O[i]      = phi_g * taugas * d_vap ;
    }


    /* Transport in gas phase (advection coef) */
    {
      /* saturation */
      double s_l       = x[I_S_L] ;
      /* Porosity */
      double phi       = x[I_Phi] ;
      /* Permeability */
      double coeff_permeability = PermeabilityCoefficient(el,phi) ;
      double k_g       = (kg_int/mu_g)*RelativePermeabilityToGas(s_l)*coeff_permeability ;
      /* Gas contents */
      double n_co2_g   = x[I_C_CO2] ;
      /* Pressures */
      double p_l       = x[I_P_L] ;
      double p_g       = x[I_P_G] ;
      double p_co2_g   = n_co2_g * RT ;
      double p_v       = VaporPressure(p_l) ;
      double p_air     = p_g - p_v - p_co2_g ;
      /* Mass densities */
      double rho_co2_g = M_CO2 * n_co2_g ;
      double rho_h2o_g = x[I_RHO_H2O_g] ;
      double rho_air   = MassDensityOfDryAir(p_air) ;
      double rho_g     = rho_air + rho_h2o_g + rho_co2_g ;
      
      KD_G[i]          = rho_g * k_g ;
      
      KC_CO2_G[i]      = n_co2_g / rho_g ;
      KC_H2O_G[i]      = rho_h2o_g / rho_g ;
    }


    /* Liquid tortuosity */
    {
      /* saturation */
      double s_l     = x[I_S_L] ;
      /* Porosity */
      double phi    = x[I_Phi] ;
      /* tortuosity liquid */
      double tauliq =  TortuosityToLiquid(phi,s_l) ;
      
      TORTUOSITY[i] = tauliq ;
    }


    /* Concentrations */
    {
      double* c = HardenedCementChemistry_GetAqueousConcentration(hcc) ;
      int nbc = HardenedCementChemistry_NbOfConcentrations ;
      int j ;
    
      for(j = 0 ; j < nbc ; j++) {
        CONCENTRATION(i)[j] = c[j] ;
      }
    }
  }
  
  return(0) ;
}



double* ComputeVariableFluxes(Element_t* el,double** u,int i,int j)
{
  double* grdij = dVariables ;

  /* Gradients */
  {
    int nn = Element_GetNbOfNodes(el) ;
    FVM_t* fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    double dij  = dist[nn*i + j] ;
    
    {
      double* xi = Variables[i] ;
      double* xj = Variables[j] ;
      int k ;
    
      for(k = 0 ; k < NbOfVariables ; k++)  {
        grdij[k] = (xj[k] - xi[k]) / dij ;
      }
    }

    {
      double* g = CementSolutionDiffusion_GetGradient(csd) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
      double* muj = CementSolutionDiffusion_GetPotentialAtPoint(csd,j) ;
      int n = CementSolutionDiffusion_NbOfConcentrations ;
      int k ;
      
      for(k = 0 ; k < n ; k++) {
        g[k] = (muj[k] - mui[k]) / dij ;
      }
    }
  }

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,grdij,i,j) ;
    
    return(w) ;
  }
}



double* ComputeFluxes(Element_t* el,double* grdij,int i,int j)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes[i] ;
  
  
  /* Transport in liquid phase */
  {
    /* Diffusion in the cement solution */
    {
      /* Gradients */
      double* g = CementSolutionDiffusion_GetGradient(csd) ;
      int n = CementSolutionDiffusion_NbOfConcentrations ;
      double* ci = CONCENTRATION(i) ;
      double* cj = CONCENTRATION(j) ;
      double tortuosity = 0.5 * (TORTUOSITY[i] + TORTUOSITY[j]) ;
      int k ;
      
      for(k = 0 ; k < n ; k++) {
        double rho = 0.5 * (ci[k] + cj[k]) ;
      
        g[k] *= tortuosity * rho ;
      }
      
      /* Molar diffusive fluxes */
      CementSolutionDiffusion_ComputeFluxes(csd) ;
    }
      
    {
      w[I_W_C]   = CementSolutionDiffusion_GetElementFluxOf(csd,C) ;
      w[I_W_Ca]  = CementSolutionDiffusion_GetElementFluxOf(csd,Ca) ;
      w[I_W_Si]  = CementSolutionDiffusion_GetElementFluxOf(csd,Si) ;
      w[I_W_Na]  = CementSolutionDiffusion_GetElementFluxOf(csd,Na) ;
      w[I_W_K ]  = CementSolutionDiffusion_GetElementFluxOf(csd,K) ;
      w[I_W_q ]  = CementSolutionDiffusion_GetIonCurrent(csd) ;
      w[I_W_Cl]  = CementSolutionDiffusion_GetElementFluxOf(csd,Cl) ;
    }
  
    /* Advection in the cement solution */
    {
      /* Mass flux of liquid */
      double grd_p_l = grdij[I_P_L] ;
      double kd_l    = 0.5 * (KD_L[i] + KD_L[j]) ;
      double w_l     = - kd_l * grd_p_l  ;
      
      /* Transfer terms */
      double kc_c_l   = 0.5 * (KC_C_L[i]   + KC_C_L[j]) ;
      double kc_ca_l  = 0.5 * (KC_Ca_L[i]  + KC_Ca_L[j]) ;
      double kc_si_l  = 0.5 * (KC_Si_L[i]  + KC_Si_L[j]) ;
      double kc_na_l  = 0.5 * (KC_Na_L[i]  + KC_Na_L[j]) ;
      double kc_k_l   = 0.5 * (KC_K_L[i]   + KC_K_L[j]) ;
      double kc_cl_l  = 0.5 * (KC_Cl_L[i]  + KC_Cl_L[j]) ;

      /* Mass flux */
      w[I_W_tot]   = w_l  ;
   
      /* Molar fluxes */
      w[I_W_C  ]  += kc_c_l  * w_l  ;
      w[I_W_Ca ]  += kc_ca_l * w_l  ;
      w[I_W_Si ]  += kc_si_l * w_l  ;
      w[I_W_Na ]  += kc_na_l * w_l  ;
      w[I_W_K  ]  += kc_k_l  * w_l  ;
      w[I_W_Cl ]  += kc_cl_l * w_l  ;
    }
  }
  
  /* Transport in gas phase */
  {
    /* Advection and diffusion in gas phase */
    {
      /* Mass flux of gas */
#ifdef E_Air
      double grd_p_g = grdij[I_P_G] ;
      double kd_g    = 0.5 * (KD_G[i] + KD_G[j]) ;
      double w_g     = - kd_g * grd_p_g  ;
#else
      /* Actually kg_int -> infinity and w_g is undetermined ! */
      double w_g     = 0  ;
#endif
      
      
      /* Molar flux of CO2 */
      double grd_co2 = grdij[I_C_CO2] ;
      double kf_co2  = 0.5 * (KF_CO2[i]   + KF_CO2[j]) ;
      double j_co2_g = - kf_co2 * grd_co2  ;
#ifdef E_Air
      double kc_co2  = 0.5 * (KC_CO2_G[i]   + KC_CO2_G[j]) ;
      double w_co2_g =   kc_co2 * w_g + j_co2_g  ;
#else
      double w_co2_g =   j_co2_g  ;
#endif
      
      /* Mass flux of water vapor */
      double grd_h2o = grdij[I_RHO_H2O_g] ;
      double kf_h2o  = 0.5 * (KF_H2O[i]   + KF_H2O[j]) ;
      double j_h2o_g = - kf_h2o * grd_h2o  ;
#ifdef E_Air
      double kc_h2o  = 0.5 * (KC_H2O_G[i]   + KC_H2O_G[j]) ;
      double w_h2o_g =   kc_h2o * w_g + j_h2o_g  ;
#else
      double w_h2o_g =   j_h2o_g  ;
#endif
      
      /* Mass flux of dry air (if E_Air == 0 this is wrong since w_g is undetermined) */
      double w_air   =   w_g - w_h2o_g - M_CO2 * w_co2_g ;
  
      w[I_W_Air]   =  w_air ;
      w[I_W_tot]  +=  w_h2o_g + M_CO2 * w_co2_g ;
      w[I_W_C  ]  +=  w_co2_g ;
    }
  }
    
  return(w) ;
}




int TangentCoefficients(Element_t* el,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    dec = NEQ*NEQ ;
  double dui[NEQ] ;
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < NEQ ; i++) {
    dui[i] =  1.e-2 * ObVal_GetValue(obval + i) ;
  }

  
  dui[E_Carbon   ] =  1.e-4 * ObVal_GetValue(obval + E_Carbon) ;
  dui[E_Sodium   ] =  1.e-3 * ObVal_GetValue(obval + E_Sodium) ;
  dui[E_Potassium] =  1.e-3 * ObVal_GetValue(obval + E_Potassium) ;
  dui[E_Calcium  ] =  1.e-4 * ObVal_GetValue(obval + E_Calcium) ;
  dui[E_Silicon  ] =  1.e-4 * ObVal_GetValue(obval + E_Silicon) ;
  dui[E_mass     ] =  1.e-4 * ObVal_GetValue(obval + E_mass) ;
  dui[E_charge   ] =  1.e+0 * ObVal_GetValue(obval + E_charge) ;
#ifdef E_eneutral
  dui[E_eneutral ] =  1.e-2 * ObVal_GetValue(obval + E_eneutral) ;
#endif
#ifdef E_Chlorine
  dui[E_Chlorine ] =  1.e-3 * ObVal_GetValue(obval + E_Chlorine) ;
#endif
#ifdef E_Air
  dui[E_Air      ] =  1.e-4 * ObVal_GetValue(obval + E_Air) ;
#endif

  
  
  for(i = 0 ; i < nn ; i++) {
    double* xi = ComputeVariables(el,u,u_n,f_n,0,dt,i) ;
    double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    int k ;
    
    if(!xi) return(-1) ;
    
    HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
    
    /* Derivation wrt LogC_CO2 -> relative value */
    #ifdef U_C_CO2
    dui[E_Carbon   ] =  1.e-4 * ObVal_GetRelativeValue(obval + E_Carbon,Un_Calcium(i)) ;
    #endif
    
    /* Derivation wrt LogC_Na -> relative value */
    #ifdef U_C_Na
    dui[E_Sodium   ] =  1.e-3 * ObVal_GetRelativeValue(obval + E_Sodium,Un_Sodium(i)) ;
    #endif
    
    /* Derivation wrt LogC_K -> relative value */
    #ifdef U_C_K
    dui[E_Potassium] =  1.e-3 * ObVal_GetRelativeValue(obval + E_Potassium,Un_Potassium(i)) ;
    #endif
    
    #if defined (U_ZN_Ca_S)
    dui[E_Calcium  ] =  1.e-4 * ObVal_GetAbsoluteValue(obval + E_Calcium,Un_Calcium(i)) ;
    dui[E_Calcium  ] *= ((xi[E_Calcium] > Un_Calcium(i)) ? 1 : -1) ;
    #elif defined (U_LogS_CH)
    dui[E_Calcium  ] =  1.e-4 * ObVal_GetAbsoluteValue(obval + E_Calcium,Un_Calcium(i)) ;
    #endif
    
    dui[E_Silicon  ] =  1.e-4 * ObVal_GetAbsoluteValue(obval + E_Silicon,Un_Silicon(i)) ;
    dui[E_Silicon  ] *= ((xi[E_Silicon] > Un_Silicon(i)) ? 1 : -1) ; 
    
    #ifdef E_eneutral
      #if defined (U_C_OH)
      /* Derivation wrt LogC_OH -> relative value */
      dui[E_eneutral ] =  1.e-2 * ObVal_GetRelativeValue(obval + E_eneutral,C_OHn(i)) ;
      #elif defined (U_Z_OH)
      dui[E_eneutral ] =  1.e-3 * ObVal_GetAbsoluteValue(obval + E_eneutral,Z_OHn(i)) * C_OHn(i) * fabs(dLogC_OHdZ_OH(Z_OHn(i))) ;
      #endif
    #endif
    
    /* Derivation wrt LogC_Cl -> relative value */
    #ifdef E_Chlorine
      #ifdef U_C_Cl
      dui[E_Chlorine ] =  1.e-3 * ObVal_GetRelativeValue(obval + E_Chlorine,C_Cln(i)) ;
      #endif
    #endif
    
    
    for(k = 0 ; k < NEQ ; k++) {
      double  dui_k = dui[k] ;
      double* dxi = ComputeVariableDerivatives(el,0,dt,xi,dui_k,k) ;
      
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
    
        cii[E_Carbon*NEQ    + k] = dxi[I_N_C] ;
        cii[E_Calcium*NEQ   + k] = dxi[I_N_Ca] ;
        cii[E_Sodium*NEQ    + k] = dxi[I_N_Na] ;
        cii[E_Silicon*NEQ   + k] = dxi[I_N_Si] ;
        cii[E_Potassium*NEQ + k] = dxi[I_N_K] ;
        cii[E_mass*NEQ      + k] = dxi[I_Mass] ;
#ifdef E_eneutral
        cii[E_eneutral*NEQ  + k] = dxi[I_N_Q] ;
#endif
#ifdef E_Chlorine
        cii[E_Chlorine*NEQ  + k] = dxi[I_N_Cl] ;
#endif
#ifdef E_Air
        cii[E_Air*NEQ       + k] = dxi[I_M_Air] ;
#endif
      }

      /* Transfer terms from node i to node j: d(wij)/d(ui_k) */
      {
        int j ;
        
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            
            {
              double* g = CementSolutionDiffusion_GetGradient(csd) ;
              double* muj = CementSolutionDiffusion_GetPotentialAtPoint(csd,j) ;
              int n = CementSolutionDiffusion_NbOfConcentrations ;
              int l ;
    
              /* On output ComputeVariableDerivatives has computed 
               * mui + d(mui) (through hcc). Then it is copied into muj */
              HardenedCementChemistry_CopyChemicalPotential(hcc,muj) ;

              /* The derivatives d(mui)/d(ui_k) */
              for(l = 0 ; l < n ; l++) {
                g[l] = (muj[l] - mui[l]) / dui_k ;
              }
            }
            
            {
              double* cij = c + (i*nn + j)*NEQ*NEQ ;
              double dij  = dist[nn*i + j] ;
              double dtdij = dt/dij ;
              double* dw = ComputeFluxes(el,dxi,i,j) ;
        
              cij[E_Carbon*NEQ    + k] = - dtdij*dw[I_W_C] ;
              cij[E_Calcium*NEQ   + k] = - dtdij*dw[I_W_Ca] ;
              cij[E_Sodium*NEQ    + k] = - dtdij*dw[I_W_Na] ;
              cij[E_Silicon*NEQ   + k] = - dtdij*dw[I_W_Si] ;
              cij[E_Potassium*NEQ + k] = - dtdij*dw[I_W_K] ;
              cij[E_mass*NEQ      + k] = - dtdij*dw[I_W_tot] ;
              //cij[E_charge*NEQ    + k] = - dtdij*dw[I_W_q] ;
              cij[E_charge*NEQ    + k] = - dw[I_W_q]/dij ;
#ifdef E_Chlorine
              cij[E_Chlorine*NEQ  + k] = - dtdij*dw[I_W_Cl] ;
#endif
#ifdef E_Air
              cij[E_Air*NEQ       + k] = - dtdij*dw[I_W_Air] ;
#endif
            }
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int n)
{
  double* v0 = Element_GetConstantTerm(el) ;
  double* x   = Variables[n] ;
  double* x_n = Variables_n(x) ;
  
  /* Primary Variables */
  x[E_Carbon   ] = LogC_CO2(n) ;
  x[E_Sodium   ] = LogC_Na(n) ;
  x[E_Potassium] = LogC_K(n) ;
  x[E_Calcium  ] = U_Calcium(n) ;
  x[E_Silicon  ] = U_Silicon(n) ;
  x[E_mass     ] = P_L(n) ;
  x[E_charge   ] = PSI(n) ;
#ifdef E_eneutral
  x[E_eneutral ] = LogC_OH(n) ;
#endif
#ifdef E_Chlorine
  x[E_Chlorine ] = LogC_Cl(n) ;
#endif
#ifdef E_Air
  x[E_Air      ] = P_G(n) ;
#endif
  
  /* Needed variables to compute secondary components */
  x_n[I_N_CH]  = N_CHn(n) ;
  x_n[I_N_CC]  = N_CCn(n) ;
  x_n[I_V_S0]  = V_S0(n) ;
  x_n[I_C_OH]  = C_OHn(n) ;
  
  {
    int k = ComputeSecondaryVariables(el,t,dt,x_n,x) ;
    
    if(k < 0) return(NULL) ;
  }
  
  return(x) ;
}


double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double du_i,int i)
{
  double* x_n = Variables_n(x) ;
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += du_i ;
  
  {
    int k = ComputeSecondaryVariables(el,t,dt,x_n,dx) ;
    
    if(k < 0) return(NULL) ;
  }
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= du_i ;
  }

  return(dx) ;
}



int  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  double logc_co2   = x[E_Carbon ] ;
  double u_calcium  = x[E_Calcium] ;
  double u_silicon  = x[E_Silicon] ;
  double p_l        = x[E_mass   ] ;
#ifdef E_Air
  double p_g        = x[E_Air    ] ; ;
#else
  double p_g        = p_g0 ;
#endif
#ifdef E_Chlorine
  double logc_cl    = x[E_Chlorine] ;
  double c_cl       = pow(10,logc_cl) ;
#else
  double logc_cl    = -99 ;
  double c_cl       = 1.e-99 ;
#endif
  
  
  /* Liquid components */
  double c_co2      = pow(10,logc_co2) ;
  double logc_co2aq = log(k_h) + logc_co2 ;
  //double c_co2aq    = k_h*c_co2 ;
    
  /* Solve cement chemistry */
  {
    double logc_na  = x[E_Sodium   ] ;
    double logc_k   = x[E_Potassium] ;
    //double logc_co2aq = log10(c_co2aq) ;
#ifdef E_eneutral
    double logc_oh  = x[E_eneutral] ;
#else
    double logc_oh  = log10(x_n[I_C_OH]) ;
#endif
    //double logc_oh  = log10(c_oh) ;
    double psi      = x[E_charge] ;

    #if defined (U_ZN_Ca_S)
    {
      double si_ch_cc = Log10SaturationIndexOfCcH(u_calcium) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CH_CC,si_ch_cc) ;
    }
    #elif defined (U_LogS_CH)
    {
      double si_ch = Log10SaturationIndexOfCH(u_calcium) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CH,si_ch) ;
    }
    #endif
        
    {
      double si_csh = Log10SaturationIndexOfCSH(u_silicon) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CSH,si_csh) ;
    }
  
    HardenedCementChemistry_SetInput(hcc,LogC_CO2,logc_co2aq) ;
    HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
    HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
    HardenedCementChemistry_SetInput(hcc,LogC_OH,logc_oh) ;
    HardenedCementChemistry_GetElectricPotential(hcc) = psi ;
    
    HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = c_cl ;
    HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,Cl) = logc_cl ;
  
    HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_CO2_H2O) ;

#ifndef E_eneutral
    {
      int k = HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      
      if(k < 0) return(-1) ;
    }
#endif
  }
  
  
  
  /* Backup */
  
  double c_q_l  = HardenedCementChemistry_GetLiquidChargeDensity(hcc) ;
  
  double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
  
  double c_c_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,C) ;
  double c_ca_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Ca) ;
  double c_na_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Na) ;
  double c_k_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
  double c_si_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Si) ;
  double c_cl_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Cl) ;
  
  double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
  double s_cc   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CC) ;
       
    
  /* Solid contents */
  /* ... as components: CH, CC, CSH */
  double n_chn      = x_n[I_N_CH] ;
  double n_ccn      = x_n[I_N_CC] ;
  double n_ch       = CHSolidContent(u_calcium,n_chn,n_ccn,s_ch,s_cc,dt) ;
  double n_cc       = CCSolidContent(u_calcium,n_chn,n_ccn,s_ch,s_cc,dt) ;
  double n_csh      = CSHSolidContent(u_silicon) ;
  
  /* ... as elements: C, Ca, Si, Cl */
  double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
  double n_si_s     = n_csh ;
  double n_ca_s     = n_ch + n_cc + x_csh * n_csh ;
  double n_c_s      = n_cc ;
  double n_cl_s     = n_csh * AdsorbedChloridePerUnitMoleOfCSH(c_cl,x_csh) ;
  
  /* ... as mass */
  double z_csh      = HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc) ;
  double m_csh      = MolarMassOfCSH(x_csh,z_csh) * n_csh ;
  double m_ch       = M_CaOH2 * n_ch ;
  double m_cc       = M_CaCO3 * n_cc ;
  double m_cl_s     = (M_Cl - M_OH) * n_cl_s ;
  double m_s        = m_ch + m_cc + m_csh + m_cl_s ;
  
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(x_csh) ;
  double v_s        = V_CH * n_ch + V_CC * n_cc + v_csh * n_csh ;
  
  
  /* Porosity */
  double v_s0     = x_n[I_V_S0] ;
  double phi_th   = phi0 + v_s0 - v_s ;
  double phi      = MAX(phi_th,phi_min) ;
  
  
  /* Pressures */
  double p_c      = p_g - p_l ;
  double p_v      = VaporPressure(p_l) ;
  double p_co2    = c_co2 * RT ;
  double p_air    = p_g - p_v - p_co2 ;
  
  
  /* Saturation */
  double s_l      = SaturationDegree(p_c) ;
  double s_g      = 1 - s_l ;
  
  
  /* Liquid contents */
  double phi_l    = phi*s_l ;
  /* ... as elements: C, Ca, Si */
  double n_c_l  = phi_l*c_c_l ;
  double n_ca_l = phi_l*c_ca_l ;
  double n_na_l = phi_l*c_na_l ;
  double n_k_l  = phi_l*c_k_l ;
  double n_si_l = phi_l*c_si_l ;
  double n_cl_l = phi_l*c_cl_l ;
  /* ... as charge */
  //double n_q_l  = phi_l*c_q_l ;
  /* ... as mass */
  double m_l    = phi_l*rho_l ;
       
       
  /* Gas contents */
  double phi_g  = phi * s_g ;
  /* ... as elements */
  double n_c_g  = phi_g * c_co2 ;
  /* ... as densities */
  double rho_co2_g   = M_CO2 * c_co2 ;
  double rho_h2o_g   = MassDensityOfWaterVapor(p_v) ;
  double rho_air_g   = M_AIR * p_air / RT ;
  double rho_noair_g = rho_h2o_g + rho_co2_g ;
  //double rho_g       = rho_air_g + rho_noair_g ;
  /* ... as masses */
  double m_air_g   = phi_g * rho_air_g ;
  double m_noair_g = phi_g * rho_noair_g ;
  //double m_g       = phi_g * rho_g ;


  /* Back up */
  

  /* Gas components */
  x[I_P_G       ] = p_g ;
  x[I_C_CO2     ] = c_co2 ;
  x[I_RHO_H2O_g ] = rho_h2o_g ;
  
  /* Liquid components */
  x[I_P_L       ] = p_l ;
  x[I_S_L       ] = s_l ;
  
  /* Solid components */
  x[I_N_CH    ] = n_ch ;
  x[I_V_S     ] = v_s ;
  x[I_N_Si_S  ] = n_si_s ;
  x[I_N_Ca_S  ] = n_ca_s ;
  x[I_N_CC    ] = n_cc ;
  x[I_N_CSH   ] = n_csh ;
  x[I_V_CSH   ] = v_csh ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  /* Element contents */
  x[I_N_C ]  = n_c_l  + n_c_s  + n_c_g ;
  x[I_N_Ca]  = n_ca_l + n_ca_s ;
  x[I_N_Na]  = n_na_l ; 
  x[I_N_K ]  = n_k_l  ;
  x[I_N_Si]  = n_si_l + n_si_s ;
  x[I_N_Cl]  = n_cl_l + n_cl_s ;
  
  /* Mass of dry air */
  x[I_M_Air] = m_air_g ;
  
  /* Total mass */
  x[I_Mass]  = m_noair_g + m_l + m_s ;
  
  /* Charge density */
  //x[I_N_Q]   = n_q_l ;
  x[I_N_Q]   = c_q_l ;
  
  /* Hydroxide ion concentration */
  x[I_C_OH]  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
  
    
  return(0) ;
}



int concentrations_oh_na_k(double c_co2,double u_calcium,double u_silicon,double c_cl,double c_na_tot,double c_k_tot)
{
/* Solve a set of 3 equations:
 * 1. Electroneutralilty
 * 2. Mass balance of Na
 * 3. Mass balance of K
 * Unknowns: c_oh, c_na, c_k.
 * On input, c_na_tot and c_k_tot are the total contents of Na and K
 */
  
  /* Initialization */
  double c_na = c_na_tot ;
  double c_k  = c_k_tot ;
  double c_oh0 = c_na + c_k ;
  double c_oh = c_oh0 ;
  
  /* c_na_tot =  c_na * (A_Na + B_Na*c_oh + C_Na*c_oh*c_oh) */
  //double A_Na = 1 ;
  //double B_Na = k_naoh/k_e + k_nahco3*k_h*c_co2/k_1 ;
  //double C_Na = k_naco3*k_h*c_co2/(k_1*k_e) ;

  /* c_k_tot =  c_k * (A_K + B_K*c_oh) */
  //double A_K = 1 ;
  //double B_K = k_koh/k_e  ;
  
  double err,tol = 1.e-8 ;
  

  /* Solve cement chemistry */
  {
    double c_co2aq    = k_h*c_co2 ;
    double logc_co2aq = log10(c_co2aq) ;
    double logc_cl    = log10(c_cl) ;

    #if defined (U_ZN_Ca_S)
    {
      double si_ch_cc = Log10SaturationIndexOfCcH(u_calcium) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CH_CC,si_ch_cc) ;
    }
    #elif defined (U_LogS_CH)
    {
      double si_ch = Log10SaturationIndexOfCH(u_calcium) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CH,si_ch) ;
    }
    #endif
        
    {
      double si_csh = Log10SaturationIndexOfCSH(u_silicon) ;
          
      HardenedCementChemistry_SetInput(hcc,SI_CSH,si_csh) ;
    }
  
    HardenedCementChemistry_SetInput(hcc,LogC_CO2,logc_co2aq) ;
    HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = c_cl ;
    HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,Cl) = logc_cl ;
  }
  
  int i = 0 ;
    
  
  do {
    double dc_oh = - c_oh ;
    double logc_na    = log10(c_na) ;
    double logc_k     = log10(c_k) ;
    
    HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
    HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
    HardenedCementChemistry_SetInput(hcc,LogC_OH,-7) ;
  
    HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_CO2_H2O) ;

    {
      int k = HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      
      if(k < 0) return(-1) ;
    }
    
    {
      double c_na_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Na) ;
      
      c_na *= c_na_tot/c_na_l ;
    }
    
    //c_na = c_na_tot/(A_Na + B_Na*c_oh + C_Na*c_oh*c_oh) ;
    
    {
      double c_k_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
      
      c_k *= c_k_tot/c_k_l ;
    }
    
    //c_k  = c_k_tot/(A_K + B_K*c_oh) ;
    
    c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
    
    dc_oh += c_oh ;
    
    err = fabs(dc_oh/c_oh) ;
    
    if(i++ > 20) {
      printf("c_na_tot = %e\n",c_na_tot) ;
      printf("c_na     = %e\n",c_na) ;
      printf("c_k_tot  = %e\n",c_k_tot) ;
      printf("c_k      = %e\n",c_k) ;
      printf("c_oh0    = %e\n",c_oh0) ;;
      printf("c_oh     = %e\n",c_oh) ;
      Message_Direct("concentrations_oh_na_k : non convergence") ;
      return(-1) ;
    }

  } while(err > tol || c_oh < 0) ;
  
  /*
  {
    printf("\n") ;
    printf("c_oh = %e \n", c_oh) ;
    printf("c_na = %e \n", c_na) ;
    printf("c_k  = %e \n", c_k) ;
    printf("c_na(kcc) = %e \n", HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Na)) ;
    printf("c_k(hcc)  = %e \n", HardenedCementChemistry_GetAqueousConcentrationOf(hcc,K)) ;
  }
  */
  
  return(0) ;
}


double PermeabilityCoefficient_KozenyCarman(Element_t* el,double phi)
/* Kozeny-Carman model */
{
  double coeff_permeability ;
  
  {
    double kozeny_carman  = (phi > 0) ? pow(phi/phi0,3.)*pow(((1 - phi0)/(1 - phi)),2.) : 0 ;
	
    coeff_permeability = kozeny_carman ;
  }
  
  return(coeff_permeability) ;
}



double PermeabilityCoefficient_VermaPruess(Element_t* el,double phi)
/* Ref:
 * A. Verma and K. Pruess,
 * Thermohydrological Conditions and Silica Redistribution Near High-Level
 * Nuclear Wastes Emplaced in Saturated Geological Formations,
 * Journal of Geophysical Research, 93(B2) 1159-1173, 1988
 * frac  = fractionnal length of pore bodies (0.8) 
 * phi_r = fraction of initial porosity (phi/phi0) at which permeability is 0 
 */
{
  double coeff_permeability ;
  
  {
    double phi_c = phi0 * phi_r ;
    //double w = 1 + (1/frac)/(1/phi_r - 1) ;
    double w = 1 + (phi_r/frac)/(1 - phi_r) ;
    double t = (phi - phi_c)/(phi0 - phi_c) ;
    double verma_pruess = (t > 0) ? t*t*(1 - frac + (frac/(w*w)))/(1 - frac + frac*(pow(t/(t + w - 1),2.))) : 0 ;

    coeff_permeability = verma_pruess ;
  }
  
  return(coeff_permeability) ;
}




double TortuosityToLiquid_OhJang(double phi,double s_l)
/* Ref:
 * Byung Hwan Oh, Seung Yup Jang, 
 * Prediction of diffusivity of concrete based on simple analytic equations, 
 * Cement and Concrete Research 34 (2004) 463 - 480.
 * tau = (m_p + sqrt(m_p**2 + phi_c/(1 - phi_c) * (Ds/D0)**(1/n)))**n
 * m_p = 0.5 * ((phi_cap - phi_c) + (Ds/D0)**(1/n) * (1 - phi_c - phi_cap)) / (1 - phi_c)
 */
{
  double phi_cap = (phi > 0) ? 0.5 * phi : 0  ;
  double phi_c = 0.17 ;          /* Percolation capilar porosity */
  double n     = 2.7 ;           /* OPC n  = 2.7  --------  Fly ash n  = 4.5 */
  double ds    = 1.e-4 ;         /* OPC ds = 1e-4 --------  Fly ash ds = 5e-5 */
  double dsn   = pow(ds,1/n) ;
  double m_phi = 0.5 * ((phi_cap - phi_c) + dsn * (1 - phi_c - phi_cap)) / (1 - phi_c) ;
  double tausat =  pow(m_phi + sqrt(m_phi*m_phi + dsn * phi_c/(1 - phi_c)),n) ;
  
  double tau =  tausat * pow(s_l,4.5) ;
    
  return(tau) ;
}




double TortuosityToLiquid_BazantNajjar(double phi,double s_l)
/* Ref:
 * Z. P. BAZANT, L.J. NAJJAR,
 * Nonlinear water diffusion in nonsaturated concrete,
 * Materiaux et constructions, 5(25), 1972.
 */
{
  double iff = 0.00029 * exp(9.95 * phi) ;
  double tausat = (iff < 1) ? iff : 1 ;
  double tau    = tausat / (1 + 625*pow((1 - s_l),4)) ;
    
  return(tau) ;
}



double TortuosityToLiquid_Xie(double phi,double s_l)
{
  double vca = 0.392 ;
  double vfa = 0.284 ;
  double vpa = 0.324  ;
  double fac = (1.-2.1*vca)/(1.-0.65*vfa)*vpa ;
  double hydrationdegree = 0.95 ;
  double wcratio = 0.5 ;
  double phiused = phi/vpa - 0.19*hydrationdegree/(wcratio + 0.32);
  double fporo = (phiused > 0.18) ? 0.001+0.07*phiused*phiused+1.8*(phiused-0.18)*(phiused-0.18) : 0.001+0.07*phiused*phiused ;
  double coefsa = pow(s_l,6.) ; 
  double iff = fac * fporo * coefsa ;
      
  return(iff) ;
}



double TortuosityToGas(double phi,double s_l)
{
  double s_g    = 1 - s_l ;
  double tausat = (phi > 0) ? pow(phi,1.74) : 0 ;
  double tau    = (s_g > 0) ? tausat * pow(s_g,3.20) : 0 ;

  return(tau) ;
}




double CHSolidContent_kin1(double n_chn,double s_ch,double dt)
{
  double av         = 1 - n_chn/n_ch0 ;
  double dn1sdt     = a_2*dn1_caoh2sdt(av,c_2) ;
  double dn_chsdt   = dn1sdt*log(s_ch) ; /* Kinetics */
  double n_ch_ki    = MAX(n_chn + dt*dn_chsdt , 0.) ;
  
  return(n_ch_ki) ;
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



double saturationdegree(double pc,double pc3,Curve_t* curve)
/* Saturation degree: regularization around 1 */
{
  double* x = Curve_GetXRange(curve) ;
  double* y = Curve_GetYValue(curve) ;
  double pc1 = x[0] ;
  double sl ;
  
  if(pc >= pc3 || pc1 >= pc3) {
    sl = Curve_ComputeValue(curve,pc) ;
  } else {
    double sl1 = y[0] ;
    double sl3 = Curve_ComputeValue(curve,pc3) ;
    
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  
  return(sl) ;
}
