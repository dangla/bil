/* General features of the model:
 * Curves for CSH:
 *   - C/S ratio
 *   - H/S ratio
 *   - Molar Volume
 * Alkalis (as sodium and potassium compounds)
 * Dissolution kinetics for CH based on spherical crystal 
 * coated by a calcite layer.
 * Dissolution and continuous decalcification of CSH
 * Precipitation/Dissolution of CC
 * Use of Zeta unknowns for Calcium and Silicon
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

#define TITLE   "Carbonation Of CBM with gas pressure (2018)"
#define AUTHORS "Morandeau-Thiery-Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ    	  (8)


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Nb of terms */
#define NVE    	  ((12 + CementSolutionDiffusion_NbOfConcentrations)*NN)
#define NVI       (8*NN*NN + 3*NN)
#define NV0       (2)


/* Indices of equations */
#define E_C       (0)
#define E_q       (1)
#define E_mass    (2)
#define E_Ca      (3)
#define E_Na      (5)
#define E_K       (6)
#define E_Si      (4)
#define E_Air     (7)



/* Indices of unknowns */
#define U_C_CO2   (0)
#define U_P_L     (2)
#define U_ZN_Ca_S (3)
#define U_PSI     (1)
#define U_C_Na    (5)
#define U_C_K     (6)
#define U_ZN_Si_S (4)
#define U_P_G     (7)


#define NOLOG_U   1
#define LOG_U     2
#define Ln10      Math_Ln10
#define U_CO2     LOG_U
#define U_Na      LOG_U
#define U_K       LOG_U


/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)





/* Names of nodal unknowns */
#if (U_CO2 == LOG_U)
  #define LogC_CO2(n)   (UNKNOWN(n,U_C_CO2))
  #define LogC_CO2n(n)  (UNKNOWNn(n,U_C_CO2))
  #define C_CO2(n)      (pow(10,UNKNOWN(n,U_C_CO2)))
  #define C_CO2n(n)     (pow(10,UNKNOWNn(n,U_C_CO2)))
#else
  #define C_CO2(n)      (UNKNOWN(n,U_C_CO2))
  #define C_CO2n(n)     (UNKNOWNn(n,U_C_CO2))
  #define LogC_CO2(n)   (log10(UNKNOWN(n,U_C_CO2)))
  #define LogC_CO2n(n)  (log10(UNKNOWNn(n,U_C_CO2)))
#endif

#define ZN_Si_S(n)    (UNKNOWN(n,U_ZN_Si_S))
#define ZN_Si_Sn(n)   (UNKNOWNn(n,U_ZN_Si_S))

#define ZN_Ca_S(n)    (UNKNOWN(n,U_ZN_Ca_S))
#define ZN_Ca_Sn(n)   (UNKNOWNn(n,U_ZN_Ca_S))

#define P_L(n)        (UNKNOWN(n,U_P_L))
#define P_Ln(n)       (UNKNOWNn(n,U_P_L))

#define P_G(n)        (UNKNOWN(n,U_P_G))
#define P_Gn(n)       (UNKNOWNn(n,U_P_G))

#define PSI(n)        (UNKNOWN(n,U_PSI))
#define PSIn(n)       (UNKNOWNn(n,U_PSI))

#if (U_Na == LOG_U)
  #define LogC_Na(n)    (UNKNOWN(n,U_C_Na))
  #define LogC_Nan(n)   (UNKNOWNn(n,U_C_Na))
  #define C_Na(n)       (pow(10,UNKNOWN(n,U_C_Na)))
  #define C_Nan(n)      (pow(10,UNKNOWNn(n,U_C_Na)))
#else
  #define C_Na(n)       (UNKNOWN(n,U_C_Na))
  #define C_Nan(n)      (UNKNOWNn(n,U_C_Na))
  #define LogC_Na(n)    (log10(UNKNOWN(n,U_C_Na)))
  #define LogC_Nan(n)   (log10(UNKNOWNn(n,U_C_Na)))
#endif

#if (U_K == LOG_U)
  #define LogC_K(n)     (UNKNOWN(n,U_C_K))
  #define LogC_Kn(n)    (UNKNOWNn(n,U_C_K))
  #define C_K(n)        (pow(10,UNKNOWN(n,U_C_K)))
  #define C_Kn(n)       (pow(10,UNKNOWNn(n,U_C_K)))
#else
  #define C_K(n)        (UNKNOWN(n,U_C_K))
  #define C_Kn(n)       (UNKNOWNn(n,U_C_K))
  #define LogC_K(n)     (log10(UNKNOWN(n,U_C_K)))
  #define LogC_Kn(n)    (log10(UNKNOWNn(n,U_C_K)))
#endif



/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define NW_C          (f)
#define NW_Cn         (f_n)
#define N_C(i)        MassAndFlux(NW_C,i,i)
#define N_Cn(i)       MassAndFlux(NW_Cn,i,i)
#define W_C(i,j)      MassAndFlux(NW_C,i,j)

#define NW_q          (f   + NN*NN)
#define NW_qn         (f_n + NN*NN)
#define N_q(i)        MassAndFlux(NW_q,i,i)
#define N_qn(i)       MassAndFlux(NW_qn,i,i)
#define W_q(i,j)      MassAndFlux(NW_q,i,j)

#define MW_tot        (f   + 2*NN*NN)
#define MW_totn       (f_n + 2*NN*NN)
#define M_tot(i)      MassAndFlux(MW_tot,i,i)
#define M_totn(i)     MassAndFlux(MW_totn,i,i)
#define W_tot(i,j)    MassAndFlux(MW_tot,i,j)

#define NW_Ca         (f   + 3*NN*NN)
#define NW_Can        (f_n + 3*NN*NN)
#define N_Ca(i)       MassAndFlux(NW_Ca,i,i)
#define N_Can(i)      MassAndFlux(NW_Can,i,i)
#define W_Ca(i,j)     MassAndFlux(NW_Ca,i,j)

#define NW_Na         (f   + 4*NN*NN)
#define NW_Nan        (f_n + 4*NN*NN)
#define N_Na(i)       MassAndFlux(NW_Na,i,i)
#define N_Nan(i)      MassAndFlux(NW_Nan,i,i)
#define W_Na(i,j)     MassAndFlux(NW_Na,i,j)

#define NW_K          (f   + 5*NN*NN)
#define NW_Kn         (f_n + 5*NN*NN)
#define N_K(i)        MassAndFlux(NW_K,i,i)
#define N_Kn(i)       MassAndFlux(NW_Kn,i,i)
#define W_K(i,j)      MassAndFlux(NW_K,i,j)

#define NW_Si         (f   + 6*NN*NN)
#define NW_Sin        (f_n + 6*NN*NN)
#define N_Si(i)       MassAndFlux(NW_Si,i,i)
#define N_Sin(i)      MassAndFlux(NW_Sin,i,i)
#define W_Si(i,j)     MassAndFlux(NW_Si,i,j)

#define MW_Air         (f   + 7*NN*NN)
#define MW_Airn        (f_n + 7*NN*NN)
#define M_Air(i)       MassAndFlux(MW_Air,i,i)
#define M_Airn(i)      MassAndFlux(MW_Airn,i,i)
#define W_Air(i,j)     MassAndFlux(MW_Air,i,j)

#define N_CH(i)       (f   + 8*NN*NN)[i]
#define N_CHn(i)      (f_n + 8*NN*NN)[i]

#define N_CC(i)       (f   + 8*NN*NN + NN)[i]
#define N_CCn(i)      (f_n + 8*NN*NN + NN)[i]

#define C_OH(i)       (f   + 8*NN*NN + 2*NN)[i]
#define C_OHn(i)      (f_n + 8*NN*NN + 2*NN)[i]




/* Names used for explicit terms */
#define TransferCoefficient(va,n)  ((va) + (n)*NN)

#define KF_CO2          TransferCoefficient(va,0)

#define KD_L            TransferCoefficient(va,1)
#define KC_C_L          TransferCoefficient(va,2)
#define KC_Ca_L         TransferCoefficient(va,3)
#define KC_Na_L         TransferCoefficient(va,4)
#define KC_K_L          TransferCoefficient(va,5)
#define KC_Si_L         TransferCoefficient(va,6)

#define KF_H2O          TransferCoefficient(va,7)

#define TORTUOSITY      TransferCoefficient(va,8)

#define KC_CO2          TransferCoefficient(va,9)
#define KC_H2O          TransferCoefficient(va,10)

#define KD_G            TransferCoefficient(va,11)

#define CONCENTRATION(i)  (TransferCoefficient(va,12) + (i)*CementSolutionDiffusion_NbOfConcentrations)



/* Names used for constant terms */
#define V_S0(n)         (v0[(0+n)])




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
#define M_H2O          MolarMassOfMolecule(H2O)
/* Molar volume of liquid water */
#define V_H2O          (18 * cm3)
/* Mass density */
#define MassDensityOfWaterVapor(p_v)   (M_H2O*(p_v)/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(p_l)          (exp(V_H2O/RT*((p_l) - p_l0)))
#define VaporPressure(p_l)             (p_v0*RelativeHumidity(p_l))
//#define LiquidPressure(hr)             (p_l0 + RT/V_H2O*(log(hr)))



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
//#define SaturationDegree(p)              (Curve_ComputeValue(SATURATION_CURVE,p))
#define SaturationDegree(p)              (saturation(p,p_c3,SATURATION_CURVE))
#define RELATIVEPERMLIQ_CURVE            (relativepermliqcurve)
#define RelativePermeabilityToLiquid(s)  (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,s))
#define RELATIVEPERMGAS_CURVE            (relativepermgascurve)
#define RelativePermeabilityToGas(s)     (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,s))
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



/* Calcium Hydroxide (Portlandite) Properties (CH)
 * ----------------------------------------------- */
#define M_CaOH2        MolarMassOfMolecule(CaO2H2)
/* Molar volume of CH solid (dm3/mole) */
#define V_CH           (33 * cm3)
//#define CHSolidContent(zn_ca_s)        CalciumContentInCH(zn_ca_s)
#define CHSolidContent(zn_ca_s,n_chn,s_ch,s_cc,dt) \
        (((s_cc) > (s_ch)) ? CHSolidContent_kin1(n_chn,s_ch,dt) : CalciumContentInCHAndCC(zn_ca_s))



/* Calcium Carbonate (Calcite) Properties (CC)
 * ------------------------------------------- */
#define M_CaCO3        MolarMassOfMolecule(CaCO3)
/* Molar volume of CC (dm3/mole) */
#define V_CC           (37 * cm3)
//#define CCSolidContent_kin(n,s,dt)     MAX((n + dt*r_cc*(s - 1)),0.)
//#define CCSolidContent(n,s,dt)         CCSolidContent_kin(n,s,dt)
#define CCSolidContent(zn_ca_s,n_chn,s_ch,s_cc,dt) \
        (CalciumContentInCHAndCC(zn_ca_s) - CHSolidContent(zn_ca_s,n_chn,s_ch,s_cc,dt))


/* Element contents in solid phases  */
#define n_ca_ref                           (n_ch0)
#define n_si_ref                           (n_csh0)
#define CalciumContentInCHAndCC(zn_ca_s)   (n_ca_ref*MAX(zn_ca_s,0.))
#define SiliconContentInCSH(zn_si_s)       (n_si_ref*MAX(zn_si_s,0.))
#define CalciumContentInCH(zn_ca_s)        (n_ca_ref*MAX(zn_ca_s,0.))



/* Access to Properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



static int     pm(const char* s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,double,int) ;
static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;

static double  dn1_caoh2sdt(double,double) ;
static double  CHSolidContent_kin1(double,double,double) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,double**,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;

static int     TangentCoefficients(Element_t*,double,double*) ;


static void    ComputePhysicoChemicalProperties(double) ;

static void    concentrations_oh_na_k(double,double,double,double,double) ;

static double  PermeabilityCoefficient_KozenyCarman(Element_t*,double) ;
static double  PermeabilityCoefficient_VermaPruess(Element_t*,double) ;
static double  TortuosityToLiquid_OhJang(double,double) ;
static double  TortuosityToLiquid_BazantNajjar(double,double) ;
static double  TortuosityToGas(double,double) ;

static double  saturation(double,double,Curve_t*) ;


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
static double a_2,c_2 ;
static double n_ch0,n_csh0,c_na0,c_k0 ;

static double p_c3 ;

static double p_g0 ;
static double p_l0 ;
static double p_v0 ;

static double d_co2 ;
static double d_vap ;

static double mu_l ;
static double mu_g ;

static double RT ;

static CementSolutionDiffusion_t* csd = NULL ;
static HardenedCementChemistry_t* hcc = NULL ;



#include "PhysicalConstant.h"
#include "AtmosphericPressure.h"
#include "WaterViscosity.h"
#include "AirViscosity.h"
#include "WaterVaporPressure.h"
#include "DiffusionCoefficientOfMoleculeInAir.h"


void ComputePhysicoChemicalProperties(double TK)
{

  /* Diffusion Coefficient Of Molecules In Air (dm2/s) */
  d_co2      = DiffusionCoefficientOfMoleculeInAir(CO2,TK) ;
  d_vap      = DiffusionCoefficientOfMoleculeInAir(H2O,TK) ;
  
  /* Viscosities */
  mu_l    = WaterViscosity(TK) ;
  mu_g    = AirViscosity(TK) ;
  
  /* Water vapor pressure */
  p_v0    = WaterVaporPressure(TK) ;
  
  /* Reference pressures */
  p_l0    = AtmosphericPressure ;
  p_g0    = AtmosphericPressure ;
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
}


#define NbOfVariables    (31)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;
#define  SecondaryVariables(x)  (x + NEQ)
#define  dSecondaryVariables    (dVariables + NEQ)


#define I_C_CO2        (30)
#define I_P_L          (8)

#define I_N_C          (9)
#define I_N_Ca         (10)
#define I_N_Si         (11)
#define I_N_K          (12)
#define I_N_Na         (13)
#define I_Mass         (14)
#define I_N_Q          (15)

#define I_N_Si_S       (16)
#define I_N_Ca_S       (17)

#define I_N_CH         (18)
#define I_N_CHn        (19)

#define I_N_CC         (20)
#define I_N_CCn        (21)

#define I_V_S          (22)
#define I_V_S0         (23)

#define I_Phi          (24)

#define I_V_CSH        (25)

#define I_C_OHn        (26)

#define I_RHO_H2O_g    (27)

#define I_M_Air        (28)

#define I_P_G          (29)


#define NbOfVariableFluxes    (8)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;

#define I_W_C           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_Na          (3)
#define I_W_K           (4)
#define I_W_tot         (5)
#define I_W_q           (6)
#define I_W_Air         (7)


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
  else if(strcmp(s,"kg_int") == 0)        return (16) ;
  else if(strcmp(s,"p_c3") == 0)          return (17) ;
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
  p_c3     = GetProperty("p_c3") ;
  
  saturationcurve  = Element_FindCurve(el,"s_l") ;
  relativepermliqcurve  = Element_FindCurve(el,"kl_r") ;
  relativepermgascurve  = Element_FindCurve(el,"kg_r") ;
  molarvolumeofcshcurve  = Element_FindCurve(el,"v_csh") ;
}


int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_C   ,"carbon") ;
  Model_CopyNameOfEquation(model,E_q   ,"charge") ;
  Model_CopyNameOfEquation(model,E_mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Ca  ,"calcium") ;
  Model_CopyNameOfEquation(model,E_Na  ,"sodium") ;
  Model_CopyNameOfEquation(model,E_K   ,"potassium") ;
  Model_CopyNameOfEquation(model,E_Si  ,"silicon") ;
  Model_CopyNameOfEquation(model,E_Air ,"air") ;
  
  
#if (U_CO2 == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_CO2 ,"logc_co2") ;
#else
  Model_CopyNameOfUnknown(model,U_C_CO2 ,"c_co2") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Si_S,"z_si") ;
  Model_CopyNameOfUnknown(model,U_P_L    ,"p_l") ;
  Model_CopyNameOfUnknown(model,U_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,U_PSI    ,"psi") ;
#if (U_Na == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_Na   ,"logc_na") ;
#else
  Model_CopyNameOfUnknown(model,U_C_Na   ,"c_na") ;
#endif
#if (U_K == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_K    ,"logc_k") ;
#else
  Model_CopyNameOfUnknown(model,U_C_K    ,"c_k") ;
#endif
  Model_CopyNameOfUnknown(model,U_P_G    ,"p_g") ;
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 17 ;
  
  InternationalSystemOfUnits_UseAsLength("decimeter") ;
  InternationalSystemOfUnits_UseAsMass("hectogram") ;

  Material_ScanProperties(mat,datafile,pm) ;
    
  /* Default initialization */
  {
    kg_int = Material_GetProperty(mat)[pm("kg_int")] ;
    kl_int = Material_GetProperty(mat)[pm("kl_int")] ;

    if(kg_int == 0.) {
      Material_GetProperty(mat)[pm("kg_int")] = kl_int ;
    }
  }
  
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
  
  {
    phi_r = Material_GetProperty(mat)[pm("phi_r")] ;
    if(phi_r == 0) phi_r = 0.7 ;
    
    Material_GetProperty(mat)[pm("phi_r")] = phi_r ;
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
        arret("ReadMatProp: no kl_r - s_l curve") ;
      }

      if((i = Curves_FindCurveIndex(curves,"kg_r")) < 0) {
        arret("ReadMatProp: no kg_r - s_l curve") ;
      }

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
  printf("\t- Total mass balance     (mass)\n") ;
  printf("\t- Charge balance         (charge)\n") ;
  printf("\t- Mass balance of air    (air)\n") ;
  
  printf("\n") ;
  printf("The 7 primary unknowns are:\n") ;
  printf("\t- Liquid pressure                  (p_l)\n") ;
  printf("\t- Gas pressure                     (p_g)\n") ;
  printf("\t- Electric potential               (psi) \n") ;
  printf("\t- Carbon dioxide gas concentration (c_co2 or logc_co2)\n") ;
  printf("\t- Potassium concentration          (c_k or logc_k)\n") ;
  printf("\t- Sodium concentration             (c_na or logc_na)\n") ;
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
  printf("\t mass      : hg !\n") ;
  printf("\t pressure  : Pa !\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;


  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;
  fprintf(ficd,"kl_int = 1.4e-17   # Intrinsic permeability (dm2)\n") ;
  fprintf(ficd,"kg_int = 1.4e-17   # Intrinsic permeability (dm2)\n") ;
  fprintf(ficd,"N_CH = 3.9        # Initial content in Ca(OH)2 (mol/L)\n") ;
  fprintf(ficd,"Radius_CH = 40.e-5  # Portlandite crystal radius \n") ;
  fprintf(ficd,"N_CSH = 2.4        # Initial content in CSH (mol/L)\n") ;
  fprintf(ficd,"C_Na = 0.019      # Total content in Na (mol/L)\n") ;
  fprintf(ficd,"C_K  = 0.012      # Total content in K  (mol/L)\n") ;
  fprintf(ficd,"D = 7.e-15        # Diffusion coef in CC (dm/mol/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Kinetic coef 2 (dm/mol/s)\n") ;
  fprintf(ficd,"frac = 0.8        # Fractionnal length of pore bodies\n") ;
  fprintf(ficd,"phi_r = 0.7       # Porosity for which permeability vanishes\n") ;
  fprintf(ficd,"Curves = my_file  # File name: p_c S_l k_rl C/S H/S V_csh\n") ;  

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = (Element_IsSubmanifold(el)) ? 0 : NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
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
      double c_na       = c_na_tot ;
      double c_k        = c_k_tot ;
      double c_co2      = C_CO2(i) ;
      double zn_ca_s    = ZN_Ca_S(i) ;
      double zn_si_s    = ZN_Si_S(i) ;

      /* Compute the concentrations of alkalis Na and K */
      concentrations_oh_na_k(c_co2,zn_ca_s,zn_si_s,c_na_tot,c_k_tot) ;
  
      c_na = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Na) ;
      c_k  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,K) ;


#if (U_Na == LOG_U)
      LogC_Na(i)  = log10(c_na) ;
#else
      C_Na(i)     = c_na ;
#endif
#if (U_K == LOG_U)
      LogC_K(i)   = log10(c_k) ;
#else
      C_K(i)      = c_k ;
#endif
      
      /* pH */
      C_OH(i) = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
      
      /* Solid contents */
      {
        double s_ch       = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
        double s_cc       = HardenedCementChemistry_GetSaturationIndexOf(hcc,CC) ;
        double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
        double n_ch       = (s_cc > s_ch) ? 0       : n_ch_cc ;
        double n_cc       = n_ch_cc - n_ch ;
        double n_si_s     = SiliconContentInCSH(zn_si_s) ;
        double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
        double v_csh      = MolarVolumeOfCSH(x_csh) ;
        double v_s0       = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;
        
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
      double* x   = ComputeVariables(el,u,f,0,0,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
    
      /* Back up */
      N_C(i)   = x[I_N_C] ;
      N_Ca(i)  = x[I_N_Ca] ;
      N_Na(i)  = x[I_N_Na] ;
      N_Si(i)  = x[I_N_Si] ;
      N_K(i)   = x[I_N_K] ; 
      M_tot(i) = x[I_Mass] ;
      N_q(i)   = x[I_N_Q] ;
      M_Air(i) = x[I_M_Air] ;

      /* Solid contents */
      N_CH(i) = x[I_N_CH] ;
      N_CC(i) = x[I_N_CC] ;
      
      /* pH */
      C_OH(i) = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
    }
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el,u,f) ;


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
        W_Air(i,j)   = w[I_W_Air] ;
        
        W_C(j,i)     = - w[I_W_C] ;
        W_Ca(j,i)    = - w[I_W_Ca] ;
        W_Na(j,i)    = - w[I_W_Na] ;
        W_Si(j,i)    = - w[I_W_Si] ;
        W_q(j,i)     = - w[I_W_q] ;
        W_K(j,i)     = - w[I_W_K] ;
        W_tot(j,i)   = - w[I_W_tot] ;
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
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Molar contents */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x   = ComputeVariables(el,u,f_n,t,dt,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
    
      /* Back up */
      N_C(i)   = x[I_N_C] ;
      N_Ca(i)  = x[I_N_Ca] ;
      N_Na(i)  = x[I_N_Na] ;
      N_Si(i)  = x[I_N_Si] ;
      N_K(i)   = x[I_N_K] ; 
      M_tot(i) = x[I_Mass] ;
      N_q(i)   = x[I_N_Q] ;
      M_Air(i) = x[I_M_Air] ;

      /* Solid contents */
      N_CH(i) = x[I_N_CH] ;
      N_CC(i) = x[I_N_CC] ;
      
      /* pH */
      C_OH(i) = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;

      {
        double c_co2      = x[I_C_CO2] ;
        
        double c_h2o = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,H2O) ;
        double c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
        double c_ca= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Ca) ;
        double c_na= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Na) ;
        double c_k = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,K) ;
        
        double n_si_s     = x[I_N_Si_S] ;
        double n_ch       = x[I_N_CH] ;
      
        if(c_co2 < 0 || c_oh <= 0 || c_h2o <= 0 || c_na < 0 || c_k < 0 || c_ca < 0 || n_si_s < 0. || n_ch < 0.) {
          double x0 = Element_GetNodeCoordinate(el,i)[0] ;
          double n_cc       = x[I_N_CC] ;
          double c_naoh    	= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaOH) ;
          double c_nahco3  	= HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaHCO3) ;
          double c_naco3 	  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,NaCO3) ;
          printf("\n") ;
          printf("en x     = %e\n",x0) ;
          printf("c_co2    = %e\n",c_co2) ;
          printf("c_oh     = %e\n",c_oh) ;
          printf("c_h2o    = %e\n",c_h2o) ;
          printf("n_cc     = %e\n",n_cc) ;
          printf("c_na     = %e\n",c_na) ;
          printf("c_k      = %e\n",c_k) ;
          printf("c_ca     = %e\n",c_ca) ;
          printf("n_si_s   = %e\n",n_si_s) ;
          printf("c_naoh   = %e\n",c_naoh) ;
          printf("c_nahco3 = %e\n",c_nahco3) ;
          printf("c_naco3  = %e\n",c_naco3) ;
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
        W_Air(i,j)   = w[I_W_Air] ;
        
        W_C(j,i)     = - w[I_W_C] ;
        W_Ca(j,i)    = - w[I_W_Ca] ;
        W_Na(j,i)    = - w[I_W_Na] ;
        W_Si(j,i)    = - w[I_W_Si] ;
        W_q(j,i)     = - w[I_W_q] ;
        W_K(j,i)     = - w[I_W_K] ;
        W_tot(j,i)   = - w[I_W_tot] ;
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
  
  TangentCoefficients(el,dt,c) ;
  {
    double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }


/* On output TangentCoefficients has computed the derivatives wrt
 * C_CO2, C_Na and C_K
 * (see ComputeVariables and ComputeVariablesDerivatives). */
 
#if (U_CO2 == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_CO2)     *= Ln10*C_CO2(0) ;
      K(i,U_C_CO2+NEQ) *= Ln10*C_CO2(1) ;
    }
  }
#endif

#if (U_Na == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_Na)     *= Ln10*C_Na(0) ;
      K(i,U_C_Na+NEQ) *= Ln10*C_Na(1) ;
    }
  }
#endif

#if (U_K == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_K)     *= Ln10*C_K(0) ;
      K(i,U_C_K+NEQ) *= Ln10*C_K(1) ;
    }
  }
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
      R(i,E_C) -= r1[i] ;
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
          g[i*nn + j] = dt * W_q(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_q) -= r1[i] ;
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
      R(i,E_Ca) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element Na: (N_Na - N_Nan) + dt * div(W_Na) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Na,NW_Nan,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Na) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element K: (N_K - N_Kn) + dt * div(W_K) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_K,NW_Kn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_K) -= r1[i] ;
    }
  }

  /*
    Conservation of element Si: (N_Si - N_Sin) + dt * div(W_Si) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Si,NW_Sin,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Si) -= r1[i] ;
    }
  }
  
  /*
    Conservation of dry air mass: (M_Air - M_Airn) + dt * div(W_Air) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,MW_Air,MW_Airn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Air) -= r1[i] ;
    }
  }

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nso = 58 ;
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
    double* x = ComputeVariables(el,u,f,t,0,j) ;
    
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
    {
      double p_l        = x[I_P_L] ;
      double p_g        = x[I_P_G] ;
      double p_c        = p_g - p_l ;
      double s_l        = SaturationDegree(p_c) ;
      
      Result_Store(r + i++,&s_l,"saturation",1) ;
    }
    
    Result_Store(r + i++,x + I_Phi,"porosity",1) ;
    
    /* Concentration in gas phase */
    Result_Store(r + i++,x + I_C_CO2,"c_co2",1) ;
    
    /* Element concentrations in liquid phase */
    Result_Store(r + i++,ptEC(Ca ),"c_ca_l",1) ;
    Result_Store(r + i++,ptEC(Si ),"c_si_l",1) ;
    Result_Store(r + i++,ptEC(Na ),"c_na_l",1) ;
    Result_Store(r + i++,ptEC(K  ),"c_k_l" ,1) ;
    Result_Store(r + i++,ptEC(C  ),"c_c_l" ,1) ;
    
    /* Portlandite */
    Result_Store(r + i++,x + I_N_CH,"n_CH",1) ;
    Result_Store(r + i++,ptS(CH),"s_ch",1) ;
    
    /* C-S-H */
    Result_Store(r + i++,x + I_N_Si_S,"n_CSH",1) ;
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
    
    /* Total element contents */
    Result_Store(r + i++,x + I_N_Ca,"n_Ca",1) ;
    Result_Store(r + i++,x + I_N_Si,"n_Si",1) ;
    Result_Store(r + i++,x + I_N_Na,"n_Na",1) ;
    Result_Store(r + i++,x + I_N_K ,"n_K" ,1) ;
    Result_Store(r + i++,x + I_N_C ,"n_C" ,1) ;
    
    
    /* Miscellaneous */
    {
      double CsurS = x[I_N_Ca_S]/x[I_N_Si_S] ;
      
      Result_Store(r + i++,&CsurS,"CsurS",1) ;
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
    
    /* Added by A. Morandeau */

    /* Transferts */
    {
      double p_l        = x[I_P_L] ;
      double p_g        = x[I_P_G] ;
      double p_c        = p_g - p_l ;
      double phi        = x[I_Phi] ;
      double s_l        = SaturationDegree(p_c) ;
      double coeff_permeability = PermeabilityCoefficient(el,phi) ;
	    double k_l  = (kl_int/mu_l)*RelativePermeabilityToLiquid(s_l)*coeff_permeability ;
      
      Result_Store(r + i++,&k_l,"k_l",1) ;
      Result_Store(r + i++,&coeff_permeability,"verma-pruess",1) ;
    }
    {
      double n_chn      = x[I_N_CHn] ;
      double av         = 1 - n_chn/n_ch0 ;
      double dn1sdt     = a_2*dn1_caoh2sdt(av,c_2) ;
      double s_ch       = *ptS(CH) ;
      double dn_chsdt   = dn1sdt*log(s_ch) ;
      double coeff_dnCH = log(s_ch) ;
    
      Result_Store(r + i++,&dn_chsdt,"dn_chsdt",1) ;
      Result_Store(r + i++,&dn1sdt,"dn1sdt",1) ;
      Result_Store(r + i++,&coeff_dnCH,"coeff_dnCH",1) ;
    }
    {
      double v_solide_csh   = x[I_V_CSH] * x[I_N_Si_S] ;
      
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
      double p_co2      = x[I_C_CO2] * RT ;
      double p_air      = p_g - p_v - p_co2 ;
      
      Result_Store(r + i++,&p_air,"air pressure",1) ;
      Result_Store(r + i++,&p_v,"vapor pressure",1) ;
      Result_Store(r + i++,&p_co2,"CO2 pressure",1) ;
      Result_Store(r + i++,&p_g,"gas pressure",1) ;
    }
  }
  
  
  if(i != nso) arret("ComputeOutputs") ;
  return(nso) ;
}


void ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Termes explicites (va)  */
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    i ; 

  /* initialization */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  /*
    Transfer coefficients
  */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x = ComputeVariables(el,u,f,0,0,i) ;
    /* pressures */
    double p_l     = x[I_P_L] ;
    double p_g     = x[I_P_G] ;
    double p_c     = p_g - p_l ;
    /* saturation */
    double s_l     = SaturationDegree(p_c) ;

    /* Porosity */
    double phi        = x[I_Phi] ;
	
    /* Permeabilities */
    double coeff_permeability = PermeabilityCoefficient(el,phi) ;
    double k_l  = (kl_int/mu_l)*RelativePermeabilityToLiquid(s_l)*coeff_permeability ;
    double k_g  = (kg_int/mu_g)*RelativePermeabilityToGas(s_l)*coeff_permeability ;
    
    
    /* Transport in gas phase (diffusion coef) */
    {
      double s_g   = 1 - s_l ;
      double phi_g = phi * s_g ;
      /* tortuosity gas */
      double taugas  = TortuosityToGas(phi,s_l) ;
      
      KF_CO2[i]      = phi_g * taugas * d_co2 ;
      KF_H2O[i]      = phi_g * taugas * d_vap ;
    }
      
    /* Transport in gas phase (advection coef) */
    {
      double n_co2_g   = x[I_C_CO2] ;
      double p_co2_g   = n_co2_g * RT ;
      double p_v       = VaporPressure(p_l) ;
      double p_air     = p_g - p_v - p_co2_g ;
      double rho_co2_g = M_CO2 * n_co2_g ;
      double rho_h2o_g = x[I_RHO_H2O_g] ;
      double rho_air   = MassDensityOfDryAir(p_air) ;
      double rho_g     = rho_air + rho_h2o_g + rho_co2_g ;
      
      KD_G[i]        = rho_g * k_g ;
      
      KC_CO2[i]      = n_co2_g / rho_g ;
      KC_H2O[i]      = rho_h2o_g / rho_g ;
    }
    
    /* Liquid tortuosity */
    {
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
    
    /* Transport in liquid phase */
    {
  
      double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
      double c_c_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,C) ;
      double c_ca_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Ca) ;
      double c_na_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Na) ;
      double c_k_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
      double c_si_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Si) ;
    
      KD_L[i]    = rho_l * k_l ;
            
      KC_C_L[i]  = c_c_l  / rho_l ;
      KC_Ca_L[i] = c_ca_l / rho_l ;
      KC_Na_L[i] = c_na_l / rho_l ;
      KC_K_L[i]  = c_k_l  / rho_l ;
      KC_Si_L[i] = c_si_l / rho_l ;
    }
  }
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
    /* Mass flux of liquid */
    double grd_p_l = grdij[I_P_L] ;
    double kd_l    = 0.5 * (KD_L[i] + KD_L[j]) ;
    double w_l     = - kd_l * grd_p_l  ;

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
    }
  
    /* Advection in the cement solution */
    {
      /* Transfer terms */
      double kc_c_l   = 0.5 * (KC_C_L[i]   + KC_C_L[j]) ;
      double kc_ca_l  = 0.5 * (KC_Ca_L[i]  + KC_Ca_L[j]) ;
      double kc_si_l  = 0.5 * (KC_Si_L[i]  + KC_Si_L[j]) ;
      double kc_na_l  = 0.5 * (KC_Na_L[i]  + KC_Na_L[j]) ;
      double kc_k_l   = 0.5 * (KC_K_L[i]   + KC_K_L[j]) ;
    
      w[I_W_tot]   = w_l  ;
   
      /* Molar fluxes */
      w[I_W_C  ]  += kc_c_l  * w_l  ;
      w[I_W_Ca ]  += kc_ca_l * w_l  ;
      w[I_W_Si ]  += kc_si_l * w_l  ;
      w[I_W_Na ]  += kc_na_l * w_l  ;
      w[I_W_K  ]  += kc_k_l  * w_l  ;
    }
  }
  
  /* Transport in gas phase */
  {
    /* Advection and diffusion in gas phase */
    {
      /* Mass flux of gas */
      double grd_p_g = grdij[I_P_G] ;
      double kd_g    = 0.5 * (KD_G[i]     + KD_G[j]) ;
      double w_g     = - kd_g * grd_p_g  ;
      
      /* Molar flux of CO2 */
      double grd_co2 = grdij[I_C_CO2] ;
      double kf_co2  = 0.5 * (KF_CO2[i]   + KF_CO2[j]) ;
      double j_co2_g = - kf_co2 * grd_co2  ;
      double kc_co2  = 0.5 * (KC_CO2[i]   + KC_CO2[j]) ;
      //double kc_co2  = (w_g > 0) ? KC_CO2[MIN(i,j)] : KC_CO2[MAX(i,j)] ;
      double w_co2_g =   kc_co2 * w_g + j_co2_g  ;
      
      /* Mass flux of water vapor */
      double grd_h2o = grdij[I_RHO_H2O_g] ;
      double kf_h2o  = 0.5 * (KF_H2O[i]   + KF_H2O[j]) ;
      double j_h2o_g = - kf_h2o * grd_h2o  ;
      double kc_h2o  = 0.5 * (KC_H2O[i]   + KC_H2O[j]) ;
      //double kc_h2o  = (w_g > 0) ? KC_H2O[MIN(i,j)] : KC_H2O[MAX(i,j)] ;
      double w_h2o_g =   kc_h2o * w_g + j_h2o_g  ;
      
      /* Mass flux of dry air */
      double w_air   =   w_g - w_h2o_g - M_CO2 * w_co2_g ;
  
      w[I_W_Air]   =  w_air ;
      w[I_W_tot]  +=  w_g - w_air ;
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

  dui[U_C_CO2  ] =  1.e-4 * ObVal_GetValue(obval + U_C_CO2) ;
  dui[U_C_Na   ] =  1.e-3 * ObVal_GetValue(obval + U_C_Na) ;
  dui[U_C_K    ] =  1.e-3 * ObVal_GetValue(obval + U_C_K) ;
  dui[U_ZN_Ca_S] =  1.e-4 * ObVal_GetValue(obval + U_ZN_Ca_S) ;
  dui[U_ZN_Si_S] =  1.e-4 * ObVal_GetValue(obval + U_ZN_Si_S) ;
  dui[U_P_L    ] =  1.e-4 * ObVal_GetValue(obval + U_P_L) ;
  dui[U_PSI    ] =  1.e+0 * ObVal_GetValue(obval + U_PSI) ;
  dui[U_P_G    ] =  1.e-4 * ObVal_GetValue(obval + U_P_G) ;
  
  
  for(i = 0 ; i < nn ; i++) {
    double* xi = ComputeVariables(el,u,f_n,0,dt,i) ;
    double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    int k ;
    
    HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;

    dui[U_C_CO2  ] =  1.e-4 * ObVal_GetValue(obval + U_C_CO2) ;
    dui[U_C_Na   ] =  1.e-3 * ObVal_GetValue(obval + U_C_Na) ;
    dui[U_C_K    ] =  1.e-3 * ObVal_GetValue(obval + U_C_K) ;
    dui[U_ZN_Ca_S] =  1.e-4 * ObVal_GetValue(obval + U_ZN_Ca_S) ;
    dui[U_ZN_Si_S] =  1.e-4 * ObVal_GetValue(obval + U_ZN_Si_S) ;
    
    dui[U_ZN_Si_S] *= ((xi[U_ZN_Si_S] > ZN_Si_Sn(i)) ? 1 : -1) ; 
    dui[U_ZN_Ca_S] *= ((xi[U_ZN_Ca_S] > ZN_Ca_Sn(i)) ? 1 : -1) ;
    
    #if (U_CO2 == LOG_U)
    dui[U_C_CO2  ] *=  C_CO2n(i) ;
    #endif
    #if (U_Na == LOG_U)
    dui[U_C_Na   ] *=  C_Nan(i) ;
    #endif
    #if (U_K == LOG_U)
    dui[U_C_K    ] *=  C_Kn(i) ;
    #endif
    
    for(k = 0 ; k < NEQ ; k++) {
      double  dui_k = dui[k] ;
      double* dxi = ComputeVariableDerivatives(el,0,dt,xi,dui_k,k) ;
      
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
    
        cii[E_C*NEQ    + k] = dxi[I_N_C] ;
        cii[E_Ca*NEQ   + k] = dxi[I_N_Ca] ;
        cii[E_Na*NEQ   + k] = dxi[I_N_Na] ;
        cii[E_Si*NEQ   + k] = dxi[I_N_Si] ;
        cii[E_K*NEQ    + k] = dxi[I_N_K] ;
        cii[E_mass*NEQ + k] = dxi[I_Mass] ;
        cii[E_Air*NEQ  + k] = dxi[I_M_Air] ;
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
        
              cij[E_C*NEQ    + k] = - dtdij*dw[I_W_C] ;
              cij[E_Ca*NEQ   + k] = - dtdij*dw[I_W_Ca] ;
              cij[E_Na*NEQ   + k] = - dtdij*dw[I_W_Na] ;
              cij[E_Si*NEQ   + k] = - dtdij*dw[I_W_Si] ;
              cij[E_K*NEQ    + k] = - dtdij*dw[I_W_K] ;
              cij[E_mass*NEQ + k] = - dtdij*dw[I_W_tot] ;
              cij[E_q*NEQ    + k] = - dtdij*dw[I_W_q] ;
              cij[E_Air*NEQ  + k] = - dtdij*dw[I_W_Air] ;
            }
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double* f_n,double t,double dt,int n)
{
  double* v0 = Element_GetConstantTerm(el) ;
  double* x = Variables[n] ;
  
  /* Primary Variables */
  x[U_C_CO2  ] = C_CO2(n) ;
  x[U_C_Na   ] = C_Na(n) ;
  x[U_C_K    ] = C_K(n) ;
  x[U_ZN_Ca_S] = ZN_Ca_S(n) ;
  x[U_ZN_Si_S] = ZN_Si_S(n) ;
  x[U_P_L    ] = P_L(n) ;
  x[U_PSI    ] = PSI(n) ;
  x[U_P_G    ] = P_G(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn]  = N_CHn(n) ;
  x[I_N_CCn]  = N_CCn(n) ;
  x[I_V_S0 ]  = V_S0(n) ;
  x[I_C_OHn]  = C_OHn(n) ;
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}


double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double du_i,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += du_i ;
  
  ComputeSecondaryVariables(el,t,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= du_i ;
  }

  return(dx) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  double c_co2      = x[U_C_CO2  ] ;
  double zn_ca_s    = x[U_ZN_Ca_S] ;
  double zn_si_s    = x[U_ZN_Si_S] ;
  double p_l        = x[U_P_L    ] ;
  double p_g        = x[U_P_G    ] ;
  
  
  /* Liquid components */
  double c_co2aq    = k_h*c_co2 ;
    
  /* Solve cement chemistry */
  {
    double c_na       = x[U_C_Na   ] ;
    double c_k        = x[U_C_K    ] ;
    double logc_co2aq = log10(c_co2aq) ;
    double logc_na    = log10(c_na) ;
    double logc_k     = log10(c_k) ;
    double c_oh       = x[I_C_OHn  ] ;
    double logc_oh    = log10(c_oh) ;
    double psi        = x[U_PSI] ;
  
    HardenedCementChemistry_SetInput(hcc,SI_CH_CC,MIN(zn_ca_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,SI_CSH,MIN(zn_si_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,LogC_CO2,logc_co2aq) ;
    HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
    HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
    HardenedCementChemistry_SetInput(hcc,LogC_OH,logc_oh) ;
    HardenedCementChemistry_GetElectricPotential(hcc) = psi ;
  
    HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_CO2_H2O) ;
      
    HardenedCementChemistry_SolveElectroneutrality(hcc) ;
  }
  
  
  
  /* Backup */
  
  double c_q  = HardenedCementChemistry_GetLiquidChargeDensity(hcc) ;
  
  double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
  
  double c_c_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,C) ;
  double c_ca_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Ca) ;
  double c_na_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Na) ;
  double c_k_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
  double c_si_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Si) ;
  
  double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
  double s_cc   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CC) ;
       
    
  /* Solid contents */
  /* ... as components: CH, CC, CSH */
  double n_ch_cc    = CalciumContentInCHAndCC(zn_ca_s) ;
  double n_chn      = x[I_N_CHn] ;
  double n_ch_ki    = CHSolidContent_kin1(n_chn,s_ch,dt) ;
  //double n_ch       = (s_cc > s_ch) ? n_ch_ki : n_ch_cc ;
  double n_ch       = CHSolidContent(zn_ca_s,n_chn,s_ch,s_cc,dt) ;
  //double n_cc       = n_ch_cc - n_ch ;
  double n_cc       = CCSolidContent(zn_ca_s,n_chn,s_ch,s_cc,dt) ;
  
  /* ... as elements: C, Ca, Si */
  //double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
  double n_si_s     = SiliconContentInCSH(zn_si_s) ;
  double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  double n_c_s      = n_cc ;
  
  /* ... as mass */
  //double z_csh      = WaterSiliconRatioInCSH(s_ch) ;
  double z_csh      = HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc) ;
  double m_csh      = MolarMassOfCSH(x_csh,z_csh) ;
  double m_s        = M_CaOH2*n_ch + M_CaCO3*n_cc + m_csh*n_si_s ;
  
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(x_csh) ;
  double v_s        = V_CH*n_ch + V_CC*n_cc + v_csh*n_si_s ;
  
  
  /* Porosity */
  double v_s0     = x[I_V_S0] ;
  double phi_th   = phi0 + v_s0 - v_s ;
  double phi      = MAX(phi_th,phi_min) ;
  
  
  /* Pressures */
  double p_c      = p_g - p_l ;
  double p_v      = VaporPressure(p_l) ;
  double p_co2    = c_co2 * RT ;
  double p_air    = p_g - p_v - p_co2 ;
  
  
  /* Saturations */
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
  /* ... as mass */
  double m_l    = phi_l*rho_l ;
       
       
  /* Gas contents */
  double phi_g  = phi * s_g ;
  /* ... as elements */
  double n_c_g  = phi_g * c_co2 ;
  /* ... as densities */
  double rho_co2_g = M_CO2 * c_co2 ;
  double rho_h2o_g = MassDensityOfWaterVapor(p_v) ;
  double rho_air   = M_AIR * p_air / RT ;
  double rho_g     = rho_air + rho_h2o_g + rho_co2_g ;
  /* ... as masses */
  double m_air     = phi_g * rho_air ;
  double m_g       = phi_g * rho_g ;


  /* Back up */
  

  /* Gas components */
  x[I_P_G       ] = p_g ;
  x[I_C_CO2     ] = c_co2 ;
  x[I_RHO_H2O_g ] = rho_h2o_g ;
  
  /* Liquid components */
  x[I_P_L       ] = p_l ;
  
  /* Solid components */
  x[I_N_CH    ] = n_ch ;
  x[I_V_S     ] = v_s ;
  x[I_N_Si_S  ] = n_si_s ;
  x[I_N_Ca_S  ] = n_ca_s ;
  x[I_N_CC    ] = n_cc ;
  x[I_V_CSH   ] = v_csh ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  /* Element contents */
  x[I_N_C ]  = n_c_l  + n_c_s  + n_c_g ;
  x[I_N_Ca]  = n_ca_l + n_ca_s ;
  x[I_N_Na]  = n_na_l ; 
  x[I_N_K ]  = n_k_l  ;
  x[I_N_Si]  = n_si_l + n_si_s ;
  
  /* Mass of dry air */
  x[I_M_Air] = m_air ;
  
  /* Total mass */
  x[I_Mass]  = m_g + m_l + m_s - m_air ;
  
  /* Charge density */
  x[I_N_Q]   = c_q ;
  
    
  return ;
}



void concentrations_oh_na_k(double c_co2,double zn_ca_s,double zn_si_s,double c_na_tot,double c_k_tot)
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
  
    HardenedCementChemistry_SetInput(hcc,SI_CH_CC,MIN(zn_ca_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,SI_CSH,MIN(zn_si_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,LogC_CO2,logc_co2aq) ;
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
      
    HardenedCementChemistry_SolveElectroneutrality(hcc) ;
    
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
      arret("concentrations_oh_na_k : non convergence") ;
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
    double w = 1 + (1/frac)/(1/phi_r - 1) ;
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
  double phi_c = 0.17 ;         /* Percolation capilar porosity */
  double n     = 2.7 ; 		      /* OPC n  = 2.7  --------  Fly ash n  = 4.5 */
  double ds    = 1.e-4 ;	      /* OPC ds = 1e-4 --------  Fly ash ds = 5e-5 */
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



double saturation(double pc,double pc3,Curve_t* curve)
/* Degre de saturation regularise autour de 1 */
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
