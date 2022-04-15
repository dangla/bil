/* General features of the model:
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* The Finite Volume Method */
#include "FEM.h"

/* Cement chemistry */
#include "HardenedCementChemistry.h"
#include "CementSolutionDiffusion.h"

/* Damage model */
#include "Damage.h"


#define TEMPERATURE   (293)

#define TITLE   "Internal/External sulfate attack of concrete (2020)" 
#define AUTHORS "Gu-Dangla"

#include "PredefinedMethods.h"



/* Indices of equations/unknowns */
enum {
  E_Sulfur,
  E_charge,
  E_Calcium,
  E_Potassium,
  E_Aluminium,
  /* Uncomment/comment the two next lines to consider/suppress electroneutrality */
  E_eneutral,
  #define E_eneutral E_eneutral
  E_kinetics,
  #define E_kinetics E_kinetics
  /* Uncomment/comment the two next lines to consider/suppress mechanics */
  //E_Mech,
  //#define E_Mech E_Mech
  E_Last
} ;



/* Nb of equations */
#define NbOfEquations      (E_Last + dim - 1)
#define NEQ                NbOfEquations





/* Generic names of nodal unknowns */
#define U_Sulfur     E_Sulfur
#define U_charge     E_charge
#define U_Calcium    E_Calcium
#define U_Potassium  E_Potassium
#define U_Aluminium  E_Aluminium
#ifdef E_eneutral
#define U_eneutral   E_eneutral
#endif
#ifdef E_kinetics
#define U_kinetics   E_kinetics
#endif
#ifdef E_Mech
#define U_Mech       E_Mech
#endif




/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */

/* Sulfur: unknown either C_SO4 or C_H2SO4 or LogC_SO4 or LogC_H2SO4 */
//#define U_C_SO4         U_Sulfur
#define U_LogC_SO4      U_Sulfur
//#define U_C_H2SO4       U_Sulfur
//#define U_LogC_H2SO4    U_Sulfur

/* charge: */
#define U_PSI       U_charge

/* Calcium:
 * - U_ZN_Ca_S: dissolution kinetics of CH; Cc at equilibrium 
 * - U_LogS_CH: dissolution kinetics of CH; precipitation kinetics of Cc */
#define U_ZN_Ca_S   U_Calcium
//#define U_LogS_CH     U_Calcium

/* Potassium: unknown either C or logC */
//#define U_C_K       U_Potassium
#define U_LogC_K    U_Potassium

/* Aluminium:
 * - U_ZN_Al_S: dissolution kinetics of AH3; Cc at equilibrium
 * - U_LogS_CH: dissolution kinetics of CH; precipitation kinetics of Cc */
//#define U_C_Al       U_Aluminium
//#define U_LogC_Al    U_Aluminium
#define U_ZN_Al_S    U_Aluminium

/* Electroneutrality: unknown either C_OH, logC_OH or Z_OH = C_H - C_OH */
#ifdef E_eneutral
//#define U_C_OH      U_eneutral
#define U_LogC_OH   U_eneutral
//#define U_Z_OH      U_eneutral
#endif

/* Crystal growth rate at the pore wall interface */
#ifdef E_kinetics
#define U_P_C      U_kinetics
#endif

/* Mechanics */
#ifdef E_Mech
#define U_Dis   U_Mech
#endif







/* Compiling options */
#define EXPLICIT  1
#define IMPLICIT  2
#define ELECTRONEUTRALITY   IMPLICIT



/* Nb of terms per point */
#define NVI   (54)    /*  nb of implicit terms per point */
#define NVE   (1 + CementSolutionDiffusion_NbOfConcentrations)     /*  nb of explicit terms per point */
#define NV0   (1)     /*  nb of constant terms per point */


/* We define some names for implicit terms (vi must be used as pointer below) */
#define NW_S        (vi)
#define NW_Sn       (vi_n)
#define N_S         (NW_S)[0]
#define N_Sn        (NW_Sn)[0]
#define W_S         (NW_S + 1)

#define NW_q        (vi + 4)
#define NW_qn       (vi_n + 4)
#define N_q         (NW_q)[0]
#define N_qn        (NW_qn)[0]
#define W_q         (NW_q + 1)

#define NW_Ca       (vi + 8)
#define NW_Can      (vi_n + 8)
#define N_Ca        (NW_Ca)[0]
#define N_Can       (NW_Can)[0]
#define W_Ca        (NW_Ca + 1)

#define NW_K        (vi + 12)
#define NW_Kn       (vi_n + 12)
#define N_K         (NW_K)[0]
#define N_Kn        (NW_Kn)[0]
#define W_K         (NW_K + 1)

#define NW_Si       (vi + 16)
#define NW_Sin      (vi_n + 16)
#define N_Si        (NW_Si)[0]
#define N_Sin       (NW_Sin)[0]
#define W_Si        (NW_Si + 1)

#define NW_Al       (vi + 20)
#define NW_Aln      (vi_n + 20)
#define N_Al        (NW_Al)[0]
#define N_Aln       (NW_Aln)[0]
#define W_Al        (NW_Al + 1)

#define NW_Cl       (vi + 24)
#define NW_Cln      (vi_n + 24)
#define N_Cl        (NW_Cl)[0]
#define N_Cln       (NW_Cln)[0]
#define W_Cl        (NW_Cl + 1)

#define Stress      (vi + 28) /* this a 3D tensor */
#define Stress_n    (vi_n + 28) /* this a 3D tensor */

#define N_CH        (vi   + 37)[0]
#define N_CHn       (vi_n + 37)[0]

#define N_CSH2      (vi   + 38)[0]
#define N_CSH2n     (vi_n + 38)[0]

#define N_AH3       (vi   + 39)[0]
#define N_AH3n      (vi_n + 39)[0]

#define N_AFm       (vi   + 40)[0]
#define N_AFmn      (vi_n + 40)[0]

#define N_AFt       (vi   + 41)[0]
#define N_AFtn      (vi_n + 41)[0]

#define N_C3AH6     (vi   + 42)[0]
#define N_C3AH6n    (vi_n + 42)[0]

#define PHI         (vi   + 43)[0]
#define PHIn        (vi_n + 43)[0]

#define PHI_C       (vi   + 44)[0]
#define PHI_Cn      (vi_n + 44)[0]


#ifndef E_eneutral
  #define C_OH      (vi   + 45)[0]
  #define C_OHn     (vi_n + 45)[0]
#endif

#define PoreRadius  (vi   + 46)[0]
#define PoreRadiusn (vi_n + 46)[0]

#define Beta_p      (vi   + 47)[0]
#define Beta_pn     (vi_n + 47)[0]

#define Straind     (vi   + 48)[0]
#define Straind_n   (vi_n + 48)[0]

#define VarPHI_C    (vi   + 49)[0]
#define VarPHI_Cn   (vi_n + 49)[0]

#define Hardv       (vi   + 50)[0]
#define Hardv_n     (vi_n + 50)[0]

#define Damage      (vi   + 51)[0]
#define Damage_n    (vi_n + 51)[0]

#define Crit        (vi   + 52)[0]
#define Crit_n      (vi_n + 52)[0]

#define GrowthRate   (vi   + 53)[0]
#define GrowthRate_n (vi_n + 53)[0]




/* We define some names for explicit terms (ve must be used as pointer below) */
#define TORTUOSITY       (ve)[0]

#define CONCENTRATION    (ve + 1) // CementSolutionDiffusion_NbOfConcentrations




/* We define some names for constant terms (v0 must be used as pointer below) */
#define V_Cem0           (v0)[0]





/* Math constants */
#define Ln10      Math_Ln10




#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define m     (InternationalSystemOfUnits_OneMeter)
#define m3    (m*m*m)
#define dm    (0.1*m)
#define cm    (0.01*m)
#define nm    (1.e-9*m)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)
#define Pa    (InternationalSystemOfUnits_OnePascal)
#define MPa   (1.e6*Pa)
#define GPa   (1.e9*Pa)
#define J     (Pa*m3)
#define mol   (InternationalSystemOfUnits_OneMole)




/* Material properties
 * ------------------- */
#define SATURATION_CURVE           (satcurve)
#define LiquidSaturationDegree(r)  (Curve_ComputeValue(SATURATION_CURVE,r))
#define dLiquidSaturationDegree(r) (Curve_ComputeDerivative(SATURATION_CURVE,r))
#define PoreEntryRadiusMax         (Curve_GetXRange(SATURATION_CURVE)[1])
#define DamageCoefficient(straind) \
        (1 - strain0/straind*exp(-(straind - strain0)/strainf))



/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CSH2 = Calcium Sulfate Dihydrate (Gypsum)
  CSH  = Calcium Silicates Hydrate
  SH   = Amorphous Silica Gel
*/


/* Calcium Silicate Hydrate Properties (C-S-H)
 * ------------------------------------------- */
//#define MOLARVOLUMEOFCSH_CURVE           (Element_GetCurve(el) + 2)
//#define MolarVolumeOfCSH(q)    (Curve_ComputeValue(MOLARVOLUMEOFCSH_CURVE,q))
#define V_CSH                  (78 * cm3)
#define V_SH                   (43 * cm3)
#define MolarVolumeOfCSH(x)    ((x)/1.7*V_CSH + (1 - (x)/1.7)*V_SH)
#define CSHSolidContent(zn_si_s)       SiliconContentInCSH(zn_si_s)



/* Calcium Hydroxide (Portlandite) Properties (CH)
 * ----------------------------------------------- */
/* Molar volume of CH solid */
#define V_CH       (33 * cm3)
#define CHSolidContent(zn_ca_s)        CalciumContentInCH(zn_ca_s)



/* Gypsum (CSH2) Properties
 * ------------------------ */
/* Molar volume of CSH2 crystal */
#define V_CSH2     (75 * cm3)
#define CSH2SolidContent_kin(n,s,dt)     MAX((n + dt*r_csh2*(s - 1)),0.)
#define CSH2SolidContent(n,s,dt)         CSH2SolidContent_kin(n,s,dt)



/* Gibbsite Properties (AH3)
 * ------------------------- */
/* Molar volume of AH3 solid */
#define V_AH3      (64.44 * cm3)
#define AH3SolidContent(zn_al_s)    (0.5*AluminiumContentInAH3(zn_al_s))



/* Monosulfoaluminate Properties (AFm = C4ASH12)
 * --------------------------------------------- */
/* Molar volume of AFm solid */
#define V_AFm      (311.26 * cm3)      /* Thermochimie (ANDRA) */
//#define AFmSolidContent(n,s,dt)     (n*pow(s,dt/t_afm))
#define AFmSolidContent(n,s,dt)     MAX((n + dt*r_afm*(s - 1)),0.)



/* Ettringite Properties (AFt = C6AS3H32)
 * -------------------------------------- */
/* Molar volume of AFt solid */
#define V_AFt      (710.32 * cm3)      /* Thermochimie (ANDRA) */
//#define AFtSolidContent(n,s,dt)     (n*pow(s,dt/t_aft))
#define AFtSolidContent(n,s,dt)     MAX((n + dt*r_aft*(s - 1)),0.)
/* Surface tension (N/m) */
#define Gamma_AFt   (0.1*Pa*m)



/* Sulfate adsorption curve 
 * ------------------------ */
#define AdsorbedSulfatePerUnitMoleOfCSH(c_so4,c_oh) \
        (alphacoef * (c_so4) / ((c_so4) + betacoef * (1.)))
//        (alphacoef * (c_so4) / ((c_so4) + betacoef * (c_oh)))



/* Crystal properties 
 * ------------------ */
#define V_Crystal       V_AFt
#define Gamma_Crystal   Gamma_AFt
/* Crystallization pressure */
#define CrystallizationPressure(beta) \
        RT/V_Crystal*log(beta)
#define dCrystallizationPressure(beta) \
        RT/V_Crystal/(beta)
/* Inverse of crystallization pressure */
#define InverseOfCrystallizationPressure(pc) \
        exp(V_Crystal/RT*pc)
        


/* Crystal growth kinetics */
#define PoreCrystalGrowthRate(s_c,beta,beta_p) \
        ((s_c) * CrystalGrowthRate(ap_AFt,dp_AFt,(beta_p)/(beta)))

#define dPoreCrystalGrowthRate(s_c,beta,beta_p) \
        ((s_c) * dCrystalGrowthRate(ap_AFt,dp_AFt,(beta_p)/(beta)) / (beta))

#define InterfaceCrystalGrowthRate(beta,beta_i) \
        (CrystalGrowthRate(ai_AFt,di_AFt,(beta_i)/(beta)))
        
#define dInterfaceCrystalGrowthRate(beta,beta_i) \
        (dCrystalGrowthRate(ai_AFt,di_AFt,(beta_i)/(beta)) / (beta))

#define CrystalGrowthRate(crys,diss,ateb) \
        ((((ateb) < 1) ? (crys) : (diss)) * (1 - (ateb)))
#define dCrystalGrowthRate(crys,diss,ateb) \
        ((((ateb) < 1) ? (crys) : (diss)) * (-1))

/* Inverse law of crystal growth kinetics */
#define InverseOfPoreCrystalGrowthRate(growth) \
        InverseOfCrystalGrowthRate(ap_AFt,dp_AFt,growth)
        
#define InverseOfCrystalGrowthRate(crys,diss,growth) \
        (1 - (growth) / (((growth) > 0) ? (crys) : (diss)))


/* Crystal - liquid interface
 * -------------------------- */
/* Interface equilibrium saturation index */
#define InterfaceEquilibriumSaturationIndex(r) \
        (exp(2*Gamma_Crystal*V_Crystal/(RT*r)))

#define dInterfaceEquilibriumSaturationIndex(r) \
        (-2*Gamma_Crystal*V_Crystal/(RT*r*r)*InterfaceEquilibriumSaturationIndex(r))

#define InverseOfInterfaceEquilibriumSaturationIndex(b) \
        (2*Gamma_Crystal*V_Crystal/(RT*log(b)))



/* Hydrogarnet Properties (C3AH6)
 * ------------------------------ */
/* Molar volume of C3AH6 solid */
#define V_C3AH6      (149.52 * cm3)
#define C3AH6SolidContent(n,s,dt)     MAX((n + dt*r_c3ah6*(s - 1)),0.)



/* Element contents in solid phases  */
//#define CalciumContentInCHAndCSH2(zn_ca_s) (n_ca_ref*MAX(zn_ca_s,0.))
#define CalciumContentInCH(zn_ca_s)        (n_ca_ref*MAX(zn_ca_s,0.))
#define SiliconContentInCSH(zn_si_s)       (n_si_ref*MAX(zn_si_s,0.))
#define AluminiumContentInAH3(zn_al_s)     (n_al_ref*MAX(zn_al_s,0.))


/* Gypsum-based porous material properties */
/* Porosity of gypsum-based porous material (-) */
#define PHI_Gyp    (0.85)
/* Molar volume of gypsum-based porous material */
#define V_Gyp      (V_CSH2/(1 - PHI_Gyp))



/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])



/* Fonctions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;


static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;
//static double* ComputeVariableFluxDerivatives(Element_t*,int,int,double) ;

static int     TangentCoefficients(Element_t*,double,double,double*) ;


static double  Radius(double,double,double,Element_t*) ;


static double TortuosityToLiquid_OhJang(double) ;
static double TortuosityToLiquid_BazantNajjar(double) ;

static void   ComputePhysicoChemicalProperties(void) ;

static double PoreWallEquilibriumSaturationIndex(double,double,double,double) ;
//static double PoreWallEquilibriumSaturationIndex_old(double,double,double,double,double,double,double) ;
static double ElasticDamageStress(double,double) ;
static double dElasticDamageStress(double,double) ;
static double DamageStrain(double,double) ;

//#define LiquidTortuosity  TortuosityToLiquid_OhJang
#define LiquidTortuosity  TortuosityToLiquid_BazantNajjar




#define ComputeFunctionGradients(...)     Damage_ComputeFunctionGradients(damage,__VA_ARGS__)
#define ReturnMapping(...)                Damage_ReturnMapping(damage,__VA_ARGS__)
#define CopyStiffnessTensor(...)          Damage_CopyStiffnessTensor(damage,__VA_ARGS__)
#define UpdateTangentStiffnessTensor(...) Damage_UpdateTangentStiffnessTensor(damage,__VA_ARGS__)



/* Parameters */
static double phi0 ;
static double phimin = 0.01 ;
static double r_afm,r_aft,r_c3ah6,r_csh2 ;
static double n_ca_ref,n_si_ref,n_al_ref ;
static double n_afm_0,n_aft_0,n_c3ah6_0,n_csh2_0 ;
static double ai_AFt,di_AFt ;
static double RT ;
static Curve_t* satcurve ;
static double Biot ;
//static double K_s ;
//static double G_s ;
//static double N_Biot ;
//static double G_Biot ;
static double strain0 ;
static double strainf ;
static double ap_AFt,dp_AFt ;
static double alphacoef ;
static double betacoef ;
static double r0 ;
static double* sig0 ;
static double  hardv0 ;
static double* cijkl ;
static double K_bulk ;
static Damage_t* damage ;
static int     damagemodel ;

#define  SetDamageModel(I) \
         do { \
           if(damagemodel < 0) { \
             damagemodel = I ; \
           } else if(damagemodel != I) { \
             Message_FatalError("Incompatible model") ; \
           } \
         } while(0)


#include "PhysicalConstant.h"
#include "Temperature.h"

void ComputePhysicoChemicalProperties(void)
{
  RT = PhysicalConstant_PerfectGasConstant * TEMPERATURE ;
}


static CementSolutionDiffusion_t* csd = NULL ;
static HardenedCementChemistry_t* hcc = NULL ;

enum VariableIndexes_e {
I_Dis = 0,
I_Dis2 = I_Dis + 2,
I_U_Calcium,
I_U_Sulfur,
I_U_Aluminium,
I_U_Potassium,
I_U_eneutral,
I_U_charge,
I_U_kinetics,

I_Strain,
I_Strain8 = I_Strain + 8,

I_GRD_U_Sulfur,
I_GRD_U_Sulfur2    = I_GRD_U_Sulfur    + 2,
I_GRD_U_Calcium,
I_GRD_U_Calcium2   = I_GRD_U_Calcium   + 2,
I_GRD_U_Potassium,
I_GRD_U_Potassium2 = I_GRD_U_Potassium + 2,
I_GRD_U_Aluminium,
I_GRD_U_Aluminium2 = I_GRD_U_Aluminium + 2,
I_GRD_U_charge,
I_GRD_U_charge2    = I_GRD_U_charge    + 2,
I_GRD_U_eneutral,
I_GRD_U_eneutral2  = I_GRD_U_eneutral  + 2,
I_GRD_U_kinetics,
I_GRD_U_linetics2  = I_GRD_U_kinetics  + 2,

I_N_q,
I_N_S,
I_N_Ca,
I_N_Si,
I_N_K,
I_N_Cl,
I_N_Al,

I_Stress,
I_Stress8 = I_Stress + 8,

I_W_S,
I_W_S2  = I_W_S  + 2,
I_W_Ca,
I_W_Ca2 = I_W_Ca + 2,
I_W_K,
I_W_K2  = I_W_K  + 2,
I_W_q,
I_W_q2  = I_W_q  + 2,
I_W_Al,
I_W_Al2 = I_W_Al + 2,

I_N_CH,
I_N_CSH2,
I_N_AH3,
I_N_AFm,
I_N_AFt,
I_N_C3AH6,
I_N_CSH,

I_PHI,
I_PHI_C,

I_V_CSH,
I_V_Cem,
I_V_Cem0,

I_C_OH,

I_PoreRadius,

I_S_C,

I_P_C,

I_Beta_p,

I_Straind,

I_VarPHI_C,

I_Hardv,

I_Damage,

I_Crit,

I_GrowthRate,

I_Last
} ;


#define NbOfVariables    (I_Last)
static double Variable[NbOfVariables] ;
static double Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;






int pm(const char *s)
{
       if(strcmp(s,"porosity") == 0)   return (0) ;
  else if(strcmp(s,"N_CH") == 0)       return (1) ;
  else if(strcmp(s,"N_Si") == 0)       return (2) ;
  else if(strcmp(s,"N_CSH") == 0)      return (2) ; /* synonym */
  else if(strcmp(s,"T_CH") == 0)       return (3) ;
  else if(strcmp(s,"T_CSH2") == 0)     return (4) ;
  else if(strcmp(s,"N_CSH2") == 0)     return (5) ;
  else if(strcmp(s,"N_AH3") == 0)      return (6) ;
  else if(strcmp(s,"N_AFm") == 0)      return (7) ;
  else if(strcmp(s,"N_AFt") == 0)      return (8) ;
  else if(strcmp(s,"N_C3AH6") == 0)    return (9) ;
  else if(strcmp(s,"T_AFm") == 0)      return (10) ;
  else if(strcmp(s,"T_AFt") == 0)      return (11) ;
  else if(strcmp(s,"R_AFm") == 0)      return (12) ;
  else if(strcmp(s,"R_AFt") == 0)      return (13) ;
  else if(strcmp(s,"R_C3AH6") == 0)    return (14) ;
  else if(strcmp(s,"R_CSH2") == 0)     return (15) ;
  else if(strcmp(s,"A_i") == 0)        return (16) ;
  else if(strcmp(s,"A_p") == 0)        return (17) ;
  else if(strcmp(s,"AlphaCoef") == 0)  return (18) ;
  else if(strcmp(s,"BetaCoef") == 0)   return (19) ;
  else if(strcmp(s,"r0") == 0)         return (20) ;
  else if(strcmp(s,"B_i") == 0)        return (21) ;
  else if(strcmp(s,"B_p") == 0)        return (22) ;
  else if(!strcmp(s,"poisson"))        return (23) ;
  else if(!strcmp(s,"young"))          return (24) ;
  else if(strcmp(s,"Biot") == 0)       return (25) ;
  else if(!strcmp(s,"sig0"))           return (26) ;
  else if(!strncmp(s,"sig0_",5)) {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(26 + 3*i + j) ;
    
    /* Pay attention below! 
     * This is the location of the initial hardening variable (hardv0)
     * (eg the initial pre-consolidation pressure). */
  } else if(!strcmp(s,"hardv0")) {
    return(35) ;
    
    /* Model 1: Mazars */
  } else if(!strcmp(s,"max_elastic_strain")) {
    SetDamageModel(1) ;
    return(35) ;
  } else if(!strcmp(s,"A_c")) {
    SetDamageModel(1) ;
    return(36) ;
  } else if(!strcmp(s,"A_t")) {
    SetDamageModel(1) ;
    return(37) ;
  } else if(!strcmp(s,"B_c")) {
    SetDamageModel(1) ;
    return(38) ;
  } else if(!strcmp(s,"B_t")) {
    SetDamageModel(1) ;
    return(39) ;
    
    /* Model 2: Marigo-Jirasek */
  } else if(!strcmp(s,"uniaxial_tensile_strength")) {
    SetDamageModel(2) ;
    return(35) ;
  } else if(!strcmp(s,"fracture_energy")) {
    SetDamageModel(2) ;
    return(36) ;
  } else if(!strcmp(s,"crack_band_width")) {
    SetDamageModel(2) ;
    return(37) ;
    
  } else if(strcmp(s,"Strain0") == 0)  return (35) ;
  else if(strcmp(s,"Strainf") == 0)    return (36) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  phi0      = GetProperty("porosity") ;
  n_ca_ref  = GetProperty("N_CH") ;
  n_si_ref  = GetProperty("N_CSH") ;
  n_al_ref  = GetProperty("N_AH3") ;
  n_csh2_0  = GetProperty("N_CSH2") ;
  n_afm_0   = GetProperty("N_AFm") ;
  n_aft_0   = GetProperty("N_AFt") ;
  n_c3ah6_0 = GetProperty("N_C3AH6") ;
  r_afm     = GetProperty("R_AFm") ;
  r_aft     = GetProperty("R_AFt") ;
  r_c3ah6   = GetProperty("R_C3AH6") ;
  r_csh2    = GetProperty("R_CSH2") ;
  ai_AFt    = GetProperty("A_i") ;
  strain0   = GetProperty("Strain0") ;
  strainf   = GetProperty("Strainf") ;
  ap_AFt    = GetProperty("A_p") ;
  alphacoef = GetProperty("AlphaCoef") ;
  betacoef  = GetProperty("BetaCoef") ;
  r0        = GetProperty("r0") ;
  Biot      = GetProperty("Biot") ;
  //G_s       = GetProperty("G_s") ;
  //K_s       = K_bulk/(1-Biot) ;
  //N_Biot    = K_s/(Biot - phi0) ;
  //G_Biot    = G_s/(0.75*phi0) ;
  di_AFt    = GetProperty("B_i") ;
  dp_AFt    = GetProperty("B_p") ;
  sig0      = &GetProperty("sig0") ;
  //hardv0    = GetProperty("hardv0") ;
  
  damage = Element_FindMaterialData(el,Damage_t,"Damage") ;
  {
    double young = GetProperty("young") ;
    double poisson = GetProperty("poisson") ;
    
    K_bulk = young / (3 - 6*poisson) ;
    cijkl   = Damage_GetStiffnessTensor(damage) ;
    hardv0  = Damage_GetHardeningVariable(damage)[0] ;
  }
  
  satcurve  = Element_FindCurve(el,"S_r") ;
}


int SetModelProp(Model_t* model)
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_Sulfur, "sulfur") ;
  Model_CopyNameOfEquation(model,E_Calcium,"calcium") ;
  Model_CopyNameOfEquation(model,E_charge, "charge") ;
  Model_CopyNameOfEquation(model,E_Potassium, "potassium") ;
  Model_CopyNameOfEquation(model,E_Aluminium,"aluminium") ;
#ifdef E_eneutral
  Model_CopyNameOfEquation(model,E_eneutral,"electroneutrality") ;
#endif
#ifdef E_kinetics
  Model_CopyNameOfEquation(model,E_kinetics,"kinetics") ;
#endif
#ifdef E_Mech
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
#endif


  /** Names of the main (nodal) unknowns */
#if   defined (U_LogC_H2SO4)
  Model_CopyNameOfUnknown(model,E_Sulfur,"logc_h2so4") ;
#elif defined (U_C_H2SO4)
  Model_CopyNameOfUnknown(model,E_Sulfur,"c_h2so4") ;
#elif defined (U_LogC_SO4)
  Model_CopyNameOfUnknown(model,E_Sulfur,"logc_so4") ;
#elif defined (U_C_SO4)
  Model_CopyNameOfUnknown(model,E_Sulfur,"c_so4") ;
#endif

  Model_CopyNameOfUnknown(model,E_Calcium,"z_ca") ;
  Model_CopyNameOfUnknown(model,E_charge,    "psi") ;
  
#if   defined (U_LogC_K)
  Model_CopyNameOfUnknown(model,E_Potassium,    "logc_k") ;
#elif defined (U_C_K)
  Model_CopyNameOfUnknown(model,E_Potassium,    "c_k") ;
#endif

  Model_CopyNameOfUnknown(model,E_Aluminium,"z_al") ;
  
#ifdef E_eneutral
  #if   defined (U_LogC_OH)
    Model_CopyNameOfUnknown(model,E_eneutral, "logc_oh") ;
  #elif defined (U_C_OH)
    Model_CopyNameOfUnknown(model,E_eneutral, "c_oh") ;
  #endif
#endif
#ifdef E_kinetics
  Model_CopyNameOfUnknown(model,E_kinetics,"p_c") ;
#endif
#ifdef E_Mech
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,E_Mech + i,name_unk[i]) ;
  }
#endif
  
  
  Model_GetComputePropertyIndex(model) = pm ;
  Model_GetNbOfVariables(model) = NbOfVariables ;
  
  Model_GetSequentialIndexOfUnknown(model)[E_Sulfur] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_Calcium] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_Potassium] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_Aluminium] = 0 ;
  Model_GetSequentialIndexOfUnknown(model)[E_charge] = 0 ;
#ifdef E_eneutral
  Model_GetSequentialIndexOfUnknown(model)[E_eneutral] = 0 ;
#endif
#ifdef E_kinetics
  Model_GetSequentialIndexOfUnknown(model)[E_kinetics] = 0 ;
#endif
#ifdef E_Mech
  for(i = 0 ; i < dim ; i++) {
    Model_GetSequentialIndexOfUnknown(model)[E_Mech+i] = 1 ;
  }
#endif
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 40 ;

  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("N_CH")]   = 1 ;
    Material_GetProperty(mat)[pm("N_Si")]   = 1 ;
    Material_GetProperty(mat)[pm("N_AH3")]  = 1 ;
    Material_GetProperty(mat)[pm("N_CSH2")] = 0 ;
    Material_GetProperty(mat)[pm("N_AFm")]  = 0 ;
    Material_GetProperty(mat)[pm("N_AFt")]  = 0 ;
    Material_GetProperty(mat)[pm("N_C3AH6")]  = 0 ;
    Material_GetProperty(mat)[pm("R_AFm")]  = 4.6e-4 ; /* 4.6e-4 (mol/L/s) Salgues 2013 */
    Material_GetProperty(mat)[pm("R_AFt")]  = 4.6e-4 ;
    Material_GetProperty(mat)[pm("R_C3AH6")] = 1.e-10 ;
    Material_GetProperty(mat)[pm("R_CSH2")]  = 1.e-10 ;
    Material_GetProperty(mat)[pm("AlphaCoef")] = 0 ;
    Material_GetProperty(mat)[pm("BetaCoef")]  = 0 ;
    Material_GetProperty(mat)[pm("r0")] = 16*nm ;
    
    Material_GetProperty(mat)[pm("B_i")] = -1 ;
    Material_GetProperty(mat)[pm("B_p")] = -1 ;

    damagemodel = -1 ;
    
    Material_ScanProperties(mat,datafile,pm) ;
    
    if(Material_GetProperty(mat)[pm("B_i")] < 0) {
      double c = Material_GetProperty(mat)[pm("A_i")] ;
      
      Material_GetProperty(mat)[pm("B_i")] = c ;
    }
    
    if(Material_GetProperty(mat)[pm("B_p")] < 0) {
      double c = Material_GetProperty(mat)[pm("A_p")] ;
      
      Material_GetProperty(mat)[pm("B_p")] = c ;
    }
  }

  {
    ComputePhysicoChemicalProperties() ;
  }
  
  /* Damage */
  {
    damage = Damage_Create() ;
      
    Material_AppendData(mat,1,damage,Damage_t,"Damage") ;
  }
  
  /* Elastic and plastic properties */
  {
    Elasticity_t* elasty = Damage_GetElasticity(damage) ;
    
    {
      /* Elasticity */
      {
        double young   = Material_GetPropertyValue(mat,"young") ;
        double poisson = Material_GetPropertyValue(mat,"poisson") ;
        
        Elasticity_SetToIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson) ;
      
        {
          double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
          Elasticity_ComputeStiffnessTensor(elasty,c) ;
        }
      }
    
      
      /* Mazars */
      if(damagemodel == 1) {
        double maxelastrain = Material_GetPropertyValue(mat,"max_elastic_strain") ;
        double A_c = Material_GetPropertyValue(mat,"A_c") ;
        double A_t = Material_GetPropertyValue(mat,"A_t") ;
        double B_c = Material_GetPropertyValue(mat,"B_c") ;
        double B_t = Material_GetPropertyValue(mat,"B_t") ;
        
        Damage_SetToMazars(damage) ;
        Damage_SetParameters(damage,maxelastrain,A_c,A_t,B_c,B_t) ;
      
      /* Marigo-Jirasek */
      } else if(damagemodel == 2) {
        double ft = Material_GetPropertyValue(mat,"uniaxial_tensile_strength") ;
        double Gf = Material_GetPropertyValue(mat,"fracture_energy") ;
        double w  = Material_GetPropertyValue(mat,"crack_band_width") ;
        
        Damage_SetToMarigoJirasek(damage) ;
        Damage_SetParameters(damage,ft,Gf,w) ;
        
      } else {
        Message_FatalError("Unknown model") ;
      }
      
    }

    #if 0
    {      
      Elasticity_PrintStiffnessTensor(elasty) ;
    }
    #endif
  }

  {
    if(!csd) csd = CementSolutionDiffusion_Create() ;
    if(!hcc) hcc = HardenedCementChemistry_Create() ;
    
    HardenedCementChemistry_SetRoomTemperature(hcc,TEMPERATURE) ;
    
    CementSolutionDiffusion_SetRoomTemperature(csd,TEMPERATURE) ;
  
    {
      Curves_t* curves = Material_GetCurves(mat) ;
      int i ;

      if((i = Curves_FindCurveIndex(curves,"S_r")) < 0) {
        arret("ReadMatProp: no cumulative pore volume fraction") ;
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
  }

  return(NbOfProp) ;
}



int PrintModelProp(Model_t* model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  int dim = Model_GetDimension(model) ;
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  
  printf("The equations are:\n") ;
  printf("\t- Conservation of mole content of Sulfur      (sulfur)\n") ;
  printf("\t- Conservation of charge ............         (charge)\n") ;
  printf("\t- Conservation of mole content of Calcium     (calcium)\n") ;
  printf("\t- Conservation of mole content of alkali      (potassium)\n") ;
  printf("\t- Conservation of mole content of Aluminium   (aluminium)\n") ;
#ifdef E_eneutral
  printf("\t- Electroneutrality                           (electroneutrality)\n") ;
#endif
#ifdef E_kinetics
  printf("\t- Crystal growth kinetics                     (kinetics)\n") ;
#endif
  printf("\t- Equilibrium equation                        (meca1,meca_2,meca_3)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
#if   defined (U_LogC_H2SO4)
  printf("\t- Sulfuric acid concentration     (logc_h2so4)\n") ;
#elif defined (U_C_H2SO4)
  printf("\t- Sulfuric acid concentration     (c_h2so4)\n") ;
#elif defined (U_LogC_SO4)
  printf("\t- Sulfate concentration           (logc_so4)\n") ;
#elif defined (U_C_SO4)
  printf("\t- Sulfate concentration           (c_so4)\n") ;
#endif
  printf("\t- Electric potential              (psi)\n") ;
  printf("\t- Zeta unknown for calcium        (z_ca)\n") ;
  printf("\t- Potassium concentration         (c_k)\n") ;
  printf("\t- Zeta unknown for aluminium      (z_al)\n") ;
#ifdef E_eneutral
  #if defined (U_LogC_OH)
  printf("\t- Hydroxide ion concentration     (logc_oh)\n") ;
  #elif defined (U_C_OH)
  printf("\t- Hydroxide ion concentration     (c_oh)\n") ;
  #endif
#endif
#ifdef E_kinetics
  printf("\t- Crystallization pressure        (p_c)\n") ;
#endif
  printf("\t- The displacement vector         (u_1,u_2,u_3)\n") ;

  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length : dm !\n") ;
  printf("\t time   : s !\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;
  fprintf(ficd,"N_CH   = 6.1      # CH mole content (moles/L)\n") ;
  fprintf(ficd,"N_K    = 0.4      # K mole content  (moles/L)\n") ;
  fprintf(ficd,"N_AH3  = 0.4      # Al mole content (moles/L)\n") ;
  fprintf(ficd,"N_AFm  = 0.1      # AFm mole content (moles/L)\n") ;
  fprintf(ficd,"N_AFt  = 0.4      # AFt mole content (moles/L)\n") ;
  fprintf(ficd,"Curves = file     # Pore volume fraction curve:  r  S_r\n") ;
  fprintf(ficd,"Curves = solid    # File name: S_CH  X_CSH  Z_CSH  S_SH\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    int    i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
#undef R
}



int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* vi0 = Element_GetImplicitTerm(el) ;
  double* ve0 = Element_GetExplicitTerm(el) ;
  double* v00 = Element_GetConstantTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  {
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      double* vi  = vi0 + p*NVI ;
      double* v0  = v00 + p*NV0 ;
      double zn_ca_s    = FEM_ComputeUnknown(fem,u,intfct,p,U_Calcium) ;
      double zn_si_s    = 1 ;
      double zn_al_s    = FEM_ComputeUnknown(fem,u,intfct,p,U_Aluminium) ;
  
      /* Solve cement chemistry */
      {
        double u_sulfur   = FEM_ComputeUnknown(fem,u,intfct,p,U_Sulfur) ;
        double logc_na    = -99 ;
        double logc_k     = FEM_ComputeUnknown(fem,u,intfct,p,U_Potassium) ;
        #ifdef E_eneutral
          double u_eneutral  = FEM_ComputeUnknown(fem,u,intfct,p,U_eneutral) ;
          #if defined (U_LogC_OH)
            double logc_oh = u_eneutral ;
          #elif defined (U_C_OH)
            double logc_oh = log10(u_eneutral) ;
          #endif
        #else
          double logc_oh = -7 ;
        #endif
  
        HardenedCementChemistry_SetInput(hcc,SI_CH,MIN(zn_ca_s,0)) ;
        HardenedCementChemistry_SetInput(hcc,SI_CSH,MIN(zn_si_s,0)) ;
        HardenedCementChemistry_SetInput(hcc,SI_AH3,MIN(zn_al_s,0)) ;
        #if defined (U_C_H2SO4) || defined (U_LogC_H2SO4)
          HardenedCementChemistry_SetInput(hcc,LogC_H2SO4,u_sulfur) ;
        #elif defined (U_C_SO4) || defined (U_LogC_SO4)
          HardenedCementChemistry_SetInput(hcc,LogC_SO4,u_sulfur) ;
        #endif
        HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
        HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
        HardenedCementChemistry_SetInput(hcc,LogC_OH,logc_oh) ;
    
        HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = 1.e-99 ;
        HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,Cl) = -99 ;
  
        HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O) ;
      
        HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      }
    
  
      {
        /* Liquid components 
         * ----------------- */
        double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
        //double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
        //double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
        //double c_oh   = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
    
        /* Solid contents 
         * -------------- */
        /* ... as components: CH, CSH2, CSH, AH3, AFm, AFt, C3AH6 */
        double n_ch       = CHSolidContent(zn_ca_s) ;
        double n_csh2     = n_csh2_0 ;
        double n_ah3      = AH3SolidContent(zn_al_s) ;
        double n_afm      = n_afm_0 ;
        double n_aft      = n_aft_0 ;
        double n_c3ah6    = n_c3ah6_0 ;
        double n_csh      = CSHSolidContent(zn_si_s) ;
        /* ... as volume */
        double v_csh      = MolarVolumeOfCSH(x_csh) ;
        double v_cem      = V_CH*n_ch + v_csh*n_csh + V_AH3*n_ah3 + V_AFm*n_afm + V_AFt*n_aft + V_C3AH6*n_c3ah6 + V_CSH2*n_csh2 ;

        /* Porosity */
        double phi_c = phi0 ;
        //double phi   = phi_c - V_CSH2*n_csh2 ;
        double phi   = phi_c ;
    

        /* Back up what is needed to compute components */
        N_CH    = n_ch ;
        N_CSH2  = n_csh2 ;
        N_AFm   = n_afm ;
        N_AFt   = n_aft ;
        N_C3AH6 = n_c3ah6 ;
        PHI     = phi ;
        PHI_C   = phi_c ;

        #ifndef E_eneutral
          C_OH    = c_oh ;
        #endif

        V_Cem0  = v_cem ;
      }
    
      {
        PoreRadius    = PoreEntryRadiusMax ;
      
        Beta_p        = 1 ;
        Straind       = strain0 ;
        VarPHI_C      = 0 ;
      }
        
      {
        int i ;
          
        for(i = 0 ; i < 9 ; i++) Stress[i]  = 0 ;
        
        Hardv  = hardv0 ;
        Damage = 0 ;
      }
    }
  }
  
  
  /* Compute here vi0 */
  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,vi0,0,0,p) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,p) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
      
      /* storage in vi */
      {
        double* vi  = vi0 + p*NVI ;
        int    i ;

        /* Molar contents */
        N_S  = x[I_N_S] ;
        N_Ca = x[I_N_Ca] ;
        N_Si = x[I_N_Si] ;
        N_K  = x[I_N_K] ;
        N_Cl = x[I_N_Cl] ;
        N_Al = x[I_N_Al] ;
    
        for(i = 0 ; i < 3 ; i++) W_S[i]  = x[I_W_S  + i] ;
        for(i = 0 ; i < 3 ; i++) W_Ca[i] = x[I_W_Ca + i] ;
        for(i = 0 ; i < 3 ; i++) W_K[i]  = x[I_W_K  + i] ;
        for(i = 0 ; i < 3 ; i++) W_Al[i] = x[I_W_Al + i] ;
        
        for(i = 0 ; i < 9 ; i++) Stress[i]  = x[I_Stress + i] ;
    
        /* Charge density */
        N_q  = x[I_N_q] ;

    
        N_CH    = x[I_N_CH] ;
        N_CSH2  = x[I_N_CSH2] ;
        N_AFm   = x[I_N_AFm] ;
        N_AFt   = x[I_N_AFt] ;
        N_C3AH6 = x[I_N_C3AH6] ;
        PHI     = x[I_PHI] ;
        PHI_C   = x[I_PHI_C] ;
      
        #ifndef E_eneutral
          C_OH    = x[I_C_OH] ;
        #endif
    
        PoreRadius    = x[I_PoreRadius] ;
      
        Beta_p        = x[I_Beta_p] ;
        Straind       = x[I_Straind] ;
        VarPHI_C      = x[I_VarPHI_C] ;
        GrowthRate    = x[I_GrowthRate] ;
        
          
        for(i = 0 ; i < 9 ; i++) Stress[i] = x[I_Stress + i] ;
          
        Hardv  = x[I_Hardv] ;
        Damage = x[I_Damage] ;
        Crit   = x[I_Crit] ;
      }
    }
  }
  
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the (current) explicit terms.
 *  IMPORTANT: if needed use only the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /*
    Transfer coefficient
  */
  
  //ComputeTransferCoefficients(el,u,f) ;
  

  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the (current) implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi0  = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  
  /* Compute here vi (with the help of vi_n if needed) */
  /* Loop on integration points */
  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vi_n,t,dt,p) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,p) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;
      
      /* storage in vi */
      {
        double* vi  = vi0 + p*NVI ;
        int    i ;

        /* Molar contents */
        N_S  = x[I_N_S] ;
        N_Ca = x[I_N_Ca] ;
        N_Si = x[I_N_Si] ;
        N_K  = x[I_N_K] ;
        N_Cl = x[I_N_Cl] ;
        N_Al = x[I_N_Al] ;
    
        for(i = 0 ; i < 3 ; i++) W_S[i]  = x[I_W_S  + i] ;
        for(i = 0 ; i < 3 ; i++) W_Ca[i] = x[I_W_Ca + i] ;
        for(i = 0 ; i < 3 ; i++) W_K[i]  = x[I_W_K  + i] ;
        for(i = 0 ; i < 3 ; i++) W_Al[i] = x[I_W_Al + i] ;
    
        /* Charge density */
        N_q  = x[I_N_q] ;

    
        N_CH    = x[I_N_CH] ;
        N_CSH2  = x[I_N_CSH2] ;
        N_AFm   = x[I_N_AFm] ;
        N_AFt   = x[I_N_AFt] ;
        N_C3AH6 = x[I_N_C3AH6] ;
        PHI     = x[I_PHI] ;
        PHI_C   = x[I_PHI_C] ;
      
        #ifndef E_eneutral
          C_OH    = x[I_C_OH] ;
        #endif
    
        PoreRadius    = x[I_PoreRadius] ;
      
        Beta_p        = x[I_Beta_p] ;
        Straind       = x[I_Straind] ;
        VarPHI_C      = x[I_VarPHI_C] ;
        GrowthRate    = x[I_GrowthRate] ;

        for(i = 0 ; i < 9 ; i++) Stress[i] = x[I_Stress + i] ;
          
        Hardv  = x[I_Hardv] ;
        Damage = x[I_Damage] ;
        Crit   = x[I_Crit] ;
      }
      
      #if 0
        printf("\n") ;
        printf("Damage:\n") ;
        printf("d: % e\n",x[I_Damage]) ;
        printf("Crystallization pressure:\n") ;
        printf("P_C: % e\n",x[I_P_C]) ;
        Math_PrintStressTensor(x+I_Strain) ;
        printf("Strain:\n") ;
        Math_PrintStressTensor(x+I_Strain) ;
        printf("Stress:\n") ;
        Math_PrintStressTensor(x+I_Stress) ;
      #endif


      {
        int test = 0 ;
        /*
        if(x[I_C_H2SO4] < 0.) test = -1 ;
        if(x[I_N_Ca_S]  < 0.) test = -1 ;
        if(x[I_N_Si_S]  < 0.) test = -1 ;
        if(x[I_N_Al_S]  < 0.) test = -1 ;
        */
        /*
        if(x[I_PHI]     < 0.) test = -1 ;
        */
      
        if(test < 0) {
          double c_h2so4 = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,H2SO4) ;
          double c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
          double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
          double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
          double s_ah3  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AH3) ;
          double s_afm  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFm) ;
          double s_aft  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFt) ;
          double s_c3ah6  = HardenedCementChemistry_GetSaturationIndexOf(hcc,C3AH6) ;


          double* h = IntFct_GetFunctionAtPoint(intfct,p) ;
          double* coor = Element_ComputeCoordinateVector(el,h) ;
          double xx = coor[0] ;
        
          printf("x         = %e\n",xx) ;
          printf("phi       = %e\n",x[I_PHI]) ;
          printf("phi_c     = %e\n",x[I_PHI_C]) ;
          printf("c_h2so4   = %e\n",c_h2so4) ;
          printf("n_ch      = %e\n",x[I_N_CH]) ;
          printf("n_csh2    = %e\n",x[I_N_CSH2]) ;
          printf("n_csh     = %e\n",x[I_N_CSH]) ;
          printf("n_ah3     = %e\n",x[I_N_AH3]) ;
          printf("n_afm     = %e\n",x[I_N_AFm]) ;
          printf("n_aft     = %e\n",x[I_N_AFt]) ;
          printf("n_c3ah6   = %e\n",x[I_N_C3AH6]) ;
          printf("s_ch      = %e\n",s_ch) ;
          printf("s_csh2    = %e\n",s_csh2) ;
          printf("s_ah3     = %e\n",s_ah3) ;
          printf("s_afm     = %e\n",s_afm) ;
          printf("s_aft     = %e\n",s_aft) ;
          printf("s_c3ah6   = %e\n",s_c3ah6) ;
          printf("zn_ca_s   = %e\n",x[I_U_Calcium]) ;
          printf("zn_al_s   = %e\n",x[I_U_Aluminium]) ;
          printf("c_oh      = %e\n",c_oh) ;
          return(-1) ;
        }
      }
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  
  
  /* Initialization */
  {
    int    i ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;
  }


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  
  /*
  ** Poromechanics matrix
  */
  #ifdef E_Mech
  {
    int ndif = NEQ - dim ;
    int n = 81 + ndif*9 + ndif*(9 + ndif) ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = TangentCoefficients(el,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,ndif,E_Mech) ;
    
    {
      int i ;
      
      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    }
  }
  #endif
  
  /*
  ** Conduction Matrix
  */
  #if 0
  {
    int ndif = NEQ - dim ;
    int n = 9*ndif*ndif ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTransferCoefficients(el,dt,c) ;
    double* kc_mass = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    double* kc_salt = FEM_ComputeConductionMatrix(fem,intfct,c+9*4,dec) ;
    double* kc_the  = FEM_ComputeConductionMatrix(fem,intfct,c+9*8,dec) ;
    int    i ;
  
    for(i = 0 ; i < nn ; i++) {
      int    j ;
      
      for(j = 0 ; j < nn ; j++) {
        K(E_Mass + i*NEQ,U_Mass + j*NEQ) += dt*kc_mass[i*nn + j] ;
        K(E_Salt + i*NEQ,U_Salt + j*NEQ) += dt*kc_salt[i*nn + j] ;
        K(E_The  + i*NEQ,U_The  + j*NEQ) += dt*kc_the[i*nn + j] ;
      }
    }
  }
  #endif
  


#if defined (U_C_H2SO4)
    {
      double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
      for(i = 0 ; i < 2*NEQ ; i++){
        K(i,E_Sulfur)     /= Ln10*C_H2SO4(0) ;
        K(i,E_Sulfur+NEQ) /= Ln10*C_H2SO4(1) ;
      }
    }
#elif defined (U_C_SO4)
    {
      double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
      for(i = 0 ; i < 2*NEQ ; i++){
        K(i,E_Sulfur)     /= Ln10*C_SO4(0) ;
        K(i,E_Sulfur+NEQ) /= Ln10*C_SO4(1) ;
      }
    }
#endif


#if defined (U_C_K)
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
  #endif
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* vi   = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = 0. ;
  }

  if(Element_IsSubmanifold(el)) return(0) ;

  
  /* Compute here the residu R(n,i) */
  
#ifdef E_Mech
  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  {
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,Stress,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  #if 0
  {
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Mech + dim - 1) -= -rbf[i] ;
    }
    
  }
  #endif
#endif
  
  
  /* 2. Conservation of sulfur */
  {
    double* ra = FEM_ComputeMassBalanceEquationResidu(fem,intfct,&N_S,&N_Sn,W_S,dt,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Sulfur) -= ra[i] ;
    }
  }
  
  
  /* 3. Conservation of calcium */
  {
    double* ra = FEM_ComputeMassBalanceEquationResidu(fem,intfct,&N_Ca,&N_Can,W_Ca,dt,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Calcium) -= ra[i] ;
    }
  }
  
  
  /* 4. Conservation of potassium */
  {
    double* ra = FEM_ComputeMassBalanceEquationResidu(fem,intfct,&N_K,&N_Kn,W_K,dt,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Potassium) -= ra[i] ;
    }
  }
  
  
  /* 5. Conservation of aluminium */
  {
    double* ra = FEM_ComputeMassBalanceEquationResidu(fem,intfct,&N_Al,&N_Aln,W_Al,dt,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_Aluminium) -= ra[i] ;
    }
  }
  
  
  /* 6. Conservation of charge: div(W_q) = 0 */
  {
    double* ra = FEM_ComputeFluxResidu(fem,intfct,W_q,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_charge) -= dt*ra[i] ;
    }
  }
  
  
#ifdef E_eneutral
  /* Electroneutrality */
  {
    double* ra = FEM_ComputeBodyForceResidu(fem,intfct,&N_q,NVI) ;
    
    {
      int i ;
    
      for(i = 0 ; i < nn ; i++) R(i,E_eneutral) -= ra[i] ;
    }
  }
#endif


#ifdef E_kinetics
  /* Kinetics */
  {
    double g1[IntFct_MaxNbOfIntPoints] ;
    int i ;
    
    for(i = 0 ; i < np ; i++ , vi += NVI, vi_n += NVI) {
      g1[i] = VarPHI_C - VarPHI_Cn - dt * GrowthRate ;
    }
    
    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
      
      for(i = 0 ; i < nn ; i++) R(i,E_kinetics) -= ra[i] ;
    }

  }
#endif
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/* Les valeurs exploitees (s) */
{
  int    NbOfOutputs = 64 ;
  double* ve0  = Element_GetExplicitTerm(el) ;
  double* vi0  = Element_GetCurrentImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  //int dim = Element_GetDimensionOfSpace(el) ;
  //FEM_t* fem = FEM_GetInstance(el) ;
  double zero = 0 ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Initialization */
  {
    int i ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }
  
  /*
    Input data
  */
  GetProperties(el) ;



  /* output quantities */
  {
    int    i ;
    /* Interpolation functions at s */
    //double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    //int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    int p = 0 ;
    /* Variables */
    double* x = ComputeVariables(el,u,u,vi0,t,0,p) ;
    
    /* Concentrations */
#define ptC(CPD)   &(HardenedCementChemistry_GetAqueousConcentrationOf(hcc,CPD))
#define ptE(CPD)   &(HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,CPD))
#define ptS(CPD)   &(HardenedCementChemistry_GetSaturationIndexOf(hcc,CPD))
#define ptPSI      &(HardenedCementChemistry_GetElectricPotential(hcc))
#define ptX_CSH    &(HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc))
#define ptZ_CSH    &(HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc))

    
    /* fluxes */
    #if 1
    double* vi   = vi0 + p*NVI ;
    double* w_s  = W_S ;
    double* w_ca = W_Ca ;
    double* w_al = W_Al ;
    double* w_si = W_Si ;
    double* w_q  = W_q ;
    /* stresses */
    double* sig = Stress ;
    double hardv = Hardv ;
    /* Damage */
    double d = Damage ;
    double crit = Crit ;
    #endif
    #if 0
    double w_s[3]  = {0,0,0} ;
    double w_ca[3] = {0,0,0} ;
    double w_al[3] = {0,0,0} ;
    double w_si[3] = {0,0,0} ;
    double w_q[3]  = {0,0,0} ;
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    /* Damage */
    double d = 0 ;
    double crit = 0 ;
    
    /* Averaging */
    for(i = 0 ; i < np ; i++) {
      double* vi  = vi0 + i*NVI ;
      double* ve  = ve0 + i*NVE ;
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_s[j]  += W_S[j] /np ;
      for(j = 0 ; j < 3 ; j++) w_ca[j] += W_Ca[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_al[j] += W_Al[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_si[j] += W_Si[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_q[j]  += W_q[j] /np ;
      for(j = 0 ; j < 9 ; j++) sig[j]  += Stress[j] /np ;
      
      hardv += Hardv/np ;
      crit  += Crit/np ;
      d     += Damage/np ;
    }
    #endif
    
    i = 0 ;
    {
      double ph        = - log10(*(ptC(H ))) ;
      
      Result_Store(r + i++,&ph,"ph",1) ;
    }
/*
    Result_Store(r + i++,ptC(OH),"c_oh",1) ;
    Result_Store(r + i++,ptC(H ),"c_h",1) ;
*/
    Result_Store(r + i++,ptC(Ca  ),"c_ca",1) ;
    Result_Store(r + i++,ptC(CaOH),"c_caoh",1) ;
    
    Result_Store(r + i++,ptC(H2SiO4),"c_h2sio4",1) ;
    Result_Store(r + i++,ptC(H3SiO4),"c_h3sio4",1) ;
    Result_Store(r + i++,ptC(H4SiO4),"c_h4sio4",1) ;
    
    Result_Store(r + i++,ptC(CaH2SiO4),"c_cah2sio4",1) ;
    Result_Store(r + i++,ptC(CaH3SiO4),"c_cah3sio4",1) ;
    
    Result_Store(r + i++,ptC(H2SO4),"c_h2so4",1) ;
    Result_Store(r + i++,ptC(HSO4 ),"c_hso4",1) ;
    Result_Store(r + i++,ptC(SO4  ),"c_so4",1) ;
    
    Result_Store(r + i++,ptC(CaSO4  ),"c_caso4aq",1) ;
    Result_Store(r + i++,ptC(CaHSO4 ),"c_cahso4",1) ;
    
    Result_Store(r + i++,ptC(Al    ),"c_al",1) ;
    Result_Store(r + i++,ptC(AlO4H4),"c_alo4h4",1) ;

    Result_Store(r + i++,ptC(K  ),"c_k",1) ;
    Result_Store(r + i++,ptC(KOH),"c_koh",1) ;
    
    Result_Store(r + i++,ptE(Ca),"C_Ca",1) ;
    Result_Store(r + i++,ptE(Si),"C_Si",1) ;
    Result_Store(r + i++,ptE(S ),"C_S",1) ;
    Result_Store(r + i++,ptE(Al),"C_Al",1) ;
    Result_Store(r + i++,ptE(K ),"C_K",1) ;
    
    Result_Store(r + i++,ptS(CH   ),"s_ch",1) ;
    Result_Store(r + i++,ptS(CSH2 ),"s_csh2",1) ;
    Result_Store(r + i++,ptS(AH3  ),"s_ah3",1) ;
    Result_Store(r + i++,ptS(AFm  ),"s_afm",1) ;
    Result_Store(r + i++,ptS(AFt  ),"s_aft",1) ;
    Result_Store(r + i++,ptS(C3AH6),"s_c3ah6",1) ;
    
    Result_Store(r + i++,(x + I_N_CH   ),"n_ch",1) ;
    Result_Store(r + i++,(x + I_N_CSH2 ),"n_csh2",1) ;
    Result_Store(r + i++,(x + I_N_CSH  ),"n_csh",1) ;
    Result_Store(r + i++,(x + I_N_AH3  ),"n_ah3",1) ;
    Result_Store(r + i++,(x + I_N_AFm  ),"n_afm",1) ;
    Result_Store(r + i++,(x + I_N_AFt  ),"n_aft",1) ;
    Result_Store(r + i++,(x + I_N_C3AH6),"n_c3ah6",1) ;
    {
      double n_csh = x[I_N_CSH] ;
      double n_so4ads = n_csh*AdsorbedSulfatePerUnitMoleOfCSH(*ptC(SO4),*ptC(OH)) ;
      
      Result_Store(r + i++,&n_so4ads,"n_so4^ads",1) ;
    }
    
    Result_Store(r + i++,(x + I_PHI),"porosite",1) ;
    Result_Store(r + i++,ptPSI,"potentiel_electrique",1) ;
    
    Result_Store(r + i++,(x + I_N_q),"charge",1) ;
    
    Result_Store(r + i++,(x + I_V_CSH),"V_CSH",1) ;
    Result_Store(r + i++,ptX_CSH,"C/S",1) ;
    
    Result_Store(r + i++,w_si,"W_Si",3) ;
    Result_Store(r + i++,w_ca,"W_Ca",3) ;
    Result_Store(r + i++,w_s ,"W_S",3) ;
    Result_Store(r + i++,w_al,"W_Al",3) ;
    Result_Store(r + i++,w_q ,"W_q",3) ;
    
    Result_Store(r + i++,&zero,"P_CSH2",1) ;
    Result_Store(r + i++,&zero,"Damage",1) ;
    
    Result_Store(r + i++,x + I_N_Ca,"N_Ca",1) ;
    Result_Store(r + i++,x + I_N_Si,"N_Si",1) ;
    Result_Store(r + i++,x + I_N_S ,"N_S",1) ;
    Result_Store(r + i++,x + I_N_Al,"N_Al",1) ;
    Result_Store(r + i++,x + I_N_K ,"N_K",1) ;
    Result_Store(r + i++,x + I_N_Cl,"N_Cl",1) ;
    
    Result_Store(r + i++,(x + I_S_C    ),"Saturation degree of crystal",1) ;
    
    {
      double radius = x[I_PoreRadius] ;
      double beta = InterfaceEquilibriumSaturationIndex(radius) ;
      
      Result_Store(r + i++,&beta,"Interface equilibrium saturation index of AFt",1) ;
    }
    
    {
      Result_Store(r + i++,x + I_Beta_p,"Pore wall equilibrium saturation index of AFt",1) ;
    }
    
    Result_Store(r + i++,(x + I_P_C    ),"Crystallization pressure",1) ;
    
    Result_Store(r + i++,x + I_Strain,"Strain tensor",9) ;
    Result_Store(r + i++,sig,"Stress tensor",9) ;
    Result_Store(r + i++,x + I_Dis,"Displacement vector",3) ;
    
    Result_Store(r + i++,&d    ,"Damage",1) ;
    Result_Store(r + i++,&hardv,"Hardening variable",1) ;
    Result_Store(r + i++,&crit ,"Yield function",1) ;
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}


#if 0
void ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Transfer coefficients  */
{
  double* ve = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;

  /* Initialization */
  for(i = 0 ; i < NVE ; i++) ve[i] = 0. ;


  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;
    
    /* Liquid tortuosity */
    {
      double phi    = x[I_PHI] ;
      double iff    = LiquidTortuosity(phi) ;
        
      TORTUOSITY[i] = iff ;
    }
    
    /* Concentrations */
    {
      double* c = HardenedCementChemistry_GetAqueousConcentration(hcc) ;
      int n = HardenedCementChemistry_NbOfConcentrations ;
      int j ;
    
      for(j = 0 ; j < n ; j++) {
        CONCENTRATION(i)[j] = c[j] ;
      }
    }
  }
}
#endif



int TangentCoefficients(Element_t* el,double t,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define T2(a,i,j)      ((a)[(i)*3+(j)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  double* vi0  = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double dui[NEQ] ;
  int ndif = NEQ - dim ;
  int  dec = 81 + ndif*9 + ndif*(9 + ndif) ;


  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  {
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }
  

  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < np ; p++) {
      double* vi  = vi0 + p*NVI ;
      double* c0 = c + p*dec ;
      /* Variables */
      double* x   = ComputeVariables(el,u,u_n,vi_n,t,dt,p) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,p) ;
      
      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;


      #if defined (U_C_H2SO4) || defined (U_C_SO4)
        dui[E_Sulfur] =  1.e-2*ObVal_GetValue(obval + E_Sulfur)/(Ln10*Un_Sulfur(i)) ;
      #endif
    
      #if defined (U_C_K)
        dui[E_Potassium] =  1.e-2*ObVal_GetValue(obval + E_Potassium)/(Ln10*Un_Potassium(i)) ;
      #endif
    
      #ifdef E_eneutral
        #if defined (U_C_OH)
        dui[E_eneutral] =  1.e-2*ObVal_GetValue(obval + E_eneutral)/(Ln10*Un_eneutral(i)) ;
        #endif
      #endif
      
      
      /* The derivative of equations w.r.t unknowns */
      

      /* Derivatives w.r.t strains */
      {
        /* Tangent stiffness matrix */
        {
          double* c1 = c0 ;
          
          CopyStiffnessTensor(c1) ;
      
          /* Tangent stiffness matrix */
          {
            /* Criterion */
            double crit = Crit ;
        
            if(crit >= 0.) {
              /* Strains */
              double* eps = x + I_Strain ;
              double crit1 = ComputeFunctionGradients(eps,&Damage,&Hardv) ;
              double fcg = UpdateTangentStiffnessTensor(c1) ;
          
              if(fcg < 0) return(-1) ;
            }
          }
          //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
        }
        
        /* Coupling matrices */
        {
          double deps = 1.e-6 ;
          double* dx = ComputeVariableDerivatives(el,t,dt,x,deps,I_Strain) ;
          double* c00 = c0 + 81 + ndif*9 ;
          
          /* assuming to be the same for the derivatives wrt I_Strain+(0/4/8)
           * and zero for the derivatives wrt others */
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      
      /* Derivatives w.r.t U_Sulfur */
      {
        double  dp = dui[E_Sulfur] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_Sulfur) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_Sulfur*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_Sulfur ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      
      /* Derivatives w.r.t U_Calcium */
      {
        double  dp = dui[E_Calcium] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_Calcium) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_Calcium*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_Calcium ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      
      /* Derivatives w.r.t U_Potassium */
      {
        double  dp = dui[E_Potassium] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_Potassium) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_Potassium*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_Potassium ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      
      /* Derivatives w.r.t U_Aluminium */
      {
        double  dp = dui[E_Aluminium] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_Aluminium) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_Aluminium*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
          
          #if 0
          printf("(point %d):\n",p) ;
          printf("Stress:\n") ;
          Math_PrintStressTensor(x+I_Stress) ;
          printf("Derivative of stress wrt Aluminium:\n") ;
          Math_PrintStressTensor(c1) ;
          #endif
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_Aluminium ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      
      /* Derivatives w.r.t U_eneutral */
      #ifdef E_eneutral
      {
        double  dp = dui[E_eneutral] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_eneutral) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_eneutral*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_eneutral ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      #endif
      
      /* Derivatives w.r.t U_kinetics */
      #ifdef E_kinetics
      {
        double  dp = dui[E_kinetics] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_kinetics) ;
        
        /* Mechanical coupling terms */
        {
          double* dsigdp = dx + I_Stress ;
          double* c1 = c0 + 81 + E_kinetics*9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + ndif*9 + 9 + E_kinetics ;
          
          {
            double* c1 = c00 + E_Sulfur*(9 + ndif) ;
        
            c1[0] = dx[I_N_S] ;
          }
          {
            double* c1 = c00 + E_Calcium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Ca] ;
          }
          {
            double* c1 = c00 + E_Potassium*(9 + ndif) ;
        
            c1[0] = dx[I_N_K] ;
          }
          {
            double* c1 = c00 + E_Aluminium*(9 + ndif) ;
        
            c1[0] = dx[I_N_Al] ;
          }
          #ifdef E_eneutral
          {
            double* c1 = c00 + E_eneutral*(9 + ndif) ;
        
            c1[0] = dx[I_N_q] ;
          }
          #endif
          #ifdef E_kinetics
          {
            double* c1 = c00 + E_kinetics*(9 + ndif) ;

            c1[0] = dx[I_VarPHI_C] - dt * dx[I_GrowthRate] ;
          }
          #endif
        }
      }
      #endif
    }
  }

  return(dec) ;
#undef C1
#undef B1
#undef T2
#undef T4
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
/** This locally defined function compute the intern variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables. */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ; 
  double* x   = Variable ;
  double* x_n = Variable_n ;
  
  /* Load the primary variables in x */
  {
    int    i ;
    
#ifdef E_Mech
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_Dis + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_Dis + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_Strain + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
#else
    {
      for(i = 0 ; i < 9 ; i++) {
        x[I_Strain + i] = 0 ;
      }
    }
#endif

    /* Other primary unknowns and their gradients */
    /* Sulfur */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_Sulfur) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Sulfur) ;
    
      x[I_U_Sulfur    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_Sulfur + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    /* Calcium */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_Calcium) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Calcium) ;
    
      x[I_U_Calcium    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_Calcium + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    /* Potassium */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_Potassium) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Potassium) ;
    
      x[I_U_Potassium    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_Potassium + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    /* charge */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_charge) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_charge) ;
    
      x[I_U_charge    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_charge + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    /* Aluminium */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_Aluminium) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_Aluminium) ;
    
      x[I_U_Aluminium    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_Aluminium + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
#if defined (E_eneutral)
    /* eneutral */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_eneutral) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_eneutral) ;
    
      x[I_U_eneutral    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_eneutral + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
#endif
#if defined (E_kinetics)
    /* eneutral */
    {
      double  v   = FEM_ComputeUnknown(fem,u,intfct,p,U_kinetics) ;
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_kinetics) ;
    
      x[I_U_kinetics    ] = v ;
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_kinetics + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
#endif
  }
  
  
  /* Variables at the previous time */
  {
    double* vi_n  = f_n + p*NVI ;

#ifdef E_Mech
    /* Strains, stresses at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_Mech) ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_Strain + i] = eps_n[i] ;
        x_n[I_Stress + i] = Stress_n[i] ;
      }
      
      x_n[I_Hardv]  = Hardv_n ;
      x_n[I_Damage] = Damage_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
#else
    {
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_Strain + i] = 0 ;
        x_n[I_Stress + i] = 0 ;
      }
      
      x_n[I_Hardv]  = 0 ;
      x_n[I_Damage] = 0 ;
    }
#endif
    
    /* Other variables at previous time step */
    {
      x_n[I_N_CH  ]  = N_CHn ;
      x_n[I_N_CSH2]  = N_CSH2n ;
      x_n[I_N_AFm ]  = N_AFmn ;
      x_n[I_N_AFt ]  = N_AFtn ;
      x_n[I_N_C3AH6] = N_C3AH6n ;
      x_n[I_PHI   ]  = PHIn ;
      x_n[I_PHI_C ]  = PHI_Cn ;
      #ifndef E_eneutral
        x_n[I_C_OH  ]  = C_OHn ;
      #endif
      x_n[I_PoreRadius]  = PoreRadiusn ;
      x_n[I_Beta_p]  = Beta_pn ;
      x_n[I_Straind] = Straind_n;
      x_n[I_VarPHI_C] = VarPHI_Cn ;
    }
  }

  /* Constant terms */
  {
    {
      double* v00 = Element_GetConstantTerm(el) ;
      double* v0  = v00 + p*NV0 ;
    
      x_n[I_V_Cem0 ]  = V_Cem0 ;
    }
  }
  
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  /* Primary variables */
  double zn_si_s    = 1 ;
  double zn_ca_s    = x[I_U_Calcium] ;
  double zn_al_s    = x[I_U_Aluminium] ;
  double psi        = x[I_U_charge] ;
  
  /* Solve cement chemistry */
  {
    #if defined (U_C_H2SO4) || defined (U_LogC_H2SO4)
      double logc_h2so4 = x[I_U_Sulfur] ;
    #elif defined (U_C_SO4) || defined (U_LogC_SO4)
      double logc_so4 = x[I_U_Sulfur] ;
    #endif
    double logc_na    = -99 ;
    double logc_k     = x[I_U_Potassium] ;
    #if defined (E_eneutral)
      #if defined (U_LogC_OH) && !defined (U_C_OH)
        double logc_oh    = x[I_U_eneutral] ;
      #elif defined (U_C_OH) && !defined (U_LogC_OH)
        double c_oh       = x[I_U_eneutral] ;
        double logc_oh    = log10(c_oh) ;
      #endif
    #else
      double logc_oh    = log10(x_n[I_C_OH]) ;
    #endif
  
    HardenedCementChemistry_SetInput(hcc,SI_CH,MIN(zn_ca_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,SI_CSH,MIN(zn_si_s,0)) ;
    HardenedCementChemistry_SetInput(hcc,SI_AH3,MIN(zn_al_s,0)) ;
    #if defined (U_C_H2SO4) || defined (U_LogC_H2SO4)
      HardenedCementChemistry_SetInput(hcc,LogC_H2SO4,logc_h2so4) ;
    #elif defined (U_C_SO4) || defined (U_LogC_SO4)
      HardenedCementChemistry_SetInput(hcc,LogC_SO4,logc_so4) ;
    #endif
    HardenedCementChemistry_SetInput(hcc,LogC_Na,logc_na) ;
    HardenedCementChemistry_SetInput(hcc,LogC_K,logc_k) ;
    HardenedCementChemistry_SetInput(hcc,LogC_OH,logc_oh) ;
    HardenedCementChemistry_GetElectricPotential(hcc) = psi ;
    
    HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = 1.e-99 ;
    HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,Cl) = -99 ;
  
    HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O) ;

    #ifndef E_eneutral
      #if (ELECTRONEUTRALITY == IMPLICIT)
        HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      #endif
    #endif
  }
  

  
  /* Backup */
  double c_q_l  = HardenedCementChemistry_GetLiquidChargeDensity(hcc) ;
  
  //double I = HardenedCementChemistry_GetIonicStrength(hcc) ;
  
  //double rho_l  = HardenedCementChemistry_GetLiquidMassDensity(hcc) ;
  
  double c_ca_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Ca) ;
  double c_si_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Si) ;
  double c_k_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,K) ;
  double c_s_l  = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,S) ;
  double c_al_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Al) ;
  double c_cl_l = HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,Cl) ;
  
  //double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
  //double s_sh   = HardenedCementChemistry_GetSaturationIndexOf(hcc,SH) ;
  double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
  //double s_ah3  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AH3) ;
  double s_afm  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFm) ;
  double s_aft  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFt) ;
  double s_c3ah6 = HardenedCementChemistry_GetSaturationIndexOf(hcc,C3AH6) ;
  
  double c_so4  = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,SO4) ;
  //double c_oh   = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
       
    
  /* The crystal saturation index */
  //double beta = s_csh2 ;
  double beta = s_aft ;
  //double beta = Math_Max(s_aft,s_csh2) ;
  
  /* Compute the saturation degree of the crystal phase as a function of beta */
  double r_n = x_n[I_PoreRadius] ;
  double r   = Radius(r_n,beta,dt,el) ;
  double s_l = LiquidSaturationDegree(r) ;
  double s_c = 1 - s_l ;

  /* Compute the saturation index at the pore wall, beta_p */
  double beta_pn   = x_n[I_Beta_p] ;
  double varphi_cn = x_n[I_VarPHI_C] ;
  double* strain_n = x_n + I_Strain ;
  double* strain   = x   + I_Strain ;
  double strainv_n = strain_n[0] + strain_n[4] + strain_n[8] ;
  double strainv   = strain[0] + strain[4] + strain[8] ;
  double dstrainv  = strainv - strainv_n ;
  //double straind_n = x_n[I_Straind] ;
  
  /* Kinetic law */
  /* Crystallization pressure */
  #ifdef E_kinetics
    double p_c = x[I_U_kinetics] ;
    double beta_p = InverseOfCrystallizationPressure(p_c) ;
  #else
    double beta_p = PoreWallEquilibriumSaturationIndex(beta_pn,beta,dstrainv,dt) ;
  //double beta_p = PoreWallEquilibriumSaturationIndex_old(beta_pn,varphi_cn,strainv_n,straind_n,beta,s_c,dt) ;
    double p_c    = CrystallizationPressure(beta_p) ;
  #endif
  
  /* Compute the crystal pore deformation */
  //double phicrate  = PoreCrystalGrowthRate(s_c,beta,beta_p) ;
  //double varphi_c  = varphi_cn + dt*phicrate ;
  double varphi_c  = varphi_cn + s_c*dstrainv ;
  
  
  /* Solid contents
   * -------------- */
  /* ... as components: CH, CSH2, CSH, AH3, AFm, AFt, C3AH6 */
  
  /* The crystal responsible for strains is either AFt or CSH2 (gypsum) */
  double v_crystal  =  phi0*s_c + varphi_c ;
  double n_crystal  =  v_crystal/V_Crystal ;

  double n_aftn     = x_n[I_N_AFt] ;
  //double n_aft      = AFtSolidContent(n_aftn,s_aft,dt) ;
  double n_aft      = v_crystal/V_AFt ;  
  //double n_aft      = (s_aft > s_csh2) ? n_crystal : 0 ;
  
  double n_csh2n    = x_n[I_N_CSH2] ;
  double n_csh2     = CSH2SolidContent(n_csh2n,s_csh2,dt) ;
  //double n_csh2     = v_crystal/V_CSH2 ;
  //double n_csh2     = (s_aft > s_csh2) ? 0 : n_crystal ;
  
  //double n_chn      = x[I_N_CHn] ;
  double n_ch       = CHSolidContent(zn_ca_s) ;
  double n_ah3      = AH3SolidContent(zn_al_s) ;
  double n_afmn     = x_n[I_N_AFm] ;
  double n_afm      = AFmSolidContent(n_afmn,s_afm,dt) ;
  double n_c3ah6n   = x_n[I_N_C3AH6] ;
  double n_c3ah6    = C3AH6SolidContent(n_c3ah6n,s_c3ah6,dt) ;
  double n_csh      = CSHSolidContent(zn_si_s) ;
  double n_so4ads   = n_csh * AdsorbedSulfatePerUnitMoleOfCSH(c_so4,c_oh) ;
  
  /* ... as elements: S, Ca, Si, Al */
  //double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
  double n_si_s     = n_csh ;
  double n_ca_s     = n_ch + n_csh2 + x_csh*n_csh + 4*n_afm + 6*n_aft + 3*n_c3ah6 ;
  double n_s_s      = n_csh2 + n_afm + 3*n_aft + n_so4ads ;
  double n_al_s     = 2*(n_ah3 + n_afm + n_aft + n_c3ah6) ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(x_csh) ;
  //double v_gyp      = V_Gyp*n_csh2 ;
  //double v_csh2     = V_CSH2*n_csh2 ;
  double v_cem      = V_CH*n_ch + v_csh*n_csh + V_AH3*n_ah3 + V_AFm*n_afm + V_AFt*n_aft + V_C3AH6*n_c3ah6 + V_CSH2*n_csh2 ;


  /* Porosities in unconfined conditions (no pressure) */
  double v_cem0     = x[I_V_Cem0] ;
  /* ... of concrete */
  double phi_con    = phi0 + v_cem0 - v_cem ;
  /* ... of gypsum */
  //double phi_gyp    = PHI_Gyp ;


  /* Total porosity */
  double varphi     = strainv ;
  double phi_c      = phi_con + varphi ;
  //double phi_t      = phi_con - v_csh2 ;
  double phi_t      = MAX(phi_c,phimin) ;
  

#if (U_PHI == IMPLICIT)
  double phi_l        = phi_t ;
#else
  double phi_l        = x_n[I_PHI] ;
#endif
    
    
  /* Liquid contents 
   * --------------- */
  /* ... as elements: S, Ca, Si, K, Cl, Al */
  double n_s_l  = phi_l*c_s_l ;
  double n_ca_l = phi_l*c_ca_l ;
  double n_si_l = phi_l*c_si_l ;
  double n_al_l = phi_l*c_al_l ;
  double n_k_l  = phi_l*c_k_l ;
  double n_cl_l = phi_l*c_cl_l ;
  double n_q_l  = phi_l*c_q_l ;


#ifndef E_eneutral
  #if (ELECTRONEUTRALITY == EXPLICIT)
  c_q_l = HardenedCementChemistry_SolveExplicitElectroneutrality(hcc) ;
  n_q_l = phi_l*c_q_l ;
  #endif
#endif
  

  /* Back up */
      
  /* Backup crystallization pressure */
  x[I_P_C] = p_c ;
  x[I_Beta_p] = beta_p ;

#ifdef E_Mech
  /* Stresses, damage and hardening variables */
  {
    double* eps   = x + I_Strain ;
    double* sig   = x + I_Stress ;
    
    {
      int    i ;
    
      /* Projection */
      {
        double hardv = x_n[I_Hardv] ;
        double d     = x_n[I_Damage] ;
        double crit  = ReturnMapping(eps,&d,&hardv) ;
        
        x[I_Crit]   = crit ;
        x[I_Hardv]  = hardv ;
        x[I_Damage] = d ;
      }
    
      /* Backup stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = 0 ;
      
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        #define C(i,j)  (cijkl[(i)*9+(j)])
        for(j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*eps[j] ;
        }
        #undef C
      }
        
      {
        double  d  = x[I_Damage] ;
        double bd  = Biot + d * (1 - Biot) ;
        
        sig[0] += - bd * s_c * p_c ;
        sig[4] += - bd * s_c * p_c ;
        sig[8] += - bd * s_c * p_c ;
      }
      
      //Damage_PrintStiffnessTensor(damage) ;
    }
  }
#endif


  /* Solid components */
  x[I_N_CH      ] = n_ch ;
  x[I_N_CSH2    ] = n_csh2 ;
  x[I_N_AH3     ] = n_ah3 ;
  x[I_N_AFm     ] = n_afm ;
  x[I_N_AFt     ] = n_aft ;
  x[I_N_C3AH6   ] = n_c3ah6 ;
  x[I_N_CSH     ] = n_csh ;
  
  x[I_V_CSH     ] = v_csh ;
  
  x[I_V_Cem     ] = v_cem ;
  
  
  /* Porosities */
  x[I_PHI       ] = phi_t ;
  x[I_PHI_C     ] = phi_c ;
  
  /* Saturation degree of crystal */
  x[I_S_C       ] = s_c ;
  
  /* Radius */
  x[I_PoreRadius] = r ;
  
  /* Strains */
  //x[I_Straind   ] = DamageStrain(strainv,straind_n) ;
  x[I_VarPHI_C  ] = varphi_c ;
  
  
  
  /* Element contents */
  x[I_N_S       ] = n_s_l  + n_s_s ;
  x[I_N_Ca      ] = n_ca_l + n_ca_s ;
  x[I_N_Si      ] = n_si_l + n_si_s ;
  x[I_N_K       ] = n_k_l  ;
  x[I_N_Cl      ] = n_cl_l  ;
  x[I_N_Al      ] = n_al_l + n_al_s ;

  /* Charge density */
  x[I_N_q       ] = n_q_l ;
  
  /* Crystal growth */
  x[I_GrowthRate] = PoreCrystalGrowthRate(s_c,beta,beta_p) ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dui,int i)
{
  double* dx  = dVariable ;
  double* x_n = Variable_n ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dui ;
  
  ComputeSecondaryVariables(el,t,dt,x_n,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dui ;
  }

  return(dx) ;
}



#if 0
double* ComputeVariableFluxes(Element_t* el,int i,int j)
{
  double* grdij = dVariable ;

  /* Gradients */
  {
    int nn = Element_GetNbOfNodes(el) ;
    FEM_t* fem   = FEM_GetInstance(el) ;

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
#endif


#if 0
double* ComputeFluxes(Element_t* el,double* grdij,int i,int j)
{
  double* ve = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes[i] ;
  
  /* Diffusion in the cement solution */
  {
    /* Gradients */
    {
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
    }
    /* Fluxes */
    {
      
      CementSolutionDiffusion_ComputeFluxes(csd) ;
      
      w[I_W_Ca]  = CementSolutionDiffusion_GetElementFluxOf(csd,Ca) ;
      w[I_W_Si]  = CementSolutionDiffusion_GetElementFluxOf(csd,Si) ;
      w[I_W_S ]  = CementSolutionDiffusion_GetElementFluxOf(csd,S) ;
      w[I_W_K ]  = CementSolutionDiffusion_GetElementFluxOf(csd,K) ;
      w[I_W_Cl]  = CementSolutionDiffusion_GetElementFluxOf(csd,Cl) ;
      w[I_W_Al]  = CementSolutionDiffusion_GetElementFluxOf(csd,Al) ;
      w[I_W_q ]  = CementSolutionDiffusion_GetIonCurrent(csd) ;
    }
  }
    
  return(w) ;
}
#endif




double TortuosityToLiquid_OhJang(double phi)
/* Ref:
 * Byung Hwan Oh, Seung Yup Jang, 
 * Prediction of diffusivity of concrete based on simple analytic equations, 
 * Cement and Concrete Research 34 (2004) 463 - 480.
 * tau = (m_p + sqrt(m_p**2 + phi_c/(1 - phi_c) * (Ds/D0)**(1/n)))**n
 * m_p = 0.5 * ((phi_cap - phi_c) + (Ds/D0)**(1/n) * (1 - phi_c - phi_cap)) / (1 - phi_c)
 */
{
  double phi_cap = 0.5 * phi  ;
  double phi_c   = 0.17 ;         /* Percolation capilar porosity */
  double n       = 2.7 ;          /* OPC n  = 2.7  --------  Fly ash n  = 4.5 */
  double ds      = 1.e-4 ;        /* OPC ds = 1e-4 --------  Fly ash ds = 5e-5 */
  double dsn     = pow(ds,1/n) ;
  double m_phi   = 0.5 * ((phi_cap - phi_c) + dsn * (1 - phi_c - phi_cap)) / (1 - phi_c) ;
  double tausat  = pow(m_phi + sqrt(m_phi*m_phi + dsn * phi_c/(1 - phi_c)),n) ;
  
  double tau =  tausat ;
    
  return tau ;
}




double TortuosityToLiquid_BazantNajjar(double phi)
/* Ref:
 * Z. P. BAZANT, L.J. NAJJAR,
 * Nonlinear water diffusion in nonsaturated concrete,
 * Materiaux et constructions, 5(25), 1972.
 */
{
  double iff = 2.9e-4 * exp(9.95 * phi) ;
  double tausat = (iff < 1) ? iff : 1 ;
    
  return(tausat) ;
}



double Radius(double r_n,double beta,double dt,Element_t* el)
{
  double r_max  = PoreEntryRadiusMax ;
  double r = r_n ;
  
  //if(beta < 1) return(r_max) ;
  
  {
    double s_ln     = LiquidSaturationDegree(r_n) ;
    //double beta_i_in  = InterfaceEquilibriumSaturationIndex(r_n) ;
    double beta_i_min = InterfaceEquilibriumSaturationIndex(r_max) ;
    //double beta_i_inf = (beta > beta_i_min) ? beta : beta_i_min ;
    double r_inf = (beta > beta_i_min) ? InverseOfInterfaceEquilibriumSaturationIndex(beta) : r_max ;
    //double s_linf  = LiquidSaturationDegree(r_inf) ;
    int iterations = 40 ;
    double tol = 1.e-6 ;
    int i ;
    
    if(r_n == r_inf) return(r_n) ;
    
    for(i = 0 ; i < iterations ; i++) {
      double beta_i  =  InterfaceEquilibriumSaturationIndex(r) ;
      double dbeta_i = dInterfaceEquilibriumSaturationIndex(r) ;
      //double dbeta_i = (beta_i_inf - beta_i_in)/(r_inf - r_n) ;
      
      double s_l   =  LiquidSaturationDegree(r) ;
      double ds_l  = dLiquidSaturationDegree(r) ;
      //double ds_l  = (s_linf - s_ln)/(r_inf - r_n) ;
      
      /* Kinetic law */
      /* Modified 03/06/2017 */
      double  scrate =  InterfaceCrystalGrowthRate(beta,beta_i) ;
      double dscrate = dInterfaceCrystalGrowthRate(beta,beta_i)*dbeta_i ;
      double  eq   =  s_l - s_ln + dt*scrate ;
      double deq   = ds_l        + dt*dscrate ;
      
      double dr    = (fabs(deq) > 0) ? - eq/deq : 0. ;
      
      /* The solution r should be in the range between r_n and r_inf.
       * So let us assume that, at a given iteration, an estimation r
       * has been found between r_n and r_inf. Then we look for a new 
       * estimation r + dr in the range between r_0 and r_1 by using
       * an under-relaxation technique so that dr should be given by 
       *         r + dr = a*r_1 + (1 - a)*r_0
       * with a being in the range between 0 and 1.
       * The under-relaxation technique is such that r_0 should be in 
       * the range between r_n and r and r_1 should be in the range 
       * between r_inf and r i.e
       *         r_0 = b0*r_n   + (1 - b0)*r
       *         r_1 = b1*r_inf + (1 - b1)*r
       * with b0 and b1 being in between 0 and 1. So we get
       *         dr = a*b1*(r_inf - r) + (1 - a)*b0*(r_n - r)
       * If b0=b1=b then
       *         dr = (a*(r_inf - r) + (1 - a)*(r_n - r))*b
       * The bigger b the larger the range where r is looked for, without
       * exceding the initial range defined by r_n and r_inf.
       * The value b=0.5 corresponds to half the initial range.
       */
      {
        /* The predicted a is computed from the predicted dr */
        //double b = 0.5 ;
        double b = 0.5 ;
        double a = (dr/b - r_n + r)/(r_inf - r_n) ;
       
        /* If the predicted a is < 0 then the used value is set to 0 */
        if(a < 0) a = 0 ;
      
        /* if the predicted a is > 1 then the used value is set to 1 */
        if(a > 1) a = 1 ;
        
        {
          dr = b*(a*(r_inf - r) + (1 - a)*(r_n - r)) ;
        }
        
        /*
        printf("\n") ;
        printf("a      = %e\n",a) ;
        printf("s_ln   = %e, ",s_ln) ;
        printf("s_linf = %e, ",s_linf) ;
        printf("s_l    = %e\n",s_l) ;
        printf("ds_l   = %e, ",ds_l) ;
        printf("dbeta_i= %e\n",dbeta_i) ;
        printf("r_n    = %e, ",r_n) ;
        printf("r_inf  = %e, ",r_inf) ;
        printf("r      = %e\n",r) ;
        */
      }
      
      r += dr ;
      if(fabs(dr/r_n) < tol) {
        return(r) ;
      }
    }
  }
  
  Message_Warning("Radius: not converged") ;
  Exception_Interrupt ;

  return(r) ;
}




double PoreWallEquilibriumSaturationIndex(double beta_pn,double beta,double deltastrainv,double dt)
{
  double beta_p = beta_pn ;
  
  if(dt > 0) {
    double rate = deltastrainv / dt ;
    double iogr = InverseOfPoreCrystalGrowthRate(rate) ;
    
    beta_p = beta*iogr ;
    
    beta_p = Math_Max(beta_p,1) ;
  
    #if 0
    printf("\n") ;
    printf("rate = %e ",rate) ;
    printf("iogr = %e ",iogr) ;
    printf("beta_p = %e",beta_p) ;
    #endif
  }
  
  return(beta_p) ;
}



#if 0
double PoreWallEquilibriumSaturationIndex_old(double beta_pn,double varphi_cn,double strainv_n,double straind_n,double beta,double sc,double dt)
{
  double beta_p   = beta_pn ;
  //double sl       = 1 - sc ;
  
  {
    int iterations = 40 ;
    double tol = 1.e-6 ;
    int i ;
    
    for(i = 0 ; i < iterations ; i++) {
      double phicrate  =  PoreCrystalGrowthRate(sc,beta,beta_p) ;
      double dphicrate = dPoreCrystalGrowthRate(sc,beta,beta_p) ;
      double pc        = CrystallizationPressure(beta_p) ;
      double dpc       = dCrystallizationPressure(beta_p) ;
      //double varphi_c  = varphi_cn + dt*phicrate ;
      //double dvarphi_c =             dt*dphicrate ;
      double strainv   = strainv_n + ((sc > 0) ? dt*phicrate/sc : 0) ;
      double dstrainv   = (sc > 0) ? dt*dphicrate/sc : 0 ;
      //double strainv  = (sc > 0) ?  varphi_c/sc - ( pc/N_Biot +  pc*sl/G_Biot) : 0 ;
      //double dstrainv = (sc > 0) ? dvarphi_c/sc - (dpc/N_Biot + dpc*sl/G_Biot) : 0 ;
      /* Effective stress */
      double stress  =  ElasticDamageStress(strainv,straind_n) ;
      double dstress = dElasticDamageStress(strainv,straind_n) ;
      double straind = DamageStrain(strainv,straind_n) ;
      double d       = DamageCoefficient(straind) ;
      double bd      = Biot + d * (1 - Biot) ;
      double eq      = stress - bd*sc*pc ;
      /* We don't derive bd */
      double deq     = dstress*dstrainv - bd*sc*dpc ;
      double dbeta_p = (fabs(deq) > 0) ? - eq/deq : 0. ;
      
      if(dbeta_p < - 0.5*beta_p) dbeta_p = - 0.5*beta_p ;
      
      beta_p += dbeta_p ;
          
      //if(beta_p < 1) beta_p = 1 ;
      
      if(fabs(dbeta_p/beta_pn) < tol) {
        return(beta_p) ;
      }
    }
  }
  
  Message_Direct("PoreWallEquilibriumSaturationIndex: not converged\n") ;
  Message_Direct("beta_p = %e ; beta_pn = %e\n",beta_p,beta_pn) ;
  Exception_Interrupt ;
  
  return(-1) ;
}





double DamageStrain(double strain,double straind)
{
  //double Y = 0.5*strain*K_bulk*strain ;
  //double K = 0.5*straind*K_bulk*straind ;
  //double crit = Y - K ;
  double crit = strain - straind ;
  
  if(crit > 0) {
    straind = strain ;
  }
  
  return(straind) ;
}



double ElasticDamageStress(double strain,double straind_n)
{
  double straind = DamageStrain(strain,straind_n) ;
  double d       = DamageCoefficient(straind) ;
  double Kd      = (1 - d) * K_bulk ;
  double stress  = Kd * strain ;
  
  return(stress) ;
}



double dElasticDamageStress(double strain,double straind_n)
{
  double dstrain = 1.e-4 * strain ;
  double a       = 0.5 ;
  double strain2 = strain - (1 - a) * dstrain ;
  double stress2 = ElasticDamageStress(strain2,straind_n) ;
  double strain1 = strain + a * dstrain ;
  double stress1 = ElasticDamageStress(strain1,straind_n) ;
  double dstress = (stress1 - stress2) / dstrain ;
  
  return(dstress) ;
}
#endif
