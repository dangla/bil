/* General features of the model:
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


#define TEMPERATURE   (293)

#define TITLE   "External sulfate attack of concrete (2017)" 
#define AUTHORS "Gu-Dangla"

#include "PredefinedMethods.h"

/* Nb of equations: 6/5 with/without electroneutrality */
//#define NEQ     (5)
#define NEQ     (6)


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)

/* Nb of (im/ex)plicit and constant terms */
#define NVE    	((1 + CementSolutionDiffusion_NbOfConcentrations)*NN)
#define NVI     (7*NN*NN + 13*NN)
#define NV0     (2)

/* Equation Indexes */
#define E_S     (0)
#define E_q     (1)
#define E_Ca    (2)
#define E_K     (3)
#define E_Al    (4)
/* Comment the next line to suppress electroneutrality equation */
#define E_el    (5)

/* Primary Unknown Indexes */
#define U_C_H2SO4     (0)
#define U_PSI         (1)
#define U_ZN_Ca_S     (2)
#define U_C_K         (3)
#define U_ZN_Al_S     (4)
/* Comment the next line to suppress unknown C_OH */
#define U_C_OH        (5)


/* Compiling options */
#define NOLOG_U     1
#define LOG_U       2
#define Ln10        (2.302585093)
#define U_H2SO4     LOG_U
#define U_K         LOG_U
#define U_OH        LOG_U
#define U_Cl        NOLOG_U
#define EXPLICIT  1
#define IMPLICIT  2
#define U_PHI     IMPLICIT
#define ELECTRONEUTRALITY   IMPLICIT


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWN_n(n,i)   Element_GetValueOfNodalUnknown(el,u_n,n,i)



/* Names of nodal unknowns */
#if (U_H2SO4 == LOG_U)
  #define LogC_H2SO4(n)   (UNKNOWN(n,U_C_H2SO4))
  #define LogC_H2SO4n(n)  (UNKNOWN_n(n,U_C_H2SO4))
  #define C_H2SO4(n)      (pow(10,UNKNOWN(n,U_C_H2SO4)))
  #define C_H2SO4n(n)     (pow(10,UNKNOWN_n(n,U_C_H2SO4)))
#else
  #define C_H2SO4(n)      (UNKNOWN(n,U_C_H2SO4))
  #define C_H2SO4n(n)     (UNKNOWN_n(n,U_C_H2SO4))
  #define LogC_H2SO4(n)   (log10(UNKNOWN(n,U_C_H2SO4)))
  #define LogC_H2SO4n(n)  (log10(UNKNOWN_n(n,U_C_H2SO4)))
#endif

#define ZN_Ca_S(n)   (UNKNOWN(n,U_ZN_Ca_S))
#define ZN_Ca_Sn(n)  (UNKNOWN_n(n,U_ZN_Ca_S))

#define ZN_Si_S(n)   (UNKNOWN(n,U_ZN_Si_S))
#define ZN_Si_Sn(n)  (UNKNOWN_n(n,U_ZN_Si_S))

#define PSI(n)       (UNKNOWN(n,U_PSI))
#define PSIn(n)      (UNKNOWN_n(n,U_PSI))

#if (U_K == LOG_U)
  #define LogC_K(n)  (UNKNOWN(n,U_C_K))
  #define LogC_Kn(n) (UNKNOWN_n(n,U_C_K))
  #define C_K(n)     (pow(10,UNKNOWN(n,U_C_K)))
  #define C_Kn(n)    (pow(10,UNKNOWN_n(n,U_C_K)))
#else
  #define C_K(n)     (UNKNOWN(n,U_C_K))
  #define C_Kn(n)    (UNKNOWN_n(n,U_C_K))
  #define LogC_K(n)  (log10(UNKNOWN(n,U_C_K)))
  #define LogC_Kn(n) (log10(UNKNOWN_n(n,U_C_K)))
#endif

#if (U_Cl == LOG_U)
  #define LogC_Cl(n)   (UNKNOWN(n,U_C_Cl))
  #define LogC_Cln(n)  (UNKNOWN_n(n,U_C_Cl))
  #define C_Cl(n)      (pow(10,UNKNOWN(n,U_C_Cl)))
  #define C_Cln(n)     (pow(10,UNKNOWN_n(n,U_C_Cl)))
#else
  #define C_Cl(n)      (UNKNOWN(n,U_C_Cl))
  #define C_Cln(n)     (UNKNOWN_n(n,U_C_Cl))
  #define LogC_Cl(n)   (log10(UNKNOWN(n,U_C_Cl)))
  #define LogC_Cln(n)  (log10(UNKNOWN_n(n,U_C_Cl)))
#endif

#define ZN_Al_S(n)   (UNKNOWN(n,U_ZN_Al_S))
#define ZN_Al_Sn(n)  (UNKNOWN_n(n,U_ZN_Al_S))

#ifdef U_C_OH
  #if (U_OH == LOG_U)
    #define LogC_OH(n)   (UNKNOWN(n,U_C_OH))
    #define LogC_OHn(n)  (UNKNOWN_n(n,U_C_OH))
    #define C_OH(n)      (pow(10,UNKNOWN(n,U_C_OH)))
    #define C_OHn(n)     (pow(10,UNKNOWN_n(n,U_C_OH)))
  #else
    #define C_OH(n)      (UNKNOWN(n,U_C_OH))
    #define C_OHn(n)     (UNKNOWN_n(n,U_C_OH))
    #define LogC_OH(n)   (log10(UNKNOWN(n,U_C_OH)))
    #define LogC_OHn(n)  (log10(UNKNOWN_n(n,U_C_OH)))
  #endif
#endif




/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])


#define NW_S          (f)
#define NW_Sn         (f_n)
#define N_S(i)        MassAndFlux(NW_S,i,i)
#define N_Sn(i)       MassAndFlux(NW_Sn,i,i)
#define W_S(i,j)      MassAndFlux(NW_S,i,j)

#define NW_q          (f   + NN*NN)
#define NW_qn         (f_n + NN*NN)
#define N_q(i)        MassAndFlux(NW_q,i,i)
#define N_qn(i)       MassAndFlux(NW_qn,i,i)
#define W_q(i,j)      MassAndFlux(NW_q,i,j)

#define NW_Ca         (f   + 2*NN*NN)
#define NW_Can        (f_n + 2*NN*NN)
#define N_Ca(i)       MassAndFlux(NW_Ca,i,i)
#define N_Can(i)      MassAndFlux(NW_Can,i,i)
#define W_Ca(i,j)     MassAndFlux(NW_Ca,i,j)

#define NW_K          (f   + 3*NN*NN)
#define NW_Kn         (f_n + 3*NN*NN)
#define N_K(i)        MassAndFlux(NW_K,i,i)
#define N_Kn(i)       MassAndFlux(NW_Kn,i,i)
#define W_K(i,j)      MassAndFlux(NW_K,i,j)

#define NW_Si         (f   + 4*NN*NN)
#define NW_Sin        (f_n + 4*NN*NN)
#define N_Si(i)       MassAndFlux(NW_Si,i,i)
#define N_Sin(i)      MassAndFlux(NW_Sin,i,i)
#define W_Si(i,j)     MassAndFlux(NW_Si,i,j)

#define NW_Al         (f   + 5*NN*NN)
#define NW_Aln        (f_n + 5*NN*NN)
#define N_Al(i)       MassAndFlux(NW_Al,i,i)
#define N_Aln(i)      MassAndFlux(NW_Aln,i,i)
#define W_Al(i,j)     MassAndFlux(NW_Al,i,j)

#define NW_Cl         (f   + 6*NN*NN)
#define NW_Cln        (f_n + 6*NN*NN)
#define N_Cl(i)       MassAndFlux(NW_Cl,i,i)
#define N_Cln(i)      MassAndFlux(NW_Cln,i,i)
#define W_Cl(i,j)     MassAndFlux(NW_Cl,i,j)


#define N_CH(i)       (f   + 7*NN*NN)[i]
#define N_CHn(i)      (f_n + 7*NN*NN)[i]

#define N_CSH2(i)     (f   + 7*NN*NN + NN)[i]
#define N_CSH2n(i)    (f_n + 7*NN*NN + NN)[i]

#define N_AH3(i)      (f   + 7*NN*NN + 2*NN)[i]
#define N_AH3n(i)     (f_n + 7*NN*NN + 2*NN)[i]

#define N_AFm(i)      (f   + 7*NN*NN + 3*NN)[i]
#define N_AFmn(i)     (f_n + 7*NN*NN + 3*NN)[i]

#define N_AFt(i)      (f   + 7*NN*NN + 4*NN)[i]
#define N_AFtn(i)     (f_n + 7*NN*NN + 4*NN)[i]

#define N_C3AH6(i)    (f   + 7*NN*NN + 5*NN)[i]
#define N_C3AH6n(i)   (f_n + 7*NN*NN + 5*NN)[i]

#define PHI(i)        (f   + 7*NN*NN + 6*NN)[i]
#define PHIn(i)       (f_n + 7*NN*NN + 6*NN)[i]

#define PHI_C(i)      (f   + 7*NN*NN + 7*NN)[i]
#define PHI_Cn(i)     (f_n + 7*NN*NN + 7*NN)[i]


#ifndef U_C_OH
  #define C_OH(i)     (f   + 7*NN*NN + 8*NN)[i]
  #define C_OHn(i)    (f_n + 7*NN*NN + 8*NN)[i]
  
  #define LogC_OH(n)  (log10(C_OH(n)))
  #define LogC_OHn(n) (log10(C_OHn(n)))
#endif


#define PoreRadius(i)  (f   + 7*NN*NN + 9*NN)[i]
#define PoreRadiusn(i) (f_n + 7*NN*NN + 9*NN)[i]


#define Beta_p(i)      (f   + 7*NN*NN + 10*NN)[i]
#define Beta_pn(i)     (f_n + 7*NN*NN + 10*NN)[i]

#define Strain(i)      (f   + 7*NN*NN + 11*NN)[i]
#define Strain_n(i)    (f_n + 7*NN*NN + 11*NN)[i]

#define Straind(i)     (f   + 7*NN*NN + 12*NN)[i]
#define Straind_n(i)   (f_n + 7*NN*NN + 12*NN)[i]





/* Names used for explicit terms */
#define TransferCoefficient(f,n)  ((f) + (n)*NN)

#define TORTUOSITY        TransferCoefficient(va,0)

#define CONCENTRATION(i)  (TransferCoefficient(va,1) + (i)*CementSolutionDiffusion_NbOfConcentrations)



/* Names used for constant terms */
#define V_Cem0(n)   (v0[(0+n)])



/* Shorthands of some units */
#include "InternationalSystemOfUnits.h"

#define m     (InternationalSystemOfUnits_OneMeter)
#define m3    (m*m*m)
#define dm    (0.1*m)
#define cm    (0.01*m)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)
#define Pa    (InternationalSystemOfUnits_OnePascal)
#define MPa   (1.e6*Pa)
#define GPa   (1.e9*Pa)
#define J     (Pa*m3)




/* Material properties */
//#define SATURATION_CURVE           (Element_GetCurve(el) + 4)
#define SATURATION_CURVE           (satcurve)
#define LiquidSaturationDegree(r)  (Curve_ComputeValue(SATURATION_CURVE,r))
#define dLiquidSaturationDegree(r) (Curve_ComputeDerivative(SATURATION_CURVE,r))
#define PoreEntryRadiusMax         (Curve_GetXRange(SATURATION_CURVE)[1])


/*
  Solids
  CH   = Calcium Hydroxide (Portlandite)
  CSH2 = Calcium Sulfate Dihydrate (Gypsum)
  CSH  = Calcium Silicates Hydrate
  SH   = Amorphous Silica Gel
*/

/* C-S-H Properties */
#define MOLARVOLUMEOFCSH_CURVE           (Element_GetCurve(el) + 2)
//#define MolarVolumeOfCSH(q)    (Curve_ComputeValue(MOLARVOLUMEOFCSH_CURVE,q))
#define V_CSH                  (78.*cm3)
#define V_SH                   (43.*cm3)
#define MolarVolumeOfCSH(x)    ((x)/1.7*V_CSH + (1 - (x)/1.7)*V_SH)
#define CSHSolidContent(zn_si_s)       SiliconContentInCSH(zn_si_s)


/* CH Properties */
/* Molar volume of CH solid (dm3/mole) */
#define V_CH       (33.*cm3)      /* (33.e-3) */
#define CHSolidContent(zn_ca_s)        CalciumContentInCH(zn_ca_s)


/* CSH2 Properties */
/* Molar volume of CSH2 crystal (dm3/mole) */
#define V_CSH2     (75.*cm3)      /* (75.e-3) */
#define CSH2SolidContent_kin(n,s,dt)     MAX((n + dt*r_csh2*(s - 1)),0.)
#define CSH2SolidContent(n,s,dt)         CSH2SolidContent_kin(n,s,dt)

/* AH3 Properties (Al2O6H6) */
/* Molar volume of AH3 solid (dm3/mole) */
#define V_AH3      (64.44*cm3)
#define AH3SolidContent(zn_al_s)    (0.5*AluminiumContentInAH3(zn_al_s))


/* AFm Properties ((Ca4).(Al2).(SO4).12(OH).6(H2O)) */
/* Molar volume of AFm solid (dm3/mole) */
#define V_AFm      (311.26*cm3)      /* Thermochimie (ANDRA) */
//#define AFmSolidContent(n,s,dt)     (n*pow(s,dt/t_afm))
#define AFmSolidContent(n,s,dt)     MAX((n + dt*r_afm*(s - 1)),0.)


/* AFt Properties ((Ca6).(Al2).3(SO4).12(OH).26(H2O)) */
/* Molar volume of AFt solid (dm3/mole) */
#define V_AFt      (710.32*cm3)      /* Thermochimie (ANDRA) */
//#define AFtSolidContent(n,s,dt)     (n*pow(s,dt/t_aft))
#define AFtSolidContent(n,s,dt)     MAX((n + dt*r_aft*(s - 1)),0.)
/* Surface tension (N/m) */
#define Gamma_AFt   (0.1*Pa*m)


/* Crystal properties */
#define V_Crystal       V_AFt
#define Gamma_Crystal   Gamma_AFt
/* Crystal - pore wall interface */
/* Crystallization pressure */
#define CrystallizationPressure(beta) \
        RT/V_Crystal*log(beta)
        
#define dCrystallizationPressure(beta) \
        RT/V_Crystal/(beta)

/* Crystal - liquid interface */
/* Interface equilibrium saturation index */
#define InterfaceEquilibriumSaturationIndex(r) \
        (exp(2*Gamma_Crystal*V_Crystal/(RT*r)))

#define dInterfaceEquilibriumSaturationIndex(r) \
        (-2*Gamma_Crystal*V_Crystal/(RT*r*r)*InterfaceEquilibriumSaturationIndex(r))

#define InverseOfInterfaceEquilibriumSaturationIndex(b) \
        (2*Gamma_Crystal*V_Crystal/(RT*log(b)))


/* C3AH6 Properties ((Ca3).(Al2).12(OH)) */
/* Molar volume of C3AH6 solid (dm3/mole) */
#define V_C3AH6      (149.52*cm3)
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
static void    ComputeSecondaryVariables(Element_t*,double,double,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;
//static double* ComputeVariableFluxDerivatives(Element_t*,int,int,double) ;

static int     TangentCoefficients(Element_t*,double,double,double*) ;


static double  Radius(double,double,double,Element_t*) ;


static double TortuosityOhJang(double) ;
static double TortuosityBazantNajjar(double) ;

static void   ComputePhysicoChemicalProperties(void) ;

static double PoreWallEquilibriumSaturationIndex(double,double,double,double,double,double) ;
static double ElasticDamageStress(double,double) ;
static double dElasticDamageStress(double,double) ;
static double DamageStrain(double,double) ;

//#define LiquidTortuosity  TortuosityOhJang
#define LiquidTortuosity  TortuosityBazantNajjar




/* Parameters */
static double phi0 ;
static double r_afm,r_aft,r_c3ah6,r_csh2 ;
static double n_ca_ref,n_si_ref,n_al_ref ;
static double n_afm_0,n_aft_0,n_c3ah6_0,n_csh2_0 ;
static double a_AFt ;
static double RT ;
static Curve_t* satcurve ;
static double K_bulk ;
static double strain0 ;
static double strainf ;
static double a_p ;


#include "PhysicalConstant.h"
#include "Temperature.h"

void ComputePhysicoChemicalProperties(void)
{
  RT = PhysicalConstant_PerfectGasConstant * TEMPERATURE ;
}


static CementSolutionDiffusion_t* csd = NULL ;
static HardenedCementChemistry_t* hcc = NULL ;




#define NbOfVariables    (NEQ+56)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;


#define I_ZN_Ca_S      (NEQ+0)
#define I_ZN_Si_S      (NEQ+1)
#define I_ZN_Al_S      (NEQ+2)

#define I_N_Q          (NEQ+4)
#define I_N_S          (NEQ+5)
#define I_N_Ca         (NEQ+6)
#define I_N_Si         (NEQ+7)
#define I_N_K          (NEQ+8)
#define I_N_Cl         (NEQ+9)
#define I_N_Al         (NEQ+10)

#define I_N_CH         (NEQ+11)
#define I_N_CSH2       (NEQ+12)
#define I_N_AH3        (NEQ+13)
#define I_N_AFm        (NEQ+14)
#define I_N_AFt        (NEQ+15)
#define I_N_C3AH6      (NEQ+16)
#define I_N_CSH        (NEQ+17)


#define I_V_Cem        (NEQ+21)

#define I_PHI          (NEQ+22)
#define I_PHI_C        (NEQ+23)

#define I_V_CSH        (NEQ+24)



#define I_N_CHn        (NEQ+31)
#define I_N_CSH2n      (NEQ+32)
#define I_N_AH3n       (NEQ+33)
#define I_N_AFmn       (NEQ+34)
#define I_N_AFtn       (NEQ+35)
#define I_N_C3AH6n     (NEQ+36)
#define I_N_CSHn       (NEQ+37)

#define I_V_Cem0       (NEQ+41)


#define I_PHIn         (NEQ+42)
#define I_PHI_Cn       (NEQ+43)


#define I_C_OHn        (NEQ+45)

#define I_Radius       (NEQ+46)
#define I_Radiusn      (NEQ+47)

#define I_S_C          (NEQ+48)

#define I_P_C          (NEQ+49)

#define I_Beta_p       (NEQ+50)
#define I_Beta_pn      (NEQ+51)

#define I_Strain       (NEQ+52)
#define I_Strain_n     (NEQ+53)

#define I_Straind      (NEQ+54)
#define I_Straind_n    (NEQ+55)




  
  

#define NbOfVariableFluxes    (7)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;
//static double dVariableFluxes[NbOfVariableFluxes] ;

#define I_W_S           (0)
#define I_W_Ca          (1)
#define I_W_Si          (2)
#define I_W_K           (3)
#define I_W_Cl          (4)
#define I_W_q           (5)
#define I_W_Al          (6)


int pm(const char *s)
{
  if(strcmp(s,"porosity") == 0)        return (0) ;
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
  else if(strcmp(s,"a_AFt") == 0)      return (16) ;
  else if(strcmp(s,"K_bulk") == 0)     return (17) ;
  else if(strcmp(s,"Strain0") == 0)    return (18) ;
  else if(strcmp(s,"Strainf") == 0)    return (19) ;
  else if(strcmp(s,"A_p") == 0)        return (20) ;
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
  a_AFt     = GetProperty("a_AFt") ;
  K_bulk    = GetProperty("K_bulk") ;
  strain0   = GetProperty("Strain0") ;
  strainf   = GetProperty("Strainf") ;
  a_p       = GetProperty("A_p") ;
  
  satcurve  = Element_FindCurve(el,"S_r") ;
}


int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  Model_GetNbOfVariables(model) = NbOfVariables ;
  Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  
  Model_CopyNameOfEquation(model,E_S, "sulfur") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_q, "charge") ;
  Model_CopyNameOfEquation(model,E_K, "potassium") ;
  Model_CopyNameOfEquation(model,E_Al,"aluminium") ;
#ifdef E_el
  Model_CopyNameOfEquation(model,E_el,"electroneutrality") ;
#endif

#if (U_H2SO4 == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_H2SO4,"logc_h2so4") ;
#else
  Model_CopyNameOfUnknown(model,U_C_H2SO4,"c_h2so4") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Ca_S,"z_ca") ;
  Model_CopyNameOfUnknown(model,U_PSI,    "psi") ;
#if (U_K == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_K,    "logc_k") ;
#else
  Model_CopyNameOfUnknown(model,U_C_K,    "c_k") ;
#endif
  Model_CopyNameOfUnknown(model,U_ZN_Al_S,"z_al") ;
#ifdef U_C_OH
  #if (U_OH == LOG_U)
    Model_CopyNameOfUnknown(model,U_C_OH, "logc_oh") ;
  #else
    Model_CopyNameOfUnknown(model,U_C_OH, "c_oh") ;
  #endif
#endif
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 21 ;

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

    Material_ScanProperties(mat,datafile,pm) ;
  }

  {
    ComputePhysicoChemicalProperties() ;
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



int PrintModelChar(Model_t* model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The 5/6 equations are:\n") ;
  printf("\t- Mass balance of S      (sulfur)\n") ;
  printf("\t- Charge balance         (charge)\n") ;
  printf("\t- Mass balance of Ca     (calcium)\n") ;
  printf("\t- Mass balance of K      (potassium)\n") ;
  printf("\t- Mass balance of Al     (aluminium)\n") ;
#ifdef E_el
  printf("\t- Electroneutrality      (electroneutrality)\n") ;
#endif
  
  printf("\n") ;
  printf("The 5/6 primary unknowns are:\n") ;
  printf("\t- Sulfuric acid concentration     (c_h2so4 or logc_h2so4)\n") ;
  printf("\t- Electric potential              (psi)\n") ;
  printf("\t- Zeta unknown for calcium        (z_ca)\n") ;
  printf("\t- Potassium concentration         (c_k)\n") ;
  printf("\t- Zeta unknown for aluminium      (z_al)\n") ;
#ifdef U_C_OH
  printf("\t- Hydroxide ion concentration     (c_oh or logc_oh)\n") ;
#endif

  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length : dm !\n") ;
  printf("\t time   : s !\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;
  fprintf(ficd,"N_CH  = 6.1       # CH mole content (moles/L)\n") ;
  fprintf(ficd,"N_K   = 0.4       # K mole content  (moles/L)\n") ;
  fprintf(ficd,"N_AH3  = 0.4      # Al mole content (moles/L)\n") ;
  fprintf(ficd,"N_AFm  = 0.1      # AFm mole content (moles/L)\n") ;
  fprintf(ficd,"N_AFt  = 0.4      # AFt mole content (moles/L)\n") ;
  fprintf(ficd,"Curves = file     # Pore volume fraction curve:  r  S_r\n") ;
  fprintf(ficd,"Curves = solid    # File name: S_CH  X_CSH  Z_CSH  S_SH\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}


int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
/* Initialise les variables du systeme (f,va) */ 
{
  double* f = Element_GetImplicitTerm(el) ;
  double* v0 = Element_GetConstantTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      double zn_ca_s    = ZN_Ca_S(i) ;
      double zn_si_s    = 1 ;
      double zn_al_s    = ZN_Al_S(i) ;
  
      /* Solve cement chemistry */
      {
        double logc_h2so4 = LogC_H2SO4(i) ;
        double logc_na    = -99 ;
        double logc_k     = LogC_K(i) ;
        double logc_oh    = -7 ;
        double c_cl       = 0 ;
  
        HardenedCementChemistry_GetInput(hcc,SI_Ca) = MIN(zn_ca_s,0) ;
        HardenedCementChemistry_GetInput(hcc,SI_Si) = MIN(zn_si_s,0) ;
        HardenedCementChemistry_GetInput(hcc,SI_Al) = MIN(zn_al_s,0) ;
        HardenedCementChemistry_GetInput(hcc,LogC_H2SO4) = logc_h2so4 ;
        HardenedCementChemistry_GetInput(hcc,LogC_Na)  = logc_na ;
        HardenedCementChemistry_GetInput(hcc,LogC_K)   = logc_k ;
        HardenedCementChemistry_GetInput(hcc,LogC_OH)  = logc_oh ;
    
        HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = c_cl ;
  
        HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O_2) ;
      
        HardenedCementChemistry_SolveElectroneutrality(hcc) ;
      }
    
  
      {
        /* Liquid components */
        double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
        //double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
        //double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
        double c_oh   = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
    
        /* Solid contents */
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
        double v_cem      = V_CH*n_ch + v_csh*n_csh + V_AH3*n_ah3 + V_AFm*n_afm + V_AFt*n_aft + V_C3AH6*n_c3ah6  ;
        /*
        double v_gyp      = V_Gyp*n_csh2 ;
        double v_csh2     = V_CSH2*n_csh2 ;
        */

        /* Porosity */
        double phi_c = phi0 ;
        double phi   = phi_c - V_CSH2*n_csh2 ;
    

        /* Back up what is needed to compute components */
        N_CH(i)    = n_ch ;
        N_CSH2(i)  = n_csh2 ;
        N_AFm(i)   = n_afm ;
        N_AFt(i)   = n_aft ;
        N_C3AH6(i) = n_c3ah6 ;
        PHI(i)     = phi ;
        PHI_C(i)   = phi_c ;

        #ifdef U_C_OH
          #if (U_OH == LOG_U)
            LogC_OH(i) = log10(c_oh) ;
          #else
            C_OH(i)    = c_oh ;
          #endif
        #else
          C_OH(i)    = c_oh ;
        #endif

        V_Cem0(i)  = v_cem ;
      }
    
      {
        PoreRadius(i)    = PoreEntryRadiusMax ;
      
        Beta_p(i)        = 1 ;
        Strain(i)        = 0 ;
        Straind(i)       = strain0 ;
      }

    }
  }
  
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,f,0,0,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;

    
      /* Back up */
      N_S(i)  = x[I_N_S] ;
      N_Ca(i) = x[I_N_Ca] ;
      N_Si(i) = x[I_N_Si] ;
      N_K(i)  = x[I_N_K] ;
      N_Cl(i) = x[I_N_Cl] ;
      N_Al(i) = x[I_N_Al] ;
      /* charge density */
      N_q(i)  = x[I_N_Q] ;

    
      N_CH(i)    = x[I_N_CH] ;
      N_CSH2(i)  = x[I_N_CSH2] ;
      N_AFm(i)   = x[I_N_AFm] ;
      N_AFt(i)   = x[I_N_AFt] ;
      N_C3AH6(i) = x[I_N_C3AH6] ;
      PHI(i)     = x[I_PHI] ;
      PHI_C(i)   = x[I_PHI_C] ;

      #ifndef U_C_OH
        C_OH(i)    = x[I_C_OH] ;
      #endif
    
      PoreRadius(i)    = x[I_Radius] ;
      
      Beta_p(i)        = x[I_Beta_p] ;
      Strain(i)        = x[I_Strain] ;
      Straind(i)       = x[I_Straind] ;
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
        double* w = ComputeVariableFluxes(el,i,j) ;
        

        W_S(i,j)     = w[I_W_S  ] ;
        W_Ca(i,j)    = w[I_W_Ca ] ;
        W_Si(i,j)    = w[I_W_Si ] ;
        W_q(i,j)     = w[I_W_q  ] ;
        W_K(i,j)     = w[I_W_K  ] ;
        W_Cl(i,j)    = w[I_W_Cl ] ;
        W_Al(i,j)    = w[I_W_Al ] ;

        W_S(j,i)     = - w[I_W_S  ] ;
        W_Ca(j,i)    = - w[I_W_Ca ] ;
        W_Si(j,i)    = - w[I_W_Si ] ;
        W_q(j,i)     = - w[I_W_q  ] ;
        W_K(j,i)     = - w[I_W_K  ] ;
        W_Cl(j,i)    = - w[I_W_Cl ] ;
        W_Al(j,i)    = - w[I_W_Al ] ;
      }
    }
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/* Thermes explicites (va)  */
{
  double* f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /*
    Transfer coefficient
  */
  
  ComputeTransferCoefficients(el,u,f) ;
  

  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Contenus molaires */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
      double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    
      HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;

      /* Back up */

      /* Molar contents */
      N_S(i)  = x[I_N_S] ;
      N_Ca(i) = x[I_N_Ca] ;
      N_Si(i) = x[I_N_Si] ;
      N_K(i)  = x[I_N_K] ;
      N_Cl(i) = x[I_N_Cl] ;
      N_Al(i) = x[I_N_Al] ;
    
      /* Charge density */
      N_q(i)  = x[I_N_Q] ;

    
      N_CH(i)    = x[I_N_CH] ;
      N_CSH2(i)  = x[I_N_CSH2] ;
      N_AFm(i)   = x[I_N_AFm] ;
      N_AFt(i)   = x[I_N_AFt] ;
      N_C3AH6(i) = x[I_N_C3AH6] ;
      PHI(i)     = x[I_PHI] ;
      PHI_C(i)   = x[I_PHI_C] ;
#ifndef U_C_OH
      C_OH(i)    = x[I_C_OH] ;
#endif
    
      PoreRadius(i)    = x[I_Radius] ;
      
      Beta_p(i)        = x[I_Beta_p] ;
      Strain(i)        = x[I_Strain] ;
      Straind(i)       = x[I_Straind] ;


      {
        int test = 0 ;
        /*
        if(x[I_C_H2SO4] < 0.) test = -1 ;
        if(x[I_N_Ca_S]  < 0.) test = -1 ;
        if(x[I_N_Si_S]  < 0.) test = -1 ;
        if(x[I_N_Al_S]  < 0.) test = -1 ;
        */
        if(x[I_PHI]     < 0.) test = -1 ;
      
        if(test < 0) {
          double c_h2so4 = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,H2SO4) ;
          double c_oh = HardenedCementChemistry_GetAqueousConcentrationOf(hcc,OH) ;
          double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
          double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
          double s_ah3  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AH3) ;
          double s_afm  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFm) ;
          double s_aft  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFt) ;
          double s_c3ah6  = HardenedCementChemistry_GetSaturationIndexOf(hcc,C3AH6) ;


          double xx = Element_GetNodeCoordinate(el,i)[0] ;
        
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
          printf("zn_ca_s   = %e\n",x[I_ZN_Ca_S]) ;
          printf("zn_al_s   = %e\n",x[I_ZN_Al_S]) ;
          printf("c_oh      = %e\n",c_oh) ;
          return(-1) ;
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
        double* w = ComputeVariableFluxes(el,i,j) ;
        

        W_S(i,j)     = w[I_W_S  ] ;
        W_Ca(i,j)    = w[I_W_Ca ] ;
        W_Si(i,j)    = w[I_W_Si ] ;
        W_q(i,j)     = w[I_W_q  ] ;
        W_K(i,j)     = w[I_W_K  ] ;
        W_Cl(i,j)    = w[I_W_Cl ] ;
        W_Al(i,j)    = w[I_W_Al ] ;

        W_S(j,i)     = - w[I_W_S  ] ;
        W_Ca(j,i)    = - w[I_W_Ca ] ;
        W_Si(j,i)    = - w[I_W_Si ] ;
        W_q(j,i)     = - w[I_W_q  ] ;
        W_K(j,i)     = - w[I_W_K  ] ;
        W_Cl(j,i)    = - w[I_W_Cl ] ;
        W_Al(j,i)    = - w[I_W_Al ] ;
      }
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
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
  
  {
    FVM_t* fvm = FVM_GetInstance(el) ;
    double c[4*NEQ*NEQ] ;
    
    TangentCoefficients(el,t,dt,c) ;
    
    {
      double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
      for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
    }
  }
  


#if (U_H2SO4 == NOLOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_H2SO4)     /= Ln10*C_H2SO4(0) ;
      K(i,U_C_H2SO4+NEQ) /= Ln10*C_H2SO4(1) ;
    }
  }
#endif

#if (U_K == NOLOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_K)     /= Ln10*C_K(0) ;
      K(i,U_C_K+NEQ) /= Ln10*C_K(1) ;
    }
  }
#endif
  
#ifdef U_C_OH
  #if (U_OH == NOLOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_OH)     /= Ln10*C_OH(0) ;
      K(i,U_C_OH+NEQ) /= Ln10*C_OH(1) ;
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
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Conservation of element S: (N_S - N_Sn) + dt * div(W_S) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = N_S(i) - N_Sn(i) ;
        } else {
          g[i*nn + j] = dt * W_S(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_S) -= r1[i] ;
      }
    }
  }
  
  /*
    Conservation of element Ca: (N_Ca - N_Can) + dt * div(W_Ca) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = N_Ca(i) - N_Can(i) ;
        } else {
          g[i*nn + j] = dt * W_Ca(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_Ca) -= r1[i] ;
      }
    }
  }
  
  /*
    Conservation of element K: (N_K - N_Kn) + dt * div(W_K) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = N_K(i) - N_Kn(i) ;
        } else {
          g[i*nn + j] = dt * W_K(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_K) -= r1[i] ;
      }
    }
  }
  
  /*
    Conservation of element Al: (N_Al - N_Aln) + dt * div(W_Al) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = N_Al(i) - N_Aln(i) ;
        } else {
          g[i*nn + j] = dt * W_Al(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_Al) -= r1[i] ;
      }
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
  
  
#ifdef E_el
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
        R(i,E_el) -= r1[i] ;
      }
    }
  }
#endif
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/* Les valeurs exploitees (s) */
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nso = 58 ;
  double zero = 0 ;
  double one = 1 ;
  int    i ;

  /* if(Element_IsSubmanifold(el)) return(0) ; */
  
  /*
    Input data
  */
  GetProperties(el) ;


  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }


  /* output quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double* x = ComputeVariables(el,u,u,f,t,0,j) ;
    /* Concentrations */
#define ptC(CPD)   &(HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,CPD))
#define ptS(CPD)   &(HardenedCementChemistry_GetSaturationIndexOf(hcc,CPD))
#define ptPSI      &(HardenedCementChemistry_GetElectricPotential(hcc))
#define ptX_CSH    &(HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc))
#define ptZ_CSH    &(HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc))

    i = 0 ;
    
    {
      double ph        = - *(ptC(H )) ;
      
      Result_Store(r + i++,&ph,"ph",1) ;
    }
    
    Result_Store(r + i++,ptC(OH),"c_oh",1) ;
    Result_Store(r + i++,ptC(H ),"c_h",1) ;
    
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
    
    Result_Store(r + i++,ptC(CaSO4),"c_caso4aq",1) ;
    Result_Store(r + i++,ptC(CaHSO4 ),"c_cahso4",1) ;

    Result_Store(r + i++,ptC(K  ),"c_k",1) ;
    Result_Store(r + i++,ptC(KOH),"c_koh",1) ;
    
    Result_Store(r + i++,(x + I_ZN_Ca_S),"zn_ca_s",1) ;
    Result_Store(r + i++,&one           ,"zn_si_s",1) ;
    
    Result_Store(r + i++,ptS(CH),"s_ch",1) ;
    Result_Store(r + i++,ptS(CSH2),"s_csh2",1) ;
    
    Result_Store(r + i++,(x + I_N_CH  ),"n_ch",1) ;
    Result_Store(r + i++,(x + I_N_CSH2),"n_csh2",1) ;
    Result_Store(r + i++,(x + I_N_CSH ),"n_csh",1) ;
    
    Result_Store(r + i++,(x + I_PHI),"porosite",1) ;
    Result_Store(r + i++,ptPSI,"potentiel_electrique",1) ;
    
    Result_Store(r + i++,(x + I_N_Q),"charge",1) ;
    
    Result_Store(r + i++,(x + I_V_CSH),"V_CSH",1) ;
    Result_Store(r + i++,ptX_CSH,"C/S",1) ;
    
    Result_Store(r + i++,&W_Si(0,1),"W_Si",1) ;
    Result_Store(r + i++,&W_Ca(0,1),"W_Ca",1) ;
    Result_Store(r + i++,&W_S(0,1) ,"W_S",1) ;
    
    Result_Store(r + i++,&zero,"P_CSH2",1) ;
    Result_Store(r + i++,&zero,"Damage",1) ;
    
    Result_Store(r + i++,ptC(Al   ),"c_al",1) ;
    Result_Store(r + i++,ptC(AlO4H4),"c_alo4h4",1) ;
    
    Result_Store(r + i++,(x + I_ZN_Al_S),"zn_al_s",1) ;
    
    Result_Store(r + i++,ptS(AH3  ),"s_ah3",1) ;
    Result_Store(r + i++,ptS(AFm  ),"s_afm",1) ;
    Result_Store(r + i++,ptS(AFt  ),"s_aft",1) ;
    Result_Store(r + i++,ptS(C3AH6),"s_c3ah6",1) ;
    
    Result_Store(r + i++,(x + I_N_AH3  ),"n_ah3",1) ;
    Result_Store(r + i++,(x + I_N_AFm  ),"n_afm",1) ;
    Result_Store(r + i++,(x + I_N_AFt  ),"n_aft",1) ;
    Result_Store(r + i++,(x + I_N_C3AH6),"n_c3ah6",1) ;
    
    Result_Store(r + i++,&W_Al(0,1),"W_Al",1) ;
    
    Result_Store(r + i++,&W_q(0,1),"W_q",1) ;
    
    Result_Store(r + i++,&N_Ca(j),"N_Ca",1) ;
    Result_Store(r + i++,&N_Si(j),"N_Si",1) ;
    Result_Store(r + i++,&N_S(j) ,"N_S",1) ;
    Result_Store(r + i++,&N_Al(j),"N_Al",1) ;
    Result_Store(r + i++,&N_K(j) ,"N_K",1) ;
    Result_Store(r + i++,&N_Cl(j),"N_Cl",1) ;
    
    Result_Store(r + i++,(x + I_S_C    ),"Saturation degree of crystal",1) ;
    Result_Store(r + i++,(x + I_Radius ),"Pore entry radius",1) ;
    
    {
      double beta = exp(2*Gamma_AFt*V_AFt/(RT*x[I_Radius])) ;
      Result_Store(r + i++,&beta,"Equilibrium saturation index of AFt",1) ;
    }
    
    Result_Store(r + i++,(x + I_P_C    ),"Crystallization pressure",1) ;
    
    Result_Store(r + i++,(x + I_Strain ),"Strain",1) ;
    
    
    if(i != nso) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(nso) ;
}



void ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Transfer coefficients  */
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;

  /* Initialization */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;


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



int TangentCoefficients(Element_t* el,double t,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
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
    dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
  }

  
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* xi   = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
    double* mui = CementSolutionDiffusion_GetPotentialAtPoint(csd,i) ;
    int k ;
    
    HardenedCementChemistry_CopyChemicalPotential(hcc,mui) ;


    #if (U_H2SO4 == NOLOG_U)
    dui[U_C_H2SO4] =  1.e-2*ObVal_GetValue(obval + U_C_H2SO4)/(Ln10*xi[U_C_H2SO4]) ;
    #endif
    
    #if (U_K == NOLOG_U)
    dui[U_C_K] =  1.e-2*ObVal_GetValue(obval + U_C_K)/(Ln10*xi[U_C_K]) ;
    #endif
    
    #ifdef U_C_OH
      #if (U_OH == NOLOG_U)
      dui[U_C_OH] =  1.e-2*ObVal_GetValue(obval + U_C_OH)/(Ln10*xi[U_C_OH]) ;
      #endif
    #endif
    
    for(k = 0 ; k < NEQ ; k++) {
      double  dui_k  = dui[k] ;
      double* dxi    = ComputeVariableDerivatives(el,t,dt,xi,dui_k,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
        
        cii[E_S*NEQ    + k] = dxi[I_N_S ] ;
        cii[E_Ca*NEQ   + k] = dxi[I_N_Ca] ;
        cii[E_K*NEQ    + k] = dxi[I_N_K ] ;
        cii[E_Al*NEQ   + k] = dxi[I_N_Al] ;
#ifdef E_el
        //cii[E_q*NEQ    + k] = dxi[I_N_Q] ;
        cii[E_el*NEQ   + k] = dxi[I_N_Q] ;
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
        
              cij[E_S*NEQ    + k] = - dtdij * dw[I_W_S ] ;
              cij[E_Ca*NEQ   + k] = - dtdij * dw[I_W_Ca] ;
              cij[E_K*NEQ    + k] = - dtdij * dw[I_W_K ] ;
              cij[E_Al*NEQ   + k] = - dtdij * dw[I_W_Al] ;
              cij[E_q*NEQ    + k] = - dtdij * dw[I_W_q ] ;
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
  /*
  Model_t* model = Element_GetModel(el) ;
  double* x      = Model_GetVariable(model,n) ;
  */
  double* x = Variables[n] ;
  
  /* Primary Variables */
  x[U_C_H2SO4 ] = LogC_H2SO4(n) ;
  x[U_ZN_Ca_S ] = ZN_Ca_S(n) ;
  x[U_C_K     ] = LogC_K(n) ;
  x[U_PSI     ] = PSI(n) ;
  x[U_ZN_Al_S ] = ZN_Al_S(n) ;
#ifdef U_C_OH
  x[U_C_OH    ] = LogC_OH(n) ;
#endif
  
  /* Needed variables to compute secondary components */
  x[I_N_CHn  ]  = N_CHn(n) ;
  x[I_N_CSH2n]  = N_CSH2n(n) ;
  x[I_N_AFmn ]  = N_AFmn(n) ;
  x[I_N_AFtn ]  = N_AFtn(n) ;
  x[I_N_C3AH6n] = N_C3AH6n(n) ;
  x[I_PHIn   ]  = PHIn(n) ;
  x[I_PHI_Cn ]  = PHI_Cn(n) ;
  x[I_C_OHn  ]  = C_OHn(n) ;
  x[I_Radiusn]  = PoreRadiusn(n) ;
  x[I_Beta_pn]  = Beta_pn(n) ;
  x[I_Strain_n] = Strain_n(n);
  x[I_Straind_n] = Straind_n(n);

  {
    double* v0   = Element_GetConstantTerm(el) ;
    
    x[I_V_Cem0 ]  = V_Cem0(n) ;
  }
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  /* Primary variables */
  double zn_si_s    = 1 ;
  double zn_ca_s    = x[U_ZN_Ca_S] ;
  double zn_al_s    = x[U_ZN_Al_S] ;
  double c_cl       = 0 ;
  double psi        = x[U_PSI] ;
  
  /* Solve cement chemistry */
  {
    double logc_h2so4 = x[U_C_H2SO4] ;
    double logc_na    = -99 ;
    double logc_k     = x[U_C_K] ;
#ifdef U_C_OH
    double logc_oh    = x[U_C_OH] ;
#else
    double logc_oh    = log10(x[I_C_OHn]) ;
#endif
  
    HardenedCementChemistry_GetInput(hcc,SI_Ca) = MIN(zn_ca_s,0) ;
    HardenedCementChemistry_GetInput(hcc,SI_Si) = MIN(zn_si_s,0) ;
    HardenedCementChemistry_GetInput(hcc,SI_Al) = MIN(zn_al_s,0) ;
    HardenedCementChemistry_GetInput(hcc,LogC_H2SO4) = logc_h2so4 ;
    HardenedCementChemistry_GetInput(hcc,LogC_Na)  = logc_na ;
    HardenedCementChemistry_GetInput(hcc,LogC_K)   = logc_k ;
    HardenedCementChemistry_GetInput(hcc,LogC_OH)  = logc_oh ;
    HardenedCementChemistry_GetElectricPotential(hcc) = psi ;
    
    HardenedCementChemistry_GetAqueousConcentrationOf(hcc,Cl) = c_cl ;
  
    HardenedCementChemistry_ComputeSystem(hcc,CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O_2) ;

#ifndef E_el
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
  
  //double s_ch   = HardenedCementChemistry_GetSaturationIndexOf(hcc,CH) ;
  //double s_sh   = HardenedCementChemistry_GetSaturationIndexOf(hcc,SH) ;
  double s_csh2 = HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2) ;
  //double s_ah3  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AH3) ;
  double s_afm  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFm) ;
  double s_aft  = HardenedCementChemistry_GetSaturationIndexOf(hcc,AFt) ;
  double s_c3ah6 = HardenedCementChemistry_GetSaturationIndexOf(hcc,C3AH6) ;
       
    
    
  /* Compute the crystal saturation as a function of s_aft */
  double r_n = x[I_Radiusn] ;
  double r   = Radius(r_n,s_aft,dt,el) ;
  double s_l = LiquidSaturationDegree(r) ;
  double s_c = 1 - s_l ;
  
  /* Compute the strain */
  double beta_pn   = x[I_Beta_pn] ;
  double strain_n  = x[I_Strain_n] ;
  double straind_n = x[I_Straind_n] ;
  double beta_p = PoreWallEquilibriumSaturationIndex(beta_pn,strain_n,straind_n,s_aft,s_c,dt) ;
  double strain  = strain_n + dt*a_p*(1 - beta_p/s_aft) ;
  /* Compute the crystal pore deformation */
  double varphi_c = s_c * strain ;
  
  
  /* Solid contents */
  /* ... as components: CH, CSH2, CSH, AH3, AFm, AFt, C3AH6 */
  //double n_chn      = x[I_N_CHn] ;
  double n_ch       = CHSolidContent(zn_ca_s) ;
  double n_csh2n    = x[I_N_CSH2n] ;
  double n_csh2     = CSH2SolidContent(n_csh2n,s_csh2,dt) ;
  double n_ah3      = AH3SolidContent(zn_al_s) ;
  double n_afmn     = x[I_N_AFmn] ;
  double n_afm      = AFmSolidContent(n_afmn,s_afm,dt) ;
  //double n_aftn     = x[I_N_AFtn] ;
  //double n_aft      = AFtSolidContent(n_aftn,s_aft,dt) ;
  double n_aft      = (phi0*s_c + varphi_c)/V_AFt ;
  double n_c3ah6n   = x[I_N_C3AH6n] ;
  double n_c3ah6    = C3AH6SolidContent(n_c3ah6n,s_c3ah6,dt) ;
  double n_csh      = CSHSolidContent(zn_si_s) ;
  /* ... as elements: S, Ca, Si, Al */
  //double x_csh      = CalciumSiliconRatioInCSH(s_ch) ;
  double x_csh      = HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) ;
  double n_si_s     = n_csh ;
  double n_ca_s     = n_ch + n_csh2 + x_csh*n_csh + 4*n_afm + 6*n_aft + 3*n_c3ah6 ;
  double n_s_s      = n_csh2 + n_afm + 3*n_aft ;
  double n_al_s     = 2*(n_ah3 + n_afm + n_aft + n_c3ah6) ;
  /* ... as volume */
  double v_csh      = MolarVolumeOfCSH(x_csh) ;
  double v_cem      = V_CH*n_ch + v_csh*n_csh + V_AH3*n_ah3 + V_AFm*n_afm + V_AFt*n_aft + V_C3AH6*n_c3ah6 ;
  double v_gyp      = V_Gyp*n_csh2 ;
  double v_csh2     = V_CSH2*n_csh2 ;


  /* Porosities in unconfined conditions (no pressure) */
  double v_cem0     = x[I_V_Cem0] ;
  /* ... of concrete */
  double phi_con    = phi0 + v_cem0 - v_cem ;
  /* ... of gypsum */
  double phi_gyp    = PHI_Gyp ;


  /* total porosity */
  double varphi     = strain ;
  double phi_c      = phi_con + varphi ;
  double phi_t      = phi_con - v_csh2 ;
  

#if (U_PHI == IMPLICIT)
  double phi_l        = phi_t ;
#else
  double phi_l        = x[I_PHIn] ;
#endif
    
    
  /* Liquid contents */
  /* ... as elements: S, Ca, Si, K, Cl, Al */
  double c_cl_l = c_cl ;
  double n_s_l  = phi_l*c_s_l ;
  double n_ca_l = phi_l*c_ca_l ;
  double n_si_l = phi_l*c_si_l ;
  double n_al_l = phi_l*c_al_l ;
  double n_k_l  = phi_l*c_k_l ;
  double n_cl_l = phi_l*c_cl_l ;
  double n_q_l  = phi_l*c_q_l ;


#ifndef E_el
  #if (ELECTRONEUTRALITY == EXPLICIT)
  c_q_l = HardenedCementChemistry_SolveExplicitElectroneutrality(hcc) ;
  n_q_l = phi_l*c_q_l ;
  #endif
#endif

  

  /* Back up */


  /* Solid components */
  x[I_N_CH      ] = n_ch ;
  x[I_N_CSH2    ] = n_csh2 ;
  x[I_N_AH3     ] = n_ah3 ;
  x[I_N_AFm     ] = n_afm ;
  x[I_N_AFt     ] = n_aft ;
  x[I_N_C3AH6   ] = n_c3ah6 ;
  x[I_N_CSH     ] = n_csh ;
  
  x[I_ZN_Ca_S   ] = zn_ca_s ;  
  x[I_ZN_Al_S   ] = zn_al_s ;  
  
  x[I_V_CSH     ] = v_csh ;
  
  x[I_V_Cem     ] = v_cem ;
  
  
  /* Porosities */
  x[I_PHI       ] = phi_t ;
  x[I_PHI_C     ] = phi_c ;
  
  /* Saturation degree of crystal */
  x[I_S_C       ] = s_c ;
  
  /* Crystallization pressure */
  x[I_P_C       ] = CrystallizationPressure(beta_p) ;
  
  /* Radius */
  x[I_Radius    ] = r ;
  
  /* Strains */
  x[I_Beta_p    ] = beta_p ;
  x[I_Strain    ] = strain ;
  x[I_Straind   ] = DamageStrain(strain,straind_n) ;
  
  
  
  /* Element contents */
  x[I_N_S       ] = n_s_l  + n_s_s ;
  x[I_N_Ca      ] = n_ca_l + n_ca_s ;
  x[I_N_Si      ] = n_si_l + n_si_s ;
  x[I_N_K       ] = n_k_l  ;
  x[I_N_Cl      ] = n_cl_l  ;
  x[I_N_Al      ] = n_al_l + n_al_s ;

  /* Charge density */
  x[I_N_Q       ] = n_q_l ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dui,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dui ;
  
  ComputeSecondaryVariables(el,t,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dui ;
  }

  return(dx) ;
}




double* ComputeVariableFluxes(Element_t* el,int i,int j)
{
  double* grdij = dVariables ;

  /* Gradients */
  {
    int nn = Element_GetNbOfNodes(el) ;
    FVM_t* fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    double dij  = dist[nn*i + j] ;

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



double TortuosityOhJang(double phi)
/** Liquid totuosity according to Oh and Jang, CCR2003 */
{
  double phi_cap = 0.5*phi  ;
  double phi_c   = 0.17 ;  /* Percolation capilar porosity */
    
  double n = 2.7 ; 		      /* OPC n = 2.7,        Fly ash n = 4.5 */
  double ds_norm = 1.e-4 ;	/* OPC ds_norm = 1e-4, Fly ash ds_norm = 5e-5 */
  double m_phi = 0.5*(pow(ds_norm,1/n) + phi_cap/(1-phi_c)*(1 - pow(ds_norm,1/n)) - phi_c/(1-phi_c)) ;
  double iff =  pow(m_phi + sqrt(m_phi*m_phi +  pow(ds_norm,1/n)*phi_c/(1-phi_c)),n) ;
    
  return(iff) ;
}



double TortuosityBazantNajjar(double phi)
/** Liquid totuosity according to Bazant and Najjar */
{
  double iff = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
  return(iff) ;
}



double Radius(double r_n,double s_aft,double dt,Element_t* el)
{
  double r_max  = PoreEntryRadiusMax ;
  double r = r_n ;
  
  //if(s_aft < 1) return(r_max) ;
  
  {
    double s_ln   = LiquidSaturationDegree(r_n) ;
    double beta_n = InterfaceEquilibriumSaturationIndex(r_n) ;
    double beta_min = InterfaceEquilibriumSaturationIndex(r_max) ;
    double beta_inf = (s_aft > beta_min) ? s_aft : beta_min ;
    double r_inf = (s_aft > beta_min) ? InverseOfInterfaceEquilibriumSaturationIndex(s_aft) : r_max ;
    double s_linf  = LiquidSaturationDegree(r_inf) ;
    int iterations = 40 ;
    double tol = 1.e-6 ;
    int i ;
    
    if(r_n == r_inf) return(r_n) ;
    
    for(i = 0 ; i < iterations ; i++) {
      double s_l   = LiquidSaturationDegree(r) ;
      double ds_l  = dLiquidSaturationDegree(r) ;
      //double ds_l  = (s_linf - s_ln)/(r_inf - r_n) ;
      double beta  = InterfaceEquilibriumSaturationIndex(r) ;
      double dbeta = dInterfaceEquilibriumSaturationIndex(r) ;
      //double dbeta = (beta_inf - beta_n)/(r_inf - r_n) ;
      
      /* Kinetic law */
      /* Modified 03/06/2017 */
      double eq    = s_l - s_ln + dt*a_AFt*(1 - beta/s_aft) ;
      double deq   = ds_l - dt*a_AFt*dbeta/s_aft ;
      //double eq    = s_l - s_ln + dt*a_AFt*(s_aft/beta - 1) ;
      //double deq   = ds_l - dt*a_AFt*s_aft*dbeta/beta/beta ;
      
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
        printf("dbeta  = %e\n",dbeta) ;
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




double PoreWallEquilibriumSaturationIndex(double beta_pn,double strain_n,double straind_n,double s_aft,double sc,double dt)
{
  double beta_p   = beta_pn ;
  
  {
    int iterations = 40 ;
    double tol = 1.e-6 ;
    int i ;
    
    for(i = 0 ; i < iterations ; i++) {
      double strain  = strain_n + dt*a_p*(1 - beta_p/s_aft) ;
      double dstrain = - dt*a_p/s_aft ;
      double stress  =  ElasticDamageStress(strain,straind_n) ;
      double dstress = dElasticDamageStress(strain,straind_n) ;
      double pc  =  CrystallizationPressure(beta_p) ;
      double dpc = dCrystallizationPressure(beta_p) ;
      double eq  = stress - sc*pc ;
      double deq = dstress*dstrain - sc*dpc  ;
      double dbeta_p  = (fabs(deq) > 0) ? - eq/deq : 0. ;
      
      if(dbeta_p < - 0.5*beta_p) dbeta_p = - 0.5*beta_p ;
      
      beta_p += dbeta_p ;
          
      //if(beta_p < 1) beta_p = 1 ;
      
      if(fabs(dbeta_p/beta_pn) < tol) {
        return(beta_p) ;
      }
    }
  }
  
  Message_Warning("PoreWallEquilibriumSaturationIndex: not converged") ;
  Exception_Interrupt ;
  
  return(-1) ;
}




double DamageStrain(double strain,double straind)
{
  double Y = 0.5*strain*K_bulk*strain ;
  double K = 0.5*straind*K_bulk*straind ;
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
  double Kd      = K_bulk * strain0/straind*exp(-(straind - strain0)/strainf) ;
  double stress  = Kd * strain ;
  
  return(stress) ;
}



double dElasticDamageStress(double strain,double straind_n)
{
  double dstrain = 1.e-4 * strain ;
  double strain2 = strain - dstrain ;
  double stress2 = ElasticDamageStress(strain2,straind_n) ;
  double strain1 = strain + dstrain ;
  double stress1 = ElasticDamageStress(strain1,straind_n) ;
  double dstress = (stress1 - stress2) * 0.5 / dstrain ;
  
  return(dstress) ;
}
