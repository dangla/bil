#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Message.h"
#include "Exception.h"
#include "Tools/Math.h"
#include "Temperature.h"
#include "CementSolutionChemistry.h"

#define Ln10      Math_Ln10


/* Units
 * ----- */
#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define dm    (0.1*InternationalSystemOfUnits_OneMeter)
#define cm    (0.01*InternationalSystemOfUnits_OneMeter)
#define dm2   (dm*dm)
#define dm3   (dm*dm*dm)
#define Liter dm3
#define cm3   (cm*cm*cm)
#define MPa   (1.e6*InternationalSystemOfUnits_OnePascal)
#define GPa   (1.e9*InternationalSystemOfUnits_OnePascal)
#define mol   InternationalSystemOfUnits_OneMole
#define sec   InternationalSystemOfUnits_OneSecond


static double poly4(double,double,double,double,double,double,double) ;

static void   (CementSolutionChemistry_AllocateMemory)             (CementSolutionChemistry_t*) ;
static void   (CementSolutionChemistry_UpdateChemicalConstants)    (CementSolutionChemistry_t*) ;
static double (CementSolutionChemistry_ComputeChargeDensity)       (CementSolutionChemistry_t*) ;
static double (CementSolutionChemistry_ComputeLiquidMassDensity)   (CementSolutionChemistry_t*) ;
static void   (CementSolutionChemistry_UpdateElementConcentrations)(CementSolutionChemistry_t*) ;
static void   (CementSolutionChemistry_Initialize)                 (CementSolutionChemistry_t*) ;
static void   (CementSolutionChemistry_TranslateConcentrationsIntoActivities)(CementSolutionChemistry_t* csc) ;



static double*   (CementSolutionChemistry_CreateValence)(void) ;
static void      (CementSolutionChemistry_InitializeValence)(double*) ;


static void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_1(CementSolutionChemistry_t*) ;

static double* instancevalence = NULL ;



double* CementSolutionChemistry_GetValence()
{
  if(!instancevalence) {
    instancevalence = CementSolutionChemistry_CreateValence() ;
  }
  
  return(instancevalence) ;
}


CementSolutionChemistry_t* CementSolutionChemistry_Create(const int n)
{
  CementSolutionChemistry_t* csc = (CementSolutionChemistry_t*) malloc(n*sizeof(CementSolutionChemistry_t)) ;
  
  assert(csc) ;
  
  {
    int i ;
      
    for(i = 0 ; i < n ; i++) {
      CementSolutionChemistry_t* csci = csc + i ;
    
      /* Memory allocation */
      CementSolutionChemistry_AllocateMemory(csci) ;
  
      /* Initialize the equilibrium constants */
      CementSolutionChemistry_UpdateChemicalConstants(csci) ;
  
  
      /* Initialize concentrations to zero and other related variables */
      CementSolutionChemistry_Initialize(csci) ;
    }
  }
  
  return(csc) ;
}



void CementSolutionChemistry_AllocateMemory(CementSolutionChemistry_t* csc)
{
  
  
  /* Allocation of space for the temperature */
  {
    Temperature_t* temp = Temperature_Create() ;
    
    CementSolutionChemistry_GetTemperature(csc) = temp ;
  }
  
  
  /* Allocation of space for the primary variable indexes */
  {
    size_t sz = CementSolutionChemistry_NbOfPrimaryVariables*sizeof(int) ;
    int* ind = (int*) malloc(sz) ;
    
    assert(ind) ;
    
    CementSolutionChemistry_GetPrimaryVariableIndex(csc) = ind ;
  }
  
  
  /* Allocation of space for the primary variables */
  {
    size_t sz = CementSolutionChemistry_NbOfPrimaryVariables*sizeof(double) ;
    double* var = (double*) malloc(sz) ;
    
    assert(var) ;
    
    CementSolutionChemistry_GetPrimaryVariable(csc) = var ;
  }
  
  
  /* Allocation of space for the activities */
  {
    size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
    double* a = (double*) malloc(sz) ;
    
    assert(a) ;
    
    CementSolutionChemistry_GetActivity(csc) = a ;
  }
  
  
  /* Allocation of space for the log of activities */
  {
    size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
    double* loga = (double*) malloc(sz) ;
    
    assert(loga) ;
    
    CementSolutionChemistry_GetLogActivity(csc) = loga ;
  }
  
  
  /* Allocation of space for the concentrations */
  {
    size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    CementSolutionChemistry_GetConcentration(csc) = c ;
  }
  
  
  /* Allocation of space for the log of concentrations */
  {
    size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
    double* logc = (double*) malloc(sz) ;
    
    assert(logc) ;
    
    CementSolutionChemistry_GetLogConcentration(csc) = logc ;
  }
  
  
  /* Allocation of space for the element concentrations */
  {
    size_t sz = CementSolutionChemistry_NbOfElementConcentrations*sizeof(double) ;
    double* ec = (double*) malloc(sz) ;
    
    assert(ec) ;
    
    CementSolutionChemistry_GetElementConcentration(csc) = ec ;
  }
  
  
  /* Allocation of space for other variables */
  {
    size_t sz = CementSolutionChemistry_NbOfOtherVariables*sizeof(double) ;
    double* var = (double*) malloc(sz) ;
    
    assert(var) ;
    
    CementSolutionChemistry_GetOtherVariable(csc) = var ;
  }
  
  
  /* Allocation of space for the equilibrium constants */
  {
    size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
    double* keq = (double*) malloc(sz) ;
    
    assert(keq) ;
    
    CementSolutionChemistry_GetLog10Keq(csc) = keq ;
  }
}



double* CementSolutionChemistry_CreateValence(void)
{
  size_t sz = CementSolutionChemistry_NbOfSpecies*sizeof(double) ;
  double* z = (double*) malloc(sz) ;
    
  assert(z) ;
    
  CementSolutionChemistry_InitializeValence(z) ;
  
  return(z);
}



#include "ElectricChargeOfIonInWater.h"


#define CementSolutionChemistry_GetValenceOf(Z,CPD) \
       ((Z)[CementSolutionChemistry_GetIndexOf(CPD)])


void CementSolutionChemistry_InitializeValence(double* z)
{
  int     n = CementSolutionChemistry_NbOfSpecies ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    z[i] = 0. ;
  }
  
/* Shorthands of macros */
#define Z(A)           ElectricChargeOfIonInWater(A)
#define Valence(CPD)   CementSolutionChemistry_GetValenceOf(z,CPD)
  {
    Valence(OH)       = Z(OH) ;
    Valence(H)        = Z(H) ;
  
    Valence(Ca)       = Z(Ca) ;
    Valence(CaOH)     = Z(CaOH) ;
    Valence(CaO2H2)   = Z(CaO2H2) ;
  
    Valence(H4SiO4)   = Z(H4SiO4) ;
    Valence(H3SiO4)   = Z(H3SiO4) ;
    Valence(H2SiO4)   = Z(H2SiO4) ;
  
    Valence(Na)       = Z(Na) ;
    Valence(NaOH)     = Z(NaOH) ;
  
    Valence(K)        = Z(K) ;
    Valence(KOH)      = Z(KOH) ;
  
    Valence(Cl)       = Z(Cl) ;
  
    Valence(H2CO3)    = Z(H2CO3) ;
    Valence(HCO3)     = Z(HCO3) ;
    Valence(CO3)      = Z(CO3) ;
    Valence(CO2)      = Z(CO2) ;
  
    Valence(H2SO4)    = Z(H2SO4) ;
    Valence(HSO4)     = Z(HSO4) ;
    Valence(SO4)      = Z(SO4) ;
  
    Valence(Al)       = Z(Al) ;
    Valence(AlO4H4)   = Z(AlO4H4) ;
  
    Valence(NaHCO3)   = Z(NaHCO3) ;
    Valence(NaCO3)    = Z(NaCO3) ;
  
    Valence(CaHCO3)   = Z(CaHCO3) ;
    Valence(CaCO3)    = Z(CaCO3) ;
  
    Valence(CaH2SiO4) = Z(CaH2SiO4) ;
    Valence(CaH3SiO4) = Z(CaH3SiO4) ;

    Valence(CaHSO4)   = Z(CaHSO4) ;
    Valence(CaSO4)    = Z(CaSO4) ;
  }
#undef Valence
#undef Z
}



#include "Log10EquilibriumConstantOfHomogeneousReactionInWater.h"


#define LogKeq(CPD)  CementSolutionChemistry_GetLog10EquilibriumConstant(csc,CPD)




void CementSolutionChemistry_UpdateChemicalConstants(CementSolutionChemistry_t* csc)
{
  double T = CementSolutionChemistry_GetRoomTemperature(csc) ;
  
  #define LogKr(R) Log10EquilibriumConstantOfHomogeneousReactionInWater(R,T)
  double logk_h2o      = LogKr(H2O__H_OH) ;
  
  /* Compounds of type I */
  /* Calcium compounds */
  double logk_caoh     = LogKr(CaOH__Ca_OH) ;
  double logk_caoh2    = LogKr(CaO2H2__Ca_2OH) ;
  
  /* Silicon compounds */
  double logk_h3sio4   = LogKr(H3SiO4_H2O__H4SiO4_OH) ;
  double logk_h2sio4   = LogKr(H2SiO4_H2O__H3SiO4_OH) ;
  
  /* Sodium compounds */
  double logk_naoh     = LogKr(NaOH__Na_OH) ;
  
  /* Potassium compounds */
  double logk_koh      = LogKr(KOH__K_OH) ;
  
  /* Carbon compounds */
  double logk_h2co3    = LogKr(H2CO3__CO2_H2O) ;
  double logk_hco3     = LogKr(HCO3_H2O__H2CO3_OH) ;
  double logk_co3      = LogKr(CO3_H2O__HCO3_OH) ;
  
  /* Sulfur compounds */
  double logk_h2so4    = LogKr(H2SO4__HSO4_H) ;
  double logk_hso4     = LogKr(HSO4__SO4_H) ;
  
  /* Aluminium compounds */
  double logk_alo4h4   = LogKr(AlO4H4__Al_4OH) ;
  
  
  /* Compounds of type II */
  /* Calcium-Silicon compounds */
  double logk_cah3sio4 = LogKr(CaH3SiO4__Ca_H3SiO4) ;
  double logk_cah2sio4 = LogKr(CaH2SiO4__Ca_H2SiO4) ;
  
  /* Calcium-Carbon compounds */
  double logk_cahco3   = LogKr(CaHCO3__Ca_HCO3) ;
  double logk_caco3    = LogKr(CaCO3__Ca_CO3) ;
  
  /* Sodium-Carbon compounds */
  double logk_nahco3   = LogKr(NaHCO3__Na_HCO3) ;
  double logk_naco3    = LogKr(NaCO3__Na_CO3) ;
  
  /* Calcium-Sulfur compounds */
  double logk_cahso4   = LogKr(CaHSO4__Ca_HSO4) ;
  double logk_caso4    = LogKr(CaSO4__Ca_SO4) ;
  
  #undef LogKr

  
  
  /* Backup */
  LogKeq(H2O)      = logk_h2o ;
  
  LogKeq(CaOH)     = logk_caoh ;
  LogKeq(CaO2H2)   = logk_caoh2 ;
  
  LogKeq(H3SiO4)   = logk_h3sio4 ;
  LogKeq(H2SiO4)   = logk_h2sio4 ;
  
  LogKeq(NaOH)     = logk_naoh ;
  
  LogKeq(KOH)      = logk_koh ;
  
  LogKeq(H2CO3)    = logk_h2co3 ;
  LogKeq(HCO3)     = logk_hco3 ;
  LogKeq(CO3)      = logk_co3 ;
  
  LogKeq(H2SO4)    = logk_h2so4 ;
  LogKeq(HSO4)     = logk_hso4 ;
  
  LogKeq(AlO4H4)   = logk_alo4h4 ;
  
  LogKeq(CaH3SiO4) = logk_cah3sio4 ;
  LogKeq(CaH2SiO4) = logk_cah2sio4 ;
  
  LogKeq(CaHCO3)   = logk_cahco3 ;
  LogKeq(CaCO3)    = logk_caco3 ;
  
  LogKeq(NaHCO3)   = logk_nahco3 ;
  LogKeq(NaCO3)    = logk_naco3 ;
  
  LogKeq(CaHSO4)   = logk_cahso4 ;
  LogKeq(CaSO4)    = logk_caso4 ;
}



void CementSolutionChemistry_PrintChemicalConstants(CementSolutionChemistry_t* csc)
{
  double T = CementSolutionChemistry_GetRoomTemperature(csc) ;
  
  Log10EquilibriumConstantOfHomogeneousReactionInWater_Print(T) ;
  
  printf("\n") ;
  
  fflush(stdout) ;
}



/* Shorthands of macros */
#define Input(U)              CementSolutionChemistry_GetInput(csc,U)
#define Concentration(CPD)    CementSolutionChemistry_GetConcentrationOf(csc,CPD)
#define LogConcentration(CPD) CementSolutionChemistry_GetLogConcentrationOf(csc,CPD)
#define Activity(CPD)         CementSolutionChemistry_GetActivityOf(csc,CPD)
#define LogActivity(CPD)      CementSolutionChemistry_GetLogActivityOf(csc,CPD)
#define C0_ref                (1 * mol / Liter)
#define LogC0_ref             log10(C0_ref)



void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_H2O(CementSolutionChemistry_t* csc)
{
  double logq_ch = Input(LogQ_CH) ;
  double logq_sh = Input(LogQ_SH) ;
  double loga_na = Input(LogA_Na) ;
  double loga_k  = Input(LogA_K) ;
  double loga_oh = Input(LogA_OH) ;
  
  
  /* Hypothesis of unit water activity (to be improved) */
  double loga_h2o = 0 ;
  
  /* Autoprotolysis of water */
  double logk_h2o = LogKeq(H2O) ;
  double loga_h   = logk_h2o + loga_h2o - loga_oh ;
  
  
  /* Chemical reactions involving compounds of type I. */
  
  /* Calcium compounds */
  double loga_ca    = logq_ch - 2*(loga_oh) ;
  double logk_caoh  = LogKeq(CaOH) ;
  double loga_caoh  = loga_ca + loga_oh - logk_caoh ; ;
  double logk_caoh2 = LogKeq(CaO2H2) ;
  double loga_caoh2 = logq_ch - logk_caoh2 ;
  
  /* Silicon compounds */
  double loga_h4sio4 = logq_sh ;
  double logk_h3sio4 = LogKeq(H3SiO4) ;
  double loga_h3sio4 = loga_h4sio4 + loga_oh - (logk_h3sio4 + loga_h2o) ;
  double logk_h2sio4 = LogKeq(H2SiO4) ;
  double loga_h2sio4 = loga_h3sio4 + loga_oh - (logk_h2sio4 + loga_h2o) ;
  
  /* Sodium compounds */
  double logk_naoh = LogKeq(NaOH) ;
  double loga_naoh = loga_na + loga_oh - logk_naoh ;
  
  /* Potassium compounds */
  double logk_koh = LogKeq(KOH) ;
  double loga_koh = loga_k + loga_oh - logk_koh ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-Silicon compounds */
  double logk_cah2sio4 = LogKeq(CaH2SiO4) ;
  double loga_cah2sio4 = loga_ca + loga_h2sio4 - logk_cah2sio4 ;
  double logk_cah3sio4 = LogKeq(CaH3SiO4) ;
  double loga_cah3sio4 = loga_ca + loga_h3sio4 - logk_cah3sio4 ;
  
  
  /* Backup activities */
  {
    LogActivity(H2O)       = loga_h2o ;
    LogActivity(OH)        = loga_oh ;
    LogActivity(H)         = loga_h ;
  
    LogActivity(Ca)        = loga_ca ;
    LogActivity(CaOH)      = loga_caoh ;
    LogActivity(CaO2H2)    = loga_caoh2 ;
  
    LogActivity(H4SiO4)    = loga_h4sio4 ;
    LogActivity(H3SiO4)    = loga_h3sio4 ;
    LogActivity(H2SiO4)    = loga_h2sio4 ;
  
    LogActivity(Na)        = loga_na ;
    LogActivity(NaOH)      = loga_naoh ;
  
    LogActivity(K)         = loga_k ;
    LogActivity(KOH)       = loga_koh ;
  
    LogActivity(CaH2SiO4)  = loga_cah2sio4 ;
    LogActivity(CaH3SiO4)  = loga_cah3sio4 ;
  }
  
  
  /* Backup concentrations */
  {
    /* Translate into concentrations providing ideality */
    double logc0          = LogC0_ref ;
    double logc_oh        = loga_oh       + logc0 ;
    double logc_h         = loga_h        + logc0 ;
  
    double logc_ca        = loga_ca       + logc0 ;
    double logc_caoh      = loga_caoh     + logc0 ;
    double logc_caoh2     = loga_caoh2    + logc0 ;
  
    double logc_h4sio4    = loga_h4sio4   + logc0 ;
    double logc_h3sio4    = loga_h3sio4   + logc0 ;
    double logc_h2sio4    = loga_h2sio4   + logc0 ;
  
    double logc_na        = loga_na       + logc0 ;
    double logc_naoh      = loga_naoh     + logc0 ;
  
    double logc_k         = loga_k        + logc0 ;
    double logc_koh       = loga_koh      + logc0 ;
  
    double logc_cah2sio4  = loga_cah2sio4 + logc0 ;
    double logc_cah3sio4  = loga_cah3sio4 + logc0 ;
    
    /* 1. Concentrations */
    Concentration(OH)        = pow(10,logc_oh) ;
    Concentration(H)         = pow(10,logc_h) ;
  
    Concentration(Ca)        = pow(10,logc_ca) ;
    Concentration(CaOH)      = pow(10,logc_caoh) ;
    Concentration(CaO2H2)    = pow(10,logc_caoh2) ;
  
    Concentration(H4SiO4)    = pow(10,logc_h4sio4) ;
    Concentration(H3SiO4)    = pow(10,logc_h3sio4) ;
    Concentration(H2SiO4)    = pow(10,logc_h2sio4) ;
     
    Concentration(Na)        = pow(10,logc_na) ;
    Concentration(NaOH)      = pow(10,logc_naoh) ;
  
    Concentration(K)         = pow(10,logc_k) ;
    Concentration(KOH)       = pow(10,logc_koh) ;
  
    Concentration(CaH2SiO4)  = pow(10,logc_cah2sio4) ;
    Concentration(CaH3SiO4)  = pow(10,logc_cah3sio4) ;
  
  
    /* 2. Log10 of concentrations */
    LogConcentration(OH)        = logc_oh ;
    LogConcentration(H)         = logc_h ;
  
    LogConcentration(Ca)        = logc_ca ;
    LogConcentration(CaOH)      = logc_caoh ;
    LogConcentration(CaO2H2)    = logc_caoh2 ;
  
    LogConcentration(H4SiO4)    = logc_h4sio4 ;
    LogConcentration(H3SiO4)    = logc_h3sio4 ;
    LogConcentration(H2SiO4)    = logc_h2sio4 ;
  
    LogConcentration(Na)        = logc_na ;
    LogConcentration(NaOH)      = logc_naoh ;
  
    LogConcentration(K)         = logc_k ;
    LogConcentration(KOH)       = logc_koh ;
  
    LogConcentration(CaH2SiO4)  = logc_cah2sio4 ;
    LogConcentration(CaH3SiO4)  = logc_cah3sio4 ;
  }
}



void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2_H2O(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_H2O(csc) ;
  
  /* Supplement with CO2 */
  double loga_co2  = Input(LogA_CO2) ;
  
  
  /* Water */
  double loga_h2o = LogActivity(H2O) ;
  double loga_oh  = LogActivity(OH) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Carbon compounds */
  double logk_h2co3 = LogKeq(H2CO3) ;
  double loga_h2co3 = loga_co2 + loga_h2o - logk_h2co3 ;
  double logk_hco3  = LogKeq(HCO3) ;
  double loga_hco3  = loga_h2co3 + loga_oh - (logk_hco3 + loga_h2o) ;
  double logk_co3   = LogKeq(CO3) ;
  double loga_co3   = loga_hco3  + loga_oh - (logk_co3  + loga_h2o) ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-Carbon compounds */
  double loga_ca     = LogActivity(Ca) ;
  double logk_cahco3 = LogKeq(CaHCO3) ;
  double loga_cahco3 = loga_ca + loga_hco3 - logk_cahco3 ;
  double logk_caco3  = LogKeq(CaCO3) ;
  double loga_caco3  = loga_ca + loga_co3  - logk_caco3 ;
  
  /* Sodium-Carbon compounds */
  double loga_na     = LogActivity(Na) ;
  double logk_nahco3 = LogKeq(NaHCO3) ;
  double loga_nahco3 = loga_na + loga_hco3 - logk_nahco3 ;
  double logk_naco3  = LogKeq(NaCO3) ;
  double loga_naco3  = loga_na + loga_co3  - logk_naco3 ;
  
  
  /* Backup activities */
  {
    LogActivity(H2CO3)  = loga_h2co3 ;
    LogActivity(HCO3)   = loga_hco3 ;
    LogActivity(CO3)    = loga_co3 ;
    LogActivity(CO2)    = loga_co2 ;
  
    LogActivity(CaHCO3) = loga_cahco3 ;
    LogActivity(CaCO3)  = loga_caco3 ;
  
    LogActivity(NaHCO3) = loga_nahco3 ;
    LogActivity(NaCO3)  = loga_naco3 ;
  }
  
  
  /* Backup concentrations */
  {
    /* Translate into concentrations providing ideality */
    double logc0             = LogC0_ref ;

    double logc_h2co3        = loga_h2co3       + logc0 ;
    double logc_hco3         = loga_hco3        + logc0 ;
    double logc_co3          = loga_co3         + logc0 ;
    double logc_co2          = loga_co2         + logc0 ;
  
    double logc_cahco3       = loga_cahco3      + logc0 ;
    double logc_caco3        = loga_caco3       + logc0 ;
  
    double logc_nahco3       = loga_nahco3      + logc0 ;
    double logc_naco3        = loga_naco3       + logc0 ;
    
    Concentration(H2CO3)     = pow(10,logc_h2co3) ;
    Concentration(HCO3)      = pow(10,logc_hco3) ;
    Concentration(CO3)       = pow(10,logc_co3) ;
    Concentration(CO2)       = pow(10,logc_co2) ;
  
    Concentration(CaHCO3)    = pow(10,logc_cahco3) ;
    Concentration(CaCO3)     = pow(10,logc_caco3) ;
  
    Concentration(NaHCO3)    = pow(10,logc_nahco3) ;
    Concentration(NaCO3)     = pow(10,logc_naco3) ;
  
  
    LogConcentration(H2CO3)  = logc_h2co3 ;
    LogConcentration(HCO3)   = logc_hco3 ;
    LogConcentration(CO3)    = logc_co3 ;
    LogConcentration(CO2)    = logc_co2 ;
  
    LogConcentration(CaHCO3) = logc_cahco3 ;
    LogConcentration(CaCO3)  = logc_caco3 ;
  
    LogConcentration(NaHCO3) = logc_nahco3 ;
    LogConcentration(NaCO3)  = logc_naco3 ;
  }
}






void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O(CementSolutionChemistry_t* csc)
{
  if(CementSolutionChemistry_InputIs(csc,SO3,LogA_H2SO4)) {
    CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_1(csc) ;
    return ;
  } else if(CementSolutionChemistry_InputIs(csc,SO3,LogA_SO4)) {
    CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_2(csc) ;
    return ;
  }
  
  arret("CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O") ;
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_1(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_H2O(csc) ;
  
  /* Supplement with SO3 */
  double loga_h2so4 = Input(LogA_H2SO4) ;
  

  /* Water */
  double loga_h  = LogActivity(H) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Sulfur compounds */
  double logk_h2so4 = LogKeq(H2SO4) ;
  double loga_hso4  = loga_h2so4 - loga_h + logk_h2so4 ;
  double logk_hso4  = LogKeq(HSO4) ;
  double loga_so4   = loga_hso4 - loga_h + logk_hso4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-sulfur compounds */
  double loga_ca     = LogActivity(Ca) ;
  double logk_cahso4 = LogKeq(CaHSO4) ;
  double loga_cahso4 = loga_ca + loga_hso4 - logk_cahso4 ;
  double logk_caso4  = LogKeq(CaSO4) ;
  double loga_caso4  = loga_ca + loga_so4  - logk_caso4 ;
  
  /* Sodium-sulfur compounds */
  
  
  /* Backup activities */
  {
    LogActivity(H2SO4)  = loga_h2so4 ;
    LogActivity(HSO4)   = loga_hso4 ;
    LogActivity(SO4)    = loga_so4 ;
  
    LogActivity(CaHSO4) = loga_cahso4 ;
    LogActivity(CaSO4)  = loga_caso4 ;
  }
  
  /* Backup concentrations */
  {
    /* Translate into concentrations providing ideality */
    double logc0             = LogC0_ref ;
    
    double logc_h2so4        = loga_h2so4  + logc0 ;
    double logc_hso4         = loga_hso4   + logc0 ;
    double logc_so4          = loga_so4    + logc0 ;
   
    double logc_cahso4       = loga_cahso4 + logc0 ;
    double logc_caso4        = loga_caso4  + logc0 ;
    
    Concentration(H2SO4)     = pow(10,logc_h2so4) ;
    Concentration(HSO4)      = pow(10,logc_hso4) ;
    Concentration(SO4)       = pow(10,logc_so4) ;
  
    Concentration(CaHSO4)    = pow(10,logc_cahso4) ;
    Concentration(CaSO4)     = pow(10,logc_caso4) ;
  
  
    LogConcentration(H2SO4)  = logc_h2so4 ;
    LogConcentration(HSO4)   = logc_hso4 ;
    LogConcentration(SO4)    = logc_so4 ;

    LogConcentration(CaHSO4) = logc_cahso4 ;
    LogConcentration(CaSO4)  = logc_caso4 ;
  }
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_2(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_H2O(csc) ;
  
  /* Supplement with SO3 */
  double loga_so4 = Input(LogA_SO4) ;
  

  /* Water */
  double loga_h  = LogActivity(H) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Sulfur compounds */
  double logk_h2so4 = LogKeq(H2SO4) ;
  double logk_hso4  = LogKeq(HSO4) ;
  double loga_hso4  = loga_so4  + loga_h - logk_hso4 ;
  double loga_h2so4 = loga_hso4 + loga_h - logk_h2so4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-sulfur compounds */
  double loga_ca     = LogActivity(Ca) ;
  double logk_cahso4 = LogKeq(CaHSO4) ;
  double loga_cahso4 = loga_ca + loga_hso4 - logk_cahso4 ;
  double logk_caso4  = LogKeq(CaSO4) ;
  double loga_caso4  = loga_ca + loga_so4  - logk_caso4 ;
  
  /* Sodium-sulfur compounds */
  
  
  /* Backup activities */
  {
    LogActivity(H2SO4)  = loga_h2so4 ;
    LogActivity(HSO4)   = loga_hso4 ;
    LogActivity(SO4)    = loga_so4 ;
  
    LogActivity(CaHSO4) = loga_cahso4 ;
    LogActivity(CaSO4)  = loga_caso4 ;
  }
  
  /* Backup concentrations */
  {
    /* Translate into concentrations providing ideality */
    double logc0             = LogC0_ref ;
    
    double logc_h2so4        = loga_h2so4  + logc0 ;
    double logc_hso4         = loga_hso4   + logc0 ;
    double logc_so4          = loga_so4    + logc0 ;
   
    double logc_cahso4       = loga_cahso4 + logc0 ;
    double logc_caso4        = loga_caso4  + logc0 ;
    
    Concentration(H2SO4)     = pow(10,logc_h2so4) ;
    Concentration(HSO4)      = pow(10,logc_hso4) ;
    Concentration(SO4)       = pow(10,logc_so4) ;
  
    Concentration(CaHSO4)    = pow(10,logc_cahso4) ;
    Concentration(CaSO4)     = pow(10,logc_caso4) ;
  
  
    LogConcentration(H2SO4)  = logc_h2so4 ;
    LogConcentration(HSO4)   = logc_hso4 ;
    LogConcentration(SO4)    = logc_so4 ;

    LogConcentration(CaHSO4) = logc_cahso4 ;
    LogConcentration(CaSO4)  = logc_caso4 ;
  }
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O-SO3 */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O(csc) ;
  
  /* Supplement with Al2O3 */
  double logq_ah3   = Input(LogQ_AH3) ;
  

  /* Water */
  double loga_oh  = LogActivity(OH) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Aluminium compounds */
  double loga_al     = 0.5*logq_ah3 - 3*loga_oh ;
  double logk_alo4h4 = LogKeq(AlO4H4) ;
  double loga_alo4h4 = loga_al + 4*loga_oh - logk_alo4h4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  
  /* Backup activities */
  {
    LogActivity(Al)     = loga_al ;
    LogActivity(AlO4H4) = loga_alo4h4 ;
  }
  
  
  /* Backup concentrations */
  {
    /* Translate into concentrations providing ideality */
    double logc0             = LogC0_ref ;
    
    double logc_al           = loga_al      + logc0 ;
    double logc_alo4h4       = loga_alo4h4  + logc0 ;
    
    Concentration(Al)        = pow(10,logc_al) ;
    Concentration(AlO4H4)    = pow(10,logc_alo4h4) ;
  
    LogConcentration(Al)     = logc_al ;
    LogConcentration(AlO4H4) = logc_alo4h4 ;
  }
}



void CementSolutionChemistry_CopyConcentrations(CementSolutionChemistry_t* csc,double* v)
{
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* c = CementSolutionChemistry_GetConcentration(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = c[i] ;
      }
    }
  }
}



void CementSolutionChemistry_CopyLogConcentrations(CementSolutionChemistry_t* csc,double* v)
{
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* logc = CementSolutionChemistry_GetLogConcentration(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = logc[i] ;
      }
    }
  }
}



void CementSolutionChemistry_CopyChemicalPotential(CementSolutionChemistry_t* csc,double* v)
{
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* logc = CementSolutionChemistry_GetLogConcentration(csc) ;
    double  epot = CementSolutionChemistry_GetElectricPotential(csc) ;
    double* z = CementSolutionChemistry_GetValence() ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = Ln10 * logc[i] + z[i] * epot ;
      }
    }
  }
}



void CementSolutionChemistry_UpdateSolution(CementSolutionChemistry_t* csc)
{
  /* Charge density */
  {
    double q = CementSolutionChemistry_ComputeChargeDensity(csc) ;
    
    CementSolutionChemistry_GetChargeDensity(csc) = q ;
  }
  
  /* Mass density */
  {
    double rho_l = CementSolutionChemistry_ComputeLiquidMassDensity(csc) ;
    
    CementSolutionChemistry_GetLiquidMassDensity(csc) = rho_l ;
  }
  
  /* Element concentrations */
  CementSolutionChemistry_UpdateElementConcentrations(csc) ;
}



/* Shorthands of macros */
#include "ElectricChargeOfIonInWater.h"
#define Z(A)          ElectricChargeOfIonInWater(A)



double CementSolutionChemistry_SolveElectroneutrality(CementSolutionChemistry_t* csc)
/** Solve the electroneutrality equation, SUM(z_i c_i) = 0,
 ** for c_h or c_oh, as root of a 4th order polynomial:
 ** ax^4 + bx^3 + cx^2 + dx + e = 0
 **/
{
  /* The primary variables are considered as constant:
   * q_ch, q_sh, q_ah3, c_na, c_k, c_co2, c_h2so4, c_cl.
   */
  
  /* Electroneutrality is written as  Sum z_i*c_i = 0, with
   * 
   * c_i = A_i * (c_h)**(n_i)  = c_i0 * (c_h/c_h0)**(n_i)
   * or equivalently
   * c_i = B_i * (c_oh)**(m_i) = c_i0 * (c_oh/c_oh0)**(m_i)
   * with 
   * m_i = -n_i and (c_h*c_oh) = (c_h0*c_oh0) = K_w = cst
   * 
   * c_i         (prop to)  (c_h)**(n_i)           : n_i  z_i
   * --------------------------------------------------------
   * c_h         (prop to)  c_h                    : +1   +1
   * c_oh        (prop to)  1 / c_h                : -1   -1
   * 
   * c_ca        (prop to)  q_ch / (c_oh)^2        : +2   +2
   * c_caoh      (prop to)  c_ca * c_oh            : +1   +1
   * c_caoh2     (prop to)  q_ch                   :  0    0
   * 
   * c_h4sio4    (prop to)  q_sh                   :  0    0
   * c_h3sio4    (prop to)  c_h4sio4 / c_h         : -1   -1
   * c_h2sio4    (prop to)  c_h3sio4 / c_h         : -2   -2
   * 
   * c_na        (prop to)  c_na                   :  0   +1
   * c_naoh      (prop to)  c_na * c_oh            : -1    0
   * 
   * c_k         (prop to)  c_k                    :  0   +1
   * c_koh       (prop to)  c_k * c_h              : -1    0
   * 
   * c_cl        (prop to)  c_cl                   :  0   -1
   * 
   * c_co2       (prop to)  c_co2                  :  0    0
   * c_h2co3     (prop to)  c_co2                  :  0    0
   * c_hco3      (prop to)  c_h2co3 / c_h          : -1   -1
   * c_co3       (prop to)  c_hco3 / c_h           : -2   -2
   * 
   * c_h2so4     (prop to)  c_h2so4                :  0    0
   * c_hso4      (prop to)  c_h2so4 / c_h          : -1   -1
   * c_so4       (prop to)  c_hso4 / c_h           : -2   -2
   * 
   * c_al        (prop to)  (q_ah3)^0.5 / (c_oh)^3 : +3   +3
   * c_alo4h4    (prop to)  c_al * (c_oh)^4        : -1   -1
   * 
   * c_cah2sio4  (prop to)  c_h2sio4 * c_ca        :  0    0
   * c_cah3sio4  (prop to)  c_h3sio4 * c_ca        : +1   +1
   * 
   * c_cahco3    (prop to)  c_hco3 * c_ca          : +1   +1
   * c_caco3     (prop to)  c_ca * c_co3           :  0    0
   * 
   * c_nahco3    (prop to)  c_na * c_hco3          : -1    0
   * c_naco3     (prop to)  c_na * c_co3           : -2   -1
   * 
   * c_caso4     (prop to)  c_ca * c_so4           :  0    0
   * c_cahso4    (prop to)  c_ca * c_hso4          : +1   +1
   */
   
  
  double x = 1. ;
   
  /* Load the c_i0 for i having non zero z_i */
  {
    double c_h        = Concentration(H) ;
    double c_oh       = Concentration(OH) ;
  
    double c_ca       = Concentration(Ca) ;
    double c_caoh     = Concentration(CaOH) ;
  
    double c_h3sio4   = Concentration(H3SiO4) ;
    double c_h2sio4   = Concentration(H2SiO4) ;
  
    double c_na       = Concentration(Na) ;
  
    double c_k        = Concentration(K) ;
  
    double c_hco3     = Concentration(HCO3) ;
    double c_co3      = Concentration(CO3) ;
  
    double c_hso4     = Concentration(HSO4) ;
    double c_so4      = Concentration(SO4) ;
  
    double c_al       = Concentration(Al) ;
    double c_alo4h4   = Concentration(AlO4H4) ;
  
    double c_cl       = Concentration(Cl) ;
  
    double c_cah3sio4 = Concentration(CaH3SiO4) ;
  
    double c_cahco3   = Concentration(CaHCO3) ;
  
    double c_naco3    = Concentration(NaCO3) ;
  
    double c_cahso4   = Concentration(CaHSO4) ;

    /* Compute the positive charge in a3,a2,a1,a0 for n = +3,+2,+1,0 */
    double a3 = Z(Al)*c_al ;
    double a2 = Z(Ca)*c_ca ;
    double a1 = Z(H)*c_h + Z(CaHCO3)*c_cahco3 + Z(CaH3SiO4)*c_cah3sio4 + Z(CaOH)*c_caoh + Z(CaHSO4)*c_cahso4 ;
    double a0 = Z(Na)*c_na + Z(K)*c_k ;
    /* Compute the negative charge in b0,b1,b2 for n = 0,-1,-2 */
    double b0 = Z(Cl)*c_cl ;
    double b1 = Z(OH)*c_oh + Z(HCO3)*c_hco3 + Z(H3SiO4)*c_h3sio4 + Z(HSO4)*c_hso4 + Z(AlO4H4)*c_alo4h4 ;
    double b2 = Z(CO3)*c_co3 + Z(H2SiO4)*c_h2sio4 + Z(NaCO3)*c_naco3 + Z(SO4)*c_so4 ;
  
    if(a3 == 0) {
      double loga_oh = LogActivity(OH) ;
      double loga_h  = LogActivity(H) ;
      double a_h     = pow(10,loga_h) ;
      double a_oh    = pow(10,loga_oh) ;
      
      /* Solve for x = (c_h/c_h0) as root of the 4th order polynomial */
      x = poly4(a2,a1,a0+b0,b1,b2,a_h,a_oh) ;
    
    } else {
      double loga_oh = LogActivity(OH) ;
      double loga_h  = LogActivity(H) ;
      double a_h     = pow(10,loga_h) ;
      double a_oh    = pow(10,loga_oh) ;
      double y[6] = {a3,a2,a1,a0,b1,b2} ;
      double tol = 1.e-4 ;
  
      x = poly4(a2,a1,a0+b0,b1,b2,a_h,a_oh) ;
      x = Math_PolishPolynomialEquationRoot(y,5,x,tol*x,20) ;
    }
  
  
    if(x < 0) {
      double y = 1/x ;
    
      printf("c_h   = %e\n",c_h*x) ;
      printf("c_oh  = %e\n",c_oh*y) ;
      printf("a3    = %e\n",a3) ;
      printf("a2    = %e\n",a2) ;
      printf("a1    = %e\n",a1) ;
      printf("a0    = %e\n",a0) ;
      printf("b1    = %e\n",b1) ;
      printf("b2    = %e\n",b2) ;
      /* Raise an interrupt signal instead of exit */
      Message_Warning("CementSolutionChemistry_SolveElectroneutrality: c_h < 0") ;
      Exception_Interrupt ;
      //arret("CementSolutionChemistry_SolveElectroneutrality: c_h<0") ;
    }
  }
  
  
  /* Update the c_i for i having non zero n_i
   * We multiply the c_i by x**n_i or by y**m_i */
  {
    double x2 = x*x ;
    double x3 = x2*x ;
    double y = 1/x ;
    double y2 = y*y ;
  
    Concentration(H)        *= x ;
    Concentration(OH)       *= y ;
  
    Concentration(Ca)       *= x2 ;
    Concentration(CaOH)     *= x ;
    //Concentration(CaO2H2)   *= 1 ;
  
    //Concentration(H4SiO4)   *= 1 ;
    Concentration(H3SiO4)   *= y ;
    Concentration(H2SiO4)   *= y2 ;
  
    //Concentration(Na)       *= 1 ;
    Concentration(NaOH)     *= y ;
  
    //Concentration(K)        *= 1 ;
    Concentration(KOH)      *= y ;
  
    //Concentration(CO2)      *= 1 ;
    //Concentration(H2CO3)    *= 1 ;
    Concentration(HCO3)     *= y ;
    Concentration(CO3)      *= y2 ;
  
    //Concentration(H2SO4)    *= 1 ;
    Concentration(HSO4)     *= y ;
    Concentration(SO4)      *= y2 ;
  
    Concentration(Al)       *= x3 ;
    Concentration(AlO4H4)   *= y ;
  
    //Concentration(CaH2SiO4) *= 1 ;
    Concentration(CaH3SiO4) *= x ;
    
    Concentration(CaHCO3)   *= x ;
    //Concentration(CaCO3)    *= 1 ;
  
    Concentration(NaCO3)    *= y2 ;
    Concentration(NaHCO3)   *= y ;
  
    //Concentration(CaSO4)    *= 1 ;
    Concentration(CaHSO4)   *= x ;
  }
  
  
  /* Update the logc_i for i having non zero n_i
   * We multiply the c_i by x**n_i or by y**m_i */
  {
    double logx = log10(x) ;
    double logx2 = 2*logx ;
    double logx3 = 3*logx ;
    double logy = -logx ;
    double logy2 = 2*logy ;
  
    LogConcentration(H)        += logx ;
    LogConcentration(OH)       += logy ;
  
    LogConcentration(Ca)       += logx2 ;
    LogConcentration(CaOH)     += logx ;
    //LogConcentration(CaO2H2)   += 0 ;
  
    //LogConcentration(H4SiO4)   += 0 ;
    LogConcentration(H3SiO4)   += logy ;
    LogConcentration(H2SiO4)   += logy2 ;
  
    //LogConcentration(Na)       += 0 ;
    LogConcentration(NaOH)     += logy ;
  
    //LogConcentration(K)        += 0 ;
    LogConcentration(KOH)      += logy ;
  
    //LogConcentration(CO2)      += 0 ;
    //LogConcentration(H2CO3)    += 0 ;
    LogConcentration(HCO3)     += logy ;
    LogConcentration(CO3)      += logy2 ;
  
    //LogConcentration(H2SO4)    += 0 ;
    LogConcentration(HSO4)     += logy ;
    LogConcentration(SO4)      += logy2 ;
  
    LogConcentration(Al)       += logx3 ;
    LogConcentration(AlO4H4)   += logy ;
  
    //LogConcentration(CaH2SiO4) += 0 ;
    LogConcentration(CaH3SiO4) += logx ;
    
    LogConcentration(CaHCO3)   += logx ;
    //LogConcentration(CaCO3)    += 0 ;
  
    LogConcentration(NaCO3)    += logy2 ;
    LogConcentration(NaHCO3)   += logy ;
  
    //LogConcentration(CaSO4)    += 0 ;
    LogConcentration(CaHSO4)   += logx ;
  }
  
  CementSolutionChemistry_TranslateConcentrationsIntoActivities(csc) ;
  
  CementSolutionChemistry_UpdateSolution(csc) ;
  
  return(CementSolutionChemistry_GetChargeDensity(csc)) ;
}





double CementSolutionChemistry_SolveExplicitElectroneutrality(CementSolutionChemistry_t* csc)
/** Solve the electroneutrality equation, SUM(z_i c_i) = 0,
 ** for c_h or c_oh, as root of a 2th order polynomial:
 ** ax^2 + bx + c = 0, keeping constant all other ion concentrations.
 **/
{
  double c_h  = Concentration(H) ;
  double c_oh = Concentration(OH) ;
  double q    = CementSolutionChemistry_GetChargeDensity(csc) ;
  double q0   = 0.5 * (q - c_h + c_oh) ;
  /* solve 2*q0 + c_h - c_oh = 0 */
  double kw   = c_h*c_oh ;
  
  if(q0 > 0) {
    c_oh =   q0 + sqrt(q0*q0 + kw) ;
    c_h  = kw/c_oh ;
  } else {
    c_h  = - q0 + sqrt(q0*q0 + kw) ;
    c_oh = kw/c_h ;
  }
  
  Concentration(H)  = c_h ;
  Concentration(OH) = c_oh ;
  
  LogConcentration(H)  = log10(c_h) ;
  LogConcentration(OH) = log10(c_oh) ;
    
  
  return(CementSolutionChemistry_ComputeChargeDensity(csc)) ;
}




/* Intern functions */

double CementSolutionChemistry_ComputeChargeDensity(CementSolutionChemistry_t* csc)
/** Return the charge **/
{
  double c_h        = Concentration(H) ;
  double c_oh       = Concentration(OH) ;
  
  double c_ca       = Concentration(Ca) ;
  double c_caoh     = Concentration(CaOH) ;
  
  double c_h3sio4   = Concentration(H3SiO4) ;
  double c_h2sio4   = Concentration(H2SiO4) ;
  
  double c_na       = Concentration(Na) ;
  
  double c_k        = Concentration(K) ;
  
  double c_cl       = Concentration(Cl) ;
  
  double c_hco3     = Concentration(HCO3) ;
  double c_co3      = Concentration(CO3) ;
  
  double c_hso4     = Concentration(HSO4) ;
  double c_so4      = Concentration(SO4) ;
  
  double c_al       = Concentration(Al) ;
  double c_alo4h4   = Concentration(AlO4H4) ;
  
  double c_cahso4   = Concentration(CaHSO4) ;
  
  double c_cah3sio4 = Concentration(CaH3SiO4) ;
  
  double c_naco3    = Concentration(NaCO3) ;
  
  double c_cahco3   = Concentration(CaHCO3) ;


  /* The charge */
  double q_0     = Z(H)*c_h + Z(OH)*c_oh ;
  double q_ca    = Z(Ca)*c_ca + Z(CaOH)*c_caoh ;
  double q_si    = Z(H2SiO4)*c_h2sio4 + Z(H3SiO4)*c_h3sio4 ;
  double q_na    = Z(Na)*c_na ;
  double q_k     = Z(K)*c_k ;
  double q_cl    = Z(Cl)*c_cl ;
  double q_c     = Z(CO3)*c_co3 + Z(HCO3)*c_hco3 ;
  double q_s     = Z(SO4)*c_so4 + Z(HSO4)*c_hso4 ;
  double q_al    = Z(Al)*c_al + Z(AlO4H4)*c_alo4h4 ;
  double q_ca_si = Z(CaH3SiO4)*c_cah3sio4 ;
  double q_ca_c  = Z(CaHCO3)*c_cahco3 ;
  double q_na_c  = Z(NaCO3)*c_naco3 ;
  double q_ca_s  = Z(CaHSO4)*c_cahso4 ;
  
  double q = q_0 + q_ca + q_si + q_na + q_k + q_cl + q_c + q_s + q_al \
           + q_ca_si + q_ca_c + q_na_c + q_ca_s ;
  
  CementSolutionChemistry_GetChargeDensity(csc) = q ;
           
  
  /* Ionic strength */
  double TwoI_0     = Z(H)*Z(H)*c_h + Z(OH)*Z(OH)*c_oh ;
  double TwoI_ca    = Z(Ca)*Z(Ca)*c_ca + Z(CaOH)*Z(CaOH)*c_caoh ;
  double TwoI_si    = Z(H2SiO4)*Z(H2SiO4)*c_h2sio4 + Z(H3SiO4)*Z(H3SiO4)*c_h3sio4 ;
  double TwoI_na    = Z(Na)*Z(Na)*c_na ;
  double TwoI_k     = Z(K)*Z(K)*c_k ;
  double TwoI_cl    = Z(Cl)*Z(Cl)*c_cl ;
  double TwoI_c     = Z(CO3)*Z(CO3)*c_co3 + Z(HCO3)*Z(HCO3)*c_hco3 ;
  double TwoI_s     = Z(SO4)*Z(SO4)*c_so4 + Z(HSO4)*Z(HSO4)*c_hso4 ;
  double TwoI_al    = Z(Al)*Z(Al)*c_al + Z(AlO4H4)*Z(AlO4H4)*c_alo4h4 ;
  double TwoI_ca_si = Z(CaH3SiO4)*Z(CaH3SiO4)*c_cah3sio4 ;
  double TwoI_ca_c  = Z(CaHCO3)*Z(CaHCO3)*c_cahco3 ;
  double TwoI_na_c  = Z(NaCO3)*Z(NaCO3)*c_naco3 ;
  double TwoI_ca_s  = Z(CaHSO4)*Z(CaHSO4)*c_cahso4 ;
  
  double TwoI       = TwoI_0 + TwoI_ca + TwoI_si + TwoI_na + TwoI_k + TwoI_cl + TwoI_c + TwoI_s + TwoI_al \
                    + TwoI_ca_si + TwoI_ca_c + TwoI_na_c + TwoI_ca_s ;
  
  CementSolutionChemistry_GetIonicStrength(csc) = 0.5*TwoI ;
          
  return(q) ;
}

#undef Z





/* Shorthands of macros */
#include "PartialMolarVolumeOfMoleculeInWater.h"
#define V(A)           (PartialMolarVolumeOfMoleculeInWater(A))
#include "MolarMassOfMolecule.h"
#define M(A)           (MolarMassOfMolecule(A))



double CementSolutionChemistry_ComputeLiquidMassDensity(CementSolutionChemistry_t* csc)
/** Return the liquid mass density **/
{
  double c_h        = Concentration(H) ;
  double c_oh       = Concentration(OH) ;
  
  double c_ca       = Concentration(Ca) ;
  double c_caoh     = Concentration(CaOH) ;
  double c_caoh2    = Concentration(CaO2H2) ;
  
  double c_h4sio4   = Concentration(H4SiO4) ;
  double c_h3sio4   = Concentration(H3SiO4) ;
  double c_h2sio4   = Concentration(H2SiO4) ;
  
  double c_na       = Concentration(Na) ;
  double c_naoh     = Concentration(NaOH) ;
  
  double c_k        = Concentration(K) ;
  double c_koh      = Concentration(KOH) ;
  
  double c_cl       = Concentration(Cl) ;
  
  double c_co2      = Concentration(CO2) ;
  double c_h2co3    = Concentration(H2CO3) ;
  double c_hco3     = Concentration(HCO3) ;
  double c_co3      = Concentration(CO3) ;
  
  double c_h2so4    = Concentration(H2SO4) ;
  double c_hso4     = Concentration(HSO4) ;
  double c_so4      = Concentration(SO4) ;
  
  double c_al       = Concentration(Al) ;
  double c_alo4h4   = Concentration(AlO4H4) ;
  
  double c_caso4    = Concentration(CaSO4) ;
  double c_cahso4   = Concentration(CaHSO4) ;
  
  double c_cah2sio4 = Concentration(CaH2SiO4) ;
  double c_cah3sio4 = Concentration(CaH3SiO4) ;
  
  double c_naco3    = Concentration(NaCO3) ;
  double c_nahco3   = Concentration(NaHCO3) ;
  
  double c_caco3    = Concentration(CaCO3) ;
  double c_cahco3   = Concentration(CaHCO3) ;
  
  
  /* Partial molar volumes */
  double v_h      = c_h*V(H) + c_oh*V(OH) ;
  double v_ca     = c_ca*V(Ca) + c_caoh*V(CaOH) + c_caoh2*V(CaO2H2) ;
  double v_si     = c_h3sio4*V(H3SiO4) + c_h4sio4*V(H4SiO4) + c_h2sio4*V(H2SiO4) ;
  double v_na     = c_na*V(Na) + c_naoh*V(NaOH) ;
  double v_k      = c_k*V(K) + c_koh*V(KOH) ;
  double v_cl     = c_cl*V(Cl) ;
  double v_c      = c_co2*V(CO2) + c_h2co3*V(H2CO3) + c_hco3*V(HCO3) + c_co3*V(CO3) ;
  double v_s      = c_h2so4*V(H2SO4) + c_hso4*V(HSO4) + c_so4*V(SO4) ;
  double v_al     = c_al*V(Al) + c_alo4h4*V(AlO4H4) ;
  double v_ca_si  = c_cah2sio4*V(CaH2SiO4) + c_cah3sio4*V(CaH3SiO4) ;
  double v_ca_c   = c_cahco3*V(CaHCO3) + c_caco3*V(CaCO3) ;
  double v_na_c   = c_nahco3*V(NaHCO3) + c_naco3*V(NaCO3) ;
  double v_ca_s   = c_cahso4*V(CaHSO4) + c_caso4*V(CaSO4) ;
  
  double v_all    = v_h + v_ca + v_si + v_na + v_k + v_c + v_cl + v_s + v_al \
                  + v_ca_si + v_ca_c + v_na_c + v_ca_s ;
  
  /* Water */
  double c_h2o    = (1 - v_all)/V(H2O) ;
       
  /* Update */
  Concentration(H2O)  = c_h2o ;
         
  /* Solution */
  double rho_0     = M(H)*c_h + M(OH)*c_oh + M(H2O)*c_h2o ;
  double rho_ca    = M(Ca)*c_ca + M(CaOH)*c_caoh + M(CaO2H2)*c_caoh2 ;
  double rho_si    = M(H3SiO4)*c_h3sio4 + M(H4SiO4)*c_h4sio4 + M(H2SiO4)*c_h2sio4 ;
  double rho_na    = M(Na)*c_na + M(NaOH)*c_naoh ;
  double rho_k     = M(K)*c_k + M(KOH)*c_koh ;
  double rho_cl    = M(Cl)*c_cl ;
  double rho_c     = M(CO2)*c_co2 + M(H2CO3)*c_h2co3 + M(HCO3)*c_hco3 + M(CO3)*c_co3 ;
  double rho_s     = M(H2SO4)*c_h2so4 + M(HSO4)*c_hso4 + M(SO4)*c_so4 ;
  double rho_al    = M(Al)*c_al + M(AlO4H4)*c_alo4h4 ;
  double rho_ca_si = M(CaH2SiO4)*c_cah2sio4 + M(CaH3SiO4)*c_cah3sio4 ;
  double rho_ca_c  = M(CaHCO3)*c_cahco3 + M(CaCO3)*c_caco3 ;
  double rho_na_c  = M(NaHCO3)*c_nahco3 + M(NaCO3)*c_naco3 ;
  double rho_ca_s  = M(CaHSO4)*c_cahso4 + M(CaSO4)*c_caso4 ;
  
  double rho_l = rho_0 + rho_ca + rho_si + rho_na + rho_k + rho_cl + rho_c + rho_s + rho_al \
               + rho_ca_si + rho_ca_c + rho_na_c + rho_ca_s ;
       
  /* Backup */
  CementSolutionChemistry_GetLiquidMassDensity(csc)  = rho_l ;


  return(rho_l) ;
}

#undef V
#undef M



void CementSolutionChemistry_Initialize(CementSolutionChemistry_t* csc)
{
  /* Concentrations */
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* v = CementSolutionChemistry_GetConcentration(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
  
  /* Log of concentrations */
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* v = CementSolutionChemistry_GetLogConcentration(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = - 99 ;
      }
    }
  }
  
  /* Log of activities */
  {
    int     n = CementSolutionChemistry_NbOfSpecies ;
    double* v = CementSolutionChemistry_GetLogActivity(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = - 99 ;
      }
    }
  }
  
  /* Element concentrations */
  {
    int     n = CementSolutionChemistry_NbOfElementConcentrations ;
    double* v = CementSolutionChemistry_GetElementConcentration(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
  
  /* Other variables */
  {
    int     n = CementSolutionChemistry_NbOfOtherVariables ;
    double* v = CementSolutionChemistry_GetOtherVariable(csc) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
}





void CementSolutionChemistry_UpdateElementConcentrations(CementSolutionChemistry_t* csc)
/** Update the element concentrations **/
{
  //double c_h        = Concentration(H) ;
  //double c_oh       = Concentration(OH) ;
  
  /* Calcium */
  double c_ca       = Concentration(Ca) ;
  double c_caoh     = Concentration(CaOH) ;
  double c_caoh2    = Concentration(CaO2H2) ;
  
  /* Silicon */
  double c_h4sio4   = Concentration(H4SiO4) ;
  double c_h3sio4   = Concentration(H3SiO4) ;
  double c_h2sio4   = Concentration(H2SiO4) ;
  
  /* Sodium */
  double c_na       = Concentration(Na) ;
  double c_naoh     = Concentration(NaOH) ;
  
  /* Potassium */
  double c_k        = Concentration(K) ;
  double c_koh      = Concentration(KOH) ;
  
  /* Carbon */
  double c_co2      = Concentration(CO2) ;
  double c_h2co3    = Concentration(H2CO3) ;
  double c_hco3     = Concentration(HCO3) ;
  double c_co3      = Concentration(CO3) ;
  
  /* Sulfur */
  double c_h2so4    = Concentration(H2SO4) ;
  double c_hso4     = Concentration(HSO4) ;
  double c_so4      = Concentration(SO4) ;
  
  /* Aluminium */
  double c_al       = Concentration(Al) ;
  double c_alo4h4   = Concentration(AlO4H4) ;
  
  /* Chlore */
  double c_cl       = Concentration(Cl) ;
  
  /* Compounds of type II */
  double c_cah2sio4 = Concentration(CaH2SiO4) ;
  double c_cah3sio4 = Concentration(CaH3SiO4) ;
  
  double c_naco3    = Concentration(NaCO3) ;
  double c_nahco3   = Concentration(NaHCO3) ;
  
  double c_caco3    = Concentration(CaCO3) ;
  double c_cahco3   = Concentration(CaHCO3) ;
  
  double c_caso4    = Concentration(CaSO4) ;
  double c_cahso4   = Concentration(CaHSO4) ;

  
  /* Concentration as element: C, Ca, Si ... */
  /* Compounds type I */
  double c_ca_l    = c_ca + c_caoh + c_caoh2 ;
  double c_si_l    = c_h2sio4 + c_h3sio4 + c_h4sio4 ;
  double c_na_l    = c_na + c_naoh ;
  double c_k_l     = c_k + c_koh ;
  double c_c_l     = c_co2 + c_h2co3 + c_hco3 + c_co3 ;
  double c_s_l     = c_h2so4 + c_hso4 + c_so4 ;
  double c_al_l    = c_al + c_alo4h4 ;
  double c_cl_l    = c_cl ;
  /* Compounds type II */
  double c_ca_si_l = c_cah2sio4 + c_cah3sio4 ;
  double c_ca_c_l  = c_cahco3 + c_caco3 ;
  double c_ca_s_l  = c_caso4 + c_cahso4 ;
  double c_na_c_l  = c_nahco3 + c_naco3 ;
  
  
  /* Backup */
  
  CementSolutionChemistry_GetElementConcentrationOf(csc,Ca) = c_ca_l + c_ca_si_l + c_ca_c_l + c_ca_s_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,Si) = c_si_l + c_ca_si_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,Na) = c_na_l + c_na_c_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,K)  = c_k_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,C)  = c_c_l + c_ca_c_l + c_na_c_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,S)  = c_s_l + c_ca_s_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,Al) = c_al_l ;
  CementSolutionChemistry_GetElementConcentrationOf(csc,Cl) = c_cl_l ;
}



void CementSolutionChemistry_TranslateConcentrationsIntoActivities(CementSolutionChemistry_t* csc)
{
  double logc0 = LogC0_ref ;
  int     n = CementSolutionChemistry_NbOfSpecies ;
  double* logc = CementSolutionChemistry_GetLogConcentration(csc) ;
  double* loga = CementSolutionChemistry_GetLogActivity(csc) ;
  
  {
    int i ;
      
    for(i = 0 ; i < n ; i++) {
      loga[i] = logc[i] - logc0 ;
    }
  }
}




double poly4(double a,double b,double c,double d,double e,double a_h,double a_oh)
/* Solve ax^4 + bx^3 + cx^2 + dx + e = 0 
 * for x in the range defined by x*a_h < 1 and a_oh/x < 1 (a_oh < x < 1/a_h)
 * because a_h and a_oh should be < 1 (ie 0 < pH < -logKw) 
 * The new values of a_h and a_oh are: x*a_h and a_oh/x.
 * Return the solution which is the closest to 1. */
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
    
    x = y[0] ;
    for(i = 1 ; i < n ; i++) {
      double x1 = fabs(x - 1) ;
      double y1 = fabs(y[i] - 1) ;
      
      if(y1 < x1) x = y[i] ;
    }

    /* Too constraining
    if((x*a_h > 1) || (a_oh/x > 1)) {
      printf("\n") ;
      printf("n    = %d\n",n) ;
      printf("a,b,c,d,e = %e,%e,%e,%e,%e\n",a,b,c,d,e) ;
      printf("x    = %e\n",x) ;
      printf("a_h  = %e\n",a_h) ;
      printf("a_oh = %e\n",a_oh) ;
      // Raise an interrupt signal instead of exit
      Message_Warning("poly4: a_h = %e > 1 or a_oh = %e > 1!",x*a_h,a_oh/x) ;
      Exception_Interrupt ;
    }
    */
  }
  
  y[0] = a ;
  y[1] = b ;
  y[2] = c ;
  y[3] = d ;
  y[4] = e ;
  
  x = Math_PolishPolynomialEquationRoot(y,4,x,tol*x,20) ;
  
  return(x) ;
}



#if 0
double poly4(double a,double b,double c,double d,double e,double a_h,double a_oh)
/* Solve ax^4 + bx^3 + cx^2 + dx + e = 0 
 * for x in the range defined by x*a_h < 1 and a_oh/x < 1 (a_oh < x < 1/a_h)
 * because a_h and a_oh should be < 1 (ie 0 < pH < -logKw) */
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
      x = y[i] ;
      if((x*a_h < 1) && (x > a_oh)) break ;
      //if((x*a_h < 1)) break ;
    }
    
    if(i == n) {
      printf("\n") ;
      printf("n    = %d\n",n) ;
      printf("a,b,c,d,e = %e,%e,%e,%e,%e\n",a,b,c,d,e) ;
      printf("x    = %e\n",x) ;
      printf("a_h  = %e\n",a_h) ;
      printf("a_oh = %e\n",a_oh) ;
      /* Raise an interrupt signal instead of exit */
      Message_Warning("poly4: a_h = %e > 1 or a_oh = %e > 1!",x*a_h,a_oh/x) ;
      Exception_Interrupt ;
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
#endif



#if 0
void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_1(CementSolutionChemistry_t* csc)
{
  double s_ch    = Input(S_CH) ;
  double s_sh    = Input(S_SH) ;
  double c_na    = Input(C_Na) ;
  double c_k     = Input(C_K) ;
  double c_oh    = Input(C_OH) ;
  
  /* Logs */
  double logs_ch    = log10(s_ch) ;
  double logs_sh    = log10(s_sh) ;
  double logc_na    = log10(c_na) ;
  double logc_k     = log10(c_k) ;
  
  /* Trial value before solving electroneutrality */
  double logc_oh    = (c_oh > 0) ? log10(c_oh) : -7 ;
  
  
  /* Hypothesis of unit water activity */
  double loga_h2o = 0 ;
  
  /* Autoprotolysis of water */
  double logk_h2o = LogKeq(H2O) ;
  double logc_h   = logk_h2o + loga_h2o - logc_oh ;
  
  
  /* Chemical reactions involving compounds of type I. */
  
  /* Calcium compounds */
  double logk_ch    = LogKsp(CH) ;
  double logq_ch    = logs_ch + logk_ch ;
  double logc_ca    = logq_ch - 2*(logc_oh) ;
  double logk_caoh  = LogKeq(CaOH) ;
  double logc_caoh  = logc_ca + logc_oh - logk_caoh ; ;
  double logk_caoh2 = LogKeq(CaO2H2) ;
  double logc_caoh2 = logq_ch - logk_caoh2 ;
  
  /* Silicon compounds */
  double logk_sh     = LogKsp(SH) ;
  double logq_sh     = logs_sh + logk_sh ;
  double logc_h4sio4 = logq_sh ;
  double logk_h3sio4 = LogKeq(H3SiO4) ;
  double logc_h3sio4 = logc_h4sio4 + logc_oh - (logk_h3sio4 + loga_h2o) ;
  double logk_h2sio4 = LogKeq(H2SiO4) ;
  double logc_h2sio4 = logc_h3sio4 + logc_oh - (logk_h2sio4 + loga_h2o) ;
  
  /* Sodium compounds */
  double logk_naoh = LogKeq(NaOH) ;
  double logc_naoh = logc_na + logc_oh - logk_naoh ;
  
  /* Potassium compounds */
  double logk_koh = LogKeq(KOH) ;
  double logc_koh = logc_k + logc_oh - logk_koh ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-Silicon compounds */
  double logk_cah2sio4 = LogKeq(CaH2SiO4) ;
  double logc_cah2sio4 = logc_ca + logc_h2sio4 - logk_cah2sio4 ;
  double logk_cah3sio4 = LogKeq(CaH3SiO4) ;
  double logc_cah3sio4 = logc_ca + logc_h3sio4 - logk_cah3sio4 ;
  
  
  /* Set variables to zero for non used variables */
  CementSolutionChemistry_Initialize(csc) ;
  
  
  /* Backup */
  //IAP(CH) = pow(10,logq_ch) ;
  //IAP(SH) = pow(10,logq_sh) ;
  
  //Concentration(H2O)       = pow(10,logc_h2o) ;
  Concentration(OH)        = pow(10,logc_oh) ;
  Concentration(H)         = pow(10,logc_h) ;
  
  Concentration(Ca)        = pow(10,logc_ca) ;
  Concentration(CaOH)      = pow(10,logc_caoh) ;
  Concentration(CaO2H2)    = pow(10,logc_caoh2) ;
  
  Concentration(H4SiO4)    = pow(10,logc_h4sio4) ;
  Concentration(H3SiO4)    = pow(10,logc_h3sio4) ;
  Concentration(H2SiO4)    = pow(10,logc_h2sio4) ;
  
  Concentration(Na)        = pow(10,logc_na) ;
  Concentration(NaOH)      = pow(10,logc_naoh) ;
  
  Concentration(K)         = pow(10,logc_k) ;
  Concentration(KOH)       = pow(10,logc_koh) ;
  
  Concentration(CaH2SiO4)  = pow(10,logc_cah2sio4) ;
  Concentration(CaH3SiO4)  = pow(10,logc_cah3sio4) ;
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_1(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O-SO3 */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O(csc) ;
  
  /* Supplement with Al2O3 */
  double s_ah3    = Input(S_AH3) ;
  double logs_ah3 = log10(s_ah3) ;
  
  /* Autoprotolysis of water */
  double c_oh       = Concentration(OH) ;
  double c_h        = Concentration(H) ;
  double logc_oh    = log10(c_oh) ;
  double logc_h     = log10(c_h) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Aluminium compounds */
  double logk_ah3    = LogKsp(AH3) ;
  double logq_ah3    = logs_ah3 + logk_ah3 ;
  double logc_al     = 0.5*logq_ah3 - 3*(logc_oh) ;
  double logk_alo4h4 = LogKeq(AlO4H4) ;
  double logc_alo4h4 = logc_al + 4*logc_oh - logk_alo4h4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  
  /* Backup */
  {
    double c_ca  = Concentration(Ca) ;
    double c_so4 = Concentration(SO4) ;
    double logc_ca  = log10(c_ca) ;
    double logc_so4 = log10(c_so4) ;
    double loga_h2o = 0 ;
    double logq_afm = 4*logc_ca + 2*logc_al + logc_so4 + 18*loga_h2o - 12*logc_h ;
    double logq_aft = 6*logc_ca + 2*logc_al + 3*logc_so4 + 38*loga_h2o - 12*logc_h ;
    
    //IAP(AH3) = pow(10,logq_ah3) ;
    //IAP(AFm) = pow(10,logq_afm) ;
    //IAP(AFt) = pow(10,logq_aft) ;
  }
  
  Concentration(Al)        = pow(10,logc_al) ;
  Concentration(AlO4H4)    = pow(10,logc_alo4h4) ;
}



void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_2(CementSolutionChemistry_t* csc)
{
  double logq_ch = Input(LogQ_CH) ;
  double logq_sh = Input(LogQ_SH) ;
  double c_na    = Input(C_Na) ;
  double c_k     = Input(C_K) ;
  double c_oh    = Input(C_OH) ;
  
  /* Logs */
  double logc_na    = log10(c_na) ;
  double logc_k     = log10(c_k) ;
  
  /* Trial value before solving electroneutrality */
  double logc_oh    = (c_oh > 0) ? log10(c_oh) : -7 ;
  
  
  /* Hypothesis of unit water activity */
  double loga_h2o = 0 ;
  
  /* Autoprotolysis of water */
  double logk_h2o = LogKeq(H2O) ;
  double logc_h   = logk_h2o + loga_h2o - logc_oh ;
  
  
  /* Chemical reactions involving compounds of type I. */
  
  /* Calcium compounds */
  double logc_ca    = logq_ch - 2*(logc_oh) ;
  double logk_caoh  = LogKeq(CaOH) ;
  double logc_caoh  = logc_ca + logc_oh - logk_caoh ; ;
  double logk_caoh2 = LogKeq(CaO2H2) ;
  double logc_caoh2 = logq_ch - logk_caoh2 ;
  
  /* Silicon compounds */
  double logc_h4sio4 = logq_sh ;
  double logk_h3sio4 = LogKeq(H3SiO4) ;
  double logc_h3sio4 = logc_h4sio4 + logc_oh - (logk_h3sio4 + loga_h2o) ;
  double logk_h2sio4 = LogKeq(H2SiO4) ;
  double logc_h2sio4 = logc_h3sio4 + logc_oh - (logk_h2sio4 + loga_h2o) ;
  
  /* Sodium compounds */
  double logk_naoh = LogKeq(NaOH) ;
  double logc_naoh = logc_na + logc_oh - logk_naoh ;
  
  /* Potassium compounds */
  double logk_koh = LogKeq(KOH) ;
  double logc_koh = logc_k + logc_oh - logk_koh ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-Silicon compounds */
  double logk_cah2sio4 = LogKeq(CaH2SiO4) ;
  double logc_cah2sio4 = logc_ca + logc_h2sio4 - logk_cah2sio4 ;
  double logk_cah3sio4 = LogKeq(CaH3SiO4) ;
  double logc_cah3sio4 = logc_ca + logc_h3sio4 - logk_cah3sio4 ;
  
  
  /* Set variables to zero for non used variables */
  CementSolutionChemistry_Initialize(csc) ;
  
  
  /* Backup */
  //IAP(CH) = pow(10,logq_ch) ;
  //IAP(SH) = pow(10,logq_sh) ;
  
  //Concentration(H2O)       = pow(10,logc_h2o) ;
  Concentration(OH)        = pow(10,logc_oh) ;
  Concentration(H)         = pow(10,logc_h) ;
  
  Concentration(Ca)        = pow(10,logc_ca) ;
  Concentration(CaOH)      = pow(10,logc_caoh) ;
  Concentration(CaO2H2)    = pow(10,logc_caoh2) ;
  
  Concentration(H4SiO4)    = pow(10,logc_h4sio4) ;
  Concentration(H3SiO4)    = pow(10,logc_h3sio4) ;
  Concentration(H2SiO4)    = pow(10,logc_h2sio4) ;
  
  Concentration(Na)        = pow(10,logc_na) ;
  Concentration(NaOH)      = pow(10,logc_naoh) ;
  
  Concentration(K)         = pow(10,logc_k) ;
  Concentration(KOH)       = pow(10,logc_koh) ;
  
  Concentration(CaH2SiO4)  = pow(10,logc_cah2sio4) ;
  Concentration(CaH3SiO4)  = pow(10,logc_cah3sio4) ;
}



void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2_2(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_2(csc) ;
  
  /* Supplement with CO2 */
  double c_co2     = Input(C_CO2) ;
  double logc_co2  = log10(c_co2) ;
  
  
  /* Hypothesis of unit water activity */
  double loga_h2o = 0 ;
  
  /* Autoprotolysis of water */
  double c_oh       = Concentration(OH) ;
  double logc_oh    = log10(c_oh) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Carbon compounds */
  double logk_h2co3 = LogKeq(H2CO3) ;
  double logc_h2co3 = logc_co2 + loga_h2o - logk_h2co3 ;
  double logk_hco3  = LogKeq(HCO3) ;
  double logc_hco3  = logc_h2co3 + logc_oh - (logk_hco3 + loga_h2o) ;
  double logk_co3   = LogKeq(CO3) ;
  double logc_co3   = logc_hco3  + logc_oh - (logk_co3  + loga_h2o) ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-Carbon compounds */
  double c_ca        = Concentration(Ca) ;
  double logc_ca     = log10(c_ca) ;
  double logk_cahco3 = LogKeq(CaHCO3) ;
  double logc_cahco3 = logc_ca + logc_hco3 - logk_cahco3 ;
  double logk_caco3  = LogKeq(CaCO3) ;
  double logc_caco3  = logc_ca + logc_co3  - logk_caco3 ;
  
  /* Sodium-Carbon compounds */
  double c_na        = Concentration(Na) ;
  double logc_na     = log10(c_na) ;
  double logk_nahco3 = LogKeq(NaHCO3) ;
  double logc_nahco3 = logc_na + logc_hco3 - logk_nahco3 ;
  double logk_naco3  = LogKeq(NaCO3) ;
  double logc_naco3  = logc_na + logc_co3  - logk_naco3 ;
  
  
  /* Backup */
  {
    double logq_cc     = logc_ca + logc_co3 ;
    
    //IAP(CC) = pow(10,logq_cc) ;
  }
  
  Concentration(H2CO3)     = pow(10,logc_h2co3) ;
  Concentration(HCO3)      = pow(10,logc_hco3) ;
  Concentration(CO3)       = pow(10,logc_co3) ;
  Concentration(CO2)       = pow(10,logc_co2) ;
  
  Concentration(CaHCO3)    = pow(10,logc_cahco3) ;
  Concentration(CaCO3)     = pow(10,logc_caco3) ;
  
  Concentration(NaHCO3)    = pow(10,logc_nahco3) ;
  Concentration(NaCO3)     = pow(10,logc_naco3) ;
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_2(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_2(csc) ;
  
  /* Supplement with SO3 */
  double c_h2so4    = Input(C_H2SO4) ;
  double logc_h2so4 = log10(c_h2so4) ;
  
  
  /* Autoprotolysis of water */
  double c_h      = Concentration(H) ;
  double logc_h   = log10(c_h) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Sulfur compounds */
  double logk_h2so4 = LogKeq(H2SO4) ;
  double logc_hso4  = logc_h2so4 - logc_h + logk_h2so4 ;
  double logk_hso4  = LogKeq(HSO4) ;
  double logc_so4   = logc_hso4 - logc_h + logk_hso4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  /* Calcium-sulfur compounds */
  double c_ca        = Concentration(Ca) ;
  double logc_ca     = log10(c_ca) ;
  double logk_cahso4 = LogKeq(CaHSO4) ;
  double logc_cahso4 = logc_ca + logc_hso4 - logk_cahso4 ;
  double logk_caso4  = LogKeq(CaSO4) ;
  double logc_caso4  = logc_ca + logc_so4  - logk_caso4 ;
  
  /* Sodium-sulfur compounds */
  
  
  
  /* Backup */
  {
    double loga_h2o = 0 ;
    double logq_csh2 = logc_ca + logc_so4 + 2*loga_h2o ;
    
    //IAP(CSH2) = pow(10,logq_csh2) ;
  }
  
  Concentration(H2SO4)     = pow(10,logc_h2so4) ;
  Concentration(HSO4)      = pow(10,logc_hso4) ;
  Concentration(SO4)       = pow(10,logc_so4) ;
  
  Concentration(CaHSO4)    = pow(10,logc_cahso4) ;
  Concentration(CaSO4)     = pow(10,logc_caso4) ;
}





void CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_2(CementSolutionChemistry_t* csc)
{
  /* Compute the system CaO-SiO2-Na2O-K2O-SO3 */
  CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_2(csc) ;
  
  /* Supplement with Al2O3 */
  double logq_ah3   = Input(LogQ_AH3) ;
  
  /* Autoprotolysis of water */
  double c_oh       = Concentration(OH) ;
  double c_h        = Concentration(H) ;
  double logc_oh    = log10(c_oh) ;
  double logc_h     = log10(c_h) ;
  
  
  /* Chemical reactions involving compounds of type I. */

  /* Aluminium compounds */
  double logc_al     = 0.5*logq_ah3 - 3*logc_oh ;
  double logk_alo4h4 = LogKeq(AlO4H4) ;
  double logc_alo4h4 = logc_al + 4*logc_oh - logk_alo4h4 ;
  
  
  /* Chemical reactions involving compounds of type II. */
  
  
  /* Backup */
  {
    double c_ca  = Concentration(Ca) ;
    double c_so4 = Concentration(SO4) ;
    double logc_ca  = log10(c_ca) ;
    double logc_so4 = log10(c_so4) ;
    double loga_h2o = 0 ;
    double logq_afm = 4*logc_ca + 2*logc_al + logc_so4 + 18*loga_h2o - 12*logc_h ;
    double logq_aft = 6*logc_ca + 2*logc_al + 3*logc_so4 + 38*loga_h2o - 12*logc_h ;
    double logq_c3ah6 = 3*logc_ca + 2*logc_al + 12*loga_h2o - 12*logc_h ;
    
    //IAP(AH3) = pow(10,logq_ah3) ;
    //IAP(AFm) = pow(10,logq_afm) ;
    //IAP(AFt) = pow(10,logq_aft) ;
    //IAP(C3AH6) = pow(10,logq_c3ah6) ;
  }
  
  Concentration(Al)        = pow(10,logc_al) ;
  Concentration(AlO4H4)    = pow(10,logc_alo4h4) ;
}
#endif
