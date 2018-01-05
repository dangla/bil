#ifndef CEMENTSOLUTIONCHEMISTRY_H
#define CEMENTSOLUTIONCHEMISTRY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct CementSolutionChemistry_s     ; 
typedef struct CementSolutionChemistry_s     CementSolutionChemistry_t ;



extern CementSolutionChemistry_t* (CementSolutionChemistry_Create)(const int) ;

extern void   (CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_H2O)          (CementSolutionChemistry_t*) ;
extern void   (CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2_H2O)      (CementSolutionChemistry_t*) ;
extern void   (CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O)      (CementSolutionChemistry_t*) ;
extern void   (CementSolutionChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O)(CementSolutionChemistry_t*) ;
extern void   (CementSolutionChemistry_PrintChemicalConstants)(CementSolutionChemistry_t*) ;

extern double (CementSolutionChemistry_SolveElectroneutrality)         (CementSolutionChemistry_t*) ;
extern double (CementSolutionChemistry_SolveExplicitElectroneutrality) (CementSolutionChemistry_t*) ;

extern void   (CementSolutionChemistry_UpdateSolution)(CementSolutionChemistry_t*) ;

extern void   (CementSolutionChemistry_CopyConcentrations)   (CementSolutionChemistry_t*,double*) ;
extern void   (CementSolutionChemistry_CopyLogConcentrations)(CementSolutionChemistry_t*,double*) ;
extern void   (CementSolutionChemistry_CopyChemicalPotential)(CementSolutionChemistry_t*,double*) ;

extern double*   (CementSolutionChemistry_GetValence)(void) ;



#define CementSolutionChemistry_GetTemperature(CSC) \
        ((CSC)->temperature)
        
#define CementSolutionChemistry_GetPrimaryVariable(CSC) \
        ((CSC)->primaryvariable)

#define CementSolutionChemistry_GetConcentration(CSC) \
        ((CSC)->concentration)

#define CementSolutionChemistry_GetLogConcentration(CSC) \
       ((CSC)->logconcentration)

#define CementSolutionChemistry_GetActivity(CSC) \
        ((CSC)->activity)

#define CementSolutionChemistry_GetLogActivity(CSC) \
       ((CSC)->logactivity)

#define CementSolutionChemistry_GetElementConcentration(CSC) \
        ((CSC)->elementconcentration)

#define CementSolutionChemistry_GetOtherVariable(CSC) \
        ((CSC)->othervariable)

#define CementSolutionChemistry_GetLog10Keq(CSC) \
        ((CSC)->log10equilibriumconstant)

#define CementSolutionChemistry_GetElectricPotential(CSC) \
        ((CSC)->electricpotential)




/* Macros for the room temperature
 * -------------------------------*/
#define CementSolutionChemistry_GetRoomTemperature(CSC) \
        Temperature_GetRoomValue(CementSolutionChemistry_GetTemperature(CSC))
        

#define CementSolutionChemistry_SetRoomTemperature(CSC,T) \
        Temperature_SetRoomTemperature(CementSolutionChemistry_GetTemperature(CSC),T)



/* Macros for primary variables
 * ----------------------------*/
 
#define CementSolutionChemistry_NbOfPrimaryVariables  (8)

/* Different primary variables may be used */
#define CementSolutionChemistry_LogQ_CH      (0)
#define CementSolutionChemistry_LogA_Ca      (0) // not used yet
#define CementSolutionChemistry_LogQ_SH      (1)
#define CementSolutionChemistry_LogA_H4SiO4  (1) // not used yet
#define CementSolutionChemistry_LogQ_AH3     (2)
#define CementSolutionChemistry_LogA_AlO4H4  (2) // not used yet

#define CementSolutionChemistry_LogA_Na      (3)
#define CementSolutionChemistry_LogA_K       (4)
#define CementSolutionChemistry_LogA_CO2     (5)
#define CementSolutionChemistry_LogA_OH      (6)
#define CementSolutionChemistry_LogA_H2SO4   (7)

#define CementSolutionChemistry_LogC_Na      (3)
#define CementSolutionChemistry_LogC_K       (4)
#define CementSolutionChemistry_LogC_CO2     (5)
#define CementSolutionChemistry_LogC_OH      (6)
#define CementSolutionChemistry_LogC_H2SO4   (7)

/* Input */
#define CementSolutionChemistry_GetInput(CSC,U) \
       (CementSolutionChemistry_GetPrimaryVariable(CSC)[CementSolutionChemistry_##U])

       
       
/* List of compound names
 * ----------------------*/
#define CementSolutionChemistry_NbOfSpecies  (31)
 
#define CementSolutionChemistry_C_0   H2O
#define CementSolutionChemistry_C_1   H
#define CementSolutionChemistry_C_2   OH

#define CementSolutionChemistry_C_3   Ca
#define CementSolutionChemistry_C_4   CaOH
#define CementSolutionChemistry_C_5   CaO2H2

#define CementSolutionChemistry_C_6   H2SiO4
#define CementSolutionChemistry_C_7   H3SiO4
#define CementSolutionChemistry_C_8   H4SiO4

#define CementSolutionChemistry_C_9   CaH2SiO4
#define CementSolutionChemistry_C_10  CaH3SiO4

#define CementSolutionChemistry_C_11  Na
#define CementSolutionChemistry_C_12  NaOH

#define CementSolutionChemistry_C_13  K
#define CementSolutionChemistry_C_14  KOH

#define CementSolutionChemistry_C_15  H2CO3
#define CementSolutionChemistry_C_16  HCO3
#define CementSolutionChemistry_C_17  CO3
#define CementSolutionChemistry_C_18  CO2

#define CementSolutionChemistry_C_19  CaHCO3
#define CementSolutionChemistry_C_20  CaCO3

#define CementSolutionChemistry_C_21  NaHCO3
#define CementSolutionChemistry_C_22  NaCO3

#define CementSolutionChemistry_C_23  H2SO4
#define CementSolutionChemistry_C_24  HSO4
#define CementSolutionChemistry_C_25  SO4
/*
#define CementSolutionChemistry_C_XX  SO3
#define CementSolutionChemistry_C_XX  S2O3

#define CementSolutionChemistry_C_XX  H2S
#define CementSolutionChemistry_C_XX  HS
#define CementSolutionChemistry_C_XX  S
#define CementSolutionChemistry_C_XX  S0
*/

#define CementSolutionChemistry_C_26  CaHSO4
#define CementSolutionChemistry_C_27  CaSO4

#define CementSolutionChemistry_C_28  Cl

#define CementSolutionChemistry_C_29  Al
#define CementSolutionChemistry_C_30  AlO4H4


        

#define CementSolutionChemistry_ListOfCompounds (\
        H2O,H,OH \
       ,Ca,CaOH,CaO2H2 \
       ,H2SiO4,H3SiO4,H4SiO4 \
       ,CaH2SiO4,CaH3SiO4 \
       ,Na,NaOH \
       ,K,KOH \
       ,H2CO3,HCO3,CO3,CO2 \
       ,CaHCO3,CaCO3 \
       ,NaHCO3,NaCO3 \
       ,H2SO4,HSO4,SO4 \
       ,CaHSO4,CaSO4 \
       ,Cl \
       ,Al,AlO4H4 \
       )
       
//       ,H2SO4,HSO4,SO4,SO3,S2O3 \
       ,H2S,HS,S,S0 \



#include "Utils.h"

#define CementSolutionChemistry_GetIndexOf(CPD) \
        Utils_CAT(CementSolutionChemistry_A_,CPD)


#include "Algos.h"

//#define CementSolutionChemistry_ENUM \
          Tuple_SEQ(Algos_MAP(CementSolutionChemistry_ListOfCompounds,CementSolutionChemistry_GetIndexOf))


enum CementSolutionChemistry_e {
  CementSolutionChemistry_ENUM
} ;



/* Macros for the concentrations
 * -----------------------------*/

#define CementSolutionChemistry_A_H2O         (0)
#define CementSolutionChemistry_A_H           (1)
#define CementSolutionChemistry_A_OH          (2)

#define CementSolutionChemistry_A_Ca          (3)
#define CementSolutionChemistry_A_CaOH        (4)
#define CementSolutionChemistry_A_CaO2H2      (5)

#define CementSolutionChemistry_A_H2SiO4      (6)
#define CementSolutionChemistry_A_H3SiO4      (7)
#define CementSolutionChemistry_A_H4SiO4      (8)

#define CementSolutionChemistry_A_CaH2SiO4    (9)
#define CementSolutionChemistry_A_CaH3SiO4    (10)

#define CementSolutionChemistry_A_Na          (11)
#define CementSolutionChemistry_A_NaOH        (12)

#define CementSolutionChemistry_A_K           (13)
#define CementSolutionChemistry_A_KOH         (14)

#define CementSolutionChemistry_A_H2CO3       (15)
#define CementSolutionChemistry_A_HCO3        (16)
#define CementSolutionChemistry_A_CO3         (17)
#define CementSolutionChemistry_A_CO2         (18)

#define CementSolutionChemistry_A_CaHCO3      (19)
#define CementSolutionChemistry_A_CaCO3       (20)

#define CementSolutionChemistry_A_NaHCO3      (21)
#define CementSolutionChemistry_A_NaCO3       (22)

#define CementSolutionChemistry_A_H2SO4       (23)
#define CementSolutionChemistry_A_HSO4        (24)
#define CementSolutionChemistry_A_SO4         (25)
#define CementSolutionChemistry_A_SO3         (xx)
#define CementSolutionChemistry_A_S2O3        (xx)

#define CementSolutionChemistry_A_H2S         (xx)
#define CementSolutionChemistry_A_HS          (xx)
#define CementSolutionChemistry_A_S           (xx)
#define CementSolutionChemistry_A_S0          (xx)

#define CementSolutionChemistry_A_CaHSO4      (26)
#define CementSolutionChemistry_A_CaSO4       (27)

#define CementSolutionChemistry_A_Cl          (28)

#define CementSolutionChemistry_A_Al          (29)
#define CementSolutionChemistry_A_AlO4H4      (30)


#define CementSolutionChemistry_GetConcentrationOf(CSC,CPD) \
       (CementSolutionChemistry_GetConcentration(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])

       
#define CementSolutionChemistry_GetLogConcentrationOf(CSC,CPD) \
       (CementSolutionChemistry_GetLogConcentration(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])


#define CementSolutionChemistry_GetActivityOf(CSC,CPD) \
       (CementSolutionChemistry_GetActivity(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])

       
#define CementSolutionChemistry_GetLogActivityOf(CSC,CPD) \
       (CementSolutionChemistry_GetLogActivity(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])



/* Macros for other variables
 * --------------------------*/
#define CementSolutionChemistry_NbOfOtherVariables       (3)

/* 1. Liquid mass density */
#define CementSolutionChemistry_GetLiquidMassDensity(CSC) \
       (CementSolutionChemistry_GetOtherVariable(CSC)[0])


/* 2. Charge density */
#define CementSolutionChemistry_GetChargeDensity(CSC) \
       (CementSolutionChemistry_GetOtherVariable(CSC)[1])


/* 3. Ionic strength */
#define CementSolutionChemistry_GetIonicStrength(CSC) \
       (CementSolutionChemistry_GetOtherVariable(CSC)[2])



/* Macros for element concentrations
 * ---------------------------------*/
#define CementSolutionChemistry_NbOfElementConcentrations       (8)

#define CementSolutionChemistry_E_Ca          (0)
#define CementSolutionChemistry_E_Si          (1)
#define CementSolutionChemistry_E_Na          (2)
#define CementSolutionChemistry_E_K           (3)
#define CementSolutionChemistry_E_C           (4)
#define CementSolutionChemistry_E_S           (5)
#define CementSolutionChemistry_E_Al          (6)
#define CementSolutionChemistry_E_Cl          (7)

#define CementSolutionChemistry_GetElementConcentrationOf(CSC,A) \
       (CementSolutionChemistry_GetElementConcentration(CSC)[CementSolutionChemistry_E_##A])






#if 0
/* Macros for ion charges
 * ----------------------*/
/* Ion charges */
#define CementSolutionChemistry_IonCharge(CSC,CPD) \
       (CementSolutionChemistry_GetChargeVariable(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])
       
/* Charge density */
#define CementSolutionChemistry_ChargeDensity(CSC) \
       (CementSolutionChemistry_GetChargeVariable(CSC)[CementSolutionChemistry_CHARGE])
       
/* Ionic strength */
#define CementSolutionChemistry_IonicStrength(CSC) \
       (CementSolutionChemistry_GetChargeVariable(CSC)[CementSolutionChemistry_IoSth])
       
/* Element charges */
#define CementSolutionChemistry_ElementCharge(CSC,A) \
       (CementSolutionChemistry_GetChargeVariable(CSC)[CementSolutionChemistry_E_##A])


#define CementSolutionChemistry_NbOfChargeVariables   (40)
#endif





/* Macros for equilibrium constants (same indices as A_CPD)
 * -------------------------------------------------------*/
#define CementSolutionChemistry_GetLog10EquilibriumConstant(CSC,CPD) \
       (CementSolutionChemistry_GetLog10Keq(CSC)[CementSolutionChemistry_GetIndexOf(CPD)])




/* Macro for the resolution of the system
 * --------------------------------------*/
#define CementSolutionChemistry_ComputeSystem(CSC,SYS) \
       (CementSolutionChemistry_ComputeSystem_##SYS(CSC))


#include "Temperature.h"

struct CementSolutionChemistry_s {
  Temperature_t* temperature ;
//  int nbofprimaryvariables ;
//  int nbofvariables ;
  double* primaryvariable ;
  double* concentration ;
  double* logconcentration ;
  double* activity ;
  double* logactivity ;
  double* elementconcentration ;
  double* othervariable ;
  double* log10equilibriumconstant ;
  double  electricpotential ;
} ;

#endif
