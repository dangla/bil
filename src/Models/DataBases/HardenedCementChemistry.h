#ifndef HARDENEDCEMENTCHEMISTRY_H
#define HARDENEDCEMENTCHEMISTRY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct HardenedCementChemistry_s     ; 
typedef struct HardenedCementChemistry_s     HardenedCementChemistry_t ;


extern HardenedCementChemistry_t* (HardenedCementChemistry_Create)(void) ;
extern HardenedCementChemistry_t* (HardenedCementChemistry_GetInstance)(void) ;
//extern void (HardenedCementChemistry_SetTemperature)(double) ;
extern void (HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2_H2O)        (HardenedCementChemistry_t*) ;
extern void (HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O)        (HardenedCementChemistry_t*) ;
extern void (HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_2)      (HardenedCementChemistry_t*) ;
extern void (HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O)  (HardenedCementChemistry_t*) ;
extern void (HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O_2)(HardenedCementChemistry_t*) ;
extern void (HardenedCementChemistry_PrintChemicalConstants)(HardenedCementChemistry_t*) ;


/* Synonyms */
#define HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2 \
        HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_CO2_H2O
        
#define HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3 \
        HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O
        
#define HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_2 \
        HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_H2O_2
        
#define HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3 \
        HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O
        
#define HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_2 \
        HardenedCementChemistry_ComputeSystem_CaO_SiO2_Na2O_K2O_SO3_Al2O3_H2O_2
        


#define HardenedCementChemistry_CopyConcentrations(hcc,v)  \
        (CementSolutionChemistry_CopyConcentrations(HardenedCementChemistry_GetCementSolutionChemistry(hcc),v))

#define HardenedCementChemistry_CopyLogConcentrations(hcc,v)  \
        (CementSolutionChemistry_CopyLogConcentrations(HardenedCementChemistry_GetCementSolutionChemistry(hcc),v))

#define HardenedCementChemistry_CopyChemicalPotential(hcc,v)  \
        (CementSolutionChemistry_CopyChemicalPotential(HardenedCementChemistry_GetCementSolutionChemistry(hcc),v))


#include "Temperature.h"

#define HardenedCementChemistry_SetTemperature(T) \
       (Temperature_RoomValue = (T))


/* The getters for attributes */
#define HardenedCementChemistry_GetTemperature(hcc) \
        ((hcc)->temperature)
#define HardenedCementChemistry_GetPrimaryVariable(hcc) \
        ((hcc)->primaryvariable)
#define HardenedCementChemistry_GetVariable(hcc) \
        ((hcc)->variable)
#define HardenedCementChemistry_GetSaturationIndex(hcc) \
        ((hcc)->saturationindex)
#define HardenedCementChemistry_GetConstant(hcc) \
        ((hcc)->constant)
#define HardenedCementChemistry_GetLog10Ksp(hcc) \
        ((hcc)->log10solubilityproductconstant)
#define HardenedCementChemistry_GetCSHCurves(hcc) \
        ((hcc)->cshcurves)
#define HardenedCementChemistry_GetCementSolutionChemistry(hcc) \
        ((hcc)->csc)
#define HardenedCementChemistry_GetLog10SaturationIndex(hcc) \
        ((hcc)->log10saturationindex)
#define HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(hcc) \
        ((hcc)->curveofcalciumsiliconratioincsh)
#define HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(hcc) \
        ((hcc)->curveofwatersiliconratioincsh)
#define HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(hcc) \
        ((hcc)->curveofsaturationindexofsh)



/* Macro for the electric potential
 * --------------------------------*/
#define HardenedCementChemistry_GetElectricPotential(hcc) \
        (CementSolutionChemistry_GetElectricPotential(HardenedCementChemistry_GetCementSolutionChemistry(hcc)))




/* Macros for primary variables
 * ----------------------------*/
#define HardenedCementChemistry_NbOfPrimaryVariables  (8)

#define HardenedCementChemistry_SI_Ca        (0)
#define HardenedCementChemistry_SI_Si        (1)
#define HardenedCementChemistry_SI_Al        (2)

#define HardenedCementChemistry_LogA_Na      (3)
#define HardenedCementChemistry_LogA_K       (4)
#define HardenedCementChemistry_LogA_CO2     (5)
#define HardenedCementChemistry_LogA_OH      (6)
#define HardenedCementChemistry_LogA_H2SO4   (7)

#define HardenedCementChemistry_LogC_Na      (3)
#define HardenedCementChemistry_LogC_K       (4)
#define HardenedCementChemistry_LogC_CO2     (5)
#define HardenedCementChemistry_LogC_OH      (6)
#define HardenedCementChemistry_LogC_H2SO4   (7)

/* Input */
#define HardenedCementChemistry_GetInput(hcc,U) \
        (HardenedCementChemistry_GetPrimaryVariable(hcc)[HardenedCementChemistry_##U])



/* Macros for saturation indexes
 * -----------------------------*/
#define HardenedCementChemistry_NbOfSaturationIndexes  (10)

#define HardenedCementChemistry_S_CH          (0)
#define HardenedCementChemistry_S_SH          (1)
#define HardenedCementChemistry_S_CC          (2)
#define HardenedCementChemistry_S_CSH2        (3)
#define HardenedCementChemistry_S_AH3         (4)
#define HardenedCementChemistry_S_AFm         (5)
#define HardenedCementChemistry_S_AFt         (6)
#define HardenedCementChemistry_S_C3AH6       (7)
#define HardenedCementChemistry_S_C2AH8       (8)
#define HardenedCementChemistry_S_CAH10       (9)

#define HardenedCementChemistry_GetSaturationIndexOf(hcc,A) \
        (HardenedCementChemistry_GetSaturationIndex(hcc)[HardenedCementChemistry_S_##A])

#define HardenedCementChemistry_GetLog10SaturationIndexOf(hcc,A) \
        (HardenedCementChemistry_GetLog10SaturationIndex(hcc)[HardenedCementChemistry_S_##A])

/* Saturation index of CH */
#define HardenedCementChemistry_GetSaturationIndexOfCH(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,CH)
        
/* Saturation index of CC */
#define HardenedCementChemistry_GetSaturationIndexOfCC(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,CC)
        
/* Saturation index of SH */
#define HardenedCementChemistry_GetSaturationIndexOfSH(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,SH)
        
/* Saturation index of CSH2 */
#define HardenedCementChemistry_GetSaturationIndexOfCSH2(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,CSH2)

/* Saturation index of AH3 */
#define HardenedCementChemistry_GetSaturationIndexOfAH3(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,AH3)

/* Saturation index of AFm */
#define HardenedCementChemistry_GetSaturationIndexOfAFm(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,AFm)

/* Saturation index of AFt */
#define HardenedCementChemistry_GetSaturationIndexOfAFt(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,AFt)

/* Saturation index of C3AH6 */
#define HardenedCementChemistry_GetSaturationIndexOfC3AH6(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,C3AH6)

/* Saturation index of C2AH8 */
#define HardenedCementChemistry_GetSaturationIndexOfC2AH8(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,C2AH8)

/* Saturation index of CAH10 */
#define HardenedCementChemistry_GetSaturationIndexOfCAH10(hcc) \
        HardenedCementChemistry_GetSaturationIndexOf(hcc,CAH10)





/* Macros for variables
 * --------------------*/
#define HardenedCementChemistry_NbOfVariables  (3)

#define HardenedCementChemistry_X_CSH         (0)
#define HardenedCementChemistry_Z_CSH         (1)
//#define HardenedCementChemistry_V_CSH         (2)

/* Ca/Si ratio in C-S-H */
#define HardenedCementChemistry_GetCalciumSiliconRatioInCSH(hcc) \
        (HardenedCementChemistry_GetVariable(hcc)[HardenedCementChemistry_X_CSH])

/* Water/Si ratio in C-S-H */
#define HardenedCementChemistry_GetWaterSiliconRatioInCSH(hcc) \
        (HardenedCementChemistry_GetVariable(hcc)[HardenedCementChemistry_Z_CSH])

/* Molar volume of C-S-H */
/*
#define HardenedCementChemistry_GetMolarVolumeOfCSH(hcc) \
       (HardenedCementChemistry_GetVariable(hcc)[HardenedCementChemistry_V_CSH])
*/



/* Macros for Aqueous concentrations
 * ---------------------------------*/
#define HardenedCementChemistry_NbOfSpecies \
        CementSolutionChemistry_NbOfSpecies
        
#define HardenedCementChemistry_GetAqueousConcentration(hcc) \
        CementSolutionChemistry_GetConcentration(HardenedCementChemistry_GetCementSolutionChemistry(hcc))

#define HardenedCementChemistry_GetLogAqueousConcentration(hcc,CPD) \
        CementSolutionChemistry_GetLogConcentration(HardenedCementChemistry_GetCementSolutionChemistry(hcc))
        
#define HardenedCementChemistry_GetAqueousConcentrationOf(hcc,CPD) \
        CementSolutionChemistry_GetConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(hcc),CPD)

#define HardenedCementChemistry_GetLogAqueousConcentrationOf(hcc,CPD) \
        CementSolutionChemistry_GetLogConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(hcc),CPD)
        

/* Synonyms */
#define HardenedCementChemistry_NbOfConcentrations \
        HardenedCementChemistry_NbOfSpecies


/* Other macros
 * ------------*/
/* Liquid mass density */
#define HardenedCementChemistry_GetLiquidMassDensity(hcc) \
        CementSolutionChemistry_GetLiquidMassDensity(HardenedCementChemistry_GetCementSolutionChemistry(hcc))

/* Liquid charge density */
#define HardenedCementChemistry_GetLiquidChargeDensity(hcc) \
        CementSolutionChemistry_GetChargeDensity(HardenedCementChemistry_GetCementSolutionChemistry(hcc))

/* Ionic strength of aqueous phase*/
#define HardenedCementChemistry_GetIonicStrength(hcc) \
        CementSolutionChemistry_GetIonicStrength(HardenedCementChemistry_GetCementSolutionChemistry(hcc))

/* Element aqueous concentrations */
#define HardenedCementChemistry_GetElementAqueousConcentration(hcc) \
        CementSolutionChemistry_GetElementConcentration(HardenedCementChemistry_GetCementSolutionChemistry(hcc))
        
#define HardenedCementChemistry_GetElementAqueousConcentrationOf(hcc,A) \
        CementSolutionChemistry_GetElementConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(hcc),A)






/* Macros for constants
 * --------------------*/
#define HardenedCementChemistry_NbOfConstants (2)

#define HardenedCementChemistry_A_CO2_EQ      (0)
#define HardenedCementChemistry_A_H2SO4_EQ    (1)


/* Equilibrium CO2 concentration */
#define HardenedCementChemistry_GetLog10EquilibriumCO2Activity(hcc) \
        (HardenedCementChemistry_GetConstant(hcc)[HardenedCementChemistry_A_CO2_EQ])

/* Equilibrium H2SO4 concentration */
#define HardenedCementChemistry_GetLog10EquilibriumH2SO4Activity(hcc) \
        (HardenedCementChemistry_GetConstant(hcc)[HardenedCementChemistry_A_H2SO4_EQ])




/* Macros for solubility product constants
 * ---------------------------------------*/
#define HardenedCementChemistry_NbOfSolubilityProductConstants (10)

#define HardenedCementChemistry_K_CH        (0)
#define HardenedCementChemistry_K_SH        (1)
#define HardenedCementChemistry_K_CC        (2)
#define HardenedCementChemistry_K_CSH2      (3)
#define HardenedCementChemistry_K_AH3       (4)
#define HardenedCementChemistry_K_AFm       (5)
#define HardenedCementChemistry_K_AFt       (6)
#define HardenedCementChemistry_K_C3AH6     (7)
#define HardenedCementChemistry_K_C2AH8     (8)
#define HardenedCementChemistry_K_CAH10     (9)

#define HardenedCementChemistry_GetLog10SolubilityProductConstantOf(hcc,CPD) \
        (HardenedCementChemistry_GetLog10Ksp(hcc)[HardenedCementChemistry_K_##CPD])



/* Macros for the resolution of the systems
 * ----------------------------------------*/
#define HardenedCementChemistry_ComputeSystem(hcc,SYS) \
        (HardenedCementChemistry_ComputeSystem_##SYS(hcc))

#define HardenedCementChemistry_SolveElectroneutrality(hcc) \
        (CementSolutionChemistry_SolveElectroneutrality(HardenedCementChemistry_GetCementSolutionChemistry(hcc)))

#define HardenedCementChemistry_SolveExplicitElectroneutrality(hcc) \
        (CementSolutionChemistry_SolveExplicitElectroneutrality(HardenedCementChemistry_GetCementSolutionChemistry(hcc)))



#include "Curves.h"

/* Macros for the CSH curves */
#define HardenedCementChemistry_GetCSHCurve(hcc) \
        (Curves_GetCurve(HardenedCementChemistry_GetCSHCurves(hcc)))

#define HardenedCementChemistry_SetDefaultCurveOfCalciumSiliconRatioInCSH(hcc) \
        do { \
          HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(hcc) = HardenedCementChemistry_GetCSHCurve(hcc) + (0) ; \
        } while(0) 

#define HardenedCementChemistry_SetDefaultCurveOfWaterSiliconRatioInCSH(hcc) \
        do { \
          HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(hcc) = HardenedCementChemistry_GetCSHCurve(hcc) + (1) ; \
        } while(0) 

#define HardenedCementChemistry_SetDefaultCurveOfSaturationIndexOfSH(hcc) \
        do { \
          HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(hcc) = HardenedCementChemistry_GetCSHCurve(hcc) + (2) ; \
        } while(0) 

/*
#define HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(hcc) \
        (HardenedCementChemistry_GetCSHCurve(hcc) + (0))

#define HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(hcc) \
        (HardenedCementChemistry_GetCSHCurve(hcc) + (1))

#define HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(hcc) \
        (HardenedCementChemistry_GetCSHCurve(hcc) + (2))
*/


#include "CementSolutionChemistry.h"

struct HardenedCementChemistry_s {
  double temperature ;
//  int nbofprimaryvariables ;
//  int nbofvariables ;
  double* primaryvariable ;
  double* constant ;
  double* variable ;
  double* saturationindex ;
  double* log10saturationindex ;
  double* log10solubilityproductconstant ;
  CementSolutionChemistry_t* csc ;
  Curves_t* cshcurves ;
  Curve_t*  curveofcalciumsiliconratioincsh ;
  Curve_t*  curveofwatersiliconratioincsh ;
  Curve_t*  curveofsaturationindexofsh ;
} ;

#endif
