#ifndef HARDENEDCEMENTCHEMISTRY_H
#define HARDENEDCEMENTCHEMISTRY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct HardenedCementChemistry_s     ; 
typedef struct HardenedCementChemistry_s     HardenedCementChemistry_t ;


extern HardenedCementChemistry_t* (HardenedCementChemistry_Create)(void) ;
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
        


#define HardenedCementChemistry_CopyConcentrations(HCC,v)  \
        (CementSolutionChemistry_CopyConcentrations(HardenedCementChemistry_GetCementSolutionChemistry(HCC),v))

#define HardenedCementChemistry_CopyLogConcentrations(HCC,v)  \
        (CementSolutionChemistry_CopyLogConcentrations(HardenedCementChemistry_GetCementSolutionChemistry(HCC),v))

#define HardenedCementChemistry_CopyChemicalPotential(HCC,v)  \
        (CementSolutionChemistry_CopyChemicalPotential(HardenedCementChemistry_GetCementSolutionChemistry(HCC),v))


#include "Temperature.h"

/* Macros for the temperature
 * --------------------------*/
#define HardenedCementChemistry_GetRoomTemperature(HCC) \
        CementSolutionChemistry_GetRoomTemperature(HardenedCementChemistry_GetCementSolutionChemistry(HCC))

       
#define HardenedCementChemistry_SetRoomTemperature(HCC,T) \
        CementSolutionChemistry_SetRoomTemperature(HardenedCementChemistry_GetCementSolutionChemistry(HCC),T)


/* The getters for attributes */
#define HardenedCementChemistry_GetPrimaryVariable(HCC) \
        ((HCC)->primaryvariable)
#define HardenedCementChemistry_GetVariable(HCC) \
        ((HCC)->variable)
#define HardenedCementChemistry_GetSaturationIndex(HCC) \
        ((HCC)->saturationindex)
#define HardenedCementChemistry_GetConstant(HCC) \
        ((HCC)->constant)
#define HardenedCementChemistry_GetLog10Ksp(HCC) \
        ((HCC)->log10solubilityproductconstant)
#define HardenedCementChemistry_GetCSHCurves(HCC) \
        ((HCC)->cshcurves)
#define HardenedCementChemistry_GetCementSolutionChemistry(HCC) \
        ((HCC)->csc)
#define HardenedCementChemistry_GetLog10SaturationIndex(HCC) \
        ((HCC)->log10saturationindex)
#define HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(HCC) \
        ((HCC)->curveofcalciumsiliconratioincsh)
#define HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(HCC) \
        ((HCC)->curveofwatersiliconratioincsh)
#define HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(HCC) \
        ((HCC)->curveofsaturationindexofsh)



/* Macro for the electric potential
 * --------------------------------*/
#define HardenedCementChemistry_GetElectricPotential(HCC) \
        (CementSolutionChemistry_GetElectricPotential(HardenedCementChemistry_GetCementSolutionChemistry(HCC)))




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
#define HardenedCementChemistry_GetInput(HCC,U) \
        (HardenedCementChemistry_GetPrimaryVariable(HCC)[HardenedCementChemistry_##U])



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

#define HardenedCementChemistry_GetSaturationIndexOf(HCC,A) \
        (HardenedCementChemistry_GetSaturationIndex(HCC)[HardenedCementChemistry_S_##A])

#define HardenedCementChemistry_GetLog10SaturationIndexOf(HCC,A) \
        (HardenedCementChemistry_GetLog10SaturationIndex(HCC)[HardenedCementChemistry_S_##A])

/* Saturation index of CH */
#define HardenedCementChemistry_GetSaturationIndexOfCH(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,CH)
        
/* Saturation index of CC */
#define HardenedCementChemistry_GetSaturationIndexOfCC(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,CC)
        
/* Saturation index of SH */
#define HardenedCementChemistry_GetSaturationIndexOfSH(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,SH)
        
/* Saturation index of CSH2 */
#define HardenedCementChemistry_GetSaturationIndexOfCSH2(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,CSH2)

/* Saturation index of AH3 */
#define HardenedCementChemistry_GetSaturationIndexOfAH3(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,AH3)

/* Saturation index of AFm */
#define HardenedCementChemistry_GetSaturationIndexOfAFm(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,AFm)

/* Saturation index of AFt */
#define HardenedCementChemistry_GetSaturationIndexOfAFt(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,AFt)

/* Saturation index of C3AH6 */
#define HardenedCementChemistry_GetSaturationIndexOfC3AH6(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,C3AH6)

/* Saturation index of C2AH8 */
#define HardenedCementChemistry_GetSaturationIndexOfC2AH8(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,C2AH8)

/* Saturation index of CAH10 */
#define HardenedCementChemistry_GetSaturationIndexOfCAH10(HCC) \
        HardenedCementChemistry_GetSaturationIndexOf(HCC,CAH10)





/* Macros for variables
 * --------------------*/
#define HardenedCementChemistry_NbOfVariables  (3)

#define HardenedCementChemistry_X_CSH         (0)
#define HardenedCementChemistry_Z_CSH         (1)
//#define HardenedCementChemistry_V_CSH         (2)

/* Ca/Si ratio in C-S-H */
#define HardenedCementChemistry_GetCalciumSiliconRatioInCSH(HCC) \
        (HardenedCementChemistry_GetVariable(HCC)[HardenedCementChemistry_X_CSH])

/* Water/Si ratio in C-S-H */
#define HardenedCementChemistry_GetWaterSiliconRatioInCSH(HCC) \
        (HardenedCementChemistry_GetVariable(HCC)[HardenedCementChemistry_Z_CSH])

/* Molar volume of C-S-H */
/*
#define HardenedCementChemistry_GetMolarVolumeOfCSH(HCC) \
       (HardenedCementChemistry_GetVariable(HCC)[HardenedCementChemistry_V_CSH])
*/



/* Macros for Aqueous concentrations
 * ---------------------------------*/
#define HardenedCementChemistry_NbOfSpecies \
        CementSolutionChemistry_NbOfSpecies
        
#define HardenedCementChemistry_GetAqueousConcentration(HCC) \
        CementSolutionChemistry_GetConcentration(HardenedCementChemistry_GetCementSolutionChemistry(HCC))

#define HardenedCementChemistry_GetLogAqueousConcentration(HCC,CPD) \
        CementSolutionChemistry_GetLogConcentration(HardenedCementChemistry_GetCementSolutionChemistry(HCC))
        
#define HardenedCementChemistry_GetAqueousConcentrationOf(HCC,CPD) \
        CementSolutionChemistry_GetConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(HCC),CPD)

#define HardenedCementChemistry_GetLogAqueousConcentrationOf(HCC,CPD) \
        CementSolutionChemistry_GetLogConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(HCC),CPD)
        

/* Synonyms */
#define HardenedCementChemistry_NbOfConcentrations \
        HardenedCementChemistry_NbOfSpecies


/* Other macros
 * ------------*/
/* Liquid mass density */
#define HardenedCementChemistry_GetLiquidMassDensity(HCC) \
        CementSolutionChemistry_GetLiquidMassDensity(HardenedCementChemistry_GetCementSolutionChemistry(HCC))

/* Liquid charge density */
#define HardenedCementChemistry_GetLiquidChargeDensity(HCC) \
        CementSolutionChemistry_GetChargeDensity(HardenedCementChemistry_GetCementSolutionChemistry(HCC))

/* Ionic strength of aqueous phase*/
#define HardenedCementChemistry_GetIonicStrength(HCC) \
        CementSolutionChemistry_GetIonicStrength(HardenedCementChemistry_GetCementSolutionChemistry(HCC))

/* Element aqueous concentrations */
#define HardenedCementChemistry_GetElementAqueousConcentration(HCC) \
        CementSolutionChemistry_GetElementConcentration(HardenedCementChemistry_GetCementSolutionChemistry(HCC))
        
#define HardenedCementChemistry_GetElementAqueousConcentrationOf(HCC,A) \
        CementSolutionChemistry_GetElementConcentrationOf(HardenedCementChemistry_GetCementSolutionChemistry(HCC),A)






/* Macros for constants
 * --------------------*/
#define HardenedCementChemistry_NbOfConstants (2)

#define HardenedCementChemistry_A_CO2_EQ      (0)
#define HardenedCementChemistry_A_H2SO4_EQ    (1)


/* Equilibrium CO2 concentration */
#define HardenedCementChemistry_GetLog10EquilibriumCO2Activity(HCC) \
        (HardenedCementChemistry_GetConstant(HCC)[HardenedCementChemistry_A_CO2_EQ])

/* Equilibrium H2SO4 concentration */
#define HardenedCementChemistry_GetLog10EquilibriumH2SO4Activity(HCC) \
        (HardenedCementChemistry_GetConstant(HCC)[HardenedCementChemistry_A_H2SO4_EQ])




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

#define HardenedCementChemistry_GetLog10SolubilityProductConstantOf(HCC,CPD) \
        (HardenedCementChemistry_GetLog10Ksp(HCC)[HardenedCementChemistry_K_##CPD])



/* Macros for the resolution of the systems
 * ----------------------------------------*/
#define HardenedCementChemistry_ComputeSystem(HCC,SYS) \
        (HardenedCementChemistry_ComputeSystem_##SYS(HCC))

#define HardenedCementChemistry_SolveElectroneutrality(HCC) \
        (CementSolutionChemistry_SolveElectroneutrality(HardenedCementChemistry_GetCementSolutionChemistry(HCC)))

#define HardenedCementChemistry_SolveExplicitElectroneutrality(HCC) \
        (CementSolutionChemistry_SolveExplicitElectroneutrality(HardenedCementChemistry_GetCementSolutionChemistry(HCC)))



#include "Curves.h"

/* Macros for the CSH curves */
#define HardenedCementChemistry_GetCSHCurve(HCC) \
        (Curves_GetCurve(HardenedCementChemistry_GetCSHCurves(HCC)))

#define HardenedCementChemistry_SetDefaultCurveOfCalciumSiliconRatioInCSH(HCC) \
        do { \
          HardenedCementChemistry_GetCurveOfCalciumSiliconRatioInCSH(HCC) = HardenedCementChemistry_GetCSHCurve(HCC) + (0) ; \
        } while(0) 

#define HardenedCementChemistry_SetDefaultCurveOfWaterSiliconRatioInCSH(HCC) \
        do { \
          HardenedCementChemistry_GetCurveOfWaterSiliconRatioInCSH(HCC) = HardenedCementChemistry_GetCSHCurve(HCC) + (1) ; \
        } while(0) 

#define HardenedCementChemistry_SetDefaultCurveOfSaturationIndexOfSH(HCC) \
        do { \
          HardenedCementChemistry_GetCurveOfSaturationIndexOfSH(HCC) = HardenedCementChemistry_GetCSHCurve(HCC) + (2) ; \
        } while(0) 



#include "CementSolutionChemistry.h"

struct HardenedCementChemistry_s {
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
