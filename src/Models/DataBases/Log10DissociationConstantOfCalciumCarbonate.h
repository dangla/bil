#ifndef LOG10DISSOCIATIONCONSTANTOFCALCIUMCARBONATE_H
#define LOG10DISSOCIATIONCONSTANTOFCALCIUMCARBONATE_H


extern void Log10DissociationConstantOfCalciumCarbonate_Print(double) ;



#define Log10DissociationConstantOfCalciumCarbonate(R,T) (Log10DissociationConstantOfCalciumCarbonate_##R(T))


/* Refs. 
 * L. Niel PLUMMER, Eurybiades BUSENBERG The solubilities of calcite, aragonite and vaterite in CO2-H2O solutions between 0 and 90°C, and an evaluation of the aqueous model for the system CaCO3-CO2-H2O, Geochimica et Cosmochimica Acta, 46:1011-1040, 1982.
* */

/*
 * Cement chemistry notation:
 * H  = H2O  ; K  = K2O  ; N  = Na2O  ;
 * C  = CaO  ; S  = SiO2 ; A  = Al2O3 ; 
 * C' = CO2  ; S' = SO3  ; F  = Fe2O3 ;
 */

#include "TemperatureDependenceOfLog10EquilibriumConstant.h"

/* Calcium Carbonate (CC'): 
 * CC' = Ca[2+] + CO3[2-] */
 
/* Calcite */
#define Log10DissociationConstantOfCalciumCarbonate_Calcite__Ca_CO3(T) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-8.45255,-7.7993E-02,2.839319E+03,7.1595E+01,0)

/* Aragonite */
#define Log10DissociationConstantOfCalciumCarbonate_Aragonite__Ca_CO3(T) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-8.30501,-7.7993E-02,2.903293E+03,7.1595E+01,0)

/* Vaterite */
#define Log10DissociationConstantOfCalciumCarbonate_Vaterite__Ca_CO3(T) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-7.87225,-7.7993E-02,3.074688E+03,7.1595E+01,0)

#endif
