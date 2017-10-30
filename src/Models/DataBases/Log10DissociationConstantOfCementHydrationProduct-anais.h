#ifndef LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_H
#define LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_H

#define Log10DissociationConstantOfCementHydrationProduct(R,T) (Log10DissociationConstantOfCementHydrationProduct_##R(T))


/*
 * Cement chemistry notation:
 * H  = H2O  ; K  = K2O  ; N  = Na2O  ;
 * C  = CaO  ; S  = SiO2 ; A  = Al2O3 ; 
 * C' = CO2  ; S' = SO3  ; F  = Fe2O3 ;
 */

#include "TemperatureDependenceOfLog10EquilibriumConstant.h"

/* Calcium Hydroxide: 
 * CH = Ca[2+] + 2OH[-] (ref [3]) */
#define Log10DissociationConstantOfCementHydrationProduct_CH__Ca_2OH(T) (22.79937) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-11.2992307115,0,7303.68926279255,3.8839330279,0)

/* Calcium Hydroxide: 
 * CH + 2H[+] = Ca[2+] + 2H2O[-] thermodem */
#define Log10DissociationConstantOfCementHydrationProduct_CH_H__Ca_H2O(T) (22.81) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-2.8493E02,-4.4711E-02,2.1380E04,1.0421E02,-7.5425E05)


/* Silica: 
 * S + 2H2O = H4SiO4 
#define Log10DissociationConstantOfCementHydrationProduct_S_2H2O__H4SiO4(T) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-2.75894,0,-7.623113E+02,0,0)*/

/* SiO2 + OH[-] + H2O = H4SiO4 ref[3]*/ 
#define Log10DissociationConstantOfCementHydrationProduct_SiO2_OH_H2O__H4SiO4(T) (1.48) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-2.1418,0,664.055,0.6,0)


/* to be completed */

/* Al-Ettringite (AFt):  (ref [3])*/
/* C3A.(CS')3.H32 = 6Ca[2+] + 2Al(OH)4[-] + 3SO4[2-] + 4OH[-] + 26H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFt__6Ca_2AlO4H4_3SO4_4OH_26H2O(T) (-44.9085) TemperatureDependenceOfLog10EquilibriumConstant_293 (T,529.277779378424,0,-34439.99676380200,-185.338642966,0)


/* Al-Ettringite (AFt):  thermodem*/
/* C3A.(CS')3.H32 + 12H[+] = 6Ca[2+] + 2Al[3+] + 3SO4[2-]  + 38H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFt_12H__6Ca_2Al_3SO4_38H2O(T) (57.01) TemperatureDependenceOfLog10EquilibriumConstant_293 (T,-6.6746E03,-1.0474,3.7879E05,2.4266E03,-2.0519E07)


/* C3A.(CS')3.H32 = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 12OH[-] + 26H2O
#define Log10DissociationConstantOfCementHydrationProduct_AFt__6Ca_2Al_3SO4_12OH_26H2O(T)    ()*/


/* Monosulfate (AFm): (ref [3])*/
/* C3A.CS'.H12    = 4Ca[2+] + 2Al(OH)4[-] + SO4[2-]  + 4OH[-] + 6H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFm__4Ca_2AlO4H4_SO4_4OH_6H2O(T)  (-29.2628) TemperatureDependenceOfLog10EquilibriumConstant_293(T,404.92625454716,0,-21017,60772698490,-146.971980391,0)

/* Monosulfate (AFm): thermodem*/
/* C3A.CS'.H12  + 12H[+]  = 4Ca[2+] + 2Al[3+] + SO4[2-]  + 18H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFm_12H__4Ca_2Al_SO4_18H2O(T)  (73.09) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-3.5426E03,-5.7055E-01,2.1351E05,1.2875E03,-1.0329E07)


/* C3A.CS'.H12    = 4Ca[2+] + 2Al[3+] + SO4[2-]  + 12OH[-] + 6H2O 
#define Log10DissociationConstantOfCementHydrationProduct_AFm__4Ca_2Al_SO4_12OH_6H2O(T)    ()*/


/* Monocarboaluminate: 
 * C3A.CC'.H11    = 4Ca[2+] + 2Al(OH)4[-] + CO3[2-] + 11H2O */


/* Gypsum:
 * CS'H2 = Ca[2+] + SO4[2-] + 2H2O  (solubility of gypsum in water is 15 mmol/L) ref [3]*/
#define Log10DissociationConstantOfCementHydrationProduct_CSH2__Ca_SO4_2H2O(T) (-4.58147) TemperatureDependenceOfLog10EquilibriumConstant_293(T,111.52942046684,0,-5116.92285811085,-39.9882855394,0)

/* Gypsum:
 * CS'H2 = Ca[2+] + SO4[2-] + 2H2O  (solubility of gypsum in water is 15 mmol/L) thermodem*/
#define Log10DissociationConstantOfCementHydrationProduct_CSH2__Ca_SO4_2H2O(T) (-4.60) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.6202E03,-2.5723E-01,8.9151E04,5.8738E02,-5.3473E06)


/* Gibbsite: */
/* AH3  = 2Al[3+] + 6OH[-] 
#define Log10DissociationConstantOfCementHydrationProduct_AH3__2Al_6OH(T) ()*/

/* AH3  = 2Al(OH)4[-] - OH[-] * ref [3]*/
#define Log10DissociationConstantOfCementHydrationProduct_AH3_2OH__2AlO4H4(T) (0.239427) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-48.08001615409,0,1546.32359095493,17.4322960704,0)

/* AH3 + 3H[+] = Al[3+] + 3H2O * thermodem*/
#define Log10DissociationConstantOfCementHydrationProduct_AH3_3H__Al_3H2O(T) (7.74) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-4.9375E02,-8.0900E-02,2.9714E04,1.7790E02,-1.2677E06)



/* Aluminium Hydrates*/
/* C3AH6 = 3Ca[2+] + 2Al[3+] + 12OH[-] (Ref (1])
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6__3Ca_2Al_12OH(T) ()*/

/* C3AH6 = 3Ca[2+] + 2Al(OH)4[-] + 4OH[-] (ref [3])*/
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6__3Ca_2AlO4H4_4OH(T) (-20.8411) TemperatureDependenceOfLog10EquilibriumConstant_293(T,291.27318390754,0,-13720.16812376770,-107.5364633634,0)

/* C3AH6 + 12H[+]= 3Ca[2+] + 2Al[3+] + 12H2O thermodem*/
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6_12H__3Ca_2Al_12H2O(T) (80.33) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.7858E03,-2.8804E-01,1.1971E05,6.4758E02,-4.6117E06)


/* C2AH8 = 2Ca[2+] + 2Al[3+] + 10OH[-] + 3H2O (ref [3])
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6__3Ca_2AlO4H4_4OH(T) (-20.5)*/

/* C2AH8 = 2Ca[2+] + 2Al(OH)4[-] + 2OH[-] + 3H2O (ref [3]) */
#define Log10DissociationConstantOfCementHydrationProduct_C2AH8__2Ca_2AlO4H4_2OH_3H2O(T) (-13.5622) TemperatureDependenceOfLog10EquilibriumConstant_293(T,154.48530223981,0,-8994.12914254083,-55.7212631893,0)

/* C2AH8 + 10H[+]= 2Ca[2+] + 2Al[3+] + 13H2O thermodem */
#define Log10DissociationConstantOfCementHydrationProduct_C2AH8_10H__2Ca_2Al_13H2O(T) (59.72) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.7347E03,-2.3802E-01,1.1041E05,6.2322E02,-4.1858E06)

/* CAH10 = Ca[2+] + 2Al(OH)4[-] + 6H2O (ref [3]) */
#define Log10DissociationConstantOfCementHydrationProduct_CAH10__Ca_2AlO4H4_6H2O(T) (-7.50282) TemperatureDependenceOfLog10EquilibriumConstant_293(T,19.33995799190,0,-4170.78109558596,-5.1912860770,0)
	
/* C4AH13 = 4Ca[2+] + 2Al(OH)4[-] + 6OH[-] + 6H2O (ref [3])*/
#define Log10DissociationConstantOfCementHydrationProduct_C4AH13__4Ca_2AlO4H4_6OH_6H2O(T) (-25.4033) TemperatureDependenceOfLog10EquilibriumConstant_293(T,407.04608084006,0,-21037.51384617330,-146.2503502912,0)

/* C4AH13 + 14H[+] = 4Ca[2+] + 2Al[3+] +  20H2O thermodem*/
#define Log10DissociationConstantOfCementHydrationProduct_C4AH13_14H__4Ca_2Al_20H2O(T) (103.67) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-2.1675E03,-3.2686E-01,1.4525E05,7.8650E02,-5.7627E06)

/* Anhydrite (metastable in presence of water):
 * CS' = Ca[2+] + SO4[2-] (solubility of anhydrite is 24 mmol/L) (ref[3])*/
#define Log10DissociationConstantOfCementHydrationProduct_CS__Ca_SO4(T) (-4.35754)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,-131.22651283827,0,-5228.78943438100,-47.7070807540,0)

/* Anhydrite (metastable in presence of water):
 * CS' = Ca[2+] + SO4[2-] (solubility of anhydrite is 24 mmol/L) thermodem*/
#define Log10DissociationConstantOfCementHydrationProduct_CS__Ca_SO4(T) (-4.44)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.6181E03,-2.6204E-01,8.9585E04,5.8663E02,-5.3589E06)


	   
/* References
 * ----------
 *
 * [1] B. Lothenbach, F. Winnefeld, Thermodynamic modelling of the hydration of Portland cement.
 * Cement and Concrete Research, 36:209-226, 2006.
 * 
 * [2] T. Thoenen and D. Kulik. Nagra/psi chemical thermodynamic data base 01/01 for the gem-selektor (v. 2-psi) geochemical modeling code: Release 28-02-03. Technical report, PSI Technical Report TM-44-03-04 about the GEMS version of Nagra/PSI chemical thermodynamic database 01/01, 2003.
 * 
 * [3] D.Jacques, Benchmarking of the cement model and detrimental chemical reactions including temperature dependent parameters, Project near surface disposal of category A waste at Dessel, 2009.
 */
#endif
