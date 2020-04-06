#ifndef LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_Cemdata18_H
#define LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_Cemdata18_H


#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18(R,T) \
        (Log10DissociationConstantOfCementHydrationProduct_Cemdata18_##R(T))


/* TO BE UPDATED */


/*
 * Cement chemistry notation:
 * H  = H2O  ; K  = K2O  ; N  = Na2O  ;
 * C  = CaO  ; S  = SiO2 ; A  = Al2O3 ; 
 * c  = CO2  ; s  = SO3  ; F  = Fe2O3 ;
 */

#include "TemperatureDependenceOfLog10EquilibriumConstant.h"
#include "PiecewiseLinearTemperatureDependence.h"

/* Three-term approximation */
#define TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,DrH,DrS,DrCp) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,DrH,DrS,DrCp)



/* Calcium Hydroxide: 
 * CH + 2H[+] = Ca[2+] + 2H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_CH_2H__Ca_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-130156,0,32.3)


/* Silica: 
 * S = SiO2(aq) */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_S__SiO2(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,21386,0,0)
        
        
/* CSH: CxSyHz = xCH + yS + (z-x)H  (CSHQ model from Kulik,2011) */
/* TobH: x = 2/3 , y = 1 , z = 3/2 */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_TobH__2o3CH_S_1o3H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,43.35e3,65,24)
/* TobD: x = 5/6 , y = 2/3 , z = 11/6 */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_TobD__5o6CH_2o3S_11o6H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,33.68e3,45,11)
/* JenH: x = 4/3 , y = 1 , z = 13/6 */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_JenH__4o3CH_S_13o6H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,51.73e3,68,16)
/* JenD: x = 3/2 , y = 2/3 , z = 5/2 */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_JenD__3o2CH_2o3S_5o2H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,35.5e3,49,3.4)
 



/* Al-Ettringite32 (AFt32 = C3A.(Cs)3.H32) */
/* Ca6Al2(SO4)3(OH)12(H2O)26 + 4H[+] = 2AlO2[-] + 6Ca[2+] + 3S6O4[2-] + 34H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFt32_4H__2AlO2_6Ca_3S6O4_34H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-23594,134.5,-694.0)

/* Al-Ettringite30 (AFt30 = C3A.(Cs)3.H30) */
/* Ca6Al2(SO4)3(OH)12(H2O)24 + 4H[+] = 2AlO2[-] + 6Ca[2+] + 3S6O4[2-] + 32H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFt30_4H__2AlO2_6Ca_3S6O4_32H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-36633,102.3,-764.6)


/* Al-Ettringite13 (AFt13 = C3A.(Cs)3.H13) */
/* Ca6Al2(SO4)3(OH)12(H2O)7 + 4H[+] = 2AlO2[-] + 6Ca[2+] + 3S6O4[2-] + 15H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFt13_4H__2AlO2_6Ca_3S6O4_15H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-596572,-1254.6,-1364.4)


/* Al-Ettringite9 (AFt9 = C3A.(Cs)3.H9) */
/* Ca6Al2(SO4)3(OH)12(H2O)3 + 4H[+] = 2AlO2[-] + 6Ca[2+] + 3S6O4[2-] + 11H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFt9_4H__2AlO2_6Ca_3S6O4_11H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-339631,-220.3,-1505.6)



/* Monosulfoalumiate12 (AFm12 = C3A.Cs.H12 */
/* Ca4Al2SO10(H2O)12 + 4H[+] = 2AlO2[-] + 4Ca[2+] + S6O4[2-] + 14H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFm12_4H__2AlO2_4Ca_S6O4_14H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-176816,-80.2,-381.2)


/* Monosulfoalumiate14 (AFm14 = C3A.Cs.H14 */
/* Ca4Al2SO10(H2O)14 + 4H[+] = 2AlO2[-] + 4Ca[2+] + S6O4[2-] + 16H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFm14_4H__2AlO2_4Ca_S6O4_16H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-185449,-109.6,-310.6)


/* Monosulfoalumiate16 (AFm16 = C3A.Cs.H16 */
/* Ca4Al2SO10(H2O)16 + 4H[+] = 2AlO2[-] + 4Ca[2+] + S6O4[2-] + 18H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AFm16_4H__2AlO2_4Ca_S6O4_18H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-148473,16.1,-246.2)
        
        


/* Monocarbonate: 
 * C3A.Cc.H11    = 4Ca[2+] + 2Al(OH)4[-] + CO3[2-] + 11H2O */


/* Gypsum: CsH2
 * CaSO4(H2O)2 = Ca[2+] + SO4[2-] + 2H2O  (solubility of gypsum in water is 15 mmol/L) */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_CsH2__2Ca_SO4_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,-1167,-91.6,-332.5)


/* Gibbsite: 0.5AH3 */
/* Al(OH)3  = AlO2[-] + H[+] + H2O */
#define Log10DissociationConstantOfCementHydrationProduct_Cemdata18_AlO3H3__AlO2_H_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_Cemdata18(T,77269,-30.4,-66.8)





/* to be completed */

/* Al-Katoite (Hydrogrossular): */
/* C3AH6 = 3Ca[2+] + 2Al[3+] + 12OH[-] (Ref (1])*/
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6__3Ca_2Al_12OH(T) \
        (xx)

/* C3AH6 = 3Ca[2+] + 2Al(OH)4[-] + 4OH[-] (ref [4]) */
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6__3Ca_2AlO4H4_4OH(T) \
        (-20.5)

/* C3AH6 + 12H[+] = 3Ca[2+] + 2Al[+] + 12H2O (ref [6]) */
#define Log10DissociationConstantOfCementHydrationProduct_C3AH6_12H__3Ca_2Al_12H2O(T) \
        (80.32)
//(PiecewiseLinearTemperatureDependence(T,89.6883,80.32,69.5665,59.7469,50.0831,42.4618,36.2972,31.2083))

/* C2AH8 = 2Ca[2+] + 2Al(OH)4[-] + 2OH[-] + 3H2O (ref [4]) */
#define Log10DissociationConstantOfCementHydrationProduct_C2AH8__2Ca_2AlO4H4_2OH_3H2O(T) \
        (-13.56)

/* CAH10 = Ca[2+] + 2Al(OH)4[-] + 6H2O (ref [4]) */
#define Log10DissociationConstantOfCementHydrationProduct_CAH10__Ca_2AlO4H4_6H2O(T) \
        (-7.6)


/* Anhydrite (metastable in presence of water):
 * Cs = Ca[2+] + SO4[2-] (solubility of anhydrite is 24 mmol/L) */
#define Log10DissociationConstantOfCementHydrationProduct_CS__Ca_SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-3.239577517,0,0,0,0) /* To be checked */


/* Magnesium carbonate (magnesite):
 * MgCO3 */
 
 
/* Hydrogranet
 * mixture of C3A.S(3-x).H(2x) and C3F.S(3-x).H(2x) (with 0<x<3) 
 * i.e. C3.Ay.F(1-y).S(3-x).H(2x) */





/* References:
 * 
 * B. Lothenbach, D.A Kulik, T. Matschei, M. Balonis, L. Baquerizo, B. Dilnesa, G.D. Miron, R.J. Myers,
 * Cemdata18: A chemical thermodynamic database for hydrated Portland cements and alkali-activated materials,
 * Cement and Concrete Research 115 (2019) 472–506.
 * 
 * D.A. Kulik,
 * Improving the structural consistency of C-S-H solid solution thermodynamic models,
 * Cement and Concrete Research 41 (2011) 477-495.
 */
 

#endif
