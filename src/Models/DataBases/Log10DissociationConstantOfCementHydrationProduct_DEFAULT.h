#ifndef LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_DEFAULT_H
#define LOG10DISSOCIATIONCONSTANTOFCEMENTHYDRATIONPRODUCT_DEFAULT_H


/*
 * Cement chemistry notation:
 * H  = H2O  ; K  = K2O  ; N  = Na2O  ;
 * C  = CaO  ; S  = SiO2 ; A  = Al2O3 ; 
 * c  = CO2  ; s  = SO3  ; F  = Fe2O3 ;
 */

#include "TemperatureDependenceOfLog10EquilibriumConstant.h"
#include "PiecewiseLinearTemperatureDependence.h"



/* Calcium Hydroxide: 
 * CH = Ca[2+] + 2OH[-] */
#define Log10DissociationConstantOfCementHydrationProduct_CH__Ca_2OH(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-5.14721,-1.0139684E-01,3.3421032E+04,2.04489E+02,-2.239338E+06)


/* Silica: 
 * S + 2H2O = H4SiO4 */
#define Log10DissociationConstantOfCementHydrationProduct_S_2H2O__H4SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-2.75894,0,-7.623113E+02,0,0)


/* to be completed */

/* Al-Ettringite (AFt):  (ref [4]) */
/* C3A.(Cs)3.H32 = 6Ca[2+] + 2Al(OH)4[-] + 3SO4[2-] + 4OH[-] + 26H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFt__6Ca_2AlO4H4_3SO4_4OH_26H2O(T) \
        (-44.9)

/* C3A.(Cs)3.H32 = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 12OH[-] + 26H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFt__6Ca_2Al_3SO4_12OH_26H2O(T) \
        (xx)

/* C3A.(Cs)3.H32 + 12H[+] = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 38H2O (ref. [5]) */
#define Log10DissociationConstantOfCementHydrationProduct_AFt_12H__6Ca_2Al_3SO4_38H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-576.29,0,44841.95,195.,0)


/* Monosulfoalumiate (AFm): (ref [4]) */
/* C3A.Cs.H12    = 4Ca[2+] + 2Al(OH)4[-] + SO4[2-]  + 4OH[-] + 6H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFm__4Ca_2AlO4H4_SO4_4OH_6H2O(T) \
        (-29.26)

/* C3A.Cs.H12    = 4Ca[2+] + 2Al[3+] + SO4[2-]  + 12OH[-] + 6H2O */
#define Log10DissociationConstantOfCementHydrationProduct_AFm__4Ca_2Al_SO4_12OH_6H2O(T) \
        (xx)

/* C3A.Cs.H12 + 12H[+]  = 4Ca[2+] + 2Al[3+] + SO4[2-] + 18H2O (ref. [5]) */
#define Log10DissociationConstantOfCementHydrationProduct_AFm_12H__4Ca_2Al_SO4_18H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,59.77,0,24477.17,-28.09,0)
        

/* Friedel's salt (ref [7]) */
/* C3AH10.CaCl2 = C4ACl2H10 = Ca4.Al2.Cl2.12(OH)·4(H2O) 
 * Ca4Al2Cl2(OH)12·4H2O = 4Ca[2+] + 2Al(OH)4[-] + 2Cl[-] + 4OH[-] + 4H2O
 * */
#define Log10DissociationConstantOfCementHydrationProduct_FriedelSalt__4Ca_2AlO4H4_2Cl_4OH_4H2O(T) \
        (-27.27)
#define Log10DissociationConstantOfCementHydrationProduct_C3AH10CaCl2__4Ca_2AlO4H4_2Cl_4OH_4H2O(T) \
        (-27.27)
        
        
/* Kuzel's salt */
/* C3AH12.0.5(Cs).0.5(CaCl2) = Ca4.Al2.Cl(SO4)0.5.12(OH)·6(H2O) 
 * Ca4Al2Cl(SO4)0.5(OH)12·6H2O = 4Ca[2+] + 2Al(OH)4[-] + Cl[-] + 0.5SO4[2-] + 4OH[-] + 6H2O*/
#define Log10DissociationConstantOfCementHydrationProduct_KuzelSalt__4Ca_2AlO4H4_Cl_halfSO4_4OH_6H2O(T) \
        (-28.53)
#define Log10DissociationConstantOfCementHydrationProduct_2KuzelSalt__8Ca_4AlO4H4_2Cl_SO4_8OH_12H2O(T) \
        (-57.06)


/* Monocarbonate: 
 * C3A.Cc.H11    = 4Ca[2+] + 2Al(OH)4[-] + CO3[2-] + 11H2O */


/* Gypsum:
 * CsH2 = Ca[2+] + SO4[2-] + 2H2O  (solubility of gypsum in water is 15 mmol/L) */
#define Log10DissociationConstantOfCementHydrationProduct_CSH2__Ca_SO4_2H2O(T) \
        (-4.58)


/* Gibbsite: */
/* AH3  = 2Al[3+] + 6OH[-] */
#define Log10DissociationConstantOfCementHydrationProduct_AH3__2Al_6OH(T) \
        (-68)

/* AH3 + 2OH[-]  = 2Al(OH)4[-] (ref [4]) */
#define Log10DissociationConstantOfCementHydrationProduct_AH3_2OH__2AlO4H4(T) \
        (-1.12)


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
 * [1] B. Lotenbach, L. Pelletier-Chaignat, F. Winnefeld, 
 * Stability in the sytem CaO-Al2O3-H2O. 
 * Cement and Concrete Research,42:1621-1634,2012.
 * 
 * [2] T. Thoenen and D. Kulik. 
 * Nagra/psi chemical thermodynamic data base 01/01 for the gem-selektor (v. 2-psi) geochemical modeling code: Release 28-02-03. Technical report, PSI Technical Report TM-44-03-04 about the GEMS version of Nagra/PSI chemical thermodynamic database 01/01, 2003.
 * 
 * [3] B. Lothenbach and F. Winnefeld, 
 * Thermodynamic modelling of the hydration of Portland cement. 
 * Cement and Concrete Research 36:209-226,2006.
 * 
 * [4] Cemdata14
 * 
 * [5] TRAN Van Quan, Contribution à la compréhension des mécanismes de dépassivation des armatures d'un béton exposé à l'eau de mer: théorie et modélisation thermodynamique, PhdThesis, Univ. Nantes, 2016.
 * 
 * [6] Thermochimie, ANDRA
 * 
 * [7] Lothenbach, Kulik, Matschei, Balonis, Baquerizo, Dilnesa, Miron, Myers
 * Cemdata18: A chemical thermodynamic database for hydrated Portland cements and 
 * alkali-activated materials, 
 * CCR 115 (2019) 472-506. 
 */
 

#endif
