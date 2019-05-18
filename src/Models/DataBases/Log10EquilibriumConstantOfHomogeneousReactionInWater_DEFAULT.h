#ifndef LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_DEFAULT_H
#define LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_DEFAULT_H



#include "TemperatureDependenceOfLog10EquilibriumConstant.h"
#include "PiecewiseLinearTemperatureDependence.h"

#define LOGK(...)   (PiecewiseLinearTemperatureDependence(__VA_ARGS__))

#if 0
#define LINEAR(T,T1,T2,L1,L2)    ((((T)-(T1))*(L2) + ((T2)-(T))*(L1))/((T2)-(T1)))

#define LOGK(T,L1,L2,L3,L4,L5,L6,L7,L8)    LOGKC(((T)-273.15),L1,L2,L3,L4,L5,L6,L7,L8)
#define LOGKC(T,L1,L2,L3,L4,L5,L6,L7,L8)   (((T) < 0  ) ? (L1)                    : LOGKGT0(T,L1,L2,L3,L4,L5,L6,L7,L8))
#define LOGKGT0(T,L1,L2,L3,L4,L5,L6,L7,L8) (((T) < 25 ) ? LINEAR(T,0,25,L1,L2)    : LOGKGT25(T,L2,L3,L4,L5,L6,L7,L8))
#define LOGKGT25(T,L2,L3,L4,L5,L6,L7,L8)   (((T) < 60 ) ? LINEAR(T,25,60,L2,L3)   : LOGKGT60(T,L3,L4,L5,L6,L7,L8))
#define LOGKGT60(T,L3,L4,L5,L6,L7,L8)      (((T) < 100) ? LINEAR(T,60,100,L3,L4)  : LOGKGT100(T,L4,L5,L6,L7,L8))
#define LOGKGT100(T,L4,L5,L6,L7,L8)        (((T) < 150) ? LINEAR(T,100,150,L4,L5) : LOGKGT150(T,L5,L6,L7,L8))
#define LOGKGT150(T,L5,L6,L7,L8)           (((T) < 200) ? LINEAR(T,150,200,L5,L6) : LOGKGT200(T,L6,L7,L8))
#define LOGKGT200(T,L6,L7,L8)              (((T) < 250) ? LINEAR(T,200,250,L6,L7) : LOGKGT250(T,L7,L8))
#define LOGKGT250(T,L7,L8)                 (((T) < 300) ? LINEAR(T,250,300,L7,L8) : (L8))
#endif



/* Autoprotolysis of water */

/* H2O = H[+] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2O__H_OH(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-14.1733,-5.069842E-02,1.332300E+04,1.022445E+02,-1.119669E+06)


/* 
 * We define compounds of type I and II as follows: 
 * 
 * Compounds of type I are monoatomic ions and compounds formed 
 * with only one element other than Hydrogen and Oxygen.
 * 
 * Compounds of type II are formed with only two different 
 * elements other than Hydrogen and Oxygen. 
 */
 
 

/* 1. Chemical reactions involving compounds of type I. */


/* Calcium compounds */

/* Ca[2+] + OH[-] = CaOH[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_OH__CaOH(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.39331,5.069842E-02,-1.332300E+04,-1.022445E+02,1.119669E+06)

/* CaOH[+] = Ca[2+] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaOH__Ca_OH(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.39331,-5.069842E-02,1.332300E+04,1.022445E+02,-1.119669E+06)

/* Ca(OH)2[0] = Ca[2+] + 2OH[-] (to be removed because Ca(OH)2[0] seems not be there) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaO2H2__Ca_2OH(T) \
        (0)

/* Ca[2+] + H2O = CaOH[+] + H[+] */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H2O__CaOH_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,1.3130E02,2.1418E-02,-1.0190E04,-4.8225E01,2.7033E05)

/* CaOH[+] + H[+] = Ca[2+] + H2O */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaOH_H__Ca_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-1.3130E02,-2.1418E-02,1.0190E04,4.8225E01,-2.7033E05)


/* Silicon compounds */

/* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H3SiO4_OH__H2SiO4_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,0.691186,5.069842E-02,-1.5903388E+04,-1.022445E+02,1.119669E+06)

/* H3SiO4[-] + H2O = H4SiO4[0] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H3SiO4_H2O__H4SiO4_OH(T) \
        (-4.19)

/* H2SiO4[2-] + H2O = H3SiO4[-] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SiO4_H2O__H3SiO4_OH(T) \
        (-0.67)

/* H4SiO4 = H3SiO4[-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__H3SiO4_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-9.88883,0,-1.337205E+03,0,0)


/* Sodium compounds */

/* NaOH[0] = Na[+] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaOH__Na_OH(T) \
        (0.18)

/* NaOH[0] + H[+] = Na[+] + H2O */ /*  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaOH_H__Na_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-40.0224,0,4902.5816,15.3,0)

/* Na[+] + H2O = NaOH[0] + H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_H2O__NaOH_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,40.0224,0,-4902.5816,-15.3,0)


/* Potassium compounds */

/* KOH[0] = K[+] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KOH__K_OH(T) \
        (0.46)

/*KOH[0] + H[+] = K[+] + H2O */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KOH_H__K_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-55.77230945322,0,0,0,0)

/* K[+] + H2O = KOH[0] + H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_K_H2O__KOH_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,55.77230945322,0,0,0,0)



/* Carbon compounds:
 * Carbone dioxide (dissolved): CO2(aq) = CO2[0]
 * Carbonic acid: H2CO3[0]
 * Bicarbonate ion (or hydrogen carbonate ion): HCO3[-]
 * Carbonate ion: CO3[2-] 
 * CO2[0] + H2CO3[0] noted as H2CO3[*] */

/* CO2[0] + H2O = H2CO3[0] */   /* Hydration reaction of CO2 */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO2_H2O__H2CO3(T) \
        (-2.77)

/* H2CO3[0] = CO2[0] + H2O */   /* Hydration reaction of CO2 */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__CO2_H2O(T) \
        (2.77)

/* H2CO3[0] = HCO3[-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__HCO3_H(T) \
        (-3.6)

/* HCO3[-] + H2O = H2CO3[0] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H2O__H2CO3_OH(T) \
        (-10.4)

/* H2CO3[*] = HCO3[-] + H[+] */
/*
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__HCO3_H(T) \
          (-6.3)
*/

/* HCO3[-] + H2O = H2CO3[*] + OH[-] */
/*
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H2O__H2CO3_OH(T) \
          (-7.648)
*/

/* HCO3[-] = CO2[0] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3__CO2_OH(T) \
        (-7.61713)

/* HCO3[-] + H[+] = CO2[0] + H2O */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H__CO2_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,6.38287,6.091964E-02,-2.183437E+04,-1.268339E+02,1.684915E+06)

/* CO2[0] + H2O = HCO3[-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO2_H2O__HCO3_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-6.38287,-6.091964E-02,2.183437E+04,1.268339E+02,-1.684915E+06)

/* HCO3[-] = CO3[2-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3__H_CO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-10.3771,-3.252849E-02,5.151790E+03,3.892561E+01,-5.637139E+05)

/* CO3[2-] + H2O = HCO3[-] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO3_H2O__HCO3_OH(T) \
        (-3.67)


/* Sulfur compounds */ /* from EQ3/6 */

/* H2S[0] = HS[-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S__HS_H(T) \
        (-LOGK(T,7.4159,6.9877,6.6467,6.4827,6.496,6.6831,7.0225,7.5536))

/* HS[-] = S[2-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HS__S_H(T) \
        (-LOGK(T,13.7100,12.9351,12.0082,11.1018,10.1202,9.2545,8.4250,7.5568))

/* H2SO4[0] = HSO4[-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SO4__HSO4_H(T) \
        (6)

/* HSO4[-] = SO4[2-] + H[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO4__SO4_H(T) \
        (-LOGK(T,1.7193,1.9791,2.4371,3.0002,3.7234,4.4683,5.2633,6.1799))

/* HSO3[-] = SO3[2-] + H[+] */   /* Ref chess.tdb */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO3__SO3_H(T) \
        (-LOGK(T,7.1301,7.2054,7.4642,7.8631,8.4488,9.11,9.8607,10.7638))

/* SO3[2-] + 2H[+] = H2SO3 */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO3_H__H2SO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,1.2948E03,2.1816E-01,-7.3030E04,-4.6772E02,4.5780E06)

/* H2SO3 = SO3[2-] + 2H[+] */ /*  thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SO3__SO3_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-1.2948E03,-2.1816E-01,7.3030E04,4.6772E02,-4.5780E06)

/* S2O3[2-] + 2H[+] = H2S2O3 */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_S2O3_H__H2S2O3(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,1.4979E03,2.3814E-01,-8.4049E04,-5.4207E02,5.0380E06))

/* H2S2O3 = S2O3[2-] + 2H[+]  thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S2O3__S2O3_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-1.4979E03,-2.3814E-01,8.4049E04,5.4207E02,-5.0380E06)

/* S2O4[2-] + 2H[+] = H2S2O4 */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_S2O4_H__H2S2O4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,1.5238E03,2.4188E-01,-8.5504E04,-5.5134E02,5.1466E06)

/* H2S2O4 = S2O4[2-] + 2H[+] */ /*  thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S2O4__S2O4_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-1.5238E03,-2.4188E-01,8.5504E04,5.5134E02,-5.1466E06)


/* Aluminium compounds */

/* Al(OH)[2+] = Al[3+] + OH[-] */ /* (ref [1]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlOH__Al_OH(T) \
        (-14 + 4.957)

/* Al(OH)2[+] = Al[3+] + 2OH[-] */ /* (ref [1]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2H2__Al_2OH(T) \
        (-28 + 10.594)

/* Al(OH)3[0] = Al[3+] + 3OH[-] */ /* (ref [1]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO3H3__Al_3OH(T) \
        (-42 + 16.432)

/* Al(OH)4[-] = Al[3+] + 4OH[-] */ /* (ref [1]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4__Al_4OH(T) \
        (-56 + 22.879)


/* Al[3+] + 2H2O = Al(OH)2[+] + 2H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_2H2O__AlO2H2_2H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(32.5363,0,-6492.5404,-8.6,0)

/* Al(OH)2[+] + 2H[+] = Al[3+] + 2H2O */ /*  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2H2_2H__Al_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-32.5363,0,6492.5404,8.6,0)

/* Al[3+] + 4H2O = Al(OH)4[-] + 4H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_4H2O__AlO4H4_4H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(32.9602,0,-10341.6037,-8.5,0)

/* Al(OH)4[-] + 4H[+] = Al[3+] + 4H2O */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4_4H__Al_4H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-32.9602,0,10341.6037,8.5,0)

/* Al[3+] + 3H2O = Al(OH)3 + 3H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlO3H3_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(91.6640,0,-11712.9732,-27.8,0))

/* Al(OH)3 + 3H[+] = Al[3+] + 3H2O */ /*  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO3H3__Al_3H2O_3H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-91.6640,0,11712.9732,27.8,0)

/* Al[3+] + H2O = AlOH[2+] + H[+] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlOH_H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-34.5703,0,-869.7066,13.1,0)

/* AlOH[2+] + H[+] = Al[3+] + H2O */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlOH_H__Al_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(34.5703,0,869.7066,-13.1,0)

/* Al[3+] + 2H2O = AlO2[-] + 4H[+] */ /* thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_2H2O__AlO2_4H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-1.7805E02,-2.6890E02,1.8672E03,6.6833E01,-7.5044E05)

/* AlO2[-] + 4H[+] = Al[3+] + 2H2O */ /* thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2_4H__Al_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(1.7805E02,2.6890E02,-1.8672E03,-6.6833E01,7.5044E05)

/* Al[3+] + 2H2O = HAlO2 + 3H[+] */ /* thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_2H2O__HAlO2_3H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(3.4561E02,6.0311E-02,-2.5787E04,-1.2376E02,1.1292E04)

/* HAlO2 + 3H[+] = Al[3+] + 2H2O */ /* thermodem*/ 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HAlO2_3H__Al_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(-3.4561E02,-6.0311E-02,2.5787E04,1.2376E02,-1.1292E04)



/* 2. Chemical reactions involving compounds of type II. */
 
 
/* Calcium-Silicon compounds */
     
/* Ca[2+] + H2SiO4[2-] = CaH2SiO4 */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H2SiO4__CaH2SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,4.6,0,0,0,0)
     
/* CaH2SiO4 = Ca[2+] + H2SiO4[2-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH2SiO4__Ca_H2SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-4.6,0,0,0,0)
       
/* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H3SiO4__CaH3SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.2,0,0,0,0)
     
/* CaH3SiO4[+] = Ca[2+] + H3SiO4[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH3SiO4__Ca_H3SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.2,0,0,0,0)

 
/* Calcium-Carbon compounds */

/* Ca[2+] + HCO3[-]    = CaHCO3[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HCO3__CaHCO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.06756,3.129400E-01,-3.476505E+04,-4.787820E+02,0) 

/* CaHCO3[+] = Ca[2+] + HCO3[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaHCO3__Ca_HCO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.06756,-3.129400E-01,3.476505E+04,4.787820E+02,0) 

/* Ca[2+] + CO3[2-]    = CaCO3[0] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_CO3__CaCO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,3.1836,-2.9944401E-01,3.551275E+04,4.8581799E+02,0) 

/* CaCO3[0] = Ca[2+] + CO3[2-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaCO3__Ca_CO3(T) \
        TemperatureDependenceOfLog10EquilibriumConstant_293(T,-3.1836,2.9944401E-01,-3.551275E+04,-4.8581799E+02,0)


 
/* Sodium-Carbon compounds */

/* NaHCO3[0] = Na[+] + HCO3[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaHCO3__Na_HCO3(T) \
        (0.25) 

/* NaCO3[-] = Na[+] + CO3[2-] */ /* (ref [4]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaCO3__Na_CO3(T) \
        (-1.27) 


 
/* Sodium-Sulfur compounds */

/* NaSO4[-] = Na[+] + SO4[2-] */ /* (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaSO4__Na_SO4(T)  (-0.7) 
/* NaSO4[-] = Na[+] + SO4[2-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaSO4__Na_SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,67.9470,0,-2908.6947,-23.8,0)

/*  Na[+] + SO4[2-] =  NaSO4[-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_SO4__NaSO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-67.9470,0,2908.6947,23.8,0)
 

 
/* Potassium-Sulfur compounds */

/* KSO4[-] = K[+] + SO4[2-] */ /* (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KSO4__K_SO4(T)    (-0.85) 
/* KSO4[-] = K[+] + SO4[2-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KSO4__K_SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,72.9765,0,-3150.0583,-25.6,0)

/*  K[+] + SO4[2-] = KSO4[-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_K_SO4__KSO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-72.9765,0,3150.0583,25.6,0)



/* Calcium-Sulfur compounds */

/* Ca[2+] + HS[-]      = CaHS[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HS__CaHS(T) \
        (1.10585)

/* Ca[2+] + S[2-]      = CaS[0] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_S__CaS(T) \
        (3.544068)

/* Ca[2+] + SO4[2-] + H[+]    = CaHSO4[+]  (phreeqc) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4_H__CaHSO4(T) \
        (3.068)

/* Ca[2+] + HSO4[-]    = CaHSO4[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HSO4__CaHSO4(T) \
        (1.10585)

/* CaHSO4[+] = Ca[2+] + HSO4[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaHSO4__Ca_HSO4(T) \
        (-1.10585)

/* Ca[2+] + SO4[2-]    = CaSO4[0] (wrong values which should be changed with those of ref [3] below)
 *                                (kept because used in Yuan's models) */
/*
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4__CaSO4(T) \
        (3.146128036)
*/

/* CaSO4[0] = Ca[2+] + SO4[2-] */
/*
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T) \
        (-3.146128036)
*/

/* CaSO4[0] = Ca[2+] + SO4[2-] */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T)  (-LOGK(T,2.0713,2.1111,2.2647,2.5111,2.9101,3.4328,4.1424,5.1853))

/* Ca[2+] + SO4[2-] = CaSO4[0] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4__CaSO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-64.2534,0,2770.0342,23.1,0)

/* CaSO4[0] = Ca[2+] + SO4[2-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,64.2534,0,-2770.0342,-23.1,0)



/* Aluminium-Silicon compounds */

/* AlH3SiO4[2+] = Al[3+] + H3SiO4[-] (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlH3SiO4__Al_H3SiO4(T)          (-7.4)

/* Al(OH)4H4SiO4[-] + H2O = Al(OH)4[-] + H4SiO4[0] */ /* (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4H4SiO4__AlO4H4_H4SiO4(T)  (-3.6)


/* AlH3SiO4[2+] = Al[3+] + H3SiO4[-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlH3SiO4__Al_H3SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,-2206.3355,0,0)

/* Al[3+] + H3SiO4[-] = AlH3SiO4[2+] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H3SiO4__AlH3SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,2206.3355,0,0)

/* Al(OH)4H4SiO4[-] + H2O = Al(OH)4[-] + H4SiO4[0] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4H4SiO4_H2O__AlO4H4_H4SiO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,-1073.3948,0,0)

/*  Al(OH)4[-] + H4SiO4[0] = Al(OH)4H4SiO4[-] + H2O */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4_H4SiO4__AlO4H4H4SiO4_H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,1073.3948,0,0)



/* Aluminium-Sulfur compounds */

/* AlSO4[+] = Al[3+] + SO4[2-] */ /* (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlSO4__Al_SO4(T)    (-3.9)

/* Al(SO4)2[-] = Al[3+] + 2SO4[2-] */ /* (ref [1]) */
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlS2O8__Al_2SO4(T)  (-5.9)


/* AlSO4[+] = Al[3+] + SO4[2-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlSO4__Al_SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,59.7476,0,-2048.6150,-22.9,0)

/* Al[3+] + SO4[2-] = AlSO4[+] */ /*  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_SO4__AlSO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-59.7476,0,2048.6150,22.9,0)

/* Al[3+] + 2SO4[2-] = Al(SO4)2[-] */ /* (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_2SO4_AlS2O8(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-129.3763,0,5505.6880,47.2,0)

/* Al(SO4)2[-] = Al[3+] + 2SO4[2-] */ /*  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlS2O8__Al_2SO4(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,129.3763,0,-5505.6880,-47.2,0)



/* Potassium-Aluminium compounds */

/* Al[3+] + K[+] + 2H2O = KAlO2 + 4H[+] */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_K_2H2O__KAlO2_4H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,6.4899E02,9.8198E-02,-4.4681E04,-2.3171E02,1.8414E06)

/* KAlO2 + 4H[+] = Al[3+] + K[+] + 2H2O */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KAlO2_4H__Al_K_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-6.4899E02,-9.8198E-02,4.4681E04,2.3171E02,-1.8414E06)



/* Sodium-Aluminium compounds */

/* Al[3+] + Na[+] + 2H2O = NaAlO2 + 4H[+] */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_Na_2H2O__NaAlO2_4H(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,7.0420E02,1.1134E-01,-4.748704,-2.5313E02,2.1869E06)

/* NaAlO2 + 4H[+] = Al[3+] + Na[+] + 2H2O */ /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaAlO2_4H__Al_Na_2H2O(T) \
        TemperatureDependenceOfLog10EquilibriumConstant(T,-7.0420E02,-1.1134E-01,4.748704,2.5313E02,-2.1869E06)



/* References
 * ----------
 *
 * [1] B. Lothenbach, F. Winnefeld, Thermodynamic modelling of the hydration of Portland cement.
 * Cement and Concrete Research, 36:209-226, 2006.
 * 
 * [2] T. Thoenen and D. Kulik. Nagra/psi chemical thermodynamic data base 01/01 for the gem-selektor (v. 2-psi) geochemical modeling code: Release 28-02-03. Technical report, PSI Technical Report TM-44-03-04 about the GEMS version of Nagra/PSI chemical thermodynamic database 01/01, 2003.
 * 
 * [3] D.Jacques, Benchmarking of the cement model and detrimental chemical reactions including temperature dependent parameters, Project near surface disposal of category A waste at Dessel, 2009.
 * 
 * [4] W. Hummel, U. Berner, E. Curti, F.J. Pearson, T. Thoenen, Nagra/PSI Chemical Thermodynamic Data Base 01/01, Universal Publishers/u-PUBLISH.com, USA, also published as Nagra Technical Report NTB 02 - 16, Nagra, Wettingen, Switzerland, 2002.
 */


#endif
