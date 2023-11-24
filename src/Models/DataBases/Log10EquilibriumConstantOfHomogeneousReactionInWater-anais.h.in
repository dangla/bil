#ifndef LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_H
#define LOG10EQUILIBRIUMCONSTANTOFHOMOGENEOUSREACTIONINWATER_H


#define Log10EquilibriumConstantOfHomogeneousReactionInWater(R,T) (Log10EquilibriumConstantOfHomogeneousReactionInWater_##R(T))



#include "TemperatureDependenceOfLog10EquilibriumConstant.h"

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



/* Autoprotolysis of water */

/* H2O = H[+] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2O__H_OH(T)     TemperatureDependenceOfLog10EquilibriumConstant_293(T,-14.1733,-5.069842E-02,1.332300E+04,1.022445E+02,-1.119669E+06)


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
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_OH__CaOH(T)    TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.39331,5.069842E-02,-1.332300E+04,-1.022445E+02,1.119669E+06)

/* CaOH[+] = Ca[2+] + OH[-]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaOH__Ca_OH(T)    TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.39331,-5.069842E-02,1.332300E+04,1.022445E+02,-1.119669E+06)

/* Ca[2+] + H2O = CaOH[+] + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H2O__CaOH_H(T) TemperatureDependenceOfLog10EquilibriumConstant(T,1.3130E02,2.1418E-02,-1.0190E04,-4.8225E01,2.7033E05)

/* CaOH[+] + H[+] = Ca[2+] + H2O * thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaOH_H__Ca_H2O(T) TemperatureDependenceOfLog10EquilibriumConstant(T,-1.3130E02,-2.1418E-02,1.0190E04,4.8225E01,-2.7033E05)

/* Ca(OH)2[0] = Ca[2+] + 2OH[-] no data*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaO2H2__Ca_2OH(T)      (0)

/* Ca[2+] + 2OH[-] = Ca(OH)2[0] no data*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_2OH__CaO2H2(T)      (0)


/* Silicon compounds */

/* H3SiO4[-] + OH[-] = H2SiO4[2-] + H2O */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H3SiO4_OH__H2SiO4_H2O(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,0.691186,5.069842E-02,-1.5903388E+04,-1.022445E+02,1.119669E+06)

/* H2SiO4[2-] + H2O = H3SiO4[-] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SiO4_H2O__H3SiO4_OH(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,-0.691186,-5.069842E-02,1.5903388E+04,1.022445E+02,-1.119669E+06)

/* H3SiO4[-] + H2O = H4SiO4[0] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H3SiO4_H2O__H4SiO4_OH(T)   (-4.19)

/* H2SiO4[2-] + H2O = H3SiO4[-] + OH[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SiO4_H2O__H3SiO4_OH(T)   (-0.67)

/* H4SiO4 = H3SiO4[-] + H[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__H3SiO4_H(T) TemperatureDependenceOfLog10EquilibriumConstant_293(T,-9.88883,0,-1.337205E+03,0,0)*/

/*  H3SiO4[-] + H[+] = H4SiO4 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__H3SiO4_H(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,9.88883,0,1.337205E+03,0,0)*/

/* H4SiO4 = H3SiO4[-] + H[+] ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__H3SiO4_H(T) TemperatureDependenceOfLog10EquilibriumConstant(T,67.7063,0,-4741.9918,-10.8137,0)

/*  H3SiO4[-] + H[+] = H4SiO4 ref [3]*/ 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H3SiO4_H__H4SiO4(T) TemperatureDependenceOfLog10EquilibriumConstant(T,-67.7063,0,4741.9918,10.8137,0)

/*2H4SiO4 = Si2O2(OH)5[-] + H2O + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_2H4SiO4__Si2O2(OH)5_H_H2O(T) TemperatureDependenceOfLog10EquilibriumConstant(T,-8.50,0,0,0,0)

/*Si2O2(OH)5[-] + H2O + H[+] = 2H4SiO4  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Si2O2(OH)5_H2O_H__2H4SiO4(T) TemperatureDependenceOfLog10EquilibriumConstant(T,8.50,0,0,0,0)

/*H4SiO4 = HSiO3[-] + H2O + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__HSiO3_H_H2O(T) TemperatureDependenceOfLog10EquilibriumConstant(T,-6.0858E02,-1.0072E-01,3.2203E04,2.2007E02,-2.1128E06)

/*HSiO3[-] + H2O + H[+] = H4SiO4 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSiO3_H_H2O__H4SiO4(T) TemperatureDependenceOfLog10EquilibriumConstant(T,6.0858E02,1.0072E-01,-3.2203E04,-2.2007E02,2.1128E06)

/*H4SiO4 = H2SiO4[2-] + 2H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H4SiO4__H2SiO4_H(T) TemperatureDependenceOfLog10EquilibriumConstant(T,3.5337E02,3.9023E-02,-2.2815E04,-1.2975E02,8.2773E05)

/*H2SiO4[2-] + 2H[+] = H4SiO4 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SiO4_H__H4SiO4(T) TemperatureDependenceOfLog10EquilibriumConstant(T,-3.5337E02,-3.9023E-02,2.2815E04,1.2975E02,-8.2773E05)

/* Sodium compounds */

/* NaOH[0] = Na[+] + OH[-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaOH__Na_OH(T)  (0.18) */

/* NaOH[0] + H[+] = Na[+] + H2O  (ref [3])*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaOH_H__Na_H2O(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,-40.0224,0,4902.5816,15.3,0)

/* Na[+] + H2O = NaOH[0] + H[+] (ref [3])*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_H2O__NaOH_H(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,40.0224,0,-4902.5816,-15.3,0)

/* NaOH[0] + H[+] = Na[+] + H2O  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaOH_H__Na_H2O(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,-5.6335E02,-8.5076E-02,3.4108E04,2.0592E02,-1.8192E06)

/* Na[+] + H2O = NaOH[0] + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_H2O__NaOH_H(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,5.6335E02,8.5076E-02,-3.4108E04,-2.0592E02,1.8192E06)



/* Potassium compounds */

/* KOH[0] = K[+] + OH[-]
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KOH__K_OH(T)  (0.46)*/

/*KOH[0] + H[+] = K[+] + H2O ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KOH_H__K_OH_H2O(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,-55.77230945322,0,0,0,0)

/* K[+] + H2O = KOH[0] + H[+] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H_H2O__KOH_H(T)     TemperatureDependenceOfLog10EquilibriumConstant(T,55.77230945322,0,0,0,0)

/*KOH[0] + H[+] = K[+] + H2O thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KOH__K_OH(T)        TemperatureDependenceOfLog10EquilibriumConstant(T,-1.4239E02,-1.6361E-02,1.1313E04,5.1782E01,-3.8639E05)

/* K[+] + H2O = KOH[0] + H[+]thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_K_H2O__KOH_H(T)     TemperatureDependenceOfLog10EquilibriumConstant(T,1.4239E02,1.6361E-02,-1.1313E04,-5.1782E01,3.8639E05)



/* Carbon compounds:
 * Carbone dioxide (dissolved): CO2(aq) = CO2[0]
 * Carbonic acid: H2CO3[0]
 * Bicarbonate ion (or hydrogen carbonate ion): HCO3[-]
 * Carbonate ion: CO3[2-] 
 * CO2[0] + H2CO3[0] noted as H2CO3[*] */

/* CO2[0] + H2O = H2CO3[0] */   /* Hydration reaction of CO2 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO2_H2O__H2CO3(T)    (-2.77)

/* H2CO3[0] = CO2[0] + H2O */   /* Hydration reaction of CO2 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__CO2_H2O(T)    (2.77)

/* H2CO3[0] = HCO3[-] + H[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__HCO3_H(T)     (-3.6)

/* HCO3[-] + H2O = H2CO3[0] + OH[-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H2O__H2CO3_OH(T)   (-10.4)

/* H2CO3[*] = HCO3[-] + H[+] 
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2CO3__HCO3_H(T)     (-6.3)

/* HCO3[-] + H2O = H2CO3[*] + OH[-] 
//#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H2O__H2CO3_OH(T)   (-7.648)

/* HCO3[-] = CO2[0] + OH[-]
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3__CO2_OH(T)   (-7.61713)

/* HCO3[-] + H[+] = CO2[0] + H2O 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3_H__CO2_H2O(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,6.38287,6.091964E-02,-2.183437E+04,-1.268339E+02,1.684915E+06)

/* CO2[0] + H2O = HCO3[-] + H[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO2_H2O__HCO3_H(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,-6.38287,-6.091964E-02,2.183437E+04,1.268339E+02,-1.684915E+06)

/* HCO3[-] = CO3[2-] + H[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HCO3__H_CO3(T)   TemperatureDependenceOfLog10EquilibriumConstant_293(T,-10.3771,-3.252849E-02,5.151790E+03,3.892561E+01,-5.637139E+05)

/* CO3[2-] + H2O = HCO3[-] + OH[-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CO3_H2O__HCO3_OH(T)   (-3.67)
*/

/* Sulfur compounds */ /* from EQ3/6 */

/* H2S[0] = HS[-] + H[+]
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S__HS_H(T)    (-LOGK(T,7.4159,6.9877,6.6467,6.4827,6.496,6.6831,7.0225,7.5536))*/

/*H2S[0] = HS[-] + H[+]* ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S__HS_H(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,92.54543936710,0,-5444.49694476074,-32.8459524374,0))
	   
/* HS[-] + H[+] = H2S[0] * ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HS_H__H2S(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,-92.54543936710,0,5444.49694476074,32.8459524374,0))

/*H2S[0] = HS[-] + H[+]* thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S__HS_H(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,-7.4841E02,-1.1982E-01,4.1347E04,2.7073E02,-2.7055E06))
	   
/* HS[-] + H[+] = H2S[0] * thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HS_H__H2S(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,7.4841E02,1.1982E-01,-4.1347E04,-2.7073E02,2.7055E06))

/* HS[-] = S[2-] + H[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HS__S_H(T)      (-LOGK(T,13.7100,12.9351,12.0082,11.1018,10.1202,9.2545,8.4250,7.5568))*/

/* HS[-] = S[2-] + H[+] ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HS__S_H(T)      (TemperatureDependenceOfLog10EquilibriumConstant(T,148.4964,0,-20433.9698,-57.0,0))

/*  S[2-] + H[+] = HS[-] ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_S_H__HS(T)      (TemperatureDependenceOfLog10EquilibriumConstant(T,-148.4964,0,20433.9698,57.0,0))
 
 /* H2SO4[0] = HSO4[-] + H[+] no data*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SO4__HSO4_H(T)    (6)

/* HSO4[-] + H[+] = H2SO4[0] no data*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO4_H__H2SO4(T)    (6)

/* HSO4[-] = SO4[2-] + H[+]
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO4__SO4_H(T)  (-LOGK(T,1.7193,1.9791,2.4371,3.0002,3.7234,4.4683,5.2633,6.1799))*/

/* HSO4[-] = SO4[2-] + H[+] ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO4__SO4_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,95.4496,0,-3428.4424,-34.7,0))

/*  SO4[2-] + H[+] = HSO4[-] ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO4_H__HSO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-95.4496,0,3428.4424,34.7,0))

/* HSO4[-] = SO4[2-] + H[+] thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO4__SO4_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-8.1698E02,-1.2950E-01,7.3030E04,4.6772E02,-4.5780E06))

/*  SO4[2-] + H[+] = HSO4[-] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO4_H__HSO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,8.1698E02,1.2950E-01,-7.3030E04,-4.6772E02,4.5780E06))


/* HSO3[-] = SO3[2-] + H[+] */   /* Ref chess.tdb 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO3__SO3_H(T) (-7.22)  (-LOGK(T,7.1301,7.2054,7.4642,7.8631,8.4488,9.11,9.8607,10.7638))*/

/* HSO3[-] = SO3[2-] + H[+] */   /* Ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO3__SO3_H(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,87.5922,0,-3812.5522,-33.1,0))

/* SO3[2-] + H[+] = HSO3[-] */   /* Ref [3] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO3_H__HSO3(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-87.5922,0,3812.5522,33.1,0))

/* HSO3[-] = SO3[2-] + H[+] */   /* thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HSO3__SO3_H(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,-8.1038E02,-1.3068E-01,4.5360E04,2.9174E02,-2.8320E06))

/* SO3[2-] + H[+] = HSO3[-] */   /*  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO3_H__HSO3(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,8.1038E02,1.3068E-01,-4.5360E04,-2.9174E02,2.8320E06))

/*SO3[2-] + 2H[+] = H2SO3 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO3_H__H2SO3(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,1.2948E03,2.1816E-01,-7.3030E04,-4.6772E02,4.5780E06))

/*H2SO3 = SO3[2-] + 2H[+]  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2SO3__SO3_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-1.2948E03,-2.1816E-01,7.3030E04,4.6772E02,-4.5780E06))

/*SO3[2-] + SO4[2-] = 0.5O2 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_SO3_SO4__O2(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,9.6719E01,1.4161E-02,-2.0795E04,-3.3793E01,5.1632E05))

/*0.5O2 = SO3[2-] + SO4[2-]  thermodem*/
#define Log10Equilibriu_mConstantOfHomogeneousReactionInWater_O2__SO3_SO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-9.6719E01,-1.4161E-02,2.0795E04,3.3793E01,-5.1632E05))

/*S2O3[2-] + 2H[+] = H2S2O3 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_S2O3_H__H2S2O3(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,1.4979E03,2.3814E-01,-8.4049E04,-5.4207E02,5.0380E06))

/*H2S2O3 = S2O3[2-] + 2H[+]  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S2O3__S2O3_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-1.4979E03,-2.3814E-01,8.4049E04,5.4207E02,-5.0380E06))

/*S2O4[2-] + 2H[+] = H2S2O4 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_S2O4_H__H2S2O4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,1.5238E03,2.4188E-01,-8.5504E04,-5.5134E02,5.1466E06))

/*H2S2O4 = S2O4[2-] + 2H[+]  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_H2S2O4__S2O4_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-1.5238E03,-2.4188E-01,8.5504E04,5.5134E02,-5.1466E06))

/* Aluminium compounds */

/* Al(OH)[2+] = Al[3+] + OH[-] (ref [1]) 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlOH__Al_OH(T)    (-14 + 4.957)

/* Al(OH)2[+] = Al[3+] + 2OH[-] (ref [1]) 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2H2__Al_2OH(T)  (-28 + 10.594)

/* Al(OH)3[0] = Al[3+] + 3OH[-] (ref [1]) 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO3H3__Al_3OH(T)  (-42 + 16.432)

/* Al(OH)4[-] = Al[3+] + 4OH[-] (ref [1]) 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4__Al_4OH(T)  (-56 + 22.879)*/


/*Al[3+] + 2H2O = Al(OH)2[+] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlO2H2(T)  (TemperatureDependenceOfLog10EquilibriumConstant(32.5363,0,-6492.5404,-8.6,0))

/*Al(OH)2[+] = Al[3+] + 2H2O  ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2H2__Al_2H2O(T) (TemperatureDependenceOfLog10EquilibriumConstant(-32.5363,0,6492.5404,8.6,0))
 
/*Al[3+] + 2H2O = Al(OH)2[+] + 2H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlO2H2_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(3.0243E02,5.4152E-02,-1.1506E0+4,-5.8055E01,6.0767E05))

/*Al(OH)2[+] + 2H[+] = Al[3+] + 2H2O thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2H2_H__Al_H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(-3.0243E02,-5.4152E-02,1.1506E0+4,5.8055E01,-6.0767E05))

/*Al[3+] + 4H2O = Al(OH)4[-] + 4H[+] ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_4H2O__AlO4H4_4H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(32.9602,0,-10341.6037,-8.5,0)

/* Al(OH)4[-] + 4H[+] = Al[3+] + 4H2O ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4_4H__Al_4H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(-32.9602,0,10341.6037,8.5,0)

/*Al[3+] + 3H2O = Al(OH)3 + 3H[+] ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlO3H3_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(91.6640,0,-11712.9732,-27.8,0))

/* Al(OH)3 + 3H[+] = Al[3+] + 3H2O  ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO3H3__Al_3H2O_3H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(-91.6640,0,11712.9732,27.8,0))

/*Al[3+] + H2O = AlOH[2+] + H[+] ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlOH2_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(-34.5703,0,-869.7066,13.1,0))

/* AlOH[2+] + H[+] = Al[3+] + H2O ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlOH2_H__Al_H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(34.5703,0,869.7066,-13.1,0))

/*Al[3+] + H2O = Al(OH)[2+] + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlOH_H(T) (TemperatureDependenceOfLog10EquilibriumConstant(1.6146E02,3.0185E-02,-1.1506E04,-5.8055E01,6.0767E05))

/*Al(OH)[2+] + H[+] = Al[3+] + H2O  thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlOH2_H__Al_H2O(T) (TemperatureDependenceOfLog10EquilibriumConstant(-1.6146E02,-3.0185E-02,1.1506E04,5.8055E01,-6.0767E05))

/*Al[3+] + 2H2O = AlO2[-] + 4H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__AlO2_H(T) (TemperatureDependenceOfLog10EquilibriumConstant(-1.7805E02,-2.6890E02,1.8672E03,6.6833E01,-7.5044E05))

/*AlO2[-] + 4H[+] = Al[3+] + 2H2O thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO2_H__Al_H2O(T) (TemperatureDependenceOfLog10EquilibriumConstant(1.7805E02,2.6890E02,-1.8672E03,-6.6833E01,7.5044E05))

/*Al[3+] + 2H2O = HAlO2 + 3H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H2O__HAlO2_H(T) (TemperatureDependenceOfLog10EquilibriumConstant(3.4561E02,6.0311E-02,-2.5787E04,-1.2376E02,1.1292E04))

/*HAlO2 + 3H[+] = Al[3+] + 2H2O thermodem*/ 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_HAlO2_H__Al_H2O(T) (TemperatureDependenceOfLog10EquilibriumConstant(-3.4561E02,-6.0311E-02,2.5787E04,1.2376E02,-1.1292E04))

/* 2. Chemical reactions involving compounds of type II. */
 
 
/* Calcium-Silicon compounds */
     
/* Ca[2+] + H2SiO4[2-] = CaH2SiO4 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H2SiO4__CaH2SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,4.6,0,0,0,0)
     
/* CaH2SiO4 = Ca[2+] + H2SiO4[2-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH2SiO4__Ca_H2SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,-4.6,0,0,0,0)
       
/* Ca[2+] + H3SiO4[-]  = CaH3SiO4[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H3SiO4__CaH3SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.2,0,0,0,0)
     
/* CaH3SiO4[+] = Ca[2+] + H3SiO4[-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH3SiO4__Ca_H3SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.2,0,0,0,0)*/

/* Ca[2+] + H2SiO4[2-] = CaH2SiO4  ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H2SiO4__CaH2SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,1371.6455,0,0)

/*  CaH2SiO4 = Ca[2+] + H2SiO4[2-]   ref[3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH2SiO4__Ca_H2SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,-1371.6455,0,0)
 
/* CaH3SiO4[+] = Ca[2+] + H3SiO4[-] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaH3SiO4__Ca_H3SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,88.1066,0,-3919.1159,-30.8,0)

/* Ca[2+] + H3SiO4[-] = CaH3SiO4[+]  ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_H3SiO4__CaH3SiO4(T)  TemperatureDependenceOfLog10EquilibriumConstant(T,-88.1066,0,3919.1159,30.8,0)


/* Calcium-Carbon compounds 

 Ca[2+] + HCO3[-]    = CaHCO3[+] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HCO3__CaHCO3(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,1.06756,3.129400E-01,-3.476505E+04,-4.787820E+02,0) 

CaHCO3[+] = Ca[2+] + HCO3[-]
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaHCO3__Ca_HCO3(T)  TemperatureDependenceOfLog10EquilibriumConstant_293(T,-1.06756,-3.129400E-01,3.476505E+04,4.787820E+02,0) 

 Ca[2+] + CO3[2-]    = CaCO3[0] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_CO3__CaCO3(T)    TemperatureDependenceOfLog10EquilibriumConstant_293(T,3.1836,-2.9944401E-01,3.551275E+04,4.8581799E+02,0) 

 CaCO3[0] = Ca[2+] + CO3[2-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaCO3__Ca_CO3(T)    TemperatureDependenceOfLog10EquilibriumConstant_293(T,-3.1836,2.9944401E-01,-3.551275E+04,-4.8581799E+02,0) 
*/


/*Sodium-Carbon compounds 

NaHCO3[0] = Na[+] + HCO3[-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaHCO3__Na_HCO3(T)  (0.25) 

NaCO3[-] = Na[+] + CO3[2-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaCO3__Na_CO3(T)  (-1.27) 

*/
 
/* Sodium-Sulfur compounds */
/* NaSO4[-] = Na[+] + SO4[2-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaSO4__Na_SO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,67.9470,0,-2908.6947,-23.8,0))

/*  Na[+] + SO4[2-] =  NaSO4[-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_SO4__NaSO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,-67.9470,0,2908.6947,23.8,0))

/*  Na[+] + SO4[2-] =  NaSO4[-] thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Na_SO4__NaSO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,9.3588E02,1.4439E-01,-5.3023E04,-3.3840E02,3.3064E06))

/*   NaSO4[-] = Na[+] + SO4[2-] thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaSO4__Na_SO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,-9.3588E02,-1.4439E-01,5.3023E04,3.3840E02,-3.3064E06))
 
 
/* Potassium-Sulfur compounds */
/* KSO4[-] = K[+] + SO4[2-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KSO4__K_SO4(T)     (TemperatureDependenceOfLog10EquilibriumConstant(T,72.9765,0,-3150.0583,-25.6,0))

/*  K[+] + SO4[2-] = KSO4[-](ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KSO4__K_SO4(T)     (TemperatureDependenceOfLog10EquilibriumConstant(T,-72.9765,0,3150.0583,25.6,0))

/* KSO4[-] = K[+] + SO4[2-] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KSO4__K_SO4(T)     (TemperatureDependenceOfLog10EquilibriumConstant(T,-9.1525E02,-1.4349E-01,5.1254E04,3.3152E02,-3.1178E06))

/*  K[+] + SO4[2-] = KSO4[-]thermodem */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_K_SO4__KSO4(T)     (TemperatureDependenceOfLog10EquilibriumConstant(T,9.1525E02,1.4349E-01,-5.1254E04,-3.3152E02,3.1178E06))


/* Calcium-Sulfur compounds */

/* Ca[2+] + HS[-]      = CaHS[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HS__CaHS(T)  (1.10585)

/* Ca[2+] + S[2-]      = CaS[0] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_S__CaS(T)  (3.544068)

/* Ca[2+] + HSO4[-]    = CaHSO4[+] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_HSO4__CaHSO4(T)  (1.10585)

/* CaHSO4[+] = Ca[2+] + HSO4[-] */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaHSO4__Ca_HSO4(T)  (-1.10585)

/* Ca[2+] + SO4[2-]    = CaSO4[0] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4__CaSO4(T)  (3.146128)*/

/*Ca[2+] + SO4[2-] = CaSO4[0] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4__CaSO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,-64.2534,0,2770.0342,23.1418,0))

/*CaSO4[0] = Ca[2+] + SO4[2-] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,64.2534,0,-2770.0342,-23.1418,0))

/*Ca[2+] + SO4[2-] = CaSO4[0] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Ca_SO4__CaSO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,1.72034E03,2.6573E-01,-9.4255E04,-6.2356E02,5.4973E06))

/*CaSO4[0] = Ca[2+] + SO4[2-] ref [3]*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,-1.72034E03,-2.6573E-01,9.4255E04,6.2356E02,-5.4973E06))



/* CaSO4[0] = Ca[2+] + SO4[2-] 
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_CaSO4__Ca_SO4(T)  (-LOGK(T,2.0713,2.1111,2.2647,2.5111,2.9101,3.4328,4.1424,5.1853))*/

/* Aluminium-Silicon compounds */

/* AlH3SiO4[2+] = Al[3+] + H3SiO4[-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlH3SiO4__Al_H3SiO4(T)          (TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,-2206.3355,0,0))

/* Al[3+] + H3SiO4[-] = AlH3SiO4[2+] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H3SiO4__AlH3SiO4(T)          (TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,2206.3355,0,0))

/* Al(OH)4H4SiO4[-] + H2O = Al(OH)4[-] + H4SiO4[0] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4H4SiO4_H2O__AlO4H4_H4SiO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,-1073.3948,0,0))

/*  Al(OH)4[-] + H4SiO4[0] = Al(OH)4H4SiO4[-] + H2O (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlO4H4_H4SiO4__AlO4H4H4SiO4_H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,0,0,1073.3948,0,0))

/*Al[3+] + H4SiO4 = AlH3SiO4[2+] + H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_H4SiO4__AlH3SiO4_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-3.4331E02,-2.3469E-02,2.5866E04,1.1700E02,-2.5178E06))

/*AlH3SiO4[2+] + H[+] = Al[3+] + H4SiO4 thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlH3SiO4_H__Al_H4SiO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,3.4331E02,2.3469E-02,-2.5866E04,-1.1700E02,2.5178E06))


/* Aluminium-Sulfur compounds */

/* AlSO4[+] = Al[3+] + SO4[2-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlSO4__Al_SO4(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,59.7476,0,-2048.6150,-22.9,0))

/* Al[3+] + SO4[2-] = AlSO4[+]  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_SO4__AlSO4(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,-59.7476,0,2048.6150,22.9,0))

/*Al[3+] + SO4[2-] = AlSO4[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_SO4__AlSO4(T)   (TemperatureDependenceOfLog10EquilibriumConstant(T,2.3193E03,3.6143E-01,-1.3494E05,-8.3586E02,8.6189E06))

/* AlSO4[+] = Al[3+] + SO4[2-] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlSO4__Al_SO4(T)    (TemperatureDependenceOfLog10EquilibriumConstant(T,-2.3193E03,-3.6143E-01,1.3494E05,8.3586E02,-8.6189E06))

/* Al[3+] + 2SO4[2-] = Al(SO4)2[-] (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_2SO4_AlS2O8(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-129.3763,0,5505.6880,47.2,0))

/* Al(SO4)2[-] = Al[3+] + 2SO4[2-]  (ref [3]) */
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_AlS2O8__Al_2SO4(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,129.3763,0,-5505.6880,-47.2,0))

/*Al[3+] + K[+] + 2H2O = KAlO2 + 4H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_K_H2O__KAlO2_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,6.4899E02,9.8198E-02,-4.4681E04,-2.3171E02,1.8414E06))

/* KAlO2 + 4H[+] = Al[3+] + K[+] + 2H2O thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_KAlO2_H__Al_K_H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-6.4899E02,-9.8198E-02,4.4681E04,2.3171E02,-1.8414E06))

/*Al[3+] + Na[+] + 2H2O = NaAlO2 + 4H[+] thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_Al_Na_H2O__NaAlO2_H(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,7.0420E02,1.1134E-01,-4.748704,-2.5313E02,2.1869E06))

/* NaAlO2 + 4H[+] = Al[3+] + Na[+] + 2H2O thermodem*/
#define Log10EquilibriumConstantOfHomogeneousReactionInWater_NaAlO2_H__Al_Na_H2O(T)  (TemperatureDependenceOfLog10EquilibriumConstant(T,-7.0420E02,-1.1134E-01,4.748704,2.5313E02,-2.1869E06))



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
