#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Exception.h"
#include "Log10EquilibriumConstantOfHomogeneousReactionInWater.h"



#define LogKr(R) Log10EquilibriumConstantOfHomogeneousReactionInWater(R,T)


void Log10EquilibriumConstantOfHomogeneousReactionInWater_Print(double T)
{
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
  
  
  
  
#define REACTITLE(...) \
  {\
    int c1 = 68 ;\
    int c2 = c1+11 ;\
    int n = printf(__VA_ARGS__) ;\
    while(n < c1) n += printf(" ") ;\
    n += printf("|  Log(k)") ;\
    while(n < c2) n += printf(" ") ;\
    printf("\n") ;\
    n = 0 ;\
    while(n < c1) n += printf("-") ;\
    n += printf("|") ;\
    while(n < c2) n += printf("-") ;\
    printf("\n") ;\
  }

#define PREACT(R,LogK) \
  {\
    double logk = LogK ;\
    int c1 = 68 ;\
    int c2 = c1+11 ;\
    int n = printf(R) ;\
    while(n < c1) n += printf(" ") ;\
    n += printf("| % g",logk) ;\
    while(n < c2) n += printf(" ") ;\
    printf("\n") ;\
  }

  

  printf("\n") ;
  
  REACTITLE("Homogeneous reaction at T = %g",T) ;
  
  printf("Homogeneous reactions involving compounds of type I\n") ;
  printf("---------------------------------------------------\n") ;
  printf("\n") ;
  {
  printf("Water\n") ;
  PREACT("H2O                 = H[+] + OH[-]",logk_h2o) ;
  }
  
  printf("\n") ;
  
  {
  printf("Calcium compounds\n") ;
  PREACT("CaOH[+]             = Ca[2+] + OH[-]",logk_caoh) ;
  PREACT("Ca[2+] + H2O        = CaOH[+] + H[+]",logk_h2o - logk_caoh) ;
  PREACT("Ca(OH)2[0]          = Ca[2+] + 2OH[-]",logk_caoh2) ;
  }
  
  printf("\n") ;
  
  {
  printf("Silicon compounds\n") ;
  PREACT("H3SiO4[-]  + H2O    = H4SiO4[0] + OH[-]",logk_h3sio4) ;
  PREACT("H2SiO4[2-] + H2O    = H3SiO4[-] + OH[-]",logk_h2sio4) ;
  PREACT("H2SiO4[2-] + H[+]   = H3SiO4[-]",logk_h2sio4 - logk_h2o) ;
  PREACT("H2SiO4[2-] + 2H[+]  = H4SiO4[0]",logk_h3sio4 + logk_h2sio4 - (2*logk_h2o)) ;
  }
  
  printf("\n") ;
  
  {
  printf("Sodium compounds\n") ;
  PREACT("NaOH[0]             = Na[+] + OH[-]",logk_naoh) ;
  printf("\n") ;
  printf("Potassium compounds\n") ;
  PREACT("KOH[0]              = K[+] + OH[-]",logk_koh) ;
  }
  
  printf("\n") ;
  
  {
  printf("Carbon compounds\n") ;
  PREACT("H2CO3[0]            = CO2[0] + H2O",logk_h2co3) ;
  PREACT("HCO3[-] + H2O       = H2CO3[0] + OH[-]",logk_hco3) ;
  PREACT("CO3[2-] + H2O       = HCO3[-] + OH[-]",logk_co3) ;
  }
  
  printf("\n") ;
  
  {
  printf("Sulfur compounds\n") ;
  PREACT("H2SO4[0]            = HSO4[-] + H[+]",logk_h2so4) ;
  PREACT("HSO4[-]             = SO4[2-] + H[+]",logk_hso4) ;
  }
  
  printf("\n") ;
  
  {
  printf("Aluminium compounds\n") ;
  PREACT("AlO4H4[-]           = Al[3+] + 4OH[-]",logk_alo4h4) ;
  }
  
  printf("\n") ;
  
  printf("Homogeneous reactions involving compounds of type II\n") ;
  printf("----------------------------------------------------\n") ;
  
  printf("\n") ;
  
  {
  printf("Calcium-silicon compounds\n") ;
  PREACT("CaH3SiO4[+]         = Ca[2+] + H3SiO4[-]",logk_cah3sio4) ;
  PREACT("CaH2SiO4[0]         = Ca[2+] + H2SiO4[2-]",logk_cah2sio4) ;
  }
  
  printf("\n") ;
  
  {
  printf("Calcium-carbon compounds\n") ;
  PREACT("CaHCO3[+]           = Ca[2+] + HCO3[-]",logk_cahco3) ;
  PREACT("CaCO3[0]            = Ca[2+] + CO3[2-]",logk_caco3) ;
  }
  
  printf("\n") ;
  
  {
  printf("Sodium-carbon compounds\n") ;
  PREACT("NaHCO3[0]           = Na[+] + HCO3[-]",logk_nahco3) ;
  PREACT("NaCO3[-]            = Na[+] + CO3[2-]",logk_naco3) ;
  }
  
  printf("\n") ;
  
  {
  printf("Calcium-sulfur compounds\n") ;
  PREACT("CaHSO4[+]           = Ca[2+] + HSO4[-]",logk_cahso4) ;
  PREACT("CaSO4[0]            = Ca[2+] + SO4[2-]",logk_caso4) ;
  }
  
  printf("\n") ;
  
  fflush(stdout) ;
  
  
#undef PREACT
#undef REACTITLE
}
