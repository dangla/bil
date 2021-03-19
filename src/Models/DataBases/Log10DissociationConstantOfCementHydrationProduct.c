#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Exception.h"
#include "Log10EquilibriumConstantOfHomogeneousReactionInWater.h"
#include "Log10DissociationConstantOfCementHydrationProduct.h"



void Log10DissociationConstantOfCementHydrationProduct_Print(double T)
{
  /* Some other constants */
  
  #define LogKr(R) Log10EquilibriumConstantOfHomogeneousReactionInWater(R,T)
  double logk_h2o      = LogKr(H2O__H_OH) ;
  #undef LogKr
  
#define REACTITLE(...) \
  {\
    int c1 = 72 ;\
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
    int c1 = 72 ;\
    int c2 = c1+11 ;\
    int n = printf(R) ;\
    while(n < c1) n += printf(" ") ;\
    n += printf("| % g",logk) ;\
    while(n < c2) n += printf(" ") ;\
    printf("\n") ;\
  }
  
  printf("\n") ;
  
  
  /* Cement hydrates dissociation reactions */
  {
    #define LogKd(R) \
            Log10DissociationConstantOfCementHydrationProduct(R,T)
        
    double logk_ch       = LogKd(CH__Ca_2OH) ;
    double logk_sh       = LogKd(S_2H2O__H4SiO4) ;
    double logk_csh2     = LogKd(CSH2__Ca_SO4_2H2O) ;
    double logk_ah3      = LogKd(AH3__2Al_6OH) ;
    double logk_ah3_bis  = LogKd(AH3_2OH__2AlO4H4) ;
    double logk_afm      = LogKd(AFm_12H__4Ca_2Al_SO4_18H2O) ;
    double logk_aft      = LogKd(AFt_12H__6Ca_2Al_3SO4_38H2O) ;
    double logk_aft_bis  = LogKd(AFt__6Ca_2AlO4H4_3SO4_4OH_26H2O) ;
    double logk_c3ah6    = LogKd(C3AH6_12H__3Ca_2Al_12H2O) ;
    double logk_c3ah6_bis= LogKd(C3AH6__3Ca_2AlO4H4_4OH) ;
    double logk_c2ah8    = LogKd(C2AH8__2Ca_2AlO4H4_2OH_3H2O) ;
    double logk_cah10    = LogKd(CAH10__Ca_2AlO4H4_6H2O) ;
    double logk_Friedel  = LogKd(FriedelSalt__4Ca_2AlO4H4_2Cl_4OH_4H2O) ;
    double logk_Kuzel    = LogKd(KuzelSalt__4Ca_2AlO4H4_Cl_halfSO4_4OH_6H2O) ;
  
    REACTITLE("Cement hydrates dissociation reactions at T = %g",T) ;
    PREACT("CH(s)              =  Ca[2+] + 2OH[-]",logk_ch) ;
    PREACT("SH(s) + H2O        =  H4SiO4[0]",logk_sh) ;
    
    PREACT("CSH2(s)            =  Ca[2+] + SO4[2-] + 2H2O",logk_csh2) ;
    
    PREACT("AH3(s)             = 2Al[3+] + 6OH[-]",logk_ah3) ;
    PREACT("AH3(s) + 2OH[-]    = 2Al(OH)4[-]",logk_ah3_bis) ;
    
    PREACT("AFm(s) + 12H[+]    = 4Ca[2+] + 2Al[3+] + SO4[2-] + 18H2O",logk_afm) ;
    PREACT("AFm(s)             = 4Ca[2+] + 2Al[3+] + SO4[2-] + 12OH[-] + 6H2O",logk_afm + 12*logk_h2o) ;
    
    PREACT("AFt(s) + 12H[+]    = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 38H2O",logk_aft) ;
    PREACT("AFt(s)             = 6Ca[2+] + 2Al[3+] + 3SO4[2-] + 12OH[-] + 26H2O",logk_aft + 12*logk_h2o) ;
    PREACT("AFt(s)             = 6Ca[2+] + 2Al(OH)4[-] + 3SO4[2-] + 4OH[-] + 26H2O",logk_aft_bis) ;
        
    PREACT("C3AH6(s) + 12H[+]  = 3Ca[2+] + 2Al[3+] + 12H2O",logk_c3ah6) ;
    PREACT("C3AH6(s)           = 3Ca[2+] + 2Al[3+] + 12OH[-]",logk_c3ah6 + 12*logk_h2o) ;
    PREACT("C3AH6(s)           = 3Ca[2+] + 2Al(OH)4[-] + 4OH[-]",logk_c3ah6_bis) ;

    PREACT("C2AH8(s)           = 2Ca[2+] + 2Al(OH)4[-] + 2OH[-] + 3H2O",logk_c2ah8) ;
    PREACT("CAH10(s)           =  Ca[2+] + 2Al(OH)4[-] + 6H2O",logk_cah10) ;
    
    PREACT("Friedel's salt(s)  = 4Ca[2+] + 2Al(OH)4[-] + 2Cl[-] + 4OH[-] + 4H2O",logk_Friedel) ;
    PREACT("Kuzel's salt(s)    = Friedel's salt(l) - Cl[-] + 0.5SO4[2-] + 2H2O", logk_Kuzel) ; //4Ca[2+] + 2Al(OH)4[-] + Cl[-] + 0.5SO4[2-] + 4OH[-] + 6H2O",logk_Kuzel) ;
    
    #undef LogKd
  }
  
  printf("\n") ;
  
  fflush(stdout) ;
  
  
#undef PREACT
#undef REACTITLE
}
