#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Exception.h"
#include "Log10DissociationConstantOfCalciumCarbonate.h"



void Log10DissociationConstantOfCalciumCarbonate_Print(double T)
{

#define REACTITLE(...) \
  {\
    int c1 = 52 ;\
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
    int c1 = 52 ;\
    int c2 = c1+11 ;\
    int n = printf(R) ;\
    while(n < c1) n += printf(" ") ;\
    n += printf("| % g",logk) ;\
    while(n < c2) n += printf(" ") ;\
    printf("\n") ;\
  }
  
  printf("\n") ;
  
  
  /* Calcite dissociation reaction */
  {
    #define LogKd(R) \
            Log10DissociationConstantOfCalciumCarbonate(R,T)
            
    double logk_calcite      = LogKd(Calcite__Ca_CO3) ;
    double logk_aragonite    = LogKd(Aragonite__Ca_CO3) ;
    double logk_vaterite     = LogKd(Vaterite__Ca_CO3) ;
    
    REACTITLE("Calcium carbonate dissociation reaction at T = %g",T) ;
    PREACT("Calcite(s)         =  Ca[2+] + CO3[2-]",logk_calcite) ;
    PREACT("Aragonite(s)       =  Ca[2+] + CO3[2-]",logk_aragonite) ;
    PREACT("Vaterite(s)        =  Ca[2+] + CO3[2-]",logk_vaterite) ;
    
    #undef LogKd
  }
  
  printf("\n") ;
  
  fflush(stdout) ;
  
  
#undef PREACT
#undef REACTITLE
}
