#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Curves.h"
#include "BilPath.h"
#include "Log10ActivityOfWaterInBrine.h"


static Log10ActivityOfWaterInBrine_t* instancelogaw = NULL ;

static Log10ActivityOfWaterInBrine_t* Log10ActivityOfWaterInBrine_Create(void) ;
static int pm(const char*) ;



Log10ActivityOfWaterInBrine_t* Log10ActivityOfWaterInBrine_Create(void)
{
  Log10ActivityOfWaterInBrine_t* logaw = (Log10ActivityOfWaterInBrine_t*) malloc(sizeof(Log10ActivityOfWaterInBrine_t)) ;
  
  if(!logaw) {
    arret("Log10ActivityOfWaterInBrine_Create(1)") ;
  }
  
  {
    /* The "XXXXXX" is here to save space for the salt name (see below) */
    char line[] = "Curves = " BIL_PATH "/src/" __FILE__ "XXXXXX" ;
    
    /* Change the suffix ".c" in ".\0" */
    //strncpy(strrchr(line,'c'),"\0",1) ;
    
    //Message_Direct("This current file is %s\n",__FILE__) ;
    
    {
      int nbofcurves = 2 ;
      Curves_t* curves = Curves_Create(nbofcurves) ;
      int nbofsalts = 0 ;
      
      /* NaCl */
      {
        char salt[] = ".NaCl" ;
        strcpy(strrchr(line,'.'),salt) ;
        
        nbofsalts += Curves_ReadCurves(curves,line) ;
      }
      
      /* Na2SO4 */
      {
        char salt[] = ".Na2SO4" ;
        strcpy(strrchr(line,'.'),salt) ;
        
        nbofsalts += Curves_ReadCurves(curves,line) ;
      }
      
      if(nbofsalts > nbofcurves) {
        arret("Log10ActivityOfWaterInBrine_Create(2)") ;
      }
      
      Log10ActivityOfWaterInBrine_GetCurves(logaw) = curves ;
      Log10ActivityOfWaterInBrine_GetNbOfSalt(logaw) = nbofsalts ;
    }
  }
  
  return(logaw) ;
}


double Log10ActivityOfWaterInBrine(const char* salt,double rho)
{
  if(!instancelogaw) {
    instancelogaw = Log10ActivityOfWaterInBrine_Create() ;
  }
  
  {
    Curves_t* curves = Log10ActivityOfWaterInBrine_GetCurves(instancelogaw) ;
    Curve_t*  curve  = Curves_GetCurve(curves) ;
    int nbofsalts = Log10ActivityOfWaterInBrine_GetNbOfSalt(instancelogaw) ;
    int salt_index = pm(salt) ;
    
    if(salt_index < 0 || salt_index > nbofsalts) {
      arret("Log10ActivityOfWaterInBrine(1): index out of range") ;
    }

    double log10aw = Curve_ComputeValue(curve + salt_index,rho) ;
  
    return(log10aw) ;
  }
}


int pm(const char* s)
{
       if(!strcmp(s,"NaCl"))     return (0) ;
  else if(!strcmp(s,"Na2SO4"))   return (1) ;
  else return (-1) ;
}
