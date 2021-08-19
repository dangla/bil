#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Message.h"
#include "Mry.h"
#include "Curves.h"
#include "BilPath.h"
#include "Temperature.h"
#include "Utils.h"
#include "LogActivityOfWaterInBrine.h"


/* Units
 * ----- */
#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define Meter   (InternationalSystemOfUnits_OneMeter)
#define M3      (Meter*Meter*Meter)
#define dm      (0.1*Meter)
#define cm      (0.01*Meter)
#define dm2     (dm*dm)
#define dm3     (dm*dm*dm)
#define Liter   dm3
#define cm3     (cm*cm*cm)
#define Pascal  (InternationalSystemOfUnits_OnePascal)
#define MPa     (1.e6*Pascal)
#define GPa     (1.e9*Pascal)
#define Mol     InternationalSystemOfUnits_OneMole
#define sec     InternationalSystemOfUnits_OneSecond


static LogActivityOfWaterInBrine_t* instancelogaw = NULL ;

static int pm(const char*) ;



LogActivityOfWaterInBrine_t* LogActivityOfWaterInBrine_Create(void)
{
  LogActivityOfWaterInBrine_t* logaw = (LogActivityOfWaterInBrine_t*) Mry_New(LogActivityOfWaterInBrine_t) ;
  
  
  /* Allocation of space for the temperature */
  #if 0
  {
    Temperature_t* temp = Temperature_Create() ;
    
    LogActivityOfWaterInBrine_GetTemperature(logaw) = temp ;
  }
  #endif
  
  
  /* The curves */
  {
    /* The beginning of the command */
    #define BASEFILE  Utils_STR(__BASEFILE__)
    char cmdbeg[] = "Curves = " BIL_PATH "/src/Models/DataBases/" BASEFILE "." ;
    
    //Message_Direct("This current file is %s\n",BASEFILE) ;
    #undef BASEFILE
    
    {
      int nbofsalts = 3 ;
      Curves_t* curves = Curves_Create(2*nbofsalts) ;
      int nbofcurves = 0 ;
      
      /* Build 2 curves at T = 253 K and T = 353 K */
      #define CMD(SALT,CMAX) \
              " concentration = Range{c1 = 0 , c2 = "#CMAX" , n = 1000} logaw253 = LinLee(1){ T = 253 , salt = "#SALT" } logaw353 = LinLee(1){ T = 353 , salt = "#SALT" }"
      
      /* NaCl */
      {
        char salt[] = "NaCl" ;
        char cmd[] = CMD(NaCl,6.e3) ;
        int len = strlen(cmdbeg) + strlen(salt) + strlen(cmd) ;
        char line[len] ;
        
        strcpy(line,cmdbeg) ;
        strcat(line,salt) ;
        strcat(line,cmd) ;
        
        nbofcurves += Curves_ReadCurves(curves,line) ;
      }
      
      /* CaCl2 */
      {
        char salt[] = "CaCl2" ;
        char cmd[] = CMD(CaCl2,6.e3) ;
        int len = strlen(cmdbeg) + strlen(salt) + strlen(cmd) ;
        char line[len] ;
        
        strcpy(line,cmdbeg) ;
        strcat(line,salt) ;
        strcat(line,cmd) ;
        
        nbofcurves += Curves_ReadCurves(curves,line) ;
      }
      
      /* Na2SO4 */
      {
        char salt[] = "Na2SO4" ;
        char cmd[] = CMD(Na2SO4,6.e3) ;
        int len = strlen(cmdbeg) + strlen(salt) + strlen(cmd) ;
        char line[len] ;
        
        strcpy(line,cmdbeg) ;
        strcat(line,salt) ;
        strcat(line,cmd) ;
        
        nbofcurves += Curves_ReadCurves(curves,line) ;
      }
      #undef CMD
      
      if(nbofcurves > 2*nbofsalts) {
        arret("LogActivityOfWaterInBrine_Create(2)") ;
      }
      
      LogActivityOfWaterInBrine_GetCurves(logaw) = curves ;
      LogActivityOfWaterInBrine_GetNbOfSalt(logaw) = nbofsalts ;
    }
  }
  
  return(logaw) ;
}


double LogActivityOfWaterInBrine(double c,double T,const char* salt)
{
  if(!instancelogaw) {
    instancelogaw = LogActivityOfWaterInBrine_Create() ;
  }
  
  {
    int nbofsalts = LogActivityOfWaterInBrine_GetNbOfSalt(instancelogaw) ;
    int salt_index = pm(salt) ;
    
    if(salt_index < 0 || salt_index >= nbofsalts) {
      arret("LogActivityOfWaterInBrine(1): index out of range") ;
    }

    {
      Curve_t* curve = LogActivityOfWaterInBrine_GetCurve(instancelogaw) ;
      double cs = c*Mol/M3 ;
      double logaw253 = Curve_ComputeValue(curve + 2*salt_index,cs) ;
      double logaw353 = Curve_ComputeValue(curve + 2*salt_index + 1,cs) ;
      double logaw = 0.01*(logaw353*(T-253) + logaw253*(353-T)) ;
  
      return(logaw) ;
    }
  }
}


int pm(const char* s)
{
       if(!strcmp(s,"NaCl"))     return (0) ;
  else if(!strcmp(s,"CaCl"))     return (1) ;
  else if(!strcmp(s,"Na2SO4"))   return (2) ;
  else return (-1) ;
}


