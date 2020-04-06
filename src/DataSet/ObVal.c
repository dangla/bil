#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "String.h"
#include "Mry.h"
#include "ObVal.h"




ObVal_t*  ObVal_New(void)
{
  ObVal_t* obval = (ObVal_t*) Mry_New(ObVal_t) ;
  
  
  /* Allocation of space for the name of unknown */
  {
    char* name = (char*) Mry_New(char[ObVal_MaxLengthOfKeyWord]) ;
  
    ObVal_GetNameOfUnknown(obval) = name ;
  }
  
  /* Initialization */
  ObVal_GetValue(obval) = -1 ; /* arbitrary negative */
  ObVal_SetTypeToAbsolute(obval) ;
  ObVal_GetRelaxationFactor(obval) = 1 ;
  
  return(obval) ;
}




void  ObVal_Scan(ObVal_t* obval,DataFile_t* datafile)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;

    
  /* Unknown and value */
  {
    char name[ObVal_MaxLengthOfKeyWord] ;
    double v ;
    int n = String_Scan(line," %s = %le",name,&v) ;
        
    if(n) {
      strcpy(ObVal_GetNameOfUnknown(obval),name) ;
      ObVal_GetValue(obval) = v ;
    } else {
      arret("ObVal_Scan: no unknown") ;
    }
      
    if(strlen(name) > ObVal_MaxLengthOfKeyWord-1)  {
      arret("ObVal_Scan: too long name of unknown") ;
    }
        
    if(ObVal_GetValue(obval) < 0.) {
      arret("ObVal_Scan: negative value") ;
    }
  }
  
  
  line = String_GetAdvancedPosition ;
  
  
  /* Absolute or relative? */
  {
    char type[ObVal_MaxLengthOfKeyWord] ;
    int n = String_Scan(line," %s",type) ;
      
    /* Set default type */
    ObVal_SetTypeToAbsolute(obval) ;
          
    /* Read the type "absolute" */
    if(String_Is(type,"Absolute")) {
      ObVal_SetTypeToAbsolute(obval) ;
    }
          
    /* Read the type "relative" */
    if(String_Is(type,"Relative")) {
      ObVal_SetTypeToRelative(obval) ;
    }
  }
  
  
  /* Relaxation factor (if any) */
  {
    double r ;
    int n = String_FindAndScanExp(line,"Relax",","," = %le",&r) ;
    
    if(n) {
      ObVal_GetRelaxationFactor(obval) = r ;
    }
  }
}
