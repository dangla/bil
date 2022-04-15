#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "InternationalSystemOfUnits.h"
#include "Mry.h"
#include "String_.h"
#include "Unit.h"




Unit_t* (Unit_New)(void)
{
  Unit_t* unit = (Unit_t*) Mry_New(Unit_t) ;
    
    
  /* Allocation of space for the name of unit */
  {
    char* name = (char*) Mry_New(char[Unit_MaxLengthOfKeyWord]) ;
  
    Unit_GetName(unit) = name ;
  }
  
  return(unit) ;
}



void (Unit_Delete)(void* self)
{
  Unit_t* unit = (Unit_t*) self ;
  
  {
    char* name = Unit_GetName(unit) ;
    
    if(name) {
      free(name) ;
      Unit_GetName(unit) = NULL ;
    }
  }
}



int (Unit_Scan)(Unit_t* unit,DataFile_t* datafile)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;

  
    
  /* Length */
  {
    char name[Unit_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Length",","," = %s",name) ;
        
    if(n) {
      strcpy(Unit_GetName(unit),name) ;
      
      InternationalSystemOfUnits_UseAsLength(name) ;
      
      return(1) ;
    }
  }
  
  
  
  /* Time */
  {
    char name[Unit_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Time",","," = %s",name) ;
        
    if(n) {
      strcpy(Unit_GetName(unit),name) ;
      
      InternationalSystemOfUnits_UseAsTime(name) ;
      
      return(1) ;
    }
  }
  
  
  
  /* Mass */
  {
    char name[Unit_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Mass",","," = %s",name) ;
        
    if(n) {
      strcpy(Unit_GetName(unit),name) ;
      
      InternationalSystemOfUnits_UseAsMass(name) ;
      
      return(1) ;
    }
  }
  
  return(0) ;
}
