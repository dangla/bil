#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "InternationalSystemOfUnits.h"
#include "Mry.h"
#include "Units.h"


Units_t* (Units_New)(void)
{
  Units_t* units = (Units_t*) Mry_New(Units_t) ;
  
  
  Units_GetNbOfUnits(units) = 0 ;
  
  
  {
    Unit_t* unit = (Unit_t*) Mry_New(Unit_t[Units_MaxNbOfUnits]) ;
    int i ;
  
    for(i = 0 ; i < Units_MaxNbOfUnits ; i++) {
      Unit_t* u = Unit_New() ;
      
      unit[i] = u[0] ;
      free(u) ;
    }
    
    Units_GetUnit(units) = unit ;
  }
  
  return(units) ;
}



void (Units_Delete)(void* self)
{
  Units_t* units = (Units_t*) self ;

  {
    int n = Units_MaxNbOfUnits ;
    Unit_t* unit = Units_GetUnit(units) ;
    
    Mry_Delete(unit,n,Unit_Delete) ;
    free(unit) ;
    Units_GetUnit(units) = NULL ;
  }
}



#if 0
Units_t* (Units_Create)(DataFile_t* datafile)
{
  Units_t* units = (Units_t*) Mry_New(Units_t) ;
  
  
  {
    int n = DataFile_CountNbOfKeyWords(datafile,"UNITS,Units",",") ;
    int i ;


    DataFile_OpenFile(datafile,"r") ;
  
  
    for(i = 0 ; i < n ; i++) {
      int again = 1 ;
  
      DataFile_SetFilePositionAfterKey(datafile,"UNITS,Units",",",i + 1) ;
  
      Message_Direct("Enter in %s","Units") ;
      Message_Direct("\n") ;
  

      do {
        char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
        char  unit[Unit_MaxLengthOfKeyWord] ;

        if(!strncmp(line,"Length",1)) {
          char* pline = strchr(line,'=') + 1 ;
          
          sscanf(pline," %s",unit) ;
      
          InternationalSystemOfUnits_UseAsLength(unit) ;

        } else if(!strncmp(line,"Time",1)) {
          char* pline = strchr(line,'=') + 1 ;
          
          sscanf(pline," %s",unit) ;
      
          InternationalSystemOfUnits_UseAsTime(unit) ;

        } else if(!strncmp(line,"Mass",1)) {
          char* pline = strchr(line,'=') + 1 ;
          
          sscanf(pline," %s",unit) ;
      
          InternationalSystemOfUnits_UseAsMass(unit) ;
      
        } else {
          
          again = 0 ;
          
        }
        
      } while(again) ;
    }
  
    DataFile_CloseFile(datafile) ;
  }
  
  return(units) ;
}
#endif



#if 1
Units_t* (Units_Create)(DataFile_t* datafile)
{
  Units_t* units = (Units_t*) Units_New() ;
  
  
  {
    char* filecontent = DataFile_GetFileContent(datafile) ;
    char* c  = String_FindToken(filecontent,"UNITS,Units",",") ;
    
    if(!c) return(units) ;
  
      
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
  
    Message_Direct("Enter in %s","Units") ;
    Message_Direct("\n") ;
    
    {
      Unit_t* unit = Units_GetUnit(units) ;
      int n = 0 ;
      
      while(Unit_Scan(unit + n,datafile)) n++ ;
      
      if(n > Units_MaxNbOfUnits) {
        arret("Units_Create: too many units") ;
      }
    }
  
  }
  
  return(units) ;
}
#endif
