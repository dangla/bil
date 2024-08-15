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

