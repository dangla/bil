#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "InternationalSystemOfUnits.h"
#include "Units.h"



Units_t* (Units_Create)(DataFile_t* datafile)
{
  Units_t* units = (Units_t*) malloc(sizeof(Units_t)) ;
  
  if(!units) arret("Units_Create") ;
  
  
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
