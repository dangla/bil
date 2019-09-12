#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Message.h"
#include "Mry.h"
#include "String.h"
#include "BCond.h"




BCond_t* BCond_New(void)
{
  BCond_t* bcond = (BCond_t*) Mry_New(BCond_t) ;
    
    
  /* Allocation of space for the name of unknowns */
  {
    char* name = (char*) Mry_New(char[BCond_MaxLengthOfKeyWord]) ;
  
    BCond_GetNameOfUnknown(bcond) = name ;
  }
    
    
  /* Allocation of space for the name of equations */
  {
    char* name = (char*) Mry_New(char[BCond_MaxLengthOfKeyWord]) ;
    
    BCond_GetNameOfEquation(bcond) = name ;
  }

  return(bcond) ;
}



void BCond_Scan(BCond_t* bcond,DataFile_t* datafile)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  
  /* Region */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Reg",","," = %d",&i) ;
    
    if(n) {
      BCond_GetRegionIndex(bcond) = i ;
    } else {
      arret("BCond_Scan: no region") ;
    }
  }
    
    
  /* Unknown */
  {
    char name[BCond_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Unk,Inc",","," = %s",name) ;
        
    if(n) {
      strcpy(BCond_GetNameOfUnknown(bcond),name) ;
    } else {
      arret("BCond_Scan: no unknown") ;
    }
      
    if(strlen(name) > BCond_MaxLengthOfKeyWord-1)  {
      arret("BCond_Scan: too long name of unknown") ;
    }
    
    if(isdigit(BCond_GetNameOfUnknown(bcond)[0])) {
      if(BCond_GetNameOfUnknown(bcond)[0] < '1') {
        arret("BCond_Scan: non positive unknown") ;
      }
    }
  }
    
    
  /* Field */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Field,Champ",","," = %d",&i) ;
        
    if(n) {
      BCond_GetFieldIndex(bcond) = i - 1 ;
    } else {
      arret("BCond_Scan: no field") ;
    }
  }
    
    
  /* Function */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Func,Fonc",","," = %d",&i) ;
        
    if(n) {
      BCond_GetFunctionIndex(bcond) = i - 1 ;
    } else {
      arret("BCond_Scan: no function") ;
    }
  }
    
    
  /* Equation (not mandatory) */
  {
    char name[BCond_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Equ",","," = %s",name) ;
      
    if(n) {
      strcpy(BCond_GetNameOfEquation(bcond),name) ;
    } else {
      strcpy(BCond_GetNameOfEquation(bcond)," ") ;
    }
      
    if(strlen(name) > BCond_MaxLengthOfKeyWord-1)  {
      arret("BCond_Scan: too long name of equation") ;
    }
        
    if(isdigit(BCond_GetNameOfEquation(bcond)[0])) {
      if(BCond_GetNameOfEquation(bcond)[0] < '1') {
        arret("BCond_Scan: non positive equation") ;
      }
    }
  }
}
