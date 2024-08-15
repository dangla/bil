#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Message.h"
#include "Mry.h"
#include "String_.h"
#include "Fields.h"
#include "Functions.h"
#include "ICond.h"



ICond_t* (ICond_New)(void)
{
  ICond_t* icond = (ICond_t*) Mry_New(ICond_t) ;
    
    
  /* Allocation of space for the name of unknowns */
  {
    char* name = (char*) Mry_New(char[ICond_MaxLengthOfKeyWord]) ;
  
    ICond_GetNameOfUnknown(icond) = name ;
  }
    
    
  /* Allocation of space for the name of files of nodal values */
  {
    char* filename = (char*) Mry_New(char[ICond_MaxLengthOfFileName]) ;
  
    ICond_GetFileNameOfNodalValues(icond) = filename ;
    ICond_GetFileNameOfNodalValues(icond)[0] = '\0' ;
  }
  
  
  /* Allocation of space for the region name */
  {
    char* name = (char*) Mry_New(char[ICond_MaxLengthOfRegionName]) ;
    
    ICond_GetRegionName(icond) = name ;
  }
  
  return(icond) ;
}



void (ICond_Delete)(void* self)
{
  ICond_t* icond = (ICond_t*) self ;
  
  {
    char* name = ICond_GetNameOfUnknown(icond) ;
    
    if(name) {
      free(name) ;
    }
  }
  
  {
    char* name = ICond_GetFileNameOfNodalValues(icond) ;
    
    if(name) {
      free(name) ;
    }
  }
  
  {
    char* name = ICond_GetRegionName(icond) ;
    
    if(name) {
      free(name) ;
    }
  }
}



void (ICond_Scan)(ICond_t* icond,DataFile_t* datafile)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  
  /* Region */
  {
    char name[ICond_MaxLengthOfRegionName] ;
    int n = String_FindAndScanExp(line,"Reg",","," = %s",name) ;
    //int i ;
    //int n = String_FindAndScanExp(line,"Reg",","," = %d",&i) ;
    
    if(n) {
      ICond_GetRegionTag(icond) = atoi(name) ;
      strncpy(ICond_GetRegionName(icond),name,ICond_MaxLengthOfRegionName)  ;
    } else {
      arret("ICond_Scan: no region") ;
    }
  }
    
    
  /* Unknown */
  {
    char name[ICond_MaxLengthOfKeyWord] ;
    int n = String_FindAndScanExp(line,"Unk,Inc",","," = %s",name) ;
        
    if(n) {
      strcpy(ICond_GetNameOfUnknown(icond),name) ;
    } else {
      arret("ICond_Scan: no unknown") ;
    }
      
    if(strlen(name) > ICond_MaxLengthOfKeyWord-1)  {
      arret("ICond_Scan: too long name of unknown") ;
    }
    
    if(isdigit(ICond_GetNameOfUnknown(icond)[0])) {
      if(ICond_GetNameOfUnknown(icond)[0] < '1') {
        arret("ICond_Scan: non positive unknown") ;
      }
    }
  }
    
    
  /* Field */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Field,Champ",","," = %d",&i) ;
    
    ICond_GetFieldIndex(icond) = -1 ;
    ICond_GetField(icond) = NULL ;
        
    if(n) {
      Fields_t* fields = ICond_GetFields(icond) ;
      int  n_fields = Fields_GetNbOfFields(fields) ;
      int ifld = i - 1 ;
      
      ICond_GetFieldIndex(icond) = ifld ;
      
      if(ifld < 0) {
        
        ICond_GetField(icond) = NULL ;
        
      } else if(ifld < n_fields) {
        Field_t* field = Fields_GetField(fields) ;
        
        ICond_GetField(icond) = field + ifld ;
          
      } else {
        
        arret("ICond_Scan: field out of range") ;
          
      }
    }
  }
    
    
  /* File */
  {
    char name[ICond_MaxLengthOfFileName] ;
    int n = String_FindAndScanExp(line,"File,Fichier",","," = %s",name) ;
        
    if(n) {
      
      if(strlen(name) > ICond_MaxLengthOfFileName-1)  {
        arret("ICond_Scan: name too long") ;
      }
      
      strcpy(ICond_GetFileNameOfNodalValues(icond),name) ;
    }
  }
    
    
  /* Function (not mandatory) */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Func,Fonc",","," = %d",&i) ;
    
    ICond_GetFunctionIndex(icond) = -1 ;
    ICond_GetFunction(icond) = NULL ;
        
    if(n) {
      Functions_t* functions = ICond_GetFunctions(icond) ;
      int n_functions = Functions_GetNbOfFunctions(functions) ;
      int ifct = i - 1 ;
      
      ICond_GetFunctionIndex(icond) = ifct ;
      
      if(ifct < 0) {
        
        ICond_GetFunction(icond) = NULL ;
        
      } else if(ifct < n_functions) {
        Function_t* fct = Functions_GetFunction(functions) ;
        
        ICond_GetFunction(icond) = fct + ifct ;
        
      } else {
        
        arret("ICond_Scan: function out of range") ;
        
      }
    }
  }
}
