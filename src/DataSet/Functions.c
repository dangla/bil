#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "String_.h"
#include "Mry.h"
#include "Functions.h"



Functions_t* (Functions_New)(const int n_fncts)
{
  Functions_t* functions   = (Functions_t*) Mry_New(Functions_t) ;
  
  
  Functions_GetNbOfFunctions(functions) = n_fncts ;
  Functions_GetFunction(functions) = NULL ;

  if(n_fncts > 0) {
    Function_t* function = (Function_t*) Mry_New(Function_t[n_fncts]) ;
    
    Functions_GetFunction(functions) = function ;
  }
  
  return(functions) ;
}



Functions_t* (Functions_Create)(DataFile_t* datafile)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"FONC,FUNC,Functions",",") ;
  int n_fncts = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Functions_t* functions = Functions_New(n_fncts) ;
  
  
  Message_Direct("Enter in %s","Functions") ;
  Message_Direct("\n") ;
  
  
  if(n_fncts <= 0) {
    return(functions) ;
  }


  {
    int i_fn ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i_fn = 0 ; i_fn < n_fncts ; i_fn++) {
      Function_t* function = Functions_GetFunction(functions) + i_fn ;
  
  
      Message_Direct("Enter in %s %d","Function",i_fn + 1) ;
      Message_Direct("\n") ;
      
      i_fn += Function_Scan(function,datafile) - 1 ;
      
      if(i_fn >= n_fncts) {
        arret("Functions_Create: too many functions") ;
      }
    }
  }
  
  return(functions) ;
}



void (Functions_Delete)(void* self)
{
  Functions_t* functions = (Functions_t*) self ;
  
  {
    int n_functions = Functions_GetNbOfFunctions(functions) ;
    
    if(n_functions > 0) {
      Function_t* function = Functions_GetFunction(functions) ;
      
      if(function) {
        int i ;
      
        for(i = 0 ; i < n_functions ; i++) {
          Function_t* functioni = function + i ;
        
          Function_Delete(functioni) ;
        }
      
        free(function) ;
      }
    }
  }
}
